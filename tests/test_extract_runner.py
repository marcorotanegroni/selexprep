"""End-to-end tests for ``selexprep.extract.runner.run_extract``."""

from __future__ import annotations

import gzip
import json
import shutil
from pathlib import Path

import pytest

from selexprep.extract.runner import run_extract
from selexprep.library.adapters import reverse_complement
from selexprep.library.report import LibraryReport

PRIMER_5P = "GGTAATACGACTCACTATAGGG"  # T7, 22 nt
PRIMER_3P = "CCATGCATGCATGCATGCAT"  # 20 nt


pytestmark = pytest.mark.skipif(shutil.which("cutadapt") is None, reason="cutadapt not on PATH")


# ---------------------------------------------------------------------------
# Synthetic FASTQ + LibraryReport builders
# ---------------------------------------------------------------------------


def _write_fastq_gz(path: Path, seqs: list[str]) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        for i, s in enumerate(seqs):
            fh.write(f"@read_{i}\n{s}\n+\n{'I' * len(s)}\n")


def _make_round_fastq(
    tmp_path: Path,
    round_no: int,
    primer_5p: str | None,
    primer_3p: str | None,
    *,
    n: int = 20,
    random_len: int = 30,
) -> Path:
    bases = "ACGT"
    seqs = []
    for i in range(n):
        rand = "".join(bases[(i * 7 + j * 13) % 4] for j in range(random_len))
        seqs.append((primer_5p or "") + rand + (primer_3p or ""))
    path = tmp_path / f"round_{round_no:02d}.fastq.gz"
    _write_fastq_gz(path, seqs)
    return path


def _make_library_report(
    *,
    primer_5p: str | None = PRIMER_5P,
    primer_3p: str | None = PRIMER_3P,
    extraction_mode: str = "BOTH_PRIMERS_SINGLE_READ",
    full_insert_recovered: bool = True,
    required_action: str = "NONE",
    orientation: str = "FORWARD",
    status: str = "HIGH",
    failure_reason: str | None = None,
) -> LibraryReport:
    """Construct a minimal LibraryReport for runner tests.

    Lots of fields are set to neutral defaults; the test only asserts on
    behavior driven by the fields that vary across scenarios.
    """
    return LibraryReport(
        primer_5p=primer_5p,
        primer_3p=primer_3p,
        variants_5p=[],
        variants_3p=[],
        known_adapter_hits={},
        extraction_mode=extraction_mode,  # type: ignore[arg-type]
        full_insert_recovered=full_insert_recovered,
        read_source="R1",
        required_action=required_action,  # type: ignore[arg-type]
        orientation=orientation,  # type: ignore[arg-type]
        n_length_mode=30,
        n_length_distribution={30: 100},
        n_length_confidence=1.0,
        match_rate_5p=0.95,
        match_rate_3p=0.95,
        position_consistency_5p=0.95,
        position_consistency_3p=0.95,
        read_fraction_used_for_inference=1.0,
        sampling_seed=42,
        confidence=0.85,
        status=status,  # type: ignore[arg-type]
        failure_reason=failure_reason,
    )


# ---------------------------------------------------------------------------
# Per-extraction_mode happy paths
# ---------------------------------------------------------------------------


def test_run_extract_both_primers_single_read(tmp_path: Path) -> None:
    fq = _make_round_fastq(tmp_path, 0, PRIMER_5P, PRIMER_3P)
    lr = _make_library_report(extraction_mode="BOTH_PRIMERS_SINGLE_READ")
    outdir = tmp_path / "out"

    result = run_extract(lr, [fq], outdir, round_map={fq.name: 0})

    assert not result.skipped
    expected = outdir / "round_00" / "extracted.fasta.gz"
    assert expected.exists()
    assert any(p == expected for p in result.outputs)


def test_run_extract_five_prime_only(tmp_path: Path) -> None:
    fq = _make_round_fastq(tmp_path, 0, PRIMER_5P, primer_3p=None)
    lr = _make_library_report(
        primer_3p=None,
        extraction_mode="FIVE_PRIME_ONLY",
        full_insert_recovered=False,
    )
    outdir = tmp_path / "out"

    result = run_extract(lr, [fq], outdir, round_map={fq.name: 0})

    assert not result.skipped
    expected = outdir / "round_00" / "partial_5p_extracted.fasta.gz"
    assert expected.exists()


def test_run_extract_three_prime_only(tmp_path: Path) -> None:
    fq = _make_round_fastq(tmp_path, 0, primer_5p=None, primer_3p=PRIMER_3P)
    lr = _make_library_report(
        primer_5p=None,
        extraction_mode="THREE_PRIME_ONLY",
        full_insert_recovered=False,
    )
    outdir = tmp_path / "out"

    result = run_extract(lr, [fq], outdir, round_map={fq.name: 0})

    assert not result.skipped
    expected = outdir / "round_00" / "partial_3p_extracted.fasta.gz"
    assert expected.exists()


def test_run_extract_paired_end_split_primers(tmp_path: Path) -> None:
    r1 = _make_round_fastq(tmp_path, 0, PRIMER_5P, primer_3p=None, random_len=80)
    r2_5p = reverse_complement(PRIMER_3P)
    r2_path = tmp_path / "round_00_R2.fastq.gz"
    bases = "ACGT"
    seqs = []
    for i in range(20):
        rand = "".join(bases[(i * 11 + j * 17) % 4] for j in range(80))
        seqs.append(r2_5p + rand)
    _write_fastq_gz(r2_path, seqs)

    lr = _make_library_report(
        extraction_mode="PAIRED_END_SPLIT_PRIMERS",
        full_insert_recovered=False,
        required_action="READ_MERGING_RECOMMENDED",
    )
    outdir = tmp_path / "out"

    result = run_extract(
        lr,
        [r1],
        outdir,
        round_map={r1.name: 0, r2_path.name: 0},
        paired_r2_inputs={0: [r2_path]},
    )

    assert not result.skipped
    r1_out = outdir / "round_00" / "partial_5p_extracted_R1.fasta.gz"
    r2_out = outdir / "round_00" / "partial_3p_extracted_R2.fasta.gz"
    assert r1_out.exists() and r2_out.exists()


# ---------------------------------------------------------------------------
# Refusal paths
# ---------------------------------------------------------------------------


def test_run_extract_refuses_unable_status(tmp_path: Path) -> None:
    fq = _make_round_fastq(tmp_path, 0, PRIMER_5P, PRIMER_3P)
    lr = _make_library_report(
        status="UNABLE_TO_INFER",
        extraction_mode="UNABLE_TO_EXTRACT",
        failure_reason="Both match rates below threshold",
    )
    outdir = tmp_path / "out"

    result = run_extract(lr, [fq], outdir, round_map={fq.name: 0})

    assert result.skipped
    assert result.skipped_reason is not None
    # Output FASTAs are NOT produced
    assert not (outdir / "round_00" / "extracted.fasta.gz").exists()


def test_run_extract_refuses_unable_extraction_mode(tmp_path: Path) -> None:
    fq = _make_round_fastq(tmp_path, 0, PRIMER_5P, PRIMER_3P)
    lr = _make_library_report(
        status="LOW",  # status OK, but extraction_mode says no
        extraction_mode="UNABLE_TO_EXTRACT",
    )
    outdir = tmp_path / "out"

    result = run_extract(lr, [fq], outdir, round_map={fq.name: 0})

    assert result.skipped


# ---------------------------------------------------------------------------
# No-clobber behavior
# ---------------------------------------------------------------------------


def test_run_extract_no_clobber_without_rebuild(tmp_path: Path) -> None:
    fq = _make_round_fastq(tmp_path, 0, PRIMER_5P, PRIMER_3P)
    lr = _make_library_report()
    outdir = tmp_path / "out"

    # First run succeeds.
    run_extract(lr, [fq], outdir, round_map={fq.name: 0})

    # Second run without --rebuild raises.
    with pytest.raises(FileExistsError, match="Pass --rebuild"):
        run_extract(lr, [fq], outdir, round_map={fq.name: 0})


def test_run_extract_rebuild_overwrites(tmp_path: Path) -> None:
    fq = _make_round_fastq(tmp_path, 0, PRIMER_5P, PRIMER_3P)
    lr = _make_library_report()
    outdir = tmp_path / "out"

    run_extract(lr, [fq], outdir, round_map={fq.name: 0})
    # Should not raise.
    result = run_extract(lr, [fq], outdir, round_map={fq.name: 0}, rebuild=True)
    assert not result.skipped


# ---------------------------------------------------------------------------
# Strand report emission
# ---------------------------------------------------------------------------


def test_run_extract_emits_strand_report_for_mixed_orientation(tmp_path: Path) -> None:
    fq = _make_round_fastq(tmp_path, 0, PRIMER_5P, PRIMER_3P)
    lr = _make_library_report(orientation="MIXED")
    outdir = tmp_path / "out"

    result = run_extract(lr, [fq], outdir, round_map={fq.name: 0})

    assert result.strand_report_path is not None
    assert result.strand_report_path.exists()
    assert "strand_report.tsv" in result.strand_report_path.name


def test_run_extract_emits_strand_report_for_reverse_orientation(tmp_path: Path) -> None:
    # All-reverse pool
    rc_5p = reverse_complement(PRIMER_5P)
    rc_3p = reverse_complement(PRIMER_3P)
    fq = _make_round_fastq(tmp_path, 0, rc_3p, rc_5p)
    lr = _make_library_report(orientation="REVERSE")
    outdir = tmp_path / "out"

    result = run_extract(lr, [fq], outdir, round_map={fq.name: 0})

    assert result.strand_report_path is not None
    assert result.strand_report_path.exists()


def test_run_extract_no_strand_report_for_forward_orientation(tmp_path: Path) -> None:
    fq = _make_round_fastq(tmp_path, 0, PRIMER_5P, PRIMER_3P)
    lr = _make_library_report(orientation="FORWARD")
    outdir = tmp_path / "out"

    result = run_extract(lr, [fq], outdir, round_map={fq.name: 0})

    assert result.strand_report_path is None


# ---------------------------------------------------------------------------
# Trim reports JSON
# ---------------------------------------------------------------------------


def test_run_extract_emits_trim_reports_json(tmp_path: Path) -> None:
    fq = _make_round_fastq(tmp_path, 0, PRIMER_5P, PRIMER_3P)
    lr = _make_library_report()
    outdir = tmp_path / "out"

    result = run_extract(lr, [fq], outdir, round_map={fq.name: 0})

    assert result.trim_reports_path is not None
    assert result.trim_reports_path.exists()
    payload = json.loads(result.trim_reports_path.read_text())
    assert isinstance(payload, list) and len(payload) >= 1
    # Each entry contains the exact cutadapt argv used.
    assert "cutadapt_cmd" in payload[0]
    assert "n_in" in payload[0]
    assert "n_out" in payload[0]


# ---------------------------------------------------------------------------
# Multi-round happy path
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Codex pass 1 regressions (2026-05-21)
# ---------------------------------------------------------------------------


def test_run_extract_multi_fastq_same_round_aggregates(tmp_path: Path) -> None:
    """Codex Phase 3 pass 1 regression: when a round has multiple input
    FASTQs, all reads must end up in the per-round output (concatenated
    deterministically). Previously each iteration overwrote the same
    target file with the latest input → only the last input's reads
    survived."""
    # Two single-end FASTQs in the SAME round (round 0), different read
    # counts so we can detect aggregation by counting the output.
    n_a = 3
    n_b = 7
    fq_a = tmp_path / "round0_part_a.fastq.gz"
    fq_b = tmp_path / "round0_part_b.fastq.gz"
    bases = "ACGT"
    seqs_a = [
        PRIMER_5P + "".join(bases[(i * 7 + j * 13) % 4] for j in range(30)) + PRIMER_3P
        for i in range(n_a)
    ]
    seqs_b = [
        PRIMER_5P + "".join(bases[((i + 100) * 7 + j * 13) % 4] for j in range(30)) + PRIMER_3P
        for i in range(n_b)
    ]
    _write_fastq_gz(fq_a, seqs_a)
    _write_fastq_gz(fq_b, seqs_b)

    lr = _make_library_report()
    outdir = tmp_path / "out"
    result = run_extract(
        lr,
        [fq_a, fq_b],
        outdir,
        round_map={fq_a.name: 0, fq_b.name: 0},
    )

    assert not result.skipped, result.skipped_reason

    extracted = outdir / "round_00" / "extracted.fasta.gz"
    assert extracted.exists()

    with gzip.open(extracted, "rt", encoding="utf-8") as fh:
        n_records = sum(1 for line in fh if line.startswith(">"))

    # The bug would have produced n_b (7) — last input overwrote first.
    # The fix concatenates both, so we expect n_a + n_b = 10.
    assert n_records == n_a + n_b, (
        f"expected {n_a + n_b} aggregated reads (Codex Phase 3 pass 1 regression), got {n_records}"
    )


def test_run_extract_sample_sheet_paired_end_demux_rebuilds_r2_inputs(tmp_path: Path) -> None:
    """Codex Phase 3 pass 1 regression: when --sample-sheet is given and the
    LR is PAIRED_END_SPLIT_PRIMERS, the runner must rebuild paired_r2_inputs
    from the demuxed ``_2.fastq.gz`` files. Previously only the R1 files
    were collected and paired_r2_inputs stayed at the caller's value
    (typically None), so the dispatch refused with 'requires --paired-r2'."""
    from unittest.mock import patch

    # Stage a pre-demuxed directory layout that mirrors what
    # ``demux_sample_sheet`` would produce: per-round folders with
    # ``_1.fastq.gz`` (R1) and ``_2.fastq.gz`` (R2).
    outdir = tmp_path / "out"
    demux_dir = outdir / "demux"
    for r in (0, 1):
        round_dir = demux_dir / f"round_{r:02d}"
        round_dir.mkdir(parents=True)
        # R1: primer_5p + random
        r1_seqs = [
            PRIMER_5P + "ACGTACGTACGTACGTACGTACGTACGTACGT"
            "ACGTACGTACGTACGTACGTACGTACGTACGT"
            "ACGTACGTACGTACGT"
            for _ in range(20)
        ]
        _write_fastq_gz(round_dir / "srr_1.fastq.gz", r1_seqs)
        # R2: revcomp(primer_3p) + random
        rc_3p = reverse_complement(PRIMER_3P)
        r2_seqs = [
            rc_3p + "ACGTACGTACGTACGTACGTACGTACGTACGT"
            "ACGTACGTACGTACGTACGTACGTACGTACGT"
            "ACGTACGTACGTACGTACGT"
            for _ in range(20)
        ]
        _write_fastq_gz(round_dir / "srr_2.fastq.gz", r2_seqs)

    # Sample sheet path is needed (even though we mock the demux call).
    sample_sheet = tmp_path / "samples.tsv"
    sample_sheet.write_text("dummy\n", encoding="utf-8")

    lr = _make_library_report(
        extraction_mode="PAIRED_END_SPLIT_PRIMERS",
        full_insert_recovered=False,
        required_action="READ_MERGING_RECOMMENDED",
    )

    # Bypass the actual demux (the test stages its output directly).
    with patch("selexprep.extract.runner.demux_sample_sheet"):
        result = run_extract(
            lr,
            [],  # ignored — sample-sheet rebuilds fastq_inputs
            outdir,
            round_map={},  # ignored — rebuilt from demuxed paths
            sample_sheet=sample_sheet,
        )

    assert not result.skipped, result.skipped_reason
    # Both rounds should have produced R1 + R2 outputs.
    for r in (0, 1):
        assert (outdir / f"round_{r:02d}" / "partial_5p_extracted_R1.fasta.gz").exists()
        assert (outdir / f"round_{r:02d}" / "partial_3p_extracted_R2.fasta.gz").exists()


def test_run_extract_sample_sheet_input_sha256_distinct_per_round(tmp_path: Path) -> None:
    """Codex Phase 3 pass 1 follow-up: in sample-sheet mode the demuxed
    inputs share basenames across rounds (``srr_1.fastq.gz`` in every
    ``round_NN/`` folder). The manifest's ``input_sha256`` must key by
    path relative to the demux dir to keep both rounds distinct;
    basename-only keying would collapse the second round onto the first."""
    from unittest.mock import patch

    from selexprep.manifest import read_manifest_json

    outdir = tmp_path / "out"
    demux_dir = outdir / "demux"

    # Stage pre-demuxed inputs (single-end, two rounds, same basename).
    for r in (0, 1):
        round_dir = demux_dir / f"round_{r:02d}"
        round_dir.mkdir(parents=True)
        # Distinct content per round so the sha256 values themselves differ.
        seqs = [PRIMER_5P + ("A" if r == 0 else "C") * 30 + PRIMER_3P for _ in range(20)]
        _write_fastq_gz(round_dir / "srr_1.fastq.gz", seqs)

    sample_sheet = tmp_path / "samples.tsv"
    sample_sheet.write_text("dummy\n", encoding="utf-8")

    lr = _make_library_report()

    with patch("selexprep.extract.runner.demux_sample_sheet"):
        result = run_extract(lr, [], outdir, round_map={}, sample_sheet=sample_sheet)

    assert not result.skipped, result.skipped_reason
    assert result.manifest_path is not None

    m = read_manifest_json(result.manifest_path)
    # Both rounds must appear with distinct keys + distinct hashes.
    assert "round_00/srr_1.fastq.gz" in m.input_sha256
    assert "round_01/srr_1.fastq.gz" in m.input_sha256
    assert m.input_sha256["round_00/srr_1.fastq.gz"] != m.input_sha256["round_01/srr_1.fastq.gz"]


def test_run_extract_multi_round_emits_one_file_per_round(tmp_path: Path) -> None:
    fq0 = _make_round_fastq(tmp_path, 0, PRIMER_5P, PRIMER_3P)
    fq1 = _make_round_fastq(tmp_path, 1, PRIMER_5P, PRIMER_3P)
    fq2 = _make_round_fastq(tmp_path, 2, PRIMER_5P, PRIMER_3P)
    lr = _make_library_report()
    outdir = tmp_path / "out"

    result = run_extract(
        lr,
        [fq0, fq1, fq2],
        outdir,
        round_map={fq0.name: 0, fq1.name: 1, fq2.name: 2},
    )

    assert not result.skipped
    for r in range(3):
        assert (outdir / f"round_{r:02d}" / "extracted.fasta.gz").exists()
