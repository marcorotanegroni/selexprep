"""Tests for Phase 4 override-primer + extract_diff.tsv path."""

from __future__ import annotations

import gzip
import shutil
from pathlib import Path

import pytest

from selexprep.extract.runner import run_extract
from selexprep.library.report import LibraryReport

PRIMER_5P = "GGTAATACGACTCACTATAGGG"
PRIMER_3P = "CCATGCATGCATGCATGCAT"


pytestmark = pytest.mark.skipif(shutil.which("cutadapt") is None, reason="cutadapt not on PATH")


def _write_fastq_gz(path: Path, seqs: list[str]) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        for i, s in enumerate(seqs):
            fh.write(f"@read_{i}\n{s}\n+\n{'I' * len(s)}\n")


def _make_round_fastq(tmp_path: Path, round_no: int, n: int = 20) -> Path:
    bases = "ACGT"
    seqs = []
    for i in range(n):
        rand = "".join(bases[(i * 7 + j * 13) % 4] for j in range(30))
        seqs.append(PRIMER_5P + rand + PRIMER_3P)
    path = tmp_path / f"round_{round_no:02d}.fastq.gz"
    _write_fastq_gz(path, seqs)
    return path


def _make_library_report() -> LibraryReport:
    return LibraryReport(
        primer_5p=PRIMER_5P,
        primer_3p=PRIMER_3P,
        variants_5p=[],
        variants_3p=[],
        known_adapter_hits={},
        extraction_mode="BOTH_PRIMERS_SINGLE_READ",
        full_insert_recovered=True,
        read_source="R1",
        required_action="NONE",
        orientation="FORWARD",
        n_length_mode=30,
        n_length_distribution={30: 20},
        n_length_confidence=1.0,
        match_rate_5p=0.95,
        match_rate_3p=0.95,
        position_consistency_5p=0.95,
        position_consistency_3p=0.95,
        read_fraction_used_for_inference=1.0,
        sampling_seed=42,
        confidence=0.85,
        status="HIGH",
        failure_reason=None,
    )


# ---------------------------------------------------------------------------
# Override without --rebuild -> outdir/overridden/ subtree
# ---------------------------------------------------------------------------


def test_override_5p_without_rebuild_writes_overridden_subtree(tmp_path: Path) -> None:
    fq = _make_round_fastq(tmp_path, 0)
    lr = _make_library_report()
    outdir = tmp_path / "out"

    # Baseline run first so we have a baseline outdir to leave alone.
    run_extract(lr, [fq], outdir, round_map={fq.name: 0})

    # Override run — no --rebuild.
    result = run_extract(
        lr,
        [fq],
        outdir,
        round_map={fq.name: 0},
        override_primer_5p="AAAAAAAAAAAAAAA",  # junk; cutadapt will drop reads
    )

    assert not result.skipped
    # Baseline outputs untouched.
    assert (outdir / "round_00" / "extracted.fasta.gz").exists()
    assert (outdir / "selexprep_manifest.json").exists()
    # Override outputs in subdir.
    assert (outdir / "overridden" / "round_00" / "extracted.fasta.gz").exists()
    assert (outdir / "overridden" / "selexprep_manifest.json").exists()
    # No extract_diff.tsv emitted in this mode (override but no --rebuild).
    assert result.extract_diff_path is None


def test_override_3p_without_rebuild_writes_overridden_subtree(tmp_path: Path) -> None:
    fq = _make_round_fastq(tmp_path, 0)
    lr = _make_library_report()
    outdir = tmp_path / "out"

    run_extract(lr, [fq], outdir, round_map={fq.name: 0})

    result = run_extract(
        lr,
        [fq],
        outdir,
        round_map={fq.name: 0},
        override_primer_3p="TTTTTTTTTTTTTTT",
    )

    assert not result.skipped
    assert (outdir / "overridden" / "round_00" / "extracted.fasta.gz").exists()


# ---------------------------------------------------------------------------
# Override + --rebuild -> overwrite baseline, emit extract_diff.tsv
# ---------------------------------------------------------------------------


def test_override_with_rebuild_overwrites_and_emits_diff(tmp_path: Path) -> None:
    fq = _make_round_fastq(tmp_path, 0)
    lr = _make_library_report()
    outdir = tmp_path / "out"

    # Baseline.
    run_extract(lr, [fq], outdir, round_map={fq.name: 0})
    assert (outdir / "selexprep_manifest.json").exists()

    # Override + rebuild: overwrites baseline, emits diff.
    result = run_extract(
        lr,
        [fq],
        outdir,
        round_map={fq.name: 0},
        rebuild=True,
        override_primer_5p="AAAAAAAAAAAAAAA",
    )

    assert not result.skipped
    assert result.extract_diff_path is not None
    assert result.extract_diff_path.exists()
    # NOT a subtree — overwrote in place.
    assert not (outdir / "overridden").exists()


def test_extract_diff_tsv_has_correct_columns(tmp_path: Path) -> None:
    fq = _make_round_fastq(tmp_path, 0)
    lr = _make_library_report()
    outdir = tmp_path / "out"

    run_extract(lr, [fq], outdir, round_map={fq.name: 0})
    result = run_extract(
        lr,
        [fq],
        outdir,
        round_map={fq.name: 0},
        rebuild=True,
        override_primer_5p="AAAAAAAAAAAAAAA",
    )

    diff_text = result.extract_diff_path.read_text(encoding="utf-8")  # type: ignore[union-attr]
    lines = diff_text.splitlines()
    header = lines[0].split("\t")
    assert header == [
        "round",
        "primer_5p_baseline",
        "primer_5p_new",
        "primer_3p_baseline",
        "primer_3p_new",
        "n_in",
        "n_out_baseline",
        "n_out_new",
        "delta_n_out",
    ]
    # One data row for round 0.
    row = lines[1].split("\t")
    assert row[0] == "0"
    assert row[1] == PRIMER_5P  # baseline 5p
    assert row[2] == "AAAAAAAAAAAAAAA"  # new 5p


# ---------------------------------------------------------------------------
# --rebuild alone (no override) -> NO extract_diff
# ---------------------------------------------------------------------------


def test_rebuild_without_override_does_not_emit_diff(tmp_path: Path) -> None:
    fq = _make_round_fastq(tmp_path, 0)
    lr = _make_library_report()
    outdir = tmp_path / "out"

    run_extract(lr, [fq], outdir, round_map={fq.name: 0})
    result = run_extract(lr, [fq], outdir, round_map={fq.name: 0}, rebuild=True)

    assert not result.skipped
    assert result.extract_diff_path is None
    assert not (outdir / "extract_diff.tsv").exists()


# ---------------------------------------------------------------------------
# Manifest emission
# ---------------------------------------------------------------------------


def test_run_extract_emits_manifest(tmp_path: Path) -> None:
    """Every run (with or without override) writes a selexprep_manifest.json."""
    fq = _make_round_fastq(tmp_path, 0)
    lr = _make_library_report()
    outdir = tmp_path / "out"

    result = run_extract(lr, [fq], outdir, round_map={fq.name: 0})

    assert result.manifest_path is not None
    assert result.manifest_path.exists()


def test_manifest_records_override_primers_when_applied(tmp_path: Path) -> None:
    fq = _make_round_fastq(tmp_path, 0)
    lr = _make_library_report()
    outdir = tmp_path / "out"

    override_5p = "ACGTACGTACGTACG"  # valid IUPAC, distinct from PRIMER_5P
    result = run_extract(
        lr,
        [fq],
        outdir,
        round_map={fq.name: 0},
        override_primer_5p=override_5p,
    )

    # Read the override manifest from the overridden/ subtree.
    from selexprep.manifest import read_manifest_json

    assert result.manifest_path is not None
    m = read_manifest_json(result.manifest_path)
    assert m.library_report.primer_5p == override_5p
    assert m.library_report.primer_3p == PRIMER_3P  # 3p unchanged
