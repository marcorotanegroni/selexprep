"""Unit tests for selexprep.extract.demux."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from selexprep.extract.demux import (
    demux_fastq,
    demux_sample_sheet,
    read_sample_sheet,
    validate_barcodes,
)


def _write_fastq_gz(path: Path, sequences: list[str]) -> None:
    with gzip.open(path, "wt") as fh:
        for i, seq in enumerate(sequences):
            fh.write(f"@read_{i}\n{seq}\n+\n{'I' * len(seq)}\n")


def _read_fastq_gz_seqs(path: Path) -> list[str]:
    if not path.exists():
        return []
    seqs = []
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                seqs.append(line.strip())
    return seqs


# ----- validate_barcodes -----


def test_validate_barcodes_passes_distant_pair() -> None:
    validate_barcodes({"AAAAA": 1, "TTTTT": 2}, max_mismatches=1)


def test_validate_barcodes_rejects_too_close() -> None:
    # Hamming distance 1, but max_mm=1 requires ≥3
    with pytest.raises(ValueError, match="Hamming-distance"):
        validate_barcodes({"AAAAA": 1, "AAAAT": 2}, max_mismatches=1)


def test_validate_barcodes_rejects_unequal_lengths() -> None:
    with pytest.raises(ValueError, match="different lengths"):
        validate_barcodes({"AAAA": 1, "AAAAA": 2}, max_mismatches=1)


# ----- demux_fastq (single-end) -----


def test_demux_fastq_single_end_routes_by_5p_barcode(tmp_path: Path) -> None:
    fq = tmp_path / "SRR1.fastq.gz"
    _write_fastq_gz(
        fq,
        [
            "AAAAACCCCC",  # barcode AAAAA → round 1, payload CCCCC
            "AAAAAGGGGG",  # barcode AAAAA → round 1, payload GGGGG
            "TTTTTGGGGG",  # barcode TTTTT → round 2, payload GGGGG
            "GGGGGAAAAA",  # no barcode match → unassigned
        ],
    )
    out = tmp_path / "out"
    report = demux_fastq(r1_path=fq, out_dir=out, barcodes={"AAAAA": 1, "TTTTT": 2}, srr="SRR1")

    assert report.total_reads == 4
    assert report.per_round == {1: 2, 2: 1}
    assert report.unassigned_reads == 1
    assert report.unassigned_fraction == 0.25

    r1_seqs = _read_fastq_gz_seqs(out / "round_01" / "SRR1.fastq.gz")
    assert r1_seqs == ["CCCCC", "GGGGG"]  # barcode stripped by default
    r2_seqs = _read_fastq_gz_seqs(out / "round_02" / "SRR1.fastq.gz")
    assert r2_seqs == ["GGGGG"]
    un_seqs = _read_fastq_gz_seqs(out / "unassigned" / "SRR1.fastq.gz")
    assert un_seqs == ["GGGGGAAAAA"]


def test_demux_fastq_preserves_barcode_when_trim_false(tmp_path: Path) -> None:
    fq = tmp_path / "SRR2.fastq.gz"
    _write_fastq_gz(fq, ["AAAAACCCCC"])
    out = tmp_path / "out"
    demux_fastq(
        r1_path=fq,
        out_dir=out,
        barcodes={"AAAAA": 1, "TTTTT": 2},
        srr="SRR2",
        trim_barcode=False,
    )
    seqs = _read_fastq_gz_seqs(out / "round_01" / "SRR2.fastq.gz")
    assert seqs == ["AAAAACCCCC"]


# ----- demux_fastq (paired-end) -----


def test_demux_fastq_paired_end_keeps_pair_sync(tmp_path: Path) -> None:
    r1 = tmp_path / "SRR3_1.fastq.gz"
    r2 = tmp_path / "SRR3_2.fastq.gz"
    _write_fastq_gz(r1, ["AAAAACCCCC", "TTTTTGGGGG"])
    _write_fastq_gz(r2, ["TTTTTGGGGG", "AAAAACCCCC"])  # different content; same record count
    out = tmp_path / "out"

    report = demux_fastq(
        r1_path=r1,
        out_dir=out,
        r2_path=r2,
        barcodes={"AAAAA": 1, "TTTTT": 2},
        srr="SRR3",
    )

    assert report.paired is True
    assert report.per_round == {1: 1, 2: 1}

    # Round 1: R1 had barcode AAAAA, payload CCCCC (stripped); R2 mate is TTTTTGGGGG
    r1_round1 = _read_fastq_gz_seqs(out / "round_01" / "SRR3_1.fastq.gz")
    r2_round1 = _read_fastq_gz_seqs(out / "round_01" / "SRR3_2.fastq.gz")
    assert r1_round1 == ["CCCCC"]
    assert r2_round1 == ["TTTTTGGGGG"]  # R2 is never trimmed


# ----- sample sheet -----


def test_read_sample_sheet_groups_barcodes_per_srr(tmp_path: Path) -> None:
    fq = tmp_path / "data.fastq.gz"
    _write_fastq_gz(fq, [])

    sheet = tmp_path / "sheet.tsv"
    sheet.write_text(
        "srr\tr1_path\tr2_path\tround\tbarcode\n"
        f"SRR1\t{fq}\t\t1\tAAAAA\n"
        f"SRR1\t{fq}\t\t2\tTTTTT\n"
        f"SRR2\t{fq}\t\t1\tGGGGG\n",
        encoding="utf-8",
    )

    jobs = read_sample_sheet(sheet)
    assert len(jobs) == 2
    srr1 = next(j for j in jobs if j.srr == "SRR1")
    assert srr1.barcodes == {"AAAAA": 1, "TTTTT": 2}
    srr2 = next(j for j in jobs if j.srr == "SRR2")
    assert srr2.barcodes == {"GGGGG": 1}


def test_demux_output_is_deterministic_across_reruns(tmp_path: Path) -> None:
    """Two independent demux runs over the same input must produce byte-for-byte
    identical .gz files. Regression guard for the deterministic-gzip writer
    in `selexprep._io` (without it, the gzip header's mtime would differ
    across runs and break manifest SHA256 reproducibility)."""
    fq = tmp_path / "SRR1.fastq.gz"
    _write_fastq_gz(fq, ["AAAAACCCCC", "AAAAAGGGGG", "TTTTTGGGGG"])

    out_a = tmp_path / "out_a"
    out_b = tmp_path / "out_b"
    demux_fastq(r1_path=fq, out_dir=out_a, barcodes={"AAAAA": 1, "TTTTT": 2}, srr="SRR1")
    demux_fastq(r1_path=fq, out_dir=out_b, barcodes={"AAAAA": 1, "TTTTT": 2}, srr="SRR1")

    for relative in (
        "round_01/SRR1.fastq.gz",
        "round_02/SRR1.fastq.gz",
    ):
        assert (out_a / relative).read_bytes() == (out_b / relative).read_bytes(), (
            f"non-deterministic demux output at {relative}"
        )


def test_demux_sample_sheet_writes_per_job_reports(tmp_path: Path) -> None:
    fq = tmp_path / "SRR1.fastq.gz"
    _write_fastq_gz(fq, ["AAAAACCCCC", "TTTTTGGGGG"])
    sheet = tmp_path / "sheet.tsv"
    sheet.write_text(
        f"srr\tr1_path\tr2_path\tround\tbarcode\nSRR1\t{fq}\t\t1\tAAAAA\nSRR1\t{fq}\t\t2\tTTTTT\n",
        encoding="utf-8",
    )
    out_root = tmp_path / "out"
    report_dir = tmp_path / "reports"

    reports = demux_sample_sheet(sheet, out_root, report_dir=report_dir)
    assert len(reports) == 1
    assert reports[0].per_round == {1: 1, 2: 1}
    assert (report_dir / "SRR1_demux.json").exists()
