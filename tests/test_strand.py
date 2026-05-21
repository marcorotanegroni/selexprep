"""Unit tests for ``selexprep.extract.strand``."""

from __future__ import annotations

import gzip
from pathlib import Path

from selexprep.extract.strand import (
    _revcomp_read,
    detect_strand_distribution,
    reorient_fastq_gz,
    write_strand_report,
)

PRIMER_5P = "GGTAATACGACTCACTATAGGG"  # T7 promoter, 22 nt
PRIMER_3P = "CCATGCATGCATGCATGCAT"  # 20 nt


# ---------------------------------------------------------------------------
# _revcomp_read
# ---------------------------------------------------------------------------


def test_revcomp_read_basic() -> None:
    assert _revcomp_read("ATCG") == "CGAT"
    assert _revcomp_read("GGGG") == "CCCC"


def test_revcomp_read_tolerates_n() -> None:
    assert _revcomp_read("ANCG") == "CGNT"
    assert _revcomp_read("NNNN") == "NNNN"


def test_revcomp_read_involutive() -> None:
    for seq in ("ATCGATCG", "GCGCATAT", "ACGTN"):
        assert _revcomp_read(_revcomp_read(seq)) == seq


# ---------------------------------------------------------------------------
# detect_strand_distribution
# ---------------------------------------------------------------------------


def _pool(primer_5p: str | None, primer_3p: str | None, n: int = 100) -> list[str]:
    """Tiny synthetic pool: primer_5p + 30-nt random + primer_3p."""
    bases = "ACGT"
    out: list[str] = []
    for i in range(n):
        rand = "".join(bases[(i * 7 + j * 13) % 4] for j in range(30))
        out.append((primer_5p or "") + rand + (primer_3p or ""))
    return out


def test_detect_strand_distribution_all_forward() -> None:
    seqs = _pool(PRIMER_5P, PRIMER_3P, n=100)
    d = detect_strand_distribution(seqs, PRIMER_5P, PRIMER_3P)
    assert d["forward"] == 100
    assert d["reverse"] == 0


def test_detect_strand_distribution_all_reverse() -> None:
    rc_3p = _revcomp_read(PRIMER_3P)
    rc_5p = _revcomp_read(PRIMER_5P)
    seqs = _pool(rc_3p, rc_5p, n=100)
    d = detect_strand_distribution(seqs, PRIMER_5P, PRIMER_3P)
    assert d["reverse"] == 100
    assert d["forward"] == 0


def test_detect_strand_distribution_mixed() -> None:
    forward = _pool(PRIMER_5P, PRIMER_3P, n=80)
    rc_3p = _revcomp_read(PRIMER_3P)
    rc_5p = _revcomp_read(PRIMER_5P)
    reverse = _pool(rc_3p, rc_5p, n=20)
    seqs = forward + reverse
    d = detect_strand_distribution(seqs, PRIMER_5P, PRIMER_3P)
    assert d["forward"] == 80
    assert d["reverse"] == 20
    assert d["ambiguous"] == 0


def test_detect_strand_distribution_no_primers_means_all_ambiguous() -> None:
    seqs = _pool(PRIMER_5P, PRIMER_3P, n=100)
    d = detect_strand_distribution(seqs, None, None)
    assert d["ambiguous"] == 100
    assert d["forward"] == 0
    assert d["reverse"] == 0


def test_detect_strand_distribution_counts_sum_to_input_size() -> None:
    seqs = _pool(PRIMER_5P, PRIMER_3P, n=100) + ["AAAA" * 10] * 5  # 5 ambiguous noise
    d = detect_strand_distribution(seqs, PRIMER_5P, PRIMER_3P)
    assert d["forward"] + d["reverse"] + d["ambiguous"] == 105


# ---------------------------------------------------------------------------
# reorient_fastq_gz
# ---------------------------------------------------------------------------


def _write_fastq_gz(path: Path, seqs: list[str]) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        for i, s in enumerate(seqs):
            fh.write(f"@read_{i}\n{s}\n+\n{'I' * len(s)}\n")


def _read_fastq(path: Path) -> list[tuple[str, str, str]]:
    """Return (name, seq, qual) tuples from a FASTQ.gz."""
    out: list[tuple[str, str, str]] = []
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            s = fh.readline().rstrip("\n")
            _ = fh.readline()
            q = fh.readline().rstrip("\n")
            out.append((h.rstrip("\n"), s, q))
    return out


def test_reorient_fastq_gz_revcomps_each_read(tmp_path: Path) -> None:
    input_fq = tmp_path / "in.fastq.gz"
    output_fq = tmp_path / "out.fastq.gz"
    seqs = ["ATCG", "GGGG", "ANCG"]
    _write_fastq_gz(input_fq, seqs)

    n = reorient_fastq_gz(input_fq, output_fq)
    assert n == 3

    out_records = _read_fastq(output_fq)
    assert [r[1] for r in out_records] == ["CGAT", "CCCC", "CGNT"]


def test_reorient_fastq_gz_reverses_quality_string(tmp_path: Path) -> None:
    input_fq = tmp_path / "in.fastq.gz"
    output_fq = tmp_path / "out.fastq.gz"
    # Custom quality so we can verify it's reversed (not just same)
    with gzip.open(input_fq, "wt", encoding="utf-8") as fh:
        fh.write('@r0\nATCG\n+\n!"#$\n')  # qual = !"#$

    reorient_fastq_gz(input_fq, output_fq)

    out_records = _read_fastq(output_fq)
    # Sequence revcomp'd, quality reversed
    assert out_records[0][1] == "CGAT"
    assert out_records[0][2] == '$#"!'


def test_reorient_fastq_gz_output_is_deterministic(tmp_path: Path) -> None:
    """Two reruns produce byte-identical output (deterministic gzip header)."""
    from selexprep._io import sha256_file

    input_fq = tmp_path / "in.fastq.gz"
    out_a = tmp_path / "a.fastq.gz"
    out_b = tmp_path / "b.fastq.gz"
    _write_fastq_gz(input_fq, ["ATCGATCG", "GGGGAAAA", "TTTTCCCC"])

    reorient_fastq_gz(input_fq, out_a)
    reorient_fastq_gz(input_fq, out_b)

    assert sha256_file(out_a) == sha256_file(out_b)
    assert out_a.read_bytes() == out_b.read_bytes()


# ---------------------------------------------------------------------------
# write_strand_report
# ---------------------------------------------------------------------------


def test_reorient_fastq_gz_raises_on_truncated_record(tmp_path: Path) -> None:
    """Codex Phase 3 pass 1 regression: a truncated FASTQ input must raise
    rather than silently produce a partial output. Previously it logged a
    warning + broke the loop, leaving the caller no way to detect that
    extraction was incomplete."""
    import pytest

    input_fq = tmp_path / "truncated.fastq.gz"
    # Write a valid first record + a truncated second record (header only).
    with gzip.open(input_fq, "wt", encoding="utf-8") as fh:
        fh.write("@r0\nATCG\n+\nIIII\n")  # complete
        fh.write("@r1\n")  # truncated — no seq/plus/qual

    output_fq = tmp_path / "out.fastq.gz"
    with pytest.raises(ValueError, match="truncated FASTQ record"):
        reorient_fastq_gz(input_fq, output_fq)


def test_write_strand_report_emits_sorted_tsv(tmp_path: Path) -> None:
    dists = {
        2: {"forward": 30, "reverse": 5, "ambiguous": 0},
        0: {"forward": 100, "reverse": 0, "ambiguous": 0},
        1: {"forward": 80, "reverse": 10, "ambiguous": 2},
    }
    path = tmp_path / "strand_report.tsv"
    write_strand_report(dists, path)

    lines = path.read_text().splitlines()
    assert lines[0] == "round\tforward\treverse\tambiguous"
    # Rounds in ascending order
    assert lines[1].startswith("0\t100\t0\t0")
    assert lines[2].startswith("1\t80\t10\t2")
    assert lines[3].startswith("2\t30\t5\t0")


def test_write_strand_report_missing_keys_default_to_zero(tmp_path: Path) -> None:
    dists = {0: {"forward": 100}}  # no reverse/ambiguous keys
    path = tmp_path / "strand_report.tsv"
    write_strand_report(dists, path)

    line = path.read_text().splitlines()[1]
    assert line == "0\t100\t0\t0"
