"""Unit tests for ``selexprep.extract.trim`` (cutadapt subprocess wrapper)."""

from __future__ import annotations

import gzip
import shutil
from pathlib import Path

import pytest

from selexprep._io import sha256_file
from selexprep.extract.trim import (
    TrimReport,
    trim_paired_split,
    trim_single_end_3p,
    trim_single_end_5p,
    trim_single_end_linked,
)

PRIMER_5P = "GGTAATACGACTCACTATAGGG"  # T7 promoter, 22 nt
PRIMER_3P = "CCATGCATGCATGCATGCAT"  # 20 nt


pytestmark = pytest.mark.skipif(shutil.which("cutadapt") is None, reason="cutadapt not on PATH")


def _synthetic_fastq(
    path: Path,
    primer_5p: str | None,
    primer_3p: str | None,
    *,
    n: int = 50,
    random_len: int = 30,
) -> None:
    """Write a tiny FASTQ.gz with primer_5p + N-region + primer_3p reads."""
    bases = "ACGT"
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        for i in range(n):
            rand = "".join(bases[(i * 7 + j * 13) % 4] for j in range(random_len))
            seq = (primer_5p or "") + rand + (primer_3p or "")
            fh.write(f"@read_{i}\n{seq}\n+\n{'I' * len(seq)}\n")


def _read_fasta_gz(path: Path) -> list[str]:
    """Return sequence lines from a FASTA.gz (cutadapt output)."""
    out: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                continue
            if line:
                out.append(line)
    return out


# ---------------------------------------------------------------------------
# BOTH_PRIMERS_SINGLE_READ (linked adapter)
# ---------------------------------------------------------------------------


def test_trim_single_end_linked_strips_both_primers(tmp_path: Path) -> None:
    input_fq = tmp_path / "in.fastq.gz"
    output_fa = tmp_path / "extracted.fasta.gz"
    _synthetic_fastq(input_fq, PRIMER_5P, PRIMER_3P, n=20, random_len=30)

    report = trim_single_end_linked(input_fq, output_fa, primer_5p=PRIMER_5P, primer_3p=PRIMER_3P)

    assert isinstance(report, TrimReport)
    assert report.return_code == 0
    assert report.n_in == 20
    assert report.n_out == 20

    seqs = _read_fasta_gz(output_fa)
    assert len(seqs) == 20
    # Every output seq is just the random region (30 nt), neither primer present.
    for s in seqs:
        assert len(s) == 30
        assert PRIMER_5P not in s
        assert PRIMER_3P not in s


def test_trim_linked_records_argv(tmp_path: Path) -> None:
    input_fq = tmp_path / "in.fastq.gz"
    output_fa = tmp_path / "extracted.fasta.gz"
    _synthetic_fastq(input_fq, PRIMER_5P, PRIMER_3P, n=10)

    report = trim_single_end_linked(input_fq, output_fa, primer_5p=PRIMER_5P, primer_3p=PRIMER_3P)

    assert "cutadapt" in report.cutadapt_cmd
    assert "--discard-untrimmed" in report.cutadapt_cmd
    assert "--fasta" in report.cutadapt_cmd
    # Linked-adapter syntax: -g "P5...P3"
    assert f"{PRIMER_5P}...{PRIMER_3P}" in report.cutadapt_cmd


# ---------------------------------------------------------------------------
# FIVE_PRIME_ONLY
# ---------------------------------------------------------------------------


def test_trim_single_end_5p_strips_only_5p(tmp_path: Path) -> None:
    input_fq = tmp_path / "in.fastq.gz"
    output_fa = tmp_path / "partial_5p_extracted.fasta.gz"
    _synthetic_fastq(input_fq, PRIMER_5P, primer_3p=None, n=20, random_len=30)

    report = trim_single_end_5p(input_fq, output_fa, primer_5p=PRIMER_5P)
    assert report.n_out == 20

    seqs = _read_fasta_gz(output_fa)
    # 5' primer removed; the random tail is what survives.
    for s in seqs:
        assert PRIMER_5P not in s
        assert len(s) == 30


# ---------------------------------------------------------------------------
# THREE_PRIME_ONLY
# ---------------------------------------------------------------------------


def test_trim_single_end_3p_strips_only_3p(tmp_path: Path) -> None:
    input_fq = tmp_path / "in.fastq.gz"
    output_fa = tmp_path / "partial_3p_extracted.fasta.gz"
    _synthetic_fastq(input_fq, primer_5p=None, primer_3p=PRIMER_3P, n=20, random_len=30)

    report = trim_single_end_3p(input_fq, output_fa, primer_3p=PRIMER_3P)
    assert report.n_out == 20

    seqs = _read_fasta_gz(output_fa)
    for s in seqs:
        assert PRIMER_3P not in s
        assert len(s) == 30


# ---------------------------------------------------------------------------
# PAIRED_END_SPLIT_PRIMERS
# ---------------------------------------------------------------------------


def test_trim_paired_split_outputs_both_files(tmp_path: Path) -> None:
    # R1 carries 5' primer + 80 nt random tail
    r1 = tmp_path / "r1.fastq.gz"
    _synthetic_fastq(r1, PRIMER_5P, primer_3p=None, n=20, random_len=80)

    # R2 carries revcomp(3' primer) + 80 nt random tail
    from selexprep.library.adapters import reverse_complement

    primer_5p_r2 = reverse_complement(PRIMER_3P)
    r2 = tmp_path / "r2.fastq.gz"
    _synthetic_fastq(r2, primer_5p_r2, primer_3p=None, n=20, random_len=80)

    out_r1 = tmp_path / "partial_5p_extracted_R1.fasta.gz"
    out_r2 = tmp_path / "partial_3p_extracted_R2.fasta.gz"

    report = trim_paired_split(
        r1,
        r2,
        out_r1,
        out_r2,
        primer_5p_r1=PRIMER_5P,
        primer_5p_r2=primer_5p_r2,
    )

    assert report.n_out == 20
    assert len(report.output_paths) == 2
    assert out_r1.exists() and out_r2.exists()

    r1_seqs = _read_fasta_gz(out_r1)
    r2_seqs = _read_fasta_gz(out_r2)
    assert len(r1_seqs) == 20
    assert len(r2_seqs) == 20


# ---------------------------------------------------------------------------
# Deterministic output
# ---------------------------------------------------------------------------


def test_trim_output_is_deterministic_across_reruns(tmp_path: Path) -> None:
    input_fq = tmp_path / "in.fastq.gz"
    out_a = tmp_path / "a.fasta.gz"
    out_b = tmp_path / "b.fasta.gz"
    _synthetic_fastq(input_fq, PRIMER_5P, PRIMER_3P, n=30, random_len=30)

    trim_single_end_linked(input_fq, out_a, primer_5p=PRIMER_5P, primer_3p=PRIMER_3P)
    trim_single_end_linked(input_fq, out_b, primer_5p=PRIMER_5P, primer_3p=PRIMER_3P)

    # The deterministic-gzip re-pack step ensures byte-identical output.
    assert sha256_file(out_a) == sha256_file(out_b)
    assert out_a.read_bytes() == out_b.read_bytes()


def test_trim_cleans_up_intermediate_uncompressed_fa(tmp_path: Path) -> None:
    """The .fa intermediate from cutadapt should be deleted after re-packing."""
    input_fq = tmp_path / "in.fastq.gz"
    output_fa = tmp_path / "extracted.fasta.gz"
    _synthetic_fastq(input_fq, PRIMER_5P, PRIMER_3P, n=10)

    trim_single_end_linked(input_fq, output_fa, primer_5p=PRIMER_5P, primer_3p=PRIMER_3P)

    # No leftover .fa or .cutadapt.json intermediates.
    leftovers = [
        p for p in tmp_path.iterdir() if p.suffix in (".fa", ".json") or ".cutadapt." in p.name
    ]
    assert leftovers == []
