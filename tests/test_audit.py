"""Unit tests for selexprep.library.audit."""

from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd
import pytest

from selexprep.library.audit import (
    TRUSEQ_R1,
    audit_raw_fastq,
    audit_trimmed_parquet,
    positional_base_freq,
)


# ----- positional_base_freq -----


def test_positional_base_freq_uniform_random() -> None:
    """Truly random region should not flag any position as constant-suspect."""
    # All 4 bases equally represented at each position
    seqs = ["ACGT", "TGCA", "GTAC", "CATG"]
    pbfs = positional_base_freq(seqs, range(4))
    assert len(pbfs) == 4
    for pbf in pbfs:
        assert pbf.top_fraction == pytest.approx(0.25)
        assert pbf.constant_suspect is False


def test_positional_base_freq_constant_position_flagged() -> None:
    """A position where >40% of reads share the same base is constant-suspect."""
    seqs = ["AAAA", "AAGA", "AATA", "AACA"]  # position 0 is always A
    pbfs = positional_base_freq(seqs, range(4))
    assert pbfs[0].top_base == "A"
    assert pbfs[0].top_fraction == pytest.approx(1.0)
    assert pbfs[0].constant_suspect is True


def test_positional_base_freq_skips_short_reads() -> None:
    """Reads shorter than the requested position don't contribute."""
    seqs = ["ACGT", "AC"]
    pbfs = positional_base_freq(seqs, range(4))
    # Position 0 has 2 contributors (both 'A'), position 3 has 1 ('T')
    assert pbfs[0].top_fraction == 1.0
    assert pbfs[3].top_base == "T"
    assert pbfs[3].top_fraction == 1.0


def test_positional_base_freq_empty_position_returns_blank() -> None:
    """Position with no contributors returns top_base='-', top_fraction=0."""
    pbfs = positional_base_freq([], range(2))
    assert pbfs[0].top_base == "-"
    assert pbfs[0].top_fraction == 0.0


# ----- audit_raw_fastq -----


def _write_fastq_gz(path: Path, sequences: list[str]) -> None:
    with gzip.open(path, "wt") as fh:
        for i, seq in enumerate(sequences):
            fh.write(f"@read_{i}\n{seq}\n+\n{'I' * len(seq)}\n")


def test_audit_raw_fastq_samples_reads(tmp_path: Path) -> None:
    fq = tmp_path / "sample.fastq.gz"
    _write_fastq_gz(fq, ["AAACCCGGGT", "TTTAAAGGGC", "GGGCCCTTAA"])
    audit = audit_raw_fastq(fq, sample_n=10, window=10)
    assert audit.n_sampled == 3
    assert audit.mean_length == 10.0
    assert audit.first_5_reads == ["AAACCCGGGT", "TTTAAAGGGC", "GGGCCCTTAA"]
    assert len(audit.first_30nt_freq) == 10
    assert len(audit.last_30nt_freq) == 10


def test_audit_raw_fastq_3prime_alignment_handles_variable_length(tmp_path: Path) -> None:
    fq = tmp_path / "varlen.fastq.gz"
    # All reads end with 'XYZ' regardless of length — should be detected by 3'-aligned audit
    _write_fastq_gz(fq, ["AAACG" + "GAC", "TTTGGGAAAC" + "GAC", "CCC" + "GAC"])
    audit = audit_raw_fastq(fq, sample_n=10, window=3)
    # Last position (3'-most) is always 'C'; second-to-last 'A'; third-to-last 'G'
    assert audit.last_30nt_freq[0].top_base == "C"
    assert audit.last_30nt_freq[0].top_fraction == pytest.approx(1.0)
    assert audit.last_30nt_freq[1].top_base == "A"
    assert audit.last_30nt_freq[2].top_base == "G"


def test_audit_raw_fastq_empty(tmp_path: Path) -> None:
    fq = tmp_path / "empty.fastq.gz"
    _write_fastq_gz(fq, [])
    audit = audit_raw_fastq(fq)
    assert audit.n_sampled == 0
    assert audit.mean_length == 0.0


# ----- audit_trimmed_parquet -----


def _write_counts_parquet(path: Path, rows: list[tuple[str, int]]) -> None:
    df = pd.DataFrame(rows, columns=["sequence", "reads"])
    df.to_parquet(path, index=False, compression="zstd")


def test_audit_trimmed_parquet_basic(tmp_path: Path) -> None:
    p = tmp_path / "round_03.counts.parquet"
    _write_counts_parquet(
        p,
        [
            ("AAAACCCCGG", 100),  # L=10, top
            ("TTTAGGCCAA", 80),   # L=10
            ("ATCG", 5),          # L=4 (outlier)
        ],
    )
    audit = audit_trimmed_parquet(p)
    assert audit.n_unique == 3
    assert audit.n_reads == 185
    assert audit.mode_length == 10
    assert audit.top_5_sequences[0] == ("AAAACCCCGG", 100)
    assert audit.truseq_r1_unique_count == 0
    assert audit.truseq_r1_read_count == 0


def test_audit_trimmed_parquet_detects_truseq_contamination(tmp_path: Path) -> None:
    p = tmp_path / "round_04.counts.parquet"
    _write_counts_parquet(
        p,
        [
            ("ACGT" + TRUSEQ_R1 + "AAAA", 50),
            ("TTTAGGCCAA", 30),
        ],
    )
    audit = audit_trimmed_parquet(p)
    assert audit.truseq_r1_unique_count == 1
    assert audit.truseq_r1_read_count == 50


def test_audit_trimmed_parquet_empty(tmp_path: Path) -> None:
    p = tmp_path / "empty.counts.parquet"
    _write_counts_parquet(p, [])
    audit = audit_trimmed_parquet(p)
    assert audit.n_unique == 0
    assert audit.n_reads == 0
    assert audit.mode_length == 0
