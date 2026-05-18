"""Unit tests for selexprep.count.counter."""

from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd
import pytest

from selexprep.count.counter import (
    count_round,
    pair_sibling,
    pool_diversity_stats,
    reverse_complement,
)

# ----- reverse_complement -----


def test_reverse_complement_basic_dna() -> None:
    assert reverse_complement("ACGT") == "ACGT"
    assert reverse_complement("AAAA") == "TTTT"
    assert reverse_complement("GGCC") == "GGCC"


def test_reverse_complement_handles_rna_uracil() -> None:
    assert reverse_complement("ACGU") == "ACGT"


def test_reverse_complement_preserves_case_in_lookup() -> None:
    assert reverse_complement("acgt") == "acgt"


def test_reverse_complement_handles_n() -> None:
    assert reverse_complement("ANNT") == "ANNT"


# ----- pair_sibling -----


def test_pair_sibling_returns_r2_when_present(tmp_path: Path) -> None:
    r1 = tmp_path / "SRR1_1.fastq.gz"
    r2 = tmp_path / "SRR1_2.fastq.gz"
    r1.touch()
    r2.touch()
    assert pair_sibling(r1) == r2


def test_pair_sibling_returns_none_for_single_end(tmp_path: Path) -> None:
    se = tmp_path / "SRR1.fastq.gz"
    se.touch()
    assert pair_sibling(se) is None


def test_pair_sibling_returns_none_when_r2_missing(tmp_path: Path) -> None:
    r1 = tmp_path / "SRR1_1.fastq.gz"
    r1.touch()
    assert pair_sibling(r1) is None


# ----- pool_diversity_stats -----


def test_pool_diversity_uniform_distribution_max_entropy() -> None:
    df = pd.DataFrame({"sequence": ["A", "C", "G", "T"], "reads": [1, 1, 1, 1]})
    shannon, top1 = pool_diversity_stats(df, n_reads=4)
    assert shannon == pytest.approx(2.0)  # log2(4) = 2 for uniform 4-way
    assert top1 == pytest.approx(0.25)


def test_pool_diversity_collapsed_distribution_low_entropy() -> None:
    df = pd.DataFrame({"sequence": ["X"], "reads": [100]})
    shannon, top1 = pool_diversity_stats(df, n_reads=100)
    assert shannon == pytest.approx(0.0)
    assert top1 == pytest.approx(1.0)


def test_pool_diversity_empty_returns_zeros() -> None:
    df = pd.DataFrame({"sequence": [], "reads": []})
    shannon, top1 = pool_diversity_stats(df, n_reads=0)
    assert shannon == 0.0
    assert top1 == 0.0


# ----- count_round end-to-end smoke -----


def _write_fastq_gz(path: Path, sequences: list[str]) -> None:
    """Write a minimal FASTQ.gz with each sequence as a single read."""
    with gzip.open(path, "wt") as fh:
        for i, seq in enumerate(sequences):
            fh.write(f"@read_{i}\n{seq}\n+\n{'I' * len(seq)}\n")


def test_count_round_minimal_no_primers(tmp_path: Path) -> None:
    fq = tmp_path / "round_00.fastq.gz"
    _write_fastq_gz(fq, ["ACGT", "ACGT", "ACGT", "GGGG", "TTTT"])
    out = tmp_path / "round_00.counts.parquet"

    stats = count_round(fq, out)

    assert out.exists()
    assert stats["n_reads"] == 5
    assert stats["n_unique"] == 3
    assert stats["top_seq_reads"] == 3
    assert stats["top1_frac"] == pytest.approx(0.6)
    assert stats["primer_unknown"] is True
    assert stats["primer_trim_applied"] is False

    df = pd.read_parquet(out)
    assert list(df.columns) == ["sequence", "reads", "rank", "rpm"]
    assert df.loc[0, "sequence"] == "ACGT"
    assert df.loc[0, "reads"] == 3


def test_count_round_empty_fastq(tmp_path: Path) -> None:
    fq = tmp_path / "empty.fastq.gz"
    _write_fastq_gz(fq, [])
    out = tmp_path / "empty.counts.parquet"

    stats = count_round(fq, out)

    assert stats["n_reads"] == 0
    assert stats["n_unique"] == 0
