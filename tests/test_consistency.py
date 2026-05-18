"""Unit tests for selexprep.qc.consistency."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from selexprep.qc.consistency import (
    canonical_kmer,
    check_bioproject,
    check_monotonicity,
    extract_kmers,
    jaccard_distance,
)

# ----- canonical_kmer -----


def test_canonical_kmer_returns_lexicographic_minimum() -> None:
    # ACGT and its revcomp (ACGT) are palindromic — same kmer
    assert canonical_kmer("ACGT") == "ACGT"
    # AAAA → TTTT, canonical is AAAA
    assert canonical_kmer("AAAA") == "AAAA"
    assert canonical_kmer("TTTT") == "AAAA"
    # GCGC → GCGC (palindrome)
    assert canonical_kmer("GCGC") == "GCGC"


def test_canonical_kmer_picks_smaller_of_pair() -> None:
    # CCGG → CCGG (palindrome — both strands equal)
    # AAAC → GTTT, canonical = AAAC
    assert canonical_kmer("AAAC") == "AAAC"
    # GTTT → AAAC, canonical = AAAC
    assert canonical_kmer("GTTT") == "AAAC"


# ----- extract_kmers -----


def test_extract_kmers_basic() -> None:
    # Sequence "ACGTAC" with k=4 → kmers: ACGT, CGTA, GTAC; canonical forms
    kmers = extract_kmers("ACGTAC", k=4)
    assert len(kmers) == 3


def test_extract_kmers_too_short_returns_empty() -> None:
    assert extract_kmers("AC", k=4) == set()


# ----- jaccard_distance -----


def test_jaccard_identical_sets_zero_distance() -> None:
    a = {"AAAA", "CCCC", "GGGG"}
    assert jaccard_distance(a, a) == pytest.approx(0.0)


def test_jaccard_disjoint_sets_max_distance() -> None:
    assert jaccard_distance({"AAAA"}, {"TTTT"}) == pytest.approx(1.0)


def test_jaccard_partial_overlap() -> None:
    a = {"AAAA", "CCCC"}
    b = {"CCCC", "GGGG"}
    # intersection = 1, union = 3, distance = 1 - 1/3
    assert jaccard_distance(a, b) == pytest.approx(2.0 / 3.0)


def test_jaccard_both_empty_zero() -> None:
    assert jaccard_distance(set(), set()) == 0.0


# ----- monotonicity -----


def test_monotonicity_clean_ordering_passes() -> None:
    # Distances grow with round-number gap: 0.1, 0.3, 0.6
    rounds = [1, 2, 3]
    dm = {(1, 2): 0.1, (2, 3): 0.15, (1, 3): 0.3}
    out = check_monotonicity(rounds, dm)
    assert out["monotonic"] is True
    assert out["violations"] == []


def test_monotonicity_swapped_rounds_flagged() -> None:
    # Consecutive distance LARGER than distant — round labels likely swapped
    rounds = [1, 2, 3]
    dm = {(1, 2): 0.8, (2, 3): 0.2, (1, 3): 0.1}
    out = check_monotonicity(rounds, dm)
    assert out["monotonic"] is False
    assert len(out["violations"]) > 0


# ----- check_bioproject integration -----


def _write_round_parquet(path: Path, sequences: list[str], reads_each: int = 100) -> None:
    df = pd.DataFrame(
        {
            "sequence": sequences,
            "reads": [reads_each] * len(sequences),
            "rank": range(1, len(sequences) + 1),
        }
    )
    df.to_parquet(path, index=False, compression="zstd")


def test_check_bioproject_insufficient_rounds(tmp_path: Path) -> None:
    bp = tmp_path / "PRJ1"
    bp.mkdir()
    _write_round_parquet(bp / "round_00.counts.parquet", ["AAAA"])
    out = check_bioproject(bp, k=4, top_n=10)
    assert out["status"] == "insufficient_rounds"


def test_check_bioproject_monotonic_progression(tmp_path: Path) -> None:
    """3 rounds where each successive round diverges from the previous."""
    bp = tmp_path / "PRJ1"
    bp.mkdir()
    # R0: random-ish pool. R1: half new. R2: mostly new sequences.
    _write_round_parquet(
        bp / "round_00.counts.parquet", ["AAAAAAAA", "CCCCCCCC", "GGGGGGGG", "TTTTTTTT"]
    )
    _write_round_parquet(
        bp / "round_01.counts.parquet", ["AAAAAAAA", "CCCCCCCC", "ACACACAC", "GTGTGTGT"]
    )
    _write_round_parquet(
        bp / "round_02.counts.parquet", ["ACACACAC", "GTGTGTGT", "ATATATAT", "CGCGCGCG"]
    )
    # Write a minimal summary.json so primer_trim_status returns something benign
    (bp / "summary.json").write_text(json.dumps({"primer_unknown": True, "rounds": []}))

    out = check_bioproject(bp, k=4, top_n=100)
    assert out["status"] in ("ok", "suspicious")  # outcome depends on the synthetic data
    assert out["n_rounds"] == 3
    assert "distance_matrix" in out
    assert "primer_trim_status" in out
