"""Unit tests for ``selexprep.qc.diversity``."""

from __future__ import annotations

import math

import pytest

from selexprep.qc.diversity import rarefy, shannon_entropy, top_n_coverage, unique_count

# ---------------------------------------------------------------------------
# rarefy
# ---------------------------------------------------------------------------


def test_rarefy_total_equals_target_depth() -> None:
    counts = {"A": 100, "B": 50, "C": 25}
    rarefied = rarefy(counts, depth=50, seed=42)
    assert sum(rarefied.values()) == 50


def test_rarefy_is_deterministic_for_same_seed() -> None:
    counts = {f"seq_{i}": i + 1 for i in range(20)}
    a = rarefy(counts, depth=30, seed=42)
    b = rarefy(counts, depth=30, seed=42)
    assert a == b


def test_rarefy_differs_for_different_seeds() -> None:
    counts = {f"seq_{i}": i + 1 for i in range(20)}
    a = rarefy(counts, depth=30, seed=42)
    b = rarefy(counts, depth=30, seed=99)
    # With diverse pool + small subsample, different seeds yield different samples
    assert a != b


def test_rarefy_returns_input_when_depth_exceeds_total() -> None:
    counts = {"A": 5, "B": 3}
    rarefied = rarefy(counts, depth=1000, seed=42)
    assert rarefied == {"A": 5, "B": 3}


def test_rarefy_empty_input_returns_empty() -> None:
    assert rarefy({}, depth=100) == {}


def test_rarefy_rejects_nonpositive_depth() -> None:
    with pytest.raises(ValueError, match="depth must be positive"):
        rarefy({"A": 1}, depth=0)
    with pytest.raises(ValueError, match="depth must be positive"):
        rarefy({"A": 1}, depth=-1)


def test_rarefy_drops_zero_count_sequences() -> None:
    counts = {"A": 1000, "B": 1, "C": 1}
    rarefied = rarefy(counts, depth=10, seed=42)
    # B and C have such low abundance they will probably not be sampled at depth 10
    # — at minimum, only sequences with sampled count > 0 appear in output.
    assert all(v > 0 for v in rarefied.values())


# ---------------------------------------------------------------------------
# shannon_entropy
# ---------------------------------------------------------------------------


def test_shannon_entropy_uniform_is_log2_n() -> None:
    counts = {f"s{i}": 1 for i in range(8)}
    # H = log2(8) = 3
    assert shannon_entropy(counts) == pytest.approx(3.0)


def test_shannon_entropy_single_sequence_is_zero() -> None:
    assert shannon_entropy({"only": 100}) == 0.0


def test_shannon_entropy_empty_is_zero() -> None:
    assert shannon_entropy({}) == 0.0


def test_shannon_entropy_handles_zero_counts() -> None:
    counts = {"A": 1, "B": 0, "C": 1}
    # B's zero contributes nothing; effective H = log2(2) = 1
    assert shannon_entropy(counts) == pytest.approx(1.0)


def test_shannon_entropy_concentrated_pool_lower_than_uniform() -> None:
    uniform = {f"s{i}": 1 for i in range(10)}
    concentrated = {"top": 99, "rest": 1}
    assert shannon_entropy(concentrated) < shannon_entropy(uniform)


# ---------------------------------------------------------------------------
# unique_count
# ---------------------------------------------------------------------------


def test_unique_count_counts_only_nonzero() -> None:
    assert unique_count({"A": 1, "B": 0, "C": 5}) == 2


def test_unique_count_empty_is_zero() -> None:
    assert unique_count({}) == 0


# ---------------------------------------------------------------------------
# top_n_coverage
# ---------------------------------------------------------------------------


def test_top_n_coverage_single_top_sequence() -> None:
    counts = {"top": 90, "rest1": 5, "rest2": 5}
    assert top_n_coverage(counts, n=1) == pytest.approx(0.9)


def test_top_n_coverage_full_when_n_exceeds_uniques() -> None:
    counts = {"A": 5, "B": 3, "C": 2}
    assert top_n_coverage(counts, n=100) == pytest.approx(1.0)


def test_top_n_coverage_empty_is_zero() -> None:
    assert top_n_coverage({}, n=10) == 0.0


def test_top_n_coverage_zero_n_is_zero() -> None:
    counts = {"A": 5, "B": 3}
    assert top_n_coverage(counts, n=0) == 0.0


def test_top_n_coverage_monotonic_in_n() -> None:
    counts = {f"s{i}": (10 - i) for i in range(10)}
    for n in range(1, 11):
        assert top_n_coverage(counts, n) >= top_n_coverage(counts, n - 1)


# ---------------------------------------------------------------------------
# Integration: rarefied entropy <= original entropy (depth-aware)
# ---------------------------------------------------------------------------


def test_rarefied_entropy_does_not_exceed_original_entropy() -> None:
    """Sanity check: subsampling can only equal or reduce diversity."""
    counts = {f"s{i}": (i + 1) for i in range(50)}
    original_h = shannon_entropy(counts)
    rarefied = rarefy(counts, depth=100, seed=42)
    rarefied_h = shannon_entropy(rarefied)
    # Allow small fp tolerance; rarefaction cannot increase entropy meaningfully.
    assert rarefied_h <= original_h + 1e-9
    # And entropy is non-negative.
    assert rarefied_h >= 0.0
    assert not math.isnan(rarefied_h)
