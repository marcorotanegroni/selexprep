"""Diversity helpers for QC flags + plots.

Locked plan line 339 lists ``qc/diversity.py`` (rarefaction) as a Phase 5
module. The rarefied diversity comparison is what makes
``unexpected_rarefied_diversity_increase`` (locked plan line 351) a
depth-aware flag rather than a raw-counts-confounded one.

Public API:

- :func:`rarefy` - sample without replacement to a target read depth
  (deterministic for a given seed).
- :func:`shannon_entropy` - base-2 entropy of the abundance distribution.
- :func:`unique_count` - number of non-zero sequences.
- :func:`top_n_coverage` - fraction of total reads in the top-N
  sequences.
"""

from __future__ import annotations

import math

import numpy as np


def rarefy(counts: dict[str, int], depth: int, *, seed: int = 42) -> dict[str, int]:
    """Subsample ``counts`` without replacement to ``depth`` total reads.

    Uses ``numpy`` multivariate hypergeometric sampling so the per-sequence
    counts after rarefaction sum exactly to ``depth`` (or to the original
    total if ``depth >= total``).

    Deterministic for a given ``seed``: two calls with the same input
    produce the same output. Zero-count sequences are dropped from the
    result.

    Args:
        counts: ``{sequence: reads}`` mapping.
        depth: target sample size (number of reads).
        seed: RNG seed.

    Raises:
        ValueError: when ``depth <= 0``.
    """
    if depth <= 0:
        raise ValueError(f"depth must be positive; got {depth}")
    if not counts:
        return {}

    keys = list(counts.keys())
    values = np.array([counts[k] for k in keys], dtype=np.int64)
    total = int(values.sum())
    if depth >= total:
        # Asking for >= what we have - return a copy of the input.
        return {k: int(v) for k, v in counts.items() if v > 0}

    rng = np.random.default_rng(seed)
    sampled = rng.multivariate_hypergeometric(values, depth)
    out: dict[str, int] = {}
    for k, c in zip(keys, sampled, strict=True):
        c_int = int(c)
        if c_int > 0:
            out[k] = c_int
    return out


def shannon_entropy(counts: dict[str, int]) -> float:
    """Base-2 Shannon entropy of the abundance distribution.

    Returns 0.0 for an empty pool or single-sequence pool. A uniform pool
    of N distinct sequences has entropy ``log2(N)``.
    """
    total = sum(counts.values())
    if total == 0:
        return 0.0
    h = 0.0
    for c in counts.values():
        if c > 0:
            p = c / total
            h -= p * math.log2(p)
    return h


def unique_count(counts: dict[str, int]) -> int:
    """Number of unique sequences with at least one read."""
    return sum(1 for c in counts.values() if c > 0)


def top_n_coverage(counts: dict[str, int], n: int) -> float:
    """Fraction of total reads concentrated in the top ``n`` sequences.

    Returns 0.0 on empty input; clamps to 1.0 when ``n`` exceeds the
    number of unique sequences.
    """
    if n <= 0:
        return 0.0
    total = sum(counts.values())
    if total == 0:
        return 0.0
    top = sorted(counts.values(), reverse=True)[:n]
    return sum(top) / total
