"""Empirical primer / constant-region inference from FASTQ-derived sequences.

Scans the longest consensus prefix/suffix of the top-N most abundant
sequences in the **earliest** available SELEX round (round_00 preferred).

**Why earliest:** naive pools have maximally diverse random regions, so
everything the reads share at the flanks IS the primer by construction. In
late rounds the random region has collapsed to winners, making primer
detection unreliable — you'd "detect" the winning aptamer.

**Algorithm:**
  For each candidate flank length L (from MAX_LEN down to MIN_LEN):
    take first/last L bases of each top-N sequence
    find the most common L-mer
    count sequences whose flank matches within ≤ 1 Hamming distance
    if hits / total ≥ CONFIDENCE: accept; otherwise try shorter L

**Phase-1 caveat:** this module is the algorithm port. Phase 2 wraps it in
``library/report.py`` to emit a full ``LibraryReport`` with adapter
blacklisting, cross-round persistence, extraction_mode classification, and
sampling-seed reproducibility. The bare detection result here is the
algorithmic core that the LibraryReport composer consumes.

**Ordering caveat:** if the earliest round was multiplexed and has not yet
been demuxed, this routine will detect the 5' barcode as the primer.
Callers must run demultiplexing first; ``selexprep`` will enforce this in
the CLI dispatcher.
"""

from __future__ import annotations

import logging
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


# Algorithm defaults.
# DEFAULT_TOP_N = None means "use every unique sequence in the parquet, no
# subsampling". In a naive SELEX pool every read carries the primer at the
# same position, so sampling top-N by read count discards evidence rather
# than concentrating it — the long tail of rare unique sequences also
# carries the primer and confirms the consensus. Pass an explicit int to
# `top_n` only if memory pressure forces it.
DEFAULT_TOP_N: int | None = None
DEFAULT_CONFIDENCE = 0.75
DEFAULT_MIN_LEN = 15
DEFAULT_MAX_LEN = 30
DEFAULT_MIN_SEQS_FOR_DETECTION = 500


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclass
class FlankResult:
    """One side of a primer-pair detection."""

    sequence: str | None
    length: int
    confidence: float  # fraction of sequences matching within Hamming ≤ 1


@dataclass
class PrimerDetection:
    """Algorithmic primer-flank detection result for one input pool.

    This is the Phase-1 raw detection; Phase 2 wraps it in `LibraryReport`
    with adapter blacklisting, cross-round persistence, and `extraction_mode`.
    """

    primer_5p: FlankResult
    primer_3p: FlankResult
    n_sequences_analyzed: int
    mean_seq_len: int
    estimated_random_region_len: int | None
    source: Path | None = None


# ---------------------------------------------------------------------------
# Hamming helper (early-exit when diff > 1)
# ---------------------------------------------------------------------------


def _hamming_le1(a: str, b: str) -> bool:
    """Return True iff Hamming distance between equal-length strings ≤ 1.

    Early-exits on the second mismatch — meaningfully faster for the hot
    path of comparing tens of thousands of L-mers to a candidate consensus.
    """
    if len(a) != len(b):
        return False
    diffs = 0
    for x, y in zip(a, b, strict=False):
        if x != y:
            diffs += 1
            if diffs > 1:
                return False
    return True


# ---------------------------------------------------------------------------
# Flank detection
# ---------------------------------------------------------------------------


def detect_flank(
    sequences: list[str],
    is_prefix: bool,
    max_len: int = DEFAULT_MAX_LEN,
    min_len: int = DEFAULT_MIN_LEN,
    confidence: float = DEFAULT_CONFIDENCE,
    min_seqs_for_detection: int = DEFAULT_MIN_SEQS_FOR_DETECTION,
) -> FlankResult:
    """Detect the longest consensus prefix (or suffix) across `sequences`.

    Returns the longest L (between `min_len` and `max_len`) at which at least
    `confidence` fraction of sequences share an L-mer within Hamming ≤ 1.
    Returns `FlankResult(None, 0, 0.0)` if no consensus reaches `confidence`.
    """
    for L in range(max_len, min_len - 1, -1):
        fragments = [(s[:L] if is_prefix else s[-L:]) for s in sequences if len(s) >= L]
        if len(fragments) < min_seqs_for_detection:
            continue
        most_common, _count = Counter(fragments).most_common(1)[0]
        hits = sum(1 for f in fragments if _hamming_le1(f, most_common))
        frac = hits / len(fragments)
        if frac >= confidence:
            return FlankResult(sequence=most_common, length=L, confidence=frac)
    return FlankResult(sequence=None, length=0, confidence=0.0)


# ---------------------------------------------------------------------------
# Top-level detection from a list of sequences
# ---------------------------------------------------------------------------


def detect_primers(
    sequences: list[str],
    confidence: float = DEFAULT_CONFIDENCE,
    max_len: int = DEFAULT_MAX_LEN,
    min_len: int = DEFAULT_MIN_LEN,
    min_seqs_for_detection: int = DEFAULT_MIN_SEQS_FOR_DETECTION,
    source: Path | None = None,
) -> PrimerDetection | None:
    """Detect 5' and 3' primer flanks from a list of sequences.

    Returns `None` if fewer than `min_seqs_for_detection` sequences are
    supplied. Either flank may still be ``None`` (no consensus at any
    length); the other can still be detected independently.
    """
    if len(sequences) < min_seqs_for_detection:
        return None

    common_kwargs = dict(
        max_len=max_len,
        min_len=min_len,
        confidence=confidence,
        min_seqs_for_detection=min_seqs_for_detection,
    )
    p5 = detect_flank(sequences, is_prefix=True, **common_kwargs)
    p3 = detect_flank(sequences, is_prefix=False, **common_kwargs)

    mean_len = round(pd.Series([len(s) for s in sequences]).mean())
    random_len: int | None = None
    if p5.sequence or p3.sequence:
        random_len = mean_len - p5.length - p3.length

    return PrimerDetection(
        primer_5p=p5,
        primer_3p=p3,
        n_sequences_analyzed=len(sequences),
        mean_seq_len=mean_len,
        estimated_random_region_len=random_len,
        source=source,
    )


# ---------------------------------------------------------------------------
# Wrapper: detect from a counts parquet
# ---------------------------------------------------------------------------


def detect_from_parquet(
    parquet_path: Path,
    top_n: int | None = DEFAULT_TOP_N,
    confidence: float = DEFAULT_CONFIDENCE,
    max_len: int = DEFAULT_MAX_LEN,
    min_len: int = DEFAULT_MIN_LEN,
    min_seqs_for_detection: int = DEFAULT_MIN_SEQS_FOR_DETECTION,
) -> PrimerDetection | None:
    """Detect primers from a counts parquet.

    By default (``top_n=None``) every unique sequence in the parquet is used —
    no subsampling. Pass an explicit positive integer to cap the input at the
    top-N most abundant sequences (useful only when memory is constrained).
    """
    df = pd.read_parquet(parquet_path, columns=["sequence", "reads"])
    if top_n is not None and top_n > 0:
        df = df.nlargest(top_n, "reads")
    sequences = df["sequence"].tolist()
    return detect_primers(
        sequences,
        confidence=confidence,
        max_len=max_len,
        min_len=min_len,
        min_seqs_for_detection=min_seqs_for_detection,
        source=parquet_path,
    )


def earliest_round_parquet(processed_bp_dir: Path) -> Path | None:
    """Return the earliest-numbered ``round_*.counts.parquet`` in `processed_bp_dir`.

    Naive-pool primer detection works best on the earliest round; the random
    region is maximally diverse there, so shared flanks are unambiguously
    primer-derived rather than aptamer-enrichment artefacts.
    """
    parquets = sorted(processed_bp_dir.glob("round_*.counts.parquet"))
    return parquets[0] if parquets else None
