"""Empirical primer / constant-region inference from FASTQ-derived sequences.

Scans the consensus prefix/suffix of the top-N most abundant sequences in
the **earliest** available SELEX round (round_00 preferred).

**Why earliest:** naive pools have maximally diverse random regions, so
everything the reads share at the flanks IS the primer by construction. In
late rounds the random region has collapsed to winners, making primer
detection unreliable — you'd "detect" the winning aptamer.

**Algorithm:**
  Walk inward from each read boundary and compute the positional consensus
  base + support. The flank boundary is the support cliff: high-support
  constant-region positions followed by a sustained drop into the random
  region. This avoids both old failure modes: hard-capping long T7-bearing
  constants at 30 nt, and over-extending by one random base when a longer
  candidate still passes Hamming ≤ 1 globally.

**Phase-1 caveat:** this module is the algorithm port. wraps it in
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
import random
import statistics
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from selexprep.library.adapters import (
    ADAPTER_PROBE_K,
    count_adapter_hits,
    matches_known_adapter_prefix,
    reverse_complement,
)
from selexprep.library.report import (
    LibraryReport,
    Orientation,
    ReadSource,
    _classify,
)

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
# Real public deposits include 14 nt library constants.
# Keep this conservative: shorter flanks need explicit future calibration.
DEFAULT_MIN_LEN = 14
DEFAULT_MAX_LEN = 60
DEFAULT_MIN_SEQS_FOR_DETECTION = 500
BOUNDARY_SUPPORT_FLOOR = 0.55
BOUNDARY_STRONG_POSITION_FRACTION = 0.80
BOUNDARY_DROP_RUN = 2
BOUNDARY_HIGH_SUPPORT_BASELINE = 0.90
BOUNDARY_HIGH_SUPPORT_DROP = 0.12
BOUNDARY_HIGH_SUPPORT_POST_MAX = 0.85

# Outer-edge core rescue.
#
# ``_positional_consensus`` walks inward from the *read edge*, so any constant
# technical sequence sitting between the read edge and the library constant —
# a truncated sequencing adapter, an index remnant, the T7 start-G — is absorbed
# into the called flank. That material is often poorly conserved (it sits at the
# low-quality end of the read), which makes the whole-primer ``match_rate``
# collapse and downgrades an otherwise two-sided library to one-sided extraction.
#
# The rescue recomputes ``match_rate`` against the high-support *core* of the
# flank (outer positions with weak support dropped) and is applied ONLY when the
# full-length rate already failed, so libraries that pass today take the exact
# same code path they take now. The core must still be a credible primer:
# ``CORE_MIN_SUPPORT`` keeps it strongly conserved and ``CORE_MIN_LEN`` stops a
# short, spurious consensus from being rescued into a false primer call (the
# specificity arm must keep making no call at all).
CORE_MIN_SUPPORT = 0.90
CORE_MIN_LEN = 12
# Consecutive well-supported positions that end the outer-edge scan. Measured
# on PRJEB62495: the outer run is 0.75/0.90/0.90/0.89/0.48/0.61, then the flank
# settles at 0.90/0.91/0.98 — while an isolated 0.77 sits 21 nt deep, well
# inside the genuine constant. Stopping at the first stable run keeps the cut
# on the outer artefact instead of chasing that interior dip.
CORE_STABLE_RUN = 3


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclass
class FlankResult:
    """One side of a primer-pair detection."""

    sequence: str | None
    length: int
    confidence: float  # fraction of sequences matching within Hamming ≤ 1
    # Per-position consensus support for the called flank, ordered 5'→3' like
    # ``sequence``. Used by the outer-edge core rescue; empty when no flank was
    # called. Internal detail — not part of the ``LibraryReport`` schema.
    supports: list[float] = field(default_factory=list)


@dataclass
class PrimerDetection:
    """Algorithmic primer-flank detection result for one input pool.

    This is the Phase-1 raw detection; wraps it in `LibraryReport`
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


def _positional_consensus(
    sequences: list[str],
    *,
    is_prefix: bool,
    max_len: int,
    min_seqs_for_detection: int,
) -> tuple[list[str], list[float]]:
    """Return consensus bases and support, walking inward from one read edge.

    Prefix consensus is returned left-to-right from read start. Suffix
    consensus is returned from read end inward; callers reverse the selected
    segment to recover the natural 5'→3' suffix sequence.
    """
    consensus: list[str] = []
    supports: list[float] = []
    for offset in range(max_len):
        bases: list[str] = []
        for seq in sequences:
            if len(seq) <= offset:
                continue
            bases.append(seq[offset] if is_prefix else seq[-(offset + 1)])
        if len(bases) < min_seqs_for_detection:
            break
        base, count = Counter(bases).most_common(1)[0]
        consensus.append(base)
        supports.append(count / len(bases))
    return consensus, supports


def _boundary_length(
    supports: list[float],
    *,
    min_len: int,
    confidence: float,
) -> int:
    """Choose a flank length from positional support values.

    ``confidence`` remains the strong-position threshold. The lower
    ``BOUNDARY_SUPPORT_FLOOR`` allows true constant positions with modest
    quality decay, while ``BOUNDARY_DROP_RUN`` prevents a single noisy base
    from truncating the flank. A candidate must still contain mostly strong
    positions, so low-complexity or biased random regions do not become
    primers merely by sitting above the floor.
    """
    if len(supports) < min_len:
        return 0

    boundary = len(supports)
    for i in range(min_len, len(supports) - BOUNDARY_DROP_RUN + 1):
        prev = supports[i - min_len : i]
        baseline = statistics.median(prev)
        window = supports[i : i + BOUNDARY_DROP_RUN]
        if (
            baseline >= BOUNDARY_HIGH_SUPPORT_BASELINE
            and supports[i] <= baseline - BOUNDARY_HIGH_SUPPORT_DROP
            and all(v < BOUNDARY_HIGH_SUPPORT_POST_MAX for v in window)
        ):
            boundary = i
            break

    floor = min(confidence, BOUNDARY_SUPPORT_FLOOR)
    for i in range(len(supports) - BOUNDARY_DROP_RUN + 1):
        if i >= boundary:
            break
        window = supports[i : i + BOUNDARY_DROP_RUN]
        if all(v < floor for v in window):
            boundary = i
            break

    if boundary < min_len:
        return 0

    called = supports[:boundary]
    strong = sum(1 for v in called if v >= confidence)
    if strong / len(called) < BOUNDARY_STRONG_POSITION_FRACTION:
        return 0
    return boundary


def detect_flank(
    sequences: list[str],
    is_prefix: bool,
    max_len: int = DEFAULT_MAX_LEN,
    min_len: int = DEFAULT_MIN_LEN,
    confidence: float = DEFAULT_CONFIDENCE,
    min_seqs_for_detection: int = DEFAULT_MIN_SEQS_FOR_DETECTION,
) -> FlankResult:
    """Detect the longest consensus prefix (or suffix) across `sequences`.

    Returns the boundary-cliff consensus length between ``min_len`` and
    ``max_len``. The reported confidence preserves the original FlankResult
    semantics: fraction of sequences whose called flank matches the consensus
    within Hamming ≤ 1.
    """
    consensus, supports = _positional_consensus(
        sequences,
        is_prefix=is_prefix,
        max_len=max_len,
        min_seqs_for_detection=min_seqs_for_detection,
    )
    L = _boundary_length(supports, min_len=min_len, confidence=confidence)
    if L:
        selected = consensus[:L] if is_prefix else list(reversed(consensus[:L]))
        sequence = "".join(selected)
        # ``supports`` is indexed from the read edge inward; reverse it for a
        # suffix so it lines up 5'→3' with ``sequence``.
        selected_supports = supports[:L] if is_prefix else list(reversed(supports[:L]))
        fragments = [(s[:L] if is_prefix else s[-L:]) for s in sequences if len(s) >= L]
        if len(fragments) >= min_seqs_for_detection:
            hits = sum(1 for f in fragments if _hamming_le1(f, sequence))
            return FlankResult(
                sequence=sequence,
                length=L,
                confidence=hits / len(fragments),
                supports=selected_supports,
            )
    return FlankResult(sequence=None, length=0, confidence=0.0)


def _high_support_core(
    sequence: str | None,
    supports: list[float],
    *,
    is_prefix: bool,
    min_support: float = CORE_MIN_SUPPORT,
    min_len: int = CORE_MIN_LEN,
) -> str | None:
    """Drop the noisy run at a flank's *outer* edge.

    The outer edge is the read edge: the start of a 5' flank, the end of a 3'
    flank. Positions are dropped up to and including the deepest weak one in
    that outer run, so an adapter/index remnant goes together with anything
    outside it.

    The scan stops after ``CORE_STABLE_RUN`` consecutive well-supported
    positions. Real flanks are not uniformly perfect — an isolated dip deep
    inside a genuine constant region is normal — and without the stop a single
    such dip would propose cutting almost the whole flank away. The inner edge,
    the boundary with the random region, is therefore never reached.

    Returns ``None`` when there is nothing to rescue: no flank, a clean outer
    edge, or a core shorter than ``min_len``.
    """
    if not sequence or len(supports) != len(sequence):
        return None

    # Order positions outer-edge-first: index 0 for a prefix, last index for a
    # suffix. Walk inward, remembering the deepest weak position, and stop once
    # the flank has settled into a stable run.
    ordered = supports if is_prefix else list(reversed(supports))
    cut = 0
    stable = 0
    for i, value in enumerate(ordered):
        if value < min_support:
            cut = i + 1
            stable = 0
            continue
        stable += 1
        if stable >= CORE_STABLE_RUN:
            break
    if cut == 0:
        return None

    core = sequence[cut:] if is_prefix else sequence[: len(sequence) - cut]
    if len(core) < min_len:
        return None
    return core


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


# ===========================================================================
# cross-round LibraryReport orchestration
# ===========================================================================
#
# Everything below this banner is work (the LibraryReport pipeline);
# everything above is the single-pool flank detector. The split is
# intentional: the design keeps the functions as the bare
# algorithmic core that wraps with calibration, persistence, and
# classification. callers (tests, future ad-hoc tools) keep working.
#
# calibration constants were reviewed:
# 6 confirmed, 4 revised (POSITION_CONSISTENCY_TOLERANCE,
# STATUS_HIGH_CUTOFF, COMPOSITE_WEIGHTS, COMPOSITE_WEIGHTS_NO_ROUND_MAP).
# Behavior-based tests mean any future tuning will not break the test
# suite. Search `CALIBRATION-REVIEWED` for the reviewed values;
# `CALIBRATION-TODO` for what still awaits review (qc flags,
# adapter blacklist composition).

# ---------------------------------------------------------------------------
# Calibration constants (reviewed 2026-05-20)
# ---------------------------------------------------------------------------
# CONFIRMED values keep their original default; REVISED values cite
# the new evidence in their comment. Benchmark recovery numbers
# will provide empirical ground truth for a future tuning pass.

# CALIBRATION-REVIEWED: CONFIRMED at 0.70.
# the design. AptaPLEX tracks primer errors per read but does
# not publish a dataset-level threshold (AptaSUITE import docs), so 0.70
# is a reasonable pre-benchmark default.
PRIMER_FOUND_MATCH_RATE_THRESHOLD = 0.7

# CALIBRATION-REVIEWED: CONFIRMED at 0.80.
# the design, 305. AptaPLEX supports randomized-region length
# bounds rather than exact modal fractions; 0.80 is an internal safety
# proxy for "sharply peaked" N-region length.
N_LENGTH_CONFIDENT_FRACTION = 0.8

# CALIBRATION-REVIEWED: CONFIRMED at 0.40.
# the design. Existing tools assume supplied primers and
# discard primer-failure reads rather than inferring from weak evidence,
# so refusing-to-extract at <40% on both sides is the safe default
# (AptaTools / AptaPLEX).
UNABLE_TO_EXTRACT_MATCH_RATE = 0.4

# CALIBRATION-REVIEWED: REVISED 2 → 3.
# AptaPLEX's default primer mismatch tolerance is 3; small public-data
# offset noise should not over-penalize an otherwise stable flank.
POSITION_CONSISTENCY_TOLERANCE = 3

# CALIBRATION-REVIEWED:
#   STATUS_HIGH_CUTOFF: REVISED 0.80 → 0.85 — "HIGH" should mean
#     paper-grade high-confidence, harder to reach via additive
#     secondary signals before benchmark calibration.
#   STATUS_MEDIUM_CUTOFF: CONFIRMED at 0.60 — "usable with caution"
#     boundary, without claiming benchmark-grade primer recovery.
#   STATUS_LOW_CUTOFF: CONFIRMED at 0.30 — below this → UNABLE_TO_INFER;
#     the separate <0.40 match-rate rule already blocks unsafe extraction.
STATUS_HIGH_CUTOFF = 0.85
STATUS_MEDIUM_CUTOFF = 0.60
STATUS_LOW_CUTOFF = 0.30

# CALIBRATION-REVIEWED: REVISED weights.
# Rationale: position_consistency deserves parity with raw match rate;
# cross-round persistence is the unique SELEX-specific signal and
# deserves the largest weight when available; adapter_clean is already
# enforced upstream as a blacklist so it should be a small confidence
# bonus, not a driver (Hoinka et al. 2015; AptaTRACE / AptaTools).
COMPOSITE_WEIGHTS = {
    "match_5p": 0.15,
    "match_3p": 0.15,
    "pos_5p": 0.15,
    "pos_3p": 0.15,
    "persistence": 0.25,
    "n_len": 0.10,
    "adapter_clean": 0.05,
}
# Rationale: with persistence absent and status already capped at MEDIUM
# , within-round evidence (match + position) gets
# equal parity weight; n_len and adapter_clean remain supporting signals.
COMPOSITE_WEIGHTS_NO_ROUND_MAP = {
    "match_5p": 0.225,
    "match_3p": 0.225,
    "pos_5p": 0.225,
    "pos_3p": 0.225,
    "persistence": 0.0,
    "n_len": 0.05,
    "adapter_clean": 0.05,
}

# CALIBRATION-REVIEWED: CONFIRMED.
# < 5% reads reversed → FORWARD; 5-95% → MIXED; > 95% → REVERSE.
# No published benchmark; conservative defaults that avoid overreacting
# to contamination/index bleed and require near-unanimous reverse
# evidence before auto-flipping every read.
ORIENTATION_REVERSED_FORWARD_MAX = 0.05
ORIENTATION_REVERSED_REVERSE_MIN = 0.95

# Top-K variants surfaced in the LibraryReport.
VARIANTS_TOP_K = 3
# ADAPTER_PROBE_K is imported from selexprep.library.adapters (single source
# of truth — same value used by count_adapter_hits and
# matches_known_adapter_prefix).


# ---------------------------------------------------------------------------
# Per-signal helpers
# ---------------------------------------------------------------------------


def _normalize_u_to_t(seq: str) -> str:
    """Convert RNA primer notation to DNA."""
    return seq.upper().replace("U", "T")


def _normalize_pool(seqs: list[str]) -> list[str]:
    """Uppercase + U→T per sequence; leaves non-ACGT residues alone here
    (audit/blacklist handle anomalies separately)."""
    return [_normalize_u_to_t(s) for s in seqs]


def _top_k_variants(
    seqs: list[str], length: int, *, is_prefix: bool, k: int = VARIANTS_TOP_K
) -> list[tuple[str, int]]:
    """Return the top-K most common ``length``-mers at the read's flank.

    Used to populate ``variants_5p`` / ``variants_3p`` for downstream
    review when the primary primer is ambiguous.
    """
    if length <= 0:
        return []
    fragments = [(s[:length] if is_prefix else s[-length:]) for s in seqs if len(s) >= length]
    return Counter(fragments).most_common(k)


def _position_consistency(
    seqs: list[str],
    primer: str | None,
    *,
    is_prefix: bool,
    tolerance: int = POSITION_CONSISTENCY_TOLERANCE,
) -> float:
    """Fraction of reads where ``primer`` appears within ±tolerance of the expected flank position.

    Hamming ≤ 1 is allowed (matches ``detect_flank`` semantics). Returns
    0.0 when ``primer`` is None or no reads are long enough.
    """
    if not primer:
        return 0.0
    L = len(primer)
    hits, total = 0, 0
    for s in seqs:
        if len(s) < L:
            continue
        total += 1
        if is_prefix:
            for offset in range(tolerance + 1):
                start = offset
                end = offset + L
                if end > len(s):
                    continue
                if _hamming_le1(s[start:end], primer):
                    hits += 1
                    break
        else:
            for offset in range(tolerance + 1):
                if offset == 0:
                    chunk = s[-L:]
                else:
                    if L + offset > len(s):
                        continue
                    chunk = s[-(L + offset) : -offset]
                if _hamming_le1(chunk, primer):
                    hits += 1
                    break
    return hits / total if total else 0.0


def _substring_match_rate(seqs: list[str], primer: str | None) -> float:
    """Fraction of reads where ``primer`` appears anywhere as a substring,
    with Hamming distance ≤ 1.

    Distinct from :func:`_position_consistency`, which requires the primer
    at the expected flank position ± tolerance. The two signals contribute
    independently to the composite confidence:

    - ``match_rate_*`` (this function) = primer is detectable in the read
      anywhere (loose; informative about presence).
    - ``position_consistency_*`` (:func:`_position_consistency`) = primer
      sits at the expected flank position (strict; informative about
      structural integrity).

    Returns 0.0 when ``primer`` is None or no read is long enough.
    """
    if not primer:
        return 0.0
    L = len(primer)
    hits, total = 0, 0
    for s in seqs:
        if len(s) < L:
            continue
        total += 1
        # Slide a length-L window across the read; first Hamming-≤1 hit wins.
        for i in range(len(s) - L + 1):
            if _hamming_le1(s[i : i + L], primer):
                hits += 1
                break
    return hits / total if total else 0.0


def _contains_hamming_le1(longer: str, shorter: str) -> bool:
    """Return True when ``shorter`` appears in ``longer`` with ≤1 mismatch."""
    if len(shorter) > len(longer):
        return False
    for i in range(len(longer) - len(shorter) + 1):
        if _hamming_le1(longer[i : i + len(shorter)], shorter):
            return True
    return False


def _primer_sequences_agree(a: str | None, b: str | None) -> bool:
    """Loose primer-equivalence check for paired-end split arbitration.

    Equal-length candidates may differ by one sequencing error. Different
    lengths are treated as agreeing only when the shorter candidate is
    contained in the longer one, again allowing one mismatch. This lets a
    single-read 3' call with a fuzzy ±1 boundary agree with the R2-derived
    insert 3' primer, while rejecting unrelated technical suffixes.
    """
    if not a or not b:
        return False
    if len(a) == len(b):
        return _hamming_le1(a, b)
    longer, shorter = (a, b) if len(a) > len(b) else (b, a)
    return _contains_hamming_le1(longer, shorter)


def _persistence_score(match_rates_per_round: list[float]) -> float | None:
    """Cross-round persistence as ``1 - clip(stdev/mean, 0, 1)``.

    Returns ``None`` if fewer than 2 rounds are available, OR if the
    mean rate is below 10% (in both cases persistence is NOT computable
    in a meaningful sense). The composite-confidence formula treats
    ``None`` as "not evaluable" and redistributes weight implicitly via
    the `if v is None: continue` short-circuit. Returning ``None`` here
    instead of ``0.0`` preserves the semantic distinction between "no
    signal to evaluate" and "evaluated and bad".
    """
    if len(match_rates_per_round) < 2:
        return None
    mean = statistics.mean(match_rates_per_round)
    if mean < 0.1:
        return None
    cv = statistics.stdev(match_rates_per_round) / mean
    return max(0.0, min(1.0, 1.0 - cv))


def _combine_persistence(
    p5: float | None,
    p3: float | None,
    primer_5p: str | None,
    primer_3p: str | None,
) -> float | None:
    """Combine 5'/3' persistence scores into one composite-input value.

    - Both primers + both scores: arithmetic mean.
    - One primer side missing: use the other.
    - Both missing: ``None``.
    """
    if primer_5p is None and primer_3p is None:
        return None
    if primer_5p is None:
        return p3
    if primer_3p is None:
        return p5
    if p5 is None and p3 is None:
        return None
    if p5 is None:
        return p3
    if p3 is None:
        return p5
    return (p5 + p3) / 2.0


def _n_length_stats(
    seqs: list[str], primer_5p_len: int, primer_3p_len: int
) -> tuple[int | None, dict[int, int], float]:
    """N-region length mode + distribution + peakedness confidence.

    ``n_length_confidence = mode_count / total`` — fraction of reads that
    fall in the modal-length bucket. A truly clean library is sharply
    peaked (≥ 0.8); smeared length distributions imply trimming failure or
    sequencing length variability.
    """
    counts: Counter[int] = Counter()
    for s in seqs:
        n = max(0, len(s) - primer_5p_len - primer_3p_len)
        counts[n] += 1
    if not counts:
        return None, {}, 0.0
    mode, mode_count = counts.most_common(1)[0]
    total = sum(counts.values())
    return mode, dict(counts), mode_count / total


def _detect_orientation(
    seqs: list[str], primer_5p: str | None, primer_3p: str | None
) -> Orientation:
    """Strand orientation summary from observed flank patterns.

    For each read, classify as "forward" (starts with primer_5p) or
    "reverse" (starts with revcomp(primer_3p)). The fraction of reverse
    reads vs total drives the FORWARD / MIXED / REVERSE call.
    """
    if not primer_5p and not primer_3p:
        return "FORWARD"

    rc_3p = reverse_complement(primer_3p) if primer_3p else None

    forward = 0
    reverse = 0
    for s in seqs:
        if primer_5p and len(s) >= len(primer_5p) and _hamming_le1(s[: len(primer_5p)], primer_5p):
            forward += 1
            continue
        if rc_3p and len(s) >= len(rc_3p) and _hamming_le1(s[: len(rc_3p)], rc_3p):
            reverse += 1
    total = forward + reverse
    if total == 0:
        return "FORWARD"
    reversed_fraction = reverse / total
    if reversed_fraction < ORIENTATION_REVERSED_FORWARD_MAX:
        return "FORWARD"
    if reversed_fraction > ORIENTATION_REVERSED_REVERSE_MIN:
        return "REVERSE"
    return "MIXED"


def _composite_confidence(signals: dict[str, float | None], *, has_round_map: bool) -> float:
    """Weighted sum of per-signal scores → composite confidence in [0, 1].

    ``None`` signals (e.g. persistence with single round) contribute 0.
    """
    weights = COMPOSITE_WEIGHTS if has_round_map else COMPOSITE_WEIGHTS_NO_ROUND_MAP
    total = 0.0
    for key, w in weights.items():
        v = signals.get(key)
        if v is None:
            continue
        total += w * v
    return max(0.0, min(1.0, total))


def _detect_paired_split_signals(r1_seqs: list[str], r2_seqs: list[str]) -> tuple[bool, str | None]:
    """Test for the paired-end split-primer pattern.

    Returns ``(has_paired_split, primer_3p_from_r2)``. ``has_paired_split``
    is True when R1 carries a strong 5' primer and R2 carries a strong 5'
    primer (which when reverse-complemented gives the real 3' primer of
    the insert), provided R1 does not already carry an agreeing 3' primer.
    A strong but conflicting R1 suffix is treated as technical readthrough
    rather than a reason to ignore R2.

    No overlap detection in v0.1 — read merging is v0.2; until then the
    caller treats this as ``required_action=READ_MERGING_RECOMMENDED``.
    """
    r1_det = detect_primers(r1_seqs)
    r2_det = detect_primers(r2_seqs)
    if r1_det is None or r2_det is None:
        return False, None
    r1_5p_strong = r1_det.primer_5p.confidence > PRIMER_FOUND_MATCH_RATE_THRESHOLD
    r1_3p_strong = r1_det.primer_3p.confidence > PRIMER_FOUND_MATCH_RATE_THRESHOLD
    r2_5p_strong = r2_det.primer_5p.confidence > PRIMER_FOUND_MATCH_RATE_THRESHOLD
    if not (r1_5p_strong and r2_5p_strong and r2_det.primer_5p.sequence is not None):
        return False, None

    primer_3p_from_r2 = reverse_complement(r2_det.primer_5p.sequence)
    r1_3p_agrees_with_r2 = _primer_sequences_agree(
        r1_det.primer_3p.sequence,
        primer_3p_from_r2,
    )
    if not r1_3p_strong or not r1_3p_agrees_with_r2:
        return True, primer_3p_from_r2
    return False, None


def _build_unable_report(
    *,
    read_source: ReadSource,
    sampling_seed: int,
    failure_reason: str,
    known_adapter_hits: dict[str, int] | None = None,
) -> LibraryReport:
    """Construct the ``UNABLE_TO_INFER`` short-circuit report.

    Used when input is empty, below detection floor, or both primer match
    rates fall below ``UNABLE_TO_EXTRACT_MATCH_RATE``.
    """
    return LibraryReport(
        primer_5p=None,
        primer_3p=None,
        variants_5p=[],
        variants_3p=[],
        known_adapter_hits=known_adapter_hits or {},
        extraction_mode="UNABLE_TO_EXTRACT",
        full_insert_recovered=False,
        read_source=read_source,
        required_action="MANUAL_PRIMERS_REQUIRED",
        orientation="FORWARD",
        n_length_mode=None,
        n_length_distribution={},
        n_length_confidence=0.0,
        match_rate_5p=0.0,
        match_rate_3p=0.0,
        position_consistency_5p=0.0,
        position_consistency_3p=0.0,
        read_fraction_used_for_inference=0.0,
        sampling_seed=sampling_seed,
        confidence=0.0,
        status="UNABLE_TO_INFER",
        failure_reason=failure_reason,
    )


def _subsample(
    seqs: list[str], max_reads: int | None, rng: random.Random
) -> tuple[list[str], float]:
    """Optionally subsample ``seqs`` to ``max_reads`` using ``rng``.

    Returns ``(subsampled_list, fraction_used)``. When ``max_reads`` is
    None or ≥ ``len(seqs)``, returns the input as-is with fraction=1.0.
    """
    if max_reads is None or max_reads <= 0 or len(seqs) <= max_reads:
        return seqs, 1.0
    sampled = rng.sample(seqs, max_reads)
    return sampled, max_reads / len(seqs)


# ---------------------------------------------------------------------------
# Top-level orchestrator
# ---------------------------------------------------------------------------


def compute_library_report(
    sequences_by_round: dict[int, list[str]],
    *,
    read_source: ReadSource,
    paired_mate_streams: dict[int, list[str]] | None = None,
    sampling_seed: int = 42,
    max_reads_per_round: int | None = None,
) -> LibraryReport:
    """Compute a full ``LibraryReport`` from per-round sequence pools.

    Algorithm:

    1. Adapter-blacklist scan on the earliest round (records hits; does
       NOT filter reads).
    2. Single-pool primer detection on the earliest round
       (``detect_primers``).
    3. Cross-round persistence — per-round match rates of the detected
       primer; persistence = ``1 - clip(stdev/mean, 0, 1)``.
    4. Position consistency at the expected flank with ±tolerance.
    5. Reverse-complement orientation check.
    6. N-region length distribution + peakedness confidence.
    7. Paired-end split detection if ``paired_mate_streams`` is provided.
    8. Composite confidence (weighted sum); ``_classify`` (in
       ``library/report.py``) maps signals → extraction_mode + workflow
       guidance + status.

    Args:
        sequences_by_round: mapping ``{round_number: [seq, ...]}``. Order
            does not matter; the lowest round number is used as the
            naive-pool reference for primer detection.
        read_source: which physical read(s) carry the random region.
        paired_mate_streams: optional R2 stream per round (only when
            ``read_source == "R1_AND_R2"``).
        sampling_seed: seeds the subsampling RNG so two runs with the
            same input produce the same report.
        max_reads_per_round: subsample cap (None = use all reads).
    """
    rng = random.Random(sampling_seed)

    if not sequences_by_round:
        return _build_unable_report(
            read_source=read_source,
            sampling_seed=sampling_seed,
            failure_reason="No sequences provided",
        )

    # Normalize + optionally subsample each round
    normalized: dict[int, list[str]] = {}
    fractions: list[float] = []
    for r, seqs in sequences_by_round.items():
        norm = _normalize_pool(seqs)
        sampled, frac = _subsample(norm, max_reads_per_round, rng)
        normalized[r] = sampled
        fractions.append(frac)
    read_fraction = sum(fractions) / len(fractions) if fractions else 0.0

    rounds_sorted = sorted(normalized.keys())
    earliest_round = rounds_sorted[0]
    earliest_seqs = normalized[earliest_round]

    # Adapter blacklist hits on the earliest round (most informative — naive
    # pool has the highest residual-adapter representation).
    adapter_hits = count_adapter_hits(earliest_seqs, k=ADAPTER_PROBE_K)

    # Single-pool primer detection on the earliest round.
    primer_detection = detect_primers(earliest_seqs)
    if primer_detection is None:
        return _build_unable_report(
            read_source=read_source,
            sampling_seed=sampling_seed,
            failure_reason=(
                f"Earliest round has fewer than {DEFAULT_MIN_SEQS_FOR_DETECTION} "
                "sequences — below detection floor"
            ),
            known_adapter_hits=adapter_hits,
        )

    primer_5p_seq = primer_detection.primer_5p.sequence
    primer_3p_seq = primer_detection.primer_3p.sequence

    # Drop primers that match known sequencing adapters (the design:
    # "Exclude from primer candidates"). Track which side was dropped so the
    # adapter_clean signal can faithfully report adapter-trap events even when
    # the OTHER primer survived (previously the signal was
    # 1.0 as long as ANY primer survived, hiding the trap).
    adapter_drop_5p = False
    adapter_drop_3p = False
    if matches_known_adapter_prefix(primer_5p_seq):
        logger.info("Dropping 5' primer candidate %r: matches known adapter", primer_5p_seq)
        primer_5p_seq = None
        adapter_drop_5p = True
    if matches_known_adapter_prefix(primer_3p_seq):
        logger.info("Dropping 3' primer candidate %r: matches known adapter", primer_3p_seq)
        primer_3p_seq = None
        adapter_drop_3p = True

    # Paired-end split detection (overrides 3' from R1 with revcomp of R2's 5').
    # For paired-split, all 3p signals (match_rate, position_consistency,
    # variants, per-round persistence) must be measured against R2 reads at
    # R2's 5' end using ``revcomp(primer_3p)`` as the lookup — that's where the
    # 3' adapter actually appears (previously these were
    # measured against R1's 3' end, where the primer cannot exist in a split
    # library by construction, giving a misleading match_rate_3p ≈ 0).
    has_paired_split = False
    r2_normalized_by_round: dict[int, list[str]] = {}
    primer_3p_lookup: str | None = primer_3p_seq
    threep_is_prefix = False  # default: 3p sits at the 3' end of R1
    if paired_mate_streams:
        # Normalize + subsample all R2 streams up-front so the per-round
        # persistence input has full coverage.
        for r in rounds_sorted:
            r2_norm = _normalize_pool(paired_mate_streams.get(r, []))
            r2_sampled, _ = _subsample(r2_norm, max_reads_per_round, rng)
            r2_normalized_by_round[r] = r2_sampled
        has_paired_split, primer_3p_from_r2 = _detect_paired_split_signals(
            earliest_seqs, r2_normalized_by_round.get(earliest_round, [])
        )
        if has_paired_split and primer_3p_from_r2 is not None:
            primer_3p_seq = primer_3p_from_r2
            # In split mode the 3' adapter appears as revcomp(primer_3p) at the
            # 5' end of R2.
            primer_3p_lookup = reverse_complement(primer_3p_seq)
            threep_is_prefix = True

    # 3p signal context: which reads + which orientation to use.
    # Normal mode: R1 reads, suffix (3' end).
    # Paired-split: R2 reads, prefix (5' end) using revcomp(primer_3p).
    threep_seqs_by_round: dict[int, list[str]] = (
        r2_normalized_by_round if has_paired_split else normalized
    )
    threep_seqs_earliest = threep_seqs_by_round[earliest_round]

    # Per-round POSITION-ANCHORED rates → persistence. Position-anchored
    # (not substring) because a true primer appears AT the flank in every
    # round — substring presence might survive aptamer enrichment.
    position_rates_5p_by_round = [
        _position_consistency(normalized[r], primer_5p_seq, is_prefix=True) for r in rounds_sorted
    ]
    position_rates_3p_by_round = [
        _position_consistency(threep_seqs_by_round[r], primer_3p_lookup, is_prefix=threep_is_prefix)
        for r in rounds_sorted
    ]
    persistence_5p = _persistence_score(position_rates_5p_by_round)
    persistence_3p = _persistence_score(position_rates_3p_by_round)
    persistence = _combine_persistence(persistence_5p, persistence_3p, primer_5p_seq, primer_3p_seq)

    # Earliest-round signals for the report — TWO DISTINCT MEASUREMENTS
    # (previously match_rate_* was aliased to
    # position_consistency_*, double-counting the same evidence in the
    # composite confidence):
    #   match_rate_*           = primer appears anywhere as substring (Hamming ≤ 1)
    #   position_consistency_* = primer appears at the expected flank ± tolerance
    match_rate_5p = _substring_match_rate(earliest_seqs, primer_5p_seq)
    match_rate_3p = _substring_match_rate(threep_seqs_earliest, primer_3p_lookup)

    # Outer-edge core rescue. Runs ONLY when the full-length rate already failed
    # the primer-found threshold, so any library that passes today follows the
    # exact same path it follows now. The rescued rate feeds classification;
    # ``extract`` still trims the full-length flank, because the weak outer
    # positions (adapter/index remnant, T7 start-G) really are in the reads and
    # must be removed. Rescuing therefore cannot move the random-region
    # boundary — it only decides whether the 3' side is usable at all.
    if (
        match_rate_5p <= PRIMER_FOUND_MATCH_RATE_THRESHOLD
        and primer_5p_seq is not None
        and primer_5p_seq == primer_detection.primer_5p.sequence
    ):
        core_5p = _high_support_core(
            primer_5p_seq, primer_detection.primer_5p.supports, is_prefix=True
        )
        if core_5p is not None:
            match_rate_5p = max(match_rate_5p, _substring_match_rate(earliest_seqs, core_5p))

    # Skipped for paired-split libraries: there the 3' lookup is
    # ``revcomp`` of a flank called on R2, so ``primer_detection.primer_3p``
    # (called on R1) is not the sequence whose supports we would be reading.
    if (
        not has_paired_split
        and match_rate_3p <= PRIMER_FOUND_MATCH_RATE_THRESHOLD
        and primer_3p_lookup is not None
        and primer_3p_lookup == primer_detection.primer_3p.sequence
    ):
        core_3p = _high_support_core(
            primer_3p_lookup, primer_detection.primer_3p.supports, is_prefix=False
        )
        if core_3p is not None:
            match_rate_3p = max(match_rate_3p, _substring_match_rate(threep_seqs_earliest, core_3p))

    position_consistency_5p = position_rates_5p_by_round[0] if position_rates_5p_by_round else 0.0
    position_consistency_3p = position_rates_3p_by_round[0] if position_rates_3p_by_round else 0.0

    # Variants — top-K flank fragments. For paired-split, 3p variants come
    # from R2's 5' end (using the revcomp-lookup length).
    p5_len = len(primer_5p_seq) if primer_5p_seq else 0
    p3_len = len(primer_3p_seq) if primer_3p_seq else 0
    p3_lookup_len = len(primer_3p_lookup) if primer_3p_lookup else 0
    variants_5p = _top_k_variants(earliest_seqs, p5_len, is_prefix=True) if p5_len else []
    variants_3p = (
        _top_k_variants(threep_seqs_earliest, p3_lookup_len, is_prefix=threep_is_prefix)
        if p3_lookup_len
        else []
    )

    # N-length distribution. Only meaningful in single-read modes — in
    # paired-end split the full insert spans R1+R2 and cannot be measured
    # from either read alone, so we surface no n-length signal.
    if has_paired_split:
        n_mode, n_dist, n_conf = None, {}, 0.0
    else:
        n_mode, n_dist, n_conf = _n_length_stats(earliest_seqs, p5_len, p3_len)

    # Orientation (always measured on R1 reads — R1's 5' end is where forward
    # vs reverse-strand inversion would appear).
    orientation = _detect_orientation(earliest_seqs, primer_5p_seq, primer_3p_seq)

    # Composite confidence.
    has_round_map = len(normalized) >= 2
    if not has_round_map:
        # UX: single-round / final-pool input is a legitimate
        # and common workflow (HT-SELEX is costly; many depositors
        # sequence only the final enriched pool). It is NOT refused — but
        # the strongest SELEX-specific signal (cross-round persistence) is
        # unavailable, so confidence is capped at MEDIUM (see
        # report._assign_status) and inference leans on within-round
        # signals only. Surface this explicitly so the MEDIUM ceiling is
        # understood rather than mistaken for a calibration wobble.
        logger.warning(
            "single round provided (%d round): cross-round persistence "
            "unavailable (the strongest SELEX-specific primer signal); "
            "confidence is capped at MEDIUM and inference relies on "
            "within-round signals only (primer match rate, flank "
            "position, low-entropy region, adapter blacklist). Verify the "
            "inferred primers before trusting extraction, or pass "
            "--override-primer-5p / --override-primer-3p if you already "
            "know them.",
            len(normalized),
        )
    # adapter_clean = 1.0 only if NO detected candidate was dropped as an
    # adapter (previously this was "we have at least one
    # primer", which hid the adapter trap whenever the other side survived).
    adapter_clean_signal = 0.0 if (adapter_drop_5p or adapter_drop_3p) else 1.0
    signals: dict[str, float | None] = {
        "match_5p": match_rate_5p,
        "match_3p": match_rate_3p,
        "pos_5p": position_consistency_5p,
        "pos_3p": position_consistency_3p,
        "persistence": persistence,
        "n_len": n_conf,
        "adapter_clean": adapter_clean_signal,
    }
    composite = _composite_confidence(signals, has_round_map=has_round_map)

    # Classify (locked decision table).
    classification = _classify(
        match_rate_5p=match_rate_5p,
        match_rate_3p=match_rate_3p,
        n_length_confidence=n_conf,
        has_paired_split=has_paired_split,
        paired_has_overlap=False,  # v0.2 — read merging not implemented
        has_round_map=has_round_map,
        composite_confidence=composite,
        primer_found_threshold=PRIMER_FOUND_MATCH_RATE_THRESHOLD,
        n_length_confident_threshold=N_LENGTH_CONFIDENT_FRACTION,
        unable_to_extract_threshold=UNABLE_TO_EXTRACT_MATCH_RATE,
        status_high_cutoff=STATUS_HIGH_CUTOFF,
        status_medium_cutoff=STATUS_MEDIUM_CUTOFF,
        status_low_cutoff=STATUS_LOW_CUTOFF,
    )

    return LibraryReport(
        primer_5p=primer_5p_seq,
        primer_3p=primer_3p_seq,
        variants_5p=variants_5p,
        variants_3p=variants_3p,
        known_adapter_hits=adapter_hits,
        extraction_mode=classification.extraction_mode,
        full_insert_recovered=classification.full_insert_recovered,
        read_source=read_source,
        required_action=classification.required_action,
        orientation=orientation,
        n_length_mode=n_mode,
        n_length_distribution=n_dist,
        n_length_confidence=n_conf,
        match_rate_5p=match_rate_5p,
        match_rate_3p=match_rate_3p,
        position_consistency_5p=position_consistency_5p,
        position_consistency_3p=position_consistency_3p,
        read_fraction_used_for_inference=read_fraction,
        sampling_seed=sampling_seed,
        confidence=composite,
        status=classification.status,
        failure_reason=classification.failure_reason,
    )
