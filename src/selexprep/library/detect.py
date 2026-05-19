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
import random
import statistics
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from selexprep.library.adapters import (
    KNOWN_ADAPTERS,
    KNOWN_ADAPTERS_RC,
    count_adapter_hits,
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


# ===========================================================================
# Phase 2 — cross-round LibraryReport orchestration
# ===========================================================================
#
# Everything below this banner is Phase 2 work (the LibraryReport pipeline);
# everything above is the Phase 1 single-pool flank detector. The split is
# intentional: the locked plan keeps the Phase 1 functions as the bare
# algorithmic core that Phase 2 wraps with calibration, persistence, and
# classification. Phase 1 callers (tests, future ad-hoc tools) keep working.
#
# All calibration constants below are PLACEHOLDERS pending Codex review
# (rate-limited 2026-05-19 → 2026-05-26). Behavior-based tests mean tuning
# numbers does NOT break the test suite. Search for `CALIBRATION-TODO` to
# find every pending value.

# ---------------------------------------------------------------------------
# Calibration constants — CALIBRATION-TODO
# ---------------------------------------------------------------------------

# CALIBRATION-TODO: locked plan line 302 ("> 0.7"); Codex sanity-checks
# against published HT-SELEX figures (Tolle 2014, Hoinka 2015, AptaSUITE
# 2017, EasyDIVER+ 2025).
PRIMER_FOUND_MATCH_RATE_THRESHOLD = 0.7

# CALIBRATION-TODO: locked plan lines 303, 305 (n_length_confidence > 0.8).
N_LENGTH_CONFIDENT_FRACTION = 0.8

# CALIBRATION-TODO: locked plan line 309 ("Both match rates < 0.4").
UNABLE_TO_EXTRACT_MATCH_RATE = 0.4

# CALIBRATION-TODO: not in locked plan; absorbs residual adapter readthrough
# that survived upstream trimming. Codex confirms the ±2 nt window.
POSITION_CONSISTENCY_TOLERANCE = 2

# CALIBRATION-TODO: not in locked plan; Codex sets cutoffs after benchmark
# recovery numbers are in.
STATUS_HIGH_CUTOFF = 0.80
STATUS_MEDIUM_CUTOFF = 0.60
STATUS_LOW_CUTOFF = 0.30

# CALIBRATION-TODO: not in locked plan. Weighted sum, weights sum to 1.0.
# Two regimes: with vs without a round map. Persistence weight gets
# redistributed to the other terms when no round map is available.
COMPOSITE_WEIGHTS = {
    "match_5p": 0.20,
    "match_3p": 0.20,
    "pos_5p": 0.10,
    "pos_3p": 0.10,
    "persistence": 0.20,
    "n_len": 0.10,
    "adapter_clean": 0.10,
}
COMPOSITE_WEIGHTS_NO_ROUND_MAP = {
    "match_5p": 0.30,
    "match_3p": 0.30,
    "pos_5p": 0.15,
    "pos_3p": 0.15,
    "persistence": 0.0,
    "n_len": 0.05,
    "adapter_clean": 0.05,
}

# CALIBRATION-TODO: not in locked plan. < 5% reads reversed → FORWARD,
# 5-95% → MIXED, > 95% → REVERSE.
ORIENTATION_REVERSED_FORWARD_MAX = 0.05
ORIENTATION_REVERSED_REVERSE_MIN = 0.95

# Top-K variants surfaced in the LibraryReport (locked plan line 297: K=3).
VARIANTS_TOP_K = 3

# Adapter probe length for `count_adapter_hits` substring scan.
ADAPTER_PROBE_K = 13


# ---------------------------------------------------------------------------
# Per-signal helpers
# ---------------------------------------------------------------------------


def _normalize_u_to_t(seq: str) -> str:
    """Convert RNA primer notation to DNA (locked plan line 296)."""
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


def _persistence_score(match_rates_per_round: list[float]) -> float | None:
    """Cross-round persistence as ``1 - clip(stdev/mean, 0, 1)``.

    Returns ``None`` if fewer than 2 rounds are available (caller falls
    back to no-round-map composite weights and caps ``status`` at MEDIUM).
    Returns ``0.0`` if the mean rate is below 10% (no signal to be
    persistent about).
    """
    if len(match_rates_per_round) < 2:
        return None
    mean = statistics.mean(match_rates_per_round)
    if mean < 0.1:
        return 0.0
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


def _matches_known_adapter(primer: str | None, k: int = ADAPTER_PROBE_K) -> bool:
    """True if ``primer``'s first ``k`` bases coincide with a known adapter or its RC."""
    if not primer:
        return False
    probe = primer.upper()[:k]
    return any(adapter[:k] == probe for adapter in KNOWN_ADAPTERS.values()) or any(
        adapter_rc[:k] == probe for adapter_rc in KNOWN_ADAPTERS_RC.values()
    )


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
    is True when R1 carries a strong 5' primer and weak 3' primer, AND R2
    carries a strong 5' primer (which when reverse-complemented gives the
    real 3' primer of the insert).

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
    if r1_5p_strong and not r1_3p_strong and r2_5p_strong and r2_det.primer_5p.sequence is not None:
        return True, reverse_complement(r2_det.primer_5p.sequence)
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

    Algorithm (locked plan §Phase 2, lines 287-313):

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

    # Single-pool primer detection (Phase 1 algorithm) on the earliest round.
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

    # Drop primers that match known sequencing adapters (locked plan line 291:
    # "Exclude from primer candidates").
    if _matches_known_adapter(primer_5p_seq):
        logger.info("Dropping 5' primer candidate %r: matches known adapter", primer_5p_seq)
        primer_5p_seq = None
    if _matches_known_adapter(primer_3p_seq):
        logger.info("Dropping 3' primer candidate %r: matches known adapter", primer_3p_seq)
        primer_3p_seq = None

    # Paired-end split detection (overrides 3' from R1 with revcomp of R2's 5').
    has_paired_split = False
    if paired_mate_streams:
        # Use the earliest round's R2 stream; normalize first.
        r2_earliest = _normalize_pool(paired_mate_streams.get(earliest_round, []))
        r2_sampled, _ = _subsample(r2_earliest, max_reads_per_round, rng)
        has_paired_split, primer_3p_from_r2 = _detect_paired_split_signals(
            earliest_seqs, r2_sampled
        )
        if has_paired_split and primer_3p_from_r2 is not None:
            primer_3p_seq = primer_3p_from_r2

    # Per-round match rates → persistence (uses the detected primer; does NOT
    # re-detect per round, which is what makes this a true cross-round signal).
    match_rates_5p_by_round = [
        _position_consistency(normalized[r], primer_5p_seq, is_prefix=True) for r in rounds_sorted
    ]
    match_rates_3p_by_round = [
        _position_consistency(normalized[r], primer_3p_seq, is_prefix=False) for r in rounds_sorted
    ]
    persistence_5p = _persistence_score(match_rates_5p_by_round)
    persistence_3p = _persistence_score(match_rates_3p_by_round)
    persistence = _combine_persistence(persistence_5p, persistence_3p, primer_5p_seq, primer_3p_seq)

    # Earliest-round signals for the report.
    match_rate_5p = match_rates_5p_by_round[0] if match_rates_5p_by_round else 0.0
    match_rate_3p = match_rates_3p_by_round[0] if match_rates_3p_by_round else 0.0
    position_consistency_5p = match_rate_5p
    position_consistency_3p = match_rate_3p

    # Variants — top-K prefixes/suffixes at the detected primer's length.
    p5_len = len(primer_5p_seq) if primer_5p_seq else 0
    p3_len = len(primer_3p_seq) if primer_3p_seq else 0
    variants_5p = _top_k_variants(earliest_seqs, p5_len, is_prefix=True) if p5_len else []
    variants_3p = _top_k_variants(earliest_seqs, p3_len, is_prefix=False) if p3_len else []

    # N-length distribution.
    n_mode, n_dist, n_conf = _n_length_stats(earliest_seqs, p5_len, p3_len)

    # Orientation.
    orientation = _detect_orientation(earliest_seqs, primer_5p_seq, primer_3p_seq)

    # Composite confidence.
    has_round_map = len(normalized) >= 2
    adapter_clean_signal = 1.0 if (primer_5p_seq is not None or primer_3p_seq is not None) else 0.0
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
