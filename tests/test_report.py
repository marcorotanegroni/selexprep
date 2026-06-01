"""Behavior-based tests for ``LibraryReport`` schema + inference pipeline.

Each test asserts on the **behavior** of ``compute_library_report`` —
``extraction_mode``, ``required_action``, ``status`` — never on threshold
constants. When the Codex calibration review tunes numbers in
``selexprep.library.detect``, these tests must stay green.

The test pools rely on the existing ``_synthetic_pool`` helper from
``tests/test_detect.py``; it is redefined here so the two test modules
stay independent.

One test per row of the locked classification table (plan lines 300-309)
plus the edge cases listed in the Phase 2 plan (status cap, adapter
blacklist demotion, orientation, U→T, determinism, schema smoke).
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from pydantic import ValidationError

from selexprep.library.adapters import KNOWN_ADAPTERS, reverse_complement
from selexprep.library.detect import compute_library_report
from selexprep.library.report import (
    LibraryReport,
    read_library_report_json,
    write_library_report_json,
)

# ---------------------------------------------------------------------------
# Synthetic pools — local duplicate of the helper in test_detect.py.
# ---------------------------------------------------------------------------


PRIMER_5P_T7 = "GGTAATACGACTCACTATAGGG"  # 22 nt, T7 promoter — a stock SELEX primer
PRIMER_3P_CCAT = "CCATGCATGCATGCATGCAT"  # 20 nt, arbitrary but distinct from PRIMER_5P_T7


def _synthetic_pool(
    primer_5p: str | None,
    primer_3p: str | None,
    n: int = 1000,
    random_len: int = 30,
    seed_offset: int = 0,
) -> list[str]:
    """primer_5p + random_region + primer_3p, deterministic random region.

    ``seed_offset`` shifts the random-region generator so different pools
    don't collide on identical random sequences.
    """
    bases = "ACGT"
    seqs: list[str] = []
    for i in range(n):
        rand = "".join(bases[((i + seed_offset) * 7 + j * 13) % 4] for j in range(random_len))
        seqs.append((primer_5p or "") + rand + (primer_3p or ""))
    return seqs


def _three_round_pool(
    primer_5p: str | None,
    primer_3p: str | None,
    *,
    n: int = 1000,
    random_len: int = 30,
) -> dict[int, list[str]]:
    """Same primers across 3 rounds (steady-state — what cross-round persistence expects)."""
    return {
        r: _synthetic_pool(primer_5p, primer_3p, n=n, random_len=random_len, seed_offset=r * 1000)
        for r in range(3)
    }


# ===========================================================================
# Classification-table coverage (locked plan lines 300-309)
# ===========================================================================


def test_both_primers_single_read() -> None:
    """Row 1: both 5' and 3' primers detected on the same read."""
    pools = _three_round_pool(PRIMER_5P_T7, PRIMER_3P_CCAT)
    report = compute_library_report(pools, read_source="R1")

    assert report.extraction_mode == "BOTH_PRIMERS_SINGLE_READ"
    assert report.full_insert_recovered is True
    assert report.required_action == "NONE"
    assert report.primer_5p is not None
    assert report.primer_3p is not None
    assert report.status in ("HIGH", "MEDIUM")  # threshold-agnostic


def test_five_prime_only_clean_n() -> None:
    """Row 2: only 5' primer + sharply peaked N-length."""
    pools = _three_round_pool(PRIMER_5P_T7, primer_3p=None, random_len=30)
    report = compute_library_report(pools, read_source="R1")

    assert report.extraction_mode == "FIVE_PRIME_ONLY"
    assert report.full_insert_recovered is False
    assert report.required_action == "NONE"
    assert report.primer_5p is not None
    assert report.primer_3p is None


def test_five_prime_only_smeared_n() -> None:
    """Row 3: only 5' primer detected but N-length distribution is too smeared
    to assert a clean random-region length → UNABLE_TO_EXTRACT."""
    # Mix five random_lens, 200 reads each → mode count is 200/1000 = 0.20,
    # well below the 0.8 confidence threshold.
    pools: dict[int, list[str]] = {}
    for r in range(3):
        round_seqs: list[str] = []
        for offset, rl in enumerate((20, 25, 30, 35, 40)):
            round_seqs.extend(
                _synthetic_pool(
                    PRIMER_5P_T7,
                    primer_3p=None,
                    n=200,
                    random_len=rl,
                    seed_offset=r * 1000 + offset * 100,
                )
            )
        pools[r] = round_seqs

    report = compute_library_report(pools, read_source="R1")

    assert report.extraction_mode == "UNABLE_TO_EXTRACT"
    assert report.required_action == "MANUAL_PRIMERS_REQUIRED"
    assert report.status == "UNABLE_TO_INFER"


def test_three_prime_only_clean_n() -> None:
    """Row 4: only 3' primer + sharply peaked N-length."""
    pools = _three_round_pool(primer_5p=None, primer_3p=PRIMER_3P_CCAT, random_len=30)
    report = compute_library_report(pools, read_source="R1")

    assert report.extraction_mode == "THREE_PRIME_ONLY"
    assert report.full_insert_recovered is False
    assert report.required_action == "NONE"
    assert report.primer_5p is None
    assert report.primer_3p is not None


def test_three_prime_only_smeared_n() -> None:
    """Row 5: only 3' primer detected but smeared N-length → UNABLE_TO_EXTRACT."""
    pools: dict[int, list[str]] = {}
    for r in range(3):
        round_seqs: list[str] = []
        for offset, rl in enumerate((20, 25, 30, 35, 40)):
            round_seqs.extend(
                _synthetic_pool(
                    primer_5p=None,
                    primer_3p=PRIMER_3P_CCAT,
                    n=200,
                    random_len=rl,
                    seed_offset=r * 1000 + offset * 100,
                )
            )
        pools[r] = round_seqs

    report = compute_library_report(pools, read_source="R1")

    assert report.extraction_mode == "UNABLE_TO_EXTRACT"
    assert report.required_action == "MANUAL_PRIMERS_REQUIRED"


def test_paired_end_split_no_overlap() -> None:
    """Row 6: R1 carries 5' primer, R2 carries revcomp(3' primer); no merger."""
    # R1: 5' primer + long random tail (3' primer never reaches the read end).
    r1_pools = _three_round_pool(PRIMER_5P_T7, primer_3p=None, random_len=80)
    # R2: revcomp(3' primer) at start + long random tail.
    rc_3p = reverse_complement(PRIMER_3P_CCAT)
    r2_pools = _three_round_pool(rc_3p, primer_3p=None, random_len=80)

    report = compute_library_report(r1_pools, read_source="R1_AND_R2", paired_mate_streams=r2_pools)

    assert report.extraction_mode == "PAIRED_END_SPLIT_PRIMERS"
    assert report.full_insert_recovered is False
    assert report.required_action == "READ_MERGING_RECOMMENDED"


def test_paired_end_split_prefers_r2_over_conflicting_r1_suffix() -> None:
    """Paired-end split must not be blocked by a strong technical R1 suffix.

    PRJNA883192 has R1 evidence for the 5' insert constant and R2 evidence
    for the 3' insert constant, but R1 also ends in a strong unrelated
    technical suffix. Its R2 primer core is followed by a C-biased random
    shoulder; paired-split must use the high-support R2 core, not the
    overextended shoulder, as the biologically correct source for the insert
    3' constant in this layout.
    """
    primer_5p = "ATGCCATCCTACCAAC"
    primer_3p = "GAGCTCTGAACTCGA"  # PRJNA883192 read-level R2-derived flank
    technical_r1_suffix = "TGAACTCCAGTCACCGAATAATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAAA"

    r1_pools = _three_round_pool(primer_5p, technical_r1_suffix, random_len=80)
    r2_core = reverse_complement(primer_3p)
    r2_pools: dict[int, list[str]] = {}
    for r in range(3):
        seqs: list[str] = []
        for i in range(1000):
            biased_tail = [
                "C" if i % 5 != 0 else "A",
                "C" if i % 5 < 3 else "G",
                "C" if i % 5 < 3 else "T",
                "C" if i % 7 < 4 else "A",
            ]
            random_tail = "".join("ACGT"[((i + r * 1000) * 7 + j * 13) % 4] for j in range(76))
            seqs.append(r2_core + "".join(biased_tail) + random_tail)
        r2_pools[r] = seqs

    report = compute_library_report(r1_pools, read_source="R1_AND_R2", paired_mate_streams=r2_pools)

    assert report.extraction_mode == "PAIRED_END_SPLIT_PRIMERS"
    assert report.required_action == "READ_MERGING_RECOMMENDED"
    assert report.primer_5p == primer_5p
    assert report.primer_3p == primer_3p
    assert report.match_rate_3p > 0.95
    assert report.position_consistency_3p > 0.95


def test_paired_r2_does_not_demote_single_read_when_r1_3p_agrees() -> None:
    """If R1 already contains the same 3' primer that R2 implies, keep Row 1.

    This prevents paired-end support from unnecessarily changing a full
    single-read library into ``PAIRED_END_SPLIT_PRIMERS``.
    """
    r1_pools = _three_round_pool(PRIMER_5P_T7, PRIMER_3P_CCAT, random_len=30)
    r2_pools = _three_round_pool(reverse_complement(PRIMER_3P_CCAT), primer_3p=None, random_len=80)

    report = compute_library_report(r1_pools, read_source="R1_AND_R2", paired_mate_streams=r2_pools)

    assert report.extraction_mode == "BOTH_PRIMERS_SINGLE_READ"
    assert report.full_insert_recovered is True
    assert report.required_action == "NONE"
    assert report.primer_3p == PRIMER_3P_CCAT


@pytest.mark.xfail(
    strict=True,
    reason=(
        "v0.2: read merging not yet implemented — when this test starts "
        "passing, the feature is ready, remove this marker"
    ),
)
def test_paired_end_split_with_overlap_v02() -> None:
    """Row 7 of the locked classification table: paired-end with detectable
    overlap → ``full_insert_recovered=True``, ``required_action=NONE``.

    Pinned via ``strict=True`` xfail: v0.1 cannot satisfy this assertion
    (always falls back to Row 6 — split-primer with no overlap), so the
    test currently fails as expected. When v0.2 ships read merging, this
    test will start passing → strict mode promotes XPASS to a CI failure,
    forcing the marker to be removed. Prevents the classic "xfail marker
    rotting after the feature lands" pitfall.
    """
    r1_pools = _three_round_pool(PRIMER_5P_T7, primer_3p=None, random_len=80)
    rc_3p = reverse_complement(PRIMER_3P_CCAT)
    r2_pools = _three_round_pool(rc_3p, primer_3p=None, random_len=80)

    report = compute_library_report(r1_pools, read_source="R1_AND_R2", paired_mate_streams=r2_pools)
    # v0.2 behavior:
    assert report.full_insert_recovered is True
    assert report.required_action == "NONE"


def test_both_match_rates_low_returns_unable() -> None:
    """Row 8: neither flank reaches the detection floor → UNABLE_TO_INFER."""
    # Pure-random pool (no shared flanks). detect_flank cannot find any
    # consensus prefix or suffix.
    pools = _three_round_pool(primer_5p=None, primer_3p=None, random_len=80)
    report = compute_library_report(pools, read_source="R1")

    assert report.extraction_mode == "UNABLE_TO_EXTRACT"
    assert report.required_action == "MANUAL_PRIMERS_REQUIRED"
    assert report.status == "UNABLE_TO_INFER"
    assert report.primer_5p is None
    assert report.primer_3p is None
    assert report.failure_reason is not None


# ===========================================================================
# Edge cases (Phase 2 plan)
# ===========================================================================


def test_status_capped_medium_no_round_map() -> None:
    """Single-round input → persistence is None → status never HIGH."""
    pools = {0: _synthetic_pool(PRIMER_5P_T7, PRIMER_3P_CCAT, n=1000)}
    report = compute_library_report(pools, read_source="R1")

    # With persistence dropped, status caps at MEDIUM per locked plan line 289.
    assert report.status in ("MEDIUM", "LOW", "UNABLE_TO_INFER")


def test_single_round_caps_status_and_warns(caplog: pytest.LogCaptureFixture) -> None:
    """Phase 6b.8 single-round UX: explicit anti-regression of the MEDIUM cap.

    A strong, clean single-round pool would score HIGH on within-round
    signals alone — but the cap in ``report._assign_status`` must force it
    to ≤ MEDIUM because cross-round persistence (the strongest
    SELEX-specific signal) is unavailable. This test pins ``!= "HIGH"``
    explicitly so a future refactor of ``_assign_status`` can't silently
    drop the cap, and confirms the user-facing warning is emitted.
    """
    pools = {0: _synthetic_pool(PRIMER_5P_T7, PRIMER_3P_CCAT, n=1000)}
    with caplog.at_level("WARNING", logger="selexprep.library.detect"):
        report = compute_library_report(pools, read_source="R1")

    # Hard anti-regression: single round must NEVER be HIGH.
    assert report.status != "HIGH"

    warnings = [r.getMessage() for r in caplog.records if r.levelname == "WARNING"]
    assert any(
        "single round" in m and "capped at MEDIUM" in m and "persistence" in m for m in warnings
    ), f"expected single-round MEDIUM-cap warning, got: {warnings}"


def test_adapter_blacklist_demotes_truseq_candidate() -> None:
    """If the detected 5' primer matches a known sequencing adapter,
    the candidate is dropped and ``known_adapter_hits`` reflects it."""
    truseq = KNOWN_ADAPTERS["TRUSEQ_R1"]
    # Use a 5' "primer" whose first 13 nt match TruSeq R1 exactly.
    truseq_like_primer = truseq + "AGGGGGT"  # 13 + 7 = 20 nt
    pools = _three_round_pool(truseq_like_primer, primer_3p=None, random_len=30)

    report = compute_library_report(pools, read_source="R1")

    # Adapter was present → recorded.
    assert report.known_adapter_hits["TRUSEQ_R1"] > 0
    # Detected candidate dropped — primer_5p should be None.
    assert report.primer_5p is None


def test_orientation_mixed_when_revcomp_present() -> None:
    """80/20 mix of forward + reverse orientation reads → orientation MIXED.

    Pure 50/50 is degenerate (no consensus at any flank); 80/20 keeps the
    primer detectable while still surfacing the reverse-strand minority.
    """
    rc_3p = reverse_complement(PRIMER_3P_CCAT)
    rc_5p = reverse_complement(PRIMER_5P_T7)

    def mixed_round(seed: int) -> list[str]:
        forward = _synthetic_pool(PRIMER_5P_T7, PRIMER_3P_CCAT, n=800, seed_offset=seed)
        reverse = _synthetic_pool(rc_3p, rc_5p, n=200, seed_offset=seed + 5000)
        return forward + reverse

    pools = {r: mixed_round(seed=r * 10000) for r in range(3)}
    report = compute_library_report(pools, read_source="R1")

    assert report.orientation == "MIXED"


def test_rna_primers_reported_as_dna() -> None:
    """Reads with U's (RNA notation) report DNA-only primers."""
    # T7 promoter with the natural T replaced by U at one position to
    # exercise the U→T normalizer.
    primer_5p_rna = PRIMER_5P_T7.replace("T", "U")
    pools = _three_round_pool(primer_5p_rna, primer_3p=None, random_len=30)
    report = compute_library_report(pools, read_source="R1")

    assert report.primer_5p is not None
    assert set(report.primer_5p) <= set("ACGT"), (
        f"primer_5p contains non-DNA bases: {report.primer_5p!r}"
    )


def test_deterministic_json_serialization(tmp_path: Path) -> None:
    """Two runs with identical inputs + sampling_seed produce identical JSON sha256."""
    pools = _three_round_pool(PRIMER_5P_T7, PRIMER_3P_CCAT)

    report_a = compute_library_report(pools, read_source="R1", sampling_seed=42)
    report_b = compute_library_report(pools, read_source="R1", sampling_seed=42)

    path_a = tmp_path / "report_a.json"
    path_b = tmp_path / "report_b.json"
    hash_a = write_library_report_json(report_a, path_a)
    hash_b = write_library_report_json(report_b, path_b)

    assert hash_a == hash_b
    assert path_a.read_bytes() == path_b.read_bytes()


def test_json_roundtrip_preserves_report(tmp_path: Path) -> None:
    pools = _three_round_pool(PRIMER_5P_T7, PRIMER_3P_CCAT)
    original = compute_library_report(pools, read_source="R1", sampling_seed=42)

    path = tmp_path / "report.json"
    write_library_report_json(original, path)
    loaded = read_library_report_json(path)

    assert loaded == original


def test_library_report_model_is_frozen() -> None:
    """LibraryReport is immutable; attempted mutation raises ValidationError."""
    pools = _three_round_pool(PRIMER_5P_T7, PRIMER_3P_CCAT)
    report = compute_library_report(pools, read_source="R1")

    with pytest.raises(ValidationError):
        report.confidence = 0.0  # type: ignore[misc]


def test_library_report_rejects_extra_fields() -> None:
    """Schema is closed (extra='forbid') — typos in field names fail loudly."""
    with pytest.raises(ValidationError):
        LibraryReport(
            primer_5p=None,
            primer_3p=None,
            variants_5p=[],
            variants_3p=[],
            known_adapter_hits={},
            extraction_mode="UNABLE_TO_EXTRACT",
            full_insert_recovered=False,
            read_source="UNKNOWN",
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
            sampling_seed=42,
            confidence=0.0,
            status="UNABLE_TO_INFER",
            failure_reason=None,
            mystery_field="oops",  # type: ignore[call-arg]
        )


def test_n_length_distribution_serialized_in_numeric_order(tmp_path: Path) -> None:
    """Deterministic JSON puts int-keyed dicts in numeric order, not lexical
    ("10" before "100" before "20") — readability + reproducibility."""
    pools: dict[int, list[str]] = {}
    for r in range(3):
        round_seqs: list[str] = []
        for offset, rl in enumerate((10, 20, 30, 100)):  # mix small + 3-digit lengths
            round_seqs.extend(
                _synthetic_pool(
                    PRIMER_5P_T7,
                    PRIMER_3P_CCAT,
                    n=250,
                    random_len=rl,
                    seed_offset=r * 1000 + offset * 100,
                )
            )
        pools[r] = round_seqs

    report = compute_library_report(pools, read_source="R1")
    path = tmp_path / "report.json"
    write_library_report_json(report, path)

    payload = json.loads(path.read_text())
    keys = list(payload["n_length_distribution"].keys())
    # Numeric order, not lexical: "20" comes before "100".
    assert keys == sorted(keys, key=int)


def test_empty_input_returns_unable() -> None:
    report = compute_library_report({}, read_source="R1")
    assert report.extraction_mode == "UNABLE_TO_EXTRACT"
    assert report.status == "UNABLE_TO_INFER"
    assert report.failure_reason is not None


def test_below_detection_floor_returns_unable() -> None:
    """Fewer than DEFAULT_MIN_SEQS_FOR_DETECTION (500) reads → cannot detect."""
    pools = {0: _synthetic_pool(PRIMER_5P_T7, PRIMER_3P_CCAT, n=100)}
    report = compute_library_report(pools, read_source="R1")

    assert report.extraction_mode == "UNABLE_TO_EXTRACT"
    assert report.failure_reason is not None
    assert "below detection floor" in report.failure_reason


# ===========================================================================
# Codex pass 1 regressions (2026-05-20)
# ===========================================================================
# Three bug fixes that survived to v0.1 RC and were caught in the second
# Codex review of Phase 2:
#   1. match_rate_* was aliased to position_consistency_* (double-count).
#   2. adapter_clean_signal ignored the drop flag (hid the adapter trap).
#   3. paired-split mode measured 3p signals against R1 instead of R2.
# These tests assert the post-fix behavior so the bugs cannot regress.


def test_match_rate_distinct_from_position_consistency() -> None:
    """Codex pass 1 regression #1: match_rate (substring anywhere, Hamming ≤ 1)
    and position_consistency (substring at flank ± tolerance) must NOT be
    aliased. Place the primer 10 nt past the read start; with tolerance=3
    the flank check fails but the substring check passes."""
    from selexprep.library.detect import (
        POSITION_CONSISTENCY_TOLERANCE,
        _position_consistency,
        _substring_match_rate,
    )

    primer = PRIMER_5P_T7  # 22 nt
    # Each read = 10 nt offset + primer + 4 nt tail. Primer starts at position 10.
    seqs = [f"ACGTACGTAC{primer}TTTT" for _ in range(100)]

    pos = _position_consistency(
        seqs, primer, is_prefix=True, tolerance=POSITION_CONSISTENCY_TOLERANCE
    )
    match = _substring_match_rate(seqs, primer)

    # Position-anchored fails: primer at position 10 > tolerance (3).
    assert pos < 0.05, f"position_consistency should be ~0, got {pos}"
    # Substring-anywhere passes: primer present in every read.
    assert match > 0.95, f"match_rate should be ~1, got {match}"


def test_adapter_clean_flag_demotes_confidence_when_primer_dropped() -> None:
    """Codex pass 1 regression #2: when a detected primer matches a known
    sequencing adapter and gets dropped, adapter_clean must register the
    event (0.0) even if the OTHER primer survived. Previously the signal
    was "1.0 if any primer survived", which hid the adapter trap."""
    truseq = KNOWN_ADAPTERS["TRUSEQ_R1"]  # "AGATCGGAAGAGC", 13 nt
    # Pool A: clean primers (no adapter trap).
    clean_pools = _three_round_pool(PRIMER_5P_T7, PRIMER_3P_CCAT)
    clean = compute_library_report(clean_pools, read_source="R1")

    # Pool B: 5' primer's first 13 nt match TruSeq R1 → adapter-drop fires
    # on the 5' side; 3' primer survives.
    truseq_primer = truseq + "AGGGGGT"  # 13 + 7 = 20 nt
    trap_pools = _three_round_pool(truseq_primer, PRIMER_3P_CCAT)
    trap = compute_library_report(trap_pools, read_source="R1")

    # 5' candidate dropped + recorded.
    assert trap.primer_5p is None
    assert trap.known_adapter_hits["TRUSEQ_R1"] > 0
    # The 3' primer survives — but adapter_clean still registers the trap,
    # so the trap report has strictly lower composite confidence than the
    # clean report (the bug would have produced equal confidence here).
    assert trap.confidence < clean.confidence, (
        f"trap confidence ({trap.confidence:.3f}) should be lower than clean "
        f"({clean.confidence:.3f}) — adapter drop must demote adapter_clean"
    )


def test_paired_split_match_rate_3p_reflects_r2_not_r1() -> None:
    """Codex pass 1 regression #3: in paired-end split mode, match_rate_3p
    + position_consistency_3p + variants_3p must reflect R2 evidence
    (where ``revcomp(primer_3p)`` actually appears at R2's 5' end) — NOT
    R1's 3' end, where the 3' adapter cannot exist by construction."""
    # R1: only 5' primer at start, long random tail (no 3' primer in R1).
    r1_pools = _three_round_pool(PRIMER_5P_T7, primer_3p=None, random_len=80)
    # R2: revcomp(3' primer) at start + long random tail.
    rc_3p = reverse_complement(PRIMER_3P_CCAT)
    r2_pools = _three_round_pool(rc_3p, primer_3p=None, random_len=80)

    report = compute_library_report(r1_pools, read_source="R1_AND_R2", paired_mate_streams=r2_pools)

    assert report.extraction_mode == "PAIRED_END_SPLIT_PRIMERS"
    # Bug would have made these ≈ 0 (measured against R1 where they can't
    # appear). Fix makes them ≈ 1 (measured against R2's 5' end).
    assert report.match_rate_3p > 0.5, (
        f"match_rate_3p should be high (R2-measured), got {report.match_rate_3p:.3f}"
    )
    assert report.position_consistency_3p > 0.5, (
        f"position_consistency_3p should be high (R2-measured), "
        f"got {report.position_consistency_3p:.3f}"
    )
    assert report.variants_3p, "variants_3p should be non-empty (R2's 5' fragments)"
