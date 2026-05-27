"""Unit tests for selexprep.fetch.metadata (round-parser cascade)."""

from __future__ import annotations

from pathlib import Path

import pytest

from selexprep.fetch.metadata import (
    RoundRecord,
    apply_seed_overrides,
    extract_round_count_from_abstract,
    parse_round,
    summarise_bioproject,
)

# ----- L1: structured sample_attributes -----


def test_l1_structured_round_attribute_high_confidence() -> None:
    r = parse_round("SRR1", sample_attributes={"selex_round": "3"})
    assert r.round_number == 3
    assert r.confidence == "HIGH"
    assert r.source_field == "sample_attributes"
    assert r.is_unassigned is False


def test_l1_handles_various_attribute_keys() -> None:
    for key in ["round", "Cycle", "iteration", "selection_round", "round_number"]:
        r = parse_round("SRR1", sample_attributes={key: "5"})
        assert r.round_number == 5
        assert r.confidence == "HIGH"


def test_l1_ignores_non_numeric_value() -> None:
    r = parse_round("SRR1", sample_attributes={"round": "five"})
    assert r.round_number is None
    assert r.confidence == "NONE"


# ----- L2: sample_title patterns -----


def test_l2_sample_title_round_word() -> None:
    r = parse_round("SRR1", sample_title="anti-VEGF Round 7")
    assert r.round_number == 7
    assert r.confidence == "HIGH"
    assert r.source_field == "sample_title"


def test_l2_sample_title_R_digit_boundary() -> None:
    r = parse_round("SRR1", sample_title="ThrombinSELEX R12 library")
    assert r.round_number == 12


def test_l2_sample_title_cycle() -> None:
    r = parse_round("SRR1", sample_title="Cell-SELEX cycle 5")
    assert r.round_number == 5


def test_l2_conflicting_numbers_downgrade_to_medium() -> None:
    r = parse_round("SRR1", sample_title="Round 3 / Cycle 7")
    assert r.confidence == "MEDIUM"
    assert r.round_candidates == [3, 7]
    # Phase 6b.4 audit refactor: genuine ambiguity is flagged via
    # ``is_unassigned`` (computed from ``len(round_candidates) > 1``)
    # rather than the old per-record ``needs_manual_review`` field.
    assert r.is_unassigned is True


# ----- L3: library_name / experiment_title -----


def test_l3_library_name_medium_confidence() -> None:
    r = parse_round("SRR1", library_name="VEGF_R6_library")
    assert r.round_number == 6
    assert r.confidence == "MEDIUM"
    assert r.source_field == "library_name"
    # Phase 6b.4 audit refactor: a single unambiguous L3 parse is NOT
    # unassigned. Previously the field ``needs_manual_review`` was True
    # here, which caused ``FetchPlan.none_confidence_runs`` to refuse
    # fetch on accessions whose round signals only lived in library_name.
    assert r.is_unassigned is False


def test_l3_falls_through_to_experiment_then_design() -> None:
    r = parse_round("SRR1", experiment_title="SELEX cycle 4 enrichment")
    assert r.round_number == 4
    assert r.source_field == "experiment_title"
    assert r.is_unassigned is False


# ----- Phase 6b.4 audit-pilot regression tests -----
#
# The N=30 Tier 2 audit pilot surfaced a policy bug: MEDIUM single-match
# parses (RAPT26-2R, SPa19-1R, R00_N16) were treated as unassigned and
# triggered FETCH_REFUSED on whole accessions. These tests pin the fixed
# behavior against the empirical cases Codex found during the pilot's
# per-accession fetch_metadata.json inspection.


def test_audit_pilot_rapt26_library_name_does_not_need_review() -> None:
    """PRJDB19138: library_name='RAPT26-2R' must parse to round 2 without
    being treated as unassigned. Pre-fix this caused the whole accession
    to be refused; post-fix it admits cleanly to rounds.tsv."""
    r = parse_round("DRR618526", library_name="RAPT26-2R")
    assert r.round_number == 2
    assert r.confidence == "MEDIUM"
    assert r.source_field == "library_name"
    assert r.is_unassigned is False


def test_audit_pilot_spa19_library_name_does_not_need_review() -> None:
    """PRJDB40016: library_name='SPa19-1R' → round 1, not unassigned."""
    r = parse_round("DRR895630", library_name="SPa19-1R")
    assert r.round_number == 1
    assert r.confidence == "MEDIUM"
    assert r.is_unassigned is False


def test_audit_pilot_r00_n16_sample_title_does_not_need_review() -> None:
    """PRJEB98610: sample_title='R00_N16' → round 0 (HIGH confidence —
    sample_title is L2; this row was already HIGH and is the cleanest
    sanity check that no L2 single-match record was ever unassigned)."""
    r = parse_round("ERR15669095", sample_title="R00_N16")
    assert r.round_number == 0
    assert r.confidence == "HIGH"
    assert r.source_field == "sample_title"
    assert r.is_unassigned is False


def test_audit_pilot_selex_round26_n40_does_not_need_review() -> None:
    """Synthetic regression for ``SELEX_round26_N40`` (a generic library_name
    pattern seen across multiple SELEX deposits). Single unambiguous
    parse, MEDIUM confidence by virtue of being L3, not unassigned."""
    r = parse_round("SRRX", library_name="SELEX_round26_N40")
    assert r.round_number == 26
    assert r.confidence == "MEDIUM"
    assert r.is_unassigned is False


def test_audit_pilot_genuine_ambiguity_still_unassigned() -> None:
    """The pre-existing 'Round 3 / Cycle 7' case must STILL be unassigned
    — multiple distinct parses are the genuine-ambiguity path that
    ``is_unassigned`` correctly identifies via
    ``len(round_candidates) > 1``."""
    r = parse_round("SRR1", sample_title="Round 3 / Cycle 7")
    assert r.confidence == "MEDIUM"
    assert r.round_candidates == [3, 7]
    assert r.is_unassigned is True


# ----- Phase 6b.5c — round-parser empirical pattern expansion -----
#
# Each pattern below was missing from the Phase 6b.4 cascade and caused a
# real SELEX deposit (per Codex's per-accession fetch_metadata.json
# inspection) to fall into NO-CONFIDENCE territory. Adding the patterns
# moves these to HIGH/MEDIUM-confidence assignments.


def test_pattern_rv01_sample_title_high() -> None:
    """PRJEB51212: sample_title='RV01' → round 1 (HIGH from L2)."""
    r = parse_round("SRR1", sample_title="RV01")
    assert r.round_number == 1
    assert r.confidence == "HIGH"
    assert r.is_unassigned is False


def test_pattern_rv01_library_name_medium() -> None:
    """library_name='RV01_s' (PRJEB51212 actual values like 'RV01_s')
    must still parse to round 1 (MEDIUM from L3)."""
    r = parse_round("SRR1", library_name="RV01_s")
    assert r.round_number == 1
    assert r.confidence == "MEDIUM"
    assert r.source_field == "library_name"
    assert r.is_unassigned is False


def test_pattern_rv_does_not_match_just_R_or_V_alone() -> None:
    """Sanity: pattern requires the glued ``RV`` literal; doesn't fire on
    ``R 01`` (space-separated, already handled by R_digit_boundary) or
    ``V01`` (no R prefix)."""
    # "V01" alone doesn't match RV_digit (needs R-V-digits)
    r = parse_round("SRR1", sample_title="V01_library")
    assert r.round_number is None


def test_pattern_glued_dnafoxr00_sample_title_high() -> None:
    """PRJEB62756: sample_title='DNAFOXR00' → round 0 (HIGH from L2)."""
    r = parse_round("SRR1", sample_title="DNAFOXR00")
    assert r.round_number == 0
    assert r.confidence == "HIGH"
    assert r.is_unassigned is False


def test_pattern_glued_dnafoxr12_library_name_medium() -> None:
    """library_name='DNAFOXR12' → round 12 (MEDIUM from L3)."""
    r = parse_round("SRR1", library_name="DNAFOXR12")
    assert r.round_number == 12
    assert r.confidence == "MEDIUM"
    assert r.is_unassigned is False


def test_pattern_glued_requires_3_letters_before_R() -> None:
    """The glued pattern requires at least 3 letters before R to limit
    false positives from short protein abbreviations like 'TFR1'
    (transcription factor R1, not a SELEX round)."""
    # TF is only 2 letters before R → must NOT match this pattern alone.
    # But "TFR1" still has a fallback: if NO pattern matches, returns None.
    r = parse_round("SRR1", library_name="TFR1")
    # Neither glued_word_R_digit (need ≥3 letters) nor R_digit_boundary
    # (no word boundary before R) should fire on TFR1.
    assert r.round_number is None


def test_pattern_digit_cyc_suffix_selex_7_cyc() -> None:
    """PRJNA385825: sample_title='Selex_7_cyc' → round 7."""
    r = parse_round("SRR1", sample_title="Selex_7_cyc")
    assert r.round_number == 7
    assert r.confidence == "HIGH"
    assert r.is_unassigned is False


@pytest.mark.parametrize(
    "text,expected",
    [
        ("7_cyc", 7),
        ("7 cyc", 7),
        ("7-cyc", 7),
        ("7_cycle", 7),
        ("12 cycles", 12),
        ("Selex_3_cyc", 3),
    ],
)
def test_pattern_digit_cyc_variants(text: str, expected: int) -> None:
    """The digit-before-cyc pattern handles common separator + plural variants."""
    r = parse_round("SRR1", sample_title=text)
    assert r.round_number == expected


def test_pattern_digit_cyc_does_not_match_cycle_word_digit() -> None:
    """``cycle 5`` (cycle BEFORE digit) is the pre-existing
    ``cycle_word_digit`` path; new pattern handles the reverse only."""
    # Both patterns may fire on "cycle 5 cycles 8" but each on its own
    # axis. We just check the simple "cycle 5" still works at HIGH.
    r = parse_round("SRR1", sample_title="cycle 5")
    assert r.round_number == 5
    assert r.confidence == "HIGH"


# ----- L4: abstract count is informative only -----


def test_l4_abstract_count_does_not_assign_round() -> None:
    r = parse_round("SRR1", abstract="We performed 10 rounds of selection.")
    assert r.round_number is None
    assert r.confidence == "NONE"
    assert "10 rounds" in r.parser_notes


def test_extract_round_count_from_abstract_patterns() -> None:
    assert extract_round_count_from_abstract("8 rounds of SELEX") == 8
    assert extract_round_count_from_abstract("rounds 1 through 12") == 12
    assert extract_round_count_from_abstract("R0 to R6") == 6
    assert extract_round_count_from_abstract("performed 9 selection cycles") == 9
    assert extract_round_count_from_abstract("") is None
    assert extract_round_count_from_abstract("no round info here") is None


# ----- L5: manual review dump -----


def test_l5_unknown_writes_manual_review_dump(tmp_path: Path) -> None:
    r = parse_round(
        "SRR999",
        sample_title="opaque title",
        bioproject_id="PRJNA1",
        manual_review_dir=tmp_path,
    )
    assert r.confidence == "NONE"
    assert r.is_unassigned is True  # round_number is None → unassigned
    dump = tmp_path / "PRJNA1_SRR999.txt"
    assert dump.exists()
    assert "MANUAL REVIEW REQUIRED" in dump.read_text()


# ----- Target hint extraction -----


def test_target_hint_extracted_from_capitalized_prefix() -> None:
    r = parse_round("SRR1", sample_title="Thrombin Round 4")
    assert r.target_hint == "Thrombin"


def test_target_hint_none_when_no_match() -> None:
    r = parse_round("SRR1", sample_title="round 1")
    assert r.target_hint is None


# ----- Seed override -----


def test_apply_seed_overrides_replaces_record() -> None:
    base = [
        RoundRecord(
            srr="SRR1",
            round_number=None,
            confidence="NONE",
            source_field="none",
            matched_pattern="none",
        ),
        RoundRecord(
            srr="SRR2",
            round_number=5,
            confidence="HIGH",
            source_field="sample_title",
            matched_pattern="round_word_digit",
        ),
    ]
    out = apply_seed_overrides(base, {"SRR1": 3, "SRR2": 4})
    assert out[0].round_number == 3
    assert out[0].confidence == "HIGH"
    assert out[0].source_field == "seed_override"
    assert out[1].round_number == 4
    assert out[1].source_field == "seed_override"


def test_apply_seed_overrides_passes_unmapped() -> None:
    base = [
        RoundRecord(
            srr="SRR1",
            round_number=5,
            confidence="HIGH",
            source_field="sample_title",
            matched_pattern="round_word_digit",
        )
    ]
    out = apply_seed_overrides(base, {"SRR99": 1})
    assert out[0].round_number == 5
    assert out[0].source_field == "sample_title"


# ----- Summary -----


def test_summarise_bioproject_flags_inconsistency() -> None:
    records = [
        RoundRecord(
            srr="SRR1",
            round_number=1,
            confidence="HIGH",
            source_field="sample_title",
            matched_pattern="round_word_digit",
        ),
        RoundRecord(
            srr="SRR2",
            round_number=None,
            confidence="NONE",
            source_field="none",
            matched_pattern="none",
        ),
    ]
    summary = summarise_bioproject("PRJ1", records, abstract="6 rounds of SELEX")
    assert summary.n_high == 1
    assert summary.n_none == 1
    assert summary.inconsistent_annotation is True
    assert summary.n_rounds_from_abstract == 6


def test_summarise_bioproject_clean_when_all_assigned() -> None:
    records = [
        RoundRecord(
            srr="SRR1",
            round_number=1,
            confidence="HIGH",
            source_field="sample_title",
            matched_pattern="round_word_digit",
        ),
        RoundRecord(
            srr="SRR2",
            round_number=2,
            confidence="HIGH",
            source_field="sample_title",
            matched_pattern="round_word_digit",
        ),
    ]
    summary = summarise_bioproject("PRJ1", records)
    assert summary.inconsistent_annotation is False
