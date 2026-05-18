"""Unit tests for selexprep.qc.coverage."""

from __future__ import annotations

from selexprep.qc.coverage import (
    STATUS_ALL,
    STATUS_MULTIPLEXED,
    STATUS_PARTIAL,
    STATUS_UNKNOWN,
    build_round_coverage_report,
    classify_bioproject_round_coverage,
)


def _bp(bp_id: str, **kw) -> dict:
    base = {"bioproject_id": bp_id, "include": "y"}
    base.update(kw)
    return base


def _round(srr: str, round_number: int | None, confidence: str = "HIGH") -> dict:
    return {
        "srr": srr,
        "round_number": str(round_number) if round_number is not None else "",
        "confidence": confidence,
    }


# ----- classify single BP -----


def test_classify_all_rounds_when_observed_equals_declared() -> None:
    bp = _bp("PRJ1", n_rounds_declared="3")
    samples = [
        {"srr": "SRR1", "bioproject_id": "PRJ1"},
        {"srr": "SRR2", "bioproject_id": "PRJ1"},
        {"srr": "SRR3", "bioproject_id": "PRJ1"},
    ]
    rounds_by_srr = {
        "SRR1": _round("SRR1", 1),
        "SRR2": _round("SRR2", 2),
        "SRR3": _round("SRR3", 3),
    }
    out = classify_bioproject_round_coverage(bp, samples, rounds_by_srr)
    assert out["round_coverage_status"] == STATUS_ALL
    assert out["n_rounds_assigned"] == 3
    assert out["assigned_rounds"] == [1, 2, 3]


def test_classify_partial_when_observed_fewer_than_declared() -> None:
    bp = _bp("PRJ1", n_rounds_declared="5")
    samples = [{"srr": "SRR1"}, {"srr": "SRR2"}]
    rounds_by_srr = {"SRR1": _round("SRR1", 1), "SRR2": _round("SRR2", 2)}
    out = classify_bioproject_round_coverage(bp, samples, rounds_by_srr)
    assert out["round_coverage_status"] == STATUS_PARTIAL


def test_classify_multiplexed_when_hint_present() -> None:
    bp = _bp(
        "PRJ1", n_rounds_declared="6", manual_curation_notes="A single SRR contains rounds 1-6"
    )
    samples = [{"srr": "SRR1"}]
    rounds_by_srr = {"SRR1": _round("SRR1", 1)}
    out = classify_bioproject_round_coverage(bp, samples, rounds_by_srr)
    assert out["round_coverage_status"] == STATUS_MULTIPLEXED
    assert out["multiplex_hint"] is True


def test_classify_unknown_no_public_srrs() -> None:
    bp = _bp("PRJ1", n_rounds_declared="3")
    out = classify_bioproject_round_coverage(bp, samples=[], rounds_by_srr={})
    assert out["round_coverage_status"] == STATUS_UNKNOWN
    assert out["reason"] == "no_public_srrs"


def test_classify_unknown_declared_round_count_missing() -> None:
    bp = _bp("PRJ1")  # no n_rounds_declared
    samples = [{"srr": "SRR1"}]
    rounds_by_srr = {"SRR1": _round("SRR1", 1)}
    out = classify_bioproject_round_coverage(bp, samples, rounds_by_srr)
    assert out["round_coverage_status"] == STATUS_UNKNOWN
    assert out["reason"] == "declared_round_count_missing"


def test_classify_unknown_when_unresolved_srrs_remain() -> None:
    bp = _bp("PRJ1", n_rounds_declared="5")
    samples = [{"srr": "SRR1"}, {"srr": "SRR_unknown"}]
    rounds_by_srr = {"SRR1": _round("SRR1", 1)}
    out = classify_bioproject_round_coverage(bp, samples, rounds_by_srr)
    assert out["round_coverage_status"] == STATUS_UNKNOWN
    assert "SRR_unknown" in out["unknown_round_srrs"]


# ----- aggregate report -----


def test_build_report_default_filter_uses_include_y() -> None:
    bioprojects = [
        _bp("PRJ1", n_rounds_declared="1", include="y"),
        _bp("PRJ2", n_rounds_declared="1", include="n"),
    ]
    samples = [{"srr": "SRR1", "bioproject_id": "PRJ1"}, {"srr": "SRR2", "bioproject_id": "PRJ2"}]
    rounds = [_round("SRR1", 1), _round("SRR2", 1)]

    report = build_round_coverage_report(bioprojects, samples, rounds)
    assert report["summary"]["n_bioprojects"] == 2
    assert len(report["included_bioprojects"]) == 1
    assert report["included_bioprojects"][0]["bioproject_id"] == "PRJ1"


def test_build_report_custom_filter() -> None:
    """Caller-supplied filter overrides the default include=y rule."""
    bioprojects = [
        _bp("PRJ1", n_rounds_declared="1", include="y", library_type_verification="RNA_confirmed"),
        _bp("PRJ2", n_rounds_declared="1", include="y", library_type_verification="DNA_confirmed"),
    ]
    samples = [{"srr": "SRR1", "bioproject_id": "PRJ1"}, {"srr": "SRR2", "bioproject_id": "PRJ2"}]
    rounds = [_round("SRR1", 1), _round("SRR2", 1)]

    def rna_only(bp: dict) -> bool:
        return bp.get("library_type_verification") == "RNA_confirmed"

    report = build_round_coverage_report(bioprojects, samples, rounds, include_filter=rna_only)
    assert len(report["included_bioprojects"]) == 1
    assert report["included_bioprojects"][0]["bioproject_id"] == "PRJ1"
