"""Tests for ``selexprep.benchmark.figure_b`` — the Tier 2 audit table emitter.

``figure_b.py`` now emits a Markdown summary table of the corpus-audit distributions
(``audit_metrics.json``), not a bar chart.
"""

from __future__ import annotations

import json
from pathlib import Path

from selexprep.benchmark.figure_b import _build_title, emit_audit_table


def _make_audit_payload() -> dict:
    return {
        "catalog_version": "v0.1.5-test",
        "sample_seed": 42,
        "sample_accessions_sha256": "deadbeef" * 8,
        "n_sampled": 30,
        "n_in_ground_truth_overlap": 0,
        "fetch_outcome_distribution": {"OK": 20, "FETCH_FAILED": 5, "FETCH_REFUSED": 5},
        "n_fetchable": 20,
        "n_with_library_report": 20,
        "library_report_status_distribution": {
            "HIGH": 12,
            "MEDIUM": 4,
            "LOW": 2,
            "UNABLE_TO_INFER": 2,
        },
        "extraction_mode_distribution": {
            "BOTH_PRIMERS_SINGLE_READ": 14,
            "UNABLE_TO_EXTRACT": 6,
        },
        "required_action_distribution": {
            "NONE": 14,
            "MANUAL_PRIMERS_REQUIRED": 6,
        },
        "inference_safe_failure_rate": 0.3,
        "n_inference_safe_failures": 6,
        "n_with_qc_run": 14,
        "flags_raised_histogram": {"0": 10, "1": 4},
        "per_accession": [],
    }


def test_emit_audit_table_writes_md(tmp_path: Path) -> None:
    audit = tmp_path / "audit.json"
    audit.write_text(json.dumps(_make_audit_payload()), encoding="utf-8")
    out = emit_audit_table(audit, tmp_path / "out")
    assert out.name == "table_audit.md"
    text = out.read_text(encoding="utf-8")
    assert "public-corpus audit" in text
    assert "Fetch outcomes" in text
    assert "| OK | 20 |" in text
    assert "| HIGH | 12 |" in text
    assert "Inference safe-failure rate:** 30% (6/20)" in text


def test_emit_audit_table_handles_empty(tmp_path: Path) -> None:
    """No data → still emits a table (with (no data) rows), no crash."""
    audit = tmp_path / "empty.json"
    audit.write_text(json.dumps({"n_sampled": 0}), encoding="utf-8")
    out = emit_audit_table(audit, tmp_path / "out")
    text = out.read_text(encoding="utf-8")
    assert "public-corpus audit" in text
    assert "(no data)" in text


def test_emit_audit_table_handles_unexpected_status_label(tmp_path: Path) -> None:
    """Forward-compatibility: an unanticipated category still renders (extras appended)."""
    audit = tmp_path / "audit.json"
    payload = _make_audit_payload()
    payload["fetch_outcome_distribution"]["NEW_FUTURE_STATUS"] = 3
    audit.write_text(json.dumps(payload), encoding="utf-8")
    text = emit_audit_table(audit, tmp_path / "out").read_text(encoding="utf-8")
    assert "| NEW_FUTURE_STATUS | 3 |" in text


# ---------------------------------------------------------------------------
# caption (_build_title) segments
# ---------------------------------------------------------------------------


def test_build_title_omits_eligibility_segment_when_classifier_did_not_run() -> None:
    """Pre-audit JSON (n_catalog_classified=0) → caption has no layer-1 segment."""
    title = _build_title(_make_audit_payload())
    assert "audit-eligible" not in title
    assert "selexprep public-corpus audit" in title


def test_build_title_includes_insdc_only_eligibility_segment_without_catalog() -> None:
    """Eligibility set but catalog total absent → INSDC-only segment."""
    payload = _make_audit_payload()
    payload["n_catalog_classified"] = 95
    payload["n_catalog_eligible"] = 24
    title = _build_title(payload)
    assert "24 of 95 INSDC rows audit-eligible" in title
    assert "non-INSDC passthrough" not in title
    assert "catalog total" not in title


def test_build_title_includes_full_catalog_denominator_when_present() -> None:
    """--catalog populates n_catalog_total + non-INSDC count → full segment."""
    payload = _make_audit_payload()
    payload["n_catalog_classified"] = 95
    payload["n_catalog_eligible"] = 24
    payload["n_catalog_total"] = 220
    payload["n_catalog_non_insdc_passthrough"] = 125
    title = _build_title(payload)
    assert "24 of 95 INSDC rows audit-eligible" in title
    assert "125 non-INSDC passthrough" in title
    assert "220 catalog total" in title
