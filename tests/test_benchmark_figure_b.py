"""Smoke tests for ``selexprep.benchmark.figure_b`` (4-panel Figure B).

Mirrors ``tests/test_benchmark_figure_a.py``. Like Phase 5's matplotlib
plot tests, we only check that PDF + PNG files are produced — byte
determinism is not guaranteed across matplotlib versions (locked plan
accepts this for plot files; audit_metrics.json IS the deterministic
source of truth).
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from selexprep.benchmark.figure_b import plot_figure_b


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


def test_plot_figure_b_writes_pdf_and_png(tmp_path: Path) -> None:
    audit = tmp_path / "audit.json"
    audit.write_text(json.dumps(_make_audit_payload()), encoding="utf-8")
    pdf, png = plot_figure_b(audit, tmp_path / "out")
    assert pdf.exists()
    assert png.exists()
    assert pdf.name == "figure_b.pdf"
    assert png.name == "figure_b.png"
    assert pdf.stat().st_size > 5000
    assert png.stat().st_size > 5000


def test_plot_figure_b_handles_empty_audit(tmp_path: Path) -> None:
    """No data → still emits both files with 'no data' labels in every panel."""
    audit = tmp_path / "empty.json"
    audit.write_text(json.dumps({"n_sampled": 0}), encoding="utf-8")
    pdf, png = plot_figure_b(audit, tmp_path / "out")
    assert pdf.exists()
    assert png.exists()


@pytest.mark.parametrize(
    "status_bucket",
    ["HIGH", "MEDIUM", "LOW", "UNABLE_TO_INFER"],
)
def test_plot_figure_b_supports_all_status_buckets(tmp_path: Path, status_bucket: str) -> None:
    """Each LibraryReport.status bucket appears in Panel B without errors."""
    audit = tmp_path / "audit.json"
    payload = _make_audit_payload()
    payload["library_report_status_distribution"] = {status_bucket: 5}
    payload["n_with_library_report"] = 5
    audit.write_text(json.dumps(payload), encoding="utf-8")
    pdf, png = plot_figure_b(audit, tmp_path / "out")
    assert pdf.exists()
    assert png.exists()


def test_plot_figure_b_renders_with_safe_failure_rate(tmp_path: Path) -> None:
    """The inference safe-failure overlay is the unique distinguishing metric.

    We don't assert text content in the rendered raster (matplotlib's
    PNG output is environment-dependent), but the render path must
    survive when the rate is populated — the overlay annotation block
    is on the same matplotlib axes as Panel D.
    """
    audit = tmp_path / "audit.json"
    payload = _make_audit_payload()
    payload["inference_safe_failure_rate"] = 0.42
    payload["n_inference_safe_failures"] = 4
    audit.write_text(json.dumps(payload), encoding="utf-8")
    pdf, png = plot_figure_b(audit, tmp_path / "out")
    assert pdf.exists()
    assert png.exists()


def test_plot_figure_b_handles_unexpected_status_label(tmp_path: Path) -> None:
    """Forward-compatibility: a status the canonical order doesn't anticipate
    still renders (sorted alphabetically among the extras)."""
    audit = tmp_path / "audit.json"
    payload = _make_audit_payload()
    payload["fetch_outcome_distribution"]["NEW_FUTURE_STATUS"] = 3
    audit.write_text(json.dumps(payload), encoding="utf-8")
    pdf, png = plot_figure_b(audit, tmp_path / "out")
    assert pdf.exists()
    assert png.exists()
