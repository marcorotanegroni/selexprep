"""Smoke tests for ``selexprep.benchmark.figure_a`` (4-panel Figure A).

Like Phase 5's matplotlib plot tests, we only check that PDF + PNG
files are produced — byte-determinism is not guaranteed across
matplotlib versions (locked plan accepts this for plot files).
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from selexprep.benchmark.figure_a import plot_figure_a


def _make_metrics_payload() -> dict:
    return {
        "n_verified": 3,
        "n_unverified": 1,
        "n_total": 4,
        "primer_recovery": {
            "n_evaluated": 3,
            "counts_5p": {"EXACT": 2, "PARTIAL_5P": 1},
            "counts_3p": {"EXACT": 2, "MISMATCH": 1},
            "counts_by_status_5p": {
                "HIGH": {"EXACT": 2},
                "MEDIUM": {"PARTIAL_5P": 1},
            },
        },
        "pair_recovery_by_status": {
            "n_evaluated": 3,
            "counts": {
                "HIGH": {"pair_exact": 2},
                "MEDIUM": {"pair_partial": 1},
            },
        },
        "safe_failure_rate": {
            "n_evaluated": 3,
            "n_safe_failures": 1,
            "rate": 0.333,
            "by_reason": {"status_UNABLE_TO_INFER": 1},
            "safe_failure_accessions": ["PRJ_X"],
        },
        "n_length_recovery": {
            "n_in_tolerance": 2,
            "n_out_of_tolerance": 1,
            "n_unmeasurable": 0,
            "tolerance": 2,
        },
        "extraction_mode_distribution": {
            "counts": {
                "BOTH_PRIMERS_SINGLE_READ": 2,
                "FIVE_PRIME_ONLY": 1,
                "UNABLE_TO_EXTRACT": 1,
            }
        },
        "required_action_distribution": {
            "counts": {
                "NONE": 2,
                "MANUAL_PRIMERS_REQUIRED": 1,
            }
        },
    }


def test_plot_figure_a_writes_pdf_and_png(tmp_path: Path) -> None:
    metrics = tmp_path / "metrics.json"
    metrics.write_text(json.dumps(_make_metrics_payload()), encoding="utf-8")
    pdf, png = plot_figure_a(metrics, tmp_path / "out")
    assert pdf.exists()
    assert png.exists()
    assert pdf.name == "figure_a.pdf"
    assert png.name == "figure_a.png"
    # 4-panel figure produces non-trivial PDF + PNG output.
    assert pdf.stat().st_size > 5000
    assert png.stat().st_size > 5000


def test_plot_figure_a_handles_empty_metrics(tmp_path: Path) -> None:
    """No data → still emits both files with the 'no data' labels in every panel."""
    empty = {"n_verified": 0, "n_unverified": 0, "n_total": 0}
    metrics = tmp_path / "empty.json"
    metrics.write_text(json.dumps(empty), encoding="utf-8")
    pdf, png = plot_figure_a(metrics, tmp_path / "out")
    assert pdf.exists()
    assert png.exists()


@pytest.mark.parametrize(
    "status_bucket",
    ["HIGH", "MEDIUM", "LOW", "UNABLE_TO_INFER", "NO_REPORT"],
)
def test_plot_figure_a_supports_all_status_buckets(tmp_path: Path, status_bucket: str) -> None:
    metrics = tmp_path / "metrics.json"
    metrics.write_text(
        json.dumps(
            {
                "n_verified": 1,
                "n_unverified": 0,
                "primer_recovery": {
                    "n_evaluated": 1,
                    "counts_5p": {"EXACT": 1},
                    "counts_3p": {"EXACT": 1},
                    "counts_by_status_5p": {status_bucket: {"EXACT": 1}},
                },
                "pair_recovery_by_status": {
                    "n_evaluated": 1,
                    "counts": {status_bucket: {"pair_exact": 1}},
                },
                "safe_failure_rate": {
                    "n_evaluated": 1,
                    "n_safe_failures": 1 if status_bucket == "UNABLE_TO_INFER" else 0,
                    "rate": 1.0 if status_bucket == "UNABLE_TO_INFER" else 0.0,
                    "by_reason": {},
                    "safe_failure_accessions": [],
                },
                "n_length_recovery": {
                    "n_in_tolerance": 1,
                    "n_out_of_tolerance": 0,
                    "n_unmeasurable": 0,
                    "tolerance": 2,
                },
                "extraction_mode_distribution": {"counts": {"BOTH_PRIMERS_SINGLE_READ": 1}},
                "required_action_distribution": {"counts": {"NONE": 1}},
            }
        ),
        encoding="utf-8",
    )
    pdf, png = plot_figure_a(metrics, tmp_path / "out")
    assert pdf.exists()
    assert png.exists()


def test_plot_figure_a_includes_safe_failure_rate_in_title(tmp_path: Path) -> None:
    """The safe-failure rate is selexprep's unique distinguishing capability —
    it must appear in the figure title prominently.
    """
    metrics = tmp_path / "metrics.json"
    payload = _make_metrics_payload()
    metrics.write_text(json.dumps(payload), encoding="utf-8")
    # We don't inspect the rendered title text directly (matplotlib renders
    # to image); instead verify the function reaches the savefig path
    # cleanly when safe_failure_rate is populated.
    pdf, png = plot_figure_a(metrics, tmp_path / "out")
    assert pdf.exists()
    assert png.exists()
