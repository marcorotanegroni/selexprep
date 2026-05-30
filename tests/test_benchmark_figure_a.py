"""Smoke tests for ``selexprep.benchmark.figure_a`` (two-arm Figure A).

Like Phase 5's matplotlib plot tests, we only check that PDF + PNG files
are produced — byte-determinism is not guaranteed across matplotlib
versions (locked plan accepts this for plot files). Phase 6b.10 reframes
the figure as a recovery / specificity two-arm benchmark; these tests
pin the new metrics.json shape.
"""

from __future__ import annotations

import json
from pathlib import Path

from selexprep.benchmark.figure_a import plot_figure_a


def _make_metrics_payload() -> dict:
    """A populated two-arm metrics.json (5 recovery + 3 specificity shape)."""
    return {
        "recovery_denominator": 3,
        "pair_recovery_by_status": {
            "n_evaluated": 3,
            "counts": {
                "HIGH": {"pair_exact": 2},
                "MEDIUM": {"pair_partial": 1},
            },
        },
        "multi_round_sensitivity": {
            "n_evaluated": 2,
            "counts": {"HIGH": {"pair_exact": 2}},
        },
        "primer_recovery": {
            "n_evaluated": 3,
            "counts_5p": {"EXACT": 2, "PARTIAL_5P": 1},
            "counts_3p": {"EXACT": 2, "MISMATCH": 1},
        },
        "specificity": {
            "n_evaluated": 3,
            "n_no_false_call": 3,
            "n_false_positive": 0,
            "false_positive_accessions": [],
            "per_row": [
                {"accession": a, "bucket": "no_false_call", "primer_5p": None, "primer_3p": None}
                for a in ("PRJEB28411", "PRJEB22637", "PRJNA990511")
            ],
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
                "UNABLE_TO_EXTRACT": 3,
            }
        },
        "required_action_distribution": {
            "counts": {
                "NONE": 2,
                "MANUAL_PRIMERS_REQUIRED": 3,
            }
        },
        "fetch_stats": {
            "PRJEB70964": {
                "accession": "PRJEB70964",
                "fetch_expected_runs": 27,
                "fetch_available_runs": 17,
                "fetch_missing_runs": 10,
                "runs_with_no_fastq_url": 0,
                "partial_fetch": True,
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


def test_plot_figure_a_renders_specificity_false_positive(tmp_path: Path) -> None:
    """The specificity arm renders a false-positive call (red bar path)."""
    payload = _make_metrics_payload()
    payload["specificity"] = {
        "n_evaluated": 3,
        "n_no_false_call": 2,
        "n_false_positive": 1,
        "false_positive_accessions": ["PRJEB22637"],
        "per_row": [
            {
                "accession": "PRJEB28411",
                "bucket": "no_false_call",
                "primer_5p": None,
                "primer_3p": None,
            },
            {
                "accession": "PRJEB22637",
                "bucket": "false_positive",
                "primer_5p": "ACGT",
                "primer_3p": None,
            },
            {
                "accession": "PRJNA990511",
                "bucket": "no_false_call",
                "primer_5p": None,
                "primer_3p": None,
            },
        ],
    }
    metrics = tmp_path / "metrics.json"
    metrics.write_text(json.dumps(payload), encoding="utf-8")
    pdf, png = plot_figure_a(metrics, tmp_path / "out")
    assert pdf.exists()
    assert png.exists()


def test_plot_figure_a_renders_without_fetch_stats(tmp_path: Path) -> None:
    """The partial-fetch note path is optional — absent fetch_stats still renders."""
    payload = _make_metrics_payload()
    payload.pop("fetch_stats")
    metrics = tmp_path / "metrics.json"
    metrics.write_text(json.dumps(payload), encoding="utf-8")
    pdf, png = plot_figure_a(metrics, tmp_path / "out")
    assert pdf.exists()
    assert png.exists()


def test_plot_figure_a_recovery_only_no_specificity(tmp_path: Path) -> None:
    """A recovery-only metrics set (specificity arm empty) still renders all panels."""
    payload = _make_metrics_payload()
    payload["specificity"] = {
        "n_evaluated": 0,
        "n_no_false_call": 0,
        "n_false_positive": 0,
        "false_positive_accessions": [],
        "per_row": [],
    }
    metrics = tmp_path / "metrics.json"
    metrics.write_text(json.dumps(payload), encoding="utf-8")
    pdf, png = plot_figure_a(metrics, tmp_path / "out")
    assert pdf.exists()
    assert png.exists()
