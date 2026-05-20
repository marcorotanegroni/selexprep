"""Smoke tests for ``selexprep.qc.plots``.

Verify each plot writes a non-empty PNG. Visual correctness is left to
manual inspection (PNGs aren't byte-deterministic across matplotlib
versions and are documented as informational-only outputs).
"""

from __future__ import annotations

from pathlib import Path

from selexprep.library.report import LibraryReport
from selexprep.manifest import build_manifest_from_extract_result
from selexprep.qc.plots import (
    plot_n_length_distribution,
    plot_per_round_panel,
    plot_primer_match_per_round,
    plot_read_retention,
)


def _make_library_report() -> LibraryReport:
    return LibraryReport(
        primer_5p="ACGTACGTACGTACGT",
        primer_3p="TGCATGCATGCATGCA",
        variants_5p=[],
        variants_3p=[],
        known_adapter_hits={},
        extraction_mode="BOTH_PRIMERS_SINGLE_READ",
        full_insert_recovered=True,
        read_source="R1",
        required_action="NONE",
        orientation="FORWARD",
        n_length_mode=30,
        n_length_distribution={30: 100},
        n_length_confidence=1.0,
        match_rate_5p=0.95,
        match_rate_3p=0.92,
        position_consistency_5p=0.95,
        position_consistency_3p=0.92,
        read_fraction_used_for_inference=1.0,
        sampling_seed=42,
        confidence=0.85,
        status="HIGH",
        failure_reason=None,
    )


def _counts_by_round(n_rounds: int = 3) -> dict[int, dict[str, int]]:
    out: dict[int, dict[str, int]] = {}
    for r in range(n_rounds):
        # Diversity drops across rounds (typical SELEX enrichment pattern).
        n_unique = 1000 // (r + 1)
        out[r] = {f"r{r}_seq_{i}": (i + 1) * 5 for i in range(n_unique)}
    return out


def _trim_reports_by_round(n_rounds: int = 3) -> dict[int, dict[str, int]]:
    return {r: {"n_in": 100_000 // (r + 1), "n_out": 95_000 // (r + 1)} for r in range(n_rounds)}


# ---------------------------------------------------------------------------
# Each plot writes a non-empty PNG
# ---------------------------------------------------------------------------


def test_plot_read_retention_smoke(tmp_path: Path) -> None:
    path = plot_read_retention(_trim_reports_by_round(), tmp_path)
    assert path.exists()
    assert path.stat().st_size > 1000  # PNG with bars should be > 1 KB
    assert path.suffix == ".png"


def test_plot_primer_match_per_round_smoke(tmp_path: Path) -> None:
    manifest = build_manifest_from_extract_result(
        library_report=_make_library_report(),
        input_paths=[],
        output_paths=[],
        accession=None,
        bioproject_id=None,
        runs=[],
        parameters={},
    )
    path = plot_primer_match_per_round(manifest, _trim_reports_by_round(), tmp_path)
    assert path.exists()
    assert path.stat().st_size > 1000


def test_plot_n_length_distribution_smoke(tmp_path: Path) -> None:
    path = plot_n_length_distribution(_counts_by_round(), tmp_path)
    assert path.exists()
    assert path.stat().st_size > 1000


def test_plot_per_round_panel_smoke(tmp_path: Path) -> None:
    path = plot_per_round_panel(_counts_by_round(), tmp_path)
    assert path.exists()
    assert path.stat().st_size > 1000


def test_plot_n_length_distribution_handles_empty_input(tmp_path: Path) -> None:
    """Empty counts -> placeholder figure, not an exception."""
    path = plot_n_length_distribution({}, tmp_path)
    assert path.exists()
