"""Figure B — Tier 2 public-corpus audit for the Application Note.

Mirrors :mod:`selexprep.benchmark.figure_a`'s 2x2 panel layout but
answers a different question: across a sampled subset of the public
HT-SELEX corpus, what does selexprep say? No per-row ground truth —
distributional metrics only.

Methodological correction from the Codex peer-review + user pass:
EACH panel labels its denominator in the panel subtitle so a reviewer
never has to guess the normalization. Two distinct denominators live
in the same figure:

- ``n_sampled`` — every sampled INSDC accession, used by Panel A
  (fetch outcomes). Includes rows that never produced a LibraryReport.
- ``n_with_library_report`` — only the rows where ``detect`` actually
  ran and emitted a report, used by Panels B/C/D. Mixing fetch failures
  with inference failures would inflate the inference safe-failure
  metric with ENA/network problems.

The figure title carries:

- ``n_sampled``
- ``n_fetchable`` (rows past the fetch stage)
- ``n_with_library_report`` (rows that produced a LibraryReport)
- catalog version
- sample seed

so a reviewer can read the headline numbers without digging into the
JSON. The audit JSON is the source of truth — PNG byte-output is not
guaranteed deterministic across matplotlib versions (locked plan
accepts this for plot files; same discipline as Figure A / Phase 5).
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)


# Stable, ordered buckets so the visual order doesn't reshuffle between
# runs. Mirrors :mod:`selexprep.benchmark.figure_a` conventions.
_STATUS_ORDER: tuple[str, ...] = ("HIGH", "MEDIUM", "LOW", "UNABLE_TO_INFER")
_EXTRACTION_MODE_ORDER: tuple[str, ...] = (
    "BOTH_PRIMERS_SINGLE_READ",
    "FIVE_PRIME_ONLY",
    "THREE_PRIME_ONLY",
    "PAIRED_END_SPLIT_PRIMERS",
    "UNABLE_TO_EXTRACT",
)
_REQUIRED_ACTION_ORDER: tuple[str, ...] = (
    "NONE",
    "MANUAL_PRIMERS_REQUIRED",
    "READ_MERGING_RECOMMENDED",
)
# Run-status order is loose; we sort by descending count at render time
# and only fix the order of "OK" and the well-known failure modes when
# they appear so the bar layout stays readable.
_FETCH_OUTCOME_ORDER: tuple[str, ...] = (
    "OK",
    "SKIPPED_READ_MERGING_RECOMMENDED",
    "DETECT_FAILED",
    "EXTRACT_REFUSED",
    "EXTRACT_FAILED",
    "COUNT_FAILED",
    "QC_FAILED",
    "FETCH_REFUSED",
    "FETCH_FAILED",
    "UNEXPECTED_FAILURE",
)


def _no_data_label(ax: Any, title: str) -> None:
    ax.set_title(f"{title} (no data)")
    ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)


def _ordered_categories(
    counts: dict[str, int], canonical_order: tuple[str, ...]
) -> tuple[list[str], list[int]]:
    """Return (labels, values) following canonical order, then any extras alphabetically.

    Extras can appear if the runner ever emits a status the canonical
    order doesn't anticipate (a forward-compatibility safety net).
    """
    in_order = [k for k in canonical_order if k in counts]
    extras = sorted(k for k in counts if k not in set(canonical_order))
    labels = in_order + extras
    values = [counts[k] for k in labels]
    return labels, values


def _panel_a_fetch_outcomes(ax: Any, audit: dict[str, Any]) -> None:
    """Panel A · Fetch outcomes · denominator = n_sampled."""
    counts: dict[str, int] = audit.get("fetch_outcome_distribution", {})
    n_sampled = int(audit.get("n_sampled", 0))
    if not counts:
        _no_data_label(ax, "A · Fetch outcomes")
        return
    labels, values = _ordered_categories(counts, _FETCH_OUTCOME_ORDER)
    ax.bar(labels, values)
    ax.set_title(f"A · Fetch outcomes\nN={n_sampled} sampled INSDC accessions")
    ax.set_ylabel("rows")
    ax.tick_params(axis="x", labelrotation=30)


def _panel_b_inference_confidence(ax: Any, audit: dict[str, Any]) -> None:
    """Panel B · Inference confidence · denominator = n_with_library_report."""
    counts: dict[str, int] = audit.get("library_report_status_distribution", {})
    n_lr = int(audit.get("n_with_library_report", 0))
    if not counts:
        _no_data_label(ax, "B · Inference confidence")
        return
    labels, values = _ordered_categories(counts, _STATUS_ORDER)
    ax.bar(labels, values)
    ax.set_title(
        "B · Inference confidence (LibraryReport.status)\n"
        f"N={n_lr} rows reaching LibraryReport (denominator excludes fetch failures)"
    )
    ax.set_ylabel("rows")


def _panel_c_extraction_mode(ax: Any, audit: dict[str, Any]) -> None:
    """Panel C · Extraction mode · denominator = n_with_library_report."""
    counts: dict[str, int] = audit.get("extraction_mode_distribution", {})
    n_lr = int(audit.get("n_with_library_report", 0))
    if not counts:
        _no_data_label(ax, "C · Extraction mode")
        return
    labels, values = _ordered_categories(counts, _EXTRACTION_MODE_ORDER)
    ax.barh(labels, values)
    ax.invert_yaxis()
    ax.set_title(f"C · Extraction mode (honest accounting)\nN={n_lr} rows reaching LibraryReport")
    ax.set_xlabel("rows")


def _panel_d_required_action(ax: Any, audit: dict[str, Any]) -> None:
    """Panel D · Required action + safe-failure overlay · denominator = n_with_library_report."""
    counts: dict[str, int] = audit.get("required_action_distribution", {})
    n_lr = int(audit.get("n_with_library_report", 0))
    rate = float(audit.get("inference_safe_failure_rate", 0.0))
    n_safe = int(audit.get("n_inference_safe_failures", 0))
    if not counts:
        _no_data_label(ax, "D · Required action")
        return
    labels, values = _ordered_categories(counts, _REQUIRED_ACTION_ORDER)
    ax.barh(labels, values)
    ax.invert_yaxis()
    ax.set_title(f"D · Required action (workflow guidance)\nN={n_lr} rows reaching LibraryReport")
    ax.set_xlabel("rows")
    # Inference safe-failure rate as overlay annotation — the unique
    # distinguishing metric vs known-primer pipelines, prominent in
    # the panel that's most informative about it.
    annotation = f"inference safe-failure rate: {rate:.0%} ({n_safe}/{n_lr})"
    ax.text(
        0.98,
        0.05,
        annotation,
        ha="right",
        va="bottom",
        transform=ax.transAxes,
        fontsize="small",
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "gray"},
    )


def plot_figure_b(audit_json: Path, outdir: Path) -> tuple[Path, Path]:
    """Render the four-panel Figure B.

    Parameters
    ----------
    audit_json
        Path to an ``audit_metrics.json`` produced by
        :func:`selexprep.benchmark.corpus_audit.write_audit_json`.
    outdir
        Directory to write ``figure_b.{pdf,png}`` into.

    Returns
    -------
    A tuple ``(pdf_path, png_path)``.
    """
    audit = json.loads(audit_json.read_text(encoding="utf-8"))
    outdir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(2, 2, figsize=(13, 11), dpi=150)
    _panel_a_fetch_outcomes(axes[0, 0], audit)
    _panel_b_inference_confidence(axes[0, 1], audit)
    _panel_c_extraction_mode(axes[1, 0], audit)
    _panel_d_required_action(axes[1, 1], audit)

    n_sampled = int(audit.get("n_sampled", 0))
    n_fetchable = int(audit.get("n_fetchable", 0))
    n_lr = int(audit.get("n_with_library_report", 0))
    catalog_version = audit.get("catalog_version") or "unspecified"
    sample_seed = audit.get("sample_seed", 42)
    # Phase 6b.5b layer-1: eligibility-stage counts. Empty when the
    # audit was generated by a pre-6b.5b aggregator (backward compat).
    n_catalog_classified = int(audit.get("n_catalog_classified", 0))
    n_catalog_eligible = int(audit.get("n_catalog_eligible", 0))

    title_parts = [
        "selexprep Figure B — public-corpus audit",
    ]
    if n_catalog_classified > 0:
        title_parts.append(
            f"{n_catalog_eligible} of {n_catalog_classified} catalog rows audit-eligible"
        )
    title_parts.append(
        f"N={n_sampled} sampled / {n_fetchable} fetchable / {n_lr} with LibraryReport"
    )
    title_parts.append(f"catalog {catalog_version}, seed {sample_seed}")
    fig.suptitle("  ·  ".join(title_parts), fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))

    pdf_path = outdir / "figure_b.pdf"
    png_path = outdir / "figure_b.png"
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, bbox_inches="tight")
    plt.close(fig)
    return pdf_path, png_path


# ---------------------------------------------------------------------------
# CLI entry (called by the audit Snakefile's figure_b rule)
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Render the Phase 6b.3a Figure B.")
    p.add_argument("--audit", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    args = p.parse_args(argv)
    pdf, png = plot_figure_b(args.audit, args.outdir)
    print(f"wrote {pdf}")
    print(f"wrote {png}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())


__all__ = ["main", "plot_figure_b"]
