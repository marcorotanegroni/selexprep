"""Figure A — primer-inference benchmark for the Application Note.

Scope pivot (2026-05-22): the figure tests selexprep's **unique claim**
(primer inference from accession reads + safe failure on ambiguity),
NOT comparator-tool count agreement. AptaPLEX / EasyDIVER+ both require
known primers as input; they cannot benchmark primer inference. Count
correlation moves to Phase 6c via a self-consistency check (inferred
vs ``--override-primer``).

Four panels arranged 2x2 (locked plan line 365, revised):

- **Panel A (top-left)** — Pair recovery by LibraryReport.status.
  The headline panel: of HIGH-confidence calls, what fraction
  recovered the primer pair exactly? Stacked bar of
  {pair_exact, pair_equivalent, pair_partial, pair_failed} per
  status bucket.

- **Panel B (top-right)** — 5'/3'/pair/partial recovery breakdown
  across all verified rows. Grouped bar showing per-side EXACT
  recovery, pair recovery, partial recovery, and pair_failed.

- **Panel C (bottom-left)** — N-length recovery within ±tolerance
  (locked plan line 365).

- **Panel D (bottom-right)** — Two-row honest accounting:
  ``extraction_mode`` distribution (top) +
  ``required_action`` distribution (bottom).

The figure title carries the headline number: N verified rows +
safe-failure rate.

PNG byte-output is non-deterministic across matplotlib versions
(accepted in Phase 5); the underlying ``metrics.json`` IS deterministic
and is the source of truth for downstream consumers.
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
# runs.
_STATUS_ORDER: tuple[str, ...] = ("HIGH", "MEDIUM", "LOW", "UNABLE_TO_INFER", "NO_REPORT")
_PAIR_BUCKETS: tuple[str, ...] = (
    "pair_exact",
    "pair_equivalent",
    "pair_partial",
    "pair_failed",
)
_EXTRACTION_MODE_ORDER: tuple[str, ...] = (
    "BOTH_PRIMERS_SINGLE_READ",
    "FIVE_PRIME_ONLY",
    "THREE_PRIME_ONLY",
    "PAIRED_END_SPLIT_PRIMERS",
    "UNABLE_TO_EXTRACT",
    "NO_REPORT",
)
_REQUIRED_ACTION_ORDER: tuple[str, ...] = (
    "NONE",
    "MANUAL_PRIMERS_REQUIRED",
    "READ_MERGING_RECOMMENDED",
    "NO_REPORT",
)


def _no_data_label(ax: Any, title: str) -> None:
    ax.set_title(f"{title} (no data)")
    ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)


def _panel_a_pair_recovery_by_status(ax: Any, pair_recovery: dict[str, Any]) -> None:
    """Stacked bar of pair_exact / equivalent / partial / failed per status."""
    counts: dict[str, dict[str, int]] = pair_recovery.get("counts", {})
    statuses = [s for s in _STATUS_ORDER if s in counts]
    if not statuses:
        _no_data_label(ax, "A · Pair recovery by status")
        return

    bottom = [0.0] * len(statuses)
    for bucket in _PAIR_BUCKETS:
        heights = [float(counts.get(s, {}).get(bucket, 0)) for s in statuses]
        if all(h == 0 for h in heights):
            continue
        ax.bar(statuses, heights, bottom=bottom, label=bucket)
        bottom = [b + h for b, h in zip(bottom, heights, strict=True)]

    ax.set_title("A · Pair recovery by LibraryReport.status")
    ax.set_xlabel("status")
    ax.set_ylabel("rows")
    ax.legend(loc="upper right", fontsize="x-small")


def _panel_b_breakdown(ax: Any, primer_recovery: dict[str, Any]) -> None:
    """Grouped bar: 5' EXACT / 3' EXACT / pair_exact / partial / failed across all rows.

    Computed from primer_recovery.counts_5p + counts_3p + the row-pair
    aggregator outputs the metrics report carries.
    """
    counts_5p = primer_recovery.get("counts_5p", {})
    counts_3p = primer_recovery.get("counts_3p", {})
    n_evaluated = primer_recovery.get("n_evaluated", 0)
    if n_evaluated == 0:
        _no_data_label(ax, "B · Per-side / pair recovery breakdown")
        return

    # Per-side EXACT counts.
    exact_5p = counts_5p.get("EXACT", 0)
    exact_3p = counts_3p.get("EXACT", 0)
    partial_5p = counts_5p.get("PARTIAL_5P", 0) + counts_5p.get("PARTIAL_3P", 0)
    partial_3p = counts_3p.get("PARTIAL_5P", 0) + counts_3p.get("PARTIAL_3P", 0)
    mismatch_5p = counts_5p.get("MISMATCH", 0) + counts_5p.get("IUPAC_UNSUPPORTED", 0)
    mismatch_3p = counts_3p.get("MISMATCH", 0) + counts_3p.get("IUPAC_UNSUPPORTED", 0)

    categories = ["5' EXACT", "3' EXACT", "5' partial", "3' partial", "5' miss", "3' miss"]
    values = [exact_5p, exact_3p, partial_5p, partial_3p, mismatch_5p, mismatch_3p]

    ax.bar(categories, values)
    ax.set_title(f"B · Per-side recovery breakdown (N={n_evaluated})")
    ax.set_ylabel("rows")
    ax.tick_params(axis="x", labelrotation=20)


def _panel_c_n_length_recovery(ax: Any, n_length_recovery: dict[str, Any]) -> None:
    """Bar of N-length recovery buckets."""
    in_tol = n_length_recovery.get("n_in_tolerance", 0)
    out_tol = n_length_recovery.get("n_out_of_tolerance", 0)
    unmeasurable = n_length_recovery.get("n_unmeasurable", 0)
    tolerance = n_length_recovery.get("tolerance", 2)
    total = in_tol + out_tol + unmeasurable
    if total == 0:
        _no_data_label(ax, "C · N-length recovery")
        return
    ax.bar(
        ["in tolerance", "out of tolerance", "unmeasurable"],
        [in_tol, out_tol, unmeasurable],
    )
    ax.set_title(f"C · N-length recovery (±{tolerance} nt)")
    ax.set_ylabel("rows")


def _panel_d_distributions(
    ax: Any, extraction_counts: dict[str, int], required_action_counts: dict[str, int]
) -> None:
    """Two horizontal bar groups in one axes: extraction_mode + required_action."""
    extraction_modes = [m for m in _EXTRACTION_MODE_ORDER if m in extraction_counts]
    required_actions = [a for a in _REQUIRED_ACTION_ORDER if a in required_action_counts]

    if not extraction_modes and not required_actions:
        _no_data_label(ax, "D · Extraction mode & required action")
        return

    # Render extraction_mode as a horizontal bar block, then required_action
    # below with a divider for readability.
    labels: list[str] = []
    heights: list[float] = []
    colors: list[str] = []
    palette_ext = "tab:blue"
    palette_req = "tab:orange"

    for m in extraction_modes:
        labels.append(f"ext · {m}")
        heights.append(float(extraction_counts[m]))
        colors.append(palette_ext)
    for a in required_actions:
        labels.append(f"req · {a}")
        heights.append(float(required_action_counts[a]))
        colors.append(palette_req)

    ax.barh(labels, heights, color=colors)
    ax.invert_yaxis()
    ax.set_title("D · Extraction mode & required action (honest accounting)")
    ax.set_xlabel("rows")
    # Light divider line between the two groups (legend would be redundant
    # given the prefix tags).


def plot_figure_a(metrics_json: Path, outdir: Path) -> tuple[Path, Path]:
    """Render the four-panel Figure A.

    Parameters
    ----------
    metrics_json
        Path to a ``metrics.json`` produced by
        :func:`selexprep.benchmark.metrics.write_metrics_json`.
    outdir
        Directory to write ``figure_a.{pdf,png}`` into.

    Returns
    -------
    A tuple ``(pdf_path, png_path)``.
    """
    metrics = json.loads(metrics_json.read_text(encoding="utf-8"))
    outdir.mkdir(parents=True, exist_ok=True)

    primer_recovery = metrics.get("primer_recovery", {})
    pair_recovery = metrics.get("pair_recovery_by_status", {})
    n_length_recovery = metrics.get("n_length_recovery", {})
    extraction_counts = metrics.get("extraction_mode_distribution", {}).get("counts", {})
    required_action_counts = metrics.get("required_action_distribution", {}).get("counts", {})

    fig, axes = plt.subplots(2, 2, figsize=(12, 10), dpi=150)
    _panel_a_pair_recovery_by_status(axes[0, 0], pair_recovery)
    _panel_b_breakdown(axes[0, 1], primer_recovery)
    _panel_c_n_length_recovery(axes[1, 0], n_length_recovery)
    _panel_d_distributions(axes[1, 1], extraction_counts, required_action_counts)

    # Headline: N verified + safe-failure rate.
    safe_failure = metrics.get("safe_failure_rate", {})
    n_verified = metrics.get("n_verified", 0)
    n_unverified = metrics.get("n_unverified", 0)
    rate = safe_failure.get("rate", 0.0)
    n_safe = safe_failure.get("n_safe_failures", 0)

    parts = ["selexprep Figure A — primer inference benchmark"]
    parts.append(f"N={n_verified} verified")
    if n_unverified:
        parts.append(f"{n_unverified} pending curation")
    parts.append(f"safe-failure rate: {rate:.0%} ({n_safe}/{n_verified} rows)")
    fig.suptitle("  ·  ".join(parts), fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    pdf_path = outdir / "figure_a.pdf"
    png_path = outdir / "figure_a.png"
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, bbox_inches="tight")
    plt.close(fig)
    return pdf_path, png_path


# ---------------------------------------------------------------------------
# CLI entry (called by the Snakefile's figure_a rule)
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Render the Phase 6 Figure A.")
    p.add_argument("--metrics", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    args = p.parse_args(argv)
    pdf, png = plot_figure_a(args.metrics, args.outdir)
    print(f"wrote {pdf}")
    print(f"wrote {png}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
