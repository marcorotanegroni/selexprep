"""Figure A — primer-inference benchmark for the Application Note.

Scope pivot (2026-05-22): the figure tests selexprep's **unique claim**
(primer inference from accession reads + safe failure on ambiguity),
NOT comparator-tool count agreement. AptaPLEX / EasyDIVER+ both require
known primers as input; they cannot benchmark primer inference.

Phase 6b.10 reframe (Codex pass-3, 4 passes + user-approved): Figure A is
a **two-arm sensitivity/specificity benchmark**, NOT a read-state
taxonomy. The benchmark set (``ground_truth.tsv``) carries a ``read_state``
label — ``raw_standard`` (recovery arm) or ``pre_trimmed`` (specificity
arm) — assigned from independent read-length + architecture evidence,
never from detect output. Out-of-scope deposits (nonstandard
architecture, aggregate-like) were removed to ``excluded_datasets.tsv``
and are NOT figure arms.

Four panels arranged 2x2:

- **Panel A (top-left)** — recovery arm (sensitivity). Pair recovery
  (exact / equivalent / partial / miss, kept separate) over the
  ``raw_standard`` rows, with a multi-round-only overlay so
  confidence-limited mono-round rows don't silently deflate the headline.

- **Panel B (top-right)** — specificity arm. False-positive primer calls
  on ``pre_trimmed`` deposits (target = 0) + the per-side (5'/3')
  no-call breakdown.

- **Panel C (bottom-left)** — N-length recovery within ±tolerance, on the
  ``raw_standard`` arm.

- **Panel D (bottom-right)** — honest accounting: ``extraction_mode`` +
  ``required_action`` distributions across all verified rows, plus a
  fetch_stats note for partial-fetch deposits (e.g. PRJEB70964 17/27).

The figure title carries the headline: complete/partial recovery on the
``raw_standard`` N + the specificity false-positive count. NEVER a flat
"X/11" (the set is not a random sample).

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

# Stable, ordered buckets so the visual order doesn't reshuffle between runs.
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

# Recovery semantics (Codex pass-4: keep exact and equivalent SEPARATE —
# do not collapse into a single "complete"): both sides EXACT → "exact";
# both sides matched via an equivalence rule (revcomp / U-T / barcode) →
# "equivalent"; exactly one side matched → "partial"; neither → "miss".
_RECOVERY_COLORS = {
    "exact": "tab:green",
    "equivalent": "yellowgreen",
    "partial": "tab:orange",
    "miss": "tab:red",
}
_SPECIFICITY_COLORS = {"no_false_call": "tab:green", "false_positive": "tab:red"}


def _no_data_label(ax: Any, title: str) -> None:
    ax.set_title(f"{title} (no data)")
    ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)


def _collapse_pair_counts(counts: dict[str, dict[str, int]]) -> dict[str, int]:
    """Sum status-bucketed pair recovery into exact / equivalent / partial / miss.

    ``counts`` is ``{status: {pair_exact/pair_equivalent/pair_partial/
    pair_failed: n}}`` (the ``pair_recovery_by_status`` cross-tab). The
    recovery arm cares about the pair-level outcome, not the status
    bucket, so we collapse over statuses — but keep exact and equivalent
    distinct (Codex pass-4: don't fold them into one "complete").
    """
    exact = equivalent = partial = miss = 0
    for bucket_counts in counts.values():
        exact += int(bucket_counts.get("pair_exact", 0))
        equivalent += int(bucket_counts.get("pair_equivalent", 0))
        partial += int(bucket_counts.get("pair_partial", 0))
        miss += int(bucket_counts.get("pair_failed", 0))
    return {"exact": exact, "equivalent": equivalent, "partial": partial, "miss": miss}


def _panel_a_recovery(ax: Any, pair_recovery: dict[str, Any], multi_round: dict[str, Any]) -> None:
    """Recovery arm: exact/equivalent/partial/miss pair recovery on raw_standard.

    Two stacked bars — all raw_standard rows + the multi-round-only
    subset (mono-round rows excluded) — so the overlay shows whether the
    confidence-limited single-round rows are dragging the headline.
    """
    all_counts = _collapse_pair_counts(pair_recovery.get("counts", {}))
    mr_counts = _collapse_pair_counts(multi_round.get("counts", {}))
    n_all = sum(all_counts.values())
    n_mr = sum(mr_counts.values())
    if n_all == 0:
        _no_data_label(ax, "A · Recovery (raw_standard)")
        return

    labels = [f"raw_standard\n(all, N={n_all})", f"multi-round only\n(N={n_mr})"]
    series = [all_counts, mr_counts]
    bottoms = [0.0, 0.0]
    for key in ("exact", "equivalent", "partial", "miss"):
        heights = [float(s.get(key, 0)) for s in series]
        ax.bar(labels, heights, bottom=bottoms, label=key, color=_RECOVERY_COLORS[key])
        bottoms = [b + h for b, h in zip(bottoms, heights, strict=True)]

    ax.set_title("A · Recovery arm — pair recovery (sensitivity)")
    ax.set_ylabel("rows")
    ax.legend(loc="upper right", fontsize="x-small")


def _panel_b_specificity(ax: Any, specificity: dict[str, Any]) -> None:
    """Specificity arm: false-positive primer calls on pre_trimmed (target 0).

    Grouped bar over {pair, 5' side, 3' side}: no-false-call (green,
    correct) vs false-positive (red). A pre-trimmed deposit's reads are
    the N-region only, so the correct behavior is to emit NO primer.
    """
    n_eval = int(specificity.get("n_evaluated", 0))
    if n_eval == 0:
        _no_data_label(ax, "B · Specificity (pre_trimmed)")
        return

    n_no_call = int(specificity.get("n_no_false_call", 0))
    n_fp = int(specificity.get("n_false_positive", 0))

    # Per-side breakdown from per_row (primer_5p / primer_3p None = no call).
    per_row = specificity.get("per_row", [])
    n5_emit = sum(1 for x in per_row if x.get("primer_5p") is not None)
    n3_emit = sum(1 for x in per_row if x.get("primer_3p") is not None)
    n5_no = n_eval - n5_emit
    n3_no = n_eval - n3_emit

    groups = ["pair", "5' side", "3' side"]
    no_call = [n_no_call, n5_no, n3_no]
    false_pos = [n_fp, n5_emit, n3_emit]

    x = range(len(groups))
    width = 0.38
    ax.bar(
        [i - width / 2 for i in x],
        no_call,
        width,
        label="no false call",
        color=_SPECIFICITY_COLORS["no_false_call"],
    )
    ax.bar(
        [i + width / 2 for i in x],
        false_pos,
        width,
        label="false positive",
        color=_SPECIFICITY_COLORS["false_positive"],
    )
    ax.set_xticks(list(x))
    ax.set_xticklabels(groups)
    ax.set_title(f"B · Specificity arm — {n_fp}/{n_eval} false-positive calls (target 0)")
    ax.set_ylabel("pre_trimmed deposits")
    ax.legend(loc="upper center", fontsize="x-small")


def _panel_c_n_length_recovery(ax: Any, n_length_recovery: dict[str, Any]) -> None:
    """N-length recovery buckets on the raw_standard arm."""
    in_tol = int(n_length_recovery.get("n_in_tolerance", 0))
    out_tol = int(n_length_recovery.get("n_out_of_tolerance", 0))
    unmeasurable = int(n_length_recovery.get("n_unmeasurable", 0))
    tolerance = n_length_recovery.get("tolerance", 2)
    total = in_tol + out_tol + unmeasurable
    if total == 0:
        _no_data_label(ax, "C · N-length recovery")
        return
    ax.bar(
        ["in tolerance", "out of tolerance", "unmeasurable"],
        [in_tol, out_tol, unmeasurable],
        color=["tab:green", "tab:red", "tab:gray"],
    )
    ax.set_title(f"C · N-length recovery (±{tolerance} nt, raw_standard)")
    ax.set_ylabel("rows")


def _panel_d_distributions(
    ax: Any,
    extraction_counts: dict[str, int],
    required_action_counts: dict[str, int],
    fetch_stats: dict[str, Any],
) -> None:
    """Honest accounting: extraction_mode + required_action + a fetch note."""
    extraction_modes = [m for m in _EXTRACTION_MODE_ORDER if m in extraction_counts]
    required_actions = [a for a in _REQUIRED_ACTION_ORDER if a in required_action_counts]

    # Partial-fetch note (e.g. PRJEB70964 17/27) — surfaced from fetch_stats.
    partial_notes = [
        f"{acc}: {fs.get('fetch_available_runs')}/{fs.get('fetch_expected_runs')} runs"
        for acc, fs in sorted(fetch_stats.items())
        if fs.get("partial_fetch")
    ]

    if not extraction_modes and not required_actions:
        _no_data_label(ax, "D · Extraction mode & required action")
        if partial_notes:
            ax.text(
                0.5,
                0.1,
                "partial fetch — " + "; ".join(partial_notes),
                ha="center",
                va="bottom",
                fontsize="x-small",
                transform=ax.transAxes,
            )
        return

    labels: list[str] = []
    heights: list[float] = []
    colors: list[str] = []
    for m in extraction_modes:
        labels.append(f"ext · {m}")
        heights.append(float(extraction_counts[m]))
        colors.append("tab:blue")
    for a in required_actions:
        labels.append(f"req · {a}")
        heights.append(float(required_action_counts[a]))
        colors.append("tab:orange")

    ax.barh(labels, heights, color=colors)
    ax.invert_yaxis()
    ax.set_title("D · Extraction mode & required action (honest accounting)")
    ax.set_xlabel("rows")
    if partial_notes:
        ax.text(
            0.98,
            0.02,
            "partial fetch — " + "; ".join(partial_notes),
            ha="right",
            va="bottom",
            fontsize="x-small",
            style="italic",
            transform=ax.transAxes,
        )


def plot_figure_a(metrics_json: Path, outdir: Path) -> tuple[Path, Path]:
    """Render the four-panel two-arm Figure A.

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

    pair_recovery = metrics.get("pair_recovery_by_status", {})
    multi_round = metrics.get("multi_round_sensitivity", {})
    specificity = metrics.get("specificity", {})
    n_length_recovery = metrics.get("n_length_recovery", {})
    extraction_counts = metrics.get("extraction_mode_distribution", {}).get("counts", {})
    required_action_counts = metrics.get("required_action_distribution", {}).get("counts", {})
    fetch_stats = metrics.get("fetch_stats", {})

    fig, axes = plt.subplots(2, 2, figsize=(12, 10), dpi=150)
    _panel_a_recovery(axes[0, 0], pair_recovery, multi_round)
    _panel_b_specificity(axes[0, 1], specificity)
    _panel_c_n_length_recovery(axes[1, 0], n_length_recovery)
    _panel_d_distributions(axes[1, 1], extraction_counts, required_action_counts, fetch_stats)

    # Headline: recovery (exact/equivalent/partial on raw_standard N, kept
    # separate per Codex pass-4) + specificity false-positive count. NEVER a
    # flat "X/11" — the set isn't a random sample, so a single recovery
    # fraction would over-read it.
    recovery_n = int(metrics.get("recovery_denominator", 0))
    collapsed = _collapse_pair_counts(pair_recovery.get("counts", {}))
    spec_n = int(specificity.get("n_evaluated", 0))
    spec_fp = int(specificity.get("n_false_positive", 0))

    parts = ["selexprep Figure A — primer-inference benchmark"]
    parts.append(
        f"recovery: {collapsed['exact']} exact / {collapsed['equivalent']} equiv / "
        f"{collapsed['partial']} partial of {recovery_n} raw_standard"
    )
    parts.append(f"specificity: {spec_fp}/{spec_n} false-positive calls on pre_trimmed")
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
