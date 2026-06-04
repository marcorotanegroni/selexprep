"""Four per-dataset QC plots.

All plots write PNGs to ``<outdir>/qc/``:

1. ``read_retention.png`` - bars: n_in vs n_out per round (trim retention).
2. ``primer_match_per_round.png`` - per-round n_out/n_in line +
   manifest's overall 5'/3' match rate reference lines.
3. ``n_length_distribution.png`` - per-round N-region length
   histograms (faceted).
4. ``per_round_panel.png`` - 3-subplot panel: unique-seq bars,
   Shannon entropy line, top-N coverage line.

**Determinism note**: matplotlib PNG output is NOT byte-deterministic
across versions (PNG metadata includes timestamps). These plots are
informational.

**QC artifact lifecycle**: ``flags.yaml`` and the 4 PNG plots are
emitted AFTER ``selexprep extract`` has already sealed
``selexprep_manifest.json``. The manifest's
``output_sha256`` therefore does NOT include the QC outputs by design
- they're a post-hoc QC report, not part of the immutable extract
provenance. ``flags.yaml`` IS deterministic by construction (sorted
keys, rounded floats); callers can hash it independently via
``selexprep._io.sha256_file`` if they need to record a QC checksum.
``.yaml`` is included in ``_HASHABLE_SUFFIXES`` so a future
``selexprep qc-amend`` command (v0.2) could append the YAML hash to
the manifest if that becomes useful.
"""

from __future__ import annotations

import logging
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # Non-interactive backend (no GUI dependency in tests).
import matplotlib.pyplot as plt

from selexprep.manifest import SelexprepManifestV1
from selexprep.qc.diversity import shannon_entropy, top_n_coverage, unique_count

logger = logging.getLogger(__name__)


# CALIBRATION-TODO: not in the design; top-N coverage uses top-100 by
# convention. confirmed.
TOP_N_COVERAGE_N = 100

_DPI = 150
_DEFAULT_FIGSIZE = (6, 4)


def _save_fig(fig: plt.Figure, path: Path) -> Path:
    """Common save logic: parent dir, tight bbox, DPI 150, close fig."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=_DPI, bbox_inches="tight")
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# Plot 1: read retention per stage
# ---------------------------------------------------------------------------


def plot_read_retention(
    trim_reports_by_round: dict[int, dict[str, int]],
    outdir: Path,
) -> Path:
    """Grouped bar chart: n_in vs n_out per round (cutadapt retention)."""
    rounds = sorted(trim_reports_by_round)
    n_in = [trim_reports_by_round[r]["n_in"] for r in rounds]
    n_out = [trim_reports_by_round[r]["n_out"] for r in rounds]

    fig, ax = plt.subplots(figsize=_DEFAULT_FIGSIZE)
    width = 0.4
    xs = [i for i, _ in enumerate(rounds)]
    ax.bar([x - width / 2 for x in xs], n_in, width=width, label="input")
    ax.bar([x + width / 2 for x in xs], n_out, width=width, label="extracted")
    ax.set_xticks(xs)
    ax.set_xticklabels([f"R{r}" for r in rounds])
    ax.set_xlabel("Round")
    ax.set_ylabel("Reads")
    ax.set_title("Read retention per round (cutadapt --discard-untrimmed)")
    ax.legend()
    return _save_fig(fig, outdir / "read_retention.png")


# ---------------------------------------------------------------------------
# Plot 2: primer match rate per round
# ---------------------------------------------------------------------------


def plot_primer_match_per_round(
    manifest: SelexprepManifestV1,
    trim_reports_by_round: dict[int, dict[str, int]],
    outdir: Path,
) -> Path:
    """Per-round trim retention as a proxy for primer match rate.

    v0.1 caveat: the LibraryReport stores only the earliest round's
    match_rate_5p / match_rate_3p. Per-round trim retention
    (``n_out / n_in``) is the closest proxy at QC time. The LR's
    overall rates are shown as horizontal reference lines for context.
    """
    rounds = sorted(trim_reports_by_round)
    retention = [
        (trim_reports_by_round[r]["n_out"] / trim_reports_by_round[r]["n_in"])
        if trim_reports_by_round[r]["n_in"] > 0
        else 0.0
        for r in rounds
    ]
    lr = manifest.library_report

    fig, ax = plt.subplots(figsize=_DEFAULT_FIGSIZE)
    ax.plot(rounds, retention, marker="o", label="trim retention (n_out / n_in)")
    ax.axhline(
        y=lr.match_rate_5p,
        linestyle="--",
        color="tab:blue",
        alpha=0.6,
        label=f"5' match rate (overall = {lr.match_rate_5p:.2f})",
    )
    ax.axhline(
        y=lr.match_rate_3p,
        linestyle="--",
        color="tab:orange",
        alpha=0.6,
        label=f"3' match rate (overall = {lr.match_rate_3p:.2f})",
    )
    ax.set_xlabel("Round")
    ax.set_ylabel("Fraction (0-1)")
    ax.set_ylim(0, 1.05)
    ax.set_title("Primer match / trim retention per round")
    ax.legend(loc="lower left", fontsize=8)
    return _save_fig(fig, outdir / "primer_match_per_round.png")


# ---------------------------------------------------------------------------
# Plot 3: N-region length distribution
# ---------------------------------------------------------------------------


def plot_n_length_distribution(
    counts_by_round: dict[int, dict[str, int]],
    outdir: Path,
) -> Path:
    """Per-round N-region length histograms, faceted (one subplot per round)."""
    rounds = sorted(counts_by_round)
    if not rounds:
        # Empty input - emit an empty placeholder figure for stability.
        fig, ax = plt.subplots(figsize=_DEFAULT_FIGSIZE)
        ax.text(0.5, 0.5, "(no rounds)", ha="center", va="center", transform=ax.transAxes)
        ax.set_title("N-region length distribution per round")
        return _save_fig(fig, outdir / "n_length_distribution.png")

    n = len(rounds)
    cols = min(3, n)
    rows = (n + cols - 1) // cols
    fig, axes_arr = plt.subplots(rows, cols, figsize=(cols * 4, rows * 3), squeeze=False)

    for idx, r in enumerate(rounds):
        ax = axes_arr[idx // cols][idx % cols]
        # Length distribution weighted by read count.
        length_counter: Counter[int] = Counter()
        for seq, count in counts_by_round[r].items():
            length_counter[len(seq)] += count
        if not length_counter:
            ax.text(0.5, 0.5, "(empty)", ha="center", va="center", transform=ax.transAxes)
        else:
            lengths = sorted(length_counter)
            weights = [length_counter[L] for L in lengths]
            ax.bar(lengths, weights, width=0.8)
        ax.set_xlabel("N-region length")
        ax.set_ylabel("Reads")
        ax.set_title(f"Round {r}")

    # Hide any unused subplots.
    for idx in range(n, rows * cols):
        axes_arr[idx // cols][idx % cols].axis("off")

    fig.suptitle("N-region length distribution per round")
    fig.tight_layout()
    return _save_fig(fig, outdir / "n_length_distribution.png")


# ---------------------------------------------------------------------------
# Plot 4: per-round diversity panel
# ---------------------------------------------------------------------------


def plot_per_round_panel(
    counts_by_round: dict[int, dict[str, int]],
    outdir: Path,
) -> Path:
    """3-subplot panel: unique sequences (bar), Shannon entropy (line),
    top-N coverage (line) per round."""
    rounds = sorted(counts_by_round)
    uniques = [unique_count(counts_by_round[r]) for r in rounds]
    entropies = [shannon_entropy(counts_by_round[r]) for r in rounds]
    coverage = [top_n_coverage(counts_by_round[r], TOP_N_COVERAGE_N) for r in rounds]

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4))

    ax1.bar(rounds, uniques)
    ax1.set_xlabel("Round")
    ax1.set_ylabel("Unique sequences")
    ax1.set_title("Unique count per round")
    ax1.set_xticks(rounds)

    ax2.plot(rounds, entropies, marker="o")
    ax2.set_xlabel("Round")
    ax2.set_ylabel("Shannon entropy (bits)")
    ax2.set_title("Shannon entropy per round")
    ax2.set_xticks(rounds)

    ax3.plot(rounds, coverage, marker="o", color="tab:red")
    ax3.set_xlabel("Round")
    ax3.set_ylabel(f"Top-{TOP_N_COVERAGE_N} coverage (fraction)")
    ax3.set_title(f"Top-{TOP_N_COVERAGE_N} coverage per round")
    ax3.set_ylim(0, 1.05)
    ax3.set_xticks(rounds)

    fig.tight_layout()
    return _save_fig(fig, outdir / "per_round_panel.png")
