"""QC orchestrator: manifest -> flags.yaml + 4 PNG plots.

The CLI's ``selexprep qc <manifest>`` verb is a thin wrapper around
:func:`run_qc`. The orchestrator:

1. Loads the manifest via :func:`selexprep.manifest.read_manifest_json`.
2. Auto-discovers ``round_*/counts.parquet`` under ``--counts-dir``
   (defaults to ``manifest_path.parent``).
3. Reads ``trim_reports.json`` for the read-retention plot.
4. Optionally reads ``strand_report.tsv`` for the strand-mix flag.
5. Computes the 8 suspicion flags.
6. Writes ``flags.yaml`` (deterministic, sorted by flag name).
7. Emits the 4 per-dataset PNG plots.

All filesystem expectations follow the /4 output layout - no
extra arguments needed in the common case.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from selexprep.manifest import read_manifest_json
from selexprep.qc.flags import Flag, compute_all_flags, write_flags_yaml
from selexprep.qc.plots import (
    plot_n_length_distribution,
    plot_per_round_panel,
    plot_primer_match_per_round,
    plot_read_retention,
)

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class QcResult:
    """Summary of one :func:`run_qc` invocation."""

    flags: list[Flag] = field(default_factory=list)
    flags_yaml_path: Path | None = None
    plot_paths: list[Path] = field(default_factory=list)

    @property
    def n_flags_raised(self) -> int:
        return len(self.flags)


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _load_counts_by_round(counts_dir: Path) -> dict[int, dict[str, int]]:
    """Load every ``round_*/counts.parquet`` under ``counts_dir``."""
    out: dict[int, dict[str, int]] = {}
    for parquet_path in sorted(counts_dir.glob("round_*/counts.parquet")):
        round_dir = parquet_path.parent.name
        if not round_dir.startswith("round_"):
            continue
        try:
            r = int(round_dir.split("_")[1])
        except (IndexError, ValueError):
            logger.warning("Skipping malformed round dir name: %s", round_dir)
            continue
        df = pd.read_parquet(parquet_path, columns=["sequence", "reads"])
        out[r] = dict(zip(df["sequence"], df["reads"].astype(int), strict=True))
    return out


def _load_trim_reports_by_round(trim_reports_path: Path) -> dict[int, dict[str, int]]:
    """Aggregate ``trim_reports.json`` entries by round."""
    if not trim_reports_path.exists():
        logger.info("No trim_reports.json at %s — skipping read-retention plot", trim_reports_path)
        return {}
    try:
        data = json.loads(trim_reports_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as e:
        logger.warning("Malformed trim_reports.json: %s", e)
        return {}
    by_round: dict[int, dict[str, int]] = {}
    for entry in data:
        paths = entry.get("output_paths") or []
        if not paths:
            continue
        round_dir = Path(paths[0]).parent.name
        if not round_dir.startswith("round_"):
            continue
        try:
            r = int(round_dir.split("_")[1])
        except (IndexError, ValueError):
            continue
        agg = by_round.setdefault(r, {"n_in": 0, "n_out": 0})
        agg["n_in"] += int(entry.get("n_in", 0))
        agg["n_out"] += int(entry.get("n_out", 0))
    return by_round


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def run_qc(
    manifest_path: Path,
    *,
    counts_dir: Path | None = None,
    outdir: Path | None = None,
) -> QcResult:
    """Run the QC pipeline for one dataset.

    Args:
        manifest_path: path to ``selexprep_manifest.json`` (emitted by
            ``extract``).
        counts_dir: directory containing ``round_*/counts.parquet`` files.
            Defaults to ``manifest_path.parent``.
        outdir: where ``flags.yaml`` + the 4 PNG plots are written.
            Defaults to ``manifest_path.parent / "qc"``.
    """
    manifest = read_manifest_json(manifest_path)
    counts_root = counts_dir if counts_dir is not None else manifest_path.parent
    qc_outdir = outdir if outdir is not None else manifest_path.parent / "qc"
    qc_outdir.mkdir(parents=True, exist_ok=True)

    counts_by_round = _load_counts_by_round(counts_root)
    trim_by_round = _load_trim_reports_by_round(manifest_path.parent / "trim_reports.json")

    strand_report_path: Path | None = manifest_path.parent / "strand_report.tsv"
    if not strand_report_path.exists():
        strand_report_path = None

    flags = compute_all_flags(
        manifest,
        counts_by_round,
        strand_report_path=strand_report_path,
        trim_reports_by_round=trim_by_round if trim_by_round else None,
    )
    flags_yaml_path = qc_outdir / "flags.yaml"
    write_flags_yaml(flags, flags_yaml_path)

    plot_paths: list[Path] = []
    if trim_by_round:
        plot_paths.append(plot_read_retention(trim_by_round, qc_outdir))
        plot_paths.append(plot_primer_match_per_round(manifest, trim_by_round, qc_outdir))
    plot_paths.append(plot_n_length_distribution(counts_by_round, qc_outdir))
    plot_paths.append(plot_per_round_panel(counts_by_round, qc_outdir))

    return QcResult(
        flags=flags,
        flags_yaml_path=flags_yaml_path,
        plot_paths=plot_paths,
    )
