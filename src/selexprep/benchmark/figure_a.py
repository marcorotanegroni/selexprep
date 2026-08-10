"""Tier 1 primer-recovery scorecard (Application Note **Table 1**).

Emits a per-deposit Markdown scorecard by joining the benchmark's deterministic
``metrics.json`` (per-side recovery outcomes) with the descriptors in
``ground_truth.tsv`` (chemistry / target / arm). The benchmark reports results
as a **table, not a chart**: deposits with categorical per-side outcomes are
far more legible per-row than as bars, and the table carries the per-deposit
detail a reader actually wants. ``metrics.json`` stays the machine-readable
source of truth; this is its presentation layer.
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Any

from selexprep.benchmark.metrics import load_ground_truth

logger = logging.getLogger(__name__)

# read_state (ground_truth) → benchmark arm shown in the scorecard.
_ARM = {
    "raw_standard": "recovery",
    "pre_trimmed": "specificity",
    "adapter_control": "adapter-control",
}

_COLUMNS = ("accession", "chemistry", "target", "arm", "five", "three", "note")
_HEADERS = ("Accession", "Chemistry", "Target", "Arm", "5'", "3'", "Note")


def _collapse_pair_counts(counts: dict[str, dict[str, int]]) -> dict[str, int]:
    """Sum the ``pair_recovery_by_status`` cross-tab into exact/equivalent/partial/miss."""
    out = {"exact": 0, "equivalent": 0, "partial": 0, "miss": 0}
    for bucket in counts.values():
        out["exact"] += int(bucket.get("pair_exact", 0))
        out["equivalent"] += int(bucket.get("pair_equivalent", 0))
        out["partial"] += int(bucket.get("pair_partial", 0))
        out["miss"] += int(bucket.get("pair_failed", 0))
    return out


def _side(side: Any) -> str:
    """Equivalence kind for one primer side (``EXACT`` / ``MISMATCH`` / …), or ``—``."""
    if not isinstance(side, dict):
        return "—"
    return str(side.get("equivalence_kind") or "—")


def _headline(metrics: dict[str, Any]) -> str:
    """One-line factual summary (counts only — no interpretation)."""
    rec = _collapse_pair_counts(metrics.get("pair_recovery_by_status", {}).get("counts", {}))
    denom = int(metrics.get("recovery_denominator", 0))
    spec = metrics.get("specificity", {})
    ac = metrics.get("adapter_control", {})
    parts = [
        f"recovery: {rec['exact']} exact / {rec['equivalent']} equivalent / "
        f"{rec['partial']} partial of {denom} evaluable",
        f"specificity: {int(spec.get('n_false_positive', 0))}/{int(spec.get('n_evaluated', 0))} "
        "false-positive calls",
    ]
    if int(ac.get("n_evaluated", 0)):
        parts.append(
            f"adapter-control: {int(ac.get('n_no_false_call', 0))}/{int(ac.get('n_evaluated', 0))} "
            "correct refusals"
        )
    return "; ".join(parts)


def build_scorecard(metrics: dict[str, Any], gt_rows: list[Any]) -> list[dict[str, str]]:
    """Join ``metrics.json`` outcomes with ground-truth descriptors, one row per verified deposit."""
    pairs = {p["accession"]: p for p in metrics.get("primer_recovery", {}).get("pairs", [])}
    controls: dict[str, dict[str, Any]] = {}
    for block in ("specificity", "adapter_control"):
        for r in metrics.get(block, {}).get("per_row", []):
            controls[r["accession"]] = r

    rows: list[dict[str, str]] = []
    for gt in gt_rows:
        if not gt.verified:
            continue
        arm = _ARM.get(gt.read_state, gt.read_state or "—")
        note = ""
        if gt.read_state == "raw_standard":
            p = pairs.get(gt.accession, {})
            five = _side(p.get("status_5p"))
            if p.get("score_3p", True):
                three = _side(p.get("status_3p"))
            else:
                three = "—"
                note = "3' read-resolved → scored on 5'"
        else:
            r = controls.get(gt.accession, {})
            five = "null" if r.get("primer_5p") is None else "CALL"
            three = "null" if r.get("primer_3p") is None else "CALL"
            note = "correct refusal" if five == "null" and three == "null" else "FALSE POSITIVE"
        rows.append(
            {
                "accession": gt.accession,
                "chemistry": gt.library_kind or "—",
                "target": gt.target_kind or "—",
                "arm": arm,
                "five": five,
                "three": three,
                "note": note,
            }
        )
    # Group by arm (recovery → specificity → adapter-control), then accession.
    arm_rank = {"recovery": 0, "specificity": 1, "adapter-control": 2}
    rows.sort(key=lambda r: (arm_rank.get(r["arm"], 9), r["accession"]))
    return rows


def emit_scorecard(metrics_json: Path, ground_truth: Path, outdir: Path) -> Path:
    """Write ``table_1.md`` — the per-deposit primer-recovery scorecard. Returns its path."""
    metrics = json.loads(metrics_json.read_text(encoding="utf-8"))
    gt_rows = load_ground_truth(ground_truth)
    rows = build_scorecard(metrics, gt_rows)
    outdir.mkdir(parents=True, exist_ok=True)

    lines = [
        "# Table 1 — primer-recovery scorecard",
        "",
        f"_{_headline(metrics)}._",
        "",
        "| " + " | ".join(_HEADERS) + " |",
        "|" + "|".join(["---"] * len(_HEADERS)) + "|",
    ]
    lines += ["| " + " | ".join(r[c] for c in _COLUMNS) + " |" for r in rows]
    lines.append("")

    out_path = outdir / "table_1.md"
    out_path.write_text("\n".join(lines), encoding="utf-8")
    return out_path


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Emit the Tier 1 primer-recovery scorecard (Table 1).")
    p.add_argument("--metrics", required=True, type=Path)
    p.add_argument("--ground-truth", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    args = p.parse_args(argv)
    out_path = emit_scorecard(args.metrics, args.ground_truth, args.outdir)
    print(f"wrote {out_path}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())


__all__ = ["build_scorecard", "emit_scorecard", "main"]
