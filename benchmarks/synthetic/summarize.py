"""Summaries of a semi-synthetic run set (DESIGN.md, "Reported per configuration").

Reads ``<run dir>/results.tsv`` and writes, next to it:

- ``summary.tsv``: one row per configuration x donor x constant family, with
  the outcome counts over seeds, the median, minimum and maximum of each
  per-read metric for the inferred and the oracle extraction, the called
  constants' categories and the statuses. A refusal counts as a recovery of 0
  (it recovers nothing), so the statistics are not taken over accepted runs
  only.
- ``outcome_by_status.tsv``: outcome counts per ``status``, negative controls
  apart. ``confidence`` is not a calibrated probability; the table is
  descriptive.
- ``totals.tsv``: outcome counts over the whole set, with coverage (runs not
  refused) and accuracy among accepted runs.

    python benchmarks/synthetic/summarize.py benchmarks/synthetic/runs/test
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

OUTCOMES = [
    "complete, correct",
    "partial, correct and declared",
    "refused",
    "wrong",
    "error",
    "negative control, refused",
    "negative control, wrong",
]
CORRECT = ("complete, correct", "partial, correct and declared")


def summarize(run_dir: Path) -> None:
    d = pd.read_csv(run_dir / "results.tsv", sep="\t", keep_default_na=True)
    d["family"] = d["family"].fillna("none")
    # Rows in the order of configurations.tsv, as the design lists them.
    listed = pd.read_csv(Path(__file__).resolve().parent / "configurations.tsv", sep="\t")
    order = {c: i for i, c in enumerate(listed["config_id"])}

    rows = []
    for (config, donor, family), s in d.groupby(["config", "donor", "family"], sort=False):
        row = {"config": config, "donor": donor, "family": family, "runs": len(s)}
        counts = s["outcome"].value_counts()
        row.update({o: int(counts.get(o, 0)) for o in OUTCOMES})
        # Median and range over the runs of the cell: with two or three seeds a
        # median alone can sit between a failed and a successful run.
        for source in ("inferred", "oracle"):
            for metric in ("precision", "recall", "n_recovery", "lowest_round_yield"):
                col = f"{source}_{metric}"
                values = s[col].dropna() if col in s else pd.Series(dtype=float)
                for stat, fn in (("median", "median"), ("min", "min"), ("max", "max")):
                    row[f"{col}_{stat}"] = (
                        round(float(getattr(values, fn)()), 4) if len(values) else ""
                    )
        for side in ("5p", "3p"):
            col = f"call_{side}_category"
            if col in s:
                row[col] = ";".join(f"{k}={v}" for k, v in s[col].value_counts().items())
        row["statuses"] = ";".join(f"{k}={v}" for k, v in s["status"].value_counts().items())
        rows.append(row)
    summary = pd.DataFrame(rows)
    summary["_o"] = summary["config"].map(order)
    summary = summary.sort_values(["_o", "donor", "family"]).drop(columns="_o")
    summary.to_csv(run_dir / "summary.tsv", sep="\t", index=False)

    real = d[~d["config"].str.startswith("NC")]
    table = pd.crosstab(real["status"], real["outcome"])
    table.to_csv(run_dir / "outcome_by_status.tsv", sep="\t")

    counts = d["outcome"].value_counts()
    accepted = sum(int(counts.get(o, 0)) for o in (*CORRECT, "wrong"))
    correct = sum(int(counts.get(o, 0)) for o in CORRECT)
    totals = {o: int(counts.get(o, 0)) for o in OUTCOMES}
    totals.update(
        runs=len(d),
        runs_with_constants=len(real),
        coverage=round(accepted / len(real), 4) if len(real) else "",
        accuracy_among_accepted=round(correct / accepted, 4) if accepted else "",
        # Runs extracted on both sides with both boundaries right: a share of
        # runs, not of reads; the accuracy above also counts partial runs.
        complete_correct_share=round(int(counts.get("complete, correct", 0)) / len(real), 4)
        if len(real)
        else "",
    )
    pd.DataFrame([totals]).to_csv(run_dir / "totals.tsv", sep="\t", index=False)


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        summarize(Path(arg))
