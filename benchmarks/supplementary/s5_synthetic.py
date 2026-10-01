# The generated document uses typographic primes and en dashes on purpose.
# ruff: noqa: RUF001
"""Supplementary Data S5: the semi-synthetic benchmark, generated from its results.

Every count and median is read from the run tables in
``benchmarks/synthetic/runs/`` (the single test run, the development runs and
the two post-test experiments), and the registered expectations are copied
from ``benchmarks/synthetic/DESIGN.md``, so S5 cannot fall behind the data. The
legend is written by hand and lives with the manuscript, outside the
repository (``paper/`` is ignored); without it the tables are generated alone.

    python benchmarks/supplementary/s5_synthetic.py --out paper/supplementary/Supplementary_S5_synthetic.md
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
SYN = HERE.parent / "synthetic"
RUNS = SYN / "runs"
PAPER_SUPP = HERE.parent.parent / "paper" / "supplementary"

SHORT = {
    "complete, correct": "complete",
    "partial, correct and declared": "partial",
    "refused": "refused",
    "wrong": "wrong",
    "error": "error",
}
NC = ("negative control, refused", "negative control, wrong")


def _pct(x: float | None) -> str:
    return "–" if x is None or pd.isna(x) else f"{100 * x:.1f}%"


def _counts(frame: pd.DataFrame) -> dict[str, int]:
    return {short: int((frame.outcome == long).sum()) for long, short in SHORT.items()}


def _table(header: list[str], rows: list[list[str]]) -> str:
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(lines)


def totals_line(frame: pd.DataFrame) -> str:
    with_constants = frame[~frame.outcome.isin(NC)]
    c = _counts(with_constants)
    accepted = c["complete"] + c["partial"] + c["wrong"]
    nc = frame[frame.outcome.isin(NC)]
    return (
        f"{len(frame):,} runs, {len(with_constants):,} with constants: {c['complete']:,} complete, "
        f"{c['partial']:,} partial, {c['refused']:,} refused, {c['wrong']:,} wrong, "
        f"{c['error']} pipeline errors; coverage {_pct((len(with_constants) - c['refused']) / len(with_constants))}, "
        f"accuracy among accepted runs {_pct((c['complete'] + c['partial']) / accepted)}. "
        f"Negative controls: {int((nc.outcome == NC[0]).sum())} of {len(nc)} refused."
    )


def per_configuration(test: pd.DataFrame) -> str:
    configs = list(csv.DictReader((SYN / "configurations.tsv").open(), delimiter="\t"))
    rows = []
    for c in configs:
        frame = test[test.config == c["config_id"]]
        if c["config_id"].startswith("NC"):
            refused = int((frame.outcome == NC[0]).sum())
            rows.append(
                [
                    c["config_id"],
                    c["description"],
                    len(frame),
                    "–",
                    "–",
                    f"{refused} (control)",
                    len(frame) - refused,
                    "–",
                    "–",
                ]
            )
            continue
        n = _counts(frame)
        rows.append(
            [
                c["config_id"],
                c["description"],
                len(frame),
                n["complete"],
                n["partial"],
                n["refused"],
                n["wrong"],
                _pct(frame.inferred_n_recovery.median()),
                _pct(frame.oracle_n_recovery.median()),
            ]
        )
    return _table(
        [
            "Configuration",
            "Condition",
            "Runs",
            "Complete",
            "Partial",
            "Refused",
            "Wrong",
            "N recovered (median)",
            "Oracle (median)",
        ],
        rows,
    )


def stratified(test: pd.DataFrame, column: str, label: str) -> str:
    frame = test[~test.outcome.isin(NC)]
    rows = []
    for key, g in frame.groupby(column):
        n = _counts(g)
        rows.append([key, len(g), n["complete"], n["partial"], n["refused"], n["wrong"]])
    return _table([label, "Runs", "Complete", "Partial", "Refused", "Wrong"], rows)


def by_status() -> str:
    t = pd.read_csv(RUNS / "test" / "outcome_by_status.tsv", sep="\t")
    rows = [[r.status, *(int(r[k]) for k in SHORT if k in t.columns)] for _, r in t.iterrows()]
    header = ["Status", *(SHORT[k] for k in SHORT if k in t.columns)]
    return _table(header, rows)


def expectations() -> str:
    """The expectations table of DESIGN.md, copied as registered."""
    lines = (SYN / "DESIGN.md").read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("| ID | Desirable |"))
    out = []
    for line in lines[start:]:
        if not line.startswith("|"):
            break
        out.append(line)
    return "\n".join(out)


def post_test() -> str:
    over = pd.read_csv(RUNS / "post-override" / "results.tsv", sep="\t")
    traj = pd.read_csv(RUNS / "post-trajectory" / "results.tsv", sep="\t")
    rows = []
    for config, g in over.groupby("config"):
        same = int((g.override_n_recovery == g.oracle_n_recovery).sum())
        n = _counts(g)
        rows.append(
            [
                f"P1, {config}",
                len(g),
                f"inferred: {n['complete']} complete, {n['refused']} refused, {n['wrong']} wrong",
                f"overrides equal the oracle in {same} of {len(g)}; N recovered (median) {_pct(g.override_n_recovery.median())}",
            ]
        )
    for config, g in traj.groupby("config"):
        n = _counts(g)
        rows.append(
            [
                f"P2, {config}",
                len(g),
                f"{n['complete']} complete, {n['refused']} refused, {n['wrong']} wrong",
                f"N recovered (median) {_pct(g.inferred_n_recovery.median())}",
            ]
        )
    return _table(["Experiment", "Runs", "Outcomes", "Extraction"], rows)


def build(legend: Path) -> str:
    test = pd.read_csv(RUNS / "test" / "results.tsv", sep="\t")
    dev = pd.read_csv(RUNS / "dev-2" / "results.tsv", sep="\t")
    parts = ["# Supplementary Data S5 — semi-synthetic benchmark", ""]
    if legend.is_file():
        parts += [legend.read_text().strip(), ""]
    parts += [
        "## Test set: totals",
        "",
        totals_line(test),
        "",
        "## Test set: outcome per configuration",
        "",
        "Counts over donors, constant families and seeds. *N recovered*: share of a run's reads "
        "whose extracted sequence equals the true random region, a refused run counting as 0; "
        "*Oracle*: the same reads extracted with the true constants.",
        "",
        per_configuration(test),
        "",
        "## Test set: outcome per donor and per constant family",
        "",
        stratified(test, "donor", "Donor"),
        "",
        stratified(test, "family", "Constant family"),
        "",
        "## Test set: outcome per status",
        "",
        by_status(),
        "",
        "## Expectations registered before the test",
        "",
        "Copied from `benchmarks/synthetic/DESIGN.md` as committed before any library was generated.",
        "",
        expectations(),
        "",
        "## Experiments after the test",
        "",
        post_test(),
        "",
        "## Development set",
        "",
        "Development run 2, the last before the test: " + totals_line(dev),
        "",
    ]
    return "\n".join(parts)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--legend", type=Path, default=PAPER_SUPP / "s5_legend.md")
    args = p.parse_args(argv)
    args.out.write_text(build(args.legend))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
