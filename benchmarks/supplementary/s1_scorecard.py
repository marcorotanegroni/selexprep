# The generated document uses typographic primes and en dashes on purpose.
# ruff: noqa: RUF001
"""Supplementary Data S1: the Tier-1 scorecard, generated from the benchmark's data.

Every outcome, length, yield and total is read from the metrics and yields of
the two scored versions (the pre-specified v0.4.1 and the released 0.4.5) and
from ``ground_truth.tsv``, so S1 cannot fall behind the data it reports. The
only text written by hand, the legend and the per-deposit notes, belongs to the
manuscript and lives with it, outside the repository (``paper/`` is ignored);
without them the table is generated with no legend and no notes.

    python benchmarks/supplementary/s1_scorecard.py --out paper/supplementary/Supplementary_S1_scorecard.md
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
BENCH = HERE.parent
RESULTS = BENCH / "results"
PAPER_SUPP = BENCH.parent / "paper" / "supplementary"

ARMS = {
    "raw_standard": "recovery",
    "pre_trimmed": "specificity",
    "adapter_control": "adapter control",
}
ARM_ORDER = ["recovery", "specificity", "adapter control"]


def _tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _pairs(metrics: dict) -> dict[str, tuple[str, str]]:
    return {
        p["accession"]: (p["status_5p"]["equivalence_kind"], p["status_3p"]["equivalence_kind"])
        for p in metrics["primer_recovery"]["pairs"]
    }


def outcome(five: str, three: str) -> str:
    """As ``compute_pair_recovery_by_status``: both exact, partial, or missed."""
    if five == three == "EXACT":
        return "both exact"
    informative = {"EXACT", "PARTIAL_5P", "PARTIAL_3P"}
    if five in informative or three in informative:
        return "partial"
    return "missed"


def lowest_yield(path: Path) -> dict[str, float]:
    out: dict[str, float] = {}
    for r in _tsv(path):
        value = float(r["frac_kept"])
        out[r["accession"]] = min(value, out.get(r["accession"], value))
    return out


def _pct(value: float | None) -> str:
    if value is None:
        return "–"
    return f"{100 * value:.2f}%" if value < 0.01 else f"{100 * value:.0f}%"


def _n_lengths(metrics: dict) -> dict[str, str]:
    out = {}
    for r in metrics["n_length_recovery"]["per_row"]:
        out[r["accession"]] = (
            f"{r['observed']} / {r['truth']}" if r.get("observed") is not None else "–"
        )
    return out


def recovery_totals(pairs: dict[str, tuple[str, str]], accessions: list[str]) -> dict[str, int]:
    rows = [pairs[a] for a in accessions]
    kinds = [outcome(f, t) for f, t in rows]
    return {
        "n": len(rows),
        "both exact": kinds.count("both exact"),
        "partial": kinds.count("partial"),
        "missed": kinds.count("missed"),
        "5p exact": sum(f == "EXACT" for f, _ in rows),
        "3p exact": sum(t == "EXACT" for _, t in rows),
    }


def build(
    results: Path = RESULTS,
    bench: Path = BENCH,
    legend: Path = PAPER_SUPP / "s1_legend.md",
    notes_tsv: Path = PAPER_SUPP / "s1_notes.tsv",
) -> str:
    new = json.loads((results / "metrics.json").read_text())
    old = json.loads((results / "metrics_v0.4.1.json").read_text())
    gt = {r["accession"]: r for r in _tsv(bench / "ground_truth.tsv")}
    notes = {r["accession"]: r["note"] for r in _tsv(notes_tsv)} if notes_tsv.is_file() else {}
    new_pairs, old_pairs = _pairs(new), _pairs(old)
    new_n = _n_lengths(new)
    new_yield = lowest_yield(results / "extraction_yield.tsv")
    old_yield = lowest_yield(results / "extraction_yield_v0.4.1.tsv")

    controls = {}
    for block in ("specificity", "adapter_control"):
        for r in new[block]["per_row"]:
            controls[r["accession"]] = r

    scored = [
        r
        for r in gt.values()
        if r["verified"].strip().lower() == "true" and r["read_state"] in ARMS
    ]
    scored.sort(key=lambda r: (ARM_ORDER.index(ARMS[r["read_state"]]), r["accession"]))

    def changed(before: str, after: str) -> str:
        """One cell for two versions: the value, or 'v0.4.1 → 0.4.5' where they differ."""
        return after if before == after else f"{before} → {after}"

    def library(r: dict[str, str]) -> str:
        return f"{r['library_kind']}; {r['target_kind']}"

    lines: list[str] = []
    for arm in ARM_ORDER:
        rows = [r for r in scored if ARMS[r["read_state"]] == arm]
        lines += ["", f"### {arm[0].upper()}{arm[1:]} arm ({len(rows)})", ""]
        if arm == "recovery":
            lines += [
                "| Accession | Library; target | 5′ | 3′ | Outcome | *N* obs / truth "
                "| Lowest run yield | Note |",
                "|---|---|---|---|---|---|---|---|",
            ]
        else:
            lines += ["| Accession | Library; target | Outcome | Note |", "|---|---|---|---|"]
        for r in rows:
            acc, note = r["accession"], notes.get(r["accession"], "")
            if arm == "recovery":
                (o5, o3), (n5, n3) = old_pairs[acc], new_pairs[acc]
                cells = [
                    changed(o5, n5),
                    changed(o3, n3),
                    changed(outcome(o5, o3), outcome(n5, n3)),
                    new_n.get(acc, "–"),
                    changed(_pct(old_yield.get(acc)), _pct(new_yield.get(acc))),
                ]
            else:
                c = controls[acc]
                call = c["primer_5p"] is not None or c["primer_3p"] is not None
                cells = ["call (false)" if call else "no call (correct)"]
            lines.append("| " + " | ".join([acc, library(r), *cells, note]) + " |")

    recovery = [r["accession"] for r in scored if ARMS[r["read_state"]] == "recovery"]
    t_old, t_new = recovery_totals(old_pairs, recovery), recovery_totals(new_pairs, recovery)
    nl = new["n_length_recovery"]
    spec, adapt = new["specificity"], new["adapter_control"]

    def totals_line(label: str, t: dict[str, int]) -> str:
        return (
            f"*{label}:* {t['both exact']} with both constants exact, {t['partial']} partial, "
            f"{t['missed']} missed, of {t['n']}; per constant, {t['5p exact']}/{t['n']} exact "
            f"at 5′ and {t['3p exact']}/{t['n']} at 3′."
        )

    out = [
        "# Supplementary Data S1 — Tier-1 benchmark: per-deposit primer-inference scorecard",
        "",
        legend.read_text().strip() if legend.is_file() else "",
        "",
        *lines,
        "",
        "**Totals, recovery arm.**",
        "",
        "- " + totals_line("v0.4.1, pre-specified", t_old),
        "- " + totals_line("0.4.5, after the fix (in-sample)", t_new),
        (
            f"- Random-region length (0.4.5): {nl['n_in_tolerance']} of "
            f"{nl['n_in_tolerance'] + nl['n_out_of_tolerance']} measurable deposits at the "
            f"published length, {nl['n_out_of_tolerance']} outside the tolerance, "
            f"{nl['n_unmeasurable']} not measurable."
        ),
        "",
        (
            f"**Totals, controls (0.4.5).** Specificity: {spec['n_false_positive']} calls in "
            f"{spec['n_evaluated']} deposits. Adapter control: {adapt['n_false_positive']} calls "
            f"in {adapt['n_evaluated']} deposits."
        ),
        "",
    ]

    ex = new.get("excluded_after_inference")
    if ex:
        ex_pairs = _pairs({"primer_recovery": ex["primer_recovery"]})
        ex_n = _n_lengths({"n_length_recovery": ex["n_length_recovery"]})
        refused = set(ex["safe_failure_rate"]["safe_failure_accessions"])
        t_ex = recovery_totals(ex_pairs, sorted(ex_pairs))
        out += [
            "**Sensitivity analysis: the deposit excluded after inference (0.4.5).**",
            "",
            "| Accession | 5′ | 3′ | Refused | *N* obs / truth | Note |",
            "|---|---|---|---|---|---|",
        ]
        for acc in ex["accessions"]:
            f, t = ex_pairs[acc]
            out.append(
                f"| {acc} | {f} | {t} | {'yes' if acc in refused else 'no'} | "
                f"{ex_n.get(acc, '–')} | {notes.get(acc, '')} |"
            )
        out += [
            "",
            "- "
            + totals_line("Recovery arm with it counted in", t_ex).replace(
                "*Recovery arm with it counted in:*", "*Recovery arm with it counted in (0.4.5):*"
            ),
            "",
        ]

    out += [
        "**Excluded before inference** (`benchmarks/excluded_datasets.tsv`).",
        "",
        "| Accession | Reason | Detail |",
        "|---|---|---|",
    ]
    for r in _tsv(bench / "excluded_datasets.tsv"):
        detail = notes.get(r["accession"]) or r["notes"]
        out.append(f"| {r['accession']} | {r['reason'].replace('_', ' ')} | {detail} |")
    unverified = new.get("skipped_unverified_accessions", [])
    if unverified:
        out += ["", "**Not scored.**", ""]
        out += [
            f"- {acc}: {notes.get(acc, 'verified=false in ground_truth.tsv')}" for acc in unverified
        ]
    return "\n".join(out) + "\n"


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--legend", type=Path, default=PAPER_SUPP / "s1_legend.md")
    p.add_argument("--notes", type=Path, default=PAPER_SUPP / "s1_notes.tsv")
    args = p.parse_args(argv)
    for label, path in (("legend", args.legend), ("notes", args.notes)):
        if not path.is_file():
            print(f"no {label} at {path}: the table is generated without it")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(build(legend=args.legend, notes_tsv=args.notes))
    print(f"S1 -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
