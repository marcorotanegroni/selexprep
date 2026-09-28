"""Donor screen, metadata steps (DESIGN.md, "Donor screen", steps 1, 2, 6, 7).

Uses the bundled catalogue's curated metadata and the ENA run records only;
nothing here looks at how selexprep behaves on the reads. Steps 3-5 (evidence
that the reads are the random region, technical residue) need the reads and the
documentation and are recorded afterwards in the same table.

Rules, fixed in DESIGN.md before this script was run:

1. Candidate: an INSDC deposit whose curated ``n_random`` states one length N.
   The field is free text. A single integer in it is that length; when the text
   holds no integer or several of them, the value is reviewed by hand in
   ``n_random_review.tsv``, which is authoritative for every accession it lists
   (amendment 2, item 9: "12-mer random sequence (within a 71-bp dsDNA)" states
   N = 12, and counting the 71-bp molecule as a second length was an
   interpretation error). The current donors are not candidates.
2. Length: a run is compatible when its FASTQ is available and its mean read
   length (base_count / read_count, halved for paired layouts) is within N ± 2.
6. Trajectory: compatible runs are grouped by sample title with the round token
   removed; the trajectory with the most rounds is taken (ties: lowest run
   accession), and it needs at least three rounds with at least 10,000 reads in
   each chosen round (earliest, middle, latest). The read count here is the
   archive's; retained reads are confirmed by the all-read check.
7. Independence from development donors, in two parts (amendment 2, item 8):
   ``step7_biosample`` no shared BioSample, which is mechanical, and
   ``step7_selection`` no shared selection or input library, which a shared
   BioSample does not cover and which comes from the documentation.

Amendment 1 (DESIGN.md, 2026-09-27): the trajectory is what the archive records
or the publication document, not what the parser reads from sample titles. For
every deposit that passes step 2, ``trajectory_documentation.tsv`` records that
review (source, trajectory, rounds, chosen runs, and pass / fail / pending /
inconclusive); it replaces the parser's step-6 result. The runs it names are then
re-checked against the archive as in step 2 and for the 10,000-read floor
(amendment 2, item 8), because substituting them would otherwise skip the checks
the parser's runs had passed.

Usage (network; writes the table):

    python benchmarks/synthetic/screen_donors.py --out benchmarks/synthetic/donor_screening.tsv
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
import time
from collections import defaultdict
from pathlib import Path

from selexprep.catalog.metadata import load_metadata
from selexprep.fetch.plan import build_fetch_plan

HERE = Path(__file__).resolve().parent
INSDC = re.compile(r"^(PRJNA|PRJEB|PRJDB|SRP|ERP|DRP)\d+$")
LENGTH_TOLERANCE = 2
MIN_ROUNDS = 3
MIN_READS = 10_000

COLUMNS = [
    "accession",
    "n_random_text",
    "n_stated",
    "n_source",
    "chemistry",
    "study_title",
    "step1_candidate",
    "runs_total",
    "runs_with_fastq",
    "runs_compatible",
    "step2_length",
    "trajectory",
    "trajectory_rounds",
    "chosen_rounds",
    "chosen_runs",
    "chosen_min_reads",
    "step2_chosen_runs",
    "step6_trajectory",
    "step7_biosample",
    "step7_selection",
    "documentation_source",
    "decision",
    "reason",
]


def stated_n(text: str) -> int | None:
    """The single length in a free-text ``n_random`` value, else None."""
    cleaned = re.sub(r"\d+(\.\d+)?\s*%", "", text or "")
    values = {int(v) for v in re.findall(r"\d+", cleaned)}
    return values.pop() if len(values) == 1 else None


def reviewed_n() -> dict[str, tuple[int | None, str]]:
    """The hand-reviewed length per accession, with its reason (amendment 2, item 9).

    Authoritative where it lists an accession: a length inside a longer molecule
    is a stated length, and a range or a second library is not.
    """
    with (HERE / "n_random_review.tsv").open() as fh:
        return {
            r["accession"]: (int(r["reviewed_n"]) if r["reviewed_n"] else None, r["reason"])
            for r in csv.DictReader(fh, delimiter="\t")
        }


def mean_read_length(run) -> float | None:
    if run.read_count <= 0:
        return None
    return run.base_count / run.read_count / (2 if run.paired_end else 1)


def trajectory_key(title: str, round_number: int) -> str:
    """Sample title with its round token removed (``rd4``, ``Round 04``, ``Cycle4``...).

    The token is the last occurrence of the round number preceded by a round
    word; failing that, the last standalone occurrence of the number. Titles
    carry other digits (``GATA4_5_6_rd5``), so the first occurrence is not it.
    """
    keyword = (
        rf"(?i)[\s_-]*(round|rd|cycle|cyc|(?<![a-z])r|(?<![a-z])c)[\s_-]*0*{round_number}(?!\d)"
    )
    bare = rf"(?<!\d)0*{round_number}(?!\d)"
    for pattern in (keyword, bare):
        matches = list(re.finditer(pattern, title))
        if matches:
            m = matches[-1]
            return (title[: m.start()] + title[m.end() :]).strip(" _-") or "(single trajectory)"
    return title.strip() or "(single trajectory)"


def text(value: object) -> str:
    """A metadata cell as text; missing values (None, NaN) as ''."""
    return value if isinstance(value, str) else ""


def middle(rounds: list[int]) -> int:
    return rounds[(len(rounds) - 1) // 2]


def fetch_plan(accession: str, attempts: int = 3):
    """ENA run records, retried: a transient archive error must not become a decision."""
    for attempt in range(attempts):
        try:
            return build_fetch_plan(accession)
        except Exception:
            if attempt == attempts - 1:
                raise
            time.sleep(5 * (attempt + 1))
    raise AssertionError("unreachable")


def check_chosen_runs(runs: list, srrs: list[str], n: int) -> tuple[str, object]:
    """Step 2 and the 10,000-read floor on the runs actually chosen (amendment 2).

    The documentation names the runs; they still have to be in the archive, have
    a FASTQ, be within N ± 2 and carry enough reads.
    """
    by_srr = {r.srr: r for r in runs}
    problems = []
    counts = []
    for srr in srrs:
        run = by_srr.get(srr)
        if run is None:
            problems.append(f"{srr}: not in the archive record")
            continue
        counts.append(run.read_count)
        if not run.fastq_urls:
            problems.append(f"{srr}: no FASTQ")
        mean = mean_read_length(run)
        if mean is None or abs(mean - n) > LENGTH_TOLERANCE:
            problems.append(f"{srr}: mean length {mean} outside N ± {LENGTH_TOLERANCE}")
        if run.read_count < MIN_READS:
            problems.append(f"{srr}: {run.read_count} reads")
    return ("pass" if not problems else "fail: " + "; ".join(problems)), (
        min(counts) if counts else ""
    )


def screen(accession: str, n: int, dev_samples: set[str], doc: dict | None) -> dict:
    plan = fetch_plan(accession)
    runs = plan.runs
    with_fastq = [r for r in runs if r.fastq_urls]
    compatible = []
    for r in with_fastq:
        mean = mean_read_length(r)
        if mean is not None and abs(mean - n) <= LENGTH_TOLERANCE:
            compatible.append(r)
    row = {
        "runs_total": len(runs),
        "runs_with_fastq": len(with_fastq),
        "runs_compatible": len(compatible),
        "step2_length": "pass" if len(compatible) >= MIN_ROUNDS else "fail",
        "step7_biosample": "pass"
        if not ({r.sample_accession for r in runs} & dev_samples)
        else "fail",
    }
    groups: dict[str, dict[int, list]] = defaultdict(lambda: defaultdict(list))
    for r in compatible:
        rn = r.round_record.round_number
        if rn is None or r.round_record.is_unassigned:
            continue
        groups[trajectory_key(r.sample_title, rn)][rn].append(r)
    if groups:
        name, by_round = sorted(
            groups.items(),
            key=lambda kv: (-len(kv[1]), min(r.srr for rs in kv[1].values() for r in rs)),
        )[0]
        rounds = sorted(by_round)
        chosen = [rounds[0], middle(rounds), rounds[-1]] if len(rounds) >= MIN_ROUNDS else rounds
        chosen_runs = [min(by_round[rn], key=lambda r: r.srr) for rn in chosen]
        min_reads = min(r.read_count for r in chosen_runs)
        row.update(
            trajectory=name,
            trajectory_rounds=",".join(map(str, rounds)),
            chosen_rounds=",".join(map(str, chosen)),
            chosen_runs=",".join(r.srr for r in chosen_runs),
            chosen_min_reads=min_reads,
            step6_trajectory="pass"
            if len(rounds) >= MIN_ROUNDS and min_reads >= MIN_READS
            else "fail",
        )
    else:
        row["step6_trajectory"] = "fail"
    if doc is not None:
        # Amendment 1: the documented trajectory replaces the parser's. Amendment 2:
        # its runs go through the same checks the parser's runs passed.
        srrs = [s for s in doc["chosen_runs"].split(",") if s]
        step2_chosen, min_reads = check_chosen_runs(runs, srrs, n) if srrs else ("", "")
        row.update(
            documentation_source=doc["documentation_source"],
            trajectory=doc["trajectory"],
            trajectory_rounds=doc["rounds_documented"],
            chosen_rounds=doc["chosen_rounds"],
            chosen_runs=doc["chosen_runs"],
            chosen_min_reads=min_reads,
            step2_chosen_runs=step2_chosen,
            step6_trajectory=doc["step6_documented"],
            step7_selection=doc.get("step7_selection", ""),
        )
    return row


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args(argv)

    donors = list(csv.DictReader((HERE / "donors.tsv").open(), delimiter="\t"))
    current = {d["donor"] for d in donors}
    dev_samples: set[str] = set()
    for acc in sorted({d["donor"] for d in donors if d["set"] == "dev"}):
        dev_samples |= {r.sample_accession for r in fetch_plan(acc).runs}

    documented = {
        d["accession"]: d
        for d in csv.DictReader((HERE / "trajectory_documentation.tsv").open(), delimiter="\t")
    }

    review = reviewed_n()
    meta = load_metadata()
    rows = []
    for _, m in meta.sort_values("bioproject_id").iterrows():
        acc = str(m["bioproject_id"])
        if not INSDC.match(acc):
            continue
        n_text = text(m["n_random"])
        if acc in review:
            n, why = review[acc]
            source = f"reviewed by hand (n_random_review.tsv): {why}"
        else:
            n, source = stated_n(n_text), "the single length in the curated field"
        row = dict.fromkeys(COLUMNS, "")
        row.update(
            accession=acc,
            n_random_text=n_text,
            n_stated=n if n is not None else "",
            n_source=source if n is not None else "",
            chemistry=text(m.get("chemistry")),
            study_title=text(m.get("study_title")),
        )
        if acc in current:
            row.update(step1_candidate="no", decision="excluded", reason="current donor")
        elif n is None:
            row.update(
                step1_candidate="no",
                decision="excluded",
                reason="no single stated N: "
                + (review[acc][1] if acc in review else "the curated field states none"),
            )
        else:
            row["step1_candidate"] = "yes"
            doc = documented.get(acc)
            try:
                row.update(screen(acc, n, dev_samples, doc))
            except Exception as exc:  # still failing after retries: flag, never exclude
                row.update(decision="error", reason=f"ENA record unavailable, rerun: {exc}")
            if not row["decision"] and row["step2_length"] == "fail":
                row.update(decision="excluded", reason="failed step2_length")
            elif not row["decision"] and doc is None:
                row.update(
                    decision="pending",
                    reason="trajectory documentation not reviewed yet: the parser's trajectory is "
                    "only a first pass (amendment 1, rule 1)",
                )
            elif not row["decision"] and doc["step6_documented"] in ("pending", "inconclusive"):
                row.update(decision=doc["step6_documented"], reason=doc["reason"])
            elif not row["decision"] and doc["step6_documented"] == "fail":
                row.update(decision="excluded", reason="step 6 (documented): " + doc["reason"])
            if not row["decision"]:
                # Every check must say "pass": an empty or pending value is not one.
                steps = (
                    "step2_chosen_runs",
                    "step6_trajectory",
                    "step7_biosample",
                    "step7_selection",
                )
                failed = [s for s in steps if str(row[s]).startswith("fail")]
                open_ = [s for s in steps if str(row[s]) != "pass" and s not in failed]
                if failed:
                    row.update(decision="excluded", reason="failed " + ", ".join(failed))
                elif open_:
                    row.update(decision="pending", reason="not yet passed: " + ", ".join(open_))
                else:
                    row.update(
                        decision="to all-read check",
                        reason="steps 3-5 pending (reads and documentation)",
                    )
        rows.append(row)
        print(acc, row["decision"], row["reason"], file=sys.stderr)

    with args.out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    passed = [r["accession"] for r in rows if r["decision"] == "to all-read check"]
    errors = [r["accession"] for r in rows if r["decision"] == "error"]
    print(f"{len(rows)} INSDC deposits screened; to all-read check: {passed}")
    if errors:
        print(f"archive errors, table not final - rerun: {errors}")
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
