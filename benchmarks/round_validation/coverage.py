"""Round-assignment coverage of the released version (PROTOCOL.md, Amendment 1).

Descriptive and in-sample. For every run of the INSDC deposits in the bundled
catalogue it records the round `fetch` would assign, with its confidence,
source field and notes, and reports:

- coverage: runs that receive a round and runs left unassigned
  (``RoundRecord.is_unassigned``), per deposit and in total; an assignment is
  never counted as correct;
- agreement with the curated maps in ``benchmarks/round_maps/``, on the runs
  both cover that receive a round, with the unassigned runs stated.

Network only (ENA run and sample records); no FASTQ is downloaded. Run it on
the release tag:

    python benchmarks/round_validation/coverage.py --out benchmarks/round_validation/results

Writes ``runs.tsv`` (one row per run), ``deposits.tsv`` (one row per deposit),
``curated_agreement.tsv`` (one row per run a curated map covers) and
``summary.json`` (totals, the commit and version it ran at). Exits 1 if any
deposit's archive records could not be read, so a partial table is never
mistaken for a complete one.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import re
import subprocess
import sys
import time
from pathlib import Path

from selexprep import __version__
from selexprep.catalog.metadata import load_metadata
from selexprep.fetch.plan import build_fetch_plan

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
INSDC = re.compile(r"^(PRJNA|PRJEB|PRJDB|SRP|ERP|DRP)\d+$")
RUN = re.compile(r"^([SDE]RR\d+)")

RUN_COLUMNS = [
    "accession",
    "run",
    "sample_accession",
    "round",
    "assigned",
    "confidence",
    "source_field",
    "matched_pattern",
    "parser_notes",
]


class SampleAttributesUnavailable(RuntimeError):
    """ENA's sample records could not be read, so rounds came from titles only."""


class _AttributeFailures(logging.Handler):
    """Collects the warning `build_fetch_plan` logs when the sample records fail.

    `build_fetch_plan` degrades on purpose: without the sample attributes it
    parses titles and library names only, and says so in a warning. For a
    coverage table that would be a silent partial result, so the warning is
    turned into a failure here.
    """

    def __init__(self) -> None:
        super().__init__(level=logging.WARNING)
        self.messages: list[str] = []

    def emit(self, record: logging.LogRecord) -> None:
        message = record.getMessage()
        if "sample attributes" in message and "could not be fetched" in message:
            self.messages.append(message)


def fetch_plan(accession: str, attempts: int = 3):
    """ENA records, retried: a transient archive error must not become a result.

    A failed request for the sample records counts as a failure too, although
    `build_fetch_plan` recovers from it; a request that succeeds with no round
    attribute does not.
    """
    plan_logger = logging.getLogger("selexprep.fetch.plan")
    for attempt in range(attempts):
        handler = _AttributeFailures()
        plan_logger.addHandler(handler)
        try:
            plan = build_fetch_plan(accession)
            if handler.messages:
                raise SampleAttributesUnavailable(handler.messages[-1])
            return plan
        except Exception:
            if attempt == attempts - 1:
                raise
            time.sleep(5 * (attempt + 1))
        finally:
            plan_logger.removeHandler(handler)
    raise AssertionError("unreachable")


def insdc_deposits() -> list[str]:
    meta = load_metadata()
    return sorted({str(a) for a in meta["bioproject_id"] if INSDC.match(str(a))})


def curated_maps(root: Path) -> dict[str, dict[str, int]]:
    """Accession -> run -> round, from ``<accession>.rounds.tsv``."""
    maps: dict[str, dict[str, int]] = {}
    for path in sorted(root.glob("*.rounds.tsv")):
        accession = path.name.split(".")[0]
        with path.open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                m = RUN.match(Path(row["file"]).name)
                if m:
                    maps.setdefault(accession, {})[m.group(1)] = int(row["round_number"])
    return maps


def run_rows(accession: str) -> list[dict]:
    rows = []
    for r in fetch_plan(accession).runs:
        rec = r.round_record
        assigned = not rec.is_unassigned
        rows.append(
            {
                "accession": accession,
                "run": r.srr,
                "sample_accession": r.sample_accession,
                "round": rec.round_number if assigned else "",
                "assigned": "yes" if assigned else "no",
                "confidence": rec.confidence,
                "source_field": rec.source_field,
                "matched_pattern": rec.matched_pattern,
                "parser_notes": rec.parser_notes,
            }
        )
    return rows


def agreement(
    rows: list[dict], maps: dict[str, dict[str, int]], covered: set[str] | None = None
) -> list[dict]:
    """One row per run a curated map covers: agree, disagree, or unassigned.

    ``covered`` holds the deposits whose archive records were read. Some curated
    maps belong to Tier-1 deposits outside the catalogue (adapter controls, some
    specificity deposits): their runs were never looked up, which is not the
    same as a run missing from its archive record.
    """
    by_run = {(r["accession"], r["run"]): r for r in rows}
    if covered is None:
        covered = {r["accession"] for r in rows}
    out = []
    for accession, runs in sorted(maps.items()):
        for run, curated in sorted(runs.items()):
            row = by_run.get((accession, run))
            if accession not in covered:
                result, parsed = "deposit not in the catalogue", ""
            elif row is None:
                result, parsed = "run not in the archive record", ""
            elif row["assigned"] == "no":
                result, parsed = "unassigned", ""
            else:
                parsed = row["round"]
                result = "agree" if int(parsed) == curated else "disagree"
            out.append(
                {
                    "accession": accession,
                    "run": run,
                    "curated_round": curated,
                    "parsed_round": parsed,
                    "result": result,
                }
            )
    return out


def _write(path: Path, rows: list[dict], columns: list[str]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def commit() -> dict:
    def git(*args: str) -> str:
        return subprocess.run(
            ["git", *args], cwd=REPO, capture_output=True, text=True, check=False
        ).stdout.strip()

    return {
        "commit": git("rev-parse", "HEAD"),
        "tag": git("describe", "--tags", "--exact-match") or None,
        "dirty": bool(git("status", "--porcelain", "--untracked-files=no")),
        "selexprep_version": __version__,
    }


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--maps", type=Path, default=REPO / "benchmarks" / "round_maps")
    p.add_argument(
        "--agreement-only",
        action="store_true",
        help="redo only the comparison with the curated maps, from <out>/runs.tsv (no network)",
    )
    args = p.parse_args(argv)
    if args.agreement_only:
        print(json.dumps(recompute_agreement(args.out, args.maps)["curated_maps"], indent=2))
        return 0
    args.out.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    deposits: list[dict] = []
    failed: list[str] = []
    for accession in insdc_deposits():
        try:
            these = run_rows(accession)
        except Exception as exc:  # still failing after retries: recorded, not dropped
            failed.append(accession)
            deposits.append({"accession": accession, "error": str(exc)})
            print(f"{accession}: archive records unavailable: {exc}", file=sys.stderr)
            continue
        rows.extend(these)
        n_assigned = sum(r["assigned"] == "yes" for r in these)
        deposits.append(
            {
                "accession": accession,
                "runs": len(these),
                "assigned": n_assigned,
                "unassigned": len(these) - n_assigned,
                "error": "",
            }
        )
        print(f"{accession}: {n_assigned}/{len(these)} runs assigned", file=sys.stderr)

    covered = {d["accession"] for d in deposits if not d["error"]}
    agree = agreement(rows, curated_maps(args.maps), covered)
    _write(args.out / "runs.tsv", rows, RUN_COLUMNS)
    _write(
        args.out / "deposits.tsv",
        deposits,
        ["accession", "runs", "assigned", "unassigned", "error"],
    )
    _write(args.out / "curated_agreement.tsv", agree, AGREEMENT_COLUMNS)
    n_assigned = sum(r["assigned"] == "yes" for r in rows)
    summary = {
        **commit(),
        "deposits": len(deposits),
        "deposits_unreadable": failed,
        "runs": len(rows),
        "runs_assigned": n_assigned,
        "runs_unassigned": len(rows) - n_assigned,
        "deposits_all_assigned": sum(
            1 for d in deposits if not d["error"] and d["runs"] and d["unassigned"] == 0
        ),
        "deposits_none_assigned": sum(
            1 for d in deposits if not d["error"] and d["runs"] and d["assigned"] == 0
        ),
        "assigned_by_confidence": {
            c: sum(1 for r in rows if r["assigned"] == "yes" and r["confidence"] == c)
            for c in ("HIGH", "MEDIUM", "LOW")
        },
        "curated_maps": curated_summary(agree),
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 1 if failed else 0


AGREEMENT_COLUMNS = ["accession", "run", "curated_round", "parsed_round", "result"]
OUTSIDE = "deposit not in the catalogue"


def curated_summary(agree: list[dict]) -> dict:
    """Counts over the maps of catalogue deposits; the others are listed apart."""
    inside = [a for a in agree if a["result"] != OUTSIDE]
    counts: dict[str, int] = {}
    for a in inside:
        counts[a["result"]] = counts.get(a["result"], 0) + 1
    return {
        "deposits": len({a["accession"] for a in inside}),
        "runs": len(inside),
        **counts,
        "outside_the_catalogue": sorted({a["accession"] for a in agree if a["result"] == OUTSIDE}),
    }


def recompute_agreement(out: Path, maps: Path) -> dict:
    """Rebuild the agreement from a finished run's ``runs.tsv``, without the network.

    The runs and their rounds stay those of the run that wrote them; only the
    comparison with the curated maps is redone, and ``summary.json`` records
    the commit that redid it.
    """
    with (out / "runs.tsv").open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    with (out / "deposits.tsv").open() as fh:
        covered = {d["accession"] for d in csv.DictReader(fh, delimiter="\t") if not d["error"]}
    agree = agreement(rows, curated_maps(maps), covered)
    _write(out / "curated_agreement.tsv", agree, AGREEMENT_COLUMNS)
    summary = json.loads((out / "summary.json").read_text())
    summary["curated_maps"] = curated_summary(agree)
    summary["curated_maps_recomputed_at"] = commit()
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


if __name__ == "__main__":
    raise SystemExit(main())
