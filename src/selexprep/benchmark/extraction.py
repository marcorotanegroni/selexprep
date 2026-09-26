"""Per-run extraction yield for the recovery arm of the Tier-1 benchmark.

The primer-recovery scorecard (:mod:`selexprep.benchmark.metrics`) scores what
``detect`` reports, and ``detect`` infers the constants from the earliest
round only. Whether those constants then extract the random region from the
*other* rounds is a separate question that the scorecard cannot see: a deposit
whose runs each carry their own inline tag outside the library constant
(PRJNA809588) scores a correct random-region length while ``extract`` keeps
almost nothing from the runs whose tag differs.

This module answers that question from what ``selexprep extract`` actually
wrote. For every recovery-arm deposit it reads ``extract/trim_reports.json``
(one entry per input FASTQ, with cutadapt's reads in / reads kept) and the
per-round extracted FASTA, and emits one row per input FASTQ:

``accession, extraction_mode, round, fastq, n_in, n_out, frac_kept,
output_len_mode, output_frac_len_mode, n_length_truth, note``

``output_len_mode`` is the modal length of the round's extracted sequences;
for a two-sided extraction it should equal ``n_length_truth``. Deposits whose
``extract`` refused (exit code 2, e.g. a report with no usable primers) get a
single row carrying the reason. Nothing is scored here; the table shows the
yield so that a collapse in any run is visible.
"""

from __future__ import annotations

import csv
import gzip
import json
from collections import Counter
from pathlib import Path

from selexprep.benchmark.metrics import load_ground_truth

COLUMNS = [
    "accession",
    "extraction_mode",
    "round",
    "fastq",
    "n_in",
    "n_out",
    "frac_kept",
    "output_len_mode",
    "output_frac_len_mode",
    "n_length_truth",
    "note",
]


def _length_mode(fasta_gz: Path) -> tuple[int | str, str]:
    """Modal sequence length of a gzipped FASTA and its share of the records."""
    lengths: Counter[int] = Counter()
    with gzip.open(fasta_gz, "rt") as fh:
        for line in fh:
            if not line.startswith(">"):
                lengths[len(line.rstrip("\n"))] += 1
    if not lengths:
        return "", ""
    mode, count = lengths.most_common(1)[0]
    return mode, f"{count / sum(lengths.values()):.4f}"


def _input_fastq(cutadapt_cmd: list[str]) -> str:
    """The R1 / single-end input of one cutadapt call (paired calls end in R1, R2)."""
    return Path(cutadapt_cmd[-2] if "-p" in cutadapt_cmd else cutadapt_cmd[-1]).name


def _last_line(path: Path) -> str:
    if not path.exists():
        return ""
    lines = [ln.strip() for ln in path.read_text(errors="replace").splitlines() if ln.strip()]
    return lines[-1] if lines else ""


def extraction_rows(results_dir: Path, accession: str, n_length_truth: int) -> list[dict]:
    """Rows of the yield table for one accession (see module docstring)."""
    acc_dir = results_dir / accession
    extract_dir = acc_dir / "extract"
    report_path = acc_dir / "library_report.json"
    mode = json.loads(report_path.read_text())["extraction_mode"] if report_path.exists() else ""
    base = {"accession": accession, "extraction_mode": mode, "n_length_truth": n_length_truth}

    rc_path = extract_dir / "extract.rc"
    if not rc_path.exists():
        return [{**base, "note": "extract not run"}]
    rc = rc_path.read_text().strip()
    trim_reports = extract_dir / "trim_reports.json"
    if rc != "0" or not trim_reports.exists():
        reason = _last_line(extract_dir / "extract.log") or "no trim_reports.json"
        return [{**base, "note": f"extract exit {rc}: {reason}"}]

    length_cache: dict[str, tuple[int | str, str]] = {}
    rows = []
    for entry in json.loads(trim_reports.read_text()):
        output = Path(entry["output_paths"][0])
        key = str(output)
        if key not in length_cache:
            length_cache[key] = _length_mode(output) if output.exists() else ("", "")
        len_mode, len_frac = length_cache[key]
        n_in, n_out = int(entry["n_in"]), int(entry["n_out"])
        rows.append(
            {
                **base,
                "round": output.parent.name.removeprefix("round_"),
                "fastq": _input_fastq(entry["cutadapt_cmd"]),
                "n_in": n_in,
                "n_out": n_out,
                "frac_kept": f"{n_out / n_in:.4f}" if n_in else "",
                "output_len_mode": len_mode,
                "output_frac_len_mode": len_frac,
                "note": "" if int(entry.get("return_code", 0)) == 0 else "cutadapt non-zero exit",
            }
        )
    return sorted(rows, key=lambda r: (r["round"], r["fastq"]))


def main(argv: list[str] | None = None) -> int:
    import argparse

    p = argparse.ArgumentParser(description="Per-run extraction yield of the recovery arm.")
    p.add_argument("--ground-truth", required=True, type=Path)
    p.add_argument("--results-dir", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    args = p.parse_args(argv)

    recovery = [
        row
        for row in load_ground_truth(args.ground_truth)
        if row.verified and row.read_state == "raw_standard"
    ]
    rows: list[dict] = []
    for row in sorted(recovery, key=lambda r: r.accession):
        acc_rows = extraction_rows(args.results_dir, row.accession, row.n_length_truth)
        rows.extend(acc_rows)
        kept = [float(r["frac_kept"]) for r in acc_rows if r.get("frac_kept")]
        if kept:
            low = sum(k < 0.5 for k in kept)
            print(
                f"{row.accession}: {len(kept)} FASTQs, kept {min(kept):.1%}-{max(kept):.1%} "
                f"per FASTQ, {low} below 50%"
            )
        else:
            print(f"{row.accession}: {acc_rows[0]['note']}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n", restval=""
        )
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
