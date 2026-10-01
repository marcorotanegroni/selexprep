"""Supplementary table of the parameters of `detect` and `extract`.

Every value is read from the installed package, never typed here, so the table
cannot drift from the version it documents; only the descriptions are written
by hand. A test checks that every public constant of the modules below is in
the table or excluded on purpose (``EXCLUDED``), so a new constant cannot be
left out unnoticed.

    python benchmarks/parameters/parameter_table.py --out benchmarks/parameters

Writes ``parameters.tsv`` and ``parameters.md``, with the package version.
"""

from __future__ import annotations

import argparse
import csv
from importlib.metadata import version
from pathlib import Path

from selexprep import __version__
from selexprep.extract import runner as extract_runner
from selexprep.library import adapters, detect

MODULES = {
    "selexprep.library.detect": detect,
    "selexprep.library.adapters": adapters,
    "selexprep.extract.runner": extract_runner,
}

# (stage, module, constant, what it controls)
PARAMETERS: list[tuple[str, str, str, str]] = [
    # --- Input to inference
    (
        "detect: input",
        "selexprep.library.detect",
        "DEFAULT_MIN_SEQS_FOR_DETECTION",
        "Fewest sequences in the earliest round for detect to infer anything; below it, it refuses.",
    ),
    (
        "detect: input",
        "selexprep.library.detect",
        "DEFAULT_TOP_N",
        "Sequences used for inference: None means every unique sequence, no subsampling.",
    ),
    # --- Flank calling (positional consensus from the read edge)
    (
        "detect: flank boundary",
        "selexprep.library.detect",
        "DEFAULT_MIN_LEN",
        "Shortest flank called as a library constant (nt).",
    ),
    (
        "detect: flank boundary",
        "selexprep.library.detect",
        "DEFAULT_MAX_LEN",
        "Longest flank examined from each read edge (nt).",
    ),
    (
        "detect: flank boundary",
        "selexprep.library.detect",
        "DEFAULT_CONFIDENCE",
        "Positional support for a 'strong' position of the flank.",
    ),
    (
        "detect: flank boundary",
        "selexprep.library.detect",
        "BOUNDARY_SUPPORT_FLOOR",
        "Support below which a position is no longer constant; the boundary stops at the first "
        "run of such positions.",
    ),
    (
        "detect: flank boundary",
        "selexprep.library.detect",
        "BOUNDARY_DROP_RUN",
        "Consecutive positions that must fall below a threshold to end the flank, so one noisy "
        "base does not.",
    ),
    (
        "detect: flank boundary",
        "selexprep.library.detect",
        "BOUNDARY_STRONG_POSITION_FRACTION",
        "Fraction of the called flank's positions that must be strong.",
    ),
    (
        "detect: flank boundary",
        "selexprep.library.detect",
        "BOUNDARY_HIGH_SUPPORT_BASELINE",
        "Drop test: median support of the preceding positions for the test to apply.",
    ),
    (
        "detect: flank boundary",
        "selexprep.library.detect",
        "BOUNDARY_HIGH_SUPPORT_DROP",
        "Drop test: fall in support from that median that marks the boundary.",
    ),
    (
        "detect: flank boundary",
        "selexprep.library.detect",
        "BOUNDARY_HIGH_SUPPORT_POST_MAX",
        "Drop test: support the following positions must stay below. A sequence that makes up "
        "this much of every round, the earliest included, can pass for constant (Known limits).",
    ),
    # --- Outer-edge core rescue
    (
        "detect: core rescue",
        "selexprep.library.detect",
        "CORE_MIN_SUPPORT",
        "Support of the high-support core used when the whole flank's match rate fails.",
    ),
    (
        "detect: core rescue",
        "selexprep.library.detect",
        "CORE_MIN_LEN",
        "Shortest core that can be rescued (nt).",
    ),
    (
        "detect: core rescue",
        "selexprep.library.detect",
        "CORE_STABLE_RUN",
        "Consecutive well-supported positions that end the outer-edge scan.",
    ),
    # --- Library constant shared by every round
    (
        "detect: shared constant",
        "selexprep.library.detect",
        "CONSTANT_CORE_LEN",
        "Inner bases of the flank, next to the random region, on which every round is aligned; "
        "never trimmed (nt).",
    ),
    (
        "detect: shared constant",
        "selexprep.library.detect",
        "CONSTANT_SEARCH_SLACK",
        "How far beyond the earliest round's flank the core may sit in another round (nt).",
    ),
    (
        "detect: shared constant",
        "selexprep.library.detect",
        "CONSTANT_MIN_ANCHORED_FRACTION",
        "Fraction of a round's reads that must carry the core for the round to vote.",
    ),
    (
        "detect: shared constant",
        "selexprep.library.detect",
        "CONSTANT_ANCHOR_AGREEMENT",
        "Agreement among single-copy reads needed to anchor reads that carry the core twice.",
    ),
    # --- Classification and status
    (
        "detect: classification",
        "selexprep.library.detect",
        "PRIMER_FOUND_MATCH_RATE_THRESHOLD",
        "Fraction of reads matching a flank for that side to count as found.",
    ),
    (
        "detect: classification",
        "selexprep.library.detect",
        "UNABLE_TO_EXTRACT_MATCH_RATE",
        "Match rate below which, on both sides, extraction is refused.",
    ),
    (
        "detect: classification",
        "selexprep.library.detect",
        "N_LENGTH_CONFIDENT_FRACTION",
        "Share of reads at the modal random-region length for a one-sided extraction to be allowed.",
    ),
    (
        "detect: classification",
        "selexprep.library.detect",
        "POSITION_CONSISTENCY_TOLERANCE",
        "Offset of a flank, in nt, still counted as the same position.",
    ),
    (
        "detect: status",
        "selexprep.library.detect",
        "STATUS_HIGH_CUTOFF",
        "Composite confidence for HIGH. Rates support for the inferred primers, not a verified "
        "boundary.",
    ),
    (
        "detect: status",
        "selexprep.library.detect",
        "STATUS_MEDIUM_CUTOFF",
        "Composite confidence for MEDIUM; a single round caps the status at MEDIUM.",
    ),
    (
        "detect: status",
        "selexprep.library.detect",
        "STATUS_LOW_CUTOFF",
        "Composite confidence for LOW; below it, UNABLE_TO_INFER.",
    ),
    (
        "detect: status",
        "selexprep.library.detect",
        "COMPOSITE_WEIGHTS",
        "Weights of the signals in the composite confidence, with a round map.",
    ),
    (
        "detect: status",
        "selexprep.library.detect",
        "COMPOSITE_WEIGHTS_NO_ROUND_MAP",
        "The same weights without a round map (no cross-round persistence).",
    ),
    # --- Orientation, variants, adapters
    (
        "detect: orientation",
        "selexprep.library.detect",
        "ORIENTATION_REVERSED_FORWARD_MAX",
        "Share of reverse reads below which the deposit is FORWARD.",
    ),
    (
        "detect: orientation",
        "selexprep.library.detect",
        "ORIENTATION_REVERSED_REVERSE_MIN",
        "Share of reverse reads above which it is REVERSE (between the two: MIXED).",
    ),
    (
        "detect: report",
        "selexprep.library.detect",
        "VARIANTS_TOP_K",
        "Flank variants listed in the report, per side.",
    ),
    (
        "detect: adapters",
        "selexprep.library.adapters",
        "ADAPTER_PROBE_K",
        "Leading bases of a known adapter that mark a call or a read as adapter (nt).",
    ),
    (
        "detect: adapters",
        "selexprep.library.adapters",
        "KNOWN_ADAPTERS",
        "Sequencing adapters excluded from primer candidates (names listed).",
    ),
    # --- Extraction
    (
        "extract: yield warning",
        "selexprep.extract.runner",
        "LOW_YIELD_RATIO",
        "An input keeping less than this fraction of the best input's yield is reported.",
    ),
    (
        "extract: yield warning",
        "selexprep.extract.runner",
        "LOW_YIELD_FLOOR",
        "An input keeping less than this fraction of its own reads is reported.",
    ),
]

# Public constants of MODULES that are not parameters of the method.
EXCLUDED = {
    "KNOWN_ADAPTERS_RC": "reverse complements of KNOWN_ADAPTERS, derived from it",
}

# extract passes no error rate or overlap to cutadapt: its defaults apply.
CUTADAPT_DEFAULTS = [
    ("--error-rate", "0.1", "Maximum error rate in a primer match (cutadapt default)."),
    ("--overlap", "3", "Shortest primer overlap at a read end (cutadapt default)."),
]


def _format(value: object) -> str:
    if isinstance(value, dict):
        if all(isinstance(v, str) for v in value.values()):
            return ", ".join(sorted(value))
        return ", ".join(f"{k} {v}" for k, v in value.items())
    return "None (all)" if value is None else str(value)


def rows() -> list[dict[str, str]]:
    out = []
    for stage, module, name, description in PARAMETERS:
        out.append(
            {
                "stage": stage,
                "parameter": name,
                "value": _format(getattr(MODULES[module], name)),
                "description": description,
                "source": module,
            }
        )
    cutadapt = version("cutadapt")
    for name, value, description in CUTADAPT_DEFAULTS:
        out.append(
            {
                "stage": "extract: cutadapt",
                "parameter": name,
                "value": value,
                "description": description,
                "source": f"cutadapt {cutadapt}",
            }
        )
    return out


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    table = rows()
    columns = ["stage", "parameter", "value", "description", "source"]
    with (args.out / "parameters.tsv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(table)
    lines = [
        f"# Parameters of `detect` and `extract` (selexprep {__version__})",
        "",
        "Values read from the installed package by `benchmarks/parameters/parameter_table.py`.",
        "",
        "| Stage | Parameter | Value | What it controls |",
        "|---|---|---|---|",
    ]
    lines += [
        f"| {r['stage']} | `{r['parameter']}` | {r['value']} | {r['description']} |" for r in table
    ]
    (args.out / "parameters.md").write_text("\n".join(lines) + "\n")
    print(f"{len(table)} parameters -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
