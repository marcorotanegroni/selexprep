"""Evaluator of the semi-synthetic benchmark (DESIGN.md, "Evaluation", and
amendment 2, items 1 and 5).

It reads the generator's truth (``truth.json`` and the read names of the input
FASTQs) and ``selexprep``'s output files (``library_report.json`` and the
extracted FASTA files) only. It imports nothing from ``selexprep``: the report
is read as plain JSON.

Run directory layout, written by ``run_synthetic.py``::

    truth.json  rounds.tsv  round_NN.fastq.gz  run.json
    detect/library_report.json      extract/round_NN/*.fasta.gz
    oracle/library_report.json      oracle/extract/round_NN/*.fasta.gz
"""

from __future__ import annotations

import gzip
import json
import re
from dataclasses import dataclass
from pathlib import Path

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
NAME = re.compile(r"^R(\d+)_(\d+)_s(\d+)_e(\d+)_r([01])([01])([01])_i([01])$")
REFUSED_STATUS = "UNABLE_TO_INFER"
REFUSED_MODE = "UNABLE_TO_EXTRACT"
TWO_SIDED = "BOTH_PRIMERS_SINGLE_READ"


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


@dataclass(frozen=True)
class Read:
    round: int
    seq: str
    start: int
    end: int
    two: bool
    five: bool
    three: bool
    flagged: bool

    def reversed(self) -> Read:
        """The same read as ``extract`` sees it after reverse-complementing it."""
        length = len(self.seq)
        return Read(
            self.round,
            revcomp(self.seq),
            length - self.end,
            length - self.start,
            self.two,
            self.three,
            self.five,
            self.flagged,
        )

    def expected(self, mode: str) -> str | None:
        """What a correct extraction emits for this read in ``mode``."""
        if mode == TWO_SIDED:
            return self.seq[self.start : self.end]
        if mode == "FIVE_PRIME_ONLY":
            return self.seq[self.start :]
        if mode == "THREE_PRIME_ONLY":
            return self.seq[: self.end]
        return None

    def eligible(self, mode: str) -> bool:
        return {
            TWO_SIDED: self.two,
            "FIVE_PRIME_ONLY": self.five,
            "THREE_PRIME_ONLY": self.three,
        }.get(mode, False)


def parse_name(name: str) -> Read | None:
    m = NAME.match(name)
    if not m:
        return None
    r, _, s, e, two, five, three, flag = m.groups()
    return Read(int(r), "", int(s), int(e), two == "1", five == "1", three == "1", flag == "1")


def _open(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def read_inputs(run_dir: Path, truth: dict) -> dict[str, Read]:
    reads: dict[str, Read] = {}
    for rnd in truth["rounds"]:
        with _open(run_dir / rnd["fastq"]) as fh:
            name = None
            for i, line in enumerate(fh):
                if i % 4 == 0:
                    name = line[1:].split()[0]
                elif i % 4 == 1:
                    parsed = parse_name(name or "")
                    if parsed is None:
                        raise ValueError(f"read name without truth: {name}")
                    reads[name] = Read(**{**parsed.__dict__, "seq": line.strip()})  # type: ignore[arg-type]
    return reads


def read_outputs(extract_dir: Path) -> dict[str, str]:
    """Read name -> emitted sequence, over every FASTA of every round."""
    out: dict[str, str] = {}
    if not extract_dir.is_dir():
        return out
    for path in sorted(extract_dir.glob("round_*/*.fasta.gz")):
        with gzip.open(path, "rt") as fh:
            name = None
            for line in fh:
                if line.startswith(">"):
                    name = line[1:].split()[0]
                elif name is not None:
                    out[name] = out.get(name, "") + line.strip()
    return out


def _ratio(num: int, den: int) -> float | None:
    """Zero denominators are not applicable, with or without output."""
    return round(num / den, 6) if den else None


def per_read(mode: str | None, reads: dict[str, Read], outputs: dict[str, str]) -> dict:
    """Per-read metrics for one extraction (DESIGN.md level 2 and 3; amendment 2, item 1).

    Reads whose indel placement touches a boundary of N are kept out and counted
    apart (amendment 2, item 5). The numerator of every recall lies inside its
    denominator.
    """
    emitted = correct = eligible = eligible_correct = n_present = n_recovered = 0
    flagged = flagged_emitted = 0
    per_round_in: dict[int, int] = {}
    per_round_out: dict[int, int] = {}
    for name, read in reads.items():
        per_round_in[read.round] = per_round_in.get(read.round, 0) + 1
        got = outputs.get(name)
        if got is not None:
            per_round_out[read.round] = per_round_out.get(read.round, 0) + 1
        if read.flagged:
            flagged += 1
            flagged_emitted += got is not None
            continue
        expected = read.expected(mode) if mode else None
        ok = got is not None and expected is not None and got == expected
        emitted += got is not None
        correct += ok
        if mode and read.eligible(mode):
            eligible += 1
            eligible_correct += ok
        if read.end > read.start:
            n_present += 1
            n_recovered += got is not None and got == read.seq[read.start : read.end]
    yields = [per_round_out.get(r, 0) / n for r, n in sorted(per_round_in.items())]
    return {
        "emitted": emitted,
        "precision": _ratio(correct, emitted),
        "recall": _ratio(eligible_correct, eligible),
        "eligible": eligible,
        "n_recovery": _ratio(n_recovered, n_present),
        "lowest_round_yield": round(min(yields), 6) if yields else None,
        "flagged_reads": flagged,
        "flagged_emitted": flagged_emitted,
    }


def call_correct(call: str | None, side: str, truth: dict, reverse: bool) -> bool:
    """A 5' call is a suffix of the true sequence before N (outer + constant) in
    every read of the earliest round provided; a 3' call, a prefix of the one
    after N. Error-free and untruncated, in read orientation."""
    if not call:
        return False
    earliest = truth["rounds"][0]
    before, after = earliest["pre"], earliest["post"]
    if reverse:
        before, after = [revcomp(s) for s in after], [revcomp(s) for s in before]
    if side == "5p":
        return all(s.endswith(call) for s in before)
    return all(s.startswith(call) for s in after)


def category(call: str | None, constant: str | None, side: str, correct: bool) -> str:
    """The called constant relative to the true one."""
    if not call:
        return "none"
    if constant and call == constant:
        return "exact"
    if correct and constant and len(call) > len(constant):
        return "includes outer material"
    if constant and (constant.endswith(call) if side == "5p" else constant.startswith(call)):
        return "shorter"
    return "other"


def outcome(truth: dict, run: dict, report: dict | None, sides_ok: dict[str, bool]) -> str:
    """One outcome per run, the first that applies (DESIGN.md, "Evaluation")."""
    if run.get("error"):
        return "error"
    assert report is not None
    refused = report["status"] == REFUSED_STATUS or report["extraction_mode"] == REFUSED_MODE
    if truth["negative_control"]:
        return "negative control, refused" if refused else "negative control, wrong"
    if refused:
        return "refused"
    mode = report["extraction_mode"]
    used = {
        TWO_SIDED: ("5p", "3p"),
        "FIVE_PRIME_ONLY": ("5p",),
        "THREE_PRIME_ONLY": ("3p",),
    }.get(mode, ("5p", "3p"))
    if not all(sides_ok[s] for s in used):
        return "wrong"
    return "complete, correct" if mode == TWO_SIDED else "partial, correct and declared"


def _load(path: Path) -> dict | None:
    return json.loads(path.read_text()) if path.is_file() else None


def evaluate(run_dir: Path) -> dict:
    truth = json.loads((run_dir / "truth.json").read_text())
    run = _load(run_dir / "run.json") or {}
    report = _load(run_dir / "detect" / "library_report.json")
    reads = read_inputs(run_dir, truth)
    oracle = truth["oracle"]
    row: dict = {
        "cell": truth["cell"],
        "config": truth["config"],
        "donor": truth["donor"],
        "family": truth["family"] or "",
        "seed": truth["seed"],
        "error": run.get("error", ""),
    }
    sides_ok = {"5p": False, "3p": False}
    if report is not None:
        reverse = report["orientation"] == "REVERSE"
        # The category compares the call with the constant as configured (for T4,
        # the full 3' constant of the later rounds), not with the oracle's, which
        # leaves out a round-specific base; whether the call equals the oracle's
        # constant is recorded apart.
        configured = truth["constants"]
        for side in ("5p", "3p"):
            call = report[f"primer_{side}"]
            sides_ok[side] = call_correct(call, side, truth, reverse)
            constant, target = configured[f"primer_{side}"], oracle[f"primer_{side}"]
            if reverse:
                # In the frame extract works in: the other side's constant, reversed.
                key = "primer_3p" if side == "5p" else "primer_5p"
                constant = revcomp(configured[key]) if configured[key] else None
                target = revcomp(oracle[key]) if oracle[key] else None
            row[f"call_{side}"] = call or ""
            row[f"call_{side}_correct"] = sides_ok[side]
            row[f"call_{side}_category"] = category(call, constant, side, sides_ok[side])
            row[f"call_{side}_is_oracle"] = bool(call) and call == target
        row.update(
            status=report["status"],
            confidence=report["confidence"],
            extraction_mode=report["extraction_mode"],
            orientation=report["orientation"],
        )
        framed = {k: v.reversed() for k, v in reads.items()} if reverse else reads
        refused = report["status"] == REFUSED_STATUS or report["extraction_mode"] == REFUSED_MODE
        if not refused and not run.get("error"):
            metrics = per_read(report["extraction_mode"], framed, read_outputs(run_dir / "extract"))
            row.update({f"inferred_{k}": v for k, v in metrics.items()})
    if "inferred_emitted" not in row:
        # A refusal or an error emits nothing: precision and recall are not
        # applicable (no mode), but the recovery of N is known and is zero, so a
        # summary over runs counts the loss instead of skipping it.
        metrics = per_read(None, reads, {})
        row.update({f"inferred_{k}": v for k, v in metrics.items()})
    row["outcome"] = outcome(truth, run, report, sides_ok)
    if oracle["mode"] is not None and not run.get("oracle_error"):
        metrics = per_read(oracle["mode"], reads, read_outputs(run_dir / "oracle" / "extract"))
        row.update({f"oracle_{k}": v for k, v in metrics.items()})
    row["oracle_error"] = run.get("oracle_error", "")
    return row


if __name__ == "__main__":
    import sys

    for d in sys.argv[1:]:
        print(json.dumps(evaluate(Path(d)), indent=2))
