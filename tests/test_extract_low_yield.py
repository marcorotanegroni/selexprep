"""`extract` reports inputs whose yield collapses (`_low_yield_inputs`).

The primers come from one round, so a run whose reads they do not fit loses
most of its reads while cutadapt and `extract` both finish normally. On the
Tier-1 benchmark PRJNA809588 kept under 1% of seven of its ten runs that way,
visible only in the per-round counts.
"""

from __future__ import annotations

from selexprep.extract.runner import LOW_YIELD_RATIO, _low_yield_inputs
from selexprep.extract.trim import TrimReport, input_fastq_name


def _report(fastq: str, n_in: int, n_out: int, *, paired: bool = False) -> TrimReport:
    cmd = ["cutadapt", "-g", "A...C", "-o", "x.fa"]
    cmd += (
        ["-p", "y.fa", f"/d/{fastq}_1.fastq.gz", f"/d/{fastq}_2.fastq.gz"]
        if paired
        else [f"/d/{fastq}.fastq.gz"]
    )
    return TrimReport(cutadapt_cmd=cmd, n_in=n_in, n_out=n_out, return_code=0)


def test_collapsed_inputs_are_named_with_their_counts():
    reports = [
        _report("SRR617", 1000, 985),
        _report("SRR627", 1000, 4),
        _report("SRR624", 1000, 951),
    ]
    (line,) = _low_yield_inputs(reports)
    assert line.startswith("SRR627.fastq.gz: kept 4 of 1,000 reads (0.4%")
    assert "best input 98.5%" in line


def test_healthy_spread_is_not_reported():
    reports = [_report("A", 1000, 990), _report("B", 1000, 790), _report("C", 1000, 830)]
    assert _low_yield_inputs(reports) == []


def test_threshold_is_relative_to_the_best_input():
    best = _report("BEST", 1000, 800)
    at_ratio = _report("EDGE", 1000, int(800 * LOW_YIELD_RATIO))
    assert _low_yield_inputs([best, at_ratio]) == []


def test_empty_inputs_are_ignored():
    assert _low_yield_inputs([_report("EMPTY", 0, 0)]) == []


def test_paired_calls_are_named_by_r1():
    assert input_fastq_name(_report("S", 10, 5, paired=True).cutadapt_cmd) == "S_1.fastq.gz"
