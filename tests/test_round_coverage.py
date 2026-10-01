"""Round-coverage script (benchmarks/round_validation, PROTOCOL.md Amendment 1).

No network: the archive rows are built by hand. The script reports coverage and
agreement with the curated maps; it never counts an assignment as correct.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parent.parent / "benchmarks" / "round_validation" / "coverage.py"


@pytest.fixture(scope="module")
def cov():
    spec = importlib.util.spec_from_file_location("round_coverage", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _row(accession, run, rnd):
    return {
        "accession": accession,
        "run": run,
        "round": rnd if rnd is not None else "",
        "assigned": "yes" if rnd is not None else "no",
    }


def test_curated_maps_are_keyed_by_run_accession(cov, tmp_path):
    (tmp_path / "PRJEB1.rounds.tsv").write_text(
        "file\tround_number\nERR1.fastq.gz\t0\nERR2_1.fastq.gz\t3\n"
    )
    assert cov.curated_maps(tmp_path) == {"PRJEB1": {"ERR1": 0, "ERR2": 3}}


def test_agreement_separates_agree_disagree_unassigned_and_missing(cov):
    rows = [_row("PRJEB1", "ERR1", 0), _row("PRJEB1", "ERR2", 4), _row("PRJEB1", "ERR3", None)]
    maps = {"PRJEB1": {"ERR1": 0, "ERR2": 3, "ERR3": 1, "ERR9": 2}}
    results = {a["run"]: a["result"] for a in cov.agreement(rows, maps)}
    assert results == {
        "ERR1": "agree",
        "ERR2": "disagree",
        "ERR3": "unassigned",
        "ERR9": "run not in the archive record",
    }


def test_the_catalogue_has_the_127_insdc_deposits_of_the_protocol(cov):
    assert len(cov.insdc_deposits()) == 127


def test_a_failed_sample_record_request_is_a_failure_not_a_result(cov, monkeypatch):
    """build_fetch_plan recovers from it (titles only) with a warning; the coverage
    table must not, so the warning triggers a retry and, if it persists, an error."""
    import logging

    calls = []

    def degraded(accession):
        calls.append(accession)
        logging.getLogger("selexprep.fetch.plan").warning(
            "build_fetch_plan: sample attributes of %s could not be fetched (%s); "
            "rounds are parsed from titles and library names only",
            accession,
            "timeout",
        )
        return "plan"

    monkeypatch.setattr(cov, "build_fetch_plan", degraded)
    monkeypatch.setattr(cov.time, "sleep", lambda s: None)
    with pytest.raises(cov.SampleAttributesUnavailable):
        cov.fetch_plan("PRJEB1", attempts=3)
    assert len(calls) == 3

    # A request that succeeds, attributes or not, is a result.
    monkeypatch.setattr(cov, "build_fetch_plan", lambda accession: "plan")
    assert cov.fetch_plan("PRJEB1") == "plan"


def test_a_transient_failure_then_success_is_retried(cov, monkeypatch):
    import logging

    state = {"n": 0}

    def flaky(accession):
        state["n"] += 1
        if state["n"] == 1:
            logging.getLogger("selexprep.fetch.plan").warning(
                "build_fetch_plan: sample attributes of %s could not be fetched (x)", accession
            )
        return "plan"

    monkeypatch.setattr(cov, "build_fetch_plan", flaky)
    monkeypatch.setattr(cov.time, "sleep", lambda s: None)
    assert cov.fetch_plan("PRJEB1") == "plan" and state["n"] == 2


def test_maps_of_deposits_outside_the_catalogue_are_counted_apart(cov):
    """Tier-1 maps of adapter controls were never looked up: not a missing run."""
    rows = [_row("PRJEB1", "ERR1", 0)]
    maps = {"PRJEB1": {"ERR1": 0}, "PRJNA9": {"SRR5": 1}}
    agree = cov.agreement(rows, maps, covered={"PRJEB1"})
    assert {a["run"]: a["result"] for a in agree} == {
        "ERR1": "agree",
        "SRR5": "deposit not in the catalogue",
    }
    summary = cov.curated_summary(agree)
    assert summary["deposits"] == 1 and summary["runs"] == 1 and summary["agree"] == 1
    assert summary["outside_the_catalogue"] == ["PRJNA9"]
