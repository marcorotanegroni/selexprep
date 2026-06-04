"""Unit tests for ``selexprep.fetch.plan`` (mocked ENA REST)."""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
import requests

from selexprep.fetch.inspect import ENA_FILEREPORT_URL
from selexprep.fetch.plan import (
    FetchPlan,
    FetchRun,
    build_fetch_plan,
    fastq_filenames_for_run,
    write_fetch_metadata_json,
)


def _mock_response(payload: list[dict] | None, *, status_code: int = 200) -> MagicMock:
    resp = MagicMock(spec=requests.Response)
    resp.status_code = status_code
    resp.json.return_value = payload or []
    if status_code >= 400:
        resp.raise_for_status.side_effect = requests.HTTPError(f"HTTP {status_code}")
    else:
        resp.raise_for_status.return_value = None
    return resp


def _row(
    srr: str = "SRR000001",
    *,
    study: str = "PRJEB1",
    study_title: str = "test",
    sample_title: str = "",
    library_name: str = "",
    experiment_title: str = "",
    sample_accession: str = "SAMEA1",
    fastq_ftp: str = "ftp.ena.example/SRR000001.fastq.gz",
    fastq_md5: str = "deadbeef",
    fastq_bytes: str = "1000",
    library_strategy: str = "OTHER",
    library_source: str = "OTHER",
) -> dict:
    return {
        "run_accession": srr,
        "study_accession": study,
        "study_title": study_title,
        "library_strategy": library_strategy,
        "library_source": library_source,
        "library_name": library_name,
        "experiment_title": experiment_title,
        "sample_title": sample_title,
        "sample_accession": sample_accession,
        "read_count": "1000",
        "base_count": "75000",
        "fastq_md5": fastq_md5,
        "fastq_bytes": fastq_bytes,
        "fastq_ftp": fastq_ftp,
    }


# ---------------------------------------------------------------------------
# build_fetch_plan
# ---------------------------------------------------------------------------


def test_build_fetch_plan_parses_single_end_run() -> None:
    rows = [_row(srr="SRR1", sample_title="Round 3 sample")]
    with patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response(rows)):
        plan = build_fetch_plan("SRR1")

    assert plan.accession == "SRR1"
    assert plan.bioproject_id == "PRJEB1"
    assert len(plan.runs) == 1
    run = plan.runs[0]
    assert run.srr == "SRR1"
    assert run.paired_end is False
    assert run.round_record.round_number == 3
    assert run.round_record.confidence == "HIGH"
    assert run.round_record.source_field == "sample_title"


def test_build_fetch_plan_parses_paired_end_run() -> None:
    rows = [
        _row(
            srr="SRR2",
            sample_title="cycle 5 PE",
            fastq_ftp="ftp.ena.example/SRR2_1.fastq.gz;ftp.ena.example/SRR2_2.fastq.gz",
            fastq_md5="aaa;bbb",
            fastq_bytes="100;200",
        )
    ]
    with patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response(rows)):
        plan = build_fetch_plan("SRR2")

    run = plan.runs[0]
    assert run.paired_end is True
    assert run.fastq_urls == [
        "ftp.ena.example/SRR2_1.fastq.gz",
        "ftp.ena.example/SRR2_2.fastq.gz",
    ]
    assert run.fastq_md5s == ["aaa", "bbb"]
    assert run.fastq_bytes == [100, 200]
    assert fastq_filenames_for_run(run) == ["SRR2_1.fastq.gz", "SRR2_2.fastq.gz"]
    assert run.round_record.round_number == 5


def test_build_fetch_plan_sorts_runs_by_srr() -> None:
    rows = [
        _row(srr="SRR3", sample_title="Round 2"),
        _row(srr="SRR1", sample_title="Round 0"),
        _row(srr="SRR2", sample_title="Round 1"),
    ]
    with patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response(rows)):
        plan = build_fetch_plan("PRJ")
    assert [r.srr for r in plan.runs] == ["SRR1", "SRR2", "SRR3"]


def test_build_fetch_plan_none_confidence_when_no_round_indicator() -> None:
    rows = [
        _row(
            srr="SRR9",
            sample_title="undifferentiated sample",
            library_name="just words",
            experiment_title="nothing useful here",
        )
    ]
    with patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response(rows)):
        plan = build_fetch_plan("SRR9")

    run = plan.runs[0]
    assert run.round_record.confidence == "NONE"
    assert run.round_record.is_unassigned is True  # round_number is None
    assert run.round_record.round_number is None
    assert plan.has_any_assigned_rounds is False
    assert plan.none_confidence_runs == [run]


def test_build_fetch_plan_has_any_assigned_rounds_true_when_mixed() -> None:
    rows = [
        _row(srr="SRR_OK", sample_title="Round 4"),
        _row(srr="SRR_BAD", sample_title="no marker"),
    ]
    with patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response(rows)):
        plan = build_fetch_plan("PRJX")

    assert plan.has_any_assigned_rounds is True
    none_runs = plan.none_confidence_runs
    assert len(none_runs) == 1
    assert none_runs[0].srr == "SRR_BAD"


def test_build_fetch_plan_medium_single_match_library_name_not_refused() -> None:
    """audit-pilot regression: PRJDB19138 had 5 runs all parsing
    cleanly from library_name (RAPT26-1R / RAPT26-2R / ...). Pre-fix every
    one was flagged as needs_manual_review (because base_confidence=MEDIUM
    on L3 fields), so ``none_confidence_runs`` returned all 5 and the
    fetch runner refused the whole accession.

    Post-fix: a single unambiguous parse from library_name produces a
    MEDIUM record with ``is_unassigned=False``. The runs admit cleanly
    to ``rounds.tsv`` and the accession is fetchable.
    """
    rows = [
        _row(srr="DRR618526", library_name="RAPT26-1R"),
        _row(srr="DRR618527", library_name="RAPT26-2R"),
        _row(srr="DRR618528", library_name="RAPT26-3R"),
        _row(srr="DRR618529", library_name="RAPT26-4R"),
        _row(srr="DRR618530", library_name="RAPT26-5R"),
    ]
    with patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response(rows)):
        plan = build_fetch_plan("PRJDB19138")

    # Every run is MEDIUM single-match — none should be in
    # none_confidence_runs (the refusal trigger).
    for run in plan.runs:
        assert run.round_record.confidence == "MEDIUM"
        assert run.round_record.round_number is not None
        assert run.round_record.is_unassigned is False

    assert plan.has_any_assigned_rounds is True
    assert plan.none_confidence_runs == []  # ← the bug was here


def test_build_fetch_plan_mixed_medium_single_and_genuine_ambiguity() -> None:
    """A plan with both MEDIUM-single-match (should pass) and
    MEDIUM-from-conflict (should be refused) only flags the latter."""
    rows = [
        _row(srr="SRR_OK", library_name="SELEX_round26_N40"),  # MEDIUM, single, OK
        _row(srr="SRR_AMBIG", sample_title="Round 3 / Cycle 7"),  # MEDIUM, conflict, refuse
    ]
    with patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response(rows)):
        plan = build_fetch_plan("PRJX")

    by_srr = {r.srr: r for r in plan.runs}
    assert by_srr["SRR_OK"].round_record.is_unassigned is False
    assert by_srr["SRR_AMBIG"].round_record.is_unassigned is True

    none_runs = plan.none_confidence_runs
    assert [r.srr for r in none_runs] == ["SRR_AMBIG"]


def test_build_fetch_plan_raises_on_empty_response() -> None:
    with (
        patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response([])),
        pytest.raises(ValueError, match="ENA returned no records"),
    ):
        build_fetch_plan("SRR_MISSING")


def test_build_fetch_plan_passes_timeout_through() -> None:
    rows = [_row(srr="SRR1", sample_title="Round 1")]
    fake = _mock_response(rows)
    with patch("selexprep.fetch.inspect.requests.get", return_value=fake) as mock_get:
        build_fetch_plan("SRR1", timeout_s=12)
    args, kwargs = mock_get.call_args
    assert args[0] == ENA_FILEREPORT_URL
    assert kwargs["timeout"] == 12


def test_build_fetch_plan_requests_extended_fields() -> None:
    """sample_title, library_name, experiment_title, sample_accession must be requested."""
    rows = [_row(srr="SRR1", sample_title="Round 1")]
    fake = _mock_response(rows)
    with patch("selexprep.fetch.inspect.requests.get", return_value=fake) as mock_get:
        build_fetch_plan("SRR1")
    _, kwargs = mock_get.call_args
    params = kwargs["params"]
    for field_name in ("sample_title", "library_name", "experiment_title", "sample_accession"):
        assert field_name in params["fields"], f"missing {field_name} in {params['fields']!r}"


# ---------------------------------------------------------------------------
# write_fetch_metadata_json
# ---------------------------------------------------------------------------


def test_write_fetch_metadata_json_sorted_keys(tmp_path: Path) -> None:
    rows = [
        _row(srr="SRR_B", sample_title="Round 1"),
        _row(srr="SRR_A", sample_title="Round 0"),
    ]
    with patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response(rows)):
        plan = build_fetch_plan("PRJ")

    out = tmp_path / "fetch_metadata.json"
    write_fetch_metadata_json(plan, out)

    payload = json.loads(out.read_text())
    assert list(payload.keys()) == sorted(payload.keys())
    # runs are SRR-sorted in the dataclass; preserve that in JSON output
    assert [r["srr"] for r in payload["runs"]] == ["SRR_A", "SRR_B"]


def test_write_fetch_metadata_json_round_record_included(tmp_path: Path) -> None:
    rows = [_row(srr="SRR1", sample_title="cycle 7")]
    with patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response(rows)):
        plan = build_fetch_plan("SRR1")

    out = tmp_path / "fetch_metadata.json"
    write_fetch_metadata_json(plan, out)

    payload = json.loads(out.read_text())
    run = payload["runs"][0]
    assert run["round_record"]["round_number"] == 7
    assert run["round_record"]["confidence"] == "HIGH"


def test_fastq_filenames_for_run_single_vs_paired() -> None:
    """fastq_filenames_for_run is the contract used by the resume oracle."""
    from selexprep.fetch.metadata import RoundRecord

    single = FetchRun(
        srr="SRRX",
        sample_accession="",
        sample_title="",
        library_name="",
        experiment_title="",
        read_count=0,
        base_count=0,
        fastq_urls=["ftp.ena.example/SRRX.fastq.gz"],
        fastq_md5s=["a"],
        fastq_bytes=[1],
        paired_end=False,
        round_record=RoundRecord(
            srr="SRRX",
            round_number=0,
            confidence="HIGH",
            source_field="x",
            matched_pattern="x",
        ),
    )
    paired = FetchRun(
        srr="SRRY",
        sample_accession="",
        sample_title="",
        library_name="",
        experiment_title="",
        read_count=0,
        base_count=0,
        fastq_urls=[
            "ftp.ena.example/SRRY_1.fastq.gz",
            "ftp.ena.example/SRRY_2.fastq.gz",
        ],
        fastq_md5s=["a", "b"],
        fastq_bytes=[1, 2],
        paired_end=True,
        round_record=RoundRecord(
            srr="SRRY",
            round_number=1,
            confidence="HIGH",
            source_field="x",
            matched_pattern="x",
        ),
    )
    assert fastq_filenames_for_run(single) == ["SRRX.fastq.gz"]
    assert fastq_filenames_for_run(paired) == ["SRRY_1.fastq.gz", "SRRY_2.fastq.gz"]


def test_fastq_filenames_for_run_fallback_when_urls_empty() -> None:
    """Defensive: ``fastq_urls`` empty (rare ENA edge case) → synthesize from SRR + paired_end.

    added URL-derived naming, but the synthesizing fallback
    survives so the resume oracle remains useful even when URLs are absent.
    """
    from selexprep.fetch.metadata import RoundRecord

    rr = RoundRecord(
        srr="SRR_FB", round_number=0, confidence="HIGH", source_field="x", matched_pattern="x"
    )
    single = FetchRun(
        srr="SRR_FB",
        sample_accession="",
        sample_title="",
        library_name="",
        experiment_title="",
        read_count=0,
        base_count=0,
        fastq_urls=[],
        fastq_md5s=[],
        fastq_bytes=[],
        paired_end=False,
        round_record=rr,
    )
    paired = FetchRun(
        srr="SRR_FB",
        sample_accession="",
        sample_title="",
        library_name="",
        experiment_title="",
        read_count=0,
        base_count=0,
        fastq_urls=[],
        fastq_md5s=[],
        fastq_bytes=[],
        paired_end=True,
        round_record=rr,
    )
    assert fastq_filenames_for_run(single) == ["SRR_FB.fastq.gz"]
    assert fastq_filenames_for_run(paired) == ["SRR_FB_1.fastq.gz", "SRR_FB_2.fastq.gz"]


def test_fetch_plan_is_frozen() -> None:
    import dataclasses

    plan = FetchPlan(
        accession="X",
        bioproject_id=None,
        study_title="",
        library_strategy="",
        library_source="",
        runs=[],
    )
    with pytest.raises(dataclasses.FrozenInstanceError):
        plan.accession = "Y"  # type: ignore[misc]
