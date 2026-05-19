"""Unit tests for ``selexprep.fetch.inspect`` (mocked ENA REST)."""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
import requests

from selexprep.fetch.inspect import (
    ENA_FILEREPORT_URL,
    InspectReport,
    RunFileInfo,
    inspect_accession,
    write_inspect_json,
)


def _mock_response(payload: list[dict] | None, *, status_code: int = 200) -> MagicMock:
    """Build a ``requests.get`` return value with a canned JSON body."""
    resp = MagicMock(spec=requests.Response)
    resp.status_code = status_code
    resp.json.return_value = payload or []
    if status_code >= 400:
        resp.raise_for_status.side_effect = requests.HTTPError(f"HTTP {status_code}")
    else:
        resp.raise_for_status.return_value = None
    return resp


# ---------------------------------------------------------------------------
# inspect_accession
# ---------------------------------------------------------------------------


def test_inspect_accession_parses_single_run() -> None:
    fake_rows = [
        {
            "run_accession": "SRR000001",
            "study_accession": "PRJEB00001",
            "study_title": "Test SELEX study",
            "library_strategy": "OTHER",
            "library_source": "OTHER",
            "read_count": "1000000",
            "base_count": "75000000",
            "fastq_md5": "abc123",
            "fastq_bytes": "12345",
            "fastq_ftp": "ftp.ena.example/SRR000001.fastq.gz",
        }
    ]
    with patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response(fake_rows)):
        report = inspect_accession("SRR000001")

    assert report.accession == "SRR000001"
    assert report.bioproject_id == "PRJEB00001"
    assert report.study_title == "Test SELEX study"
    assert report.library_strategy == "OTHER"
    assert report.library_source == "OTHER"
    assert len(report.runs) == 1
    run = report.runs[0]
    assert run.run_accession == "SRR000001"
    assert run.read_count == 1000000
    assert run.base_count == 75000000
    assert run.fastq_size_bytes == [12345]
    assert run.fastq_md5 == ["abc123"]


def test_inspect_accession_parses_paired_end_semicolon_lists() -> None:
    fake_rows = [
        {
            "run_accession": "SRR000002",
            "study_accession": "PRJEB00002",
            "study_title": "Paired-end study",
            "library_strategy": "OTHER",
            "library_source": "OTHER",
            "read_count": "500000",
            "base_count": "37500000",
            "fastq_md5": "abc;def",
            "fastq_bytes": "10000;11000",
            "fastq_ftp": "f1.gz;f2.gz",
        }
    ]
    with patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response(fake_rows)):
        report = inspect_accession("SRR000002")

    run = report.runs[0]
    assert run.fastq_size_bytes == [10000, 11000]
    assert run.fastq_md5 == ["abc", "def"]


def test_inspect_accession_parses_multi_run_study() -> None:
    fake_rows = [
        {
            "run_accession": f"SRR00000{i}",
            "study_accession": "PRJEB99",
            "study_title": "Multi-run study",
            "library_strategy": "OTHER",
            "library_source": "OTHER",
            "read_count": str(1000 * i),
            "base_count": str(75000 * i),
            "fastq_md5": f"hash{i}",
            "fastq_bytes": str(1000 * i),
        }
        for i in range(1, 5)
    ]
    with patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response(fake_rows)):
        report = inspect_accession("PRJEB99")

    assert len(report.runs) == 4
    assert report.runs[0].run_accession == "SRR000001"
    assert report.runs[3].run_accession == "SRR000004"
    assert report.bioproject_id == "PRJEB99"


def test_inspect_accession_tolerant_of_missing_fields() -> None:
    """Some SRA submissions omit library_strategy / library_source."""
    fake_rows = [
        {
            "run_accession": "SRR000003",
            "study_accession": "",
            "study_title": "",
            # library_strategy/library_source absent entirely
            "read_count": "0",
            "base_count": "0",
            "fastq_md5": "",
            "fastq_bytes": "",
        }
    ]
    with patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response(fake_rows)):
        report = inspect_accession("SRR000003")

    assert report.library_strategy == ""
    assert report.library_source == ""
    assert report.bioproject_id is None  # empty string -> None
    assert report.runs[0].fastq_size_bytes == []
    assert report.runs[0].fastq_md5 == []


def test_inspect_accession_raises_on_empty_response() -> None:
    with (
        patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response([])),
        pytest.raises(ValueError, match="ENA returned no records"),
    ):
        inspect_accession("SRR_DOES_NOT_EXIST")


def test_inspect_accession_propagates_http_error() -> None:
    with (
        patch(
            "selexprep.fetch.inspect.requests.get",
            return_value=_mock_response(None, status_code=503),
        ),
        pytest.raises(requests.HTTPError),
    ):
        inspect_accession("SRR000001")


def test_inspect_accession_passes_timeout_through() -> None:
    """timeout_s arg is forwarded to requests.get."""
    fake = _mock_response(
        [
            {
                "run_accession": "SRR0",
                "study_accession": "P0",
                "study_title": "",
                "library_strategy": "",
                "library_source": "",
                "read_count": "0",
                "base_count": "0",
                "fastq_md5": "",
                "fastq_bytes": "",
            }
        ]
    )
    with patch("selexprep.fetch.inspect.requests.get", return_value=fake) as mock_get:
        inspect_accession("SRR0", timeout_s=15)
    args, kwargs = mock_get.call_args
    assert args[0] == ENA_FILEREPORT_URL
    assert kwargs["timeout"] == 15


# ---------------------------------------------------------------------------
# write_inspect_json
# ---------------------------------------------------------------------------


def test_write_inspect_json_emits_sorted_keys(tmp_path: Path) -> None:
    report = InspectReport(
        accession="SRR1",
        bioproject_id="PRJ1",
        study_title="t",
        library_strategy="OTHER",
        library_source="OTHER",
        runs=[
            RunFileInfo(
                run_accession="SRR1",
                read_count=100,
                base_count=7500,
                fastq_size_bytes=[1000],
                fastq_md5=["xyz"],
            )
        ],
    )
    out = tmp_path / "inspect.json"
    write_inspect_json(report, out)

    payload = json.loads(out.read_text())
    assert payload["accession"] == "SRR1"
    assert payload["runs"][0]["run_accession"] == "SRR1"
    # Top-level keys are sorted alphabetically.
    keys = list(payload.keys())
    assert keys == sorted(keys)
