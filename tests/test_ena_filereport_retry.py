"""ENA's filereport service can answer HTTP 200 with an error message in place
of the records, intermittently for the same query. ``fetch`` must retry instead
of reading it as a run without FASTQ (TSV) or failing on malformed JSON."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest
import requests

from selexprep.fetch import download, inspect

ERROR_TSV = (
    "run_accession\tfastq_ftp\tfastq_md5\tfastq_bytes\n"
    'ERROR occurred. Not all results may have been written.Query: {"bool":{}}'
)
GOOD_TSV = (
    "run_accession\tfastq_ftp\tfastq_md5\tfastq_bytes\n"
    "SRR1\tftp.example/SRR1_1.fastq.gz;ftp.example/SRR1_2.fastq.gz\tm1;m2\t10;20\n"
)
ERROR_JSON = '[\n{"run_accession":"SRR1","fastq_ftp":"x"}\nERROR occurred. Not all results'


def _resp(text: str, payload=None) -> MagicMock:
    r = MagicMock(spec=requests.Response)
    r.status_code = 200
    r.raise_for_status.return_value = None
    r.text = text
    if payload is None:
        r.json.side_effect = ValueError("Expecting ',' delimiter")
    else:
        r.json.return_value = payload
    return r


@pytest.fixture(autouse=True)
def no_sleep(monkeypatch):
    monkeypatch.setattr(inspect.time, "sleep", lambda s: None)


def test_an_error_answer_is_retried_until_the_records_come():
    answers = [_resp(ERROR_TSV), _resp(ERROR_TSV), _resp(GOOD_TSV)]
    with patch("selexprep.fetch.inspect.requests.get", side_effect=answers) as get:
        text = inspect.get_ena_filereport({"accession": "SRR1"}, timeout_s=5, as_json=False)
    assert text == GOOD_TSV and get.call_count == 3


def test_a_persistent_error_raises_a_service_error_not_no_data():
    with (
        patch("selexprep.fetch.inspect.requests.get", return_value=_resp(ERROR_TSV)) as get,
        pytest.raises(inspect.EnaServiceError, match="archive's service, not the data"),
    ):
        inspect.get_ena_filereport({"accession": "SRR1"}, timeout_s=5, as_json=False)
    assert get.call_count == inspect.FILEREPORT_ATTEMPTS


def test_malformed_json_is_retried():
    rows = [{"run_accession": "SRR1", "fastq_ftp": "x"}]
    answers = [_resp(ERROR_JSON), _resp('[{"run_accession":"SRR1"}]', payload=rows)]
    with patch("selexprep.fetch.inspect.requests.get", side_effect=answers):
        got = inspect.query_ena_filereport("PRJ1", fields="run_accession,fastq_ftp")
    assert got == rows


def test_a_persistent_json_error_is_a_value_error_the_runner_classifies():
    """``run`` turns a ValueError from fetch into FETCH_FAILED with its message."""
    with (
        patch("selexprep.fetch.inspect.requests.get", return_value=_resp(ERROR_JSON)),
        pytest.raises(ValueError, match="filereport service answered with an error"),
    ):
        inspect.query_ena_filereport("PRJ1", fields="run_accession")


def test_download_retries_the_file_list_then_downloads(tmp_path, monkeypatch):
    answers = [_resp(ERROR_TSV), _resp(GOOD_TSV)]
    fetched = []
    monkeypatch.setattr(
        download, "stream_download", lambda url, dest, **kw: fetched.append(dest.name) or True
    )
    with patch("selexprep.fetch.inspect.requests.get", side_effect=answers):
        assert download.download_srr_ena_direct("SRR1", tmp_path)
    assert fetched == ["SRR1_1.fastq.gz", "SRR1_2.fastq.gz"]


def test_download_reports_a_service_error_and_fails_without_downloading(
    tmp_path, monkeypatch, caplog
):
    monkeypatch.setattr(download, "stream_download", lambda *a, **k: pytest.fail("no download"))
    with patch("selexprep.fetch.inspect.requests.get", return_value=_resp(ERROR_TSV)):
        assert not download.download_srr_ena_direct("SRR1", tmp_path)
    assert "archive's service" in caplog.text
    assert "no fastq_ftp" not in caplog.text
