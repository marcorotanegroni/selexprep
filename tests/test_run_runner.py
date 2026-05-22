"""Tests for ``selexprep.run.runner.run_batch`` (orchestration logic).

Heavy stages (``compute_library_report``, ``run_extract``, ``run_qc``,
``count_fasta``) are mocked so the tests verify orchestration / status
transitions, not pipeline algorithms (those are covered in
``test_detect.py``, ``test_extract_runner.py``, etc.).
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
import requests

from selexprep.extract.runner import ExtractResult
from selexprep.library.report import LibraryReport
from selexprep.qc.runner import QcResult
from selexprep.run import run_batch

# ---------------------------------------------------------------------------
# Synthetic fixtures
# ---------------------------------------------------------------------------


def _mock_response(payload: list[dict] | None, *, status_code: int = 200) -> MagicMock:
    resp = MagicMock(spec=requests.Response)
    resp.status_code = status_code
    resp.json.return_value = payload or []
    if status_code >= 400:
        resp.raise_for_status.side_effect = requests.HTTPError(f"HTTP {status_code}")
    else:
        resp.raise_for_status.return_value = None
    return resp


def _row(srr: str, *, sample_title: str, paired: bool = False) -> dict:
    ftp = (
        f"ftp.ena.example/{srr}_1.fastq.gz;ftp.ena.example/{srr}_2.fastq.gz"
        if paired
        else f"ftp.ena.example/{srr}.fastq.gz"
    )
    md5 = "aaa;bbb" if paired else "aaa"
    sizes = "100;200" if paired else "100"
    return {
        "run_accession": srr,
        "study_accession": "PRJEB1",
        "study_title": "T",
        "library_strategy": "OTHER",
        "library_source": "OTHER",
        "library_name": "",
        "experiment_title": "",
        "sample_title": sample_title,
        "sample_accession": f"SAM_{srr}",
        "read_count": "10",
        "base_count": "100",
        "fastq_md5": md5,
        "fastq_bytes": sizes,
        "fastq_ftp": ftp,
    }


def _write_fastq(path: Path, n_records: int = 800) -> None:
    """Drop a deterministic gzipped FASTQ at `path` (synthetic; primer-irrelevant).

    Sized above ``validate_fastq_gz``'s 1024-byte floor — the Codex-pass-1
    resume oracle rejects smaller files as potentially-corrupt.
    """
    import random as _random

    rng = _random.Random(0xC0DE)
    path.parent.mkdir(parents=True, exist_ok=True)
    body = []
    for i in range(n_records):
        seq = "".join(rng.choice("ACGT") for _ in range(60))
        body.append(f"@r{i}\n{seq}\n+\n{'I' * 60}\n")
    with gzip.open(path, "wb") as fh:
        fh.write("".join(body).encode())


def _stub_dl_single(srr: str, output_dir: Path, **kw) -> bool:
    _write_fastq(output_dir / f"{srr}.fastq.gz")
    return True


def _stub_dl_paired(srr: str, output_dir: Path, **kw) -> bool:
    _write_fastq(output_dir / f"{srr}_1.fastq.gz")
    _write_fastq(output_dir / f"{srr}_2.fastq.gz")
    return True


def _make_lr(**overrides: object) -> LibraryReport:
    base: dict[str, object] = {
        "primer_5p": "GGGAATACG",
        "primer_3p": "ATGCATGC",
        "variants_5p": [],
        "variants_3p": [],
        "known_adapter_hits": {},
        "extraction_mode": "BOTH_PRIMERS_SINGLE_READ",
        "full_insert_recovered": True,
        "read_source": "R1",
        "required_action": "NONE",
        "orientation": "FORWARD",
        "n_length_mode": 24,
        "n_length_distribution": {24: 1000},
        "n_length_confidence": 1.0,
        "match_rate_5p": 0.95,
        "match_rate_3p": 0.92,
        "position_consistency_5p": 0.95,
        "position_consistency_3p": 0.92,
        "read_fraction_used_for_inference": 1.0,
        "sampling_seed": 42,
        "confidence": 0.85,
        "status": "HIGH",
        "failure_reason": None,
    }
    base.update(overrides)
    return LibraryReport(**base)  # type: ignore[arg-type]


def _write_accessions_tsv(path: Path, accessions: list[str]) -> None:
    path.write_text("accession\n" + "\n".join(accessions) + "\n", encoding="utf-8")


def _stub_extract_success(*args, **kwargs) -> ExtractResult:
    """Pretend extract ran: emit manifest.json + per-round extracted.fasta.gz + trim_reports.json."""
    outdir: Path = kwargs.get("outdir") or args[2]
    round_map: dict[str, int] = kwargs["round_map"]
    rounds = sorted(set(round_map.values()))
    for r in rounds:
        round_dir = outdir / f"round_{r:02d}"
        round_dir.mkdir(parents=True, exist_ok=True)
        with gzip.open(round_dir / "extracted.fasta.gz", "wb") as fh:
            fh.write(b">r0\nACGT\n>r1\nACGT\n")
    (outdir / "selexprep_manifest.json").write_text(json.dumps({"stub": True}))
    (outdir / "trim_reports.json").write_text("[]")
    return ExtractResult(skipped_reason=None)


def _stub_qc_success(manifest_path: Path, **kw) -> QcResult:
    qc_dir = manifest_path.parent / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)
    (qc_dir / "flags.yaml").write_text("- name: example_flag\n  severity: info\n")
    return QcResult(flags=[], flags_yaml_path=qc_dir / "flags.yaml", plot_paths=[])


# ---------------------------------------------------------------------------
# TSV validation
# ---------------------------------------------------------------------------


def test_run_batch_refuses_duplicate_accessions(tmp_path: Path) -> None:
    tsv = tmp_path / "accs.tsv"
    _write_accessions_tsv(tsv, ["PRJ_A", "PRJ_A"])
    with pytest.raises(ValueError, match="duplicate"):
        run_batch(tsv, tmp_path / "out")


def test_run_batch_refuses_missing_accession_column(tmp_path: Path) -> None:
    tsv = tmp_path / "accs.tsv"
    tsv.write_text("not_accession\nPRJ_A\n", encoding="utf-8")
    with pytest.raises(ValueError, match="must have an 'accession' column"):
        run_batch(tsv, tmp_path / "out")


# ---------------------------------------------------------------------------
# Fetch refusal paths (user peer-review point 3)
# ---------------------------------------------------------------------------


def test_run_batch_records_fetch_refused_status_for_all_none_runs(tmp_path: Path) -> None:
    rows = [_row("SRR_X", sample_title="nothing")]
    tsv = tmp_path / "accs.tsv"
    _write_accessions_tsv(tsv, ["PRJ_BAD"])

    with patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response(rows)):
        report = run_batch(tsv, tmp_path / "out")

    assert len(report.rows) == 1
    row = report.rows[0]
    assert row.status == "FETCH_REFUSED"
    assert row.last_stage_completed == "fetch"
    assert "NONE-confidence" in row.notes


# ---------------------------------------------------------------------------
# Paired-end handling (user peer-review point 1)
# ---------------------------------------------------------------------------


def test_run_batch_paired_end_threads_R2_through_detect_and_extract(tmp_path: Path) -> None:
    rows = [_row("SRR_PE", sample_title="Round 1", paired=True)]
    tsv = tmp_path / "accs.tsv"
    _write_accessions_tsv(tsv, ["PRJ_PE"])

    captured: dict[str, object] = {}

    def fake_compute(seqs, *, read_source, paired_mate_streams=None, **kw):
        captured["read_source"] = read_source
        captured["paired_mate_streams_keys"] = (
            sorted(paired_mate_streams.keys()) if paired_mate_streams else None
        )
        return _make_lr(read_source=read_source)

    def fake_extract(*args, **kwargs):
        captured["paired_r2_inputs_present"] = kwargs.get("paired_r2_inputs") is not None
        captured["paired_r2_round_keys"] = (
            sorted(kwargs["paired_r2_inputs"].keys()) if kwargs.get("paired_r2_inputs") else None
        )
        return _stub_extract_success(*args, **kwargs)

    with (
        patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response(rows)),
        patch("selexprep.fetch.runner.download_srr", side_effect=_stub_dl_paired),
        patch("selexprep.run.runner.compute_library_report", side_effect=fake_compute),
        patch("selexprep.run.runner.run_extract", side_effect=fake_extract),
        patch("selexprep.run.runner.run_qc", side_effect=_stub_qc_success),
        patch("selexprep.run.runner.count_fasta") as mock_count,
    ):
        report = run_batch(tsv, tmp_path / "out")

    assert report.rows[0].status == "OK"
    assert captured["read_source"] == "R1_AND_R2"
    assert captured["paired_mate_streams_keys"] == [1]
    assert captured["paired_r2_inputs_present"] is True
    assert captured["paired_r2_round_keys"] == [1]
    # count_fasta should have been called for round 01
    assert mock_count.call_count >= 1


# ---------------------------------------------------------------------------
# Split-primer skip (user peer-review point 2)
# ---------------------------------------------------------------------------


def test_run_batch_split_primer_skips_count_and_qc(tmp_path: Path) -> None:
    rows = [_row("SRR_SP", sample_title="Round 0", paired=True)]
    tsv = tmp_path / "accs.tsv"
    _write_accessions_tsv(tsv, ["PRJ_SP"])

    def fake_compute(seqs, *, read_source, paired_mate_streams=None, **kw):
        return _make_lr(
            extraction_mode="PAIRED_END_SPLIT_PRIMERS",
            full_insert_recovered=False,
            required_action="READ_MERGING_RECOMMENDED",
            read_source="R1_AND_R2",
        )

    with (
        patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response(rows)),
        patch("selexprep.fetch.runner.download_srr", side_effect=_stub_dl_paired),
        patch("selexprep.run.runner.compute_library_report", side_effect=fake_compute),
        patch(
            "selexprep.run.runner.run_extract", side_effect=_stub_extract_success
        ) as mock_extract,
        patch("selexprep.run.runner.run_qc") as mock_qc,
        patch("selexprep.run.runner.count_fasta") as mock_count,
    ):
        report = run_batch(tsv, tmp_path / "out")

    row = report.rows[0]
    assert row.status == "SKIPPED_READ_MERGING_RECOMMENDED"
    assert row.required_action == "READ_MERGING_RECOMMENDED"
    # extract still runs (it produces partial_5p_R1 + partial_3p_R2 outputs)
    mock_extract.assert_called_once()
    # count + qc must be skipped
    mock_count.assert_not_called()
    mock_qc.assert_not_called()


# ---------------------------------------------------------------------------
# Resume oracles (user peer-review point 5)
# ---------------------------------------------------------------------------


def test_run_batch_fetch_resume_oracle_triggers_redownload_on_missing_fastq(
    tmp_path: Path,
) -> None:
    """Pre-populate fetch artefacts then delete one FASTQ; resume must re-run fetch."""
    rows = [_row("SRR_R", sample_title="Round 0")]
    tsv = tmp_path / "accs.tsv"
    _write_accessions_tsv(tsv, ["PRJ_R"])
    outdir = tmp_path / "out"

    # First pass: populate everything.
    with (
        patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response(rows)),
        patch("selexprep.fetch.runner.download_srr", side_effect=_stub_dl_single),
        patch("selexprep.run.runner.compute_library_report", return_value=_make_lr()),
        patch("selexprep.run.runner.run_extract", side_effect=_stub_extract_success),
        patch("selexprep.run.runner.run_qc", side_effect=_stub_qc_success),
        patch("selexprep.run.runner.count_fasta"),
    ):
        run_batch(tsv, outdir)

    # Now delete the FASTQ and re-run with resume.
    fastq = outdir / "PRJ_R" / "round_00" / "SRR_R.fastq.gz"
    assert fastq.exists()
    fastq.unlink()

    download_call_count = 0

    def counting_dl(srr: str, output_dir: Path, **kw):
        nonlocal download_call_count
        download_call_count += 1
        return _stub_dl_single(srr, output_dir)

    with (
        patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response(rows)),
        patch("selexprep.fetch.runner.download_srr", side_effect=counting_dl),
        patch("selexprep.run.runner.compute_library_report", return_value=_make_lr()),
        patch("selexprep.run.runner.run_extract", side_effect=_stub_extract_success),
        patch("selexprep.run.runner.run_qc", side_effect=_stub_qc_success),
        patch("selexprep.run.runner.count_fasta"),
    ):
        run_batch(tsv, outdir, resume=True)

    assert download_call_count == 1  # re-fetched the missing run
    assert fastq.exists()


def test_run_batch_resume_skips_completed_stages(tmp_path: Path) -> None:
    """With every sentinel present, no stage's heavy function is invoked."""
    rows = [_row("SRR_S", sample_title="Round 0")]
    tsv = tmp_path / "accs.tsv"
    _write_accessions_tsv(tsv, ["PRJ_S"])
    outdir = tmp_path / "out"

    def fake_count(extracted: Path, parquet: Path, *args, **kw) -> dict:
        # Drop an empty parquet so the resume oracle on pass 2 sees it.
        import pandas as pd

        pd.DataFrame({"sequence": [], "reads": [], "rank": [], "rpm": []}).to_parquet(parquet)
        return {}

    # First pass.
    with (
        patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response(rows)),
        patch("selexprep.fetch.runner.download_srr", side_effect=_stub_dl_single),
        patch("selexprep.run.runner.compute_library_report", return_value=_make_lr()),
        patch("selexprep.run.runner.run_extract", side_effect=_stub_extract_success),
        patch("selexprep.run.runner.run_qc", side_effect=_stub_qc_success),
        patch("selexprep.run.runner.count_fasta", side_effect=fake_count) as count_call,
    ):
        run_batch(tsv, outdir)
        first_count_calls = count_call.call_count

    # Second pass with --resume — heavy stages should not run again.
    with (
        patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response(rows)),
        patch("selexprep.fetch.runner.download_srr") as mock_dl,
        patch("selexprep.run.runner.compute_library_report") as mock_detect,
        patch("selexprep.run.runner.run_extract") as mock_extract,
        patch("selexprep.run.runner.run_qc") as mock_qc,
        patch("selexprep.run.runner.count_fasta") as mock_count_resume,
    ):
        run_batch(tsv, outdir, resume=True)

    mock_dl.assert_not_called()
    mock_detect.assert_not_called()
    mock_extract.assert_not_called()
    mock_qc.assert_not_called()
    mock_count_resume.assert_not_called()
    assert first_count_calls >= 1


# ---------------------------------------------------------------------------
# Error policy
# ---------------------------------------------------------------------------


def test_run_batch_continues_after_failure_by_default(tmp_path: Path) -> None:
    rows_bad = [_row("SRR_NONE", sample_title="nothing")]
    rows_good = [_row("SRR_OK", sample_title="Round 0")]

    def fake_get(url, params=None, **kw):
        acc = params["accession"]
        if acc == "PRJ_BAD":
            return _mock_response(rows_bad)
        return _mock_response(rows_good)

    tsv = tmp_path / "accs.tsv"
    _write_accessions_tsv(tsv, ["PRJ_BAD", "PRJ_GOOD"])

    with (
        patch("selexprep.fetch.inspect.requests.get", side_effect=fake_get),
        patch("selexprep.fetch.runner.download_srr", side_effect=_stub_dl_single),
        patch("selexprep.run.runner.compute_library_report", return_value=_make_lr()),
        patch("selexprep.run.runner.run_extract", side_effect=_stub_extract_success),
        patch("selexprep.run.runner.run_qc", side_effect=_stub_qc_success),
        patch("selexprep.run.runner.count_fasta"),
    ):
        report = run_batch(tsv, tmp_path / "out")

    statuses = {r.accession: r.status for r in report.rows}
    assert statuses == {"PRJ_BAD": "FETCH_REFUSED", "PRJ_GOOD": "OK"}


def test_run_batch_stop_on_error_halts_after_first_failure(tmp_path: Path) -> None:
    rows_bad = [_row("SRR_NONE", sample_title="nothing")]
    rows_good = [_row("SRR_OK", sample_title="Round 0")]

    def fake_get(url, params=None, **kw):
        acc = params["accession"]
        return _mock_response(rows_bad if acc == "PRJ_BAD" else rows_good)

    tsv = tmp_path / "accs.tsv"
    _write_accessions_tsv(tsv, ["PRJ_BAD", "PRJ_GOOD"])

    with (
        patch("selexprep.fetch.inspect.requests.get", side_effect=fake_get),
        patch("selexprep.fetch.runner.download_srr", side_effect=_stub_dl_single),
        patch("selexprep.run.runner.compute_library_report", return_value=_make_lr()),
        patch("selexprep.run.runner.run_extract", side_effect=_stub_extract_success),
        patch("selexprep.run.runner.run_qc", side_effect=_stub_qc_success),
        patch("selexprep.run.runner.count_fasta"),
    ):
        report = run_batch(tsv, tmp_path / "out", stop_on_error=True)

    assert [r.accession for r in report.rows] == ["PRJ_BAD"]


def test_run_batch_detect_failure_records_status(tmp_path: Path) -> None:
    rows = [_row("SRR_X", sample_title="Round 0")]
    tsv = tmp_path / "accs.tsv"
    _write_accessions_tsv(tsv, ["PRJ_DETECT_FAIL"])

    with (
        patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response(rows)),
        patch("selexprep.fetch.runner.download_srr", side_effect=_stub_dl_single),
        patch(
            "selexprep.run.runner.compute_library_report",
            side_effect=RuntimeError("boom"),
        ),
    ):
        report = run_batch(tsv, tmp_path / "out")

    row = report.rows[0]
    assert row.status == "DETECT_FAILED"
    assert "boom" in row.notes


# ---------------------------------------------------------------------------
# Summary writer
# ---------------------------------------------------------------------------


def test_run_batch_emits_run_summary_tsv_sorted_by_accession(tmp_path: Path) -> None:
    rows = [_row("SRR_S", sample_title="Round 0")]
    tsv = tmp_path / "accs.tsv"
    _write_accessions_tsv(tsv, ["PRJ_Z", "PRJ_A"])

    with (
        patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response(rows)),
        patch("selexprep.fetch.runner.download_srr", side_effect=_stub_dl_single),
        patch("selexprep.run.runner.compute_library_report", return_value=_make_lr()),
        patch("selexprep.run.runner.run_extract", side_effect=_stub_extract_success),
        patch("selexprep.run.runner.run_qc", side_effect=_stub_qc_success),
        patch("selexprep.run.runner.count_fasta"),
    ):
        report = run_batch(tsv, tmp_path / "out")

    assert report.summary_tsv is not None
    body = report.summary_tsv.read_text().strip().splitlines()
    header = body[0].split("\t")
    assert header[0] == "accession" and "status" in header
    data_rows = [line.split("\t")[0] for line in body[1:]]
    assert data_rows == ["PRJ_A", "PRJ_Z"]


# ---------------------------------------------------------------------------
# Manual-review separation (user peer-review point 3)
# ---------------------------------------------------------------------------


def test_run_batch_with_manual_review_keeps_rounds_tsv_clean(tmp_path: Path) -> None:
    """When --allow-manual-review is on, NONE-confidence runs land in round_unknown/
    AND in manual_review.tsv, and rounds.tsv stays clean.
    """
    rows = [
        _row("SRR_OK", sample_title="Round 4"),
        _row("SRR_MR", sample_title="nothing"),
    ]
    tsv = tmp_path / "accs.tsv"
    _write_accessions_tsv(tsv, ["PRJ_MIX"])

    with (
        patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response(rows)),
        patch("selexprep.fetch.runner.download_srr", side_effect=_stub_dl_single),
        patch("selexprep.run.runner.compute_library_report", return_value=_make_lr()),
        patch("selexprep.run.runner.run_extract", side_effect=_stub_extract_success),
        patch("selexprep.run.runner.run_qc", side_effect=_stub_qc_success),
        patch("selexprep.run.runner.count_fasta"),
    ):
        report = run_batch(tsv, tmp_path / "out", allow_manual_review=True)

    assert report.rows[0].status == "OK"
    acc_dir = tmp_path / "out" / "PRJ_MIX"
    rounds_body = sorted((acc_dir / "rounds.tsv").read_text().strip().splitlines()[1:])
    assert rounds_body == ["SRR_OK.fastq.gz\t4"]
    assert (acc_dir / "manual_review.tsv").exists()
    assert "SRR_MR" in (acc_dir / "manual_review.tsv").read_text()
    assert (acc_dir / "round_unknown" / "SRR_MR.fastq.gz").exists()


# ===========================================================================
# Codex pass 1 regression tests
# ===========================================================================


def test_run_batch_fetch_http_error_records_FETCH_FAILED_not_extract_failed(
    tmp_path: Path,
) -> None:
    """Codex pass 1 blocking fix: HTTPError from build_fetch_plan must NOT
    surface as EXTRACT_FAILED via the outer except. Stage-classify it as
    FETCH_FAILED inside the fetch block.
    """
    tsv = tmp_path / "accs.tsv"
    _write_accessions_tsv(tsv, ["PRJ_HTTP"])

    err_response = MagicMock(spec=requests.Response)
    err_response.status_code = 503
    err_response.raise_for_status.side_effect = requests.HTTPError("503 ENA down")

    with patch("selexprep.fetch.inspect.requests.get", return_value=err_response):
        report = run_batch(tsv, tmp_path / "out")

    row = report.rows[0]
    assert row.status == "FETCH_FAILED"
    assert row.last_stage_completed == "fetch"
    assert "HTTPError" in row.notes


def test_run_batch_corrupt_fastq_is_redownloaded_by_resume_oracle(
    tmp_path: Path,
) -> None:
    """Codex pass 1 blocking fix: a corrupt .fastq.gz (passes Path.exists()
    but fails gzip -t) must trigger re-download on --resume.
    """
    rows = [_row("SRR_C", sample_title="Round 0")]
    tsv = tmp_path / "accs.tsv"
    _write_accessions_tsv(tsv, ["PRJ_CORRUPT"])
    outdir = tmp_path / "out"

    # First pass: populate.
    with (
        patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response(rows)),
        patch("selexprep.fetch.runner.download_srr", side_effect=_stub_dl_single),
        patch("selexprep.run.runner.compute_library_report", return_value=_make_lr()),
        patch("selexprep.run.runner.run_extract", side_effect=_stub_extract_success),
        patch("selexprep.run.runner.run_qc", side_effect=_stub_qc_success),
        patch("selexprep.run.runner.count_fasta"),
    ):
        run_batch(tsv, outdir)

    # Corrupt the FASTQ: keep the path, garble the gzip stream.
    fastq = outdir / "PRJ_CORRUPT" / "round_00" / "SRR_C.fastq.gz"
    assert fastq.exists()
    fastq.write_bytes(b"this is not a gzip stream at all")

    dl_count = 0

    def counting_dl(srr: str, output_dir: Path, **kw) -> bool:
        nonlocal dl_count
        dl_count += 1
        # Real download_srr would overwrite — mirror that.
        for p in output_dir.glob(f"{srr}*.fastq.gz"):
            p.unlink()
        return _stub_dl_single(srr, output_dir)

    with (
        patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response(rows)),
        patch("selexprep.fetch.runner.download_srr", side_effect=counting_dl),
        patch("selexprep.run.runner.compute_library_report", return_value=_make_lr()),
        patch("selexprep.run.runner.run_extract", side_effect=_stub_extract_success),
        patch("selexprep.run.runner.run_qc", side_effect=_stub_qc_success),
        patch("selexprep.run.runner.count_fasta"),
    ):
        run_batch(tsv, outdir, resume=True)

    # The corrupt file was re-downloaded.
    assert dl_count == 1
    # And the resulting file passes validation.
    from selexprep.fetch.download import validate_fastq_gz

    assert validate_fastq_gz(fastq)


def test_run_batch_count_yaml_flags_uses_yaml_parser(tmp_path: Path) -> None:
    """Codex pass 1 non-blocking fix: a flags.yaml whose first key is
    'evidence' (sort-keys=True) must still be counted correctly on resume.
    """
    from selexprep.run.runner import _count_yaml_flags

    yaml_path = tmp_path / "flags.yaml"
    yaml_path.write_text(
        # safe_dump(sort_keys=True) would emit `- evidence:` first when
        # `evidence` is a non-empty dict (e < n < s alphabetically).
        "- evidence:\n"
        "    a: 1\n"
        "  name: flag_one\n"
        "  severity: warn\n"
        "- evidence:\n"
        "    a: 2\n"
        "  name: flag_two\n"
        "  severity: info\n",
        encoding="utf-8",
    )

    # The old line-startswith implementation would return 0; YAML-parsing
    # returns 2.
    assert _count_yaml_flags(yaml_path) == 2


def test_fastq_filenames_for_run_derives_from_urls() -> None:
    """Codex pass 1 non-blocking fix: filenames should derive from URLs
    (matching ``download_srr_ena_direct`` semantics), not be synthesized
    from the SRR + paired_end flag.
    """
    from selexprep.fetch.metadata import RoundRecord
    from selexprep.fetch.plan import FetchRun, fastq_filenames_for_run

    # Hypothetical alternate naming (e.g., a future ENA URL convention or
    # a Zenodo mirror): the function must respect what the URLs say.
    run = FetchRun(
        srr="SRR_X",
        sample_accession="",
        sample_title="",
        library_name="",
        experiment_title="",
        read_count=0,
        base_count=0,
        fastq_urls=[
            "ftp.ena.example/foo_subdir/SRR_X_unusual.fastq.gz",
        ],
        fastq_md5s=["a"],
        fastq_bytes=[1],
        paired_end=False,
        round_record=RoundRecord(
            srr="SRR_X",
            round_number=0,
            confidence="HIGH",
            source_field="x",
            matched_pattern="x",
        ),
    )
    assert fastq_filenames_for_run(run) == ["SRR_X_unusual.fastq.gz"]


def test_run_batch_unexpected_error_uses_unexpected_failure_status(tmp_path: Path) -> None:
    """Codex pass 1 blocking fix sanity check: an exception that escapes
    every per-stage try/except is classified as UNEXPECTED_FAILURE — never
    EXTRACT_FAILED (which would falsely point at the extract stage).
    """
    tsv = tmp_path / "accs.tsv"
    _write_accessions_tsv(tsv, ["PRJ_BOOM"])

    # Force the OS-level mkdir to raise (a true "unexpected" condition).
    with patch(
        "selexprep.run.runner._process_one_accession",
        side_effect=OSError("filesystem boom"),
    ):
        report = run_batch(tsv, tmp_path / "out")

    row = report.rows[0]
    assert row.status == "UNEXPECTED_FAILURE"
    assert "OSError" in row.notes
