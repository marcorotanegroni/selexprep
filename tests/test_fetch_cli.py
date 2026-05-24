"""Tests for selexprep.fetch.runner (orchestrator + ``selexprep fetch`` CLI).

The download is monkeypatched to drop deterministic fixture FASTQs so the
test suite has no network dependency.
"""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
import requests
from typer.testing import CliRunner

from selexprep.cli import app
from selexprep.fetch import run_fetch
from selexprep.fetch.runner import check_fetch_inventory

# ---------------------------------------------------------------------------
# Fixtures
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


def _row(
    srr: str,
    *,
    sample_title: str,
    library_name: str = "",
    experiment_title: str = "",
    paired: bool = False,
) -> dict:
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
        "library_name": library_name,
        "experiment_title": experiment_title,
        "sample_title": sample_title,
        "sample_accession": f"SAM_{srr}",
        "read_count": "10",
        "base_count": "100",
        "fastq_md5": md5,
        "fastq_bytes": sizes,
        "fastq_ftp": ftp,
    }


def _stub_download_writes_files(target_dir: Path, srr: str, paired: bool) -> bool:
    """A stand-in for download_srr that drops valid FASTQ.gz files.

    Each file is sized above ``validate_fastq_gz``'s 1024-byte minimum so
    the Codex-pass-1 resume oracle accepts it. 800 records of 60-base
    pseudo-random sequences (~80 KB raw, ~6 KB gzipped — well above
    threshold).
    """
    import gzip
    import random as _random

    rng = _random.Random(0xC0DE)
    target_dir.mkdir(parents=True, exist_ok=True)
    names = [f"{srr}_1.fastq.gz", f"{srr}_2.fastq.gz"] if paired else [f"{srr}.fastq.gz"]
    body_lines: list[str] = []
    for i in range(800):
        seq = "".join(rng.choice("ACGT") for _ in range(60))
        body_lines.append(f"@read_{i}\n{seq}\n+\n{'I' * 60}\n")
    payload = "".join(body_lines).encode()
    for name in names:
        with gzip.open(target_dir / name, "wb") as fh:
            fh.write(payload)
    return True


# ---------------------------------------------------------------------------
# run_fetch happy paths
# ---------------------------------------------------------------------------


def test_run_fetch_single_end_emits_rounds_tsv_and_metadata(tmp_path: Path) -> None:
    rows = [
        _row("SRR1", sample_title="Round 0"),
        _row("SRR2", sample_title="Round 1"),
    ]

    def fake_dl(srr: str, output_dir: Path, backend: str = "ena", **kw):
        return _stub_download_writes_files(output_dir, srr, paired=False)

    with (
        patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response(rows)),
        patch("selexprep.fetch.runner.download_srr", side_effect=fake_dl),
    ):
        result = run_fetch("PRJ", tmp_path)

    assert result.refused_reason is None
    assert sorted(result.downloaded_srrs) == ["SRR1", "SRR2"]
    assert (tmp_path / "round_00" / "SRR1.fastq.gz").exists()
    assert (tmp_path / "round_01" / "SRR2.fastq.gz").exists()

    rounds_tsv = (tmp_path / "rounds.tsv").read_text()
    lines = [ln.strip() for ln in rounds_tsv.strip().splitlines()]
    assert lines[0] == "file\tround_number"
    body = sorted(lines[1:])
    assert body == ["SRR1.fastq.gz\t0", "SRR2.fastq.gz\t1"]

    payload = json.loads((tmp_path / "fetch_metadata.json").read_text())
    assert payload["accession"] == "PRJ"
    assert {r["srr"] for r in payload["runs"]} == {"SRR1", "SRR2"}


def test_run_fetch_paired_end_emits_both_R1_R2(tmp_path: Path) -> None:
    rows = [_row("SRR3", sample_title="Round 2", paired=True)]

    def fake_dl(srr: str, output_dir: Path, backend: str = "ena", **kw):
        return _stub_download_writes_files(output_dir, srr, paired=True)

    with (
        patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response(rows)),
        patch("selexprep.fetch.runner.download_srr", side_effect=fake_dl),
    ):
        result = run_fetch("PRJ_PE", tmp_path)

    assert result.refused_reason is None
    assert (tmp_path / "round_02" / "SRR3_1.fastq.gz").exists()
    assert (tmp_path / "round_02" / "SRR3_2.fastq.gz").exists()

    rounds_lines = sorted((tmp_path / "rounds.tsv").read_text().strip().splitlines()[1:])
    assert rounds_lines == ["SRR3_1.fastq.gz\t2", "SRR3_2.fastq.gz\t2"]


# ---------------------------------------------------------------------------
# Refusal paths (user peer-review point 3)
# ---------------------------------------------------------------------------


def test_run_fetch_refuses_when_all_runs_none_confidence(tmp_path: Path) -> None:
    rows = [_row("SRR_X", sample_title="just text"), _row("SRR_Y", sample_title="also nothing")]
    with (
        patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response(rows)),
        patch("selexprep.fetch.runner.download_srr") as mock_dl,
    ):
        result = run_fetch("PRJ_ALL_NONE", tmp_path, allow_manual_review=True)

    # All-unassigned fails even with --allow-manual-review (no rounds to
    # anchor). Phase 6b.4 audit-pilot fix: the refusal message now says
    # "unassigned" rather than "NONE-confidence" because the old wording
    # was misleading once MEDIUM-single-match records stopped being
    # lumped into the refusal bucket.
    assert result.refused_reason is not None
    assert "unassigned" in result.refused_reason
    mock_dl.assert_not_called()
    assert not (tmp_path / "rounds.tsv").exists()


def test_run_fetch_refuses_none_confidence_without_manual_review_flag(tmp_path: Path) -> None:
    rows = [
        _row("SRR_OK", sample_title="Round 5"),
        _row("SRR_BAD", sample_title="nothing"),
    ]
    with (
        patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response(rows)),
        patch("selexprep.fetch.runner.download_srr") as mock_dl,
    ):
        result = run_fetch("PRJ_MIX", tmp_path)

    assert result.refused_reason is not None
    assert "SRR_BAD" in result.refused_reason
    assert "--allow-manual-review" in result.refused_reason
    mock_dl.assert_not_called()


def test_run_fetch_allow_manual_review_keeps_rounds_tsv_clean(tmp_path: Path) -> None:
    """NONE-confidence run goes to round_unknown/ + manual_review.tsv; rounds.tsv excludes it."""
    rows = [
        _row("SRR_OK", sample_title="Round 5"),
        _row("SRR_MR", sample_title="nothing"),
    ]

    def fake_dl(srr: str, output_dir: Path, backend: str = "ena", **kw):
        return _stub_download_writes_files(output_dir, srr, paired=False)

    with (
        patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response(rows)),
        patch("selexprep.fetch.runner.download_srr", side_effect=fake_dl),
    ):
        result = run_fetch("PRJ_MIX", tmp_path, allow_manual_review=True)

    assert result.refused_reason is None
    assert (tmp_path / "round_05" / "SRR_OK.fastq.gz").exists()
    assert (tmp_path / "round_unknown" / "SRR_MR.fastq.gz").exists()

    rounds_body = sorted((tmp_path / "rounds.tsv").read_text().strip().splitlines()[1:])
    assert rounds_body == ["SRR_OK.fastq.gz\t5"]  # SRR_MR excluded

    mr_body = (tmp_path / "manual_review.tsv").read_text().strip().splitlines()
    assert mr_body[0].startswith("file\tsrr\t")
    assert any("SRR_MR" in line for line in mr_body[1:])


# ---------------------------------------------------------------------------
# Resume oracle (user peer-review point 5)
# ---------------------------------------------------------------------------


def test_check_fetch_inventory_detects_missing_fastq(tmp_path: Path) -> None:
    rows = [_row("SRR_RESUME", sample_title="Round 0", paired=True)]

    def fake_dl(srr: str, output_dir: Path, backend: str = "ena", **kw):
        return _stub_download_writes_files(output_dir, srr, paired=True)

    with (
        patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response(rows)),
        patch("selexprep.fetch.runner.download_srr", side_effect=fake_dl),
    ):
        result = run_fetch("PRJ_RESUME", tmp_path)

    # All files present → empty missing list.
    assert check_fetch_inventory(result.plan, tmp_path) == []

    # Delete one R2 → resume oracle should flag it.
    (tmp_path / "round_00" / "SRR_RESUME_2.fastq.gz").unlink()
    missing = check_fetch_inventory(result.plan, tmp_path)
    assert missing == ["SRR_RESUME_2.fastq.gz"]


def test_run_fetch_skips_already_present_fastqs(tmp_path: Path) -> None:
    rows = [_row("SRR_SKIP", sample_title="Round 0")]
    # Pre-populate the file.
    target_dir = tmp_path / "round_00"
    _stub_download_writes_files(target_dir, "SRR_SKIP", paired=False)

    def fake_dl(srr: str, output_dir: Path, backend: str = "ena", **kw):
        raise AssertionError("download_srr should not be called when file already present")

    with (
        patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response(rows)),
        patch("selexprep.fetch.runner.download_srr", side_effect=fake_dl),
    ):
        result = run_fetch("PRJ_SKIP", tmp_path)

    assert result.downloaded_srrs == []
    assert result.skipped_srrs == ["SRR_SKIP"]


# ---------------------------------------------------------------------------
# CLI smoke tests
# ---------------------------------------------------------------------------


def test_cli_fetch_dry_run_no_download(tmp_path: Path) -> None:
    rows = [_row("SRR1", sample_title="Round 0")]
    with (
        patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response(rows)),
        patch("selexprep.fetch.runner.download_srr") as mock_dl,
    ):
        result = CliRunner().invoke(
            app, ["fetch", "PRJ_DRY", "--outdir", str(tmp_path), "--dry-run"]
        )

    assert result.exit_code == 0, result.output
    mock_dl.assert_not_called()
    assert (tmp_path / "fetch_metadata.json").exists()
    assert (tmp_path / "rounds.tsv").exists()


def test_cli_fetch_invalid_backend_exits_non_zero(tmp_path: Path) -> None:
    result = CliRunner().invoke(
        app, ["fetch", "PRJ", "--outdir", str(tmp_path), "--backend", "bogus"]
    )
    assert result.exit_code == 2
    assert "invalid --backend" in result.output


def test_cli_fetch_refusal_exits_with_code_2(tmp_path: Path) -> None:
    rows = [_row("SRR_NONE", sample_title="nothing useful")]
    with patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response(rows)):
        result = CliRunner().invoke(app, ["fetch", "PRJ", "--outdir", str(tmp_path)])
    assert result.exit_code == 2
    # Phase 6b.4 audit-pilot fix: refusal message says "unassigned" now,
    # not "NONE-confidence". See tests/test_metadata.py for the
    # MEDIUM-single-match regression suite that motivated the wording
    # change.
    assert "unassigned" in result.output


def test_cli_fetch_propagates_value_error(tmp_path: Path) -> None:
    with patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response([])):
        result = CliRunner().invoke(app, ["fetch", "SRR_MISSING", "--outdir", str(tmp_path)])
    assert result.exit_code == 2
    assert "ENA returned no records" in result.output


@pytest.mark.parametrize("backend", ["auto", "ena", "kingfisher", "sra"])
def test_cli_fetch_accepts_all_documented_backends(backend: str, tmp_path: Path) -> None:
    """Document that the four backend names are accepted by the CLI parser."""
    rows = [_row("SRR1", sample_title="Round 0")]

    def fake_dl(srr: str, output_dir: Path, backend: str = "ena", **kw):
        return _stub_download_writes_files(output_dir, srr, paired=False)

    with (
        patch("selexprep.fetch.inspect.requests.get", return_value=_mock_response(rows)),
        patch("selexprep.fetch.runner.download_srr", side_effect=fake_dl),
    ):
        result = CliRunner().invoke(
            app,
            ["fetch", "PRJ", "--outdir", str(tmp_path), "--backend", backend],
        )
    assert result.exit_code == 0, result.output
