"""Unit tests for selexprep.fetch.download (offline only — network calls are
not exercised here; mark @pytest.mark.network in follow-up integration tests)."""

from __future__ import annotations

import hashlib
from pathlib import Path
from unittest import mock

from selexprep.fetch.download import (
    download_srr,
    kingfisher_available,
    md5_file,
    needs_manual_review,
    round_dir,
    safe_dir_name,
    sratoolkit_available,
    srr_present,
    validate_fastq_gz,
    zenodo_expected_md5,
)

# ----- safe_dir_name -----


def test_safe_dir_name_replaces_colons() -> None:
    assert safe_dir_name("zenodo:1234567") == "zenodo_1234567"


def test_safe_dir_name_passes_clean_bp() -> None:
    assert safe_dir_name("PRJNA315881") == "PRJNA315881"


# ----- round_dir -----


def test_round_dir_single_target(tmp_path: Path) -> None:
    out = round_dir(tmp_path, "PRJ1", 3)
    assert out == tmp_path / "PRJ1" / "round_03"


def test_round_dir_multi_target(tmp_path: Path) -> None:
    out = round_dir(tmp_path, "PRJ1", 3, target_hint="VEGF")
    assert out == tmp_path / "PRJ1" / "VEGF" / "round_03"


def test_round_dir_unknown_round(tmp_path: Path) -> None:
    out = round_dir(tmp_path, "PRJ1", None)
    assert out == tmp_path / "PRJ1" / "round_unknown"


def test_round_dir_zero_padded() -> None:
    out = round_dir(Path("/r"), "PRJ1", 5)
    assert out.name == "round_05"


def test_round_dir_handles_colons_in_bpid(tmp_path: Path) -> None:
    out = round_dir(tmp_path, "zenodo:123", 1)
    assert "zenodo_123" in str(out)


# ----- needs_manual_review -----


def test_needs_manual_review_truthy() -> None:
    for raw in ("true", "TRUE", "1", "yes", "y"):
        assert needs_manual_review({"needs_manual_review": raw})


def test_needs_manual_review_falsy() -> None:
    for raw in ("false", "0", "no", "", "None"):
        assert not needs_manual_review({"needs_manual_review": raw})


def test_needs_manual_review_missing_key() -> None:
    assert not needs_manual_review({})


# ----- srr_present -----


def test_srr_present_finds_single_end(tmp_path: Path) -> None:
    bp_dir = tmp_path / "PRJ1"
    bp_dir.mkdir()
    (bp_dir / "SRR1.fastq.gz").touch()
    assert srr_present(tmp_path, "PRJ1", "SRR1")


def test_srr_present_false_when_missing(tmp_path: Path) -> None:
    assert not srr_present(tmp_path, "PRJ1", "SRR_missing")


def test_srr_present_no_prefix_collision(tmp_path: Path) -> None:
    bp_dir = tmp_path / "PRJ1"
    bp_dir.mkdir()
    (bp_dir / "SRR12345.fastq.gz").touch()
    # Should NOT match SRR1234 even though SRR12345 exists
    assert not srr_present(tmp_path, "PRJ1", "SRR1234")


# ----- md5_file -----


def test_md5_file_matches_hashlib(tmp_path: Path) -> None:
    p = tmp_path / "data.bin"
    payload = b"ACGT" * 10_000
    p.write_bytes(payload)
    expected = hashlib.md5(payload).hexdigest()
    assert md5_file(p) == expected


def test_md5_file_streams_chunks(tmp_path: Path) -> None:
    """File larger than chunk size still produces the canonical MD5."""
    p = tmp_path / "big.bin"
    chunk = 1 << 20  # 1 MiB
    payload = b"X" * (chunk * 3 + 17)
    p.write_bytes(payload)
    assert md5_file(p) == hashlib.md5(payload).hexdigest()


# ----- validate_fastq_gz -----


def test_validate_fastq_gz_rejects_missing(tmp_path: Path) -> None:
    assert not validate_fastq_gz(tmp_path / "nope.fastq.gz")


def test_validate_fastq_gz_rejects_too_small(tmp_path: Path) -> None:
    p = tmp_path / "stub.fastq.gz"
    p.write_bytes(b"X" * 10)
    assert not validate_fastq_gz(p, min_size=1024)


def test_validate_fastq_gz_accepts_valid_gzip(tmp_path: Path) -> None:
    """Write a real gzipped file (any content) and confirm gzip -t passes."""
    import gzip

    p = tmp_path / "ok.fastq.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("@read_0\nACGT\n+\nIIII\n" * 100)
    assert validate_fastq_gz(p, min_size=10)


def test_validate_fastq_gz_rejects_corrupt(tmp_path: Path) -> None:
    p = tmp_path / "broken.fastq.gz"
    # Garbage bytes long enough to pass the size check but not gzip -t
    p.write_bytes(b"\xff" * 2048)
    assert not validate_fastq_gz(p)


# ----- zenodo_expected_md5 -----


def test_zenodo_md5_strips_prefix() -> None:
    assert (
        zenodo_expected_md5("md5:abc123def456abc123def456abc123ab")
        == "abc123def456abc123def456abc123ab"
    )


def test_zenodo_md5_passes_bare_hex() -> None:
    assert (
        zenodo_expected_md5("abc123def456abc123def456abc123ab")
        == "abc123def456abc123def456abc123ab"
    )


def test_zenodo_md5_rejects_non_md5() -> None:
    assert zenodo_expected_md5("sha256:longhash") is None
    assert zenodo_expected_md5(None) is None
    assert zenodo_expected_md5("") is None


# ----- Backend probes -----


def test_kingfisher_available_true_when_found() -> None:
    with mock.patch("selexprep.fetch.download.shutil.which", return_value="/path/kingfisher"):
        assert kingfisher_available()


def test_kingfisher_available_false_when_missing() -> None:
    with mock.patch("selexprep.fetch.download.shutil.which", return_value=None):
        assert not kingfisher_available()


def test_sratoolkit_available_requires_both_binaries() -> None:
    def fake_which(name: str) -> str | None:
        return "/path/" + name if name in ("prefetch", "fasterq-dump") else None

    with mock.patch("selexprep.fetch.download.shutil.which", side_effect=fake_which):
        assert sratoolkit_available()

    with mock.patch(
        "selexprep.fetch.download.shutil.which",
        side_effect=lambda n: "/path/prefetch" if n == "prefetch" else None,
    ):
        assert not sratoolkit_available()


# ----- Dispatcher order (ENA-direct first; licensing-correct default) -----


def test_download_srr_tries_ena_direct_first(tmp_path: Path) -> None:
    """Dispatcher MUST try ENA-direct first so the MIT-only default install works
    without kingfisher (GPL-3.0). Regression guard for the licensing-driven
    backend order."""
    with (
        mock.patch("selexprep.fetch.download.download_srr_ena_direct", return_value=True) as ena,
        mock.patch("selexprep.fetch.download.download_srr_kingfisher", return_value=True) as kf,
        mock.patch("selexprep.fetch.download.download_srr_sratoolkit", return_value=True) as sra,
    ):
        result = download_srr("SRR1", tmp_path)
        assert result is True
        ena.assert_called_once()
        kf.assert_not_called()
        sra.assert_not_called()


def test_download_srr_falls_back_to_kingfisher_if_ena_fails(tmp_path: Path) -> None:
    """When ENA-direct fails and kingfisher is installed, kingfisher is tried."""
    with (
        mock.patch("selexprep.fetch.download.download_srr_ena_direct", return_value=False),
        mock.patch("selexprep.fetch.download.kingfisher_available", return_value=True),
        mock.patch("selexprep.fetch.download.download_srr_kingfisher", return_value=True) as kf,
        mock.patch("selexprep.fetch.download.sratoolkit_available", return_value=False),
    ):
        assert download_srr("SRR1", tmp_path)
        kf.assert_called_once()


def test_download_srr_skips_kingfisher_when_not_available(tmp_path: Path) -> None:
    """Default install (no kingfisher binary) must still produce a clean fall-through."""
    with (
        mock.patch("selexprep.fetch.download.download_srr_ena_direct", return_value=False),
        mock.patch("selexprep.fetch.download.kingfisher_available", return_value=False),
        mock.patch("selexprep.fetch.download.sratoolkit_available", return_value=False),
    ):
        assert not download_srr("SRR1", tmp_path)


def test_download_srr_returns_false_when_all_backends_exhausted(tmp_path: Path) -> None:
    with (
        mock.patch("selexprep.fetch.download.download_srr_ena_direct", return_value=False),
        mock.patch("selexprep.fetch.download.kingfisher_available", return_value=True),
        mock.patch("selexprep.fetch.download.download_srr_kingfisher", return_value=False),
        mock.patch("selexprep.fetch.download.sratoolkit_available", return_value=True),
        mock.patch("selexprep.fetch.download.download_srr_sratoolkit", return_value=False),
    ):
        assert not download_srr("SRR1", tmp_path)
