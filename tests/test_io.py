"""Unit tests for selexprep._io (deterministic output helpers)."""

from __future__ import annotations

import gzip
from pathlib import Path

from selexprep._io import (
    open_gzip_text_deterministic,
    sha256_file,
    write_gzip_text_deterministic,
)

# ----- deterministic gzip -----


def test_deterministic_gzip_produces_bit_identical_bytes_across_runs(tmp_path: Path) -> None:
    """Two writes of the same content to DIFFERENT destination paths must
    produce byte-for-byte identical files. Both `mtime` and `FNAME` must be
    suppressed — the gzip format normally embeds both."""
    a = tmp_path / "a.fastq.gz"
    b = tmp_path / "subdir" / "b.fastq.gz"

    payload = "@read_0\nACGT\n+\nIIII\n" * 100

    with write_gzip_text_deterministic(a) as fh:
        fh.write(payload)
    with write_gzip_text_deterministic(b) as fh:
        fh.write(payload)

    # Critical: the BYTES on disk must match — not just the decompressed content.
    # This is what makes manifest SHA256 reproducible across reruns AND across
    # rename / move operations.
    assert a.read_bytes() == b.read_bytes()
    assert sha256_file(a) == sha256_file(b)


def test_deterministic_gzip_header_suppresses_mtime_and_fname(tmp_path: Path) -> None:
    """Inspect the raw gzip header bytes to confirm both mtime (bytes 4-7) and
    the FNAME flag bit (bit 3 of byte 3) are cleared."""
    p = tmp_path / "x.fastq.gz"
    with write_gzip_text_deterministic(p) as fh:
        fh.write("ACGT\n")
    raw = p.read_bytes()
    # gzip header: magic (2) + CM (1) + FLG (1) + MTIME (4) + XFL (1) + OS (1)
    assert raw[:2] == b"\x1f\x8b", "missing gzip magic"
    flg = raw[3]
    assert flg & 0b00001000 == 0, "FNAME flag must be cleared"  # FLG.FNAME
    mtime = int.from_bytes(raw[4:8], "little")
    assert mtime == 0, f"mtime must be 0, got {mtime}"


def test_open_gzip_text_deterministic_caller_must_close(tmp_path: Path) -> None:
    """Non-context-manager variant returns a writer the caller must close."""
    path = tmp_path / "x.fastq.gz"
    fh = open_gzip_text_deterministic(path)
    fh.write("hello\n")
    fh.close()

    # Round-trip read
    with gzip.open(path, "rt") as g:
        assert g.read() == "hello\n"


def test_deterministic_gzip_creates_parent_dirs(tmp_path: Path) -> None:
    nested = tmp_path / "a" / "b" / "c" / "out.fastq.gz"
    with write_gzip_text_deterministic(nested) as fh:
        fh.write("x\n")
    assert nested.exists()


# ----- sha256_file -----


def test_sha256_file_matches_hashlib(tmp_path: Path) -> None:
    import hashlib

    p = tmp_path / "data.bin"
    payload = b"ACGT" * 10_000
    p.write_bytes(payload)
    assert sha256_file(p) == hashlib.sha256(payload).hexdigest()


def test_sha256_file_streams_large_input(tmp_path: Path) -> None:
    """File larger than chunk size still produces canonical SHA256."""
    import hashlib

    p = tmp_path / "big.bin"
    chunk = 1 << 20
    payload = b"X" * (chunk * 3 + 17)
    p.write_bytes(payload)
    assert sha256_file(p) == hashlib.sha256(payload).hexdigest()
