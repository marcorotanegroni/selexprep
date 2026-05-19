"""Deterministic I/O helpers for reproducible-output guarantees.

The selexprep manifest promises bit-identical SHA256 across reruns for
FASTA/FASTQ/TSV/JSON outputs. The default ``gzip.open`` writer embeds the
current ``mtime`` in the gzip header, which breaks this guarantee silently.
``open_gzip_text_deterministic`` forces ``mtime=0`` so the output bytes are
reproducible across machines and across reruns on the same machine.

See the plan's "Reproducible-output discipline" section: every ``.gz`` write
must go through this module.
"""

from __future__ import annotations

import gzip
import hashlib
import io
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import IO


def open_gzip_text_deterministic(
    path: Path,
    encoding: str = "utf-8",
) -> io.TextIOWrapper:
    """Open a deterministic-gzip writer in text mode.

    The gzip format normally embeds two non-deterministic bits in the header:

    1. ``mtime`` — the wall-clock time of writing.
    2. ``FNAME`` — the destination filename (so two runs writing to ``a.gz``
       and ``b.gz`` produce different headers even for identical content).

    We suppress both: ``mtime=0`` and an empty ``filename=""`` (with an
    explicit ``fileobj``) prevent either field from making it into the
    header. The resulting bytes depend only on the payload, so any two
    invocations with the same content yield byte-for-byte identical output
    — which is the property the manifest's SHA256 checks rely on.

    The returned ``TextIOWrapper`` is the caller's to close. Closing it
    cascades through the inner ``GzipFile`` to the underlying raw file
    because we set ``gz.myfileobj = raw``, which is how stdlib ``gzip``
    flags file ownership for its own close path.

    Example::

        with open_gzip_text_deterministic(Path("out.fastq.gz")) as fh:
            fh.write("@read_0\\nACGT\\n+\\nIIII\\n")
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = open(path, "wb")  # noqa: SIM115 — closed via GzipFile.myfileobj cascade
    gz = gzip.GzipFile(filename="", mode="wb", mtime=0, fileobj=raw)
    # Force GzipFile.close() to close `raw` for us so the caller-visible close
    # on the TextIOWrapper produces a clean cascade with no file-handle leak.
    gz.myfileobj = raw  # type: ignore[attr-defined]
    return io.TextIOWrapper(gz, encoding=encoding, write_through=True)


@contextmanager
def write_gzip_text_deterministic(
    path: Path,
    encoding: str = "utf-8",
) -> Iterator[IO[str]]:
    """Context-managed version of :func:`open_gzip_text_deterministic`."""
    fh = open_gzip_text_deterministic(path, encoding=encoding)
    try:
        yield fh
    finally:
        fh.close()


def sha256_file(path: Path, chunk: int = 1 << 20) -> str:
    """Stream-SHA256 a file in 1 MiB chunks (memory-bounded for multi-GB inputs)."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()
