"""Cutadapt subprocess wrapper for primer-aware extraction.

the design: "cutadapt invoked as subprocess (CLI is the stable
contract, not the Python API)." Per-extraction-mode adapter flags are
applied accordingly.

**Determinism discipline.** Cutadapt's own gzip writer embeds the current
``mtime`` in the gzip header, breaking ``output_sha256`` reproducibility.
This module therefore:

1. Tells cutadapt to write **uncompressed** ``.fa`` output.
2. Re-gzips that ``.fa`` with ``_io.open_gzip_text_deterministic``
   (``mtime=0`` header) to ``<name>.fasta.gz``.
3. Deletes the intermediate ``.fa``.

Cutadapt's structured ``--json`` report supplies ``n_in`` / ``n_out``
counts without fragile stderr parsing. The exact argv is recorded in
:class:`TrimReport` so the manifest can reproduce the invocation
verbatim.

Public API:

- :func:`trim_single_end_linked` — BOTH_PRIMERS_SINGLE_READ
- :func:`trim_single_end_5p` — FIVE_PRIME_ONLY
- :func:`trim_single_end_3p` — THREE_PRIME_ONLY
- :func:`trim_paired_split` — PAIRED_END_SPLIT_PRIMERS
- :class:`TrimReport` — exact argv + read counts for each invocation
"""

from __future__ import annotations

import json
import logging
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

from selexprep._common import resolve_cutadapt
from selexprep._io import open_gzip_text_deterministic

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class TrimReport:
    """Exact cutadapt invocation + I/O counts (for manifest reproducibility)."""

    cutadapt_cmd: list[str]
    n_in: int
    n_out: int
    return_code: int
    output_paths: list[Path] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _run_cutadapt(argv: list[str], json_report_path: Path) -> tuple[int, dict]:
    """Run cutadapt; return (return_code, parsed_json_report).

    Captures stderr only on failure to keep the happy path quiet. Mirrors
    the subprocess pattern in `selexprep.count.counter._run_cutadapt`.

    ``argv[0]`` is the literal ``"cutadapt"`` (recorded verbatim in the
    manifest for portability); execution substitutes the resolved absolute
    path so the call works even when the environment isn't activated.
    """
    exe = resolve_cutadapt()
    if exe is None:
        raise RuntimeError(
            "cutadapt not found on PATH or alongside the Python interpreter. "
            "Install with `pip install cutadapt` or `uv add cutadapt`; "
            "selexprep declares it as a core dependency."
        )
    full = [exe, *argv[1:], "--json", str(json_report_path)]
    logger.info("cutadapt argv: %s", " ".join(full))
    try:
        subprocess.run(full, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        tail = (e.stderr or "")[-500:]
        logger.error("cutadapt failed (exit %d):\n%s", e.returncode, tail)
        raise RuntimeError(f"cutadapt failed: {tail}") from e

    report: dict = json.loads(json_report_path.read_text(encoding="utf-8"))
    return 0, report


def _read_counts(report: dict) -> tuple[int, int]:
    """Pull (n_in, n_out) from a cutadapt --json report.

    Cutadapt's JSON schema (5.x) puts these under ``read_counts.input`` and
    ``read_counts.output`` for single-end; for paired-end the same keys
    refer to read-pair counts.
    """
    rc = report.get("read_counts", {})
    return int(rc.get("input", 0)), int(rc.get("output", 0))


def _repack_fasta_deterministic(tmp_fa: Path, output_fa_gz: Path) -> None:
    """Re-gzip a cutadapt-produced uncompressed FASTA with mtime=0 header.

    Writes to ``output_fa_gz`` byte-for-byte deterministically (so the
    manifest's ``output_sha256`` is stable across reruns) and removes the
    intermediate ``tmp_fa``.
    """
    output_fa_gz.parent.mkdir(parents=True, exist_ok=True)
    with (
        open(tmp_fa, encoding="utf-8") as fh_in,
        open_gzip_text_deterministic(output_fa_gz) as fh_out,
    ):
        for line in fh_in:
            fh_out.write(line)
    tmp_fa.unlink()


# ---------------------------------------------------------------------------
# Public trim functions — one per extraction_mode
# ---------------------------------------------------------------------------


def trim_single_end_linked(
    input_fastq: Path,
    output_fasta_gz: Path,
    *,
    primer_5p: str,
    primer_3p: str,
) -> TrimReport:
    """Trim 5' + 3' primers in linked-adapter mode (BOTH_PRIMERS_SINGLE_READ).

    The ``-g "P5...P3"`` syntax tells cutadapt to require BOTH adapters
    in the same read, trimming the random region between them.
    ``--discard-untrimmed`` drops reads that don't carry the full
    linked-adapter pattern.
    """
    output_fasta_gz.parent.mkdir(parents=True, exist_ok=True)
    tmp_fa = output_fasta_gz.with_suffix("").with_suffix(".fa")
    tmp_json = output_fasta_gz.with_suffix(".cutadapt.json")

    argv = [
        "cutadapt",
        "--discard-untrimmed",
        "--fasta",
        "-g",
        f"{primer_5p}...{primer_3p}",
        "-o",
        str(tmp_fa),
        str(input_fastq),
    ]
    try:
        rc, report = _run_cutadapt(argv, tmp_json)
        n_in, n_out = _read_counts(report)
        _repack_fasta_deterministic(tmp_fa, output_fasta_gz)
        return TrimReport(
            cutadapt_cmd=argv,
            n_in=n_in,
            n_out=n_out,
            return_code=rc,
            output_paths=[output_fasta_gz],
        )
    finally:
        # Belt-and-suspenders: _repack_fasta_deterministic unlinks tmp_fa on
        # success; this catches the failure paths (cutadapt raise, repack
        # raise) so we don't orphan intermediate files.
        tmp_fa.unlink(missing_ok=True)
        tmp_json.unlink(missing_ok=True)


def trim_single_end_5p(
    input_fastq: Path,
    output_fasta_gz: Path,
    *,
    primer_5p: str,
) -> TrimReport:
    """Trim 5' primer only (FIVE_PRIME_ONLY).

    Emits ``partial_5p_extracted.fasta.gz`` — the filename signals to
    downstream ML pipelines that the 3' end is unbounded (random region
    length unknown per read).
    """
    output_fasta_gz.parent.mkdir(parents=True, exist_ok=True)
    tmp_fa = output_fasta_gz.with_suffix("").with_suffix(".fa")
    tmp_json = output_fasta_gz.with_suffix(".cutadapt.json")

    argv = [
        "cutadapt",
        "--discard-untrimmed",
        "--fasta",
        "-g",
        primer_5p,
        "-o",
        str(tmp_fa),
        str(input_fastq),
    ]
    try:
        rc, report = _run_cutadapt(argv, tmp_json)
        n_in, n_out = _read_counts(report)
        _repack_fasta_deterministic(tmp_fa, output_fasta_gz)
        return TrimReport(
            cutadapt_cmd=argv,
            n_in=n_in,
            n_out=n_out,
            return_code=rc,
            output_paths=[output_fasta_gz],
        )
    finally:
        tmp_fa.unlink(missing_ok=True)
        tmp_json.unlink(missing_ok=True)


def trim_single_end_3p(
    input_fastq: Path,
    output_fasta_gz: Path,
    *,
    primer_3p: str,
) -> TrimReport:
    """Trim 3' primer only (THREE_PRIME_ONLY).

    Emits ``partial_3p_extracted.fasta.gz`` — 5' end is unbounded.
    """
    output_fasta_gz.parent.mkdir(parents=True, exist_ok=True)
    tmp_fa = output_fasta_gz.with_suffix("").with_suffix(".fa")
    tmp_json = output_fasta_gz.with_suffix(".cutadapt.json")

    argv = [
        "cutadapt",
        "--discard-untrimmed",
        "--fasta",
        "-a",
        primer_3p,
        "-o",
        str(tmp_fa),
        str(input_fastq),
    ]
    try:
        rc, report = _run_cutadapt(argv, tmp_json)
        n_in, n_out = _read_counts(report)
        _repack_fasta_deterministic(tmp_fa, output_fasta_gz)
        return TrimReport(
            cutadapt_cmd=argv,
            n_in=n_in,
            n_out=n_out,
            return_code=rc,
            output_paths=[output_fasta_gz],
        )
    finally:
        tmp_fa.unlink(missing_ok=True)
        tmp_json.unlink(missing_ok=True)


def trim_paired_split(
    r1_fastq: Path,
    r2_fastq: Path,
    out_r1_fasta_gz: Path,
    out_r2_fasta_gz: Path,
    *,
    primer_5p_r1: str,
    primer_5p_r2: str,
) -> TrimReport:
    """Trim a paired-end split-primer library (PAIRED_END_SPLIT_PRIMERS).

    R1 carries the forward 5' primer; R2 carries
    ``reverse_complement(real_3p)`` at its 5' end. cutadapt's paired mode
    enforces pair-sync (both reads kept or both dropped). Two separate
    FASTA outputs are emitted — joining by read ID alone without merging
    is biologically wrong, so v0.1 keeps them split.
    """
    out_r1_fasta_gz.parent.mkdir(parents=True, exist_ok=True)
    out_r2_fasta_gz.parent.mkdir(parents=True, exist_ok=True)
    tmp_r1 = out_r1_fasta_gz.with_suffix("").with_suffix(".fa")
    tmp_r2 = out_r2_fasta_gz.with_suffix("").with_suffix(".fa")
    tmp_json = out_r1_fasta_gz.with_suffix(".cutadapt.json")

    argv = [
        "cutadapt",
        "--discard-untrimmed",
        "--fasta",
        "-g",
        primer_5p_r1,
        "-G",
        primer_5p_r2,
        "-o",
        str(tmp_r1),
        "-p",
        str(tmp_r2),
        str(r1_fastq),
        str(r2_fastq),
    ]
    try:
        rc, report = _run_cutadapt(argv, tmp_json)
        n_in, n_out = _read_counts(report)
        _repack_fasta_deterministic(tmp_r1, out_r1_fasta_gz)
        _repack_fasta_deterministic(tmp_r2, out_r2_fasta_gz)
        return TrimReport(
            cutadapt_cmd=argv,
            n_in=n_in,
            n_out=n_out,
            return_code=rc,
            output_paths=[out_r1_fasta_gz, out_r2_fasta_gz],
        )
    finally:
        tmp_r1.unlink(missing_ok=True)
        tmp_r2.unlink(missing_ok=True)
        tmp_json.unlink(missing_ok=True)
