"""Shared helpers for selexprep modules."""

from __future__ import annotations

import csv
import logging
import sys
from datetime import datetime
from pathlib import Path


# Canonical per-SRR FASTQ filename conventions produced by kingfisher /
# sra-toolkit. Single-end runs produce `{srr}.fastq.gz`; paired-end runs
# produce `{srr}_1.fastq.gz` and `{srr}_2.fastq.gz`. Match EXACTLY to avoid
# prefix collisions (SRR1234 must not match SRR12345).


def iter_srr_files(root: Path, srr: str, suffix: str = ".fastq.gz") -> list[Path]:
    """Return FASTQ files for `srr` under `root`, matching exact filenames only.

    Supports both single-end (`{srr}{suffix}`) and paired-end
    (`{srr}_1{suffix}`, `{srr}_2{suffix}`) layouts.
    """
    if not root.exists():
        return []
    valid_names = {f"{srr}{suffix}", f"{srr}_1{suffix}", f"{srr}_2{suffix}"}
    return [p for p in root.rglob(f"*{suffix}") if p.name in valid_names]


def load_csv(path: Path) -> list[dict[str, str]]:
    """Load a CSV file into a list of dicts. Returns `[]` if the file is missing."""
    try:
        with open(path, encoding="utf-8") as f:
            return list(csv.DictReader(f))
    except FileNotFoundError:
        logging.getLogger(__name__).warning("File not found: %s", path)
        return []


def parse_round_number(parquet_path: Path) -> int | str:
    """Extract the round number from a per-round parquet filename.

    Handles `.counts.parquet` and `.clusters.parquet` suffixes. Returns the
    integer round when the trailing token parses as int, otherwise returns
    the bare stem (so callers can detect and route non-round files).
    """
    stem = parquet_path.stem.replace(".counts", "").replace(".clusters", "")
    tail = stem.split("_")[-1]
    try:
        return int(tail)
    except ValueError:
        return stem


def setup_logging(
    logger_name: str,
    log_dir: Path | None = None,
    log_prefix: str | None = None,
) -> logging.Logger:
    """Configure root logging with optional rotating per-run file handler.

    Idempotent: only attaches a stream handler once. When `log_dir` is given,
    a timestamped log file is added on top of the stream handler.
    """
    root = logging.getLogger()
    if not root.handlers:
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
            handlers=[logging.StreamHandler(sys.stdout)],
        )
    if log_dir is not None:
        log_dir.mkdir(exist_ok=True)
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        prefix = log_prefix or logger_name
        fh = logging.FileHandler(log_dir / f"{prefix}_{ts}.log", encoding="utf-8")
        fh.setFormatter(
            logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s")
        )
        root.addHandler(fh)
    return logging.getLogger(logger_name)
