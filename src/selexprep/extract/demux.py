"""Barcode-driven demultiplexer for SELEX runs that pool multiple rounds.

Some SELEX submissions pack multiple rounds into a single sequencing run, with
each round tagged by a distinct 5' barcode. The de-pooling must be done
before per-round counting; otherwise the rounds collide in the count tables.

**v0.1 scope (Codex-locked):** sample-sheet / sample-supplied barcodes only.
Barcode inference (auto-discovery without prior knowledge) is v0.2 and not
implemented here.

Public API:
- ``demux_fastq()`` — split one multiplexed FASTQ (single- or paired-end)
  into per-round files in ``<out_dir>/round_NN/`` and ``<out_dir>/unassigned/``.
- ``read_sample_sheet()`` — parse a minimal TSV sample sheet into
  ``DemuxJob`` records suitable for batch invocation.
- ``validate_barcodes()`` — pre-flight check that barcodes are Hamming-distant
  enough that ``max_mismatches`` cannot alias them.
- ``demux_sample_sheet()`` — batch helper: run every job in a sample sheet.

For paired-end input, the round is decided from R1's 5' barcode and BOTH R1
and R2 are written to the destination round folder, kept in pair-sync so that
downstream paired-mode trimming works.
"""

from __future__ import annotations

import csv
import gzip
import json
import logging
from collections import Counter
from collections.abc import Iterator
from dataclasses import asdict, dataclass, field
from pathlib import Path

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclass
class DemuxJob:
    """One row of a sample sheet — what to demux."""

    srr: str
    r1_path: Path
    r2_path: Path | None
    barcodes: dict[str, int]  # barcode sequence → round number


@dataclass
class DemuxReport:
    """Outcome of one demultiplexing run."""

    srr: str
    input_r1: str
    input_r2: str | None
    paired: bool
    total_reads: int
    assigned_reads: int
    unassigned_reads: int
    per_round: dict[int, int] = field(default_factory=dict)
    unassigned_fraction: float = 0.0
    max_mismatches: int = 1
    trim_barcode: bool = True


# ---------------------------------------------------------------------------
# Barcode matching
# ---------------------------------------------------------------------------


def _hamming(a: str, b: str) -> int:
    return sum(1 for x, y in zip(a, b, strict=False) if x != y)


def _match_barcode(seq: str, barcodes: dict[str, int], max_mm: int) -> int | None:
    """Return the round number if any barcode matches the read's 5' prefix.

    First-match-wins — caller must ensure barcodes are mutually distinguishable
    via ``validate_barcodes()`` first.
    """
    for bc, rn in barcodes.items():
        if len(seq) < len(bc):
            continue
        if _hamming(seq[: len(bc)], bc) <= max_mm:
            return rn
    return None


def validate_barcodes(barcodes: dict[str, int], max_mismatches: int) -> None:
    """Raise ValueError if any pair of barcodes is too close to be distinguished.

    Barcodes must be pairwise Hamming-distance ≥ ``2 * max_mismatches + 1``
    (the standard error-correcting-code constraint) and have identical length.
    """
    items = list(barcodes.items())
    min_allowed = 2 * max_mismatches + 1
    for i, (bc_i, rn_i) in enumerate(items):
        for bc_j, rn_j in items[i + 1 :]:
            if len(bc_i) != len(bc_j):
                raise ValueError(
                    f"barcodes have different lengths: {bc_i} ({len(bc_i)}) vs "
                    f"{bc_j} ({len(bc_j)})"
                )
            d = _hamming(bc_i, bc_j)
            if d < min_allowed:
                raise ValueError(
                    f"barcodes {bc_i} (round {rn_i}) and {bc_j} (round {rn_j}) "
                    f"are Hamming-distance {d}; need ≥{min_allowed} for "
                    f"max_mismatches={max_mismatches}"
                )


# ---------------------------------------------------------------------------
# FASTQ iteration
# ---------------------------------------------------------------------------


def _iter_fastq(fq_path: Path) -> Iterator[tuple[str, str, str, str]]:
    """Yield (header, seq, plus, qual) tuples (newline-terminated) from a FASTQ.gz."""
    with gzip.open(fq_path, "rt", encoding="utf-8", errors="replace") as fh:
        while True:
            header = fh.readline()
            if not header:
                return
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            if not qual:
                return
            yield header, seq, plus, qual


def _iter_fastq_pair(
    r1_path: Path, r2_path: Path
) -> Iterator[tuple[tuple[str, str, str, str], tuple[str, str, str, str]]]:
    """Yield ((h1,s1,p1,q1), (h2,s2,p2,q2)) in lockstep.

    Raises if R1 and R2 diverge in record count (would corrupt pair sync).
    """
    with gzip.open(r1_path, "rt", encoding="utf-8", errors="replace") as f1, gzip.open(
        r2_path, "rt", encoding="utf-8", errors="replace"
    ) as f2:
        while True:
            r1 = (f1.readline(), f1.readline(), f1.readline(), f1.readline())
            r2 = (f2.readline(), f2.readline(), f2.readline(), f2.readline())
            if not r1[0] and not r2[0]:
                return
            if not r1[0] or not r2[0]:
                raise RuntimeError(
                    f"R1/R2 record count mismatch between {r1_path.name} and {r2_path.name}"
                )
            if not r1[3] or not r2[3]:
                return
            yield r1, r2


# ---------------------------------------------------------------------------
# Single-FASTQ demux
# ---------------------------------------------------------------------------


def demux_fastq(
    r1_path: Path,
    out_dir: Path,
    barcodes: dict[str, int],
    srr: str,
    r2_path: Path | None = None,
    max_mismatches: int = 1,
    trim_barcode: bool = True,
) -> DemuxReport:
    """Split a multiplexed FASTQ (SE or PE) into per-round files by 5' barcode.

    Output layout under `out_dir`:
        out_dir/round_NN/<srr>.fastq.gz                  (SE)
        out_dir/round_NN/<srr>_1.fastq.gz + _2.fastq.gz  (PE)
        out_dir/unassigned/<srr>(_1/_2).fastq.gz         — reads matching no barcode

    The barcode is stripped from R1 only (R2 doesn't carry it) when
    ``trim_barcode=True``.
    """
    validate_barcodes(barcodes, max_mismatches)
    paired = r2_path is not None

    def _open_pair(directory: Path) -> tuple:
        directory.mkdir(parents=True, exist_ok=True)
        if paired:
            w1 = gzip.open(directory / f"{srr}_1.fastq.gz", "wt", encoding="utf-8")
            w2 = gzip.open(directory / f"{srr}_2.fastq.gz", "wt", encoding="utf-8")
            return (w1, w2)
        return (gzip.open(directory / f"{srr}.fastq.gz", "wt", encoding="utf-8"), None)

    writers: dict = {}
    for rn in set(barcodes.values()):
        writers[rn] = _open_pair(out_dir / f"round_{rn:02d}")
    writers["_unassigned"] = _open_pair(out_dir / "unassigned")

    counts: Counter = Counter()
    total = 0
    try:
        if paired:
            assert r2_path is not None  # for type-checker
            for (h1, s1, p1, q1), (h2, s2, p2, q2) in _iter_fastq_pair(r1_path, r2_path):
                total += 1
                s1_stripped = s1.rstrip("\n")
                rn = _match_barcode(s1_stripped, barcodes, max_mismatches)
                key = "_unassigned" if rn is None else rn
                w1, w2 = writers[key]
                if rn is not None and trim_barcode:
                    bc_len = next(
                        len(bc) for bc in barcodes
                        if _hamming(s1_stripped[: len(bc)], bc) <= max_mismatches
                    )
                    w1.write(h1)
                    w1.write(s1_stripped[bc_len:] + "\n")
                    w1.write(p1)
                    w1.write(q1.rstrip("\n")[bc_len:] + "\n")
                else:
                    w1.write(h1); w1.write(s1); w1.write(p1); w1.write(q1)
                w2.write(h2); w2.write(s2); w2.write(p2); w2.write(q2)
                counts[key] += 1
        else:
            for header, seq, plus, qual in _iter_fastq(r1_path):
                total += 1
                seq_stripped = seq.rstrip("\n")
                rn = _match_barcode(seq_stripped, barcodes, max_mismatches)
                key = "_unassigned" if rn is None else rn
                w1, _ = writers[key]
                if rn is not None and trim_barcode:
                    bc_len = next(
                        len(bc) for bc in barcodes
                        if _hamming(seq_stripped[: len(bc)], bc) <= max_mismatches
                    )
                    w1.write(header)
                    w1.write(seq_stripped[bc_len:] + "\n")
                    w1.write(plus)
                    w1.write(qual.rstrip("\n")[bc_len:] + "\n")
                else:
                    w1.write(header); w1.write(seq); w1.write(plus); w1.write(qual)
                counts[key] += 1
    finally:
        for w1, w2 in writers.values():
            if w1 is not None:
                w1.close()
            if w2 is not None:
                w2.close()

    per_round = {int(k): v for k, v in counts.items() if k != "_unassigned"}
    return DemuxReport(
        srr=srr,
        input_r1=str(r1_path),
        input_r2=str(r2_path) if r2_path else None,
        paired=paired,
        total_reads=total,
        assigned_reads=sum(per_round.values()),
        unassigned_reads=counts["_unassigned"],
        per_round=per_round,
        unassigned_fraction=counts["_unassigned"] / total if total else 0.0,
        max_mismatches=max_mismatches,
        trim_barcode=trim_barcode,
    )


# ---------------------------------------------------------------------------
# Sample-sheet API
# ---------------------------------------------------------------------------


def read_sample_sheet(path: Path) -> list[DemuxJob]:
    """Parse a TSV sample sheet into ``DemuxJob`` records.

    Expected columns (header row):
        srr          — accession / run ID, used for output filenames
        r1_path      — path to single-end or R1 paired-end FASTQ.gz
        r2_path      — path to R2 paired-end FASTQ.gz (blank for SE)
        round        — integer round number this barcode tags
        barcode      — 5' DNA barcode sequence (ACGT only)

    Multiple rows per `srr` are merged into one job with a barcode→round dict.
    """
    rows: list[dict[str, str]] = []
    with open(path, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append({k.strip(): (v or "").strip() for k, v in row.items()})

    jobs: dict[str, DemuxJob] = {}
    for row in rows:
        srr = row.get("srr") or row.get("accession") or ""
        if not srr:
            continue
        r1 = Path(row["r1_path"])
        r2 = Path(row["r2_path"]) if row.get("r2_path") else None
        bc = row["barcode"].upper()
        rn = int(row["round"])

        if srr not in jobs:
            jobs[srr] = DemuxJob(srr=srr, r1_path=r1, r2_path=r2, barcodes={bc: rn})
        else:
            jobs[srr].barcodes[bc] = rn
    return list(jobs.values())


def demux_sample_sheet(
    sample_sheet: Path,
    out_root: Path,
    max_mismatches: int = 1,
    trim_barcode: bool = True,
    report_dir: Path | None = None,
) -> list[DemuxReport]:
    """Run ``demux_fastq`` over every job in a sample sheet.

    Per-job reports are returned and (if ``report_dir`` is given) written as
    ``{report_dir}/{srr}_demux.json``.
    """
    jobs = read_sample_sheet(sample_sheet)
    reports: list[DemuxReport] = []
    if report_dir is not None:
        report_dir.mkdir(parents=True, exist_ok=True)

    for job in jobs:
        report = demux_fastq(
            r1_path=job.r1_path,
            out_dir=out_root,
            barcodes=job.barcodes,
            srr=job.srr,
            r2_path=job.r2_path,
            max_mismatches=max_mismatches,
            trim_barcode=trim_barcode,
        )
        reports.append(report)
        if report_dir is not None:
            (report_dir / f"{job.srr}_demux.json").write_text(
                json.dumps(asdict(report), indent=2)
            )
    return reports
