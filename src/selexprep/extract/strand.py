"""Strand orientation handling for SELEX extraction.

Phase 2's ``library/detect.py`` diagnoses each library's orientation
(``FORWARD`` / ``REVERSE`` / ``MIXED``); Phase 3 acts on that diagnosis:

- ``FORWARD``: no rewriting; trim as-is.
- ``REVERSE``: revcomp every read in the input FASTQ before trimming.
  Typical case is a SELEX submission where R1 was generated from the
  antisense strand of the library.
- ``MIXED``: emit ``strand_report.tsv`` with per-round forward/reverse
  counts but do NOT auto-flip — guessing wrong on the majority would
  destroy correctly-oriented reads. cutadapt's ``--discard-untrimmed``
  naturally drops the misoriented minority.

Public API:

- :func:`detect_strand_distribution` — count forward/reverse/ambiguous
  reads in a sequence list.
- :func:`reorient_fastq_gz` — revcomp every record in a FASTQ.gz; output
  is deterministic (mtime=0 header) so the manifest SHA256 is stable.
- :func:`write_strand_report` — TSV aggregator across rounds.
"""

from __future__ import annotations

import gzip
import logging
from pathlib import Path

from selexprep._io import open_gzip_text_deterministic

logger = logging.getLogger(__name__)


# CALIBRATION-TODO: not in locked plan; Codex confirms after Phase 6
# benchmark recovery numbers. v0.1 emits aggregated per-round counts;
# v0.2 may emit per-read tags for downstream filtering tools.
STRAND_REPORT_PER_READ = False


# Watson-Crick complement that tolerates ``N`` in raw sequencer reads.
# Distinct from the strict ``library.adapters.reverse_complement`` (which
# rejects Ns by design, since it operates on curated primer strings).
_READ_COMPLEMENT = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
    "U": "A",
    "N": "N",
    "a": "t",
    "t": "a",
    "c": "g",
    "g": "c",
    "u": "a",
    "n": "n",
}


def _revcomp_read(seq: str) -> str:
    """Reverse-complement a read sequence. Tolerates N (N -> N)."""
    # Use a translation table where possible for speed on long reads.
    return "".join(_READ_COMPLEMENT.get(b, "N") for b in reversed(seq))


def _hamming_le1(a: str, b: str) -> bool:
    """Hamming distance <= 1 with early exit. Local copy to avoid the
    cross-package private import."""
    if len(a) != len(b):
        return False
    diffs = 0
    for x, y in zip(a, b, strict=False):
        if x != y:
            diffs += 1
            if diffs > 1:
                return False
    return True


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def detect_strand_distribution(
    seqs: list[str],
    primer_5p: str | None,
    primer_3p: str | None,
) -> dict[str, int]:
    """Count forward / reverse / ambiguous reads in ``seqs``.

    A read is **forward** if its 5'-end matches ``primer_5p`` (Hamming <=
    1); **reverse** if its 5'-end matches ``reverse_complement(primer_3p)``;
    **ambiguous** otherwise. Returns ``{"forward": N, "reverse": M,
    "ambiguous": K}`` — sums to ``len(seqs)``.

    When both primers are ``None`` the call returns all-ambiguous (there
    is nothing to test against).
    """
    rc_3p = _revcomp_read(primer_3p) if primer_3p else None
    counts = {"forward": 0, "reverse": 0, "ambiguous": 0}
    for s in seqs:
        matched_forward = (
            primer_5p is not None
            and len(s) >= len(primer_5p)
            and _hamming_le1(s[: len(primer_5p)], primer_5p)
        )
        if matched_forward:
            counts["forward"] += 1
            continue
        matched_reverse = (
            rc_3p is not None and len(s) >= len(rc_3p) and _hamming_le1(s[: len(rc_3p)], rc_3p)
        )
        if matched_reverse:
            counts["reverse"] += 1
        else:
            counts["ambiguous"] += 1
    return counts


def reorient_fastq_gz(input_path: Path, output_path: Path) -> int:
    """Read ``input_path`` (FASTQ.gz); revcomp every record; write
    ``output_path`` (deterministic FASTQ.gz with ``mtime=0`` header).

    Quality strings are reversed (the per-base quality of position N becomes
    the quality of the revcomp's position L-N-1; no complement on
    qualities themselves). Read names are kept verbatim — downstream tools
    rely on them for pair-sync.

    Returns the number of records rewritten.
    """
    n_records = 0
    with (
        gzip.open(input_path, "rt", encoding="utf-8", errors="replace") as fh,
        open_gzip_text_deterministic(output_path) as out,
    ):
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            if not (header and seq and plus and qual):
                # Fail loud rather than silently produce a partial output —
                # the project's "no silent miscalls" discipline applies here
                # too. Callers handle the exception (extract/runner.py
                # surfaces it via ``ExtractResult.skipped_reason``).
                raise ValueError(
                    f"reorient_fastq_gz: truncated FASTQ record at index "
                    f"{n_records} in {input_path}"
                )
            seq_stripped = seq.rstrip("\n")
            qual_stripped = qual.rstrip("\n")
            out.write(header)
            out.write(_revcomp_read(seq_stripped))
            out.write("\n")
            out.write(plus)
            out.write(qual_stripped[::-1])
            out.write("\n")
            n_records += 1
    return n_records


def write_strand_report(
    distributions_by_round: dict[int, dict[str, int]],
    path: Path,
) -> None:
    """Write per-round forward/reverse/ambiguous counts as a sorted TSV.

    Columns: ``round\\tforward\\treverse\\tambiguous``. Rounds are emitted
    in ascending numeric order. The file is plain UTF-8 (not gzipped); it
    is small and meant to be human-readable for QC review.
    """
    lines = ["round\tforward\treverse\tambiguous\n"]
    for r in sorted(distributions_by_round):
        d = distributions_by_round[r]
        lines.append(
            f"{r}\t{d.get('forward', 0)}\t{d.get('reverse', 0)}\t{d.get('ambiguous', 0)}\n"
        )
    path.write_text("".join(lines), encoding="utf-8")
