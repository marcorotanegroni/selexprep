"""Forensic audits of SELEX library structure.

Two complementary audits:

- ``audit_raw_fastq`` — samples reads from a raw FASTQ.gz (pre-extraction) and
  characterises their length, the first 5 reads verbatim, and positional base
  frequencies at the 5' start and 3'-aligned end. Surfaces constant flanks
  that primer inference must account for.

- ``audit_trimmed_parquet`` — reads a per-round counts parquet (post-extraction)
  and reports length distribution, positional base frequency over the
  modal-length subset, top sequences, and Illumina TruSeq R1 contamination
  rate. A position with >40% top base is flagged as a constant-residue
  suspect (a truly random N-region should be 22-30% per base at every
  position).

Both audits return dataclasses; presentation (logging, HTML, terminal table) is
the caller's responsibility.
"""

from __future__ import annotations

import gzip
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

# Illumina TruSeq Read 1 adapter — used as a contamination probe in
# audit_trimmed_parquet. NOT a SELEX primer; presence in trimmed reads
# indicates extraction failure / under-trimming.
TRUSEQ_R1 = "AGATCGGAAGAGC"

_CONSTANT_SUSPECT_FRACTION = 0.40
_DEFAULT_SAMPLE_N = 2000


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------


@dataclass
class PositionalBaseFreq:
    """Base frequencies at a single position across a set of sequences."""

    position: int
    base_frequencies: dict[str, float] = field(default_factory=dict)
    top_base: str = "-"
    top_fraction: float = 0.0

    @property
    def constant_suspect(self) -> bool:
        return self.top_fraction > _CONSTANT_SUSPECT_FRACTION


@dataclass
class RawReadStructureAudit:
    """Audit of raw FASTQ reads sampled pre-extraction."""

    fastq_path: Path
    n_sampled: int
    mean_length: float
    first_5_reads: list[str]
    first_30nt_freq: list[PositionalBaseFreq]
    last_30nt_freq: list[PositionalBaseFreq]


@dataclass
class TrimmedRegionAudit:
    """Audit of a per-round counts parquet (post-extraction)."""

    parquet_path: Path
    n_unique: int
    n_reads: int
    length_distribution_top5: dict[int, int]
    mode_length: int
    first_10nt_freq: list[PositionalBaseFreq]
    last_10nt_freq: list[PositionalBaseFreq]
    top_5_sequences: list[tuple[str, int]]
    truseq_r1_unique_count: int
    truseq_r1_read_count: int


# ---------------------------------------------------------------------------
# Positional base frequency
# ---------------------------------------------------------------------------


def positional_base_freq(seqs: list[str], positions: range) -> list[PositionalBaseFreq]:
    """Compute base frequencies at each position in `positions` across `seqs`.

    Out-of-range positions are skipped silently per read. Frequencies are
    normalised to the number of reads contributing to that position.
    """
    counters: list[Counter] = [Counter() for _ in positions]
    for seq in seqs:
        for slot, pos in enumerate(positions):
            if 0 <= pos < len(seq):
                counters[slot][seq[pos]] += 1

    out: list[PositionalBaseFreq] = []
    for slot, pos in enumerate(positions):
        c = counters[slot]
        total = sum(c.values())
        if total == 0:
            out.append(PositionalBaseFreq(position=pos))
            continue
        freqs = {b: c[b] / total for b in "ACGTN" if c[b] > 0}
        top_base, top_count = c.most_common(1)[0]
        out.append(
            PositionalBaseFreq(
                position=pos,
                base_frequencies=freqs,
                top_base=top_base,
                top_fraction=top_count / total,
            )
        )
    return out


def _positional_base_freq_3prime(seqs: list[str], window: int) -> list[PositionalBaseFreq]:
    """3'-aligned positional base frequencies (handles variable-length reads).

    Position 0 corresponds to the last base of each read; position `window-1`
    is the `window`-th base from the end. Useful for spotting 3' constant
    flanks that variable-length sequencing reads would otherwise mask.
    """
    counters: list[Counter] = [Counter() for _ in range(window)]
    for seq in seqs:
        for i in range(min(window, len(seq))):
            counters[i][seq[len(seq) - 1 - i]] += 1

    out: list[PositionalBaseFreq] = []
    for i in range(window):
        c = counters[i]
        total = sum(c.values())
        if total == 0:
            out.append(PositionalBaseFreq(position=-(i + 1)))
            continue
        freqs = {b: c[b] / total for b in "ACGTN" if c[b] > 0}
        top_base, top_count = c.most_common(1)[0]
        out.append(
            PositionalBaseFreq(
                position=-(i + 1),
                base_frequencies=freqs,
                top_base=top_base,
                top_fraction=top_count / total,
            )
        )
    return out


# ---------------------------------------------------------------------------
# Raw FASTQ audit (pre-extraction)
# ---------------------------------------------------------------------------


def _sample_fastq(fastq_gz: Path, n: int) -> list[str]:
    """Read up to `n` sequences from a FASTQ.gz (stdlib gzip, no deps)."""
    seqs: list[str] = []
    with gzip.open(fastq_gz, "rt", encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                seqs.append(line.strip())
                if len(seqs) >= n:
                    break
    return seqs


def audit_raw_fastq(
    fastq_gz: Path,
    sample_n: int = _DEFAULT_SAMPLE_N,
    window: int = 30,
) -> RawReadStructureAudit:
    """Sample reads from a raw FASTQ.gz and characterise their structure.

    Returns 5' positional frequencies over the first `window` positions and
    3'-aligned frequencies over the last `window` positions. The 3'-aligned
    view handles variable-length reads correctly.
    """
    seqs = _sample_fastq(fastq_gz, sample_n)
    if not seqs:
        return RawReadStructureAudit(
            fastq_path=fastq_gz,
            n_sampled=0,
            mean_length=0.0,
            first_5_reads=[],
            first_30nt_freq=[],
            last_30nt_freq=[],
        )

    return RawReadStructureAudit(
        fastq_path=fastq_gz,
        n_sampled=len(seqs),
        mean_length=sum(len(s) for s in seqs) / len(seqs),
        first_5_reads=seqs[:5],
        first_30nt_freq=positional_base_freq(seqs, range(window)),
        last_30nt_freq=_positional_base_freq_3prime(seqs, window),
    )


# ---------------------------------------------------------------------------
# Trimmed parquet audit (post-extraction)
# ---------------------------------------------------------------------------


def audit_trimmed_parquet(parquet_path: Path) -> TrimmedRegionAudit:
    """Audit a per-round counts parquet for length, base composition, and
    TruSeq R1 contamination.

    Expects columns `sequence` and `reads`. Empty parquets are handled by
    returning a zero-filled audit (mode_length=0, empty lists).
    """
    df = pd.read_parquet(parquet_path, columns=["sequence", "reads"])
    if df.empty:
        return TrimmedRegionAudit(
            parquet_path=parquet_path,
            n_unique=0,
            n_reads=0,
            length_distribution_top5={},
            mode_length=0,
            first_10nt_freq=[],
            last_10nt_freq=[],
            top_5_sequences=[],
            truseq_r1_unique_count=0,
            truseq_r1_read_count=0,
        )

    df["L"] = df["sequence"].str.len()
    total_reads = int(df["reads"].sum())
    L_dist = df.groupby("L")["reads"].sum().sort_values(ascending=False).head(5)
    mode_L = int(L_dist.index[0])

    mode_seqs = df.loc[df["L"] == mode_L, "sequence"].tolist()
    first10 = positional_base_freq(mode_seqs, range(min(10, mode_L)))
    last10 = positional_base_freq(mode_seqs, range(max(0, mode_L - 10), mode_L))

    top5 = df.nlargest(5, "reads")[["sequence", "reads"]]
    top_5_sequences = [(str(row["sequence"]), int(row["reads"])) for _, row in top5.iterrows()]

    truseq_mask = df["sequence"].str.contains(TRUSEQ_R1, regex=False)
    truseq_unique = int(truseq_mask.sum())
    truseq_reads = int(df.loc[truseq_mask, "reads"].sum())

    return TrimmedRegionAudit(
        parquet_path=parquet_path,
        n_unique=len(df),
        n_reads=total_reads,
        length_distribution_top5={int(L): int(r) for L, r in L_dist.items()},
        mode_length=mode_L,
        first_10nt_freq=first10,
        last_10nt_freq=last10,
        top_5_sequences=top_5_sequences,
        truseq_r1_unique_count=truseq_unique,
        truseq_r1_read_count=truseq_reads,
    )
