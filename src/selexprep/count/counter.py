"""Per-round SELEX sequence counting (FASTQ.gz → parquet).

One FASTQ.gz round → parquet with columns (sequence, reads, rank, rpm).

- Uses pyfastx when available (C extension, ~10x faster) and falls back to
  stdlib gzip+manual FASTQ parse otherwise.
- Primer-aware trimming via cutadapt CLI when primers are known:
  * Both 5' and 3' primers known → anchored LINKED adapter (`^P5...P3`) — read
    is kept only if BOTH primers match at the correct positions. Prevents
    mid-read false-trims that unanchored `-g` would allow.
  * Single-primer and no-primer modes are still supported.
  * `--minimum-length 15` drops empty/stub reads; `--no-indels` tightens match.
  * Paired-end reads (R1 + R2 siblings) are trimmed together, with primers
    reverse-complemented for R2. Only R1 is counted (standard SELEX practice).
- Multi-target BioProjects: when the raw dir contains per-target subdirs
  (e.g. `raw/<bp>/TNF/round_00/` and `raw/<bp>/IL-6/round_00/`), each target
  is processed into its own `processed/<bp>/<target>/` tree. Mixed layouts
  (top-level rounds AND target subdirs in the same BP) raise a hard error —
  the pipeline does not guess which belongs to which.
"""

from __future__ import annotations

import collections
import gzip
import logging
import math
import subprocess
import tempfile
from collections.abc import Iterator
from pathlib import Path

import numpy as np
import pandas as pd

from selexprep._common import resolve_cutadapt

logger = logging.getLogger(__name__)


_RESERVED_SUBDIRS = {"unassigned", "processed_external"}
_MIN_TRIMMED_LEN = 15
_MAX_TRIMMED_LEN = 200


# ---------------------------------------------------------------------------
# FASTQ reading
# ---------------------------------------------------------------------------


def _iter_sequences(fastq_gz: Path) -> Iterator[str]:
    try:
        import pyfastx

        fq = pyfastx.Fastq(str(fastq_gz), build_index=False)
        for _, seq, _ in fq:
            yield seq
    except ImportError:
        with gzip.open(fastq_gz, "rt", encoding="utf-8", errors="replace") as fh:
            while True:
                header = fh.readline()
                if not header:
                    break
                seq = fh.readline().strip()
                fh.readline()
                fh.readline()
                if seq:
                    yield seq


# ---------------------------------------------------------------------------
# Primer handling (reverse-complement for paired-end R2)
# ---------------------------------------------------------------------------


def reverse_complement(seq: str) -> str:
    """Reverse complement of a DNA/RNA sequence (U/T-tolerant)."""
    table = str.maketrans("ACGTUNacgtun", "TGCAANtgcaan")
    return seq.translate(table)[::-1]


def pair_sibling(r1: Path) -> Path | None:
    """Return the R2 sibling for an `_1.fastq.gz` file, or None if absent/unpaired."""
    name = r1.name
    if not name.endswith("_1.fastq.gz"):
        return None
    r2 = r1.parent / (name[: -len("_1.fastq.gz")] + "_2.fastq.gz")
    return r2 if r2.exists() else None


def _cutadapt_strand_args(
    primer_5p: str | None,
    primer_3p: str | None,
    strand: str,
) -> list[str]:
    """Build the -g/-a (forward) or -G/-A (reverse) args for one strand.

    Reverse strand uses reverse-complemented primers.
    """
    g_flag, a_flag = ("-g", "-a") if strand == "forward" else ("-G", "-A")
    if strand == "reverse":
        p5 = reverse_complement(primer_3p) if primer_3p else None
        p3 = reverse_complement(primer_5p) if primer_5p else None
    else:
        p5 = primer_5p
        p3 = primer_3p

    if p5 and p3:
        # Linked adapter, anchored at 5' end — requires BOTH primers, keeps only
        # the random region between them.
        return [g_flag, f"^{p5}...{p3}"]
    if p5:
        return [g_flag, f"^{p5}"]
    if p3:
        return [a_flag, p3]
    return []


def _trim_pair(
    r1: Path,
    r2: Path | None,
    primer_5p: str | None,
    primer_3p: str | None,
    tmp_dir: Path,
    cap_length: int | None = None,
) -> tuple[Path, Path | None, bool]:
    """Trim primers via cutadapt.

    Returns (trimmed_r1, trimmed_r2_or_None, was_trimmed). Falls back to the
    untrimmed inputs (and was_trimmed=False) if cutadapt is unavailable,
    primers are unknown, or trimming fails.
    """
    if not primer_5p and not primer_3p:
        return r1, r2, False
    cutadapt_exe = resolve_cutadapt()
    if cutadapt_exe is None:
        logger.warning("cutadapt not found — skipping trim for %s", r1.name)
        return r1, r2, False

    t1 = tmp_dir / f"{r1.stem}_trimmed.fastq.gz"
    cmd = [
        cutadapt_exe,
        "--quiet",
        "--discard-untrimmed",
        "--minimum-length",
        str(_MIN_TRIMMED_LEN),
        "--maximum-length",
        str(_MAX_TRIMMED_LEN),
        "--error-rate",
        "0.1",
        "--no-indels",
    ]
    cmd += _cutadapt_strand_args(primer_5p, primer_3p, "forward")

    # 5p-only mode: when reads are LONGER than expected post-strip random region
    # (i.e., a 3p constant residue + adapter readthrough survives the 5p-anchored
    # trim), cap R1 to the random_region_len. Without this cap, the parquet would
    # carry a zero-info constant tail and 3'-end sequencing errors on that tail
    # would inflate fake-uniqueness in the unique-sequence count.
    if primer_5p and not primer_3p and cap_length:
        cmd += ["--length", str(cap_length)]

    if r2 is not None:
        t2: Path | None = tmp_dir / f"{r2.stem}_trimmed.fastq.gz"
        cmd += _cutadapt_strand_args(primer_5p, primer_3p, "reverse")
        # When only one primer is known, the opposite read may not physically
        # contain the revcomp adapter (read too short to span the random region
        # back to the other constant). Cutadapt's default --pair-filter=any
        # would then discard the pair on the unknown side. Switch to filtering
        # on the read whose primer IS known.
        if primer_5p and not primer_3p:
            cmd += ["--pair-filter=first"]
        elif primer_3p and not primer_5p:
            cmd += ["--pair-filter=second"]
        cmd += ["-o", str(t1), "-p", str(t2), str(r1), str(r2)]
    else:
        t2 = None
        cmd += ["-o", str(t1), str(r1)]

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        return t1, t2, True
    except subprocess.CalledProcessError as e:
        logger.warning("cutadapt failed for %s: %s", r1.name, (e.stderr or "")[-300:])
        return r1, r2, False


# ---------------------------------------------------------------------------
# Counting
# ---------------------------------------------------------------------------


def _length_distribution_warning(df: pd.DataFrame, expected_len: int | None = None) -> bool:
    lengths = df["sequence"].str.len()
    stddev = float(lengths.std())
    if stddev > 5:
        return True
    if expected_len is not None:
        mean_len = float(lengths.mean())
        if abs(mean_len - expected_len) > 8:
            return True
    return False


def _counter_into(
    counter: collections.Counter,
    fastq_gz: Path,
    primer_5p: str | None,
    primer_3p: str | None,
    cap_length: int | None = None,
) -> bool:
    """Accumulate sequences from `fastq_gz` into `counter`.

    Returns True iff primer trimming was actually applied. R2 files are skipped
    here — they're picked up via their R1 sibling in paired-mode cutadapt.
    R1-only counting matches standard SELEX practice.
    """
    if fastq_gz.name.endswith("_2.fastq.gz"):
        return False

    r2 = pair_sibling(fastq_gz)
    with tempfile.TemporaryDirectory() as tmp:
        input_r1, _input_r2, was_trimmed = _trim_pair(
            fastq_gz,
            r2,
            primer_5p,
            primer_3p,
            Path(tmp),
            cap_length=cap_length,
        )
        for seq in _iter_sequences(input_r1):
            counter[seq] += 1
    return was_trimmed


def pool_diversity_stats(df: pd.DataFrame, n_reads: int) -> tuple[float, float]:
    """Return (shannon_entropy_bits, top1_frac).

    Shannon is a SELEX-progression signal — starts high in naive pools and
    collapses as winners dominate.
    """
    if n_reads <= 0 or df.empty:
        return 0.0, 0.0
    probs = df["reads"].to_numpy(dtype=np.float64) / n_reads
    probs = probs[probs > 0]
    if probs.size == 0:
        return 0.0, 0.0
    shannon = float(-np.sum(probs * np.log2(probs)))
    top1_frac = float(df["reads"].iloc[0] / n_reads)
    return shannon, top1_frac


def _counter_to_parquet(
    counter: collections.Counter,
    output_parquet: Path,
    expected_random_len: int | None,
) -> dict:
    n_reads = sum(counter.values())
    n_unique = len(counter)
    if n_reads == 0:
        return {
            "n_reads": 0,
            "n_unique": 0,
            "top_seq_reads": 0,
            "top_seq_rpm": 0.0,
            "top1_frac": 0.0,
            "shannon_entropy_bits": 0.0,
            "n_singletons": 0,
            "singleton_frac": 0.0,
            "n_seqs_ge_10": 0,
            "n_seqs_ge_100": 0,
            "length_distribution_warning": False,
        }

    df = pd.DataFrame(counter.items(), columns=["sequence", "reads"])
    df.sort_values("reads", ascending=False, inplace=True)
    df.reset_index(drop=True, inplace=True)
    df["rank"] = df.index + 1
    df["rpm"] = (df["reads"] / n_reads) * 1_000_000

    output_parquet.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(output_parquet, index=False, compression="zstd")

    shannon, top1_frac = pool_diversity_stats(df, n_reads)
    reads_arr = df["reads"].to_numpy()
    n_singletons = int((reads_arr == 1).sum())
    n_seqs_ge_10 = int((reads_arr >= 10).sum())
    n_seqs_ge_100 = int((reads_arr >= 100).sum())
    top_row = df.iloc[0]
    return {
        "n_reads": int(n_reads),
        "n_unique": int(n_unique),
        "top_seq_reads": int(top_row["reads"]),
        "top_seq_rpm": float(top_row["rpm"]),
        "top1_frac": top1_frac,
        "shannon_entropy_bits": shannon,
        "n_singletons": n_singletons,
        "singleton_frac": n_singletons / n_unique if n_unique else 0.0,
        "n_seqs_ge_10": n_seqs_ge_10,
        "n_seqs_ge_100": n_seqs_ge_100,
        "length_distribution_warning": _length_distribution_warning(df, expected_random_len),
    }


def _iter_fasta_sequences(fasta_gz: Path) -> Iterator[str]:
    """Yield sequence lines from a (gzipped) FASTA file.

    Phase-3 outputs are FASTA (cutadapt's ``--fasta`` mode); ``count_fasta``
    consumes them. Multi-line FASTAs are not produced by selexprep but are
    tolerated (concatenates wrapped lines per record).
    """
    current: list[str] = []
    with gzip.open(fasta_gz, "rt", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith(">"):
                if current:
                    yield "".join(current)
                    current = []
                continue
            if line:
                current.append(line)
        if current:
            yield "".join(current)


def count_fasta(
    fasta_gz: Path,
    output_parquet: Path,
    expected_random_len: int | None = None,
) -> dict:
    """Count unique sequences from a primer-stripped FASTA.gz -> per-round parquet.

     entry point: consumes the output of ``selexprep extract``
    . Output schema matches
     :func:`count_round`: columns ``sequence``, ``reads``, ``rank``, ``rpm``.

     Unlike ``count_round``, no cutadapt invocation - the input is already
     primer-stripped by ``extract``. Pure aggregation: read FASTA, build
     Counter, emit parquet.
    """
    counter: collections.Counter[str] = collections.Counter()
    for seq in _iter_fasta_sequences(fasta_gz):
        counter[seq] += 1
    return _counter_to_parquet(counter, output_parquet, expected_random_len)


def _iter_fastq_sequences(fastq_gz: Path) -> Iterator[str]:
    """Yield sequence lines (line 2 of each 4-line record) from a FASTQ(.gz).

    Companion to :func:`_iter_fasta_sequences` for the
    ``--from-pretrimmed-fastq`` opt-in path. Truncated records (incomplete
    final 4-block) raise ``ValueError`` rather than silently producing
    partial output - same policy as ``extract/strand.py:reorient_fastq_gz``.

    **v0.1 validation scope (deliberate; harden in v0.2 if this path
    becomes a publicly advertised workflow).** This iterator only checks
    record completeness (4 lines per record). It does NOT enforce:

    - Header line begins with ``@`` (a malformed FASTQ missing the
      header sigil would silently use whatever is on that line as the
      header, then yield the next line as the sequence).
    - Separator line begins with ``+``.
    - ``len(sequence) == len(quality_string)`` per record.

    These checks would catch hand-corrupted or wrongly-formatted FASTQ
    inputs. They're deferred because the ``--from-pretrimmed-fastq`` flag
    is explicitly a power-user opt-in: the user is asserting they've
    pre-trimmed the FASTQ with a real tool (AptaPLEX, EasyDIVER+,
    cutadapt), which by construction produces well-formed FASTQ. If
    selexprep promotes this path to a first-class workflow in v0.2, add
    those validators here (and surface them as ``ValueError`` matching
    the truncation policy).
    """
    opener = gzip.open if fastq_gz.suffix == ".gz" else open
    with opener(fastq_gz, "rt", encoding="utf-8", errors="replace") as fh:
        record_idx = 0
        while True:
            header = fh.readline()
            if not header:
                return
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            if not (header and seq and plus and qual):
                raise ValueError(
                    f"_iter_fastq_sequences: truncated FASTQ record at index "
                    f"{record_idx} in {fastq_gz}"
                )
            yield seq.rstrip("\n")
            record_idx += 1


def count_fastq_pretrimmed(
    fastq_gz: Path,
    output_parquet: Path,
    expected_random_len: int | None = None,
) -> dict:
    """Count unique sequences from a pre-trimmed FASTQ(.gz) -> per-round parquet.

    Opt-in path for users who have already primer-stripped their FASTQ
    via an external tool (AptaPLEX, EasyDIVER+, manual cutadapt, ...).
    The caller is responsible for the trimming state - selexprep cannot
    verify it. If primers/adapters remain in the input, unique-sequence
    counts will be inflated. See the CLI ``--from-pretrimmed-fastq``
    flag for the warning surfaced at invocation time.

    Output schema is identical to :func:`count_fasta` /
    :func:`count_round` (sequence, reads, rank, rpm).
    """
    counter: collections.Counter[str] = collections.Counter()
    for seq in _iter_fastq_sequences(fastq_gz):
        counter[seq] += 1
    return _counter_to_parquet(counter, output_parquet, expected_random_len)


def count_round(
    fastq_gz: Path,
    output_parquet: Path,
    primer_5p: str | None = None,
    primer_3p: str | None = None,
    expected_random_len: int | None = None,
) -> dict:
    """Count a single FASTQ.gz → parquet + stats dict."""
    logger.info("Counting %s", fastq_gz.name)
    counter: collections.Counter = collections.Counter()
    was_trimmed = _counter_into(
        counter,
        fastq_gz,
        primer_5p,
        primer_3p,
        cap_length=expected_random_len,
    )

    stats = _counter_to_parquet(counter, output_parquet, expected_random_len)
    stats["primer_trim_applied"] = bool(was_trimmed and stats["n_reads"] > 0)
    stats["primer_unknown"] = not (primer_5p or primer_3p)

    if stats["n_reads"] == 0:
        logger.warning("  No reads found in %s", fastq_gz)
    else:
        logger.info(
            "  %d reads → %d unique → %s",
            stats["n_reads"],
            stats["n_unique"],
            output_parquet.name,
        )
    if stats["length_distribution_warning"]:
        logger.warning(
            "  Suspicious length distribution in %s (possible adapter contamination)",
            fastq_gz.name,
        )
    return stats


def count_round_multi(
    fastq_files: list[Path],
    output_parquet: Path,
    primer_5p: str | None = None,
    primer_3p: str | None = None,
    expected_random_len: int | None = None,
) -> dict:
    """Merge counts from multiple FASTQs (same round, multi-SRR) in-memory."""
    logger.info("Counting %d FASTQs → %s", len(fastq_files), output_parquet.name)
    counter: collections.Counter = collections.Counter()
    any_trimmed = False
    for fq in fastq_files:
        if _counter_into(counter, fq, primer_5p, primer_3p, cap_length=expected_random_len):
            any_trimmed = True

    stats = _counter_to_parquet(counter, output_parquet, expected_random_len)
    stats["primer_trim_applied"] = bool(any_trimmed and stats["n_reads"] > 0)
    stats["primer_unknown"] = not (primer_5p or primer_3p)
    logger.info("  %d reads → %d unique (merged)", stats["n_reads"], stats["n_unique"])
    return stats


# ---------------------------------------------------------------------------
# Per-directory orchestrator (internal)
# ---------------------------------------------------------------------------


def _count_rounds_in_dir(
    raw_dir: Path,
    processed_dir: Path,
    primer_5p: str | None,
    primer_3p: str | None,
    expected_random_len: int | None,
) -> dict[int | str, dict]:
    """Count every `round_*/` subdir directly under `raw_dir` → `processed_dir`."""
    results: dict[int | str, dict] = {}
    round_dirs = sorted([d for d in raw_dir.glob("round_*") if d.is_dir()])

    for rd in round_dirs:
        label = rd.name
        try:
            rn: int | str = int(rd.name.split("_")[-1])
        except ValueError:
            rn = label

        # Only count R1 and single-end FASTQs; _2 siblings are handled by cutadapt
        # during paired-mode trimming. Prevents double-counting reads.
        fastq_files = [
            p for p in sorted(rd.glob("*.fastq.gz")) if not p.name.endswith("_2.fastq.gz")
        ]
        if not fastq_files:
            logger.warning("  No FASTQ files in %s", rd)
            continue

        out_parquet = (
            processed_dir / f"round_{rn:02d}.counts.parquet"
            if isinstance(rn, int)
            else processed_dir / f"{label}.counts.parquet"
        )

        if out_parquet.exists():
            logger.info("  [skip] %s already counted", out_parquet.name)
            df = pd.read_parquet(out_parquet, columns=["reads"])
            total = int(df["reads"].sum())
            reads_arr = df["reads"].to_numpy() if total else np.array([], dtype=int)
            n_unique = len(df)
            n_singletons = int((reads_arr == 1).sum()) if total else 0
            results[rn] = {
                "n_reads": total,
                "n_unique": n_unique,
                "top_seq_reads": int(df["reads"].max()) if total else 0,
                "top_seq_rpm": float(df["reads"].max() / total * 1e6) if total else 0.0,
                "top1_frac": float(df["reads"].max() / total) if total else 0.0,
                "shannon_entropy_bits": math.nan,
                "n_singletons": n_singletons,
                "singleton_frac": n_singletons / n_unique if n_unique else 0.0,
                "n_seqs_ge_10": int((reads_arr >= 10).sum()) if total else 0,
                "n_seqs_ge_100": int((reads_arr >= 100).sum()) if total else 0,
                "primer_trim_applied": None,
                "primer_unknown": not (primer_5p or primer_3p),
                "length_distribution_warning": False,
                "cached": True,
            }
            continue

        if len(fastq_files) == 1:
            results[rn] = count_round(
                fastq_files[0],
                out_parquet,
                primer_5p=primer_5p,
                primer_3p=primer_3p,
                expected_random_len=expected_random_len,
            )
        else:
            results[rn] = count_round_multi(
                fastq_files,
                out_parquet,
                primer_5p=primer_5p,
                primer_3p=primer_3p,
                expected_random_len=expected_random_len,
            )
    return results


# ---------------------------------------------------------------------------
# Public entry point — detects single- vs multi-target layout
# ---------------------------------------------------------------------------


def _detect_target_subdirs(bp_raw_dir: Path) -> list[Path]:
    """Return direct subdirs of `bp_raw_dir` that contain `round_*/` (per-target).

    Excludes reserved names (`unassigned`, `processed_external`) and top-level
    `round_*` dirs themselves.
    """
    if not bp_raw_dir.exists():
        return []
    targets = []
    for d in sorted(bp_raw_dir.iterdir()):
        if not d.is_dir():
            continue
        if d.name in _RESERVED_SUBDIRS:
            continue
        if d.name.startswith("round_"):
            continue
        if any(child.is_dir() and child.name.startswith("round_") for child in d.iterdir()):
            targets.append(d)
    return targets


def count_all_rounds(
    bp_raw_dir: Path,
    bp_processed_dir: Path,
    primer_5p: str | None = None,
    primer_3p: str | None = None,
    expected_random_len: int | None = None,
) -> dict:
    """Count every round for a BioProject.

    Return shape:
      - Single-target (flat layout, `raw/<bp>/round_NN/`):
          ``{round_number_or_label: stats_dict, ...}``
      - Multi-target (``raw/<bp>/<target>/round_NN/``):
          ``{"_multi_target": True, "targets": {target_name: {round: stats, ...}, ...}}``

    A mixed layout (both top-level rounds AND target subdirs) raises
    ``ValueError`` — the pipeline does not guess which belongs to which.
    """
    top_level_rounds = [d for d in bp_raw_dir.glob("round_*") if d.is_dir()]
    target_subdirs = _detect_target_subdirs(bp_raw_dir)

    if top_level_rounds and target_subdirs:
        raise ValueError(
            f"mixed round layout in {bp_raw_dir}: top-level round_* dirs "
            f"({[d.name for d in top_level_rounds]}) coexist with per-target "
            f"subdirs ({[d.name for d in target_subdirs]}). Restructure so "
            f"only one layout is used per BioProject."
        )

    if target_subdirs:
        logger.info(
            "multi-target layout detected in %s — processing %d targets: %s",
            bp_raw_dir,
            len(target_subdirs),
            [d.name for d in target_subdirs],
        )
        per_target: dict[str, dict[int | str, dict]] = {}
        for t_dir in target_subdirs:
            t_processed = bp_processed_dir / t_dir.name
            t_processed.mkdir(parents=True, exist_ok=True)
            per_target[t_dir.name] = _count_rounds_in_dir(
                t_dir,
                t_processed,
                primer_5p,
                primer_3p,
                expected_random_len,
            )
        return {"_multi_target": True, "targets": per_target}

    return _count_rounds_in_dir(
        bp_raw_dir,
        bp_processed_dir,
        primer_5p,
        primer_3p,
        expected_random_len,
    )
