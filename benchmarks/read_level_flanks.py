"""Independent read-level flank consensus for recovery-arm ground-truth curation.

Phase 6b.10 Troncone 2 (non-circularity, Codex pass-4). The recovery-arm
``primer_*_truth`` in ``ground_truth.tsv`` must reflect the flanks PHYSICALLY
PRESENT in the deposited reads — and they must be derived from evidence that is
**independent of ``selexprep detect``** (otherwise the benchmark is circular:
testing detect against a truth detect itself produced).

This script is that independent evidence. It is deliberately self-contained —
**no ``selexprep`` import**, pure read inspection. For each accession it samples
the raw FASTQs from a run manifest (default ``fastqs.manifest`` = R1 /
single-end; optionally ``fastqs.r2.manifest`` = R2) and computes
position-anchored base consensus for the first ``--flank-len`` bases (5' end)
and the last ``--flank-len`` bases (3' end), with per-position support. The
constant region is where support stays high; it ends where support collapses
toward 0.25 (random region). The emitted consensus flanks + the provenance row
are what you transcribe into ``ground_truth.tsv`` (read-level orientation),
citing this table.

Worked context: PRJDB9110/9111 (RaptRanker, Ishida 2020) deposit the
reverse-complement / cDNA strand with a T7 promoter prepended, so the paper's
forward RNA-template constants don't match the read strand. Expect the 5'
consensus to start with the T7 promoter ``TAATACGACTCACTATAGGG`` and the 3'
consensus to be the reverse-complement of the paper's 5' constant.

Usage (run on HPC, where the benchmark FASTQs live):

    python benchmarks/read_level_flanks.py \\
        --results-dir benchmarks/results \\
        --accession PRJDB9110 --accession PRJDB9111 \\
        --sample 200000 --flank-len 45 \\
        --out benchmarks/read_level_truth_provenance.tsv

Paired-end 3' constant (e.g. PRJNA883192, whose 3' truth lives in R2): read
the R2 manifest instead, then reverse-complement the reported R2 5' consensus
to get the insert's 3' flank — keeps the R2-derived truth reproducible by the
same tool as the R1 rows, not hand-asserted:

    python benchmarks/read_level_flanks.py \\
        --results-dir benchmarks/results --accession PRJNA883192 \\
        --manifest-name fastqs.r2.manifest --sample 0 \\
        --out benchmarks/read_level_truth_provenance.tsv

The per-position tables print to stdout; the suggested flanks + support are
appended to the ``--out`` provenance TSV (created with a header if absent).
"""

from __future__ import annotations

import argparse
import gzip
import statistics
import sys
from collections import Counter
from pathlib import Path

_PROVENANCE_HEADER = [
    "accession",
    "n_reads_sampled",
    "flank_len",
    "support_threshold",
    "mode_read_len",
    "consensus_5p_flank",
    "support_5p_min",
    "consensus_3p_flank",
    "support_3p_min",
    "fastqs",
    "note",
]


def _iter_sequences(fastq: Path, limit: int | None):
    """Yield up to ``limit`` sequence lines from a (gzipped) FASTQ."""
    opener = gzip.open if str(fastq).endswith(".gz") else open
    n = 0
    with opener(fastq, "rt") as fh:  # type: ignore[operator]
        for i, line in enumerate(fh):
            if i % 4 == 1:  # the sequence line of each 4-line record
                yield line.strip()
                n += 1
                if limit is not None and n >= limit:
                    return


def _manifest_fastqs(
    results_dir: Path, accession: str, manifest_name: str = "fastqs.manifest"
) -> list[Path]:
    """FASTQs for an accession from the named run manifest (or a glob fallback).

    Default ``fastqs.manifest`` is the primary (R1 / single-end) stream. Pass
    ``manifest_name="fastqs.r2.manifest"`` to read the R2 mates — used to
    derive a paired-end 3' constant: the reverse-complement of the R2 5'
    consensus IS the insert's 3' flank (no detect involved). This keeps R2-
    derived truth reproducible by the same tool as the R1-derived rows.
    """
    manifest = results_dir / accession / manifest_name
    if manifest.exists():
        paths = [Path(ln.strip()) for ln in manifest.read_text().splitlines() if ln.strip()]
        if paths:
            return paths
    # Fallback glob: R2 mates for an r2 manifest, primary reads otherwise.
    want_r2 = "r2" in manifest_name.lower()
    return sorted(
        p
        for p in (results_dir / accession).glob("round_*/*.fastq.gz")
        if p.name.endswith("_2.fastq.gz") == want_r2
    )


def _consensus(counters: list[Counter]) -> tuple[str, list[float]]:
    """Per-position top base + support fraction over a list of base counters."""
    bases: list[str] = []
    supports: list[float] = []
    for c in counters:
        total = sum(c.values())
        if total == 0:
            break
        base, cnt = c.most_common(1)[0]
        bases.append(base)
        supports.append(cnt / total)
    return "".join(bases), supports


def _boundary(
    supports: list[float],
    threshold: float,
    *,
    baseline_window: int = 10,
    high_baseline: float = 0.90,
    min_drop: float = 0.12,
    post_max: float = 0.85,
) -> int:
    """Position of the constant→random support cliff (= end of constant).

    Two detectors, mirroring ``selexprep.detect`` so the recipe and detect
    agree on the boundary:

    1. **High-support cliff** — a sharp, sustained drop from a strong
       (≥``high_baseline``) constant baseline. This catches a clean constant
       whose post-constant tail is *low-complexity but still above the floor*
       (e.g. a C-rich tail at ~60%): a naive floor test would run past the
       real cliff into it (the PRJNA883192 R2 case: 96%→79% at pos 15, then
       a C-rich tail at 57-61%).
    2. **Floor + sustained-drop** fallback — for noisier constants (~80-86%,
       Illumina 3'-end decay) where no ≥0.90 baseline exists; the next-position
       guard avoids truncating on a single noisy dip.

    Set ``threshold`` below the constant band and above the random band
    (~0.25-0.35 undoped): default 0.6. The per-position table is
    authoritative; this is a heuristic.
    """
    # (1) high-support baseline cliff
    for i in range(1, len(supports)):
        prev = supports[max(0, i - baseline_window) : i]
        if not prev:
            continue
        baseline = statistics.median(prev)
        nxt = supports[i + 1] if i + 1 < len(supports) else 0.0
        if (
            baseline >= high_baseline
            and supports[i] <= baseline - min_drop
            and supports[i] < post_max
            and nxt < post_max
        ):
            return i
    # (2) floor + sustained-drop
    for i, s in enumerate(supports):
        if s < threshold and (i + 1 >= len(supports) or supports[i + 1] < threshold):
            return i
    return len(supports)


def analyze(
    results_dir: Path,
    accession: str,
    sample: int,
    flank_len: int,
    threshold: float,
    manifest_name: str = "fastqs.manifest",
) -> dict | None:
    fastqs = _manifest_fastqs(results_dir, accession, manifest_name)
    if not fastqs:
        print(f"[{accession}] no FASTQs found under {results_dir / accession}", file=sys.stderr)
        return None

    five = [Counter() for _ in range(flank_len)]
    three = [Counter() for _ in range(flank_len)]  # index 0 = last base (3' end)
    lengths: Counter = Counter()
    n = 0
    for fq in fastqs:
        if not fq.exists():
            print(f"[{accession}] missing {fq}", file=sys.stderr)
            continue
        remaining = None if sample == 0 else sample - n
        for seq in _iter_sequences(fq, remaining):
            length = len(seq)
            lengths[length] += 1
            for i in range(min(flank_len, length)):
                five[i][seq[i]] += 1
                three[i][seq[length - 1 - i]] += 1
            n += 1
        if sample and n >= sample:
            break

    cons5, sup5 = _consensus(five)
    cons3_from_end, sup3 = _consensus(three)  # read 3'->5'
    b5 = _boundary(sup5, threshold)
    b3 = _boundary(sup3, threshold)
    flank5 = cons5[:b5]
    flank3 = cons3_from_end[:b3][::-1]  # reverse to read 5'->3'

    return {
        "accession": accession,
        "manifest_name": manifest_name,
        "n_reads": n,
        "flank_len": flank_len,
        "threshold": threshold,
        "mode_read_len": lengths.most_common(1)[0][0] if lengths else 0,
        "consensus_5p_flank": flank5,
        "support_5p_min": min(sup5[:b5]) if b5 else 0.0,
        "consensus_3p_flank": flank3,
        "support_3p_min": min(sup3[:b3]) if b3 else 0.0,
        "fastqs": ";".join(p.name for p in fastqs),
        "cons5": cons5,
        "sup5": sup5,
        "cons3_from_end": cons3_from_end,
        "sup3": sup3,
        "b5": b5,
        "b3": b3,
    }


def _print_tables(r: dict) -> None:
    print(
        f"\n{'=' * 70}\n{r['accession']}  (n={r['n_reads']} reads, mode_len={r['mode_read_len']})"
    )
    print(f"{'-' * 70}\n5' end (pos 0 = read start):")
    for i, (b, s) in enumerate(zip(r["cons5"], r["sup5"], strict=True)):
        mark = "  <- constant" if i < r["b5"] else ("  <- random?" if i == r["b5"] else "")
        print(f"  pos {i:>2}: {b}  {s:6.1%}{mark}")
    print("3' end (pos 0 = read END, reading inward):")
    for i, (b, s) in enumerate(zip(r["cons3_from_end"], r["sup3"], strict=True)):
        mark = "  <- constant" if i < r["b3"] else ("  <- random?" if i == r["b3"] else "")
        print(f"  -{i + 1:>2}: {b}  {s:6.1%}{mark}")
    print(
        f"\n  => read-level 5' flank: {r['consensus_5p_flank']}  "
        f"(min support {r['support_5p_min']:.1%}, {r['b5']} nt)"
    )
    print(
        f"  => read-level 3' flank: {r['consensus_3p_flank']}  "
        f"(min support {r['support_3p_min']:.1%}, {r['b3']} nt)"
    )


def _append_provenance(out: Path, r: dict) -> None:
    new = not out.exists()
    with out.open("a", encoding="utf-8") as fh:
        if new:
            fh.write("\t".join(_PROVENANCE_HEADER) + "\n")
        fh.write(
            "\t".join(
                str(x)
                for x in [
                    r["accession"],
                    r["n_reads"],
                    r["flank_len"],
                    r["threshold"],
                    r["mode_read_len"],
                    r["consensus_5p_flank"],
                    f"{r['support_5p_min']:.4f}",
                    r["consensus_3p_flank"],
                    f"{r['support_3p_min']:.4f}",
                    r["fastqs"],
                    (
                        "independent read-level consensus (no detect); R2 mates → "
                        "insert 3' = revcomp(this 5' consensus)"
                        if "r2" in str(r.get("manifest_name", "")).lower()
                        else "independent read-level consensus (no detect); "
                        "transcribe into ground_truth.tsv"
                    ),
                ]
            )
            + "\n"
        )


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--results-dir", type=Path, default=Path("benchmarks/results"))
    p.add_argument("--accession", action="append", required=True, help="repeatable")
    p.add_argument("--sample", type=int, default=200000, help="reads to sample (0 = all)")
    p.add_argument("--flank-len", type=int, default=45)
    p.add_argument(
        "--support-threshold",
        type=float,
        default=0.6,
        help="constant/random cliff threshold; set below the constant band "
        "(~0.8-0.99) and above the random band (~0.25-0.35). Default 0.6.",
    )
    p.add_argument(
        "--manifest-name",
        default="fastqs.manifest",
        help="manifest under <results>/<acc>/ to read. Default 'fastqs.manifest' "
        "(R1/single-end). Use 'fastqs.r2.manifest' for a paired-end 3' constant: "
        "the insert 3' flank = revcomp of the reported R2 5' consensus (no detect).",
    )
    p.add_argument("--out", type=Path, default=None, help="provenance TSV (appended)")
    args = p.parse_args(argv)

    rc = 0
    for acc in args.accession:
        r = analyze(
            args.results_dir,
            acc,
            args.sample,
            args.flank_len,
            args.support_threshold,
            args.manifest_name,
        )
        if r is None:
            rc = 1
            continue
        _print_tables(r)
        if args.out is not None:
            _append_provenance(args.out, r)
    if args.out is not None:
        print(f"\nprovenance -> {args.out}")
    return rc


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
