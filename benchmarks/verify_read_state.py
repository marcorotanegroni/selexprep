"""All-read verification of the read state recorded in ``read_state_evidence.tsv``.

The read state of every benchmark deposit (raw full-length vs pre-trimmed) and
the presence of its published constants were first recorded from a ``zcat``
spot-check of the leading reads. That is how PRJEB70964 came to be filed as a
raw library carrying its constants: a full look shows every run is a 35-nt
random-region-only deposit. This script replaces the spot-check with a pass
over **every read of every FASTQ** that the benchmark fetched, so the recorded
read state is a measurement rather than an impression.

It is deliberately self-contained -- **no ``selexprep`` import** -- for the same
non-circularity reason as ``read_level_flanks.py``: the evidence that decides
which arm a deposit belongs to must not come from the tool being scored. The
published constants are read from ``ground_truth.tsv``; nothing is inferred.

Per FASTQ (and per accession/mate in aggregate) it reports:

* read count and the full length distribution (mode, share at mode, min,
  median, max);
* the fraction of reads that carry the published 5' constant at the read start
  and the published 3' constant at the read end (up to one substitution), the
  same for the reverse-complement orientation (which is also where the 3'
  constant sits on an R2 mate), and the fraction carrying either constant
  exactly anywhere in the read (catches offset constants);
* the fraction carrying the TruSeq Read 1 or Nextera adapter probe (either
  strand) anywhere;
* a positional profile of the first and last ``--profile-len`` bases: the
  consensus base and its support as one digit per position (``9`` = at least
  90% of reads agree, ``2`` = 20-29%, i.e. random).

Rows are per file, so a deposit whose runs differ (different constructs pooled
under one round, a run with a different read length) is visible as such.

Usage (on HPC, where the benchmark FASTQs live):

    python benchmarks/verify_read_state.py \\
        --results-dir benchmarks/results --accession PRJEB70964 \\
        --jobs 8 --out benchmarks/read_state_full.tsv

    # every deposit in ground_truth.tsv
    python benchmarks/verify_read_state.py --jobs 8 \\
        --out benchmarks/read_state_full.tsv
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
from collections.abc import Iterator
from multiprocessing import Pool
from pathlib import Path

import numpy as np

_SCRIPT_PATH = "benchmarks/verify_read_state.py"

# Same probes as the package's adapter list, restated here so this script does
# not import selexprep (see module docstring).
_TRUSEQ_R1 = "AGATCGGAAGAGC"
_NEXTERA = "CTGTCTCTTATACACATCT"

_BASES = b"ACGT"
_MAX_LEN = 4096
_CHUNK = 200_000

_HEADER = [
    "accession",
    "mate",
    "file",
    "n_reads",
    "len_mode",
    "frac_len_mode",
    "len_min",
    "len_median",
    "len_max",
    "const_5p_len",
    "const_3p_len",
    "frac_5p_at_start",
    "frac_3p_at_end",
    "frac_rc3p_at_start",
    "frac_rc5p_at_end",
    "frac_5p_anywhere",
    "frac_3p_anywhere",
    "frac_truseq_r1",
    "frac_nextera",
    "cons_5p",
    "support_5p",
    "cons_3p",
    "support_3p",
    "script_commit",
]

_COUNT_KEYS = [
    "5p_at_start",
    "3p_at_end",
    "rc3p_at_start",
    "rc5p_at_end",
    "5p_anywhere",
    "3p_anywhere",
    "truseq_r1",
    "nextera",
]


def _revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def _one_mismatch(seq: str) -> frozenset[str]:
    """``seq`` plus every sequence at Hamming distance 1 (N counts as a base)."""
    out = {seq}
    for i, base in enumerate(seq):
        for sub in "ACGTN":
            if sub != base:
                out.add(seq[:i] + sub + seq[i + 1 :])
    return frozenset(out)


def _git_commit() -> str:
    """Short commit of this script, ``-dirty`` if it differs from HEAD."""
    try:
        rev = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True,
            text=True,
            timeout=5,
            check=False,
        ).stdout.strip()
        if not rev:
            return "unknown"
        dirty = (
            subprocess.run(
                ["git", "diff", "--quiet", "HEAD", "--", _SCRIPT_PATH], check=False
            ).returncode
            != 0
        )
        return rev + ("-dirty" if dirty else "")
    except Exception:
        return "unknown"


def _iter_sequences(fastq: Path) -> Iterator[str]:
    """Every sequence line of a FASTQ; gzip is decompressed in a child process."""
    if str(fastq).endswith(".gz"):
        tool = shutil.which("pigz") or "gzip"
        proc = subprocess.Popen([tool, "-dc", str(fastq)], stdout=subprocess.PIPE)
        assert proc.stdout is not None
        try:
            for i, line in enumerate(proc.stdout):
                if i % 4 == 1:
                    yield line.rstrip().decode("ascii").upper()
        finally:
            proc.stdout.close()
            if proc.wait() != 0:
                raise RuntimeError(f"{tool} failed on {fastq}")
    else:
        with fastq.open() as fh:
            for i, line in enumerate(fh):
                if i % 4 == 1:
                    yield line.rstrip().upper()


def _profile_add(profile: np.ndarray, chunk: list[str], width: int) -> None:
    """Add per-position A/C/G/T counts of equal-width strings into ``profile``."""
    arr = np.frombuffer("".join(chunk).encode("ascii"), dtype=np.uint8).reshape(-1, width)
    for b, code in enumerate(_BASES):
        profile[b] += (arr == code).sum(axis=0)


def scan_fastq(task: tuple[str, str, str, str, str, int]) -> dict:
    """Stream one FASTQ end to end and return its raw counts."""
    accession, mate, path, c5, c3, width = task
    fastq = Path(path)
    set_5p_start = _one_mismatch(c5) if c5 else frozenset()
    set_3p_end = _one_mismatch(c3) if c3 else frozenset()
    rc3, rc5 = _revcomp(c3), _revcomp(c5)
    set_rc3_start = _one_mismatch(rc3) if c3 else frozenset()
    set_rc5_end = _one_mismatch(rc5) if c5 else frozenset()
    l5, l3 = len(c5), len(c3)
    truseq_rc, nextera_rc = _revcomp(_TRUSEQ_R1), _revcomp(_NEXTERA)

    lengths = [0] * (_MAX_LEN + 1)
    counts = dict.fromkeys(_COUNT_KEYS, 0)
    prof5 = np.zeros((4, width), dtype=np.int64)
    prof3 = np.zeros((4, width), dtype=np.int64)
    heads: list[str] = []
    tails: list[str] = []
    n = 0
    for seq in _iter_sequences(fastq):
        n += 1
        lengths[min(len(seq), _MAX_LEN)] += 1
        if l5:
            if seq[:l5] in set_5p_start:
                counts["5p_at_start"] += 1
            if seq[-l5:] in set_rc5_end:
                counts["rc5p_at_end"] += 1
            if c5 in seq:
                counts["5p_anywhere"] += 1
        if l3:
            if seq[-l3:] in set_3p_end:
                counts["3p_at_end"] += 1
            if seq[:l3] in set_rc3_start:
                counts["rc3p_at_start"] += 1
            if c3 in seq:
                counts["3p_anywhere"] += 1
        if _TRUSEQ_R1 in seq or truseq_rc in seq:
            counts["truseq_r1"] += 1
        if _NEXTERA in seq or nextera_rc in seq:
            counts["nextera"] += 1
        heads.append(seq[:width].ljust(width, "N"))
        tails.append(seq[-width:].rjust(width, "N"))
        if len(heads) >= _CHUNK:
            _profile_add(prof5, heads, width)
            _profile_add(prof3, tails, width)
            heads.clear()
            tails.clear()
    if heads:
        _profile_add(prof5, heads, width)
        _profile_add(prof3, tails, width)

    return {
        "accession": accession,
        "mate": mate,
        "file": fastq.name,
        "n_reads": n,
        "lengths": lengths,
        "counts": counts,
        "prof5": prof5,
        "prof3": prof3,
        "const_5p_len": l5,
        "const_3p_len": l3,
    }


def _profile_strings(profile: np.ndarray) -> tuple[str, str]:
    """Consensus base and one-digit support per position (``9`` = >=90%)."""
    totals = profile.sum(axis=0)
    cons, digits = [], []
    for j in range(profile.shape[1]):
        if totals[j] == 0:
            cons.append("N")
            digits.append("-")
            continue
        b = int(profile[:, j].argmax())
        cons.append(chr(_BASES[b]))
        digits.append(str(min(9, int(10 * profile[b, j] / totals[j]))))
    return "".join(cons), "".join(digits)


def _row(r: dict, commit: str) -> dict:
    n = r["n_reads"]
    lengths = r["lengths"]
    present = [length for length, c in enumerate(lengths) if c]
    if n:
        mode = max(range(len(lengths)), key=lengths.__getitem__)
        half, cum, median = (n + 1) // 2, 0, 0
        for length, c in enumerate(lengths):
            cum += c
            if cum >= half:
                median = length
                break
    cons5, sup5 = _profile_strings(r["prof5"])
    cons3, sup3 = _profile_strings(r["prof3"])

    def frac(key: str, needs: int) -> str:
        return f"{r['counts'][key] / n:.4f}" if n and needs else ""

    return {
        "accession": r["accession"],
        "mate": r["mate"],
        "file": r["file"],
        "n_reads": n,
        "len_mode": mode if n else "",
        "frac_len_mode": f"{lengths[mode] / n:.4f}" if n else "",
        "len_min": present[0] if present else "",
        "len_median": median if n else "",
        "len_max": present[-1] if present else "",
        "const_5p_len": r["const_5p_len"],
        "const_3p_len": r["const_3p_len"],
        "frac_5p_at_start": frac("5p_at_start", r["const_5p_len"]),
        "frac_3p_at_end": frac("3p_at_end", r["const_3p_len"]),
        "frac_rc3p_at_start": frac("rc3p_at_start", r["const_3p_len"]),
        "frac_rc5p_at_end": frac("rc5p_at_end", r["const_5p_len"]),
        "frac_5p_anywhere": frac("5p_anywhere", r["const_5p_len"]),
        "frac_3p_anywhere": frac("3p_anywhere", r["const_3p_len"]),
        "frac_truseq_r1": frac("truseq_r1", 1),
        "frac_nextera": frac("nextera", 1),
        "cons_5p": cons5,
        "support_5p": sup5,
        "cons_3p": cons3,
        "support_3p": sup3,
        "script_commit": commit,
    }


def _aggregate(parts: list[dict]) -> dict:
    """Sum the raw counts of several files into one ``ALL`` record."""
    first = parts[0]
    return {
        "accession": first["accession"],
        "mate": first["mate"],
        "file": "ALL",
        "n_reads": sum(p["n_reads"] for p in parts),
        "lengths": [sum(col) for col in zip(*(p["lengths"] for p in parts), strict=True)],
        "counts": {k: sum(p["counts"][k] for p in parts) for k in _COUNT_KEYS},
        "prof5": sum(p["prof5"] for p in parts),
        "prof3": sum(p["prof3"] for p in parts),
        "const_5p_len": first["const_5p_len"],
        "const_3p_len": first["const_3p_len"],
    }


def _read_ground_truth(path: Path) -> dict[str, tuple[str, str]]:
    with path.open(newline="") as fh:
        return {
            row["accession"]: (
                (row.get("primer_5p_truth") or "").strip().upper(),
                (row.get("primer_3p_truth") or "").strip().upper(),
            )
            for row in csv.DictReader(fh, delimiter="\t")
        }


def _resolve(entry: str, acc_dir: Path) -> Path:
    """A manifest entry as written (relative to the CWD) or under the accession dir."""
    p = Path(entry)
    if p.is_absolute() or p.exists():
        return p
    return acc_dir / p.name if (acc_dir / p.name).exists() else acc_dir / p


def _fastqs(results_dir: Path, accession: str) -> dict[str, list[Path]]:
    """R1/single-end and (if present) R2 FASTQs of an accession, by mate."""
    acc_dir = results_dir / accession
    out: dict[str, list[Path]] = {}
    for mate, manifest_name in (("R1", "fastqs.manifest"), ("R2", "fastqs.r2.manifest")):
        manifest = acc_dir / manifest_name
        if manifest.exists():
            entries = [ln.strip() for ln in manifest.read_text().splitlines() if ln.strip()]
            if entries:
                out[mate] = [_resolve(e, acc_dir) for e in entries]
    if not out:
        files = sorted(acc_dir.glob("round_*/*.fastq.gz"))
        r1 = [p for p in files if not p.name.endswith("_2.fastq.gz")]
        r2 = [p for p in files if p.name.endswith("_2.fastq.gz")]
        out = {k: v for k, v in (("R1", r1), ("R2", r2)) if v}
    return out


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--results-dir", type=Path, default=Path("benchmarks/results"))
    p.add_argument("--ground-truth", type=Path, default=Path("benchmarks/ground_truth.tsv"))
    p.add_argument(
        "--accession", action="append", help="repeatable; default = every ground_truth row"
    )
    p.add_argument("--profile-len", type=int, default=30)
    p.add_argument("--jobs", type=int, default=1, help="FASTQs scanned in parallel")
    p.add_argument("--out", type=Path, required=True, help="TSV, overwritten")
    args = p.parse_args(argv)

    truth = _read_ground_truth(args.ground_truth)
    accessions = args.accession or list(truth)
    commit = _git_commit()

    tasks = []
    rc = 0
    for acc in accessions:
        c5, c3 = truth.get(acc, ("", ""))
        by_mate = _fastqs(args.results_dir, acc)
        if not by_mate:
            print(f"[{acc}] no FASTQs under {args.results_dir / acc}", file=sys.stderr)
            rc = 1
            continue
        for mate, paths in by_mate.items():
            for path in paths:
                if not path.exists():
                    print(f"[{acc}] missing {path}", file=sys.stderr)
                    rc = 1
                    continue
                tasks.append((acc, mate, str(path), c5, c3, args.profile_len))

    with Pool(max(1, args.jobs)) as pool:
        results = pool.map(scan_fastq, tasks, chunksize=1)

    rows = []
    groups: dict[tuple[str, str], list[dict]] = {}
    for r in results:
        rows.append(_row(r, commit))
        groups.setdefault((r["accession"], r["mate"]), []).append(r)
    for (acc, mate), parts in groups.items():
        agg = _row(_aggregate(parts), commit)
        rows.append(agg)
        print(
            f"{acc} {mate}: {agg['n_reads']} reads in {len(parts)} files | "
            f"length mode {agg['len_mode']} ({agg['frac_len_mode']}), "
            f"range {agg['len_min']}-{agg['len_max']} | "
            f"5' at start {agg['frac_5p_at_start'] or 'n/a'}, "
            f"3' at end {agg['frac_3p_at_end'] or 'n/a'}, "
            f"TruSeq {agg['frac_truseq_r1']}\n"
            f"    5' {agg['cons_5p']}\n       {agg['support_5p']}\n"
            f"    3' {agg['cons_3p']}\n       {agg['support_3p']}"
        )

    with args.out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=_HEADER, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"\n{len(rows)} rows -> {args.out}")
    return rc


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
