"""Read-level donor checks and pool sampling (DESIGN.md: donor screen steps 3-5,
"Which donor reads are used", amendment 2).

For every chosen run of every donor — the current donors in ``donors.tsv`` and
the candidates marked "to all-read check" in ``donor_screening.tsv`` — the run's
FASTQ is downloaded from ENA and every read is read once to record:

- the read count, the modal length and the fraction of reads within N +- 2 nt
  (the reads that are used; the rest are excluded before sampling). The
  mechanical part of step 3 — modal length equal to N and at least 70% of reads
  within N +- 2 — is reported as ``step3_reads``; the documentation part of step
  3 is decided outside this script;
- technical-residue signals at the read edges, tested at **every** length from 8
  nt up to the whole probe, anchored at the edge (amendment 2, item 6): the start
  of a known adapter (TruSeq R1/R2, Nextera, small-RNA 3') at the read end, and,
  where the donor's constants are documented, the part of the 5' constant next to
  N at the read start or the part of the 3' constant next to N at the read end.
  Each read contributes its longest match per edge, counted once per probe
  family (TruSeq R1 and R2 are one family). A signal in >= 5% of reads triggers
  an investigation, never an exclusion by itself; the diversity of the carrying
  reads (distinct / carrying) and the longest length seen are recorded. Which
  probe was used on each side is recorded too, and a constant side is marked
  ``unknown`` or ``too short`` when it was not probed, so "none" no longer stands
  for both "checked, nothing found" and "not checked" (item 7);
- descriptors, never used to exclude: distinct-read fraction, the share of the
  most frequent read, and the highest positional support among the first and
  last 10 positions of the reads within N +- 2.

It then samples 10,000 of the reads within N +- 2 without replacement, **once per
seed** (1-3 for development donors, 101-103 for test donors and candidates;
amendment 2), writes each pool one read per line and records its SHA-256.

Nothing here runs selexprep, ranks a donor or excludes one: admission is decided
afterwards and recorded in ``donor_screening.tsv``.

    python benchmarks/synthetic/check_and_sample_donors.py \\
        --cache benchmarks/synthetic/fastq_cache --pools benchmarks/synthetic/pools \\
        --out benchmarks/synthetic/donor_checks.tsv
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import random
import shutil
import sys
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
TOLERANCE = 2
SAMPLE_SIZE = 10_000
DEV_SEEDS = (1, 2, 3)
TEST_SEEDS = (101, 102, 103)
RESIDUE_K = 8
# A documented fragment shorter than the probe floor is still usable down to 6 nt,
# where a random match (1/4096) stays fifty times below the 5% trigger.
RESIDUE_K_MIN = 6
RESIDUE_FRACTION = 0.05
WITHIN_FRACTION = 0.70
PROFILE = 10
ADAPTERS = {
    "TruSeq_R1": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
    "TruSeq_R2": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
    "Nextera": "CTGTCTCTTATACACATCT",
    "smallRNA_3p": "TGGAATTCTCGGGTGCCAAGG",
}
# TruSeq R1 and R2 share their first 13 nt: a residue of 8-13 nt cannot tell them
# apart, so both count toward one signal. Counting them apart would split one
# residue between labels and could keep each part under the 5% trigger.
FAMILY = {"TruSeq_R1": "TruSeq", "TruSeq_R2": "TruSeq"}
COLUMNS = [
    "donor",
    "set",
    "role",
    "round",
    "run_accession",
    "n_stated",
    "reads",
    "modal_length",
    "within_n2",
    "within_n2_fraction",
    "excluded_fraction",
    "step3_reads",
    "probes_start",
    "probes_end",
    "residue_signals",
    "residue_max_length",
    "residue_carrier_diversity",
    "distinct_fraction",
    "top_read_share",
    "max_support_first10",
    "max_support_last10",
    "sample_seeds",
    "sample_sizes",
    "sample_sha256",
    "sample_dir",
]


class EdgeProbes:
    """Residue probes anchored at one read edge, matched at every length >= kmin.

    A read carrying a residue longer than the probe floor does not end with the
    floor-length probe — ``N + AGATCGGAA`` does not end with ``AGATCGGA`` — so
    every length is kept and indexed by the ``kmin`` nucleotides that sit at the
    edge. Each read gets its longest match, once.
    """

    def __init__(self, side: str) -> None:
        self.side = side
        self._labels: dict[str, set[str]] = defaultdict(set)
        self._index: dict[str, list[tuple[str, str]]] | None = None
        self._key_len = 0
        self.described: list[str] = []

    def add(self, name: str, sequence: str, kmin: int = RESIDUE_K) -> None:
        """``sequence`` is the adapter (read end) or the constant (read start).

        Matches are reported by probe family (``FAMILY``), so a read counts once
        toward one signal whichever adapter of the family its residue fits.
        """
        floor = min(kmin, len(sequence))
        if floor < RESIDUE_K_MIN:
            self.described.append(f"{name}={sequence}(too short, <{RESIDUE_K_MIN}nt)")
            return
        for k in range(floor, len(sequence) + 1):
            fragment = sequence[:k] if self.side == "end" else sequence[-k:]
            self._labels[fragment].add(FAMILY.get(name, name))
        # Show the part that sits at the edge: the adapter's start, the constant's end.
        if len(sequence) == floor:
            shown = sequence
        elif self.side == "end":
            shown = f"{sequence[:floor]}..."
        else:
            shown = f"...{sequence[-floor:]}"
        self.described.append(f"{name}={shown}(>={floor}nt)")
        self._index = None

    def build(self) -> None:
        # The index key is the shortest fragment length on this edge: a longer key
        # would never match a shorter probe.
        self._key_len = min(len(f) for f in self._labels) if self._labels else 0
        index: dict[str, list[tuple[str, str]]] = defaultdict(list)
        for fragment, names in self._labels.items():
            key = fragment[-self._key_len :] if self.side == "end" else fragment[: self._key_len]
            index[key].append((fragment, "/".join(sorted(names))))
        for entries in index.values():
            entries.sort(key=lambda e: -len(e[0]))
        self._index = index

    def match(self, seq: str) -> tuple[str, int] | None:
        if self._index is None:
            self.build()
        assert self._index is not None
        if not self._key_len or len(seq) < self._key_len:
            return None
        key = seq[-self._key_len :] if self.side == "end" else seq[: self._key_len]
        for fragment, label in self._index.get(key, ()):
            if seq.endswith(fragment) if self.side == "end" else seq.startswith(fragment):
                return label, len(fragment)
        return None


def chosen_runs() -> list[dict]:
    """(donor, set, N, role, round, run) for current donors and candidates."""
    runs = []
    for d in csv.DictReader((HERE / "donors.tsv").open(), delimiter="\t"):
        runs.append(
            dict(
                donor=d["donor"],
                set=d["set"],
                n=int(d["n_stated"]),
                role=d["role"],
                round=d["round"],
                run=d["run_accession"],
            )
        )
    for c in csv.DictReader((HERE / "donor_screening.tsv").open(), delimiter="\t"):
        if c["decision"] != "to all-read check":
            continue
        roles = ("earliest", "middle", "latest")
        for role, rnd, run in zip(
            roles, c["chosen_rounds"].split(","), c["chosen_runs"].split(","), strict=True
        ):
            runs.append(
                dict(
                    donor=c["accession"],
                    set="candidate",
                    n=int(c["n_stated"]),
                    role=role,
                    round=rnd,
                    run=run,
                )
            )
    return runs


def documented_constants() -> dict[str, dict[str, tuple[str, int]]]:
    """Per donor and side, the fragment next to N that documentation gives.

    The sides are independent (amendment 2, item 7): a deposit with one known
    side is probed on that side. ``ground_truth.tsv`` carries whole published
    constants; ``candidate_constants.tsv`` carries the shorter fragments the
    archive documents for the candidates, with their source.
    """
    out: dict[str, dict[str, tuple[str, int]]] = defaultdict(dict)
    with (REPO / "benchmarks" / "ground_truth.tsv").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            for side, column in (("5p", "primer_5p_truth"), ("3p", "primer_3p_truth")):
                if r[column]:
                    out[r["accession"]][side] = (r[column].upper(), RESIDUE_K)
    with (HERE / "candidate_constants.tsv").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if not r["fragment"]:
                continue
            fragment = r["fragment"].upper()
            out[r["accession"]].setdefault(r["side"], (fragment, len(fragment)))
    return dict(out)


def md5sum(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def download(run: str, cache: Path) -> Path:
    """The run's (R1 / single-end) FASTQ from ENA, cached and md5-checked."""
    url = (
        "https://www.ebi.ac.uk/ena/portal/api/filereport?accession="
        f"{run}&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv"
    )
    with urllib.request.urlopen(url, timeout=120) as resp:
        line = resp.read().decode().splitlines()[1].split("\t")
    ftp, md5 = line[1].split(";")[0], line[2].split(";")[0]
    target = cache / Path(ftp).name
    if not (target.exists() and md5sum(target) == md5):
        cache.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen("https://" + ftp, timeout=600) as resp, target.open("wb") as fh:
            shutil.copyfileobj(resp, fh)
        if md5sum(target) != md5:
            raise RuntimeError(f"md5 mismatch for {run}")
    return target


def sequences(path: Path):
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                yield line.rstrip().upper()


def read_key(seq: str) -> bytes:
    """A stable identity for a read: Python's str hash is salted per process."""
    return hashlib.blake2b(seq.encode(), digest_size=8).digest()


def max_support(profile: np.ndarray) -> float:
    totals = profile.sum(axis=0)
    ok = totals > 0
    return float((profile[:, ok].max(axis=0) / totals[ok]).max()) if ok.any() else 0.0


def check_and_sample(item: dict, constants: dict[str, tuple[str, int]], cache: Path, pools: Path):
    n = item["n"]
    start = EdgeProbes("start")
    end = EdgeProbes("end")
    for name, seq in ADAPTERS.items():
        end.add(name, seq)
    for side, probes in (("5p", start), ("3p", end)):
        if side in constants:
            fragment, kmin = constants[side]
            probes.add(f"constant_{side}", fragment, kmin=kmin)
        else:
            probes.described.append(f"constant_{side}=unknown")

    seeds = DEV_SEEDS if item["set"] == "dev" else TEST_SEEDS
    rngs = {
        seed: random.Random(
            int(
                hashlib.sha256(f"{item['donor']}:{item['run']}:{seed}".encode()).hexdigest()[:16],
                16,
            )
        )
        for seed in seeds
    }
    reservoirs: dict[int, list[str]] = {seed: [] for seed in seeds}
    lengths: Counter[int] = Counter()
    reads: Counter[bytes] = Counter()
    hits: Counter[str] = Counter()
    longest: Counter[str] = Counter()
    carriers: dict[str, set[bytes]] = defaultdict(set)
    head = np.zeros((4, PROFILE), dtype=np.int64)
    tail = np.zeros((4, PROFILE), dtype=np.int64)
    heads: list[str] = []
    tails: list[str] = []
    total = kept = 0
    for seq in sequences(download(item["run"], cache)):
        total += 1
        lengths[len(seq)] += 1
        key = read_key(seq)
        reads[key] += 1
        for probes in (start, end):
            found = probes.match(seq)
            if found:
                label, length = found
                hits[label] += 1
                longest[label] = max(longest[label], length)
                carriers[label].add(key)
        if abs(len(seq) - n) > TOLERANCE:
            continue
        kept += 1
        heads.append(seq[:PROFILE].ljust(PROFILE, "."))
        tails.append(seq[-PROFILE:].rjust(PROFILE, "."))
        if len(heads) >= 200_000:
            for prof, chunk in ((head, heads), (tail, tails)):
                arr = np.frombuffer("".join(chunk).encode(), dtype=np.uint8).reshape(-1, PROFILE)
                for b, code in enumerate(b"ACGT"):
                    prof[b] += (arr == code).sum(axis=0)
            heads.clear()
            tails.clear()
        for seed, reservoir in reservoirs.items():
            if len(reservoir) < SAMPLE_SIZE:
                reservoir.append(seq)
            else:
                j = rngs[seed].randrange(kept)
                if j < SAMPLE_SIZE:
                    reservoir[j] = seq
    for prof, chunk in ((head, heads), (tail, tails)):
        if chunk:
            arr = np.frombuffer("".join(chunk).encode(), dtype=np.uint8).reshape(-1, PROFILE)
            for b, code in enumerate(b"ACGT"):
                prof[b] += (arr == code).sum(axis=0)

    signals = {k: v / total for k, v in hits.items() if total and v / total >= RESIDUE_FRACTION}
    modal = lengths.most_common(1)[0][0] if lengths else None
    within = kept / total if total else 0.0
    failures = []
    if modal != n:
        failures.append(f"modal length {modal} != N {n}")
    if within < WITHIN_FRACTION:
        failures.append(f"{within:.1%} within N +- {TOLERANCE}")
    out_dir = pools / item["donor"]
    out_dir.mkdir(parents=True, exist_ok=True)
    digests = {}
    for seed, reservoir in reservoirs.items():
        body = "\n".join(reservoir) + "\n"
        name = f"{item['role']}_round{item['round']}_{item['run']}_seed{seed}.txt"
        (out_dir / name).write_text(body)
        digests[seed] = hashlib.sha256(body.encode()).hexdigest()
    return {
        "donor": item["donor"],
        "set": item["set"],
        "role": item["role"],
        "round": item["round"],
        "run_accession": item["run"],
        "n_stated": n,
        "reads": total,
        "modal_length": modal if modal is not None else "",
        "within_n2": kept,
        "within_n2_fraction": f"{within:.4f}",
        "excluded_fraction": f"{1 - within:.4f}" if total else "",
        "step3_reads": "pass" if not failures else "fail: " + "; ".join(failures),
        "probes_start": ";".join(start.described),
        "probes_end": ";".join(end.described),
        "residue_signals": ";".join(f"{k}={v:.4f}" for k, v in sorted(signals.items())) or "none",
        "residue_max_length": ";".join(f"{k}={longest[k]}" for k in sorted(signals)),
        "residue_carrier_diversity": ";".join(
            f"{k}={len(carriers[k]) / hits[k]:.4f}" for k in sorted(signals)
        ),
        "distinct_fraction": f"{len(reads) / total:.4f}" if total else "",
        "top_read_share": f"{reads.most_common(1)[0][1] / total:.4f}" if total else "",
        "max_support_first10": f"{max_support(head):.3f}",
        "max_support_last10": f"{max_support(tail):.3f}",
        "sample_seeds": ";".join(str(s) for s in seeds),
        "sample_sizes": ";".join(f"{s}={len(reservoirs[s])}" for s in seeds),
        "sample_sha256": ";".join(f"{s}={digests[s]}" for s in seeds),
        "sample_dir": str(out_dir.relative_to(REPO))
        if out_dir.is_relative_to(REPO)
        else str(out_dir),
    }


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--cache", type=Path, default=HERE / "fastq_cache")
    p.add_argument("--pools", type=Path, default=HERE / "pools")
    p.add_argument("--out", type=Path, default=HERE / "donor_checks.tsv")
    p.add_argument("--only", action="append", help="restrict to these donors (testing)")
    args = p.parse_args(argv)

    constants = documented_constants()
    rows = []
    for item in chosen_runs():
        if args.only and item["donor"] not in args.only:
            continue
        print(
            f"{item['donor']} {item['role']} round {item['round']} {item['run']}", file=sys.stderr
        )
        rows.append(
            check_and_sample(item, constants.get(item["donor"], {}), args.cache, args.pools)
        )
    with args.out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"{len(rows)} runs checked -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
