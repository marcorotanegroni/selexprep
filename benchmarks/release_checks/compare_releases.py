"""Run detect -> extract -> count with an earlier release and with the working
tree on the same inputs, and compare every output.

Inputs: synthetic libraries generated here (forward, reverse, mixed orientation,
sequencing errors, paired split, random reads, correct, truncated, degenerate
and wrongly oriented primer overrides), the bundled demo, and optionally the
real DL-SELEX subsets shipped with selexprepR (``--cs-dir``).

Extracted reads (decompressed) and count tables must be identical; JSON outputs
are compared field by field after removing the two output roots, and every
differing field is printed. Logs and plots are not compared.

    uv run python benchmarks/release_checks/compare_releases.py v0.4.6 /tmp/cmp \
        --cs-dir ~/selexprepR/inst/extdata > benchmarks/release_checks/v0.4.6_vs_0.4.7.txt
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
import os
import random
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
PY = Path(sys.executable)
VERSIONS: dict[str, str] = {}

P5 = "GGGAGACAAGAATAAACGCTCAA"
P3 = "TTCGACAGGAGGCTCACAACAGGC"
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


def mutate(s: str, rate: float, rng: random.Random) -> str:
    out = []
    for b in s:
        r = rng.random()
        if r < rate:
            out.append(rng.choice([x for x in "ACGT" if x != b]))
        elif r < rate * 1.1:
            out.append("N")
        else:
            out.append(b)
    return "".join(out)


def write_fq(path: Path, seqs: list[str]) -> None:
    with gzip.open(path, "wt") as fh:
        for i, s in enumerate(seqs):
            fh.write(f"@r{i}\n{s}\n+\n{'I' * len(s)}\n")


def library(rng, n_reads, *, n_len=40, err=0.005, pool=None, enrich=0.0):
    seqs = []
    for _ in range(n_reads):
        if pool and rng.random() < enrich:
            core = rng.choice(pool)
        else:
            L = n_len + rng.choice([0, 0, 0, 0, -1, 1])
            core = "".join(rng.choice("ACGT") for _ in range(L))
        seqs.append(mutate(P5 + core + P3, err, rng))
    return seqs


def make_synthetic(root: Path) -> list[dict]:
    rng = random.Random(7)
    pool = ["".join(rng.choice("ACGT") for _ in range(40)) for _ in range(30)]
    cases = []

    def case(name, rounds, *, r2=None, override=None):
        d = root / "inputs" / name
        d.mkdir(parents=True, exist_ok=True)
        files, r2files = [], []
        for i, seqs in enumerate(rounds):
            f = d / f"{name}_r{i}.fastq.gz"
            write_fq(f, seqs)
            files.append(f)
            if r2 is not None:
                g = d / f"{name}_r{i}_R2.fastq.gz"
                write_fq(g, r2[i])
                r2files.append(g)
        rm = d / "rounds.tsv"
        rows = [f"{f.name}\t{i}\n" for i, f in enumerate(files)] + [
            f"{g.name}\t{i}\n" for i, g in enumerate(r2files)
        ]
        rm.write_text("file\tround_number\n" + "".join(rows))
        cases.append(
            dict(name=name, fastq=files, r2=r2files or None, round_map=rm, override=override)
        )

    fwd = [library(rng, 4000, pool=pool, enrich=e) for e in (0.0, 0.2, 0.5)]
    case("fwd", fwd)
    case("fwd_override_ok", fwd, override=(P5, P3))
    case("fwd_override_short5", fwd, override=(P5[:-4], None))
    case("fwd_override_degenerate", fwd, override=(P5[:-1] + "R", P3))
    rev = [[rc(s) for s in r] for r in fwd]
    case("rev", rev)
    case("rev_override", rev, override=(P5, None))
    mixed = [fwd[0], [rc(s) for s in fwd[1]], fwd[2]]
    case("mixed", mixed)
    noisy = [library(rng, 4000, err=0.03, pool=pool, enrich=e) for e in (0.0, 0.3)]
    case("noisy", noisy)
    # Paired: 100-nt N, 80-nt reads, no overlap of the constants.
    frags = [library(rng, 3000, n_len=100, pool=None) for _ in range(2)]
    r1 = [[f[:80] for f in fr] for fr in frags]
    r2 = [[rc(f)[:80] for f in fr] for fr in frags]
    case("paired", r1, r2=r2)
    case("paired_override", r1, r2=r2, override=(P5, P3))
    # Paired with overlap: 40-nt N, 80-nt reads.
    frags = [library(rng, 3000, n_len=40) for _ in range(2)]
    case(
        "paired_overlap",
        [[f[:80] for f in fr] for fr in frags],
        r2=[[rc(f)[:80] for f in fr] for fr in frags],
    )
    # Random reads: detect refuses.
    case(
        "random",
        [["".join(rng.choice("ACGT") for _ in range(87)) for _ in range(3000)] for _ in range(2)],
    )
    return cases


CS_DIR: Path | None = None


def real_cases(root: Path) -> list[dict]:
    cases = []
    demo = REPO / "examples" / "demo"
    files = [demo / "demo_R1.fastq.gz", demo / "demo_R3.fastq.gz"]
    cases.append(
        dict(name="demo", fastq=files, r2=None, round_map=demo / "rounds.tsv", override=None)
    )
    if CS_DIR is None:
        return cases
    ext = CS_DIR
    files = [ext / f"cs_round_0{i}.fastq.gz" for i in (3, 5, 7)]
    rm = root / "inputs" / "cs_rounds.tsv"
    rm.parent.mkdir(parents=True, exist_ok=True)
    rm.write_text(
        "file\tround_number\n"
        + "".join(f"{f.name}\t{r}\n" for f, r in zip(files, (3, 5, 7), strict=True))
    )
    cases.append(dict(name="cs", fastq=files, r2=None, round_map=rm, override=None))
    return cases


def run(version: str, args: list[str], log: Path) -> int:
    env = dict(os.environ, PYTHONPATH=VERSIONS[version])
    env["PATH"] = f"{PY.parent}{os.pathsep}{env['PATH']}"
    with log.open("w") as fh:
        return subprocess.run(
            [str(PY), "-m", "selexprep.cli", *args], stdout=fh, stderr=subprocess.STDOUT, env=env
        ).returncode


def run_case(case: dict, version: str, out: Path, report_from: Path | None) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    common = ["--round-map", str(case["round_map"])]
    r2 = []
    for f in case["r2"] or []:
        r2 += ["--paired-r2", str(f)]
    rcs = {}
    if report_from is None:
        rcs["detect"] = run(
            version,
            ["detect", *map(str, case["fastq"]), *common, *r2, "--outdir", str(out / "detect")],
            out / "detect.log",
        )
        report = out / "detect" / "library_report.json"
    else:
        report = report_from
    ov = []
    if case["override"]:
        p5, p3 = case["override"]
        if p5:
            ov += ["--override-primer-5p", p5]
        if p3:
            ov += ["--override-primer-3p", p3]
    if report.is_file():
        rcs["extract"] = run(
            version,
            [
                "extract",
                *map(str, case["fastq"]),
                "--library-report",
                str(report),
                *common,
                *r2,
                *ov,
                "--outdir",
                str(out / "extract"),
            ],
            out / "extract.log",
        )
    for fa in (
        sorted((out / "extract").glob("round_*/extracted.fasta.gz"))
        if (out / "extract").is_dir()
        else []
    ):
        rnd = fa.parent.name.split("_")[1]
        rcs[f"count_{rnd}"] = run(
            version,
            ["count", str(fa), "--round", rnd, "--outdir", str(out / "count")],
            out / f"count_{rnd}.log",
        )
    return rcs


def fasta_text(p: Path) -> str:
    with gzip.open(p, "rt") as fh:
        return fh.read()


def norm_json(obj, root_a: str, root_b: str):
    s = json.dumps(obj, sort_keys=True).replace(root_a, "ROOT").replace(root_b, "ROOT")
    return json.loads(s)


def diff_dicts(a, b, path=""):
    out = []
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            out += diff_dicts(a.get(k, "<missing>"), b.get(k, "<missing>"), f"{path}.{k}")
    elif isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        for i, (x, y) in enumerate(zip(a, b, strict=True)):
            out += diff_dicts(x, y, f"{path}[{i}]")
    elif a != b:
        out.append((path, a, b))
    return out


VOLATILE = (
    "created_at",
    "timestamp",
    "seconds",
    "elapsed",
    "selexprep_version",
    "version",
    "git",
    "duration",
    "date",
    "runtime",
)


def compare(case: dict, a: Path, b: Path) -> list[str]:
    msgs = []
    for p in sorted(a.rglob("*")):
        if p.is_dir() or p.suffix in (".log", ".png"):
            continue
        q = b / p.relative_to(a)
        rel = p.relative_to(a)
        if not q.exists():
            msgs.append(f"  MISSING in new: {rel}")
            continue
        if p.name.endswith(".fasta.gz"):
            if fasta_text(p) != fasta_text(q):
                msgs.append(f"  FASTA DIFFERS: {rel}")
        elif p.suffix == ".parquet":
            da, db = pd.read_parquet(p), pd.read_parquet(q)
            if not da.equals(db):
                msgs.append(f"  COUNTS DIFFER: {rel}")
        elif p.suffix in (".json", ".yaml", ".yml", ".tsv", ".txt"):
            if p.suffix == ".json":
                ja = norm_json(json.loads(p.read_text()), str(a), str(b))
                jb = norm_json(json.loads(q.read_text()), str(a), str(b))
                for path, x, y in diff_dicts(ja, jb):
                    if any(v in path.lower() for v in VOLATILE):
                        continue
                    msgs.append(f"  {rel}{path}: {str(x)[:90]!s} -> {str(y)[:90]!s}")
            else:
                ta = p.read_text().replace(str(a), "ROOT")
                tb = q.read_text().replace(str(b), "ROOT")
                if ta != tb:
                    msgs.append(f"  TEXT DIFFERS: {rel}")
        else:
            if p.read_bytes() != q.read_bytes():
                msgs.append(f"  BYTES DIFFER: {rel}")
    for q in sorted(b.rglob("*")):
        if q.is_file() and q.suffix not in (".log", ".png") and not (a / q.relative_to(b)).exists():
            msgs.append(f"  EXTRA in new: {q.relative_to(b)}")
    return msgs


def main() -> None:
    global CS_DIR
    ap = argparse.ArgumentParser()
    ap.add_argument("old_tag")
    ap.add_argument("workdir", type=Path)
    ap.add_argument("--cs-dir", type=Path)
    args = ap.parse_args()
    CS_DIR = args.cs_dir
    root = args.workdir
    if root.exists():
        shutil.rmtree(root)
    old_src = root / "old_src"
    old_src.mkdir(parents=True)
    tar = subprocess.run(
        ["git", "-C", str(REPO), "archive", args.old_tag, "src"], check=True, capture_output=True
    ).stdout
    with tarfile.open(fileobj=io.BytesIO(tar)) as tf:
        tf.extractall(old_src, filter="data")
    VERSIONS["old"] = str(old_src / "src")
    VERSIONS["new"] = str(REPO / "src")
    cases = make_synthetic(root) + real_cases(root)
    # Real data, also with override of the detected primers (filled in after detect).
    for case in cases:
        name = case["name"]
        res = {}
        for v in VERSIONS:
            res[v] = run_case(case, v, root / "runs" / name / v, None)
        same_rc = res["old"] == res["new"]
        msgs = compare(case, root / "runs" / name / "old", root / "runs" / name / "new")
        print(
            f"[{name}] rc old={res['old']} new={res['new']}{'' if same_rc else '  <-- RC DIFFER'}"
        )
        for m in msgs:
            print(m)
        if not msgs:
            print("  identical")
        if name in ("demo", "cs"):
            rep = json.loads(
                (root / "runs" / name / "old" / "detect" / "library_report.json").read_text()
            )
            if rep.get("primer_5p"):
                ov_case = dict(
                    case,
                    name=f"{name}_override_detected",
                    override=(rep["primer_5p"], rep.get("primer_3p")),
                )
                report = root / "runs" / name / "old" / "detect" / "library_report.json"
                r = {
                    v: run_case(ov_case, v, root / "runs" / ov_case["name"] / v, report)
                    for v in VERSIONS
                }
                msgs = compare(
                    ov_case,
                    root / "runs" / ov_case["name"] / "old",
                    root / "runs" / ov_case["name"] / "new",
                )
                print(f"[{ov_case['name']}] rc old={r['old']} new={r['new']}")
                for m in msgs:
                    print(m)
                if not msgs:
                    print("  identical")


if __name__ == "__main__":
    main()
