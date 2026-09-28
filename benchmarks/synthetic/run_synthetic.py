"""Harness of the semi-synthetic benchmark (DESIGN.md, "Configurations", "The
oracle", "Freezing and changes").

For every cell of the chosen set it generates the libraries (``simulate.py``),
runs ``selexprep detect`` and ``selexprep extract`` through the command line, as
a user would, extracts the same reads with the oracle report, and evaluates the
run (``evaluate.py``). A cell whose ``result.json`` exists is not run again, so
an interrupted job resumes. Bulky files (FASTQs, extracted reads) are deleted
after evaluation unless ``--keep``; they are regenerable from the seed.

    python benchmarks/synthetic/run_synthetic.py --set dev \\
        --pools benchmarks/synthetic/pools --out benchmarks/synthetic/runs/dev --jobs 8

Writes ``<out>/results.tsv`` (one row per run) and ``<out>/manifest.json`` (the
commit and state of the code that ran).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import time
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
TIMEOUT_S = 1800


def _module(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


simulate = _module("simulate")
evaluate = _module("evaluate")


def _env() -> dict[str, str]:
    # cutadapt and selexprep live next to the interpreter in the project's venv.
    env = dict(os.environ)
    env["PATH"] = f"{Path(sys.executable).parent}{os.pathsep}{env.get('PATH', '')}"
    return env


def _selexprep(args: list[str], log: Path) -> tuple[int, float]:
    t0 = time.monotonic()
    with log.open("w") as fh:
        try:
            rc = subprocess.run(
                [sys.executable, "-m", "selexprep.cli", *args],
                stdout=fh,
                stderr=subprocess.STDOUT,
                env=_env(),
                timeout=TIMEOUT_S,
            ).returncode
        except subprocess.TimeoutExpired:
            rc = -1
    return rc, round(time.monotonic() - t0, 2)


def write_oracle_report(oracle: dict, path: Path) -> None:
    """The configuration's truth as a ``LibraryReport`` (public schema)."""
    from selexprep.library.report import LibraryReport, write_library_report_json

    two_sided = oracle["mode"] == "BOTH_PRIMERS_SINGLE_READ"
    report = LibraryReport(
        primer_5p=oracle["primer_5p"],
        primer_3p=oracle["primer_3p"],
        variants_5p=[],
        variants_3p=[],
        known_adapter_hits={},
        extraction_mode=oracle["mode"],
        full_insert_recovered=two_sided,
        read_source="R1",
        required_action="NONE",
        orientation="FORWARD",
        n_length_mode=None,
        n_length_distribution={},
        n_length_confidence=1.0,
        match_rate_5p=1.0 if oracle["primer_5p"] else 0.0,
        match_rate_3p=1.0 if oracle["primer_3p"] else 0.0,
        position_consistency_5p=1.0 if oracle["primer_5p"] else 0.0,
        position_consistency_3p=1.0 if oracle["primer_3p"] else 0.0,
        read_fraction_used_for_inference=1.0,
        sampling_seed=0,
        confidence=1.0,
        status="HIGH",
        failure_reason=None,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    write_library_report_json(report, path)


def run_cell(cell, pools: Path, out: Path, keep: bool) -> dict:
    run_dir = out / "runs" / cell.cell_id
    result = run_dir / "result.json"
    if result.is_file():
        return json.loads(result.read_text())
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True)
    truth = simulate.generate(cell, pools, run_dir)
    fastqs = [str(run_dir / r["fastq"]) for r in truth["rounds"]]
    rounds = str(run_dir / "rounds.tsv")
    run: dict = {}

    rc, run["detect_seconds"] = _selexprep(
        [
            "detect",
            *fastqs,
            "--round-map",
            rounds,
            "--outdir",
            str(run_dir / "detect"),
        ],
        run_dir / "detect.log",
    )
    report_path = run_dir / "detect" / "library_report.json"
    if rc != 0 or not report_path.is_file():
        run["error"] = f"detect exited {rc}" if rc != -1 else "detect timed out"
    else:
        report = json.loads(report_path.read_text())
        refused = (
            report["status"] == "UNABLE_TO_INFER"
            or report["extraction_mode"] == "UNABLE_TO_EXTRACT"
        )
        if not refused:
            rc, run["extract_seconds"] = _selexprep(
                [
                    "extract",
                    *fastqs,
                    "--library-report",
                    str(report_path),
                    "--round-map",
                    rounds,
                    "--outdir",
                    str(run_dir / "extract"),
                ],
                run_dir / "extract.log",
            )
            if rc != 0:
                run["error"] = f"extract exited {rc}" if rc != -1 else "extract timed out"

    if truth["oracle"]["mode"] is not None:
        oracle_report = run_dir / "oracle" / "library_report.json"
        write_oracle_report(truth["oracle"], oracle_report)
        rc, run["oracle_seconds"] = _selexprep(
            [
                "extract",
                *fastqs,
                "--library-report",
                str(oracle_report),
                "--round-map",
                rounds,
                "--outdir",
                str(run_dir / "oracle" / "extract"),
            ],
            run_dir / "oracle.log",
        )
        if rc != 0:
            run["oracle_error"] = f"oracle extract exited {rc}"

    (run_dir / "run.json").write_text(json.dumps(run, indent=2, sort_keys=True) + "\n")
    row = evaluate.evaluate(run_dir)
    row.update({k: v for k, v in run.items() if k.endswith("_seconds")})
    result.write_text(json.dumps(row, indent=2, sort_keys=True) + "\n")
    if not keep:
        for path in run_dir.glob("round_*.fastq.gz"):
            path.unlink()
        for sub in (run_dir / "extract", run_dir / "oracle" / "extract"):
            shutil.rmtree(sub, ignore_errors=True)
    return row


def _star(args):
    return run_cell(*args)


def verify_pools(which: str, pools: Path) -> str:
    """Check every donor pool the set uses against the SHA-256 in ``donors.tsv``.

    Returns one digest over them; raises ``ValueError`` naming the first pool
    that is missing or differs. The pools are part of the experiment's identity:
    results made from other pools must not be resumed as these.
    """
    seeds = simulate.DEV_SEEDS if which == "dev" else simulate.TEST_SEEDS
    digests = []
    with (HERE / "donors.tsv").open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["set"] != which:
                continue
            expected = dict(pair.split("=") for pair in row["sample_sha256"].split(";"))
            for seed in seeds:
                path = simulate.pool_path(
                    pools, row["donor"], row["role"], int(row["round"]), row["run_accession"], seed
                )
                if not path.is_file():
                    raise ValueError(f"missing pool {path}")
                got = hashlib.sha256(path.read_bytes()).hexdigest()
                if got != expected.get(str(seed)):
                    raise ValueError(f"pool {path} differs from donors.tsv")
                digests.append(got)
    return hashlib.sha256("".join(digests).encode()).hexdigest()


def manifest(which: str, pools_sha256: str) -> dict:
    def git(*args: str) -> str:
        return subprocess.run(
            ["git", *args], cwd=REPO, capture_output=True, text=True, check=False
        ).stdout.strip()

    # The tracked files that differ from the commit, and a hash of how they
    # differ: the same file edited twice without a commit is another experiment.
    modified = [
        line[3:] for line in git("status", "--porcelain", "--untracked-files=no").splitlines()
    ]
    return {
        "set": which,
        "commit": git("rev-parse", "HEAD"),
        "dirty": bool(modified),
        "modified": modified,
        "modified_sha256": hashlib.sha256(git("diff", "HEAD").encode()).hexdigest(),
        "pools_sha256": pools_sha256,
        "python": sys.version.split()[0],
    }


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--set", choices=("dev", "test"), required=True)
    p.add_argument("--pools", type=Path, default=HERE / "pools")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--jobs", type=int, default=1)
    p.add_argument("--config", action="append", help="only these configurations")
    p.add_argument("--donor", action="append", help="only these donors")
    p.add_argument("--keep", action="store_true", help="keep FASTQs and extracted reads")
    args = p.parse_args(argv)

    todo = [
        c
        for c in simulate.cells(args.set)
        if (not args.config or c.config in args.config)
        and (not args.donor or c.donor in args.donor)
    ]
    try:
        pools_sha256 = verify_pools(args.set, args.pools)
    except ValueError as exc:
        print(f"donor pools: {exc}", file=sys.stderr)
        return 2
    current = manifest(args.set, pools_sha256)
    if args.set == "test" and current["dirty"]:
        # The test set runs once, at a recorded commit: nothing uncommitted.
        print(f"test set needs a clean tree; modified: {current['modified']}", file=sys.stderr)
        return 2
    args.out.mkdir(parents=True, exist_ok=True)
    stored_path = args.out / "manifest.json"
    if stored_path.is_file():
        # Resuming reuses every stored result.json, so it is allowed only for the
        # same experiment: same set, commit, uncommitted changes and pools.
        # Otherwise old results would be attributed to new code or data. The
        # first manifest is kept.
        stored = json.loads(stored_path.read_text())
        keys = ("set", "commit", "modified", "modified_sha256", "pools_sha256")
        if any(stored.get(k) != current[k] for k in keys):
            diff = {k: (stored.get(k), current[k]) for k in keys if stored.get(k) != current[k]}
            print(
                f"{args.out} holds results of another experiment {diff}; use a new --out",
                file=sys.stderr,
            )
            return 2
    else:
        stored_path.write_text(json.dumps(current, indent=2) + "\n")
    print(f"{len(todo)} runs ({args.set})", file=sys.stderr)
    jobs = [(c, args.pools, args.out, args.keep) for c in todo]
    rows = []
    with Pool(args.jobs) as pool:
        for i, row in enumerate(pool.imap_unordered(_star, jobs), 1):
            rows.append(row)
            print(f"[{i}/{len(jobs)}] {row['cell']}: {row['outcome']}", file=sys.stderr)
    rows.sort(key=lambda r: r["cell"])
    columns = list(dict.fromkeys(k for r in rows for k in r))
    with (args.out / "results.tsv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"{len(rows)} runs -> {args.out / 'results.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
