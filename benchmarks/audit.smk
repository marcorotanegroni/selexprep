# selexprep Phase 6b.3a — Tier 2 corpus-audit pipeline
#
# Scaffolding only: this Snakefile defines the DAG but is NOT executed
# in CI. The real-data HPC run that lands ``audit_results/`` artifacts
# (``audit_accessions.tsv`` + ``audit_metrics.json`` + ``figure_b.{pdf,png}``)
# is the 6b.4 follow-up commit with no code changes.
#
# DAG (locked plan + Codex peer-review):
#
#   rule sample_corpus  → audit_accessions.tsv (+ .manifest.json sidecar)
#   rule run_corpus     → run_summary.tsv (via ``selexprep run --resume``)
#   rule aggregate_audit → audit_metrics.json
#   rule figure_b       → figure_b.{pdf,png}
#
# Methodological correction folded into ``rule aggregate_audit``:
# ``inference_safe_failure_rate`` is computed ONLY among rows with a
# LibraryReport (denominator = ``n_with_library_report``), NEVER mixed
# with fetch failures. The aggregator enforces this; the Snakefile
# just passes the reproducibility envelope (catalog version + seed +
# sample sha256) through the sidecar manifest emitted by
# ``rule sample_corpus``.
#
# Reproducibility:
#
#   snakemake -s benchmarks/audit.smk --cores 1
#
# Or dry-run:
#
#   snakemake -s benchmarks/audit.smk --cores 1 --dry-run
#
# Configuration knobs (override at the CLI via ``--config key=value``):
#
#   --config n_sample=30      sample size (default 30)
#   --config seed=42          rng seed (default 42)

from pathlib import Path

# Anchor paths to the Snakefile directory (Codex peer-review fix): otherwise
# ``snakemake -s benchmarks/audit.smk`` from the repo root looks for
# ``ground_truth.tsv`` in the repo root and dies with a MissingInputException.
# ``workflow.basedir`` is the standard Snakemake idiom for "the directory of
# the main Snakefile" — stable across snakemake 6+. Output paths are absolute
# so the run lands in ``benchmarks/audit_results/`` regardless of cwd.
_BENCHMARKS_DIR = Path(workflow.basedir).resolve()
OUTROOT = str(_BENCHMARKS_DIR / "audit_results")
GROUND_TRUTH = str(_BENCHMARKS_DIR / "ground_truth.tsv")
N_SAMPLE = int(config.get("n_sample", 30))
SEED = int(config.get("seed", 42))


rule all:
    input:
        OUTROOT + "/figure_b.pdf",
        OUTROOT + "/figure_b.png",
        OUTROOT + "/audit_metrics.json",


rule sample_corpus:
    output:
        accessions=OUTROOT + "/audit_accessions.tsv",
        manifest=OUTROOT + "/audit_accessions.manifest.json",
    params:
        n=N_SAMPLE,
        seed=SEED,
        ground_truth=GROUND_TRUTH,
    shell:
        "python -m selexprep.benchmark.corpus_audit sample "
        "--n {params.n} --seed {params.seed} "
        "--exclude-ground-truth {params.ground_truth} "
        "--out {output.accessions}"


rule run_corpus:
    input:
        accessions=OUTROOT + "/audit_accessions.tsv",
    output:
        summary=OUTROOT + "/run_summary.tsv",
    params:
        runs_dir=OUTROOT + "/runs",
    shell:
        # --resume so a re-run picks up where the previous one left off
        # (matches the audit's "interruptible HPC job" use case).
        "selexprep run {input.accessions} --outdir {params.runs_dir} --resume "
        "&& cp {params.runs_dir}/run_summary.tsv {output.summary}"


rule aggregate_audit:
    input:
        summary=OUTROOT + "/run_summary.tsv",
        ground_truth=GROUND_TRUTH,
        manifest=OUTROOT + "/audit_accessions.manifest.json",
    output:
        OUTROOT + "/audit_metrics.json",
    shell:
        "python -m selexprep.benchmark.corpus_audit aggregate "
        "--run-summary {input.summary} "
        "--ground-truth {input.ground_truth} "
        "--sample-manifest {input.manifest} "
        "--out {output}"


rule figure_b:
    input:
        OUTROOT + "/audit_metrics.json",
    output:
        pdf=OUTROOT + "/figure_b.pdf",
        png=OUTROOT + "/figure_b.png",
    shell:
        "python -m selexprep.benchmark.figure_b "
        "--audit {input} --outdir {OUTROOT}"
