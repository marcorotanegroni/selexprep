# selexprep Tier 2 corpus-audit pipeline
#
# Scaffolding only: this Snakefile defines the DAG but is NOT executed
# in CI. The real-data HPC run that lands ``audit_results/`` artifacts
# (``audit_accessions.tsv`` + ``audit_metrics.json`` + ``table_audit.md``)
# is run separately on the HPC cluster.
#
# DAG (eligibility layer):
#
#   rule classify_catalog → eligibility.tsv (per-accession audit-eligibility
#                            classification via ENA fetch per row)
#   rule sample_corpus    → audit_accessions.tsv (samples only from
#                            ELIGIBLE_HT_SELEX_ROUNDS rows) + .manifest.json
#   rule run_corpus       → run_summary.tsv (via ``selexprep run --resume``)
#   rule aggregate_audit  → audit_metrics.json (includes
#                            catalog_classification_distribution from the
#                            full catalog classification)
#   rule figure_b         → table_audit.md
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

# Anchor paths to the Snakefile directory (path fix): otherwise
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
# Catalog snapshot path (resolved at import time from the installed
# package; the classifier reads the same file).
import selexprep.catalog.reader as _cat_reader
CATALOG_CSV = str(_cat_reader.catalog_path())


rule all:
    input:
        OUTROOT + "/table_audit.md",
        OUTROOT + "/audit_metrics.json",


rule classify_catalog:
    """classify every INSDC catalog row before sampling.

    Hits ENA once per accession to build a FetchPlan, then applies
    ``selexprep.benchmark.eligibility.classify_plan``. Only
    ``ELIGIBLE_HT_SELEX_ROUNDS`` rows feed ``sample_corpus`` below;
    the other buckets (NON_SELEX_ASSAY, NO_ROUND_STRUCTURE,
    MIXED_PROJECT_NEEDS_GROUPING, FETCH_DEAD) are counted and reported
    in ``audit_metrics.json``'s ``catalog_classification_distribution``.

    ~200 ENA queries; takes a couple of minutes. ``--limit`` available
    via config knob for smoke runs.
    """
    output:
        eligibility=OUTROOT + "/eligibility.tsv",
    params:
        catalog=CATALOG_CSV,
        limit=int(config.get("classify_limit", 0)),
    shell:
        "python -m selexprep.benchmark.eligibility classify-catalog "
        "--catalog {params.catalog} "
        "--out {output.eligibility} "
        + (" --limit {params.limit}" if int(config.get("classify_limit", 0)) > 0 else "")


rule sample_corpus:
    input:
        eligibility=OUTROOT + "/eligibility.tsv",
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
        "--eligibility {input.eligibility} "
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
        eligibility=OUTROOT + "/eligibility.tsv",
    output:
        OUTROOT + "/audit_metrics.json",
    params:
        # pass the catalog snapshot so the aggregator can
        # populate ``n_catalog_total`` + ``n_catalog_non_insdc_passthrough``
        # — surfacing the full-catalog denominator in the audit JSON +
        # audit-table caption (the eligibility classifier only sees INSDC
        # rows, so without this segment a reviewer reads "X of N
        # audit-eligible" as "X of all catalog rows").
        catalog=CATALOG_CSV,
    shell:
        "python -m selexprep.benchmark.corpus_audit aggregate "
        "--run-summary {input.summary} "
        "--ground-truth {input.ground_truth} "
        "--sample-manifest {input.manifest} "
        "--eligibility {input.eligibility} "
        "--catalog {params.catalog} "
        "--out {output}"


rule figure_b:
    input:
        OUTROOT + "/audit_metrics.json",
    output:
        OUTROOT + "/table_audit.md",
    shell:
        "python -m selexprep.benchmark.figure_b "
        "--audit {input} --outdir {OUTROOT}"
