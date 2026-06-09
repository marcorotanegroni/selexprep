# Examples

A suggested path for a first-time user, in order:

### 1. See it run — offline, no setup

[`01_offline_toy_pipeline.ipynb`](01_offline_toy_pipeline.ipynb) runs the full
`detect → extract → count → qc` pipeline on a tiny synthetic library — **no
network, no accession, no HPC**. It's the fastest way to see what each stage
consumes and produces, and the annotated walkthrough of the by-hand path you'd
use on your own local FASTQs.

> The demo library is *synthetic* (planted constants), so it proves the pipeline
> runs deterministically end-to-end — not how accurate inference is on messy real
> reads. For that evidence, see the Tier 1 scorecard under [`../benchmarks/`](../benchmarks/).

### 2. Run it for real — one public accession

The normal use is a single command. See the
[project quickstart](../README.md#quick-start); in short:

```bash
printf 'accession\nPRJNA615076\n' > accessions.tsv
selexprep run accessions.tsv --outdir out
```

Outputs land in `out/<accession>/`: `library_report.json`,
`round_NN/counts.parquet`, `qc/flags.yaml`, plus the corpus-level
`out/run_summary.tsv`. [`run_public_accession.sh`](run_public_accession.sh) is
exactly this as a runnable script (`bash examples/run_public_accession.sh`) —
it fetches from ENA, so it needs network.

### 3. Read the output — interpret a LibraryReport

[`02_library_report_interpretation.ipynb`](02_library_report_interpretation.ipynb)
walks through the `LibraryReport` fields (status, extraction_mode,
required_action, …) on cached real outputs under `data/` — high-confidence,
partial, and safe-failure cases. Read this once and the JSON in your own `out/`
folders stops being opaque.
