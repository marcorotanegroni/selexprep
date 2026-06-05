# Examples

Normal use is a single command — see the [project README](../README.md#quick-start)
for the full quickstart. In short:

```bash
printf 'accession\nPRJNA615076\n' > accessions.tsv
selexprep run accessions.tsv --outdir out --resume
```

Then look in `out/<accession>/`: `library_report.json`, `round_NN/counts.parquet`,
`qc/flags.yaml`, plus the corpus-level `out/run_summary.tsv`.

## What's here

- **[`run_public_accession.sh`](run_public_accession.sh)** — the command above as a
  runnable script (`bash examples/run_public_accession.sh`). Fetches from ENA, so it
  needs network.
- **[`01_offline_toy_pipeline.ipynb`](01_offline_toy_pipeline.ipynb)** — offline
  reproducibility demo: the full `detect → extract → count → qc` pipeline on
  synthetic data, no network. A demo/verification, *not* a usage guide.
- **[`02_library_report_interpretation.ipynb`](02_library_report_interpretation.ipynb)**
  — how to read a `LibraryReport`, using cached real outputs under `data/`.
