# Quick start

Process a single public SELEX accession from raw reads to QC — without supplying
primers. Prefer a fully offline, hands-on version? See the [examples](examples.md):
`01_offline_toy_pipeline.ipynb` runs every stage on synthetic data with no network.

## Install

```bash
pip install selexprep        # pulls in cutadapt
```

Working from a clone? `uv pip install -e .` from the repo root.

## Preview, then run

```bash
# Preview metadata only — no download
selexprep inspect PRJNA615076

# One accession, end to end
printf "accession\nPRJNA615076\n" > accessions.tsv
selexprep run accessions.tsv --outdir out
```

`run` chains `fetch → detect → extract → count → qc` for each accession and
writes everything under `out/<accession>/`:

- `library_report.json` — inferred primers + `status` / `extraction_mode` /
  `required_action` (see the [LibraryReport reference](library-report.md)).
- `round_NN/extracted.fasta.gz` + `counts.parquet` — per round (round numbers are
  whatever the deposit used; PRJNA615076's are 0, 2, 3, …, 11).
- `selexprep_manifest.json` — reproducibility manifest (dependency versions,
  deterministic hashes, captured parameters).
- `qc/` — four diagnostic PNGs + `flags.yaml`.

## Running the stages yourself

`run` exists because per-deposit details vary — notably **paired-end** layout
(R1 passed positionally, R2 via `--paired-r2`) and per-round file naming. To
drive the individual commands (`detect`, `extract`, `count`, `qc`) on a small
local dataset, see [`01_offline_toy_pipeline.ipynb`](examples.md); for every flag and the
local `--round-map` contract, see the [CLI reference](cli.md).

## When inference is uncertain

If `detect` returns `status: UNABLE_TO_INFER`, `extract` refuses rather than
guessing — supply `--override-primer-5p` / `--override-primer-3p` to proceed.
See [Limitations](limitations.md) for what v0.1 does and does not do.
