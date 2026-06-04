# CLI reference

`selexprep` is a Typer app. The core pipeline is
`inspect → fetch → detect → extract → count → qc`; `run` batches it and
`catalog` browses the bundled discovery catalog. Every command supports
`--help`; `selexprep --version` prints the version.

## `inspect` — preview an accession (no download)

```bash
selexprep inspect <ACCESSION> [--outdir DIR] [--timeout-s 30]
```

ENA/SRA/DDBJ metadata only: round count, `library_strategy`, BioProject, file
sizes. With `--outdir`, also writes `inspect.json`.

## `fetch` — download FASTQ + round map

```bash
selexprep fetch <ACCESSION> --outdir OUT [--backend ena|auto] [--allow-manual-review] [--dry-run] [--timeout-s 30]
```

Writes `OUT/rounds.tsv` (the trusted round→file map consumed by
`detect`/`extract`), per-round FASTQs under `OUT/round_NN/`, and
`OUT/fetch_metadata.json`. Default backend `ena` (fail-fast, reproducible);
`auto` adds the kingfisher → sra-toolkit fallback (GPL-3.0, opt-in). Runs whose
round assignment is NONE-confidence are refused unless `--allow-manual-review`
(then they land in `round_unknown/` and `manual_review.tsv`, never in
`rounds.tsv`).

## `detect` — infer primers / library structure

```bash
selexprep detect FASTQ... --round-map rounds.tsv --outdir OUT [--paired-r2 R2...] [--sampling-seed 42] [--max-reads-per-round N]
```

Emits `OUT/library_report.json`. `--round-map` is required for local FASTQs
(cross-round persistence is a core signal); for fetched data, pass the
`rounds.tsv` that `fetch` wrote. `--paired-r2` enables split-primer paired-end
detection. See the [LibraryReport reference](library-report.md).

## `extract` — trim constants + extract the random region

```bash
selexprep extract FASTQ... --library-report report.json --round-map rounds.tsv --outdir OUT [--paired-r2 R2...] [--override-primer-5p SEQ] [--override-primer-3p SEQ] [--rebuild]
```

Trims via `cutadapt` and writes `round_NN/extracted.fasta.gz` (full inserts) or
`partial_{5p,3p}_extracted*.fasta.gz` (single-side / split), plus
`selexprep_manifest.json`. If `detect` returned `UNABLE_TO_INFER`, `extract`
refuses unless you supply `--override-primer-*`. `--rebuild` overwrites and
emits `extract_diff.tsv`.

## `count` — per-round unique sequences

```bash
selexprep count EXTRACTED.fasta.gz --round R1 --outdir OUT [--from-pretrimmed-fastq]
```

Writes `OUT/round_NN/counts.parquet` (raw counts + RPM). `--from-pretrimmed-fastq`
opts into counting an externally-trimmed FASTQ (warns); by default FASTQ input is
rejected to keep the `extract → count` contract unambiguous.

## `qc` — plots + suspicion flags

```bash
selexprep qc MANIFEST.json [--counts-dir DIR] [--outdir DIR]
```

Reads `selexprep_manifest.json` + `round_*/counts.parquet`; writes four PNGs
(read retention, primer match per round, N-length distribution, per-round panel)
and a depth-aware `flags.yaml`.

## `run` — batch driver

```bash
selexprep run accessions.tsv --outdir OUT [--resume] [--stop-on-error] [--backend ena|auto]
```

Runs the full chain per accession (TSV with an `accession` column), emitting
`OUT/<accession>/…` plus corpus-level QC. `--resume` skips completed accessions;
by default per-accession failures are logged and skipped (`--stop-on-error` to
halt). See the [batch tutorial](batch.md).

## `catalog` — browse the bundled discovery catalog

```bash
selexprep catalog list                 # list public-SELEX bioprojects
selexprep catalog show <ACCESSION>     # full detail for one
selexprep catalog version              # catalog snapshot id
selexprep catalog refresh              # re-run the broad ENA discovery queries
```
