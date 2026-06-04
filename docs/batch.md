# Batch processing

`selexprep run` processes many accessions in one pass — emitting per-dataset
outputs plus corpus-level QC — and is resumable.

## Input: an accessions TSV

```text
accession
PRJNA615076
PRJEB28411
```

One accession per line, under an `accession` header column.

## Run

```bash
selexprep run accessions.tsv --outdir out --resume
```

For each accession `selexprep` runs the full `fetch → detect → extract → count →
qc` chain (the round map is auto-populated from ENA/SRA metadata) and writes
`out/<accession>/`. `--resume` skips accessions already completed, so a killed
run restarts where it left off.

## Outputs

- `out/<accession>/` — per-dataset `library_report.json`, extracted FASTA(s),
  per-round counts, `selexprep_manifest.json`, and a `qc/` folder.
- Corpus-level QC (e.g. a last-round uniques histogram across datasets) in the
  batch root.

Each dataset carries its own `LibraryReport`, so a `HIGH` dataset and an
`UNABLE_TO_INFER` one can sit side by side in the same batch — the latter is
flagged, never silently miscalled. See the
[LibraryReport reference](library-report.md) for the fields.
