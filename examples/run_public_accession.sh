#!/usr/bin/env bash
# Run selexprep end-to-end on one public SELEX accession.
# Fetches FASTQs from ENA, so it needs network access.
# Usage:  bash examples/run_public_accession.sh
set -euo pipefail

demo="$(dirname "$0")/demo_public"
mkdir -p "$demo"
printf 'accession\nPRJNA615076\n' >"$demo/accessions.tsv"

selexprep run "$demo/accessions.tsv" --outdir "$demo/out" --resume

echo
echo "Done. Key outputs to look at:"
echo "  $demo/out/PRJNA615076/library_report.json     # inferred primers + status"
echo "  $demo/out/PRJNA615076/round_*/counts.parquet  # per-round unique-sequence counts"
echo "  $demo/out/PRJNA615076/qc/flags.yaml           # depth-aware QC flags"
echo "  $demo/out/run_summary.tsv                     # corpus-level summary"
