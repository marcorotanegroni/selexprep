#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --cpus-per-task=8
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCH --output=verify_read_state-%j.log
#SBATCH --partition=dgx5highpriority

# All-read verification of the benchmark read states (see the docstring of
# benchmarks/verify_read_state.py). Submit from the repository root:
#
#   sbatch benchmarks/verify_read_state.sh                          # every ground_truth row
#   sbatch benchmarks/verify_read_state.sh --accession PRJEB70964   # one deposit
#
# PY must point at an interpreter with numpy (the selexprep environment has it).
set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
PY="${PY:-python3}"
OUT="${OUT:-benchmarks/read_state_full.tsv}"

"$PY" benchmarks/verify_read_state.py \
    --results-dir benchmarks/results \
    --ground-truth benchmarks/ground_truth.tsv \
    --jobs "${SLURM_CPUS_PER_TASK:-8}" \
    --out "$OUT" \
    "$@"
