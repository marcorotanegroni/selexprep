#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --output=synthetic-%j.log
#SBATCH --partition=dgx5highpriority

# Semi-synthetic benchmark (DESIGN.md, "Freezing and changes"). Submit from the
# repository root after the donor pools exist; SET is dev (default) or test:
#
#   sbatch benchmarks/synthetic/run_synthetic.sh
#   sbatch --export=ALL,SET=test benchmarks/synthetic/run_synthetic.sh
#
# Resumable: a run whose result.json exists is not repeated. Writes
# benchmarks/synthetic/runs/$SET/results.tsv and manifest.json (committed); the
# per-run directories stay on the HPC. The test set is run once, after the freeze.
set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
# --frozen: use uv.lock as committed, without rewriting it, so the tree stays clean.
python3 -m uv run --frozen --extra bench python benchmarks/synthetic/run_synthetic.py \
    --set "${SET:-dev}" \
    --pools benchmarks/synthetic/pools \
    --out "benchmarks/synthetic/runs/${SET:-dev}" \
    --jobs "${SLURM_CPUS_PER_TASK:-1}"
echo FINITO
