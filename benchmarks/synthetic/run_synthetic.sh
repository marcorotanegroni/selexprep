#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --output=synthetic_dev-%j.log
#SBATCH --partition=dgx5highpriority

# Semi-synthetic benchmark, development set (DESIGN.md, "Freezing and changes",
# step 4). Submit from the repository root after the donor pools exist:
#
#   sbatch benchmarks/synthetic/run_synthetic.sh
#
# Resumable: a run whose result.json exists is not repeated. Writes
# benchmarks/synthetic/runs/dev/results.tsv and manifest.json (committed); the
# per-run directories stay on the HPC. The test set is run once, after the
# freeze, with --set test and its own output directory.
set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
python3 -m uv run --extra bench python benchmarks/synthetic/run_synthetic.py \
    --set dev \
    --pools benchmarks/synthetic/pools \
    --out benchmarks/synthetic/runs/dev \
    --jobs "${SLURM_CPUS_PER_TASK:-1}"
echo FINITO
