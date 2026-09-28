#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --cpus-per-task=2
#SBATCH --time=06:00:00
#SBATCH --mem=24G
#SBATCH --output=donor_checks-%j.log
#SBATCH --partition=dgx5highpriority

# Read-level donor checks and pool sampling for the semi-synthetic benchmark
# (see the docstring of check_and_sample_donors.py). Submit from the repository
# root:
#
#   sbatch benchmarks/synthetic/check_and_sample_donors.sh
#
# Writes benchmarks/synthetic/donor_checks.tsv (committed) and the sampled pools
# under benchmarks/synthetic/pools/ (kept on the HPC; regenerable from ENA).
set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
python3 -m uv run --extra bench python benchmarks/synthetic/check_and_sample_donors.py \
    --cache benchmarks/synthetic/fastq_cache \
    --pools benchmarks/synthetic/pools \
    --out benchmarks/synthetic/donor_checks.tsv
echo FINITO
