#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00
#SBATCH --mem=8G
#SBATCH --output=round_coverage-%j.log
#SBATCH --partition=dgx5highpriority

# Round-assignment coverage of the released version (PROTOCOL.md, Amendment 1).
# Submit from the repository root, checked out at the release tag:
#
#   sbatch benchmarks/round_validation/coverage.sh
#
# Network only (ENA run and sample records). Writes
# benchmarks/round_validation/results/ (committed).
set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
# --frozen: use uv.lock as committed, without rewriting it, so the tree stays clean.
python3 -m uv run --frozen python benchmarks/round_validation/coverage.py \
    --out benchmarks/round_validation/results
echo FINITO
