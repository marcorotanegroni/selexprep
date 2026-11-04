#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --output=regenerate_reports-%j.log
#SBATCH --partition=dgx5highpriority

# Regenerates the cached LibraryReports that
# examples/02_library_report_interpretation.ipynb reads (examples/data/), with
# the checked-out version. Submit from the repository root:
#
#   sbatch examples/regenerate_reports.sh
#
# Uses the Tier-1 FASTQs under benchmarks/results/<accession>/ when they are
# there; otherwise fetches them first (network). detect is called as the
# Tier-1 Snakefile calls it (R1 positionally, R2 via --paired-r2, the fetched
# rounds.tsv).
set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
SP="python3 -m uv run --frozen selexprep"
$SP --version

for acc in PRJNA615076 PRJEB62495 PRJEB70964 PRJNA1395820 PRJNA809588; do
    src="benchmarks/results/$acc"
    if [ ! -s "$src/rounds.tsv" ]; then
        echo "== $acc: fetching"
        $SP fetch "$acc" --outdir "$src"
    fi
    find "$src" -name '*.fastq.gz' ! -name '*_2.fastq.gz' | sort > "$src/fastqs.manifest"
    find "$src" -name '*_2.fastq.gz' | sort > "$src/fastqs.r2.manifest"
    r2_args=()
    while IFS= read -r fq; do r2_args+=(--paired-r2 "$fq"); done < "$src/fastqs.r2.manifest"
    out="$(mktemp -d)"
    echo "== $acc: detect"
    # shellcheck disable=SC2046
    $SP detect $(cat "$src/fastqs.manifest") ${r2_args[@]+"${r2_args[@]}"} --round-map "$src/rounds.tsv" --outdir "$out"
    cp "$out/library_report.json" "examples/data/$acc.library_report.json"
    rm -rf "$out"
done
echo FINITO
