#!/usr/bin/env bash
# Root owns resource checks and serial dependencies. No nested submissions.
#SBATCH --job-name=ground_certify
#SBATCH --account=acf-utk0011
#SBATCH --partition=short
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3500M
#SBATCH --time=00:15:00
set -euo pipefail
: "${VALIDATION_ROOT:?Set VALIDATION_ROOT}"
: "${CERTIFICATION_TRANCHE:?Set absolute staged tranche directory}"
(cd "$CERTIFICATION_TRANCHE" && sha256sum -c SHA256SUMS)
module load anaconda3/2024.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 JULIA_NUM_THREADS=1
export JULIA_DEPOT_PATH="$VALIDATION_ROOT/depot"
export RENYI_ROOT="$VALIDATION_ROOT/source/renyi_disjoint_corrected"
CERTIFICATION_OUTPUT="$VALIDATION_ROOT/results/ground_certification_${SLURM_JOB_ID}"
mkdir "$CERTIFICATION_OUTPUT"
/usr/bin/time -v -o "$CERTIFICATION_OUTPUT/oracle.resources.txt" \
 "$VALIDATION_ROOT/quspin-venv/bin/python" "$CERTIFICATION_TRANCHE/oracle_quspin.py" \
 --manifest "$CERTIFICATION_TRANCHE/manifest.toml" --output "$CERTIFICATION_OUTPUT/oracle.json"
/usr/bin/time -v -o "$CERTIFICATION_OUTPUT/target.resources.txt" \
 "$VALIDATION_ROOT/runtime/julia-1.12.1/bin/julia" --project="$RENYI_ROOT" \
 "$CERTIFICATION_TRANCHE/target.jl" --manifest "$CERTIFICATION_TRANCHE/manifest.toml" --output "$CERTIFICATION_OUTPUT/target.toml"
"$VALIDATION_ROOT/quspin-venv/bin/python" "$CERTIFICATION_TRANCHE/compare.py" \
 --manifest "$CERTIFICATION_TRANCHE/manifest.toml" --oracle "$CERTIFICATION_OUTPUT/oracle.json" \
 --target "$CERTIFICATION_OUTPUT/target.toml" --output "$CERTIFICATION_OUTPUT/comparison.json"
