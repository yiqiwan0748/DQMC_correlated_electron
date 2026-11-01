#!/usr/bin/env bash
# Root owns submission/serial dependencies. No nested submissions.
#SBATCH --job-name=spectral_ed
#SBATCH --account=acf-utk0011
#SBATCH --partition=short
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3500M
#SBATCH --time=00:10:00
set -euo pipefail
: "${VALIDATION_ROOT:?Set VALIDATION_ROOT}"
: "${SPECTRAL_TRANCHE:?Set absolute staged tranche directory}"
module load anaconda3/2024.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 JULIA_NUM_THREADS=1
export JULIA_DEPOT_PATH="$VALIDATION_ROOT/depot"
export RENYI_ROOT="$VALIDATION_ROOT/source/renyi_disjoint_corrected"
SPECTRAL_OUTPUT="$VALIDATION_ROOT/results/spectral_dynamics_${SLURM_JOB_ID}"
mkdir "$SPECTRAL_OUTPUT"
/usr/bin/time -v -o "$SPECTRAL_OUTPUT/oracle.resources.txt" \
 "$VALIDATION_ROOT/quspin-venv/bin/python" "$SPECTRAL_TRANCHE/oracle_quspin.py" \
 --manifest "$SPECTRAL_TRANCHE/manifest.toml" --output "$SPECTRAL_OUTPUT/oracle.json"
/usr/bin/time -v -o "$SPECTRAL_OUTPUT/target.resources.txt" \
 "$VALIDATION_ROOT/runtime/julia-1.12.1/bin/julia" --project="$RENYI_ROOT" \
 "$SPECTRAL_TRANCHE/target_ed.jl" --manifest "$SPECTRAL_TRANCHE/manifest.toml" --output "$SPECTRAL_OUTPUT/target.toml"
"$VALIDATION_ROOT/quspin-venv/bin/python" "$SPECTRAL_TRANCHE/compare.py" \
 --manifest "$SPECTRAL_TRANCHE/manifest.toml" --oracle "$SPECTRAL_OUTPUT/oracle.json" \
 --target "$SPECTRAL_OUTPUT/target.toml" --output "$SPECTRAL_OUTPUT/comparison.json"
