#!/usr/bin/env bash
# Root launches only after inspecting the measured pilot and resource receipt.
#SBATCH --job-name=tenpy_full
#SBATCH --account=acf-utk0011
#SBATCH --partition=short
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3500M
#SBATCH --time=00:20:00
set -euo pipefail
: "${VALIDATION_ROOT:?Set VALIDATION_ROOT}"
: "${TENPY_TRANCHE:?Set absolute deployed tranche directory}"
module load anaconda3/2024.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 JULIA_NUM_THREADS=1
export JULIA_DEPOT_PATH="$VALIDATION_ROOT/depot"
export RENYI_ROOT="$VALIDATION_ROOT/source/renyi_corrected"
TENPY_PYTHON="$VALIDATION_ROOT/tenpy-venv/bin/python"
TENPY_OUTPUT="$VALIDATION_ROOT/results/tenpy_full_${SLURM_JOB_ID}"
mkdir "$TENPY_OUTPUT"
/usr/bin/time -v -o "$TENPY_OUTPUT/oracle.resources.txt" timeout 480 \
  "$TENPY_PYTHON" "$TENPY_TRANCHE/oracle_tenpy.py" --manifest "$TENPY_TRANCHE/manifest.toml" --output "$TENPY_OUTPUT/oracle.json"
/usr/bin/time -v -o "$TENPY_OUTPUT/target.resources.txt" timeout 600 \
  "$VALIDATION_ROOT/runtime/julia-1.12.1/bin/julia" --project="$RENYI_ROOT" \
  "$TENPY_TRANCHE/target_itensor.jl" --manifest "$TENPY_TRANCHE/manifest.toml" --output "$TENPY_OUTPUT/target.toml"
"$TENPY_PYTHON" "$TENPY_TRANCHE/compare.py" --manifest "$TENPY_TRANCHE/manifest.toml" \
  --oracle "$TENPY_OUTPUT/oracle.json" --target "$TENPY_OUTPUT/target.toml" --output "$TENPY_OUTPUT/comparison.json"
