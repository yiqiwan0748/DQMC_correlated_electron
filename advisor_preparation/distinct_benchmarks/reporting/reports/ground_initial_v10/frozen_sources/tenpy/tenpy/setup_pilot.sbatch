#!/usr/bin/env bash
# Root submits with --dependency=afterany:6457686. This file never submits jobs.
#SBATCH --job-name=tenpy_pilot
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
: "${TENPY_TRANCHE:?Set absolute deployed tranche directory}"
module load anaconda3/2024.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 JULIA_NUM_THREADS=1
export JULIA_DEPOT_PATH="$VALIDATION_ROOT/depot"
export RENYI_ROOT="$VALIDATION_ROOT/source/renyi_corrected"
TENPY_VENV="$VALIDATION_ROOT/tenpy-venv"
TENPY_OUTPUT="$VALIDATION_ROOT/results/tenpy_pilot_${SLURM_JOB_ID}"
mkdir "$TENPY_OUTPUT"
if [[ ! -d "$TENPY_VENV" ]]; then python -m venv "$TENPY_VENV"; fi
"$TENPY_VENV/bin/python" -m pip install --only-binary=:all: --disable-pip-version-check -r "$TENPY_TRANCHE/requirements.txt"
"$TENPY_VENV/bin/python" -m pip freeze > "$TENPY_OUTPUT/pip_freeze.txt"
/usr/bin/time -v -o "$TENPY_OUTPUT/oracle.resources.txt" timeout 240 \
  "$TENPY_VENV/bin/python" "$TENPY_TRANCHE/oracle_tenpy.py" --manifest "$TENPY_TRANCHE/manifest.toml" --output "$TENPY_OUTPUT/oracle.json" --pilot
/usr/bin/time -v -o "$TENPY_OUTPUT/target.resources.txt" timeout 240 \
  "$VALIDATION_ROOT/runtime/julia-1.12.1/bin/julia" --project="$RENYI_ROOT" \
  "$TENPY_TRANCHE/target_itensor.jl" --manifest "$TENPY_TRANCHE/manifest.toml" --output "$TENPY_OUTPUT/target.toml" --pilot
"$TENPY_VENV/bin/python" "$TENPY_TRANCHE/compare.py" --manifest "$TENPY_TRANCHE/manifest.toml" \
  --oracle "$TENPY_OUTPUT/oracle.json" --target "$TENPY_OUTPUT/target.toml" --output "$TENPY_OUTPUT/comparison.json"
