#!/usr/bin/env bash
# Submission and serial dependency are owned by root. No nested submissions.
#SBATCH --job-name=thermal_bridges
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
: "${THERMAL_BRIDGES_TRANCHE:?Set absolute deployed tranche path}"
module load anaconda3/2024.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 JULIA_NUM_THREADS=1
export JULIA_DEPOT_PATH="$VALIDATION_ROOT/depot"
export RENYI_ROOT="$VALIDATION_ROOT/source/renyi_disjoint_corrected"
BRIDGES_OUTPUT="$VALIDATION_ROOT/results/thermal_bridges_${SLURM_JOB_ID}"
mkdir "$BRIDGES_OUTPUT"
/usr/bin/time -v -o "$BRIDGES_OUTPUT/oracle.resources.txt" \
 "$VALIDATION_ROOT/quspin-venv/bin/python" "$THERMAL_BRIDGES_TRANCHE/oracle_quspin.py" \
 --manifest "$THERMAL_BRIDGES_TRANCHE/manifest.toml" --output "$BRIDGES_OUTPUT/oracle.json"
/usr/bin/time -v -o "$BRIDGES_OUTPUT/target.resources.txt" \
 "$VALIDATION_ROOT/runtime/julia-1.12.1/bin/julia" --project="$RENYI_ROOT" \
 "$THERMAL_BRIDGES_TRANCHE/target_thermal.jl" --manifest "$THERMAL_BRIDGES_TRANCHE/manifest.toml" --output "$BRIDGES_OUTPUT/target.toml"
"$VALIDATION_ROOT/quspin-venv/bin/python" "$THERMAL_BRIDGES_TRANCHE/compare.py" \
 --manifest "$THERMAL_BRIDGES_TRANCHE/manifest.toml" --oracle "$BRIDGES_OUTPUT/oracle.json" \
 --target "$BRIDGES_OUTPUT/target.toml" --output "$BRIDGES_OUTPUT/comparison.json"
