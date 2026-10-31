#!/usr/bin/env bash
#SBATCH --job-name=cv_validate_baseline
#SBATCH --account=acf-utk0011
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:45:00
set -euo pipefail
: "${VALIDATION_ROOT:?Set VALIDATION_ROOT at submission}"
cd "$VALIDATION_ROOT/source"
export JULIA_DEPOT_PATH="$VALIDATION_ROOT/depot"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export GKSwstype=100
julia_bin="$VALIDATION_ROOT/runtime/julia-1.12.1/bin/julia"
"$julia_bin" --project=renyi renyi/test/runtests.jl > "$VALIDATION_ROOT/logs/renyi_baseline.log" 2>&1
"$julia_bin" --project=renyi renyi/analysis/validate_grover_rank1_update.jl > "$VALIDATION_ROOT/logs/grover_baseline.log" 2>&1
"$julia_bin" --project=alm isaac/alm_baseline.jl > "$VALIDATION_ROOT/logs/alm_baseline.log" 2>&1
date -u > "$VALIDATION_ROOT/logs/BASELINE_PASS"
