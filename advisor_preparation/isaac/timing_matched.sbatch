#!/usr/bin/env bash
#SBATCH --job-name=cv_matched_timing
#SBATCH --account=acf-utk0011
#SBATCH --partition=short
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
set -euo pipefail
: "${VALIDATION_ROOT:?Set VALIDATION_ROOT}"
export JULIA_DEPOT_PATH="$VALIDATION_ROOT/depot"
export JULIA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export GKSwstype=100
export RENYI_ROOT="$VALIDATION_ROOT/source/renyi_corrected"
cd "$VALIDATION_ROOT/source"
export TIMING_SIZES=8 TIMING_REPEATS=3 TIMING_ROUNDS=2
export TIMING_OUTPUT="$VALIDATION_ROOT/results/matched_timing_pilot_20261022.toml"
"$VALIDATION_ROOT/runtime/julia-1.12.1/bin/julia" --project=renyi_corrected isaac/diagnostics/benchmark_matched_lazy_grover.jl
export TIMING_SIZES=8,16,32,64,128 TIMING_REPEATS=7 TIMING_ROUNDS=12
export TIMING_OUTPUT="$VALIDATION_ROOT/results/matched_timing_full_20261022.toml"
"$VALIDATION_ROOT/runtime/julia-1.12.1/bin/julia" --project=renyi_corrected isaac/diagnostics/benchmark_matched_lazy_grover.jl
