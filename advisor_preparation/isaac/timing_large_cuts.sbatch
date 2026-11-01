#!/usr/bin/env bash
#SBATCH --job-name=cv_large_cut_timing
#SBATCH --account=acf-utk0011
#SBATCH --partition=short
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
set -euo pipefail
: "${VALIDATION_ROOT:?Set VALIDATION_ROOT}"
export JULIA_DEPOT_PATH="$VALIDATION_ROOT/depot"
export JULIA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export GKSwstype=100
export RENYI_ROOT="$VALIDATION_ROOT/source/renyi_corrected"
cd "$VALIDATION_ROOT/source"
export TIMING_SIZES=64 TIMING_REPEATS=7 TIMING_ROUNDS=12
export TIMING_NA=48 TIMING_OUTPUT="$VALIDATION_ROOT/results/matched_timing_cut48_20261022.toml"
"$VALIDATION_ROOT/runtime/julia-1.12.1/bin/julia" --project=renyi_corrected isaac/diagnostics/benchmark_matched_large_cuts.jl
export TIMING_NA=64 TIMING_OUTPUT="$VALIDATION_ROOT/results/matched_timing_cut64_20261022.toml"
"$VALIDATION_ROOT/runtime/julia-1.12.1/bin/julia" --project=renyi_corrected isaac/diagnostics/benchmark_matched_large_cuts.jl
export VALIDATION_ATTEMPT=external_atomic_corrected_16384
export VALIDATION_ID_PREFIX=ATOM- VALIDATION_ID_LIST=ATOM-05
export VALIDATION_MC_MEASUREMENTS=16384
bash isaac/external_corrected.sbatch
"$VALIDATION_ROOT/runtime/julia-1.12.1/bin/julia" --project=renyi_corrected isaac/diagnostics/plot_matched_timing.jl "$VALIDATION_ROOT/results/matched_timing_full_20261022.toml" "$VALIDATION_ROOT/results/matched_timing_scaling_20261022"
