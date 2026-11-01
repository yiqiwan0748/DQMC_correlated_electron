#!/usr/bin/env bash
#SBATCH --job-name=cv_external_validate
#SBATCH --account=acf-utk0011
#SBATCH --partition=short
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:45:00
set -euo pipefail
: "${VALIDATION_ROOT:?Set VALIDATION_ROOT}"
: "${VALIDATION_ATTEMPT:?Set unique VALIDATION_ATTEMPT}"
export JULIA_DEPOT_PATH="$VALIDATION_ROOT/depot"
export JULIA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export GKSwstype=100
export RENYI_ROOT="$VALIDATION_ROOT/source/renyi"
export VALIDATION_ORACLES="$VALIDATION_ROOT/results/quspin_attempt02/oracle_values.toml"
export VALIDATION_OUTPUT="$VALIDATION_ROOT/results/$VALIDATION_ATTEMPT.toml"
export VALIDATION_RUN_DQMC="${VALIDATION_RUN_DQMC:-0}"
export VALIDATION_MC_MEASUREMENTS="${VALIDATION_MC_MEASUREMENTS:-4096}"
cd "$VALIDATION_ROOT/source"
"$VALIDATION_ROOT/runtime/julia-1.12.1/bin/julia" --project=renyi validation_design/run_external_validation.jl
