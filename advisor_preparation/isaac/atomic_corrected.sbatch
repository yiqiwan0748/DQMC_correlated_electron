#!/usr/bin/env bash
#SBATCH --job-name=cv_atomic_corrected
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
"$VALIDATION_ROOT/runtime/julia-1.12.1/bin/julia" --project=renyi_corrected renyi_corrected/test/runtests.jl
export ATOMIC_DIAGNOSTIC_OUTPUT="$VALIDATION_ROOT/results/atomic_corrected_kernel_diagnostic.toml"
"$VALIDATION_ROOT/runtime/julia-1.12.1/bin/julia" --project=renyi_corrected isaac/diagnostics/diagnose_atomic_odd_sweep.jl --kernel
export VALIDATION_ATTEMPT=external_atomic_corrected_4096
export VALIDATION_ID_PREFIX=ATOM-
bash isaac/external_corrected.sbatch
