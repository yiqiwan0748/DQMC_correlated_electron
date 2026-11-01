#!/usr/bin/env bash
#SBATCH --job-name=cv_validate_setup
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
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1
export JULIA_PKG_PRECOMPILE_AUTO=0
julia_bin="$VALIDATION_ROOT/runtime/julia-1.12.1/bin/julia"
sha256sum -c SHA256SUMS > "$VALIDATION_ROOT/logs/source_verification.txt"
"$julia_bin" --project=renyi -e 'using Pkg; Pkg.instantiate(); Pkg.precompile(); using HubbardRenyiDQMC, ITensors, ITensorMPS; println("RENYI_ENVIRONMENT_READY")'
"$julia_bin" --project=alm -e 'using Pkg; Pkg.instantiate(); Pkg.precompile(); using AltermagnetDQMCTrial; println("ALM_ENVIRONMENT_READY")'
date -u > "$VALIDATION_ROOT/logs/SETUP_PASS"
