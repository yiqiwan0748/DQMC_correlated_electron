#!/usr/bin/env bash
#SBATCH --job-name=cv_external_ED
#SBATCH --account=acf-utk0011
#SBATCH --partition=short
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
set -euo pipefail
: "${VALIDATION_ROOT:?Set VALIDATION_ROOT}"
export JULIA_DEPOT_PATH="$VALIDATION_ROOT/stdlib-depot"
export JULIA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export RENYI_ROOT="$VALIDATION_ROOT/source/renyi"
export VALIDATION_ORACLES="$VALIDATION_ROOT/results/quspin_attempt02/oracle_values.toml"
export VALIDATION_OUTPUT="$VALIDATION_ROOT/results/ed_reference_attempt01.toml"
cd "$VALIDATION_ROOT/source"
"$VALIDATION_ROOT/runtime/julia-1.12.1/bin/julia" --startup-file=no validation_design/run_ed_reference_checks.jl
