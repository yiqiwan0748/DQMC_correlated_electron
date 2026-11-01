#!/usr/bin/env bash
#SBATCH --job-name=cv_X06_pilot
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
export JULIA_DEPOT_PATH="$VALIDATION_ROOT/depot"
export JULIA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export RENYI_ROOT="$VALIDATION_ROOT/source/renyi_corrected"
export COVERAGE_ORACLE="$VALIDATION_ROOT/results/quspin_attempt02/oracle_values.toml"
export COVERAGE_MODE=pilot COVERAGE_FIRST=1 COVERAGE_LAST=4
cd "$VALIDATION_ROOT/source"
attempt="$VALIDATION_ROOT/results/distinct_X06_pilot01"
test ! -e "$attempt"
mkdir -p "$attempt"
export COVERAGE_OUTPUT="$attempt/pilot.toml"
"$VALIDATION_ROOT/runtime/julia-1.12.1/bin/julia" --startup-file=no --project="$RENYI_ROOT" distinct_benchmarks/dqmc_coverage/run_coverage.jl
"$VALIDATION_ROOT/quspin-venv/bin/python" distinct_benchmarks/dqmc_coverage/summarize_coverage.py --mode pilot --inputs "$attempt/pilot.toml" --output "$attempt/summary.json"
