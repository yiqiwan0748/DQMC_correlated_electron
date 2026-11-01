#!/usr/bin/env bash
#SBATCH --job-name=cv_X06_final
#SBATCH --account=acf-utk0011
#SBATCH --partition=short
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3500M
#SBATCH --time=02:00:00
set -euo pipefail
: "${VALIDATION_ROOT:?Set VALIDATION_ROOT}"
: "${COVERAGE_FIRST:?Set first frozen replicate}"
: "${COVERAGE_LAST:?Set last frozen replicate}"
: "${COVERAGE_PILOT_SUMMARY:?Set actual pilot summary JSON}"
"$VALIDATION_ROOT/quspin-venv/bin/python" -c 'import json,sys; x=json.load(open(sys.argv[1])); assert x["status"]=="PILOT_ONLY" and x["raw_audit"]=="PASS"' "$COVERAGE_PILOT_SUMMARY"
export JULIA_DEPOT_PATH="$VALIDATION_ROOT/depot"
export JULIA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export RENYI_ROOT="$VALIDATION_ROOT/source/renyi_corrected"
export COVERAGE_ORACLE="$VALIDATION_ROOT/results/quspin_attempt02/oracle_values.toml"
export COVERAGE_MODE=final
cd "$VALIDATION_ROOT/source"
mkdir -p "$VALIDATION_ROOT/results/distinct_X06_final01"
export COVERAGE_OUTPUT="$VALIDATION_ROOT/results/distinct_X06_final01/shard_${COVERAGE_FIRST}_${COVERAGE_LAST}.toml"
"$VALIDATION_ROOT/runtime/julia-1.12.1/bin/julia" --startup-file=no --project="$RENYI_ROOT" distinct_benchmarks/dqmc_coverage/run_coverage.jl
