#!/usr/bin/env bash
#SBATCH --job-name=cv_journal_v2
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
export JULIA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export JULIA_DEPOT_PATH="$VALIDATION_ROOT/depot"
export DISTINCT_ALM_ROOT="$VALIDATION_ROOT/source/distinct_benchmarks/frozen_alm"
attempt="$VALIDATION_ROOT/results/distinct_journal_attempt02"
test ! -e "$attempt"
mkdir "$attempt"
export DISTINCT_JOURNAL_OUTPUT="$attempt/journal.toml"
"$VALIDATION_ROOT/runtime/julia-1.12.1/bin/julia" --startup-file=no "$VALIDATION_ROOT/source/distinct_benchmarks/journal/run_raghu_v2.jl"
