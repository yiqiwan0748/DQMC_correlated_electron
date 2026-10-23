#!/usr/bin/env bash
#SBATCH --job-name=cv_distinct_tranche01
#SBATCH --account=acf-utk0011
#SBATCH --partition=short
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3500M
#SBATCH --time=00:20:00
set -euo pipefail
: "${VALIDATION_ROOT:?Set VALIDATION_ROOT}"
export JULIA_DEPOT_PATH="$VALIDATION_ROOT/depot"
export JULIA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export GKSwstype=100
export RENYI_ROOT="$VALIDATION_ROOT/source/renyi_corrected"
export DISTINCT_ALM_ROOT="$VALIDATION_ROOT/source/distinct_benchmarks/frozen_alm"
cd "$VALIDATION_ROOT/source"
attempt="$VALIDATION_ROOT/results/distinct_tranche01_attempt01"
test ! -e "$attempt"
mkdir -p "$attempt"
export DISTINCT_JOURNAL_OUTPUT="$attempt/journal.toml"
# Continue independent families after a scientific FAIL; preserve every status.
set +e
export ANALYTIC_OUTPUT="$attempt/analytic.toml"
"$VALIDATION_ROOT/runtime/julia-1.12.1/bin/julia" --startup-file=no --project=renyi_corrected distinct_benchmarks/analytic/run_analytic_tranche.jl
analytic_exit=$?
"$VALIDATION_ROOT/runtime/julia-1.12.1/bin/julia" --startup-file=no --project=renyi_corrected distinct_benchmarks/journal/run_raghu.jl
journal_exit=$?
"$VALIDATION_ROOT/quspin-venv/bin/python" distinct_benchmarks/packages/quspin_distinct.py --output-dir "$attempt/oracle"
oracle_exit=$?
export DISTINCT_TARGET_OUTPUT="$attempt/target.toml"
"$VALIDATION_ROOT/runtime/julia-1.12.1/bin/julia" --startup-file=no --project=renyi_corrected distinct_benchmarks/packages/user_ed_distinct.jl
target_exit=$?
"$VALIDATION_ROOT/quspin-venv/bin/python" distinct_benchmarks/packages/compare.py --oracle-dir "$attempt/oracle" --target "$attempt/target.toml" --output "$attempt/packages.json"
compare_exit=$?
printf 'analytic_exit=%s\njournal_exit=%s\noracle_exit=%s\ntarget_exit=%s\ncompare_exit=%s\n' "$analytic_exit" "$journal_exit" "$oracle_exit" "$target_exit" "$compare_exit" > "$attempt/process_status.txt"
exit 0
