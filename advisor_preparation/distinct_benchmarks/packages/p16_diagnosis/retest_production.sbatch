#!/usr/bin/env bash
#SBATCH --job-name=cv_disjoint_repair
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
: "${P16_PROOF_JSON:?Set successful independent proof JSON}"
: "${P16_ISOLATED_TOML:?Set successful isolated candidate result TOML}"
export JULIA_DEPOT_PATH="$VALIDATION_ROOT/depot"
export JULIA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export GKSwstype=100
cd "$VALIDATION_ROOT/source"
diag=distinct_benchmarks/packages/p16_diagnosis
attempt="$VALIDATION_ROOT/results/distinct_packages_disjoint_repair_attempt01"
"$VALIDATION_ROOT/quspin-venv/bin/python" "$diag/prepare_retest.py" --root "$VALIDATION_ROOT" --proof "$P16_PROOF_JSON" --isolated "$P16_ISOLATED_TOML" --attempt "$attempt"
export RENYI_ROOT="$VALIDATION_ROOT/source/renyi_disjoint_corrected"
export DISTINCT_TARGET_OUTPUT="$attempt/target.toml"
"$VALIDATION_ROOT/runtime/julia-1.12.1/bin/julia" --startup-file=no --project="$RENYI_ROOT" distinct_benchmarks/packages/user_ed_distinct.jl
"$VALIDATION_ROOT/quspin-venv/bin/python" distinct_benchmarks/packages/compare.py --oracle-dir "$VALIDATION_ROOT/results/distinct_tranche01_attempt01/oracle" --target "$attempt/target.toml" --output "$attempt/packages.json"
export P16_SUBSETS_OUTPUT="$attempt/production_subsets.toml"
"$VALIDATION_ROOT/runtime/julia-1.12.1/bin/julia" --startup-file=no --project="$RENYI_ROOT" "$diag/recheck_production_subsets.jl"
"$VALIDATION_ROOT/quspin-venv/bin/python" "$diag/prepare_retest.py" --audit --root "$VALIDATION_ROOT" --proof "$P16_PROOF_JSON" --attempt "$attempt"
