#!/usr/bin/env bash
#SBATCH --job-name=cv_p16_sign_proof
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
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1
export JULIA_DEPOT_PATH="$VALIDATION_ROOT/depot"
export RENYI_ROOT="$VALIDATION_ROOT/source/renyi_corrected"
driver="$VALIDATION_ROOT/source/distinct_benchmarks/packages/p16_diagnosis"
attempt="$VALIDATION_ROOT/results/distinct_p16_diagnosis_attempt01"
test ! -e "$attempt"
mkdir "$attempt"
"$VALIDATION_ROOT/quspin-venv/bin/python" "$driver/verify_quspin_regrouping.py" --output "$attempt/permutation_proof.json"
export P16_ORACLE_TOML="$VALIDATION_ROOT/results/distinct_tranche01_attempt01/oracle/oracle.toml"
export P16_RECHECK_OUTPUT="$attempt/candidate_recheck.toml"
"$VALIDATION_ROOT/runtime/julia-1.12.1/bin/julia" --startup-file=no "$driver/recheck_user_v2.jl"
