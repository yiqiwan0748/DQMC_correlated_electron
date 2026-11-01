#!/usr/bin/env bash
#SBATCH --job-name=cv_X06_final_audit
#SBATCH --account=acf-utk0011
#SBATCH --partition=short
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3500M
#SBATCH --time=00:03:00
set -euo pipefail
: "${VALIDATION_ROOT:?Set VALIDATION_ROOT}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$VALIDATION_ROOT/quspin-venv/bin/python" \
  "$VALIDATION_ROOT/source/distinct_benchmarks/dqmc_coverage/summarize_coverage.py" \
  --mode final \
  --inputs "$VALIDATION_ROOT/results/distinct_X06_final01/shard_1_400.toml" \
  --output "$VALIDATION_ROOT/results/distinct_X06_final01/summary_attempt01.json"
