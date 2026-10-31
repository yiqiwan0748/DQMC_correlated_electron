#!/usr/bin/env bash
#SBATCH --job-name=cv_quspin_oracle
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
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
module load anaconda3/2024.06
"$VALIDATION_ROOT/quspin-venv/bin/python" "$VALIDATION_ROOT/source/external_oracles/quspin_export.py" \
  --manifest "$VALIDATION_ROOT/source/external_oracles/cases.json" \
  --output-dir "$VALIDATION_ROOT/results/quspin_attempt01"
