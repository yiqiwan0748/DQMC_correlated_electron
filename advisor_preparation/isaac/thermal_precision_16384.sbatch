#!/usr/bin/env bash
#SBATCH --job-name=cv_thermal_precision
#SBATCH --account=acf-utk0011
#SBATCH --partition=short
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
set -euo pipefail
: "${VALIDATION_ROOT:?Set VALIDATION_ROOT}"
export VALIDATION_ATTEMPT=external_thermal_corrected_16384
export VALIDATION_ID_PREFIX=QT-
export VALIDATION_ID_LIST='QT-09,QT-15,QT-16'
export VALIDATION_MC_MEASUREMENTS=16384
cd "$VALIDATION_ROOT/source"
bash isaac/external_corrected.sbatch
