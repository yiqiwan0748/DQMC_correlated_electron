#!/usr/bin/env bash
#SBATCH --job-name=cv_ALF_setup_only
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
: "${ALF_PREPARATION_ROOT:?Choose new dedicated ALF setup output directory}"
# The operator loads a coherent GNU/LP64 BLAS+LAPACK module set before submission.
# Inventory mode is the default. Build mode must be named explicitly after inventory review.
bash "$VALIDATION_ROOT/source/distinct_benchmarks/alf_preparation/setup_candidate.sh" "${ALF_SETUP_MODE:-inventory}"
