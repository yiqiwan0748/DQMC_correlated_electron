#!/usr/bin/env bash
#SBATCH --job-name=cv_ALF_LP64_setup
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
: "${ALF_LINKER_ROOT:?Choose a new dedicated linker receipt directory}"
bash "$VALIDATION_ROOT/source/distinct_benchmarks/alf_preparation/setup_v2/prepare_and_build.sh" "${ALF_SETUP_V2_MODE:-probe}"
