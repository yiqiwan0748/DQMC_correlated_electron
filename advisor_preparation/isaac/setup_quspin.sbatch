#!/usr/bin/env bash
#SBATCH --job-name=cv_quspin_setup
#SBATCH --account=acf-utk0011
#SBATCH --partition=short
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:25:00
set -euo pipefail
: "${VALIDATION_ROOT:?Set VALIDATION_ROOT}"
module load anaconda3/2024.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python -m venv "$VALIDATION_ROOT/quspin-venv"
"$VALIDATION_ROOT/quspin-venv/bin/python" -m pip install --disable-pip-version-check quspin==1.0.1
"$VALIDATION_ROOT/quspin-venv/bin/python" -m pip freeze > "$VALIDATION_ROOT/logs/quspin_requirements.txt"
"$VALIDATION_ROOT/quspin-venv/bin/python" -c 'import quspin,numpy,scipy; print("QuSpin",quspin.__version__,"NumPy",numpy.__version__,"SciPy",scipy.__version__)'
date -u > "$VALIDATION_ROOT/logs/QUSPIN_SETUP_PASS"
