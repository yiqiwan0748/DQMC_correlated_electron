#!/usr/bin/env bash
#SBATCH --job-name=pyscf_setup_v1
#SBATCH --account=acf-utk0011
#SBATCH --partition=short
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3500M
#SBATCH --time=00:10:00
#SBATCH --no-requeue
set -euo pipefail
: "${VALIDATION_ROOT:?}"
: "${PYSCF_SETUP_ROOT:?}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
module load anaconda3/2024.06
out="$VALIDATION_ROOT/results/pyscf_setup_v1_${SLURM_JOB_ID}"
mkdir "$out"
(cd "$PYSCF_SETUP_ROOT/wheels" && sha256sum -c SHA256SUMS)
test ! -e "$VALIDATION_ROOT/pyscf-venv"
/usr/bin/time -v -o "$out/install.resources.txt" bash -c '
set -euo pipefail
"$VALIDATION_ROOT/quspin-venv/bin/python" -m venv "$VALIDATION_ROOT/pyscf-venv"
"$VALIDATION_ROOT/pyscf-venv/bin/python" -m pip install --no-index --find-links="$PYSCF_SETUP_ROOT/wheels" -r "$PYSCF_SETUP_ROOT/requirements.txt"
'
"$VALIDATION_ROOT/pyscf-venv/bin/python" -m pip list --format=json > "$out/installed_packages.json"
"$VALIDATION_ROOT/pyscf-venv/bin/python" - "$out/environment.json" <<'PY'
import sys,json,platform
import numpy,scipy,h5py,pyscf
from pathlib import Path
v={x.__name__:x.__version__ for x in (numpy,scipy,h5py,pyscf)}
assert v=={'numpy':'2.4.6','scipy':'1.17.1','h5py':'3.15.1','pyscf':'2.14.0'},v
Path(sys.argv[1]).write_text(json.dumps({'status':'ENVIRONMENT_READY_NOT_SCIENTIFIC_PASS','versions':v,'python':sys.version,'platform':platform.platform(),'scientific_solver_executed':False},indent=2)+'\n')
print(json.dumps(v))
PY
