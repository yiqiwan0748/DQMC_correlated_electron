#!/usr/bin/env bash
#SBATCH --job-name=alf_adapter_v3
#SBATCH --account=acf-utk0011
#SBATCH --partition=short
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3500M
#SBATCH --time=00:22:00
#SBATCH --dependency=afterany:6458547
# Adapter eligibility only; root owns announcement/submission and resource check.
# Estimated2-6minutes after start; queue delay is separate. Zero P21 completions.
# BEGIN EXACT REVIEWED V3 PLAN BODY
set -euo pipefail
export VALIDATION_ROOT=/lustre/isaac24/scratch/ywang287/advisor-validation-20261022
R="$VALIDATION_ROOT"
archive="$R/source_releases/alf-measurement-observers-v3-source.tar.gz"
echo 'eaff11972ff34831893d51132ad805b6e415703a3dc10edc97a910dea050bf06  '"$archive" | sha256sum -c -
mkdir "$R/source/p21_measurement_proof_v3"
tar -xzf "$archive" -C "$R/source/p21_measurement_proof_v3"
observer_root="$R/source/p21_measurement_proof_v3/alf_preparation/measurement_preparation_v3"
export ALF_NATIVE_ROOT="$R/external/alf_build_attempt01/alf_build/ALF"
export P21_LINKER_COMPAT="$R/external/alf_build_attempt01/compat"
export P21_USER_BASE="$R/source/renyi_disjoint_corrected"
export P21_PROOF_ROOT="$R/results/p21_adapter_proof_v3_attempt01"
export P21_PYTHON="$R/quspin-venv/bin/python"
export P21_JULIA="$R/runtime/julia-1.12.1/bin/julia"
export JULIA_DEPOT_PATH="$R/depot"
test "$(cat "$R/external/alf_build_attempt01/status.txt")" = ALF_BUILT_NOT_BENCHMARKED
test "$(cat "$R/external/alf_build_attempt01/alf_build/status.txt")" = SETUP_BUILT_NOT_BENCHMARKED
test -x "$ALF_NATIVE_ROOT/Prog/ALF.out"
(cd "$ALF_NATIVE_ROOT" && sha256sum -c ../binary_sha256.txt)
python3 - "$observer_root" "$ALF_NATIVE_ROOT" <<'PY'
import hashlib,json,sys
from pathlib import Path
b,n=map(Path,sys.argv[1:])
for rel,want in json.loads((b/'pre_execution_receipt.json').read_text())['files'].items():
    assert hashlib.sha256((b/rel).read_bytes()).hexdigest()==want,rel
for row in json.loads((b/'upstream_files.json').read_text())['files']:
    assert hashlib.sha256((n/row['path']).read_bytes()).hexdigest()==row['sha256'],row['path']
PY
if ! type module >/dev/null 2>&1; then source /etc/profile.d/modules.sh; fi
module load gcc/10.2.0 openblas/0.3.12
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 JULIA_NUM_THREADS=1
provider=/spack/spack-0.16.3/apps/linux-rhel8-cascadelake/gcc-10.2.0/openblas-0.3.12-c4x3eon3lkfddumq466n3n6b5gpqdbvh/lib/libopenblas.so
python3 "$observer_root/check_native_linkage.py" "$ALF_NATIVE_ROOT/Prog/ALF.out" "$provider" "$R/results/p21_native_linkage_v3_attempt01.json"
bash "$observer_root/run_adapter_proofs.sh"
