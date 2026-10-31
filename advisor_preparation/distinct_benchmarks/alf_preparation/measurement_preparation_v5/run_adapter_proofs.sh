#!/usr/bin/env bash
# Candidate execution command. No submission or external statistical benchmark.
# Native fixed-configuration propagation plus32 diagnostic user update calls.
set -euo pipefail
: "${ALF_NATIVE_ROOT:?Existing pinned native ALF build tree}"
: "${P21_USER_BASE:?Existing corrected user source tree}"
: "${P21_PROOF_ROOT:?Choose a NEW result/build directory}"
: "${P21_LINKER_COMPAT:?Existing isolated LP64 libblas/liblapack compatibility directory}"
: "${P21_PYTHON:?Pinned QuSpin-venv Python with NumPy/SciPy}"
: "${P21_JULIA:?Existing Julia1.12.1 executable}"
own=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
test ! -e "$P21_PROOF_ROOT"
case "$P21_PROOF_ROOT" in *[[:space:]]*) exit 2;; esac
test "$(git -C "$ALF_NATIVE_ROOT" rev-parse HEAD)" = f4c500cc5d7bcf089400c63e4e75382c558a4ed3
if ! type module >/dev/null 2>&1; then source /etc/profile.d/modules.sh; fi
module load gcc/10.2.0 openblas/0.3.12
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 JULIA_NUM_THREADS=1
provider_dir=/spack/spack-0.16.3/apps/linux-rhel8-cascadelake/gcc-10.2.0/openblas-0.3.12-c4x3eon3lkfddumq466n3n6b5gpqdbvh/lib
provider=$(readlink -f "$provider_dir/libopenblas.so")
test "$(readlink -f "$P21_LINKER_COMPAT/libblas.so")" = "$provider"
test "$(readlink -f "$P21_LINKER_COMPAT/liblapack.so")" = "$provider"
export LIBRARY_PATH="$P21_LINKER_COMPAT${LIBRARY_PATH:+:$LIBRARY_PATH}"
export LD_LIBRARY_PATH="$provider_dir${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
mkdir -p "$P21_PROOF_ROOT"
python3 "$own/prepare_observers.py" --alf-root "$ALF_NATIVE_ROOT" --user-root "$P21_USER_BASE" --output "$P21_PROOF_ROOT/overlay" > "$P21_PROOF_ROOT/overlay.log"
# The recipe must produce the independently reviewed exact overlay hashes.
python3 - "$own/overlay_v1/overlay_receipt.json" "$P21_PROOF_ROOT/overlay/overlay_receipt.json" <<'PY'
import json,sys
a,b=[json.load(open(p)) for p in sys.argv[1:]]
assert a['alf_commit']==b['alf_commit'] and a['generator_sha256']==b['generator_sha256']
assert a['files']==b['files']
PY
cp -a "$ALF_NATIVE_ROOT" "$P21_PROOF_ROOT/ALF_observer"
cp "$P21_PROOF_ROOT/overlay/alf/Prog/main.F90" "$P21_PROOF_ROOT/ALF_observer/Prog/main.F90"
cp "$P21_PROOF_ROOT/overlay/alf/Prog/Hamiltonians/Hamiltonian_Hubbard_smod.F90" "$P21_PROOF_ROOT/ALF_observer/Prog/Hamiltonians/Hamiltonian_Hubbard_smod.F90"
cd "$P21_PROOF_ROOT/ALF_observer"
export ALF_FLAGS_EXT='-fno-fast-math -DOBS_LEGACY'
set +u
source configure.sh GNU noMPI STAB3 NO-INTERACTIVE NO-FALLBACK > "$P21_PROOF_ROOT/configure.log" 2>&1
set -u
[[ "$ALF_FLAGS_PROG" != *-DMPI* && "$ALF_FLAGS_PROG" != *-DHDF5* ]]
[[ "$ALF_FLAGS_PROG" == *-DOBS_LEGACY* ]]
printf 'ALF_FC=%s\nALF_LIB=%s\nALF_FLAGS_PROG=%s\n' "$ALF_FC" "$ALF_LIB" "$ALF_FLAGS_PROG" > "$P21_PROOF_ROOT/effective_build.txt"
# v2 provenance repair: clear all copied native objects in the NEW tree.
# Upstream dependencies do not track flag changes, so a clean rebuild is required.
timeout --signal=TERM --kill-after=5s 60s make -j1 cleanall > "$P21_PROOF_ROOT/clean.log" 2>&1
# v3: upstream cleanall omits this one known generated Fortran module.
# Remove ONLY that untracked compiler artifact in the NEW copied tree, with SHA log.
python3 - "$P21_PROOF_ROOT/ALF_observer" "$P21_PROOF_ROOT/known_compiler_module_cleanup.json" "$ALF_NATIVE_ROOT" > "$P21_PROOF_ROOT/known_compiler_module_cleanup.log" 2>&1 <<'PY_KNOWN_MODULE'
from pathlib import Path
import hashlib,json,subprocess,sys
root=Path(sys.argv[1]).resolve(strict=True);output=Path(sys.argv[2]);native=Path(sys.argv[3]).resolve(strict=True)
assert root.name=='ALF_observer' and root!=native,'Cleanup must be inside the NEW observer copy'
assert not output.exists(),'Use a new cleanup receipt'
rel='Prog/opt_time_dependent_mod.mod';target=root/rel
assert not target.is_symlink(),'Known compiler artifact must not be a symlink'
tracked=subprocess.run(['git','-C',str(root),'ls-files','--error-unmatch','--',rel],capture_output=True,text=True,timeout=10)
assert tracked.returncode==1 and not tracked.stdout.strip(),'Do not remove a tracked ALF source file'
receipt={'schema':'p21_known_compiler_module_cleanup_v1','relative_path':rel,'scope':'NEW copied ALF_observer only',
         'source_tree':str(root),'preserved_native_tree':str(native),'tracked_source':False,'strict_stale_assertion_unchanged':True}
if target.exists():
 assert target.is_file(),'Known compiler artifact is not a regular file'
 data=target.read_bytes();receipt.update(bytes=len(data),sha256_before_removal=hashlib.sha256(data).hexdigest())
 target.unlink();receipt['status']='KNOWN_UNTRACKED_COMPILER_MODULE_REMOVED'
else:
 receipt['status']='KNOWN_COMPILER_MODULE_ALREADY_ABSENT'
output.write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
PY_KNOWN_MODULE
python3 - "$P21_PROOF_ROOT/ALF_observer" <<'PY_CLEAN'
from pathlib import Path
import sys
root=Path(sys.argv[1]);stale=[]
for sub in ('Prog','Libraries','Analysis'):
 for path in (root/sub).rglob('*'):
  if path.is_file() and path.suffix.lower() in ('.o','.mod','.smod','.a'):
   stale.append(str(path.relative_to(root)))
assert not stale,('Clean target left reusable objects/archives',stale)
assert not (root/'Prog/ALF.out').exists(),'Clean target left copied ALF binary'
print('COPIED_OBJECTS_MODULES_ARCHIVES_REMOVED_PASS')
PY_CLEAN
timeout --signal=TERM --kill-after=15s 15m make -j1 program > "$P21_PROOF_ROOT/build.log" 2>&1
binary="$P21_PROOF_ROOT/ALF_observer/Prog/ALF.out"
python3 "$own/check_native_linkage.py" "$binary" "$provider" "$P21_PROOF_ROOT/observer_linkage_audit.json" > "$P21_PROOF_ROOT/linkage.log"
mkdir "$P21_PROOF_ROOT/fixed_configuration"
cp "$own/parameters_proof" "$P21_PROOF_ROOT/fixed_configuration/parameters"
# v5: native serial ALF reads one positive integer before Hamiltonian setup.
# Fixed proof fields are still overwritten before UDV/Green construction.
test ! -e "$P21_PROOF_ROOT/fixed_configuration/seeds"
cp "$own/seeds_proof" "$P21_PROOF_ROOT/fixed_configuration/seeds"
(cd "$P21_PROOF_ROOT/fixed_configuration" && sha256sum seeds) > "$P21_PROOF_ROOT/native_startup_seed_sha256.txt"
cd "$P21_PROOF_ROOT/fixed_configuration"
P21_PROOF_ONLY=1 timeout --signal=TERM --kill-after=5s 90s "$binary" > "$P21_PROOF_ROOT/native_fixed_configuration.log" 2>&1
"$P21_PYTHON" "$own/check_alf_proof.py" --run-dir "$P21_PROOF_ROOT/fixed_configuration" --alf-root "$P21_PROOF_ROOT/ALF_observer" --overlay-receipt "$P21_PROOF_ROOT/overlay/overlay_receipt.json" --binary "$binary" --output "$P21_PROOF_ROOT/alf_adapter_proof.json" > "$P21_PROOF_ROOT/alf_adapter_proof.log" 2>&1
export P21_USER_OVERLAY="$P21_PROOF_ROOT/overlay/user"
export P21_USER_PROOF_OUTPUT="$P21_PROOF_ROOT/user_adapter_proof.toml"
timeout --signal=TERM --kill-after=10s 180s "$P21_JULIA" --startup-file=no --project="$P21_USER_BASE" "$own/user_buffer_probe.jl" > "$P21_PROOF_ROOT/user_adapter_proof.log" 2>&1
printf '%s\n' 'ADAPTER_PROOFS_PASS_P21_EXTERNAL_BENCHMARK_NOT_RUN' > "$P21_PROOF_ROOT/status.txt"
