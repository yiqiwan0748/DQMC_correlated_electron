#!/usr/bin/env bash
# New wrapper around the preserved v1 candidate. Default is only the ABI probe.
set -euo pipefail
: "${ALF_LINKER_ROOT:?Choose a NEW dedicated receipt directory}"
case "$ALF_LINKER_ROOT" in *[[:space:]]*) exit 2;; esac
case "${1:-probe}" in probe|build) mode=${1:-probe};; *) exit 2;; esac
test ! -e "$ALF_LINKER_ROOT"
mkdir -p "$ALF_LINKER_ROOT"
own=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
if ! type module >/dev/null 2>&1; then
  if [[ -r /etc/profile.d/modules.sh ]]; then source /etc/profile.d/modules.sh; fi
fi
type module >/dev/null 2>&1
module load gcc/10.2.0 openblas/0.3.12
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
openblas_dir=/spack/spack-0.16.3/apps/linux-rhel8-cascadelake/gcc-10.2.0/openblas-0.3.12-c4x3eon3lkfddumq466n3n6b5gpqdbvh/lib
provider=$(readlink -f "$openblas_dir/libopenblas.so")
test -f "$provider"
{
  date -u
  uname -a
  if command -v lscpu >/dev/null 2>&1; then lscpu; fi
  module -t list 2>&1
  module show gcc/10.2.0 2>&1
  module show openblas/0.3.12 2>&1
  command -v gfortran
  gfortran --version
  gfortran -print-search-dirs
  for name in libblas.so liblapack.so libopenblas.so; do gfortran -print-file-name="$name"; done
  ls -l "$openblas_dir"
  if command -v ldconfig >/dev/null 2>&1; then ldconfig -p | awk '/libblas|liblapack|libopenblas|libgfortran/'; fi
  printf 'PROVIDER=%s\n' "$provider"
} > "$ALF_LINKER_ROOT/linker_inventory.txt"
test "$(gfortran -dumpfullversion)" = "10.2.0"
python3 "$own/check_openblas.py" "$provider" "$ALF_LINKER_ROOT/openblas_metadata.json" > "$ALF_LINKER_ROOT/openblas_metadata.log" 2>&1
# GNU ALF configure checks -llapack -lblas. Use isolated aliases to this SAME
# verified provider if Spack offers only -lopenblas; no system files are changed.
mkdir "$ALF_LINKER_ROOT/compat"
ln -s "$provider" "$ALF_LINKER_ROOT/compat/libblas.so"
ln -s "$provider" "$ALF_LINKER_ROOT/compat/liblapack.so"
export LIBRARY_PATH="$ALF_LINKER_ROOT/compat${LIBRARY_PATH:+:$LIBRARY_PATH}"
export LD_LIBRARY_PATH="$openblas_dir${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
printf 'LIBRARY_PATH=%s\nLD_LIBRARY_PATH=%s\n' "$LIBRARY_PATH" "$LD_LIBRARY_PATH" > "$ALF_LINKER_ROOT/linker_environment.txt"
timeout --signal=TERM --kill-after=5s 60s gfortran -std=f2008 -O0 -g -fcheck=all -ffree-line-length-none "$own/lp64_probe.f90" -o "$ALF_LINKER_ROOT/lp64_probe" -L"$ALF_LINKER_ROOT/compat" -Wl,-rpath,"$openblas_dir" -llapack -lblas -fopenmp > "$ALF_LINKER_ROOT/probe_compile.log" 2>&1
ldd "$ALF_LINKER_ROOT/lp64_probe" > "$ALF_LINKER_ROOT/probe_linkage.txt"
python3 - "$ALF_LINKER_ROOT/probe_linkage.txt" "$provider" <<'PY'
from pathlib import Path
import sys
lines=Path(sys.argv[1]).read_text().splitlines();wanted=Path(sys.argv[2]).resolve();providers=[]
for line in lines:
 if any(name in line for name in ['libopenblas','liblapack','libblas']):
  assert '=>' in line and 'not found' not in line,line
  path=Path(line.split('=>',1)[1].strip().split()[0]).resolve();providers.append(path)
assert providers and all(p==wanted for p in providers),(providers,wanted)
PY
timeout --signal=TERM --kill-after=2s 10s "$ALF_LINKER_ROOT/lp64_probe" > "$ALF_LINKER_ROOT/probe_run.log" 2>&1
awk '/^LP64_BLAS_LAPACK_PROBE_PASS$/ {found=1} END {exit !found}' "$ALF_LINKER_ROOT/probe_run.log"
sha256sum "$own/lp64_probe.f90" "$own/check_openblas.py" "$own/prepare_and_build.sh" "$ALF_LINKER_ROOT/lp64_probe" > "$ALF_LINKER_ROOT/probe_sha256.txt"
printf '%s\n' 'LP64_SETUP_PROBE_PASS_NO_ALF_SCIENCE' > "$ALF_LINKER_ROOT/status.txt"
if [[ "$mode" == build ]]; then
  export ALF_PREPARATION_ROOT="$ALF_LINKER_ROOT/alf_build"
  python3 - "$own/../setup_candidate.sh" <<'PY_CHECK'
from pathlib import Path
import hashlib,sys
assert hashlib.sha256(Path(sys.argv[1]).read_bytes()).hexdigest()=="915bad50ae59d6cbe26ee072e2e92d0c1ec07e087f41e5eff35b7585f24282f3",'Preserved v1 setup source differs'
PY_CHECK
  bash "$own/../setup_candidate.sh" build
  ldd "$ALF_PREPARATION_ROOT/ALF/Prog/ALF.out" > "$ALF_LINKER_ROOT/alf_linkage.txt"
  printf '%s\n' 'ALF_BUILT_NOT_BENCHMARKED' > "$ALF_LINKER_ROOT/status.txt"
fi
