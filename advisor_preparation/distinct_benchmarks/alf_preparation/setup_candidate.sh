#!/usr/bin/env bash
# Candidate only. No ALF.out execution, no physics samples, no job submission.
set -euo pipefail
: "${ALF_PREPARATION_ROOT:?Choose a NEW dedicated directory on ISAAC}"
case "$ALF_PREPARATION_ROOT" in *[[:space:]]*) printf '%s\n' 'ALF build paths must not contain whitespace' >&2; exit 2;; esac
case "${1:-inventory}" in inventory|build) mode=${1:-inventory};; *) exit 2;; esac
if [[ -e "$ALF_PREPARATION_ROOT" ]]; then
  printf '%s\n' 'Choose a new ALF_PREPARATION_ROOT; existing receipts are preserved' >&2
  exit 2
fi
mkdir -p "$ALF_PREPARATION_ROOT"
{
  date -u
  uname -a
  if type module >/dev/null 2>&1; then
    module -t list 2>&1
    module -t avail GCC GCCcore OpenBLAS BLAS LAPACK 2>&1
  fi
  for command in git make python3 gfortran timeout; do
    command -v "$command" || true
    "$command" --version 2>&1 | head -n 2 || true
  done
} > "$ALF_PREPARATION_ROOT/environment_inventory.txt"
if [[ "$mode" == inventory ]]; then
  printf '%s\n' 'INVENTORY_ONLY_NO_INSTALL' > "$ALF_PREPARATION_ROOT/status.txt"
  exit 0
fi
for command in git make python3 gfortran timeout; do command -v "$command" >/dev/null; done
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
alf_pin=f4c500cc5d7bcf089400c63e4e75382c558a4ed3
alf_src="$ALF_PREPARATION_ROOT/ALF"
timeout --signal=TERM --kill-after=10s 180s git clone --single-branch --branch ALF-2.6 --no-checkout https://github.com/ALF-QMC/ALF.git "$alf_src" > "$ALF_PREPARATION_ROOT/download.log" 2>&1
git -C "$alf_src" checkout --detach "$alf_pin" > "$ALF_PREPARATION_ROOT/checkout.log" 2>&1
test "$(git -C "$alf_src" rev-parse HEAD)" = "$alf_pin"
git -C "$alf_src" status --porcelain=v1 > "$ALF_PREPARATION_ROOT/source_status_before_build.txt"
test ! -s "$ALF_PREPARATION_ROOT/source_status_before_build.txt"
git -C "$alf_src" ls-tree -r HEAD > "$ALF_PREPARATION_ROOT/source_git_tree.txt"
cd "$alf_src"
export ALF_FLAGS_EXT=-fno-fast-math
# Upstream configure is not written for nounset; restore strict mode afterwards.
set +u
source configure.sh GNU noMPI STAB3 NO-INTERACTIVE NO-FALLBACK > "$ALF_PREPARATION_ROOT/configure.log" 2>&1
set -u
{
  printf 'ALF_FC=%s\nALF_LIB=%s\nALF_FLAGS_PROG=%s\nALF_FLAGS_MODULES=%s\n' "$ALF_FC" "$ALF_LIB" "$ALF_FLAGS_PROG" "$ALF_FLAGS_MODULES"
  gfortran --version
} > "$ALF_PREPARATION_ROOT/effective_build.txt"
timeout --signal=TERM --kill-after=15s 15m make -j1 program ana > "$ALF_PREPARATION_ROOT/build.log" 2>&1
test -x Prog/ALF.out
test -x Analysis/ana.out
sha256sum Prog/ALF.out Analysis/ana.out > "$ALF_PREPARATION_ROOT/binary_sha256.txt"
ldd Prog/ALF.out > "$ALF_PREPARATION_ROOT/linkage.txt"
# Record any generated tracked-file changes, instead of assuming a pristine build tree.
git status --porcelain=v1 > "$ALF_PREPARATION_ROOT/source_status_after_build.txt"
printf '%s\n' 'SETUP_BUILT_NOT_BENCHMARKED' > "$ALF_PREPARATION_ROOT/status.txt"
