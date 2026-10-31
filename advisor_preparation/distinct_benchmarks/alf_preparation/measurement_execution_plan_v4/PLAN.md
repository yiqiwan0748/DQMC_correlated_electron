# P21 adapter-proof execution plan v4 — wait for root's resource confirmation

**P21 remains NOT_RUN_AS_DEFINED.** Execute revision4 only after independent source review and root's fresh resource-lane confirmation. V2job6458532 stopped on a stale compiler module; v3job6458622 cleaned correctly but stopped on a proof-main reference to private Obs_scal. V4 uses public ham%Pr_obs and a proof-only Hubbard print override; native sampled printing delegates unchanged to Pr_obs_base. The supported legacy writer/parser, full clean rebuild, targeted known-module cleanup, linkage checks, propagation, both proof drivers and all tolerances are unchanged. Native and v1–v3 trees/receipts are preserved. This is eligibility only, adding zero scientific completions.

The observer source is frozen in `measurement_preparation_v4/pre_execution_receipt.json`. The self-contained archive in this directory contains exactly those29 files plus that receipt, rooted at `alf_preparation/measurement_preparation_v4/`:

| Item | Exact value |
|---|---|
| Local bundle | `alf-measurement-observers-v4-source.tar.gz` |
| Bytes | 109,802 |
| SHA256 | `961f9700e6ac41efd2ad778f073cd0ae8dd75a82bee56d85c505a0145633eb7e` |
| ISAAC campaign root `R` | `/lustre/isaac24/scratch/ywang287/advisor-validation-20261022` |
| Proposed staged archive | `$R/source_releases/alf-measurement-observers-v4-source.tar.gz` |
| New extraction root | `$R/source/p21_measurement_proof_v4` |
| New result/build root | `$R/results/p21_adapter_proof_v4_attempt01` |
| Existing native ALF tree | `$R/external/alf_build_attempt01/alf_build/ALF` |
| Existing LP64 compatibility directory | `$R/external/alf_build_attempt01/compat` |
| User source | `$R/source/renyi_disjoint_corrected` (its kernel must have frozen SHA `b72bbbf…ada4be`) |

Submission receipt6458258 identifies `external/alf_build_attempt01`. The CV agent's read-only completion check reported exit0, elapsed1m11s, outer marker `ALF_BUILT_NOT_BENCHMARKED` and nested marker `SETUP_BUILT_NOT_BENCHMARKED`. These are setup results only. The commands below require both markers, executable and recorded binary hashes again; they also check its actual linked provider before proceeding.

Root should stage the archive and `bundle_receipt.json`, then run the following **only after confirming a free resource lane**, inside a one-CPU,3500MB,22-minute allocation. No second setup/science process is spawned by the wrapper. Existing source trees and results are preserved; all new paths must be unused.

```bash
set -euo pipefail
export VALIDATION_ROOT=/lustre/isaac24/scratch/ywang287/advisor-validation-20261022
R="$VALIDATION_ROOT"
archive="$R/source_releases/alf-measurement-observers-v4-source.tar.gz"
echo '961f9700e6ac41efd2ad778f073cd0ae8dd75a82bee56d85c505a0145633eb7e  '"$archive" | sha256sum -c -
mkdir "$R/source/p21_measurement_proof_v4"
tar -xzf "$archive" -C "$R/source/p21_measurement_proof_v4"
observer_root="$R/source/p21_measurement_proof_v4/alf_preparation/measurement_preparation_v4"
export ALF_NATIVE_ROOT="$R/external/alf_build_attempt01/alf_build/ALF"
export P21_LINKER_COMPAT="$R/external/alf_build_attempt01/compat"
export P21_USER_BASE="$R/source/renyi_disjoint_corrected"
export P21_PROOF_ROOT="$R/results/p21_adapter_proof_v4_attempt01"
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
python3 "$observer_root/check_native_linkage.py" "$ALF_NATIVE_ROOT/Prog/ALF.out" "$provider" "$R/results/p21_native_linkage_v4_attempt01.json"
bash "$observer_root/run_adapter_proofs.sh"
```

The native prerequisites are GNU10.2.0 plus the LP64 OpenBLAS0.3.12 library already probed by6458231, one thread, serial/noHDF5/STAB3, and pinned ALF commit `f4c500cc5d7bcf089400c63e4e75382c558a4ed3`. The new build makes legacy vector output explicit and appends `-fno-fast-math`. Every BLAS/LAPACK SONAME must resolve to that same provider. The original and observer binaries are separately hashed and linked-provider checked. The existing QuSpin environment supplies NumPy/SciPy for deterministic proof algebra only; it is not substituted for ALF as the P21 external reference.

The wrapper's hard timeouts are60seconds for cleaning and15minutes for the one-process full rebuild,90seconds for the native fixed-configuration traversal, and180seconds for the user eligibility driver. It requests no allocation itself. The previous native full build took1m11s, but the new proof workload has no measured runtime yet;22minutes is the overall cap, not a promised duration.

Expected evidence is `overlay/overlay_receipt.json`, configure/clean/build logs, known_compiler_module_cleanup.json with the removed artifact SHA/size, and its cleanup log and effective flags, `observer_linkage_audit.json`, raw native `fixed_configuration/P21_*` files, `alf_adapter_proof.json`, and `user_adapter_proof.toml`. Successful completion writes `ADAPTER_PROOFS_PASS_P21_EXTERNAL_BENCHMARK_NOT_RUN`. The native proof must report zero MC sweeps; the user proof makes32 diagnostic updates and requires exact source/state/RNG replay, complete time grids, and every contemporaneous pre/post buffer identity within the frozen tolerances. Any missing output, compile/API failure, incomplete grid, convention mismatch, or failed check preserves the attempt as unresolved and stops progression.

After these eligibility proofs pass, a separate fixed-budget statistical protocol and independent seeds are still required before root submits an interacting ALF-versus-user sampling comparison. All1408 real components must remain visible. No proof-only success may be relabeled as a P21 scientific pass.

Revision4 diagnosis and exact native-bin grammar are in `measurement_preparation_v4/COMPILE_FAILURE_AND_FIX.md`; the newly retained pinned Hamiltonian_main_mod.F90 reference proves private-storage/public-dispatch scope. Compilation has not been demonstrated for this revision. Independently review the30-member source closure before any root submission. The old statistical candidate referencesv3; sampling will require a fresh closure referencing the successful proof revision and actual sampling implementation, without editing its frozen predecessor.
