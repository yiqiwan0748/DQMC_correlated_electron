This versioned setup wrapper uses the **actual ISAAC module inventory**: `gcc/10.2.0` and `openblas/0.3.12`, built together for RHEL8. The initial inventory job6458092 instead found default system gfortran8.5, make4.2.1 and Python3.8.8 with no loaded modules. The initial setup files and their preparation receipt remain unchanged; this directory adds an explicit compatible-toolchain candidate. No execution or installation by the preparation agent is implied.

The wrapper defaults to **probe only**. It loads those two modules, records compiler/search paths, module definitions, node CPU and library files, and queries the actual OpenBLAS library with Python's standard-library ctypes. The [OpenBLAS0.3.12 source](https://github.com/OpenMathLib/OpenBLAS/blob/v0.3.12/driver/others/openblas_get_config.c) includes `USE64BITINT` in its configuration string for an ILP64 build and `NO_LAPACK` when LAPACK is omitted. Either flag rejects this candidate. The exact library SHA256, returned configuration, required DDOT/DGEMM/ZGEMM/ZHEEV symbols and one-thread setting are recorded.

The source provider identified by the actual module is

```
/spack/spack-0.16.3/apps/linux-rhel8-cascadelake/gcc-10.2.0/openblas-0.3.12-c4x3eon3lkfddumq466n3n6b5gpqdbvh/lib/libopenblas.so
```

ALF's GNU configure script explicitly probes `-llapack -lblas`. If the Spack module supplies only `libopenblas.so`, changing ALF_LIB after configuration would be too late to satisfy that probe. The wrapper creates **local compatibility symlinks** `compat/liblapack.so` and `compat/libblas.so` pointing to the *same metadata-verified OpenBLAS provider*. LIBRARY_PATH exposes them to the normal compiler linker; the module library path supplies the actual runtime SONAME. No system library, upstream source or prior environment receipt is modified. Dynamic linkage must show that same provider for all BLAS/LAPACK entries.

The standalone2×2 Fortran program checks default32-bit INTEGER, real and complex matrix multiplication, DDOT, both Hermitian eigenvalues, INFO, eigenpair residuals and orthogonality. Zero-high-word padded integer inputs remain bounded if an unexpected64-bit routine is encountered; an INFO canary detects an8-byte LAPACK write. The prior metadata gate supplies the decisive OpenBLAS integer-width check. Compile and runtime are capped at60seconds and10seconds respectively. This is a compiler/ABI setup probe, not a physical benchmark or P21 completion.

After staging, the exact root-owned command for the probe is:

```bash
sbatch --time=00:03:00 \
  --export=ALL,VALIDATION_ROOT="$VALIDATION_ROOT",ALF_LINKER_ROOT="$VALIDATION_ROOT/external/alf_linker_probe01",ALF_SETUP_V2_MODE=probe \
  "$VALIDATION_ROOT/source/distinct_benchmarks/alf_preparation/setup_v2/probe_and_build_candidate.sbatch"
```

A successful receipt reports `LP64_SETUP_PROBE_PASS_NO_ALF_SCIENCE`. It has run no ALF code. The operator may then choose a **different new output directory** and `ALF_SETUP_V2_MODE=build` with the20-minute sbatch cap. Build mode repeats the verified probe, checks the preserved v1 setup script SHA, then calls that original candidate to clone the immutable ALF2.6 commit and perform its15-minute capped build. A completed build reports `ALF_BUILT_NOT_BENCHMARKED`; it still never runs ALF.out. One CPU and3500MB are the upper resource limits in both modes.

The native continuous-HS observable adapters, user unequal-time observer and slice-order proof remain separate scientific prerequisites. P21 must remain NOT_RUN_AS_DEFINED after a successful compiler probe or build.
