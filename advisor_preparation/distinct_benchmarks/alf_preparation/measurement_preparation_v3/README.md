# P21 native unequal-time measurement preparation — revision3

**P21 is NOT_RUN_AS_DEFINED.** Revision3 preserves v1 and the failed v2 attempt6458532. After cleanall, it hashes/logs and removes only the known untracked compiler artifact Prog/opt_time_dependent_mod.mod in the NEW copied ALF tree, then retains the strict no-stale-object/module/archive/binary assertion and full flagged rebuild. The v2 attempt stopped before physics; no scientific gate is relaxed. The overlays and proof commands are prepared and syntax checked. Native ALF compilation, both adapter eligibility proofs, and an independently frozen statistical protocol are still required. The completed LP64 setup probe establishes the available ABI/provider; it is not an interacting ALF benchmark.

The native ALF overlay exports both physical spin Green matrices at all eleven times, retaining all 704 complex entries. It reads the actual buffers passed to ALF's existing `ObserT`; `Tau_m`, stabilization, hopping, interaction, and sampling kernels remain unchanged. The user overlay adds callbacks around the actual ordinary replica-update path, retaining maintained buffers **before** stabilization. See [CORRESPONDENCE.md](CORRESPONDENCE.md) for the operator, slice order, exact Gaussian/Hirsch identity, raw-bin layout, and moving-cut endpoint treatment.

`overlay_v1/` contains the complete modified source files and minimal diffs. `prepare_observers.py` regenerates precisely those overlays from pinned originals in a new directory, verifying their hashes without editing the input trees. `upstream_reference/` preserves seven inspected ALF source files byte-for-byte at the pinned Git blobs; `upstream_files.json` records their upstream blob SHA, SHA256, sizes and source URLs. Original GPL and attribution/origin notices are preserved; `license.gpl` and `license.additional` accompany this clearly marked modification. This is not an official ALF version. The previous preparation and setup_v2 files are unchanged.

The two eligibility proofs are deliberately bounded:

- Native ALF uses one deterministic nonuniform HS configuration and **zero Monte Carlo sweeps**. It exports actual `PROPR`, `PROPRM1`, and hopping matrices, invokes the actual native `Tau_m` traversal, then stops before sampling. `check_alf_proof.py` verifies all entries, endpoints, exact interaction identity, full time counts, and the distinguishable `VK` order.
- The user proof makes 32 diagnostic update calls (sixteen forward/reverse cycles), checking every maintained pre/post matrix against that callback's current HS configuration. Original and observed trajectories must match exactly, including fields, Greens, acceptance, determinants/signs, LDR factors, direction, and next RNG value. `UserMatrixObserver.jl` also demonstrates the primary raw paired-bin export without replacing buffers with recomputed matrices.

No scientific submission is performed by these files. Root execution, once the native build is available and the setup lane is free, can use the following command **inside a one-CPU,3500MB,22-minute allocation**. Supply the two actual setup paths rather than guessing an uncompleted build location:

```bash
export VALIDATION_ROOT=/lustre/isaac24/scratch/ywang287/advisor-validation-20261022
export ALF_NATIVE_ROOT=/actual/completed/alf_build/ALF
export P21_LINKER_COMPAT=/actual/completed/linker_setup/compat
export P21_USER_BASE="$VALIDATION_ROOT/source/renyi_corrected"
export P21_PROOF_ROOT="$VALIDATION_ROOT/results/p21_adapter_proof_attempt01"
export P21_PYTHON="$VALIDATION_ROOT/quspin-venv/bin/python"
export P21_JULIA="$VALIDATION_ROOT/runtime/julia-1.12.1/bin/julia"
export JULIA_DEPOT_PATH="$VALIDATION_ROOT/depot"
bash "$VALIDATION_ROOT/source/distinct_benchmarks/alf_preparation/measurement_preparation_v3/run_adapter_proofs.sh"
```

The wrapper copies the original native build to a new tree, applies only the reviewed observers, configures serial GNU/STAB3/legacy-vector output, cleans only that new tree and rebuilds all objects with one make process, and checks that all BLAS/LAPACK names resolve to the one pinned OpenBLAS provider. It then runs the two eligibility probes with fixed timeouts. Existing source, binaries, results and frozen scripts are preserved. Output directories must be new; a compiler/API or convention mismatch stops the attempt and preserves its logs. The checker requires SciPy/NumPy only from the already installed pinned QuSpin environment. Source and license copies may be published; compiled binaries and compatibility-library symlinks need not be copied into the public record.

Even `ADAPTER_PROOFS_PASS_P21_EXTERNAL_BENCHMARK_NOT_RUN` is only an eligibility result. Later science must compare independently sampled native ALF and actual-user outputs with frozen disjoint seeds, fixed burn/bin budgets, all 1408 real components, quality checks, and simultaneous equivalence intervals. The candidate eight-chain/family-alpha.01/margin.02 settings in the earlier plan are **not** a completed statistical freeze or a scientific pass. This preparation covers the declared real-hopping ordinary replica-update path; general orbital mixing, complex hopping and Grover-fused paths require further coverage.
