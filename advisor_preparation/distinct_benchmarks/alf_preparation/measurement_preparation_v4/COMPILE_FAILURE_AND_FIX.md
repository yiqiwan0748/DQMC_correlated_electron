# Revision3 compiler failure and revision4 print-interface repair

P21 remains **NOT_RUN_AS_DEFINED**. ISAAC job6458622 ended FAILED/2:0 after48seconds, with88112KiB maximum batch RSS. Its clean guard passed and library/program prerequisites compiled, then compilation stopped at `main.F90:804`: `Obs_scal` had no implicit type. No observer executable, fixed-configuration proof, user-buffer proof, or scientific sample was produced. [Original text receipts](../../results/alf_adapter_proof_v3_6458622/transfer_hashes.json) preserve nine explicitly selected files with byte counts and remote/local SHA256 matches; the scheduler/hash receipt and all original log text are retained. Native ALF and preparations1–3 are untouched.

The retrieved build log contains exactly **one Fortran Error diagnostic**, followed by the two ordinary make error exits. Its failing line was `Call Print_bin_Vec(Obs_scal(I),Group_Comm)` in the proof branch of our additive main overlay. The clean marker `COPIED_OBJECTS_MODULES_ARCHIVES_REMOVED_PASS` and the targeted compiler-module cleanup receipt show that revision2's stale-object setup problem was resolved; it is not the cause of this failure.

## Pinned source evidence and exact fix

At ALF commit `f4c500cc5d7bcf089400c63e4e75382c558a4ed3`, [Hamiltonian_main_mod.F90](https://github.com/ALF-QMC/ALF/blob/f4c500cc5d7bcf089400c63e4e75382c558a4ed3/Prog/Hamiltonian_main_mod.F90) declares module-default `private` and keeps `Obs_scal` private, except for an unrelated `__PGI` conditional. Public `ham` has a `Pr_obs` type-bound procedure; `Pr_obs_base` prints its private observables within their valid parent-module scope. The Hubbard submodule already accesses those parent-private observables through Fortran submodule host association. This exact additional reference is preserved byte-for-byte at Git blob `13ec1174618fcc50e2607ca3227710760d9619e5` in `upstream_reference/`, including GPL/origin notices.

Revision4 removes the proof main's direct `Use Observables`/`Obs_scal` access. It calls the public `ham%Pr_obs(Ltau)` interface instead. The additive Hubbard override `Pr_obs => P21_Pr_obs` has the same `nopass`/integer-input signature as its native base. In proof mode it checks dimensions and positive measurement counts, then calls the unchanged native writer for only P21 vectors5–7. Those are measured by the actual native `Tau_m` traversal. Other observables have zero measurements in this deliberately fixed-configuration traversal, so printing every allocated observable through the unmodified base would divide some of them by zero. In ordinary sampled mode the override delegates directly to `Pr_obs_base(Ltau)`, preserving ALF's original complete printing behavior. It never makes private storage public, defines `__PGI`, substitutes independently recomputed matrices, or changes sampling/propagation/stabilization/HS algebra.

## Legacy output is supported; its format is unchanged

The [pinned native observables module](https://github.com/ALF-QMC/ALF/blob/f4c500cc5d7bcf089400c63e4e75382c558a4ed3/Prog/observables_mod.F90) explicitly defines `OBS_LEGACY` locally when neither HDF5 nor that macro is selected (lines45–47). The successful native build used GNU10.2.0/noMPI/noHDF5/STAB3 and appended `-fno-fast-math`; revision3 used the same settings plus the redundant explicit `-DOBS_LEGACY`. The actual compiler successfully compiled `observables_mod.F90` with that flag before it reached main. Thus the observed failure is a visibility error, not evidence of an obsolete output mode. Revision4 preserves those flags and the exact existing writer/parser. No output-format migration or changed proof tolerance is needed.

Native `Print_bin_Vec` divides the accumulated complex numerator vector and real sign by `Obs%N`. Each appended `_scal` line is:

```text
I10(width+1) [width repetitions of " (" E25.17E3 "," E25.17E3 ")"] E26.17E3(mean_sign)
```

The leading integer counts the complex-vector entries plus one real sign. Each complex pair is a sign-weighted numerator **mean**, not a ratio already divided by the sign. `P21_GT0_scal` and `P21_G0T_scal` each have352 complex entries and leading marker353. `P21_Counts_scal` has11 complex raw-count entries, marker12, print divisor1 and final sign1. In the one-configuration proof each count is exactly1. The existing checker requires exactly one line per vector, parses all complex pairs, checks marker/counts/sign, reshapes tau-major/spin/row/column, and divides the two matrix-vector means by the sign. It continues to require the full time grid, both orders/endpoints, native `VK` propagation and all unchanged tolerances. Statistical sampling later must restore paired sums using raw counts; no per-bin-ratio averaging is introduced.

## Limits of this repair record

The repair is generated and statically reviewed; it has **not** compiled or executed in revision4. There is no local Fortran compiler. A new one-CPU3500MB bounded root-controlled attempt must still rebuild the fresh copied tree and pass both original eligibility proofs. It may uncover later compiler/API/convention issues that the stopped revision3 build could not reach. All prior failures remain separate receipts; no proof-only result establishes P21 science.
