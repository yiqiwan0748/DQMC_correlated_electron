# Validated working package

`HubbardRenyiDQMC/` contains the current working source. Exact historical and
intermediate source versions remain in `../validation_source/`; numerical
claims refer to the hashes in their own execution receipts.

The October validation found and corrected two independent defects:

- The isolated-atom systematic sampling scan could be reducible. Independent
  50% holds before each local proposal restore the intended sampling kernel;
  the exact finite-field diagnosis and corrected external comparisons are in
  `../external_oracles/ATOMIC_SAMPLER_DIAGNOSIS.md`.
- Spatial ED entropy for a noncontiguous subsystem omitted within-spin
  fermionic permutation signs. The minimal repair preserves the sparse
  contraction and adds the missing occupied B-before-A crossing parities.
  The Hamiltonian and Lanczos solver are unchanged. All seven registered
  QuSpin comparisons and all 16 four-site spatial subsets passed through the
  actual repaired source on ISAAC job 6457716. The independent permutation
  proof, original failure, exact patch and raw retest results are under
  `../distinct_benchmarks/packages/p16_diagnosis/` and
  `../distinct_benchmarks/results/packages_disjoint_repair_attempt01/`.

The repaired ED source SHA256 is
`f8f94cfc13e11bb8b5e64e8ac8f4ade5ed3a13b3976c83b5ab5db89d86ba5f28`.
The original ED/thermal analysis pair is also preserved separately under
`../validation_source/ed_before_disjoint_repair/` so the original tranche can
be audited independently of later working-source changes.

These finite fixtures do not certify every model, observable, numerical regime
or unsupported measurement API. The distinct-benchmark ledger records method
coverage and unexecuted questions explicitly.
