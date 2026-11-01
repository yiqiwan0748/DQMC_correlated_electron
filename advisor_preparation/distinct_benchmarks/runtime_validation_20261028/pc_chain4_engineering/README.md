# Four-site Hubbard MPS measurement check

**Saved-record audit: PASS, 19,094 checks and zero failures.** This is one excluded engineering arm; the full X08 benchmark remains `NOT_RUN`, with zero new catalogue credit. The audit summary is in [AUDIT_SUMMARY.json](AUDIT_SUMMARY.json).

This check compares measurements of the same saved matrix product states against an independently written exact finite-data calculation. It uses a deliberately small bond dimension, chi = 2. The original execution recorded consistency checks on the saved representations and observables; the saved-record audit authenticates those records and their documented certificates. This does not establish ground-state accuracy or convergence of the optimizer.

The model is the open four-site Hubbard chain with hopping t = 1, uncentered on-site interaction U = 4, chemical potential 0, and two up and two down electrons. Starting from alternating up/down occupations, the library optimizer performs 12 registered sweeps. Snapshots are retained after sweeps 3, 6 and 12. All controls, including the seed and noise schedule, are in [protocol.json](protocol.json).

| Saved sweep | Same-state energy | Raw state norm squared | S2, cut 1 | S2, cut 2 | S2, cut 3 | Optimizer final-bond Ritz value |
| --- | --- | --- | --- | --- | --- | --- |
| 3 | 0.0 | 0.9999999999999998 | 0.6931471805599453 | -0.0 | -0.0 | See the separate Ritz table |
| 6 | 0.0 | 1.0 | 0.6931471805599453 | -0.0 | -0.0 | See the separate Ritz table |
| 12 | 0.0 | 1.0 | 0.6931471805599453 | -0.0 | -0.0 | -0.8284271247461902 |

These are stored observations, copied without recomputing energies or entropies. The Ritz diagnostic belongs to a local optimization step before the retained MPS truncation. Its value must not be substituted for the expectation value of the saved post-truncation state. Their difference is material in this coarse arm. The recorded truncation diagnostics also remain visible, rather than being used to assert convergence. Signed `-0.0` values are preserved.

The exact reference decodes each saved binary64 coefficient as an exact dyadic rational. It independently constructs the eight-mode fermion algebra, contracts the saved tensors, checks all 256 state amplitudes, and evaluates state norms, Hamiltonian actions, spin correlations and reduced-state purities. Natural logarithms are enclosed by a finite rational series with an explicit tail bound. The reference checks the represented finite states; it supplies no independent ground-state solution. The execution recorded 2,799 passing checks, 17 state records and 21 entropy-cut records. Those are the original execution outcomes now covered by the completed saved-record audit. That audit independently rechecked recorded rational inequalities, finite logarithm enclosures, callback/scalar associations and source/byte coverage; it did not recompute the underlying physical residuals or purities. No reference replay or scientific run was repeated to prepare this packet.

The user supplied the Hubbard MPO builder and `renyi2_bond_entropy` in [the original analysis source](code/compute_s2_itensor_dmrg.jl). ITensorMPS supplied the DMRG optimizer and native observer. The measurement, capture and checkpoint adapters and [the independent exact replay](code/EngineeringCompare.jl) were newly written for this campaign. This record makes no authored-optimizer or performance-speedup claim. The complete bounded execution used one CPU, took approximately 144.66 seconds and reached about 1.62 GB of memory; these are whole-execution figures, not solver-only timings.

The [data dictionary](DATA_DICTIONARY.md) explains the tables. The [provenance record](PROVENANCE.json) records the exact copied code hashes, projection rules and limits. Raw spin components, local observables, all 72 bond diagnostics, all 12 stopping callbacks and the full recorded algebraic replay remain available. No values were clipped, averaged, symmetrized, phase fitted or converted to absolute values.

The supplied scientific JSON retains the numerical tensors, amplitudes, operator witnesses and observation records needed to inspect the finite-data reference. Execution infrastructure and opaque serialized snapshots are private. The original user source and unchanged exact-reference module are included, but the original measurement/capture runner is not included and there is no tested standalone public replay command. This is a scientific data projection, not a complete reproducibility release. No captured code was executed, binary snapshot deserialized, or scientific calculation repeated while preparing it.

Remaining work for the full benchmark includes an independent ground-state reference, all registered convergence arms and the discarded-weight/ground-state-bias questions. This one arm does not resolve those questions.

The unchanged `publication_status` field in [protocol.json](protocol.json) records the earlier preparation stage. Current audit status is in [AUDIT_SUMMARY.json](AUDIT_SUMMARY.json). Two earlier audit attempts stopped on audit-schema assumptions and remain privately retained; the final audit corrected those interpretations without changing or rerunning the original scientific code or data.
