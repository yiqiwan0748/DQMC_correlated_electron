# Boundary spin response: full P17 comparison

**PASS — one distinct scientific benchmark.** A weak magnetic field at one end of a six-site Hubbard chain produces the same spatial spin response in the user Hubbard-MPO calculation, independent TeNPy DMRG and QuSpin exact diagonalization. All **6,578 frozen gates pass**, with zero failures.

The chain has six particles, open ends, hopping 1 and repulsion 4. The edge term is `-h Sz`, with `Sz=(n_up-n_down)/2`. The fixed fields are 0, ±0.02, ±0.01 and ±0.005. Both MPS backends use two truncation arms and two initial states: 56 states in total. QuSpin provides complete eigensystems for all 49 field/spin-sector combinations. These variants cover P17 and do not add catalogue ideas.

| Comparison | Required checks | Largest recorded difference including its radius | Frozen limit |
| --- | ---: | ---: | ---: |
| Spatial spin profiles |1176|3.607e-13|2e-7|
| Central-difference spin response |504|6.824e-11|1e-4|
| Finite-field response versus zero-field spectral response |144|3.591e-4|1e-3|

The full run completed in 69 min 05 s on one CPU with 3500 MB requested memory. This is a verification timing, with no production-solver speedup claim. ITensorMPS 0.4.1 and TeNPy 1.1.1 supply their respective optimizers; QuSpin 1.0.1 supplies the independent exact reference. The user code supplies the Hubbard MPO. Boundary-field, state-export, contraction and certificate adapters are new validation code.

Independent review passed 807 capture checks, 963 saved-metadata checks, 9688 saved-runtime checks and 17002 scientific receipt checks. All eight before/after runtime observations and complete auditor source inventories were checked. The [scientific projection](full_comparison_scientific_projection.json) preserves every original gate, all 56 state certificates, all 49 reference-sector certificates and all six zero-field spectral responses. Its original pending-review/zero-credit fields are retained verbatim; the separate [reconciliation summary](scientific_audit_summary.json) records the final decision. [Machine-readable results](attempt_result.json) retain exact maxima and limits.

The uncertainty budgets are conditional floating arithmetic, not directed-rounding interval arithmetic. Saved-receipt audits check source/runtime association and recorded controls, configurations and necessary scalar relationships; they do not repeat every Hamiltonian, tensor contraction or eigensystem calculation. Boolean-only intermediates, elementwise comparisons and the ground-distance denominator have the explicit proof limits listed in the reconciliation summary. All 22 negative subtractive-variance diagnostics remain unchanged; the positive H-residual certificate is separate. This validates one finite-cluster ground-state response and makes no thermodynamic-order or whole-package reliability claim.

The [audit history](audit_history.json) retains the configuration-attribution and capture-reader audit failures with their narrow repairs. The [earlier saved comparison component](../p17_profiled_comparison_complete_v2_20261027/README.md), earlier 420-second timeout and original global-pilot failures remain separate and unchanged. Runtime/source/synthetic checks add zero catalogue ideas. Only P17 gains one completion; the new ledger has 50 PASS, 2 UNRESOLVED and 48 NOT_RUN.
