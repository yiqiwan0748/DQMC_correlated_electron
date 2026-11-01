# Full interaction-quench curve: numerical comparison complete, workflow unresolved

**11,692 reported numerical gates passed; the overall workflow failed its final output census. Full P13 remains NOT_RUN and adds zero catalogue credit.** The ledger stays **50 PASS / 0 FAIL / 2 UNRESOLVED / 48 NOT_RUN**.

This calculation follows how entanglement changes after suddenly weakening electron repulsion in a small Hubbard chain. It compares the user's Hamiltonian and entropy-measurement route with an independently constructed TeNPy model, initial state, evolution and entropy calculation. The saved comparison covers the complete registered time curve and a growth slope over a fixed window.

The open chain has six sites, three up-spin and three down-spin electrons, hopping one, and zero chemical potential. Repulsion changes from eight to four. The initial state is the ground state in that fixed spin/particle sector, and the second Renyi entropy is measured across the central three-site cut. Both routes use complete dense spectra for direct evolution at each of 33 times from zero to four, spaced by one eighth. This spacing is an observation grid, with no chained time integration. Both registered MPS conversion controls are retained.

| Time | Saved native primary entropy center | Saved TeNPy primary entropy center |
| --- | ---: | ---: |
| 0 | 0.7934674920595309 | 0.7934674920595287 |
| 1 | 1.149203680594738 | 1.1492036805947414 |
| 2 | 1.0405550160669603 | 1.0405550160669594 |
| 3 | 1.1886687142947603 | 1.1886687142947596 |
| 4 | 0.9441889525934604 | 0.9441889525934585 |

These are stored centers. The full data retain their uncertainty radii and all registered acceptance gates; center agreement alone is not the acceptance criterion. `time_series.csv` copies all 132 saved backend/control/time entropy records, including radii and actual bond ranks.

The four saved growth records use the predetermined nine points from time 0.5 through 1.5, with an intercept and exact rational weights. The native primary slope is 0.0883793939449672 with radius 3.908370727918892e-9; the TeNPy primary slope is 0.08837939394497095 with radius 5.70672623643643e-10. These are finite-window slopes, with no universal or asymptotic rate claim. The JSON retains both controls for both routes.

All 11,692 original scientific gate occurrences, their order, limits and reported values are preserved. Coverage includes four model certificates, two initial-state certificates, 132 state certificates and four fixed-window growth records. The separate saved-association audit recorded 26,498 passing checks and zero failed association checks; an independent saved-receipt review passed 34 focused checks. These audits check recorded source/data/runtime/input/limit relationships and do not independently repeat the physical certificates.

The numerical report originally says `PASS_FULL_PENDING_INDEPENDENT_AUDIT`; that value is retained unchanged in the scientific projection. The scheduler outcome is `FAILED 2:0`, and the parent workflow remains `UNRESOLVED` because its final output census exceeded its time allowance after the three numerical roles returned successfully. Scheduler success, error-free parent completion and the original output-census publication remain **three unmet workflow predicates**. A separate complete byte-preservation capture retained the original outputs and empty directories. It does not retroactively pass those completion predicates. `WORKFLOW_STATUS.json` records this distinction explicitly.

The native route here uses the user's Hubbard MPO builder and Renyi-entropy measurement API, with new dense-spectrum, complex-state conversion and validation adapters. ITensors/ITensorMPS and TeNPy provide the underlying tensor and numerical routines. This record does not establish production DMRG time evolution, a new user-developed optimizer, a speedup, or whole-package reliability.

The eight numerical sources and frozen protocol are byte-identical to the [previous public source record](../../runtime_validation_20261028/p13_quench_pilot_rayleigh_recomparison/source/). `SOURCE_REFERENCES.json` gives their paths and hashes. Historical prospective status/resource wording in that protocol is retained as source history; this record states the actual completed-comparison and failed-workflow outcomes. The [original precision-unresolved pilot and its 15 failed gates](../../runtime_validation_20261028/p13_quench_pilot_precision_unresolved/README.md) and subsequent versioned records remain preserved.

The certificate uses conditional floating error accounting, including the recorded extended-precision assumptions; it is not directed-rounding interval arithmetic. Public binary64 scalar exports do not independently reconstruct every internal extended-precision enclosure. The packet omits raw producer tensors/eigensystems, dependency locks, execution drivers and private execution qualification, so it is not a complete standalone rerun. See `AUDIT_SCOPE.md` for the evidence boundary.

Run `python3 verify_packet.py` to check public hashes, finite data, recorded statuses and the copied CSV fields. With the repository available, `python3 verify_packet.py --repository-root /path/to/repository` also verifies the eight referenced numerical files. The checker performs no scientific package import, physical computation or network action.
