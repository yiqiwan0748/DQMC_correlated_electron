# Repaired interaction-quench pilot: workflow completed

**All 8,889 reported numerical gates passed, and the original three-role workflow and output census completed successfully.** This five-time qualification pilot adds zero catalogue credit. Full P13 remains **NOT_RUN**; the ledger stays **50 PASS / 0 FAIL / 2 UNRESOLVED / 48 NOT_RUN**.

This calculation follows the early change in quantum entanglement after suddenly weakening electron repulsion in a small Hubbard chain. It compares the user's Hamiltonian and Renyi-entropy measurement route with an independently constructed TeNPy model, initial state, evolution and entropy calculation.

The open chain has six sites, three up-spin and three down-spin electrons, hopping one and zero chemical potential. Repulsion changes from eight to four. The initial state is the ground state in that fixed spin/particle sector, and the second Renyi entropy is measured across the central three-site cut. Both routes use complete dense spectra for direct evolution at five times from zero through 0.5, spaced by one eighth. This is an observation grid, with no chained time integration. Both registered MPS conversion controls are retained.

| Time | Saved native primary entropy center | Saved TeNPy primary entropy center |
| --- | ---: | ---: |
| 0 | 0.7934674920595309 | 0.7934674920595287 |
| 0.125 | 0.837014112067983 | 0.8370141120679822 |
| 0.25 | 0.9427179354452171 | 0.9427179354452155 |
| 0.375 | 1.05399738616138 | 1.053997386161378 |
| 0.5 | 1.125863399063166 | 1.1258633990631646 |

These are stored centers. Agreement of centers alone is not the acceptance criterion: the full data retain radii, values and limits for all 8,889 ordered scientific gate records. `time_series.csv` copies all 20 backend/control/time entropy records, including radii and actual bond ranks. The lossless gzip JSON also retains four model certificates and two initial-state certificates.

The original pilot comparison status remains `PASS_PILOT_ZERO_CREDIT_PENDING_INDEPENDENT_AUDIT` in the scientific data. The original workflow completed with no errors and published its complete output census. A saved-association audit passed 9,321 checks with zero failures; an independent focused review of that saved receipt passed 152 checks with zero failures. These inspect recorded source/data/runtime/input/limit relationships. They do not independently repeat the physical calculations.

Five observation times do not cover the registered full curve or its fixed nine-point growth window. Growth records remain absent, and their original NOT_RUN statuses are preserved. The [earlier full-curve workflow failure](../p13_full_curve_workflow_unresolved/README.md), including its failed final output census, remains unchanged. The [original precision-unresolved pilot and its 15 failed gates](../../runtime_validation_20261028/p13_quench_pilot_precision_unresolved/README.md) remain preserved as well.

The native route uses the user's Hubbard MPO builder and Renyi-entropy API with disclosed dense-spectrum, complex-state conversion and validation adapters. ITensors/ITensorMPS and independent TeNPy provide underlying tensor and numerical routines. This record establishes neither production DMRG time evolution, a new user-developed optimizer, a speedup, nor whole-package reliability.

The eight referenced numerical files are byte-identical to the [previous public source record](../../runtime_validation_20261028/p13_quench_pilot_rayleigh_recomparison/source/). Historical prospective wording in that frozen protocol remains source history; this record states the actual pilot outcome. The bounds use conditional floating error accounting and recorded extended-precision assumptions, rather than directed-rounding interval arithmetic. Public binary64 exports do not reconstruct every internal extended-precision enclosure.

Raw producer tensors/eigensystems, dependency locks, execution drivers and private execution records are omitted, so this packet is not a complete standalone reproduction. See `AUDIT_SCOPE.md`. Run `python3 verify_packet.py` for public byte/status/CSV consistency; `--repository-root` additionally checks the eight source references. The checker performs no physical calculation, scientific import or network action.
