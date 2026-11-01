# Full Hubbard-quench entanglement reconciled

**51 PASS / 0 FAIL / 2 UNRESOLVED / 47 NOT_RUN** across exactly 100 distinct questions. Only P13 changes from the immutable [50-question checkpoint](../full_boundary_reconciliation_v18/README.md). All 99 other row objects, all 50 prior PASS and both unresolved histories remain exact.

P13 follows second Renyi entropy after repulsion changes from eight to four in an open six-site Hubbard chain. The complete native–TeNPy comparison passes all 11,692 frozen gate occurrences across 132 states at 33 times and four fixed nine-point growth fits. The [scientific packet](../../../runtime_validation_20261101/p13_full_quench_audited/README.md) retains the complete curves, certificates, exact raw-state audit values and source references. All 16 published files were read back and verified before this leaf was created.

This adds one finite-cluster question, using disclosed dense-spectrum/MPS adapters and conditional floating error accounting. It does not establish a user-native real-time evolution solver, speedup, thermodynamic scaling or whole-package reliability. All earlier failures remain preserved. [Closure audit](audit.json).

Run `python3 verify_leaf.py --repository-root /path/to/repository` for the read-only saved-ledger/packet consistency check. It runs no physical calculation.
The verifier refuses optimized Python modes so its validation checks cannot be disabled.
