# Distinct scientific benchmark ledger

**22/100 PASS; 1 FAIL; 2 UNRESOLVED; 75 NOT_RUN.**

Exactly 100 catalogue IDs; one status per distinct idea. Parameters, matrix elements, fixture/check/component counts are not additional ideas.

ED/thermal-ED analysis and production free-band construction only for this tranche. No production DQMC sampling or DMRG pass is inferred.

Every raw input is copied byte-for-byte under `inputs/` and hashed in `input_manifest.json`. Inputs and the catalogue are never modified. PASS requires the frozen-source provenance and a completed scientific comparison. FAIL and incomplete/unresolved outcomes remain visible. A11 is summarized by its actual L1 metric, never the analytical driver’s empty generic maximum-absolute field.

| Category | PASS | FAIL | UNRESOLVED | NOT_RUN |
|---|---:|---:|---:|---:|
| Analytical | 16 | 0 | 0 | 9 |
| External-reference error control | 0 | 0 | 0 | 25 |
| Independent package | 6 | 1 | 0 | 18 |
| Published journal target | 0 | 0 | 2 | 23 |

## Executed and unresolved ideas

| ID | Status | Target scope | Metric-specific maxima |
|---|---|---|---|
| A01 | PASS | user ED and/or thermal ED analysis functions plus explicit analytical-benchmark adapters | maximum_absolute: 5.55112e-16 |
| A02 | PASS | user ED and/or thermal ED analysis functions plus explicit analytical-benchmark adapters | maximum_absolute: 1.66533e-16 |
| A03 | PASS | user ED and/or thermal ED analysis functions plus explicit analytical-benchmark adapters | maximum_absolute: 7.10543e-15 |
| A04 | PASS | user ED and/or thermal ED analysis functions plus explicit analytical-benchmark adapters | maximum_absolute: 1.28786e-14 |
| A05 | PASS | user ED and/or thermal ED analysis functions plus explicit analytical-benchmark adapters | maximum_absolute: 7.10543e-15 |
| A06 | PASS | user ED and/or thermal ED analysis functions plus explicit analytical-benchmark adapters | maximum_absolute: 3.41949e-14 |
| A07 | PASS | user ED and/or thermal ED analysis functions plus explicit analytical-benchmark adapters | maximum_absolute: 2.22045e-15 |
| A08 | PASS | user ED and/or thermal ED analysis functions plus explicit analytical-benchmark adapters | maximum_absolute: 1.89848e-14 |
| A09 | PASS | user ED and/or thermal ED analysis functions plus explicit analytical-benchmark adapters | maximum_absolute: 1.65961e-12 |
| A10 | PASS | user ED and/or thermal ED analysis functions plus explicit analytical-benchmark adapters | maximum_absolute: 6.05072e-15 |
| A11 | PASS | user ED and/or thermal ED analysis functions plus explicit analytical-benchmark adapters | L1: 1.44869e-14 |
| A12 | PASS | user ED and/or thermal ED analysis functions plus explicit analytical-benchmark adapters | maximum_absolute: 4.44089e-16 |
| A13 | PASS | user ED and/or thermal ED analysis functions plus explicit analytical-benchmark adapters | maximum_absolute: 3.46945e-18 |
| A14 | PASS | user ED and/or thermal ED analysis functions plus explicit analytical-benchmark adapters | maximum_absolute: 8.88178e-15 |
| A22 | PASS | user ED and/or thermal ED analysis functions plus explicit analytical-benchmark adapters | maximum_absolute: 4.44089e-16 |
| A23 | PASS | user ED and/or thermal ED analysis functions plus explicit analytical-benchmark adapters | maximum_absolute: 1.0001e-14 |
| P01 | PASS | user ED state/Hamiltonian plus new observable adapter versus independent QuSpin | maximum_absolute: 1.06581e-14 |
| P02 | PASS | user ED state/Hamiltonian plus new observable adapter versus independent QuSpin | maximum_absolute: 5.32907e-15 |
| P03 | PASS | user ED state/Hamiltonian plus new observable adapter versus independent QuSpin | maximum_absolute: 7.21645e-15 |
| P04 | PASS | user ED state/Hamiltonian plus new observable adapter versus independent QuSpin | maximum_absolute: 1.19349e-15 |
| P05 | PASS | user ED state/Hamiltonian plus new observable adapter versus independent QuSpin | maximum_absolute: 9.99201e-16 |
| P06 | PASS | user ED state/Hamiltonian plus new observable adapter versus independent QuSpin | maximum_absolute: 3.44169e-15 |
| P16 | FAIL | user ED state/Hamiltonian plus new observable adapter versus independent QuSpin | maximum_absolute: 0.0472809 |
| J20 | UNRESOLVED | Completion not established | No completed numeric comparison |
| J21 | UNRESOLVED | Completion not established | No completed numeric comparison |

The complete 100-row ledger and all explanations are in [ledger_100.csv](ledger_100.csv) and [ledger_100.json](ledger_100.json). Fixture and scalar-component totals are descriptive evidence, not pass-count multipliers. Existing user analysis functions and new observable adapters are explicitly identified in the JSON; only J20/J21 exercise the original production band/real-space builders.

## Provenance

- Frozen release: `tranche01`
- Frozen archive SHA256: `a5c95efd46268165a9556a7a56371651257f2e66da63850d2ba35a9bcec9d09c`
- Catalogue SHA256: `59a780c5994db50a56da260fc51ae639f22b674b3bc24f878cde1a6e7bc31818`
- Reporting tool SHA256: `e131a57d9f8e3c2c1ef3951b3b8d18a9049fe2688023263c87bd9f285465d49e`

Audit notes:

- Malformed journal; exact bytes preserved: Invalid value (at end of document)
