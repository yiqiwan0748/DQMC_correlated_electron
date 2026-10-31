# Saved boundary-spin-response comparison completed; full benchmark awaits reconciliation

The saved-input comparison now completes its frozen pilot grid: the native and TeNPy DMRG outputs agree with the QuSpin reference under all **1,393 required gates**, with **zero FAIL gates**. It checks the response to a weak boundary magnetic field, including the spin profiles and linear-response comparisons from the saved states. This is a **complete saved component pending full reconciliation**, with **zero new catalogue credit**. Full production P17 remains **NOT_RUN**, unadmitted, with resource limits unset.

The complete census contains 277 reference gates, 6 spectral-response gates, 513 native gates, 429 TeNPy gates and 168 cross-state gates. All 12 state certificates completed: 6 native and 6 TeNPy. The analysis also completed 21 reference-sector certificates, 6 zero-field spectral responses and qualification of 6 unique operator payloads. The cross-state total comprises 108 profile-pair, 24 linear-window and 36 response-pair gates. These are internal checks of one component, not 1,393 distinct benchmark ideas.

No producer was rerun: this attempt reuses the three original solver outputs and preserves the comparison mathematics, precision/error rules, physical gates and five-million-call ceiling. Its prospective comparison budget was 1,500 seconds within a 30-minute, one-CPU, 3500 MB allocation. The comparison child returned zero after 720.097 seconds, and the parent finished in 727.048 seconds. The scheduler records **COMPLETED / 0:0**, with 728 seconds elapsed. The prior 420-second profiling timeout remains preserved in [the earlier packet](../p17_profiled_comparison_v1_20261027/README.md), including its 690 partial PASS gates and incomplete status. Earlier global-pilot and diagnostic failures remain unchanged.

The saved profile contains **1,758 events** and **419,247 arithmetic-function calls**. These calls operate on arrays of different sizes; they are not hardware floating-point-operation counts.

| Selected completed verification scope | Intervals | Wall seconds | Self-CPU seconds |
| --- | ---: | ---: | ---: |
| Runtime inventory hashing | 1 | 88.412 | 1.949 |
| All reference work | 1 | 12.755 | 12.338 |
| Represented-MPO contractions | 56 | 339.136 | 335.999 |
| State reconstruction | 12 | 258.867 | 257.591 |

The state-reconstruction intervals include 129.040 seconds of extended-MPS contraction and 129.263 seconds of binary64-export-budget contraction. Those subintervals are **nested** within state reconstruction. The MPO contractions are nested within operator qualification; reference, operator, state and cross-state scopes are themselves included within the unchanged comparison-run interval. Do not add parents and children. Selected scopes do not exhaust total time. No scope was interrupted in this completed attempt.

These are verification timings, with instrumentation included. No matched baseline, instrumentation-overhead estimate, exclusive low-level cause or production-solver speedup is established. The certificate uses **conditional floating-point error accounting**, not directed-rounding interval arithmetic.

An independent saved-evidence audit passed **251 checks**, with zero failures. A separate saved-runtime-metadata audit passed **29 checks**, reparsing 173 mapped segments, 33 files and 191 loaded modules in each saved before/after report. These audits checked the saved identities and exact component census; they did not recompute quantum results or freshly inspect remote providers. An initial audit-harness namespace assumption was corrected for this versioned capture; that audit-only failure is preserved separately, and no physical comparison or acceptance gate changed.

`attempt_result.json` records the audited outcome and hash provenance; `observed_profile.json` records selected inclusive timings, actual call counters and completed component coverage. Raw runtime inventories, process details, large arrays and execution-control records are excluded. The complete component must still be reconciled into a separately reviewed full production workflow. The ledger remains **49 PASS / 0 FAIL / 2 UNRESOLVED / 49 NOT_RUN**, with P17 **NOT_RUN**.
