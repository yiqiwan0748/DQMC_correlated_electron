# Measured comparison cost; complete comparison still unresolved

The new profiling attempt reused the three original solver outputs and preserved the original comparison mathematics, precision/error rules, physical gates and five-million-call ceiling. It reached the **420-second comparison cap**. The parent and child report **UNRESOLVED**, and Slurm records **FAILED / 2:0**. The comparison did not finish; it adds **zero catalogue credit**.

It saved **1,098 timing events**, actual arithmetic-function counters and the partial scientific report. There are **690 recorded PASS gates and zero FAIL gates**:277 reference,6 spectral-response and407 native. Four native states completed (zero and negative pilot field, both starts); no TeNPy state completed. All21 reference-sector certificates and6 zero-field responses were saved. The108 profile-pair,24 linear-window and36 response-pair gates remain pending. Producer-file completeness is separate from completion of this analysis.

Selected measured verification scopes:

| Scope | Completed scopes | Wall seconds | Self-CPU seconds |
| --- | ---: | ---: | ---: |
| Runtime inventory hashing | 1 | 69.238 | 4.536 |
| All reference work | 1 | 11.780 | 11.693 |
| Represented-MPO contractions | 35 | 203.579 | 202.291 |
| State reconstruction | 4 | 112.435 | 111.688 |

The four state-reconstruction intervals contain56.033seconds of extended-MPS contraction and56.197seconds of binary64-export-budget contraction. Those subintervals are **nested** and must not be added again to their parent. All35 MPO contractions closed normally. Interruption occurred during the fifth native state's extended-MPS contraction scope,9.623wall/9.537CPU seconds and26 tracked calls after that scope's start marker. Its operator/sector/row context had correctly been cleared.

The saved counter total is **256,666 arithmetic-function calls**, including249,480 in the completed MPO scopes. These calls operate on differently sized arrays; they are not hardware floating-point-operation counts. The counters and phase observations do not establish a speedup, an exclusive low-level cause or the overhead of instrumentation. These are verification timings, not production DMRG solver timings. No producer was rerun.

The 15-minute allocation used one CPU and3500MB. The parent finished in425.003seconds. Cleanup sent only SIGINT, took2.788seconds and left no live owned group processes. The child returned2 after saving an explicitly unresolved partial report. The reached time cap cannot be converted to PASS by successful cleanup.

An independent saved-evidence audit passed **243 checks**, with zero failures; a separate before/after runtime-metadata audit passed **29 checks**. The latter reparsed173 mapped segments,33 files and191 loaded modules in each saved report. Both recorded runtime reports are PASS. The audit verified the16-file analysis identity separately from the unchanged70-file producer source, all three original inputs,52 preserved records and the actual release/review bindings. It did not recompute quantum results, execute another runtime or freshly hash/restat remote providers.

`attempt_result.json` records audited outcome and identity hashes; `observed_profile.json` records selected timings, actual calls and partial coverage. Raw maps, inventories, commands, approval bodies, private paths and large arrays are excluded. The certificate remains conditional floating-point error accounting, not directed-rounding interval arithmetic.

All earlier global-pilot and diagnostic timeouts remain preserved. There is no complete comparison PASS or full-production admission. The ledger stays **49 PASS / 0 FAIL / 2 UNRESOLVED / 49 NOT_RUN**, with P17 **NOT_RUN**. A future time budget must account for the unmeasured remaining package states and operators rather than treating this partial profile as a complete runtime estimate.
