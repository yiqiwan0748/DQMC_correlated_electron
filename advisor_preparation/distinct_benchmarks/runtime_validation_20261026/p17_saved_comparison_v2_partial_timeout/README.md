# Timeout preserves partial comparison evidence

The v2 saved-input diagnostic reached the same 420-second run cap and preserved the unchanged comparator's partial report during bounded cleanup. Its parent remains **UNRESOLVED**, with child status **TIME_CAP_EXCEEDED** and Slurm **FAILED / 2:0**. The child report's raw `RUNNING` label means **partial data saved on interruption**; it does not describe an active or passed job.

The actual cleanup sent only SIGINT, completed in `1.372` seconds and left no live child-group processes. Child return code was `130`; measured child time including cleanup was `421.393` seconds and parent time was `425.695` seconds. The declared cleanup policy allowed SIGINT30/TERM3/KILL2 within a 35-second cleanup reserve, with a separate 10-second finalization reserve. Reaching the run cap irrevocably leaves the parent unresolved, even if cleanup can save more evidence.

The partial report contains **829 recorded PASS gates and zero recorded FAIL gates**:

| Recorded gate prefix | Count |
| --- | ---: |
| QuSpin reference | 277 |
| Zero-field spectral response | 6 |
| Native | 513 |
| TeNPy | 33 |

Its completed state list contains **six native cases** across all three pilot fields and both starts, and **zero completed TeNPy states**. All 21 reference-sector certificate records and six zero-field spectral responses were saved. The last TeNPy gates cover local operators and represented MPO Hamiltonians in sectors0,1,2. Planned later states and cross-package profile/response gates must not be inferred from this partial census. Actual arithmetic-call counts were not saved; the unchanged source retains its five-million-call ceiling.

The saved stderr records interrupt delivery during `compare.py`'s `contract_mpo` at line149, through multiplication in `certificate.py` at lines247,165,81 and then NumPy reduction code. This locates the reported interruption; it does not prove that this stack consumed all or most elapsed time. `partial_evidence.json` includes only stack basenames, functions, line numbers and source hashes, without private paths or raw traceback text.

Both actual before/after runtime reports record PASS. The independent audit reparsed 173 mapped segments and 189 loaded modules in each against frozen source/provider pins. It also verified source-bound traceback snippets, exact saved inputs, signal cleanup, the new seven-file diagnostic, all 70 upstream source files, three original producers, 52 preserved records and both approvals. Its **239 saved-evidence checks** passed with zero failed checks. This audit did not recompute quantum results or certificates, hash/restat remote providers anew or execute another runtime.

The comparator, certificate, numerical gates and producer files remain unchanged. No producer was rerun. `attempt_result.json` records the outcome and identity hashes; `partial_evidence.json` records the observed partial census and sanitized stack. Raw maps, installed-file inventories, commands, paths, approval bodies, raw traceback and large physical arrays are excluded.

Recorded partial passes do not establish a completed physical comparison. The original V8 and saved-comparison V1 failures remain preserved; this attempt grants no whole-pilot pass, full-production admission or catalogue completion. The ledger remains **49 PASS / 0 FAIL / 2 UNRESOLVED / 49 NOT_RUN**, with P17 **NOT_RUN**.
