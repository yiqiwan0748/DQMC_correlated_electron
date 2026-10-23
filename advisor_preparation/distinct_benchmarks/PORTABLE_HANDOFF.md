# Portable scientific handoff: 23 October 2026

The current [immutable scientific ledger](reporting/reports/pc_reconciliation_v12/README.md) records **42 PASS, 2 UNRESOLVED and 56 NOT_RUN** across 100 distinct ideas. SHA256: `976272652f95647c133eed9c6261edeedfcf0833b38cf8bc44095398474d6338`. It preserves all 95 unrelated rows and all 39 prior PASS rows from the [previous checkpoint](reporting/reports/ground_initial_v10/README.md). The user explicitly resumed work on October 23; automatic follow-up remains paused.

Clone draft PR 9's branch and record the actual checked-out commit:

```bash
git clone --branch codex/research-record-validation-20261022 https://github.com/yiqiwan0748/DQMC_correlated_electron.git
cd DQMC_correlated_electron
git rev-parse HEAD
```

The [scientific continuation packets](pc_continuation_20261023/README.md) provide compact exact text data, complete source/member hashes, original independent mathematical auditors, frozen comparator programs and portable replay wrappers. They do not require private handoff ZIPs, controller policies, credentials, email, contacts or CV data. Extraction rejects links, path traversal, duplicates and altered bytes. Use your own authorized SSH setup if inspecting ISAAC; keys and host control sockets are excluded.

| Question | Reconciled scientific state | Principal original jobs |
| --- | --- | --- |
| X09 | PASS at unchanged gates with frozen sixteen-sweep precision; original eight-sweep failed gates retained |6460656; original 6458708 |
| X13 | PASS across complete finite charged thresholds/operator visibility; new adapter through unchanged native Hamiltonian |6462017 |
| X16 | PASS for all finite binding intervals under declared conditional binary64 budget; pilot excluded |6462471; pilot 6462342 |
| P21 | UNRESOLVED: two individual-chain ESS failures and 376 complex precision failures across 704 complex entries |6461010; proof 6460659 and pilot 6460727 excluded from credit |
| X01 | UNRESOLVED at selected stage 2; quartic entropy halfwidth misses fixed precision; final stage 3 incomplete |6460771; stage 2 parts 6460990/91/92/94 and analysis 6460996; stage 3 parts 6462170–73 and analysis 6462174 |

P21 uses actual native ALF continuous-Gaussian fields versus actual maintained user Hirsch buffers, matched finite-step V_K order. Its frozen production paired-chain/complex-disk protocol explicitly supersedes the older unexecuted Welch/pilot candidate. Its 704 complex entries are checks within one question, not 704 ideas. X01's complete first/second stages, failed original analysis, analysis-only repair, incomplete final native outputs and cancellations remain visible; no pooling, fourth look or replacement is applied. The repeated X12 result adds no second completion.

Use Python 3.12 with NumPy 2.4.6 and SciPy 1.17.1 for the [four read-only packet replays](pc_continuation_20261023/README.md). They reconstruct full numerical gates without native solvers, Monte Carlo, ALF execution or submission. A successful integrity audit leaves P21/X01 scientifically UNRESOLVED. Setup-only PySCF job 6465769 contributes zero scientific completions. Do not infer a current live queue or new source release from this archival snapshot; inspect fresh scheduler receipts and freeze any newly authorized scientific run separately.

The October 22 deadline was fulfilled. Original jobs 6458825, 6458826 and 6458828 were cancelled before starting, while L32 pilot 6458711 completed and independently passed 177 pilot checks with zero catalogue credit. The [initial stop receipt](submissions/user_stop_deadline_20261022.json), [final stop receipt](submissions/user_stop_final_state_20261022.json), [scheduler accounting](submissions/user_stop_final_accounting_20261022.txt) and [L32 audit](reporting/reports/bethe_l32_audit_6458711_portable_v1/README.md) remain immutable. [handoff_state.json](handoff_state.json) retains that historical stop snapshot separately from today's audited ledger.
