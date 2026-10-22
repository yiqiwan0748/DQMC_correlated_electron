# Portable handoff

**Stop request: 22 October 2026, 17:00 EDT (21:00 UTC), regardless of outcomes.** Automatic follow-up is paused. Jobs 6458825, 6458826 and 6458828 were cancelled before starting at 16:22:12 EDT. L32 pilot 6458711 completed at 16:22:48 EDT after 21m56s, exit 0 and MaxRSS 1,735,276 KiB. Its pilot quality PASS has been [independently audited](reporting/reports/bethe_l32_audit_6458711_portable_v1/README.md) with 177 checks and adds no catalogue completion. All campaign jobs are stopped before the deadline. The [initial stop receipt](submissions/user_stop_deadline_20261022.json) preserves the original submissions and administrative override; the [final stop receipt](submissions/user_stop_final_state_20261022.json) and [scheduler accounting](submissions/user_stop_final_accounting_20261022.txt) confirm no campaign job remains pending or running. No scientific gate changed. **No automatic resume or new scientific submission is authorized.**

## Open the same record on another computer

```bash
git clone --branch codex/research-record-validation-20261022 https://github.com/yiqiwan0748/DQMC_correlated_electron.git
cd DQMC_correlated_electron
git rev-parse HEAD
```

Record the actual cloned commit. This is [draft PR 9](https://github.com/yiqiwan0748/DQMC_correlated_electron/pull/9). Use the new computer's own authorized ISAAC SSH setup; computer-specific keys and control sockets are not part of the handoff. Scientific submission receipts identify the existing campaign paths, so inspect those before choosing a scratch/source directory.

## Scientific state and original evidence

The immutable [ground_initial_v10 ledger](reporting/reports/ground_initial_v10/ledger_100.json) records exactly 100 distinct IDs: **39 PASS, 1 UNRESOLVED, 60 NOT_RUN**. SHA256: `abee1b7af4263979880d24afc679d1553567e2141caf16c41888c7c46cf2e228`. All 38 prior PASS rows and 98 other row objects are unchanged. The initial 100 parameter comparisons/188 solver arms are a separate count.

X12 passes its interacting canonical-dimer prescribed-state comparison. X09 remains unresolved at the new certification adapter's eight-sweep budget versus the production default of 14. Full Hamiltonian/MPO checks agree with QuSpin; three starts miss optimization-precision gates. [Original job 6458708 raw states and audits](results/ground_certification_6458708/) retain 15 failed detailed gates and unchanged thresholds. Do not overwrite that attempt with a later retest.

| Stage | Current verified status | Frozen source and receipts |
|---|---|---|
|L32 Bethe/DMRG pilot 6458711|Completed 16:22:48 EDT; independently audited pilot quality PASS; zero catalogue completions|[source](source_releases/bethe-l32-pilot-v3-freeze.json), [submission](submissions/bethe_l32_pilot_6458711.json)|
|Precision retest 6458825|Cancelled before start; NOT_RUN_USER_STOP|[revision 2](source_releases/ground-certification-precision-v2-freeze.json), [submission](submissions/ground_certification_precision_v2_6458825.json)|
|ALF eligibility 6458826|Cancelled before start; P21 still NOT_RUN|[revision 5](source_releases/alf-measurement-proof-v5-submission-freeze.json), [submission](submissions/alf_measurement_proof_v5_6458826.json)|
|X01 timing pilot 6458828|Cancelled before start; zero completions|[pilot](source_releases/continuum-bias-timing-pilot-v1-freeze.json), [submission](submissions/continuum_bias_timing_pilot_6458828.json)|

Cancelled unrun stages add no scientific pass or failure. The original precision, ALF and timing protocols remain available. ALF revision 4 compiled/linked and passed its provider check, then failed on a missing startup seed; [those text receipts](results/alf_adapter_proof_v4_6458709/failure_audit.json) remain intact. Revision 5 only supplies that input. Successful eligibility would still not complete P21 sampling.

## Read-only reconciliation before any explicit resume

```bash
sacct -j6458711,6458825,6458826,6458828 --format=JobID,State,ExitCode,Elapsed,MaxRSS,AllocCPUS,ReqMem -P
squeue -j6458711,6458825,6458826,6458828
```

[handoff_state.json](handoff_state.json) holds full source-archive hashes and the stop policy. The L32 audit verifies all 32 bound inputs, six sector/precision solves and byte-identical frozen summary replay. Its finite-size pilot is not a thermodynamic benchmark or a catalogue completion. Preserve every cancelled, failed or partial attempt. After a new explicit human instruction, choose a new recorded output directory, verify existing source/runtime hashes and freeze changed launch paths. Existing extraction/output guards deliberately prevent overwrites. Fresh resource admission includes unrelated jobs and the previous one CPU, 3500 MB/two-total-running-job policy.

Start with the [ground_initial_v10 replay guide](reporting/GROUND_INITIAL_V10_REPRODUCE.md). It uses Python 3.11+ standard-library postprocessing and the five published user source/environment snapshots, with no solver execution or Julia installation. Earlier reports and failures are retained. Solver reconstruction uses [ISAAC_REPRODUCE](../isaac/ISAAC_REPRODUCE.md) and the linked frozen source releases. The existing runtime is Julia 1.12.1, QuSpin 1.0.1/NumPy 2.4.6/SciPy 1.17.1 and ITensors 0.9.30/ITensorMPS 0.4.1; verify the frozen Project/Manifest rather than upgrading libraries. No CV, contacts or private continuation files are included.
