# L32-only Bethe DMRG feasibility pilot

Release **bethe_l32_pilot_v3**, prepared and frozen before execution. The earlier L8/L12 and L24 releases, all raw results and both L24 audit versions remain intact. No L32 numerical physics, job submission, L48 calculation or thermodynamic extrapolation was performed to prepare this release. Root owns the live-queue/resource check, staging and any later submission.

This release measures resource and finite-chain convergence feasibility. It cannot complete A18, A19 or A20. The native target retains `scientific_completion_allowed=false`, and the unchanged scientific summarizer branches always report zero distinct scientific passes and PILOT_ONLY/UNRESOLVED for those IDs. The catalogue statuses stay NOT_RUN.

## Why an L32 pilot is now eligible

All six L24 cases in ISAAC job6458547 passed the original quality gates. The portable independent audit passes168 checks and freshly replays the exact frozen summarizer from its source archive; its output matches the original summary byte-for-byte. The audit's two original temporary inputs are now permanent, clearly named sidecars in the L24 raw directory. The first auditor and first report are separately preserved.

The resource continuation criteria were published before the L24 run in the [earlier pilot audit at commit66950c7](https://github.com/yiqiwan0748/DQMC_correlated_electron/blob/66950c791f7f939ca176e5427be62cd8dfca8348/advisor_preparation/distinct_benchmarks/reporting/reports/bethe_pilot_audit_6458253/README.md), dated2026-10-22T19:12:58Z. The exact published bytes and their Git blob identity are preserved, rather than silently using a later README revision.

| Engineering condition | L24 measurement | Result |
|---|---:|---|
| Job elapsed time at most30min |991sec =16m31s|Meets condition|
| Process peak RSS at most2300MiB |1,484,548KiB =1449.75390625MiB|Meets condition|

These are previously declared resource continuation gates, separate from solver tolerances. `engineering_prerequisite_audit.json` is an exact copy of the new portable audit; `engineering_gate_anchor_66950c7.md` is the exact published reference. Their hashes are bound by `pilot.toml` and the source checksum closure. The wrapper checks them before launching either numerical process. This authorizes preparation of a bounded feasibility pilot, not automatic submission or a bulk scientific conclusion.

## Exact fixture and unchanged numerical protocol

Use one32-site open Hubbard chain, hopping1, interaction4, physical chemical potential0, and the actual existing user `hubbard_mpo` with ITensorMPS as the optimizer. The physical Hamiltonian is uncentered. There is no ED-space allocation or Monte Carlo sampler call.

| Sector | Nup | Ndown | Total N | 2Sz |
|---|---:|---:|---:|---:|
| half |16|16|32|0|
| minus |16|15|31|1|
| plus |17|16|33|1|

The two arms are exactly the L24 arms:

| Arm | Bond cap | Sweeps | Cutoff | Variance/site limit | Final truncation-error limit |
|---|---:|---:|---:|---:|---:|
| coarse |256|18|1e-12|1e-6|1e-8|
| refined |512|24|1e-13|1e-8|1e-10|

There are six solves. The generic even-length initializer, seeded shuffle, QN-preserving random MPS, low-bond schedule, noise and local eigensolver tolerance are unchanged. Seed bases841001/841002 and the unchanged formula `base+1000*L+10*sector_index` give half873011/873012, minus873021/873022 and plus873031/873032. These differ from earlier length-specific seeds. Every completed arm immediately writes a new TOML receipt with its full local profiles and truncation history; a slow, failed or missing arm is retained and is never replaced.

`pilot_target.jl` and `oracle_lieb_wu.py` are byte-identical to the L24 and original pilot files. The original user MPO source, Project/Manifest files and runtime versions are pinned to the same hashes. All quality gates, per-arm thresholds, oracle quadrature settings, spin-sector rules and physical conventions are unchanged. The summarizer differs only in one planning sentence describing L32 rather than L24; after normalizing that string, its AST is identical. No numerical acceptance branch changed.

The independent oracle remains the same U4/t1 thermodynamic Lieb–Wu quadrature with the same cutoffs, integration checks and derivative check. Its primary sources and derivations remain in the [original protocol](../bethe_ansatz/README.md). Finite-chain energy/site, charge gap `Eplus+Eminus−2Ehalf` and double occupancy are diagnostics relative to the bulk integrals. They are not required to equal the infinite-chain values. Particle-hole identity `Eplus−Eminus=U`, norm/particle counts, contracted energy, variance, truncation and coarse/refined observable stability retain all prior thresholds.

## Resources and limits

Request one CPU,3500MiB and a60-minute Slurm wall limit. The target command has a55-minute timeout; the independent tiny quadrature and metadata use the remaining margin. The engineering runtime range is **20–50minutes after start**, excluding queue delay. It is empirical and unmeasured at L32. Applying only the fixed-chi length factor4/3 to the actual L24 job time16m31s gives about22m01s before an additional guard. Effective QN blocks, entanglement, retained bonds, convergence, variance contractions, cold compilation and filesystem waits can change the scaling substantially; the range and cap are not guarantees. The L24 memory measurement establishes headroom for this bounded test, not a rigorous L32 memory bound.

`pilot.sbatch` writes a separate `results/bethe_l32_pilot_JOBID` directory and retains `/usr/bin/time -v` receipts. Root supplies `VALIDATION_ROOT` and `BETHE_TRANCHE`; source-root selection remains `source/renyi_disjoint_corrected`. No dependency, nested submission or automatic follow-on run is embedded. Root must keep the user's total running scientific jobs within the existing concurrency limit and state a useful EDT checkpoint after examining the current queue.

A memory limit, timeout, incomplete case or failed quality gate is a measured feasibility limitation. Preserve partial arm receipts and logs, then assess a separately versioned resource/precision proposal without loosening existing gates. No larger pilot is authorized by this release. In particular, **L48 is not included or queued**.

No polynomial fit, thermodynamic error envelope or held-out-size validation is implemented here. The original full proposal's training/held-out sizes and scientific choices remain unchanged and unexecuted. L32 feasibility cannot replace those checks or create a scientific benchmark completion. Any reuse of bridge data in a future full release must be declared before that full analysis.

## Review and exact changes

`protocol_changes.diff` and `protocol_changes.json` record every change from L24: chain length/release/predecessor identity, resource range/caps, frozen engineering prerequisites, one prose sentence and the wrapper's pre-execution metadata checks. All executable scientific target/oracle code stays identical. Static review parses Python AST, Julia syntax, TOML and shell syntax; checks exact predecessor/target/oracle/source identities; verifies all physics/arm/gate/oracle settings remain equal; and derives the six requested sectors without evaluating a numerical physics routine. `SHA256SUMS` closes this release.
