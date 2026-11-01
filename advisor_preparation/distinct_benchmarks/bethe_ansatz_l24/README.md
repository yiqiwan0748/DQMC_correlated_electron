# L24-only Bethe DMRG resource/quality bridge

Release: **bethe_l24_pilot_v2 — FROZEN, NOT RUN**. This is a new pilot release; the original `bethe_ansatz/` release and job 6458253 receipts remain immutable. No local physics, job submission, L32 run, full extrapolation, or source-kernel change was performed in preparing it. Root owns staging, the live-queue check, and submission.

This bridge measures whether the proposed larger bond cap is feasible before any thermodynamic programme. It cannot complete A18, A19, or A20. The target reports `scientific_completion_allowed=false`; the summarizer always reports zero distinct scientific passes and UNRESOLVED for those three IDs, even if its separate pilot-quality field passes. The canonical catalogue statuses remain NOT_RUN.

## Exact bounded fixture

One open chain of 24 sites, hopping 1, interaction 4, physical mu 0. It uses the same uncentered physical Hubbard Hamiltonian and actual user MPO as the audited L8/L12 pilot. All three electron sectors are retained:

| Sector | Nup | Ndown | Total N | 2Sz |
|---|---:|---:|---:|---:|
| half |12|12|24|0|
| minus |12|11|23|1|
| plus |13|12|25|1|

The existing generic initializer changes one Dn product-state site to Emp for removal or UpDn for addition, then performs the same seeded shuffle and QN-conserving random-MPS construction. No fixed-particle-sector ED space is allocated. The user source's `hubbard_mpo(...; boundary="obc")` and ITensorMPS supply the actual Hamiltonian and optimizer.

| Arm | Maximum bond dimension | Sweeps | Cutoff | Variance/site limit | Final truncation-error limit |
|---|---:|---:|---:|---:|---:|
| coarse |256|18|1e-12|1e-6|1e-8|
| refined |512|24|1e-13|1e-8|1e-10|

There are **six solves**. Noise, the low-bond ramp, local eigensolver tolerance, initialization algorithm and seed formula are unchanged. Seed bases remain 841001/841002; the `1000*L+10*sector_index` offsets give disjoint L24 seeds 865011/865012, 865021/865022 and 865031/865032. Each completed solve is immediately saved to a new per-arm TOML receipt; no failed or slow arm may be dropped or replaced. Aggregate and resource receipts are written into a new job-specific directory.

## What changed and what stayed identical

`protocol_changes.diff` records the exact manifest, wrapper and one-line summarizer changes. `protocol_changes.json` lists every changed manifest field.

- The only chain length is now 24, replacing pilot sizes 8/12.
- Coarse settings change from chi128/12sweeps/cutoff1e-11 to chi256/18sweeps/cutoff1e-12. Refined settings change from chi256/18sweeps to chi512/24sweeps; its cutoff remains1e-13.
- The Slurm job name and result-directory prefix identify this separate L24 release. The wall cap changes from 15 to 45 minutes, and the target command cap from 780 to 2400 seconds. Resource metadata records the preliminary runtime range and predecessor identity.
- The Python summarizer changes **only** its inherited `planning_limit` sentence so it correctly describes an L24-only pilot. Its AST is otherwise identical after replacing that string; all arithmetic and acceptance branches remain unchanged.

`pilot_target.jl` and `oracle_lieb_wu.py` are **byte-identical** copies of the frozen v1 files. The existing Julia driver reads lengths, arms and sweep counts from the manifest and only requires even L. The existing summarizer forms its expected solve set from the manifest, so `[24] × three sectors × two arms` requires all six cases without a hardcoded case count. Its numerical checks are unchanged.

Every scalar quality tolerance, every per-arm variance/truncation threshold, all seed bases, all oracle quadrature parameters and versions, all physical conventions, and all target/source/environment hashes are identical to v1. The exact user MPO source SHA is `ff1a0becd7305420afa097291ea033c4ba1b12222eae078f9eed87d0c5c2838f`. The remote root remains `source/renyi_disjoint_corrected`; the adapter does not call the ED or Monte Carlo samplers. ITensorMPS remains responsible for the optimizer; this adapter supplies initialization, observables and receipts.

The independent Lieb–Wu quadrature is rerun under the new manifest identity with the same fixed U4/t1 formulas and numerical gates. Its primary reference and full derivations remain in the preserved [v1 protocol](../bethe_ansatz/README.md). Oracle integration, fermion-sector conventions and finite-size offsets have not been redefined.

## Unchanged pilot acceptance and limits

The checks require correct Nup/Ndown and normalization, reported versus contracted energy agreement, the original variance and truncation limits, particle-hole identity `Eplus-Eminus=U`, and the original coarse/refined stability margins for energy/site, mean double occupancy and `Eplus+Eminus-2*Ehalf`. Local profiles and complete sweep records remain required. All three sectors and both arms must be present. A source/version mismatch, timeout, missing case, or failed quality gate remains unresolved. No tolerance was widened.

The source does not infer thermodynamic equality from these finite-chain values. No bulk fit, model selection, uncertainty envelope or held-out-size check is performed. The original full scientific proposal, including its training sizes and held-out 56/64 checks, remains unchanged and unexecuted. Bridge data cannot justify changing its scientific choices after seeing results. L32 is neither included nor automatically queued.

## Resources and parent handoff

`pilot.sbatch` requests one CPU, 3500M and 45 minutes; the scientific target has a 40-minute timeout. Root supplies `VALIDATION_ROOT` and `BETHE_TRANCHE`, pointing the latter to this new release. Outputs use `results/bethe_l24_pilot_${SLURM_JOB_ID}`. The existing separate QuSpin environment and pinned Julia 1.12.1 runtime are used; no environment install or nested submission is included.

The preliminary engineering range is 15–40 minutes **after start**, not a measured L24 runtime. The earlier job took 6m03s with substantial cold/first-use overhead, while its L12 refined solves took about 19 seconds each at cap 256. Doubling both length and bond cap can increase work sharply, and chi 512 memory has not been measured. The hard 3500M limit remains authoritative. Scheduler wait is separate. Root should state an EDT check point before submission and retain partial receipts if a hard limit is reached.

Before a larger bridge or full programme, inspect maximum RSS, actual retained bond dimensions, all six timings, variance/truncation diagnostics, symmetry checks and arm differences. A full thermodynamic execution requires a separate implementation, source release and controlled extrapolation review. No success of this resource pilot can establish a scientific catalogue pass by itself.
