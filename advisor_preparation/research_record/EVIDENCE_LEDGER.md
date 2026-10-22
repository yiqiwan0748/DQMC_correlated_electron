# Evidence ledger

Status date: 22 October 2026. Historical runs are reported from their saved data and audit records. The assembler did not rerun their simulations, reconstruct large ED vectors, or certify every solver. SHA-256 integrity and the listed scalar arithmetic are independently checked by `verify_evidence.py`.

## R1 — Five subsystem entropy benchmarks

- **Quantity and setup:** Second Rényi entropy `−log Tr(ρ_A²)`, natural logarithm. Periodic half-filled 1D Hubbard chain, `L=10,t=1,U=2,μ=0,β=40,Δτ=0.05,n_stab=10,Ninc=16`; four root seeds `[1234,11234,21234,31234]`, 500 burn-in sweeps, 4,000 measurements and 40 bins per node. Contiguous `A=1:nA`.
- **Source:** `source_snapshot/HubbardRenyiDQMC/results/processed/BLOCKSIZE_U2_beta40_Ninc16_summary.csv`; detailed method and checkpoint precision in `docs/l10_u2_beta40_blocksize_benchmark.md` within that snapshot. Authoritative table commit `88c21ae3345ad7414d86b37e777e10c35e6451f2`; original archived-worktree table SHA-256 `1f799d83d74b6ad31f21122f5bb5a40935846e449603096aa2e08f8e987fc3db`, preserved separately at `evidence/HubbardRenyiDQMC_archive_worktree/BLOCKSIZE_U2_beta40_Ninc16_summary.csv`. The canonical Git-blob copy under `source_snapshot/` has SHA-256 `80109e5dd636c2d4d296b6eda159757e99e2927f84bf45fe37e2b26fba6600be`: Git normalized one CRLF line ending. The six CSV rows and numerical values are identical.

| nA | S₂ DQMC | MC standard error | Ground-state ED | Recomputed difference / MC error |
|---:|---:|---:|---:|---:|
| 1 | 1.308138170759 | 0.001209659236 | 1.308461130863 | −0.267 |
| 2 | 1.274574334430 | 0.002188963049 | 1.273890131213 | +0.313 |
| 3 | 1.634894520051 | 0.002334522219 | 1.633561808649 | +0.571 |
| 4 | 1.524569914598 | 0.003041891872 | 1.524217285675 | +0.116 |
| 5 | 1.716441955736 | 0.003040377094 | 1.716551661433 | −0.036 |

- **Independent arithmetic:** Maximum `abs((DQMC−ED)/error)=0.5708711578`. The half-chain difference is `−0.000109705697`.
- **Limit:** ED is a ground state in `N↑=N↓=5`; DQMC is grand canonical at finite β. These are different ensembles. Quoted MC errors do not bound all thermal/Trotter bias. Different ED serializations differ by about `10⁻¹¹`, negligible on this MC scale.
- **Reproduce:** Use the exact source snapshot and parameter files, a fresh output directory, the node/grid runner and `analysis/summarize_grid.jl`; see `REPRODUCE.md`. The safe fast integrity check uses the preserved table and no simulation.

## R2 — Finite-size data through L=22

- **Setup:** Same periodic half-filled `U=2` chain, half-chain subsystem, `β=4L,Δτ=0.05,Ninc=16`; four seeds, 500 burn-in, 4,000 measurements and 40 bins per node. Each L20/L22 production grid contains 68 node/root pairs.
- **Data:** `results/processed/FSS_1D_U2_beta4L_L20_Ninc16.csv`, `FSS_1D_U2_beta4L_L22_Ninc16.csv`, and `DMRG_S2_U2_periodic.csv` in the source snapshot. Production code commit `c4c873ffa16908a6641f851cb6c9ad41a9a22c29`; promoted summaries commit `c766667ee7ee327bbf0744bfacbe01b82fd4cea0`.

| L | β | S₂ DQMC | Nominal MC error | DMRG | Difference / DQMC error |
|---:|---:|---:|---:|---:|---:|
| 20 | 80 | 1.805179361272 | 0.008209270540 | 1.808766655046 | −0.437 |
| 22 | 88 | 2.025158749526 | 0.004730458606 | 2.021914792180 | +0.686 |

- **Diagnostics:** Targeted 8,000-burn-in checks covered only L20 nodes `6,11,13` and L22 node `14` (16 node/root jobs). L20's conservative standalone error is `0.0095`; its canonical fit retains the nominal error above. L22 blocking was stable.
- **Limit:** DMRG references have no quoted uncertainty and reach a max link dimension of about 4,800. Original large-L production and long-burn raw corpora remain at documented external Windows paths; their manifests and promoted summaries are preserved. Their full raw-file hashes were not rechecked in this assembly. L20/L22 lack independent β and Δτ convergence tests.
- **Provenance/method:** `docs/one_d_hubbard_u2_fss_L4_to_L22_20260820.md` in the snapshot records all summary, source and raw-manifest hashes plus exact reproduction commands.

## R3 — Finite-size interpretation

- **Evidence:** The canonical fit input has 22 observations: 10 DQMC, 5 ED and 7 DMRG. `oneD_hubbard_fss_fit_summary_L4_to_L22_20260820.csv` has 289 fit rows.
- **Finding:** The five principal DQMC models have χ²/dof approximately `477.89,29.28,3.90,80.56,244.18`. Their fitted slopes are finite-window diagnostics, not accepted asymptotic central charges. Preserve the mod-4 branch structure and report model inadequacy.
- **Reproduce:** `analysis/build_oneD_hubbard_fss_L4_to_L22_input.jl` then `analysis/fit_oneD_hubbard_fss_with_oscillation.jl`, using fresh output paths and the snapshot's FSS environment. Exact commands are in the detailed FSS note.

## R4 — Documented Grover update optimization

- **Preserved identity:** `g=G₁G₂+(I−G₁)(I−G₂)` and `det(g+p rᵀ)/det(g)=1+rᵀg⁻¹p`, with safeguards and explicit fallbacks. This reduces the Grover determinant-ratio calculation from cubic to quadratic in subsystem size.
- **Historical validation:** `docs/grover_rank1_speedup.md` and `docs/dqmc_renyi_ee.md` report 2,318 checks, 12 same-seed direct/optimized rounds, and 36 matching historical bin files after omitting the run-name field.
- **Historical timing:** Isolated L8 kernel: `0.550769 → 0.204545 s` (documented 2.692658×). Matched L10 interior-node means: `95.9854 → 25.1035 s` (documented 3.8236×, 28 nodes). Exact timings are workload/hardware dependent.
- **Verification level:** Mathematical identity and code are inspectable. The 2026-10-22 record builder has not independently retimed the code or recovered the original full performance corpus. Treat the quoted speedups as documented historical results.
- **CV discrepancy:** No source for a 4.59× equal-/unequal-time Green-function speedup was found in the recovered Rényi tree or June source. That may be a different optimization. Its implementation, workload, baseline, timing and dates need separate evidence; do not substitute the Grover timing.

## R5 — Electronic DMRG against an exact finite-cluster reference

- **Setup:** Open 2×2 two-orbital Hubbard cluster, eight electronic sites, `U=20,N=6,t=(-1,-1.75,-0.85,-0.65)`, centered interaction. Physical Sz=0,1 and seeds 314159,271828 give four arms. Bond caps 32,64,128,256, six sweeps per cap and six extra zero-noise sweeps at 256; cutoff `10⁻¹²`. ITensors 0.9.30 and ITensorMPS 0.4.1.
- **Files:** `evidence/AltermagnetDQMCTrial/results/literature_transfer/alm_dmrg_benchmark_20260911/attempt01/` contains returned CSVs, metadata, source and summary. `.../alm_method_benchmark_20260911/dmrg_return_audit.toml` is the independent returned-CSV audit.
- **Result:** 4/4 final arms passed. Maximum final energy error `3.197442310920451e−14`, full spin-dot matrix entry error `7.346068198188505e−13`, and squared-spin error `8.881784197001252e−16`; minimum ground-subspace fidelity reported `0.9999999999999989`.
- **Meaning:** Agreement of electronic MPO/DMRG and a small ED reference, including explicit fermion-sign fixtures and contact terms. These errors describe this finite problem, with an exactly representable final bond cap.
- **Limit:** The independent audit recomputes CSV arithmetic and compares to the original ED reference; it does not reconstruct MPS or independently rerun MPO actions. The 109 MB serialized exact reference and original MPS states remain in the source repository and are not copied here. The reference SHA-256 is `ae3ca51533fe8f91b46cd9fbb8f65e06dcd91b4a160739c3020487fb58eb1c09`. No bulk order claim follows.

## R6 — Verified selected eight-cell ED singlet

- **Setup:** Tilted periodic eight-cell two-orbital cluster, vectors `(2,2),(-2,2)`, 16 orbitals, `N↑=N↓=6,U=20,t=(-1,-1.75,-0.85,-0.65)`, total K=X. Energy uses `U(n↑−1/2)(n↓−1/2)`, so uncentered energy is centered energy plus 40.
- **File:** `evidence/AltermagnetDQMCTrial/results/isaac_20261014/moreo_configurations/solve003_terminal_validation_20261015.json` records an accepted ISAAC job, sector dimension 8,015,616, source manifest and fresh native original-H/full-S² actions.
- **Result:** `E=−65.59824971805813`, `||Hψ−Eψ||=1.1365542933084897e−10`, full-S² residual `3.987806681587894e−14`.
- **Limit:** This is an excited singlet: an earlier X state with S=3 lies lower. Sector minimum and global ground state are explicitly uncertified. The receipt audits small records of fresh native operations; large remote vectors were not revalidated locally. No altermagnetic phase is proved.

## Attribution and timeline to confirm

The Rényi Project.toml names Yingqing Wang; Git metadata alone cannot assign intellectual or implementation credit. Confirm which algorithm, solver, optimization, data campaign and analysis portions were personally implemented, which relied on collaborators or AI assistance, and how JDQMCFramework, JDQMCMeasurements, ITensors and ITensorMPS should be credited. The recovered August 2026 improvements must not be silently attributed to an appointment ending in July 2026. The numerical record can retain later dated work as a separate continuation.
