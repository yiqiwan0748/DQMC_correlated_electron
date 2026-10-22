# Hubbard-model research record

Numerical evidence for advisor conversations, assembled 22 October 2026. This record preserves the calculation settings, source code, benchmark tables and limits behind each claim. Historical numerical runs are distinguished from checks performed while assembling this record.

## Rényi entropy of the 1D Hubbard chain

The Julia workflow evaluates `S₂ = −log Tr(ρ_A²)` with a two-replica, adjacent-ratio λ-grid DQMC estimator. At half filling, `t=1`, `U=2`, periodic `L=10`, `β=40`, `Δτ=0.05`, and 16 λ intervals:

- The half-chain result is **S₂ = 1.716442 ± 0.003040**, compared with ground-state ED **1.716552**.
- All five contiguous subsystem sizes `nA=1…5` lie within **0.6 Monte Carlo standard errors** of ED.
- The finite-size DQMC study reaches **L=22**. At `β=88`, its half-chain result is **2.025159 ± 0.004730**, compared with DMRG **2.021915**.
- Strong `L mod 4` corrections and poor finite-window fits prevent an asymptotic central-charge claim.

DQMC uses a finite-temperature grand-canonical ensemble; ED and DMRG use fixed-particle ground states. The DMRG table reports no uncertainty. Agreement is a cross-method consistency check, with those limits retained.

The optimized Grover determinant update has a documented **3.8236×** speedup on one matched L=10 workload, with unchanged historical bin values. This August 2026 timing differs from the **4.59× Green-function speedup** in the supplied CV, for which a source has not yet been identified. The CV claim should be supported separately.

[Benchmark table](source_snapshot/HubbardRenyiDQMC/results/processed/BLOCKSIZE_U2_beta40_Ninc16_summary.csv) · [Scientific account](source_snapshot/HubbardRenyiDQMC/docs/dqmc_renyi_ee.md) · [Reproduction](REPRODUCE.md) · [Evidence ledger](EVIDENCE_LEDGER.md)

## Related ED and electronic DMRG validation

A separate two-orbital Hubbard benchmark compared electronic DMRG against ED on an open 2×2 cluster at `U=20`, total `N=6`. Four independent sector/seed arms passed the registered energy, variance, full spin-correlation and particle-number checks. The saved maximum final energy error is **3.20×10⁻¹⁴**; this measures numerical agreement for that small reference problem.

The larger tilted eight-cell ED study verified a selected `K=X`, `S=0` eigenstate at `U=20`, with energy **−65.5982497181** and original-H residual **1.14×10⁻¹⁰**. A known lower `S=3` state makes this an excited singlet. Neither a certified ground state nor a bulk altermagnetic phase is claimed.

## Source and contribution record

The complete Rényi source was recovered from the 28 August backup, at commit `b0adedf042304a72bc1aea618abdcda1962242e9`; its exact committed tree is preserved under `source_snapshot/`. The separate backup bundle advertises a different HEAD. The old Downloads copy is an incomplete scaffold and is unsuitable as the numerical source.

The project declares Yingqing Wang as author. It depends on JDQMCFramework/JDQMCMeasurements and ITensors/ITensorMPS. Personal contributions, collaborator roles and work dates need the author's confirmation before using first-person achievement claims. The ledger supports measured outcomes; it does not infer sole authorship from a repository.

`provenance.json` records source paths, commits and SHA-256 hashes. `python3 verify_evidence.py` verifies the preserved files and recalculates the small arithmetic checks without starting a scientific simulation.
