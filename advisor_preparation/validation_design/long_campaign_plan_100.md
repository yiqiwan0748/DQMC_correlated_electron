# Supporting internal regression and extended numerical campaign design

This supporting design is excluded from the canonical 100 external scientific benchmark comparisons. It contains exactly 25 DQMC, 25 ED, 25 DMRG and 25 Rényi/workflow cases. It is a design, not an execution report. Every row starts as `NOT_RUN`. One long Julia suite can cover multiple cases, but a passing suite does not establish an unimplemented quantitative case. Do not report assertion counts as independent physics tests.

The complete August 28 HubbardRenyiDQMC archive was inspected at `/tmp/advisor_renyi_evidence_20261022/HubbardRenyiDQMC`. Its source snapshot is being preserved separately. The older Downloads package has an empty module and only scaffold tests and is not the source for this campaign. The ALM implementation was inspected at `/home/yiqiwan/code/AltermagnetDQMCTrial`; its isolated ITensor benchmark has its own pinned environment.

## What can run first

Run the ALM suite, the complete Rényi suite and the focused Grover validation in ISAAC batch allocations, then audit the registered DMRG preflight/run. Keep the source commit or archive hash, Project/Manifest hash, Julia version, BLAS/thread counts, Slurm job ID, command, resource usage and complete stdout/stderr. Write all new outputs to unique campaign directories. Preserve accepted raw data, archived estimator source and prior snapshots. A scheduler exit code zero is not a scientific pass by itself.

`existing_test` means a corresponding source assertion/fixture exists; it does not mean that this campaign executed it. `existing_capability` means the code supports the calculation but the full specification still needs an adapter or additional evidence. `planned_adapter` means implementation remains. The exact maps and numeric tolerances are in JSON/CSV. The initial inventory is: {'existing_test': 37, 'planned_adapter': 45, 'existing_capability': 18}.

## Ensembles and conventions

Use Hc = sum(K_ij c†_i c_j) + U sum_i(nup_i-1/2)(ndown_i-1/2) - mu_c N for centered DQMC. Chain ED and DMRG references use Hu = K + U sum_i D_i - mu_u N in a specified fixed sector. The grand-canonical conversion is mu_u=mu_c+U/2, Hc=Hu+UL/4. At the same numerical mu in a fixed N sector, Eu-Ec=UN/2-UL/4. Entropies of a given pure state ignore an overall energy constant, but thermal ensemble weights do not ignore a wrong chemical potential.

Finite-beta DQMC targets a grand-canonical thermal density matrix and a finite imaginary-time discretization. Fixed-N ground-state ED and ordinary DMRG target pure states. They can be compared only after matching Hamiltonian, geometry, boundary, sector/manifold, subsystem, normalization and energy reference, and separately bounding thermal and Trotter bias. The source has a useful L4,U4,beta16,dtau0.05,nA2 example: ground S2=0.7521220360017, exact thermal S2=0.7712028368091 and asymmetric finite-dtau S2=0.7689951876109. Recalculate these in this campaign before presenting them as newly verified results. A pure-state complement identity must not be applied to a thermal mixture.

All entropies use natural logarithms. For normalized Schmidt probabilities p, S2=-log(sum(p^2)); this differs from von Neumann entropy -sum(p log p). For a number-conserving free Gaussian state with restricted correlation eigenvalues nu, S2=-sum log(nu^2+(1-nu)^2), with both spin species included. For a thermal full system, purity=Z(2beta)/Z(beta)^2 for the same Hamiltonian. Explicit fermionic parity must accompany regrouping all-up/all-down coefficient bases into local site order.

The ALM DMRG reference is a fully open 2x2 physical-cell system with eight spinful orbital modes, N6, U20, hoppings (-1,-1.75,-0.85,-0.65), physical Sz0 and Sz1 sectors, and centered interaction. It is a different model and geometry from the one-dimensional Rényi chain. Keep their evidence distinct.

## Statistical acceptance

The JSON records a familywise error rate of 1%. Freeze the full number C of scalar comparisons across all stochastic cases and parameter arms before running them. Use Bonferroni alpha=0.01/C and a Student-t or validated simultaneous block-bootstrap interval. For an accuracy claim, require the entire confidence interval for estimate minus oracle to lie within the stated equivalence margin. A mere non-significant discrepancy is insufficient evidence of accuracy.

Use independent validation seeds, retain paired numerator/denominator blocks, measure autocorrelation, and require sufficient effective samples and stable blocking. Preserve within-node forward/backward covariance. Do not take log(abs(ratio)) to conceal an invalid sign, silently drop nonfinite bins, or count a noisy denominator as a valid estimator. Insufficient precision, failed convergence, unavailable dependency, missing output or exhausted budget is `UNRESOLVED`; disagreement beyond a well-resolved margin is `FAIL`. Report every case, including failed and pending cases.

Initial T0/T1/T2/T3 budgets are unmeasured planning ranges in JSON, not runtime forecasts. Pilot timings must precede a production allocation estimate. Separate queue delay, compilation time, solver time and convergence uncertainty. The catalogue builder itself performs no numerical simulation.

## Case list

| ID | Test | Source status | Group | Tier |
|---|---|---|---|---|
| DQMC-01 | Parameter and imaginary-time-grid rejection | existing_test | alm_suite | T0 |
| DQMC-02 | Four-state Hubbard-Stratonovich identity | existing_test | alm_suite | T0 |
| DQMC-03 | Kinetic Fourier assembly and periodic bond multiplicity | existing_test | alm_suite | T0 |
| DQMC-04 | Free equal-time Green matrix | existing_test | alm_suite | T0 |
| DQMC-05 | Free grand-canonical thermodynamic trace | existing_test | alm_suite | T0 |
| DQMC-06 | Fixed-field cyclic propagator at every time slice | existing_test | alm_suite | T0 |
| DQMC-07 | Local HS flip determinant ratio | existing_test | alm_suite | T0 |
| DQMC-08 | Accepted rank-one Green update | existing_test | alm_suite | T0 |
| DQMC-09 | Rejected update leaves physical state unchanged | planned_adapter | new_alm_adapter | T0 |
| DQMC-10 | Full sweep direction and stabilization consistency | existing_test | alm_suite | T0 |
| DQMC-11 | Stabilization-period independence with fixed fields | planned_adapter | new_alm_adapter | T1 |
| DQMC-12 | Centered-interaction and chemical-potential energy identities | existing_test | alm_suite | T0 |
| DQMC-13 | Onsite spin contact identity | existing_capability | alm_suite | T0 |
| DQMC-14 | Exactly-once sign reweighting and covariance | existing_test | alm_suite | T0 |
| DQMC-15 | Zero or unresolved sign denominator is withheld | existing_test | alm_suite | T0 |
| DQMC-16 | Free unequal-time Green function | existing_test | alm_suite | T0 |
| DQMC-17 | Fermionic endpoint discontinuity | existing_test | alm_suite | T0 |
| DQMC-18 | Finite spatial Fourier normalization and orbital rotation | existing_test | alm_suite | T0 |
| DQMC-19 | Single-slice Monte Carlo versus exact HS enumeration | existing_capability | alm_enum_single | T2 |
| DQMC-20 | Two-slice Monte Carlo versus exact HS enumeration | existing_capability | alm_enum_multi | T2 |
| DQMC-21 | Exhaustive local detailed balance | planned_adapter | new_alm_adapter | T1 |
| DQMC-22 | Finite-temperature ED thermodynamics at matched ensemble | planned_adapter | new_alm_adapter | T3 |
| DQMC-23 | Trotter extrapolation versus finite-dtau thermal oracle | planned_adapter | new_chain_adapter | T3 |
| DQMC-24 | Independent-chain equilibration and blocked-error calibration | planned_adapter | new_chain_adapter | T3 |
| DQMC-25 | Low-temperature numerical stability stress | planned_adapter | new_chain_adapter | T2 |
| ED-01 | Particle-sector basis completeness | planned_adapter | new_ed_adapter | T0 |
| ED-02 | Canonical fermion anticommutation | planned_adapter | new_ed_adapter | T0 |
| ED-03 | Full Hubbard dimer spectrum | existing_test | alm_suite | T0 |
| ED-04 | Dimer singlet-triplet gap | existing_test | alm_suite | T0 |
| ED-05 | Dimer ground-state double occupancy | existing_test | alm_suite | T0 |
| ED-06 | Atomic four-state grand partition function | existing_test | alm_suite | T0 |
| ED-07 | Noninteracting many-body energies | planned_adapter | new_ed_adapter | T1 |
| ED-08 | Independent Jordan-Wigner Hamiltonian action | planned_adapter | new_ed_adapter | T1 |
| ED-09 | Hermiticity of dense and matrix-free operator | planned_adapter | new_ed_adapter | T1 |
| ED-10 | Centered and uncentered energy conversion | planned_adapter | new_ed_adapter | T0 |
| ED-11 | Dense versus Lanczos low eigenvalues | planned_adapter | new_ed_adapter | T1 |
| ED-12 | Lanczos residual and nonconvergence reporting | existing_test | renyi_suite | T0 |
| ED-13 | Degenerate eigenspace rotation invariance | planned_adapter | new_ed_adapter | T1 |
| ED-14 | Particle-number conservation | planned_adapter | new_ed_adapter | T0 |
| ED-15 | Spin SU(2) symmetry | planned_adapter | new_ed_adapter | T1 |
| ED-16 | Momentum projectors and fermionic translations | existing_capability | alm_suite | T0 |
| ED-17 | Local spin contact expectation | planned_adapter | new_ed_adapter | T0 |
| ED-18 | Spin correlation sum equals total S squared | existing_capability | alm_suite | T1 |
| ED-19 | One-body density matrix representability | planned_adapter | new_ed_adapter | T0 |
| ED-20 | Hamiltonian energy from separate observables | planned_adapter | new_ed_adapter | T0 |
| ED-21 | Hellmann-Feynman interaction derivative | planned_adapter | new_ed_adapter | T1 |
| ED-22 | Grand-potential chemical-potential derivative | planned_adapter | new_ed_adapter | T1 |
| ED-23 | Lehmann spectral zeroth moment | existing_test | alm_suite | T1 |
| ED-24 | Spectral causality and imaginary-time reconstruction | existing_capability | alm_suite | T1 |
| ED-25 | Twist gauge equivalence and two-pi periodicity | existing_capability | alm_suite | T1 |
| DMRG-01 | Electron local-state and double-occupancy mapping | planned_adapter | new_dmrg_adapter | T0 |
| DMRG-02 | ED to MPS to ED coefficient round trip | existing_test | alm_dmrg_preflight | T1 |
| DMRG-03 | Nonlocal up-spin hopping fermionic sign | existing_test | alm_dmrg_preflight | T1 |
| DMRG-04 | Nonlocal down-spin hopping fermionic sign | existing_test | alm_dmrg_preflight | T1 |
| DMRG-05 | Independent real-space hopping reconstruction | existing_test | alm_dmrg_preflight | T1 |
| DMRG-06 | Random full-sector MPO action versus Fock Hamiltonian | existing_test | alm_dmrg_preflight | T1 |
| DMRG-07 | Onsite centered-interaction constant | planned_adapter | new_dmrg_adapter | T0 |
| DMRG-08 | MPS norm and gauge invariance | planned_adapter | new_dmrg_adapter | T1 |
| DMRG-09 | Particle-number QN and local density sum | existing_capability | alm_dmrg_run | T2 |
| DMRG-10 | Physical spin QN normalization | existing_capability | alm_dmrg_run | T2 |
| DMRG-11 | MPO energy versus dense coefficient contraction | existing_capability | alm_dmrg_run | T2 |
| DMRG-12 | Energy variance versus full residual norm | planned_adapter | new_dmrg_adapter | T2 |
| DMRG-13 | Ground-state energy against certified ED | existing_capability | alm_dmrg_run | T2 |
| DMRG-14 | Ground-manifold overlap against certified ED | planned_adapter | new_dmrg_adapter | T2 |
| DMRG-15 | Full local density profile against ED | existing_capability | alm_dmrg_run | T2 |
| DMRG-16 | Full double-occupancy profile against ED | existing_capability | alm_dmrg_run | T2 |
| DMRG-17 | Full spin correlation matrix against ED | existing_capability | alm_dmrg_run | T2 |
| DMRG-18 | Spin contact terms from local observables | existing_capability | alm_dmrg_run | T2 |
| DMRG-19 | Total spin from complete correlation sum | existing_capability | alm_dmrg_run | T2 |
| DMRG-20 | Bond-dimension convergence of energy and entropy | planned_adapter | new_dmrg_adapter | T2 |
| DMRG-21 | Extra-sweep convergence at fixed bond dimension | planned_adapter | new_dmrg_adapter | T2 |
| DMRG-22 | Cutoff tightening with ample bond capacity | planned_adapter | new_dmrg_adapter | T2 |
| DMRG-23 | Independent initial MPS agreement | existing_capability | alm_dmrg_run | T2 |
| DMRG-24 | Checkpoint serialization preserves the state | planned_adapter | new_dmrg_adapter | T1 |
| DMRG-25 | PBC versus APBC seam MPO and site-order permutation | planned_adapter | new_dmrg_adapter | T1 |
| RENYI-01 | Zero entropy of product states | planned_adapter | new_renyi_adapter | T0 |
| RENYI-02 | Spin-singlet Bell pair entropy | planned_adapter | new_renyi_adapter | T0 |
| RENYI-03 | Maximally entangled spinful-site entropy | planned_adapter | new_renyi_adapter | T0 |
| RENYI-04 | Maximally mixed thermal subsystem | existing_test | renyi_suite | T0 |
| RENYI-05 | Schmidt spectrum versus explicit partial trace | planned_adapter | new_renyi_adapter | T0 |
| RENYI-06 | Pure-state complement equality | planned_adapter | new_renyi_adapter | T1 |
| RENYI-07 | Empty subsystem and full thermal purity | planned_adapter | new_renyi_adapter | T1 |
| RENYI-08 | Fermionic subsystem sign and noncontiguous cut | planned_adapter | new_renyi_adapter | T0 |
| RENYI-09 | Free-fermion correlation-matrix purity | planned_adapter | new_renyi_adapter | T1 |
| RENYI-10 | Interacting ED and DMRG ground-state entropy | planned_adapter | new_renyi_adapter | T2 |
| RENYI-11 | Known finite-temperature and Trotter ensemble separation | existing_test | renyi_suite | T1 |
| RENYI-12 | Grover rank-one determinant lemma | existing_test | renyi_rank1 | T0 |
| RENYI-13 | Single-proposal Grover tracker determinant and inverse | existing_test | renyi_rank1 | T0 |
| RENYI-14 | Sequential accepted Grover updates and rebuild clock | existing_test | renyi_rank1 | T0 |
| RENYI-15 | Rejected Grover proposal preserves accepted state | existing_test | renyi_rank1 | T0 |
| RENYI-16 | Explicit fallback on singular or ill-conditioned Grover matrix | existing_test | renyi_rank1 | T0 |
| RENYI-17 | Direct and optimized same-seed Markov trajectory | existing_test | renyi_rank1 | T1 |
| RENYI-18 | Replica stream reproducibility and independence | existing_test | renyi_suite | T0 |
| RENYI-19 | Adjacent lambda-ratio telescoping and endpoint rules | existing_capability | renyi_suite | T0 |
| RENYI-20 | Bidirectional lambda consistency in Monte Carlo | planned_adapter | new_renyi_adapter | T3 |
| RENYI-21 | Small interacting entropy versus exact replicated thermal target | planned_adapter | new_renyi_adapter | T3 |
| RENYI-22 | Nonlinear entropy uncertainty and shared-node covariance | planned_adapter | new_renyi_adapter | T1 |
| RENYI-23 | Strict raw-data schema and complete grid provenance | existing_test | renyi_suite | T0 |
| RENYI-24 | Launcher resume and atomic-output integrity | existing_test | renyi_suite | T0 |
| RENYI-25 | Low-temperature three-method entropy bridge | planned_adapter | new_renyi_adapter | T3 |

## Independent references

[Grover (2013)](https://arxiv.org/abs/1307.1486) supplies the Gaussian reduced-density-matrix construction and determinant approach. [Broecker and Trebst (2014)](https://arxiv.org/abs/1404.3027) supplies an independent replica formulation and Hubbard benchmarks at zero and finite temperature. [Assaad, Lang and Parisen Toldin (2014)](https://arxiv.org/abs/1311.5851) addresses stabilization. The [ITensorMPS DMRG FAQ](https://itensor.github.io/ITensorMPS.jl/stable/faq/DMRG.html) recommends ED comparisons and convergence checks; small variance alone does not establish a ground state. Exact finite-dimensional formulas and tolerance choices in this catalogue are derived validation specifications, not numbers attributed to those papers.

## Audit record per execution

For each ID, store `status`, `source_hash`, `runner_hash`, resolved parameters and ensemble, oracle hash, measured values, absolute/relative errors, statistical interval and comparison count, ESS/blocking diagnostics, Slurm job and array index, start/end UTC, runtime/maximum RSS, output hashes and failure reason. Link all IDs covered by a shared run to that same run record. Preserve the distinction between one case, one parameter arm, one assertion and one Slurm job.
