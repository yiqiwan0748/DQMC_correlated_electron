# 100 external scientific benchmark comparisons

The current campaign contains **100 benchmark comparisons**, with repeated parameter families identified explicitly. It uses an independent established package and exact analytical reference methods. Internal unit tests, seed checks, schema checks and optimized-versus-original kernel checks do **not** count toward these 100.

| Reference family | Comparisons | Target calculations |
|---|---:|---|
| QuSpin exact canonical ground states | 64 | Both the user's ED and DMRG, on the same finite chain and sector |
| QuSpin grand-canonical finite-time-step thermal states | 16 | The user's thermal ED and interacting DQMC |
| Exact free-fermion Fourier/Gaussian formulas | 12 | The user's noninteracting DQMC Green functions and Rényi estimator |
| Exact four-state Hubbard atom | 8 | The user's thermal ED and interacting DQMC at zero hopping |

These are not 100 independent numerical algorithms or 100 independent publications. For example, eight observables/subsystems are compared across eight specified ground-state fixtures. Each case asks a separate observable/model/ensemble question; different random seeds are independent samples of one case.

The canonical definitions are [validation_100.json](validation_100.json), [validation_100.csv](validation_100.csv) and [benchmark_manifest.toml](benchmark_manifest.toml). The executable consumer is [run_external_validation.jl](run_external_validation.jl). No execution is asserted by the manifest: all rows begin as `NOT_RUN`. Actual receipts are separate TOML files from the execution host.

## What is compared

The eight ground fixtures are `(L,U,boundary) = (4,0,APBC), (4,2,PBC), (4,2,APBC), (4,4,PBC), (4,4,APBC), (4,8,PBC), (6,2,PBC), (6,4,APBC)`, all with `t=1`, `Nup=Ndown=L/2` and fixed particle number. Their observables are total uncentered ground energy, total double occupancy, longitudinal spin correlations at separations one and two, up-spin one-body correlations at separations one and two, and second Rényi entropies of the first site and the first half-chain. The independent exporter requires a unique ground state in the specified sector. The `U=0,L=4,APBC` fixture is closed shell.

The two thermal fixtures are `L=4,t=1,U=2,beta=1,dtau=0.1` and `L=4,t=1,U=4,beta=2,dtau=0.1`, both centered chemical potential zero and PBC. They use the same eight observables/subsystems. The reference is the normalized asymmetric product `[exp(-dtau V) exp(-dtau K)]^(beta/dtau)`, with the same ordering as the production DQMC. Its matrix need not be Hermitian at finite `dtau`; the entropy calculation uses `Tr(rho_A^2)`, not a Hermitian eigensolver applied to an asymmetric matrix.

The four free fixtures are `(L,beta,boundary) = (8,0.5,PBC), (8,4,APBC), (6,2,PBC), (6,1,APBC)`. Each compares total kinetic energy, the nearest-neighbor up-spin one-body correlator and half-chain Rényi entropy. The two atomic fixtures have `L=4,t=0`, `(U,beta,mu_centered)=(2,1,0)` and `(4,0.5,0.3)`. Each compares centered internal energy, total particle number, total double occupancy and two-site-subsystem Rényi entropy.

## Independent sources and formulas

QuSpin constructs and diagonalizes its own spinful fermion Hamiltonian and performs its own partial trace; it shares no implementation with the user's Julia ED or ITensor MPO. Its [spinful-fermion documentation](https://quspin.github.io/QuSpin/generated/quspin.basis.spinful_fermion_basis_1d.html) defines the operator/basis API. The package is described in [Weinberg and Bukov, SciPost Physics 7, 020 (2019)](https://doi.org/10.21468/SciPostPhys.7.2.020). The exporter and detailed provenance note are in `../external_oracles/`; the receipt fingerprints installed distribution versions and wheel records as well as exporter and input hashes.

For a free chain, the independent reference uses `epsilon_k=-2t cos(k)`, with `k=2pi n/L` for PBC and `k=2pi(n+1/2)/L` for APBC. With `f_k=1/[1+exp(beta(epsilon_k-mu))]`, the one-body correlator is the exact finite Fourier sum and total kinetic energy is `2 sum_k epsilon_k f_k`. If `nu_a` are eigenvalues of the restricted one-spin correlation matrix, spinful `S2=-2 sum_a log[nu_a^2+(1-nu_a)^2]`. This follows the correlation-matrix construction of [Peschel, Journal of Physics A 36, L205 (2003)](https://arxiv.org/abs/cond-mat/0212631), also reviewed by [Peschel and Eisler (2009)](https://arxiv.org/abs/0906.1663). At `U=0` the auxiliary-field noise and Trotter error vanish, so these are deterministic comparisons.

The atomic reference is derived directly from the four occupations of the [Hubbard Hamiltonian](https://doi.org/10.1098/rspa.1963.0204). For centered interaction, the grand energies of empty, up, down and double states are `U/4`, `-U/4-mu`, `-U/4-mu`, `U/4-2mu`. Normalize their Boltzmann weights to probabilities `p0,pup,pdown,pdouble`. For `L` independent sites, `N=L(pup+pdown+2pdouble)`, `D=L pdouble`, centered **internal** energy is `U(D-N/2+L/4)` and `S2(A)=-nA log(sum p^2)`. Internal energy excludes `-mu*N`; this distinction matters in the second atomic fixture. Since hopping is zero, the interaction and kinetic parts commute and there is no Trotter bias.

## Ensembles and energy conventions

DQMC uses `Hc=K+U sum(nup-1/2)(ndown-1/2)-mu_c N`. The fixed-sector ED and DMRG use `Hu=K+U sum D-mu_u N`. Grand-canonical equivalence requires `mu_u=mu_c+U/2`, giving `Hc=Hu+UL/4`. At the same numerical chemical potential in a fixed particle sector, `Eu-Ec=UN/2-UL/4`. The runner compares fixed-sector uncentered energies directly and converts thermal centered energies explicitly.

Finite-beta thermal states and pure ground states are separate targets. No ground-state Rényi value is substituted for a finite-temperature one. The original package's `L4,U4,beta16,dtau0.05,nA2` audit illustrates why: it records ground `S2=0.7521220360017`, exact thermal `0.7712028368091` and finite-step asymmetric thermal `0.7689951876109`. The broader campaign retains low-temperature and continuum extrapolation tests; this initial 100-case batch makes matched finite-size comparisons.

All entropies use natural logarithms. The DMRG energy and variance gates precede observable acceptance. The small-system ED target explicitly checks `norm(Hpsi-Epsi)<=1e-9`; DMRG requires finite variance in `[-1e-9,1e-7]`, norm error at most `1e-9`, matching particle number and consistent reported/contracted energy. Independent QuSpin ground-energy agreement is also required; a low variance alone would certify an eigenstate without certifying that it is the ground state.

## Run on ISAAC

The parent campaign prepares the pinned Julia environment and a separate Python environment containing QuSpin. Use compute-node batch allocations for the calculations. First generate the independent 80 QuSpin references using `../external_oracles/quspin_export.py` and its `cases.json`; the remaining 20 analytical references are computed directly from their formulas in the Julia consumer. All output paths must be new.

A deterministic staging run loads and compares ED/DMRG/thermal-ED/free-DQMC targets, explicitly leaving the interacting Monte Carlo arms unresolved:

```bash
RENYI_ROOT=/path/to/HubbardRenyiDQMC \
VALIDATION_ORACLES=/path/to/oracle_values.toml \
VALIDATION_OUTPUT=/new/campaign/deterministic_results.toml \
VALIDATION_RUN_DQMC=0 \
julia --project=/path/to/HubbardRenyiDQMC --startup-file=no \
  /path/to/validation_design/run_external_validation.jl
```

For all Monte Carlo arms, use a new output path and `VALIDATION_RUN_DQMC=1` (the default). `VALIDATION_ID_PREFIX=QT-` or `ATOM-` permits splitting independent jobs. `VALIDATION_MC_MEASUREMENTS` defaults to 4096 per independent chain and must be divisible by 64. There are eight preregistered seeds, 512 burn-in sweeps, 64 blocks and four lambda increments. Separate subsystem entropy arms have separate fixed grids. The selected physics parameters, source files, manifest, oracle file, Julia/thread settings, Slurm job and full paired bins are recorded in the receipts.

Compilation, queue delay, oracle generation and numerical runtime must be reported separately. Initial small ED/DMRG comparisons should fit a modest batch, while full fixed-budget DQMC entropy grids may require several CPU hours. These are unmeasured planning ranges, not finish promises; revise them after the first measured arm. The runner does not submit jobs or alter existing jobs.

## Quantitative acceptance and uncertainty

Ground-state ED and deterministic thermal-ED comparisons use absolute tolerance `1e-8`; analytical free comparisons use `1e-9` and atomic thermal-ED comparisons `1e-10`. DMRG uses `1e-7` for energy and `2e-5` for other observables. These tolerances are much larger than small-matrix roundoff but much smaller than the physical order-one scales; they are accuracy targets, not fitting allowances. Retain actual errors, even when they pass.

Interacting DQMC uses independent-chain means after paired block jackknife bias correction. The original campaign divided a 1% budget by 100 scalar comparisons. The corrected campaign uses alpha=0.01/400 and explicitly accounts for the original inspections in [its registered protocol](../isaac/corrected_monte_carlo_protocol.md). The runner evaluates the exact Student-t7 critical value from its elementary trigonometric tail; eight independent chain means give seven degrees of freedom. Its confidence interpretation assumes approximately normal sufficiently equilibrated independent chain means. Estimated effective independent bin count must be at least 128 across the eight chains; every chain's standard error must change by at most 50% when adjacent blocks are combined. A split-chain Rhat gate of at most1.05 on paired jackknife pseudovalues now also rejects between-chain discrepancies. These diagnostics do not prove equilibration. The broader campaign adds longer-chain, initial-state, autocorrelation and stopping-rule audits.

`PASS` requires the **entire** simultaneous confidence interval for estimate minus oracle to fit inside the stated equivalence margin: `0.05` for thermal energy and entropy, `0.03` for the remaining interacting observables (atomic entropy uses `0.05`). A wide interval, failed blocking/effective-bin gate, missing oracle, disabled run, dependency error or unconverged execution is `UNRESOLVED`. A resolved discrepancy beyond the margin is `FAIL`. The code does not silently take absolute values inside logarithms to conceal a negative sampled ratio. Failed or unresolved cases must remain visible; running longer is a new registered campaign, not a reason to overwrite the first receipt.

## Supporting work, excluded from the 100

`long_campaign_plan_100.*` preserves the earlier broader design, including deterministic regression tests and longer numerical stress/extrapolation tests. Its identifiers and source mappings are separate from this external benchmark batch. `supporting_regression_partial.jl` is an unfinished draft of kernel checks and has no execution driver. Neither contributes to the external count. `build_external_manifest.py` regenerates only the canonical external manifest; `build_long_campaign_plan.py` regenerates only the supporting plan.

## Sampler defect found during execution

The original atomic fixture with five imaginary-time slices exposed a reducible systematic Metropolis scan. Two of eight chains stayed in closed alternating-field cycles. The small reported within-chain errors and effective-bin counts did not establish correct exploration. The original receipts are diagnostic history, including their apparent passes. A working-copy patch adds independent50% idling at every local proposal in all three update paths. The immutable historical source remains available. Exact finite-state transition-graph evidence and fresh corrected-kernel comparisons are reported separately; see the [corrected protocol](../isaac/corrected_monte_carlo_protocol.md). The corrected summary must rerun all24 stochastic arms rather than reuse any original-kernel stochastic pass.
