# Distinct analytical benchmarks: ED and thermal ED

This additive tranche implements **16 distinct questions** from A01–A25. Parameter points and multiple observables inside one question do not become additional benchmark ideas. The manifest lists every grid before execution; nine unavailable or uncontrolled comparisons remain `NOT_RUN`. Neither the source package nor earlier receipts are modified.

Run on ISAAC with Julia standard libraries only:

```sh
RENYI_ROOT=/absolute/path/to/corrected/HubbardRenyiDQMC \
ANALYTIC_OUTPUT=/new/path/analytic_results.toml \
JULIA_NUM_THREADS=1 julia --startup-file=no run_analytic_tranche.jl
```

`ANALYTIC_IDS=A03,A04,...` selects a subset. Existing output paths are rejected. The driver writes a receipt after every idea, recording every fixture, measured value, independent analytical target, error, tolerance, timing and source hash. A missing API or exception is `UNRESOLVED`; an executed comparison outside tolerance is `FAIL`. A group passes only when every declared fixture and comparison completes and passes. A group left unselected remains `NOT_RUN`. No stochastic error bars or multiple-testing corrections apply to these deterministic calculations.

## The tested calculation is the user's source

The driver imports the actual `analysis/compute_s2_ed.jl` and `analysis/audit_u4_l4_thermal_ed.jl`. It obtains Hamiltonians through `build_space` and `apply_hamiltonian!`, calls the user's Lanczos solver, calls the user's thermal density routine, and uses the user's spatial fermionic partial trace and Rényi functional. Analytical targets are computed separately from finite-mode sums, exact dimer/atomic solutions, binomial counting and Hubbard theorems. The driver does not manufacture a second analytical answer and label that a package result.

New measurement adapters extract Lehmann correlators, susceptibility, spectral residues, spin classification and RDM eigenvalues from those source Hamiltonians or eigenvectors. These validate the Hamiltonian and the stated added adapter. They **do not certify a corresponding production DQMC or DMRG measurement API**. Those methods require their own executions and uncertainty checks. Sparse Jordan–Wigner annihilation operators are explicit observable definitions, not replacement Hubbard Hamiltonians.

Every Lanczos ground state is checked against the entire dense spectrum of the same tiny source Hamiltonian and its true residual. This certifies the minimum before comparison with the independent closed form; a small residual alone would not establish a ground state. Thermal computations use the exact exponential, with no finite-Trotter claim.

## Hamiltonian, ensemble and fermion conventions

ED uses `H=K+U D-mu_uncentered N`. Thermal ED uses `Hc=K+U(D-N/2+L/4)-mu_centered N`. Thus `mu_uncentered=mu_centered+U/2`; centered **internal** energy excludes the `-mu N` term. Grand partition functions include it. Both spin species and natural logarithms are retained. All grand-canonical comparisons include every particle sector. Canonical beta-zero entropy is evaluated in its actual fixed-number sector, never replaced by a grand-canonical identity state.

The imported basis orders all up modes before all down modes. Observable Jordan–Wigner signs use that order, and the source partial trace supplies the sign for grouping spatial modes. `Gplus_ij(tau)=<c_i(tau)c_j†>` is positive at `0+`, and the endpoint `beta` means `beta-minus`; antiperiodic extension outside this interval is not silently substituted.

The native L2 periodic builder counts two directed ring bonds. The dimer adapter supplies `t_input=t_physical/2`, producing exactly one physical bond with coefficient `-t_physical`; this boundary/normalization conversion is recorded in every dimer fixture. A23 removes the seam by averaging the actual PBC and APBC Hamiltonians at L3. It tests the open three-site bipartite chain, not an unavailable general-graph API. The full N=3 sector is diagonalized, so both spin components and the lowest multiplet are included.

## Independent references and acceptance

- **A01–A02:** the four atomic energies are `U/4, -U/4-mu, -U/4-mu, U/4-2mu`. Their Boltzmann weights give N, D, centered internal E and log Z. The two spin-up spectral poles are `-U/2-mu` and `U/2-mu`, with residues `1-<n_down>` and `<n_down>`. Every supported transition, both aggregated residues and six imaginary-time points are compared. Absolute tolerance is `1e-10`.
- **A03–A05:** for one physical dimer bond, `R=sqrt(U²+16t²)`, `E2=(U-R)/2`, `E1=-|t|`, `E3=U-|t|`, spin gap `(R-U)/2`, charge gap `R-2|t|`. Singlet and triplet are checked by full S², not only Sz. With total doublon weight `D=(1-U/R)/2`, the one-site probabilities are `D/2,D/2,(1-D)/2,(1-D)/2`. Energies/gaps use `1e-9`; probabilities and entropy use `1e-10`. These follow the symmetric exact solution discussed by [Carrascal et al.](https://arxiv.org/html/1502.02194).
- **A06–A10:** independent momenta are `(2pi*m+phi)/L`, with phi=0 or pi, dispersion `epsilon=-2t cos(k)` and `f=1/(1+exp(beta*(epsilon-mu)))`. Density, energy, partition function, entropy, full imaginary-time matrix, all spectral poles/residues, compressibility and the discrete connected Lindhard response are distinct physical questions. Density response includes both spins and division by L. The `q=0,omega=0` response subtracts the disconnected mean and uses the exact coincident-energy limit. Tolerances are `1e-10` for thermodynamics and `1e-9` for dynamical/response quantities. These elementary free-mode results are derived independently here; a package documentation link in the parent catalogue is context, not a claim that those formulas were imported from it.
- **A11:** every restricted free correlation eigenvalue contributes factors nu and 1-nu to the complete many-body RDM spectrum, independently for both spins. The complete sorted spectrum must have L1 error at most `1e-9`; comparing only S2 would not satisfy this question. The Gaussian RDM construction follows [Peschel](https://arxiv.org/abs/cond-mat/0212631).
- **A12:** a fixed-spin-sector local configuration with j occupied A sites has probability `choose(L-nA,N-j)/choose(L,N)`. Summing its squared probability over `choose(nA,j)` configurations gives each spin's purity; their product gives the spinful purity. The source thermal routine is actually evaluated at beta=0 on its sector Hamiltonian, embedded in full Fock space, and spatially reduced. This is direct combinatorial counting, with `1e-10` tolerance.
- **A13:** the grand-canonical beta-zero state has `S2=nA*ln4` and purity `4^(-nA)`, independently of finite U and mu. The source exponential/partial trace must reach that limit to `1e-10`.
- **A14:** with Ndown=0, the uncentered interaction vanishes. Closed-shell filled free modes give the complete one-body projector and energy for arbitrary tested U. Energies/projectors use `1e-9`, and doublon occupancy must be at most `1e-12`. A finite gap between the last occupied and first empty free modes is required before comparing one pure ground-state projector.
- **A22:** the staggered eta tower at m pairs has uncentered E=mU. Its energy residual, dense excited-eigenspace overlap and full pair density matrix are checked to `1e-9`. Diagonal pair density is m/L; off diagonal it is `(-1)^(i+j)*m*(L-m)/(L*(L-1))`. These are specified excited states, not ground states. The reference is [Yang's eta-pairing construction](https://doi.org/10.1103/PhysRevLett.63.2144).
- **A23:** the connected open L3 chain has sublattice sizes 2 and 1, U>0 and half filling. Its complete lowest spin multiplet must have S=1/2, hence S²=3/4, multiplicity two, vanishing S² variance and a commuting spin operator. The entire fixed-N sector is included. This instantiates [Lieb's theorem](https://doi.org/10.1103/PhysRevLett.62.1201), including its stated assumptions, rather than guessing S from Sz. Tolerance is `1e-9`; multiplet dimension is exact.

The largest full-Fock matrix is 256×256 (L4). L6 polarized tests use a small fixed sector. All grids are finite and contain no uncontrolled thermodynamic extrapolation. Expected elapsed time is a few minutes including compilation on one ISAAC CPU; this is an estimate until an actual receipt exists.

## Provenance and interpretation

The receipt hashes the two imported source files and every tranche file. Importing the thermal audit also imports its nested copy of `compute_s2_ed.jl`, already included in the hash closure. Only Julia standard libraries are used. The corrected sampler can be the enclosing package, but passing these ED comparisons does not independently execute or certify that sampler. Scientific pass counts are **distinct benchmark IDs**, while fixture/comparison counts are supporting detail.
