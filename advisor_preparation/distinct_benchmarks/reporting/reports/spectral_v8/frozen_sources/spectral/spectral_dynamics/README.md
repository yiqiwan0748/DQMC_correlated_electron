# Spectral and time-response tranche: P08, P09, P10

Status: **NOT RUN**. This protocol defines three distinct questions on a periodic four-site Hubbard ring with hopping 1 and repulsion 4. There are no local numerical results. Root owns deployment and serial submission. Previous source releases, results, the CV, and the user package are unchanged.

This tests the actual user's ED Hamiltonian/ground-state APIs and full-Fock thermal-density API, with **new additive observable and time adapters**. It does not establish that native DQMC real-time dynamics, DQMC unequal-time spin measurements, a native user electron-spectrum routine, or DMRG dynamics work. No Monte Carlo sampling is performed. Each finite scientific question counts once; its time points, momenta, spins, transitions, and normalization checks are not extra ideas.

## Independent source and target attribution

`oracle_quspin.py` constructs the complete fermionic Hamiltonian, density operators, and annihilation operators independently with QuSpin 1.0.1. Its real-time reference uses QuSpin `exp_op`; its sector spectra use SciPy eigensolvers on QuSpin matrices; its finite-temperature response uses the complete QuSpin spectrum. It does not import user matrices, vectors, state labels, time series, operator adapters, or previous values.

`target_ed.jl` includes the existing `analysis/audit_u4_l4_thermal_ed.jl`, which imports adjacent `analysis/compute_s2_ed.jl`. The adapter calls existing `build_space`, `apply_hamiltonian!`, `ground_state_lanczos`, `full_fock_matrices`, and `exact_thermal_density`. It materializes tiny sector matrices by applying the actual Hamiltonian to unit vectors. Spectral propagation, fermion operators between sectors, Lehmann grouping, and the two-particle spin response are newly added adapter logic; existing user APIs do not natively supply these outputs in this tranche.

Both user source hashes are checked before inclusion. The active remote source is `source/renyi_disjoint_corrected`, with ED SHA `f8f94cfc13e11bb8b5e64e8ac8f4ade5ed3a13b3976c83b5ab5db89d86ba5f28` and unchanged thermal SHA `658c7c6bf465c70463450bca093d6e7eff589ac649d1532c98750e694bc1ddb5`. QuSpin 1.0.1, NumPy 2.4.6, SciPy 1.17.1, and Julia 1.12.1 are pinned. The Julia adapter uses standard libraries only.

The primary API references are official QuSpin documentation for [spinful fermions](https://quspin.github.io/QuSpin/generated/quspin.basis.spinful_fermion_basis_1d.html), [Hamiltonians](https://quspin.github.io/QuSpin/generated/quspin.operators.hamiltonian.html), [exponential action](https://quspin.github.io/QuSpin/generated/quspin.operators.exp_op.html), and the [Hubbard dynamics example](https://quspin.github.io/QuSpin/examples/example6.html). The [spinful operator source](https://quspin.github.io/QuSpin/modules/quspin/basis/basis_1d/fermion.html) maps the spin separator to the appropriate modes; the [exponential-action source](https://quspin.github.io/QuSpin/modules/quspin/operators/exp_op_core.html) shows the `Ns × time_count` output layout and its use of SciPy exponential multiplication. These were read during preparation; operator and evolution conventions are frozen before execution.

## Shared physical conventions

The physical Hamiltonian is

\[
H=-t\sum_{\langle ij\rangle,\sigma}(c^\dagger_{i\sigma}c_{j\sigma}+\mathrm{h.c.})+U\sum_i n_{i\uparrow}n_{i\downarrow}.
\]

The ring has four undirected bonds, each with both hopping directions exactly once. Real-time dynamics use this uncentered H at physical mu zero. Momentum phases use zero-based position `j=0,1,2,3`, although displayed sites are one-based. The full-Fock thermal source instead uses `Gc=H-mu_physical*N+U*L/4`, where `mu_physical=mu_centered+U/2`. The target constructs physical H from the source at centered mu zero, adding `U*N/2-U*L/4`. At nonzero thermal physical mu it obtains a separate density from the source at `mu_centered=mu_physical-U/2`; it does not mistake grand-generator eigenvalues for physical fixed-number energies.

## P08: charge-imbalance quench

The initial normalized pure state has two electrons of each spin: one doublon at sites 2 and 4, with sites 1 and 3 empty. No preparation Hamiltonian, imaginary-time cooling, or disorder averaging is implied. At time zero this product state is evolved by the uniform interacting H above.

The full time grid is 25 points from 0 through 6 in units of inverse hopping, with step 0.25. Observables are every site's total density, total doublon number, and

\[
I(t)=\frac{N_{\rm even}(t)-N_{\rm odd}(t)}{N},\qquad N=4.
\]

Even and odd refer to one-based site numbers. Thus `I(0)=+1`, total doublon number is initially 2, and the initial density vector is `[0,2,0,2]`. Doublon **number**, rather than number divided by L, is compared. Energy, norm, and both spin-resolved electron counts are conserved.

The target materializes the actual user 36-dimensional `(Nup,Ndown)=(2,2)` Hamiltonian and propagates through its full eigendecomposition. The independent reference applies QuSpin `exp_op` to the full 256-dimensional state and also compares that action with its own independent dense spectral propagation. There is no finite time step in either physical propagator; the displayed time spacing is an observation grid. Eigensystem residuals, unitary conservation, exponential-action accuracy, positivity, the exact initial profile, and a nonconstant imbalance are eligibility gates. Every time point must agree within `1e-8`.

## P09: interacting electron poles and residues

The initial state is the canonical `(2,2)` ground state of physical H. The target explicitly calls the actual user Lanczos API; its residual and agreement with the dense ground-state projector are checked. A unique ground state in that sector must be resolved. Both spin species and all four momenta `k=2*pi*m/4` are retained, using

\[
c_{k\sigma}=L^{-1/2}\sum_j e^{-ikj}c_{j\sigma}.
\]

For spin up, the complete addition/removal sectors are `(3,2)` and `(1,2)`; for spin down they are `(2,3)` and `(2,1)`. Each sector has dimension 24. The orbital annihilator in the target's spin-grouped Fock basis includes the parity of occupied preceding modes; a down-spin operator also includes the parity of **all** occupied up-spin modes. Adjoint operators supply creation. Same-mode canonical anticommutation is checked over the full Fock space for every spin/momentum. Both cross-spin anticommutators, `{c_up,c_down}` and `{c_up,c_down^dagger}`, are also checked on the full Fock space at frozen momentum indices 0 and 1. The latter checks detect a missing down-spin parity even though the initial ground state has even up-spin number; same-mode checks and this state's spectra alone would miss that error.

Frequencies are measured about physical mu 2, with uncentered physical energies:

\[
\omega_m^+=E_m(N+1)-E_0(N)-\mu,\qquad
\omega_m^-=E_0(N)-E_m(N-1)-\mu.
\]

Addition lies at positive frequency and removal at negative frequency for this frozen half-filled fixture. Residues are squared addition/removal matrix elements. There is no broadening, interpolation, or fitted line shape.

Individual eigenvectors within a degenerate multiplet have arbitrary gauge. All energies are sorted, and a cluster contains energies within `1e-8` of its **first** member; its residue is the sum of the complete multiplet's weights. The arithmetic mean energy is its representative. The comparison uses clustered frequencies, aggregate weights, and multiplicities. All complete-sector energy groups are retained, including groups with zero spectral weight; zero-weight groups are not described as visible electron poles. This avoids threshold-dependent deletion and eigenvector-by-eigenvector comparisons. The maximum cluster span, all 24 target states, nonnegative weights, and complete channel grid are checked.

For each spin/momentum, total addition plus removal weight is 1, removal weight equals `<c_k^dagger c_k>`, and addition weight equals the norm of `c_k^dagger|0>`. The first moment is checked against both the direct grand-generator matrix expression and the homogeneous Hubbard identity

\[
\int \omega A_\sigma(k,\omega)d\omega=-2t\cos k+U\langle N_{\bar\sigma}\rangle/L-\mu.
\]

The second moment is checked independently by applying `G-g0` to the addition/removal vectors and summing their squared norms. These checks fix the Fourier normalization, removal-frequency sign, energy zero, and completeness. The entire unbroadened spectrum must agree with the independent construction within `1e-8` for frequencies and aggregate residues.

## P10: finite-temperature imaginary-time spin response

This is a **continuum-time full grand-canonical** thermal trace, not a finite-Trotter target or a fixed-spin canonical trace. Inverse temperatures are 0.5, 2, and 6. Physical mu 1.3 corresponds to centered mu -0.7 in the existing source. Every momentum `q=2*pi*m/4` and every fraction `tau/beta=0,1/8,...,1` is retained: 108 grid points in total.

Use the physical spin operator, including its factor 1/2,

\[
S_q^z=L^{-1/2}\sum_j e^{-iqj}(n_{j\uparrow}-n_{j\downarrow})/2,
\]

and the connected correlator

\[
C(q,\tau)=\frac{\operatorname{Tr}[e^{-(\beta-\tau)G}S_q^z e^{-\tau G}S_{-q}^z]}Z
-|\langle S_q^z\rangle|^2,\qquad G=H-\mu N.
\]

The Fourier factor makes this equivalent to dividing a correlator of **unnormalized** spin Fourier sums by L. A second division by L would be wrong. No tau integral or extra beta factor is applied; this is the imaginary-time response function itself, not its integrated static susceptibility. Zero-field spin means vanish in the full trace, and the connected subtraction is nevertheless evaluated explicitly, including at q0.

Each backend evaluates the complete Lehmann sum with stable shifted Boltzmann factors. Its nonnegative spectral terms are invariant under degenerate-eigenbasis rotations, so individual thermal matrix elements are not compared. The thermal density is additionally checked against a normalized direct matrix exponential. The source's common centered constant cancels in the shifted normalized sum.

The gates require nonnegative response/spectral weights to numerical tolerance, correct equal-time endpoints, the **bosonic** KMS relation `C(q,tau)=C(-q,beta-tau)`, and a tau-independent uniform-spin response because total Sz is conserved. Both endpoints equal the direct `Tr(rho*S_q*S_-q)` minus the disconnected term, since these longitudinal operators commute. They are not fermionic antiperiodic Green functions. A nonzero imaginary-time variation at nonzero momentum is required so the fixture tests actual dynamics.

The equal-time sum rule is

\[
\sum_q C(q,0)=\sum_j[\langle (S_j^z)^2\rangle-\langle S_j^z\rangle^2]
=\langle N\rangle/4-\langle D_{\rm total}\rangle/2
\]

in this zero-field trace. It checks the spin-1/2 normalization and momentum normalization independently. Every beta, q, and tau point must agree deterministically within `1e-8`; no stochastic confidence band is invented for exact thermal traces.

## Acceptance, resources, and limits

`manifest.toml` fixes all grids, sources, versions, thresholds, and conventions. `SHA256SUMS` closes the complete local protocol. `compare.py` is standard-library postprocessing. It checks complete grids with no duplicates, exact receipt/source identities, finite values, all frozen eligibility conditions, and the cross-backend observables. Missing inputs, wrong provenance, incomplete spectra/grids, or unresolved normalization/eigensystem gates produce UNRESOLVED. Observable disagreement after eligibility gates produces FAIL. PASS applies only to this finite fixture and the specified adapter. Previous failures and source amendments must be retained by the parent ledger.

`run.sbatch` uses the existing separate QuSpin environment and pinned Julia runtime, one CPU, 3500 MB and a ten-minute unmeasured cap. It runs oracle, target, and comparison sequentially into a new job-specific output directory, records `/usr/bin/time -v` resources, and submits no further jobs. Root supplies `VALIDATION_ROOT` and `SPECTRAL_TRANCHE`. No CV or private personal files are needed for staging. A timeout or scientific execution error cannot establish a pass. Syntax parsing and independent read-only review are recorded separately; preparing this release executed no physics.
