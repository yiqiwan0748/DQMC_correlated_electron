# Thermal ensemble bridges: X02, X03, P07

Status: **NOT RUN**. This frozen tranche contains three distinct scientific questions. No numerical result is claimed by this protocol. Root owns deployment, the serial Slurm dependency, execution, and the canonical catalogue. The cap is one CPU, 3500 MB, ten minutes; this is an unmeasured cap, not a runtime prediction. The full Fock space has 256 states.

The fixture is a periodic four-site Hubbard ring, hopping 1, interaction 4, and a contiguous two-site subsystem. All grids, thresholds, versions, and source hashes are in `manifest.toml`; `SHA256SUMS` closes the executable source. Output paths are new, job-specific directories. Existing raw data and previous receipts are not overwritten.

## Independence and attribution

`oracle_quspin.py` independently constructs the entire physical Hamiltonian and particle-number operators with QuSpin 1.0.1. It performs SciPy full diagonalizations and QuSpin fermionic partial traces. It does not import the user's Julia matrices, eigenvectors, densities, RDMs, or prior benchmark results.

`target_thermal.jl` calls the existing user functions `full_fock_matrices`, `exact_thermal_density`, `renyi2`, and `fixed_sector_ground_state` in `analysis/audit_u4_l4_thermal_ed.jl`. The additive adapter introduces ensemble projectors, chemical-potential matching, response moments, and the cooling bound; it is not a new thermal-ED solver. That original module imports adjacent `compute_s2_ed.jl`; both source files are hash-checked **before** inclusion. The latter is the separately tested disjoint-subsystem production repair, SHA `f8f94cfc13e11bb8b5e64e8ac8f4ade5ed3a13b3976c83b5ab5db89d86ba5f28`. This tranche uses contiguous cuts and does not itself retest disjoint cuts.

The source dependency is `source/renyi_disjoint_corrected` on ISAAC. Neither the earlier immutable source nor the frozen TeNPy tranche is changed. Execution uses Julia 1.12.1 and the already separate QuSpin environment: QuSpin 1.0.1, NumPy 2.4.6, SciPy 1.17.1. The target requires only Julia standard libraries.

Primary references for the independent operator and trace APIs are the official [spinful fermion basis documentation](https://quspin.github.io/QuSpin/generated/quspin.basis.spinful_fermion_basis_1d.html), the official [Hamiltonian documentation](https://quspin.github.io/QuSpin/generated/quspin.operators.hamiltonian.html), and Weinberg and Bukov, [QuSpin Part II, SciPost Physics 7, 020 (2019)](https://arxiv.org/abs/1804.06782). QuSpin's `Nf=(Nup,Ndown)` fixes both spin-resolved counts; our canonical projector below instead fixes their sum. The oracle therefore starts from the full basis with no `Nf` restriction.

## Physical and centered chemical potentials

Define the physical Hamiltonian and total electron number by

\[
H_{\rm phys}=K+U\sum_i n_{i\uparrow}n_{i\downarrow},\qquad N=\sum_{i,\sigma}n_{i\sigma}.
\]

The existing source uses the centered grand generator

\[
G_c=K+U\sum_i(n_{i\uparrow}-\tfrac12)(n_{i\downarrow}-\tfrac12)-\mu_cN
=H_{\rm phys}-(\mu_c+U/2)N+UL/4.
\]

Thus `mu_physical=mu_centered+U/2`. The additive constant drops out of normalized Gibbs densities, but must be removed when reporting physical energies. The target builds `Hphys` once from the source matrix at **centered mu zero**, adding `U*N/2-UL/4`. For every requested physical chemical potential it separately obtains the thermal density from the original source at `mu_centered=mu_physical-U/2`. It never adds only `U*N/2` to a generator that still contains a nonzero centered chemical-potential term.

## X02: when is the thermal entropy close to the ground-state value?

At physical chemical potential 2, inverse temperatures are 0.5, 1, 2, 4, 8, 16, 32, 64, 128. The oracle first checks that the global grand-canonical ground state is unique, has total number 4 with two electrons of each spin, and has an eigenvector residual at most `1e-10`. Physical mu must lie at least `1e-6` inside the adjacent total-number charge-gap interval. The global first excitation must be separated by more than `1e-8`. The target also compares its existing fixed-spin ground-state routine with this global ground state.

Let the total excited-state Gibbs probability be `epsilon=sum(p[excited])`, evaluated directly rather than subtracting the ground weight from one. For a unique ground state,

\[
\rho_\beta=(1-\epsilon)|0\rangle\langle0|+\epsilon\tau,
\quad \|\rho_{A,\beta}-\rho_{A,0}\|_1\le2\epsilon.
\]

Contractivity of the partial trace and `||rhoA_beta+rhoA_0||_infinity<=2` imply the conservative bound

\[
|p_\beta-p_0|\le4\epsilon,\qquad p=\operatorname{Tr}\rho_A^2.
\]

If `B=4*epsilon<p0`, the second Renyi entropy in natural logarithms obeys

\[
|S_{2,\beta}-S_{2,0}|\le-\log(1-B/p_0).
\]

The lower-purity side gives the larger absolute logarithmic change. The bound is marked invalid when `B>=p0`; `-1` is only a finite serialization sentinel for an invalid bound. No monotonicity of the entropy is assumed.

A declared cooling point must satisfy this rigorous bound at most `0.001`, measured absolute entropy bias at most `0.001`, and leakage outside total number 4 at most `0.0001`. Both backends must provide at least one qualifying point. Statistical Monte Carlo uncertainty is absent here. This finite four-site certificate does **not** certify the historical ten-site, inverse-temperature-40 calculation in the research record.

## X03: canonical fixed total number versus a matched-mean grand ensemble

At inverse temperatures 0.5, 1, 2, 4, 8, the canonical trace projects onto **total** number 3 and includes all spin projections. Its dimension is `binomial(8,3)=56`. A fixed two-up/one-down sector has dimension 24 and is explicitly reported separately; it cannot replace the canonical ensemble.

Each backend independently bisects physical mu in `[-8,8]` until the bracket width is at most `1e-12`, then verifies mean number 3 within `1e-9`. The susceptibility `d<N>/dmu=beta*Var(N)` must exceed `1e-8` to resolve the root. The GCE and canonical energy per site, subsystem purity, and entropy are each compared against the corresponding independent oracle. Their differences quantify a finite-system ensemble bias; they are **not** required to vanish. The canonical weight of the fixed two-up/one-down sector is recorded to make the spin distinction visible.

## P07: differentiate physical energy at fixed physical chemical potential

Physical mu is held at 0.7, hence centered mu at -1.3, at inverse temperatures 0.5, 1, 2, 4, 8. With `Gphys=Hphys-mu_physical*N`, direct differentiation of the normalized trace gives

\[
C_\mu=\left.\frac{d\langle H_{\rm phys}\rangle}{dT}\right|_{\mu_{\rm phys}}
=\beta^2\operatorname{Cov}(H_{\rm phys},G_{\rm phys})
=\beta^2[\operatorname{Var}(H_{\rm phys})-\mu_{\rm phys}\operatorname{Cov}(H_{\rm phys},N)].
\]

By comparison, `beta^2*Var(Gphys)` differentiates the mean grand generator, a different observable. Neither expression should be called the fixed-number heat capacity. The fixture must separate the two expressions by at least `1e-5` at one declared temperature, or remains unresolved.

Both backends independently check `C_mu` against the five-point centered temperature derivative

\[
\frac{E(T-2h)-8E(T-h)+8E(T+h)-E(T+2h)}{12h},\qquad h=0.001T.
\]

Physical mu remains fixed at all four displaced temperatures; no density-matching root is used in P07. The derivative tolerance is `1e-7`; the independent covariance identity and cross-backend moments use `1e-8`.

## Acceptance, execution, and limits

The comparator checks source and manifest receipts, pinned versions, finite numeric values, all declared grid points, density/RDM trace errors at most `1e-10`, Hamiltonian Hermiticity and number conservation. Each question yields one status: failed eligibility, absent receipts, wrong provenance, or unresolved spectral/root conditions produce `UNRESOLVED`; a resolved observable disagreement over `1e-8` produces `FAIL`; agreement after the question's gates produces `PASS`. The exact entropy-bound identity allows `2e-10` for floating-point roundoff. No passed energy arm is counted as a separate idea.

`run.sbatch` executes oracle, target, and comparator sequentially, without nested submissions. Root supplies `VALIDATION_ROOT` and `THERMAL_BRIDGES_TRANCHE` and chooses the serial dependency. `/usr/bin/time -v` records each scientific command's resource usage. No environment install is included and no simulation has been run as part of preparing these files. Syntax parsing and read-only mathematical/source review are recorded in `static_review.json`.
