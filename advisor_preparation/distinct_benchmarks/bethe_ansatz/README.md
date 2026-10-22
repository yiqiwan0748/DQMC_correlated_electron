# Bethe-ansatz reference pilot: A18, A19, A20

**Prepared, not run.** Only the bounded pilot is executable. It cannot pass a thermodynamic catalogue question. No local scientific computation, job submission, source change, or change to an earlier tranche was performed. Root owns review, staging, and serial scheduling.

The three questions are the half-filled bulk ground-state energy (A18), the full addition/removal charge gap (A19), and the bulk double occupancy (A20). They are distinct observables, even though the same finite-chain ground states supply them. Sizes, sectors, seeds, and convergence arms never count as additional questions.

## Primary reference and conventions

The reference is Lieb and Wu, [The one-dimensional Hubbard model: A reminiscence](https://arxiv.org/pdf/cond-mat/0207529), specifically Eq. 54 and Eqs. 64–65, on PDF pages 21 and 25. This is an author account of the analytical solution, published as Physica A 321 (2003), DOI [10.1016/S0378-4371(02)01785-5](https://doi.org/10.1016/S0378-4371(02)01785-5). The original discovery is [Phys. Rev. Lett. 20, 1445 (1968)](https://doi.org/10.1103/PhysRevLett.20.1445). The formulas below restore the hopping scale t from the paper's t=1 convention.

All target energies refer to the uncentered, zero-chemical-potential Hamiltonian

\[
H=-t\sum_{i,\sigma}(c^\dagger_{i\sigma}c_{i+1,\sigma}+\mathrm{h.c.})+U\sum_i n_{i\uparrow}n_{i\downarrow}.
\]

The user MPO already has this convention. No centered interaction or grand-generator shift is introduced. At half filling, a centered interaction would shift energy per site by `-U/4` and its U derivative by `-1/4`; confusing these would invalidate A18/A20. The finite addition/removal gap cancels a common constant and a term linear in total number, but the pilot's particle-hole check is explicitly for physical mu zero.

The independent oracle uses thermodynamic integrals. It is **not** a finite-size Lieb–Wu equation solver, and its numbers are not exact finite OBC-chain energies.

## Independent quadrature oracle

Writing `f(x)=1/(1+exp(U*x/(2*t)))`, the exact bulk energy per site is

\[
e_\infty=-4t\int_0^\infty \frac{J_0(x)J_1(x)}{x}f(x)\,dx.
\]

The paper gives `mu_minus=2t-4t*integral[J1(x)*f(x)/x]` and `mu_plus+mu_minus=U`. Consequently the **full** charge gap, not a single-holon half-gap, is

\[
\Delta_c=\mu_+-\mu_-=U-4t+8t\int_0^\infty \frac{J_1(x)}x f(x)\,dx.
\]

Hellmann–Feynman for the uncentered Hamiltonian and differentiation under this convergent integral give

\[
d_\infty=\frac{\partial e_\infty}{\partial U}
=2\int_0^\infty J_0(x)J_1(x) f(x)[1-f(x)]\,dx
=\lim_{L\to\infty}\frac1L\sum_i\langle n_{i\uparrow}n_{i\downarrow}\rangle.
\]

The factor 2 follows from the derivative of the Fermi denominator, `df/dU=-x*f*(1-f)/(2*t)`. The oracle additionally checks this analytic derivative against a five-point U derivative of the energy integral, with predeclared step 0.01.

`oracle_lieb_wu.py` imports no target matrices, MPS, user code, or prior results. It uses pinned NumPy 2.4.6 and SciPy 1.17.1, adaptive QUADPACK on `[0,40]`, and separately generated Gauss–Legendre rules of orders 128 and 256. The formulas use stable logistic functions. The removable `J1(x)/x` limit at zero is `1/2`.

For `a=U/(2*t)>0` and cutoff `X`, `|J_n(x)|<=1` and `f(x)<=exp(-a*x)` give absolute tail bounds:

\[
|R_e|\le\frac{4t e^{-aX}}{aX},\qquad
|R_{\Delta}|\le\frac{8t e^{-aX}}{aX},\qquad
|R_d|\le\frac{2e^{-aX}}a.
\]

The estimated numerical uncertainty is the maximum of QUADPACK's error estimate, disagreement with either fixed Gauss rule, and a `1e-12` roundoff floor, plus the analytic tail bound. It must not exceed `1e-9`; an integration warning is an execution failure. The tail bounds are rigorous, while floating-point quadrature/error estimates are not an interval-arithmetic proof. The fixed U4/t1 fixture avoids the severe subtraction problem in the weak-coupling charge gap. U near zero, attractive U, or t=0 are unsupported by this release.

The numerical APIs are described in the official [SciPy quad documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad.html) and [Gauss–Legendre documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.roots_legendre.html). The installed versions, rather than the current documentation version, are enforced in the receipts.

## Frozen pilot and actual target attribution

`pilot.toml` fixes open L8 and L12 chains at U4,t1 and physical mu0. At each size, half-filled, one-electron-removed, and one-electron-added states are obtained with independently seeded coarse/refined DMRG arms: chi 128/256, 12/18 sweeps, and cutoffs `1e-11`/`1e-13`. This is 12 finite-chain solves. The low-bond sweeps and noise ramp are fixed before execution. Initial random-MPS seed and complete sweep truncation diagnostics are recorded. Each completed solve is saved immediately to a new per-arm receipt; a timeout retains completed receipts.

For even `L=2m`, the sectors are

| Sector | Total N | Nup | Ndown | 2Sz |
|---|---:|---:|---:|---:|
| half | L | m | m |0|
| minus | L-1 | m | m-1 |1|
| plus | L+1 | m+1 | m |1|

The sectors have the minimum allowed absolute spin projection. In an SU(2)-invariant Hamiltonian they contain representatives of every allowed total-spin multiplet at the specified total number; no odd-particle Sz0 state is requested. Starting from alternating Up/Dn labels, removal changes one Dn to Emp and addition changes one Dn to UpDn before a seeded shuffle. The charge gap is `Eplus+Eminus-2*Ehalf` at the **same chain length**, never energies from different-length lattices. On this bipartite OBC chain, particle-hole symmetry combined with a spin flip requires `Eplus-Eminus=U`; this is a quality check, not the external gap oracle.

`pilot_target.jl` includes the actual existing `analysis/compute_s2_itensor_dmrg.jl` and calls its `hubbard_mpo`. It follows the user's QN-conserving random-MPS/DMRG/observer workflow; ITensorMPS supplies the optimizer. The additive adapter supplies N±1 initialization, local double occupancy, contracted-energy variance, and resource receipts. It does not claim a new optimizer or a modification of the user kernel. Source SHA, Project.toml SHA, Manifest.toml SHA, Julia 1.12.1, ITensors 0.9.30, and ITensorMPS 0.4.1 are enforced before execution. The remote source root is `source/renyi_disjoint_corrected`; the previously repaired ED source is not used by this MPO pilot.

The APIs match the official [ITensorMPS DMRG examples](https://docs.itensor.org/ITensorMPS/stable/examples/DMRG.html) and the user's source. The adapter reports direct `Nupdn` expectation values rather than estimating target double occupancy by a numerical derivative of target energies.

Quality gates check normalization, exact target electron/spin counts, reported versus contracted energy, variance per site, final truncation error, particle-hole energies, and coarse/refined stability of energy, gap, and double occupancy. Operator normalization is checked through full local profiles and their means. All thresholds are in `pilot.toml`; a quality failure remains unresolved and cannot be hidden by dropping an arm. Finite-chain offsets from the bulk oracle are diagnostic values and are never equality gates.

`summarize_pilot.py` always emits **zero scientific completions** and UNRESOLVED for A18/A19/A20, even when its separate pilot-quality field passes. Both source hashes and full case coverage are checked. The pilot allocation cap is one CPU, 3500 MB, 15 minutes, with a 13-minute target-command cap. These are unmeasured limits, not predicted requirements. Root must inspect timing, maximum bond dimensions, truncation/variance diagnostics and `/usr/bin/time -v` memory receipts before a larger-size allocation. There is no full-size runtime bound from L8/L12 alone.

## Predeclared full-size proposal; separate execution release required

`full_proposal.toml` records a prospective scientific design before the pilot: U4,t1 OBC, training sizes 16, 24, 32, 40, 48; held-out sizes 56 and 64; all three sectors and two independently seeded convergence arms. Pilot sizes 8/12 are excluded from all bulk fits. The full executor and allocation are **not** supplied or authorized by this release. A bridge-size timing/quality pilot and memory estimate are needed before choosing final chi/sweep caps and creating a separately reviewed full execution release.

Boundary contributions produce `1/L` terms in chain averages. The primary model is quadratic in `32/L`; a cubic model and dropping the smallest training size test sensitivity. Neither held-out size is used to fit or select the model. Every model must predict both held-out results within the observable-specific tolerances. Full-chain mean double occupancy, including boundaries, is extrapolated; a central-site value cannot be silently substituted. The charge gap has gapless-spin finite-size corrections, so a polynomial ansatz may fail; it is not a theorem about its remainder.

The proposed uncertainty is a deterministic empirical envelope containing model-intercept spread, leave-one-training-size-out intercept spread, training residuals, held-out prediction errors, propagated DMRG arm uncertainty, and oracle numerical uncertainty. A conservative proposed implementation is the sum of the maximum of the first four terms plus the propagated arm term and oracle error. All least-squares models use the training data only. This is **not** a statistical confidence interval or a rigorous thermodynamic error bar. Ground-state DMRG values are not independent Monte Carlo samples.

The proposed final equivalence test requires `abs(primary_intercept-oracle)+empirical_uncertainty <= margin`, and separately limits that uncertainty. Margins are 0.0002 for energy/site and double occupancy, 0.03 for the full charge gap. The maximum empirical uncertainties and held-out tolerances are separately predeclared in `full_proposal.toml`. Full uncertainty propagation, especially amplification of per-sector energy errors in the charge-gap difference, must be implemented and independently reviewed before freezing the full runner. Missing sizes, poor DMRG convergence, failed predictions, unstable extrapolation, or insufficient precision yield UNRESOLVED. No oracle-guided size deletion, changed polynomial, or widened threshold is allowed after seeing full-size data. Any extension needs a new protocol and preserved original attempt.

Even a successful future full release will validate a specified observable at one interaction and the stated empirical extrapolation precision, not every Hubbard parameter or all DMRG calculations. No finite-size computation proves the thermodynamic Bethe ansatz.
