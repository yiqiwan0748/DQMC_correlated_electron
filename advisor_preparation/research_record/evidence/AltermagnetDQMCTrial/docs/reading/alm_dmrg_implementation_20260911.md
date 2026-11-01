# Independent electronic DMRG implementation — 11 September 2026

The include-safe module and standalone runner are
`experiments/alm_dmrg_20260911/run_benchmark.jl`. They implement the
[registered exact-locking milestone](alm_exact_locking_plan_20260911.md),
using ITensor electronic MPS optimization independently of the previous
constrained-path solver. Implementation and package setup are not themselves
scientific validation; the saved numerical status determines the outcome.

## Fixed model, environment and scope

There are four fully open physical cells and eight spinful orbital modes,
ordered `(dxz,dyz)` within each cell. The model has six electrons, U=20,
and hoppings `(-1,-1.75,-.85,-.65)`. The only interaction is
`U*(nup-1/2)*(ndown-1/2)` on each orbital mode. The MPO is assembled from
the exact-reference one-spin K, including both hopping directions and both
spins. An independent bond assembly from the registered Bloch coefficients
checks K before use. Periodic wrapping is absent, so t4 survives.

The isolated environment pins ITensorMPS0.4.1 and ITensors0.9.30; their
tagged Project files declare Julia1.10-compatible1.x versions. Both package
versions are checked at runtime. No production dependency, production solver,
or frozen one-orbital reference is modified. Sources:
[ITensorMPS Project](https://raw.githubusercontent.com/ITensor/ITensorMPS.jl/v0.4.1/Project.toml),
[ITensors Project](https://raw.githubusercontent.com/ITensor/ITensors.jl/v0.9.30/Project.toml).

The runner requires one Julia thread and sets one BLAS thread. The detached
remote service enforces the registered4GiB,150%CPU,two-hour resource cap;
the runner does not replace that service-level limit with its own timer.
There is no implicit larger lattice, extra trial search, or adaptive sweep grid.

## Electron and ED ordering

An ITensor `Electron` site has states `(Emp,Up,Dn,UpDn)`, retaining double
occupancy. It conserves total Nf and Sz, equivalent here to conserving
Nup and Ndown. Its QN named `Sz` stores twice physical Sz; the operator
`"Sz"` has physical eigenvalues±1/2. No total-S=0 or separate orbital
population constraint is imposed. The current QN API is Abelian. Sources:
[Electron site documentation](https://itensor.github.io/ITensorMPS.jl/stable/IncludedSiteTypes.html),
[QN tutorial](https://itensor.github.io/ITensorMPS.jl/stable/tutorials/QN_DMRG.html).

The ED vector is up-major/down-minor, with all ascending up creation
operators preceding all ascending down creation operators. The local
Electron tensor instead interleaves modes, with a doubly occupied local
state `cup†*cdown†|0>`. Conversion multiplies each coefficient by
`(-1)^inversions`, counting occupied pairs `down_i,up_j` with `i<j`.
Both tensor-to-ED and ED-to-tensor conversions apply this parity. Full
coefficient tensors are split into MPS at cutoff0 and maxdim256, sufficient
for the maximum Schmidt rank of eight local-dimension4 sites.

Hopping uses `Cdagup/Cup` and `Cdagdn/Cdn` OpSum terms with the signed K
coefficient directly. ITensor supplies fermion ordering and Jordan–Wigner
strings. No manual string or extra overall hopping minus sign is added.
The onsite terms are `U*Nupdn-U/2*Ntot+U/4*Id`; the centered energy differs
from the uncentered one by−20 at the registered N6,U20,M8. Sources:
[OpSum API](https://itensor.github.io/ITensorMPS.jl/stable/OpSum.html),
[fermion sorting implementation](https://github.com/ITensor/ITensorMPS.jl/blob/712f18f675bbda7d0beaf82618a2d62bc6e77adb/src/opsum_to_mpo/opsum_to_mpo_generic.jl).

## Preflight actually performed

The runner accepts the plain-array
`results/literature_transfer/alm_exact_locking_20260911/dmrg_reference_target.jls`.
It checks geometry, model parameters, particle counts, mask uniqueness,
Hermiticity, ground-basis orthonormality, and eigenpair residuals.
The actual references have a unique Sz0 singlet and a unique Sz1 triplet
minimum; no selected-vector versus degenerate-ensemble ambiguity arises.

MPO action is checked on **fourteen vectors**, seven in each fixed sector:

- Three independently randomized vectors with full support over every ED
  basis configuration, using preflight RNG seeds94871 and94872 by sector.
- Two nontrivial product vectors, one containing a doubly occupied mode.
- One up-spin and one down-spin product fixture containing a nonlocal hop
  across an odd number of occupied same-spin modes. Its negative fermionic
  sign is checked analytically against the corresponding dense-H entry.

Each vector passes ED-to-QN-tensor-to-MPS and inverse round-trip checks.
The electronic MPO is applied with direct (`naive`) contraction, cutoff0,
and maxdim256; its recovered coefficients are compared with dense `H*v`.
Maximum component error must be≤2e−10, relative action error≤2e−11, and
norm leakage≤1e−8. K and the explicit fermionic-sign matrix entries use2e−12.
These are randomized full-sector action checks plus explicit fixtures;
**they are not a comparison of every Hamiltonian basis column**.

`--preflight-only` stops after those checks. `--run` performs the same
preflight before any DMRG stage and stops on a preflight exception. The
registered execution uses the latter mode, avoiding a redundant separate
preflight campaign. The script had a successful local syntax-only parse
before deployment; ITensor execution is assessed remotely.

## Frozen scientific schedule and outputs

There are exactly four arms: physical Sz0 and Sz1, each with seeds314159
and271828. Their product configurations differ, and an explicit Xoshiro RNG
randomizes each into a QN MPS of initial link dimension4. They optimize
independently. Within each arm, states continue through caps32,64,128,256
with six sweeps per cap, then six extra zero-noise sweeps at256.

Cutoff is1e−12. The local eigensolver uses tolerance1e−11,
Krylov dimension10, and maximum8 rebuild iterations. Initial-cap noise is
`[1e-6,1e-7,1e-8,0,0,0]`; subsequent cap stages use
`[1e-8,1e-9,0,0,0,0]`; the last extra stage is all zero. These are DMRG
optimization controls, not fields added to the Hamiltonian. Sources:
[DMRG API](https://itensor.github.io/ITensorMPS.jl/stable/DMRG.html),
[random MPS implementation](https://raw.githubusercontent.com/ITensor/ITensorMPS.jl/v0.4.1/src/mps.jl).

Every output directory must be new. Saved files include stage and
initialization-agreement CSVs, every spin matrix and local observable,
per-stage logs and checkpoints, four final MPSs in `finalstates.jls`, a
summary TOML, metadata, journal, and final status. The runner copies its
source, Project and Manifest into the output and records their hashes,
package entry-source hashes, versions, and reference hash. An exception
preserves its traceback and latest available state. A failed final gate
does not launch more optimization.

## Observables, normalization and acceptance

The measured full matrix is
`spin_dot=Czz+(Cpm+Cmp)/2`, where the three correlation matrices contain
physical `Sz*Sz`, `S+*S-`, and `S-*S+`. Its sum is total S². Contact terms
are included and independently compared with
`3/4*(density-2*double_occupancy)`. Sources:
[correlation API](https://itensor.github.io/ITensorMPS.jl/stable/MPSandMPO.html),
[variance formula](https://itensor.github.io/ITensorMPS.jl/stable/faq/DMRG.html).

For orbital-fast coefficients `a=(1,-1,1,-1,...)`, the runner's
`Sodd_scalar=a'*spin_dot*a/(3*Ncell)` is a **rotationally normalized**
scalar. The exact study's scalar normalized only by Ncell is three times
this field. `Sodd_z=a'*Czz*a/Ncell` agrees with the rotationally normalized
scalar in an isotropic singlet, but need not agree in the Sz1 state.
All registered matrix gates compare full `spin_dot` directly, without a
factor-of-three conversion. `Delta=sum(a.*<Sz>)/Ncell` is also recorded.

Primary final gates are absolute energy error≤1e−8, energy variance≤1e−8,
each spin-dot entry error≤1e−6, total S² error≤1e−6, and particle/Sz
errors≤1e−10. The moment-subtraction variance retains its raw sign and
allows only a roundoff-scale negative value within1e−8; an independent
ED-vector residual variance must also pass. Additional internal checks
cover energy from dense-H action, contact identities, and sector leakage.
Both initializations must agree within twice the single-reference tolerances.

Ground-subspace fidelity, density and double-occupancy errors are reported;
they are not silently promoted into additional preregistered physics gates.
Variance checks identify eigenstates, not necessarily the ground state;
the exact reference energy and fidelity address that distinction here.
Success validates this finite-cluster solver/observable calculation. It
does not settle the bulk magnetic competition or the U20 coexistence
candidate near density1.425: this benchmark has density1.5.
