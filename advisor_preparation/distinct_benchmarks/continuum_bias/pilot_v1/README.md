# X01 excluded timing pilot, release v1

This release implements **only a timing and numerical-eligibility pilot**. X01 remains **NOT_RUN**, with **zero scientific completions**, whatever the pilot's measurements or eligibility status. There is no continuum fit, statistical-equivalence decision, production run or automatic permission to extend the budget. The separate `../PROPOSAL.md` retains the proposed three precision stages, two-observable margins and fixed fit checks for review before any full campaign.

## Frozen fixture and native paths

Use a four-site periodic Hubbard chain with hopping 1, repulsion 4, inverse temperature 1 and centered chemical potential 0. This is the full 256-state grand-canonical ensemble; half filling means mean particle number4, rather than a fixed spin/particle sector. Physical chemical potential is 2. Report **total physical internal energy**, `Hphys=K+4*sum(nup*ndown)`, and natural-log second Renyi entropy of sites1–2, including both spin species. This is finite-temperature subsystem entropy.

The pilot uses time steps 1/4 and 1/32, giving 4 and 32 slices, two independent bundles per time step, 512 burn-in rounds and 1024 measured rounds. Retain 64 bins of 16 measurements. Stabilize every 2 slices and use one complete update round per measurement.

Each bundle has one independent scalar-energy stream and four native two-replica entropy streams. The left lambda nodes are0,1/4,1/2,3/4, each with forward increment 1/4. Thus the full pilot has 4 scalar calls and 16 native entropy calls, with 36 actual replica/stream initializations. `seeds.toml` fixes every seed before execution; `generate_seeds.py` reproduces them with integer hashing and the unchanged native UInt64 splitmix derivation. Every actual initialization seed is positive and unique. The whole namespace is excluded from every future production stage, fit and confidence interval.

`target_pilot.jl` hashes the six native source/environment files in the manifest **before including them**. It calls actual `initializer` and `update_replica!` for energy, and unchanged `run_lambda_node_simulation` from `run/run_lambda_node.jl` for entropy. That native runner invokes `update_replica_rank1_fused!`; its lambda-zero path uses ordinary updates. No sampler, update kernel or production file is modified. New work is the scalar Wick/Green measurement adapter, raw-bin retention, timing and receipt checks. For one spin, `C=I-transpose(G)`; energy is `tr(K*(Cup+Cdown))+U*sum(diag(Cup).*diag(Cdown))`.

**Cumulative stability gate:** carry the returned `deltaG` into every next energy update, including burn-in. The native entropy runner likewise carries its maximum through burn-in, both replicas and measurements. Save its actual return `r.deltaG`, rather than a final-sweep surrogate. All 20 cumulative maxima must be finite, nonnegative and at most `1e-8`. The scalar cumulative bin history must be monotone, include the burn-in endpoint and end at the recorded maximum.

## Independent continuum and finite-step references

The oracle constructs the Hamiltonian and number operators independently with QuSpin 1.0.1 and the complete spinful basis, using pinned NumPy 2.4.6 and SciPy 1.17.1. It never imports the user's ED or DQMC code. The continuum density is normalized `exp(-beta*G)` with

`G=K+V=Hphys-2*N+4*I`,

`V=4*sum((nup-1/2)*(ndown-1/2))`.

Eigenstate and matrix-exponential continuum densities must agree. References are prepared for all four proposed time steps 1/4,1/8,1/16,1/32, even though this excluded pilot samples only the endpoints. The scientific endpoint of any later campaign remains the continuum density.

The source fixes asymmetric, exact-kinetic propagators. For the **HS-summed onsite-interaction transfer operator**, set `T=exp(-hV)exp(-hK)`, `D=exp(-hV/2)` and `S=D exp(-hK)D`. Then `T=D S D^-1` and normalized `rhoT=D rhoS D^-1`. Since the onsite even interaction splits into A and B, the spatial partial trace gives `rhoT,A=DA rhoS,A DA^-1`. Therefore their purities and S2 are exactly equal. This is not a single-configuration HS identity or an argument for nonlocal interactions.

The oracle explicitly checks full-state transfer/density similarity, direct integer transfer powers, reduced-state similarity, purity/S2 invariance, trace normalization, symmetric positivity, eigensolver residuals and mean half filling. It builds `VA_local=TrB(VA)/dim(B)` using QuSpin's linear public partial trace, avoiding a guessed reduced-basis order. **Do not Hermitize the asymmetric reduced operator**: native purity is `Tr(rhoA*rhoA)`, whose similarity to the positive symmetric state is checked. All declared identity errors must be at most `1e-10`; symmetric-state/transfer eigenvalues must be at least `-1e-12`. No raw asymmetric-eigenvalue positivity assumption is substituted.

For the stated real Hamiltonian, the symmetric transfer product has an even effective-Hamiltonian expansion near zero step. The purity identity gives leading O(h²) entropy bias. Also `Tr(Hphys*rhoT)=Tr(D^-1 Hphys D*rhoS)`; odd nested commutators of real symmetric V and Hphys are antisymmetric and trace to zero against symmetric rhoS. Thus physical internal energy has an even leading expansion too. These explicit derivations support the later proposal; this timing release performs no extrapolation.

## Pilot receipt analysis and failure handling

Require all 4 scalar rows, 16 entropy-node rows and 20 matching individual receipts, exact seeds/time steps/lambda values/cuts, finite 64-bin arrays, positive half-filled phase bins and source/runtime closure. Missing cases, a timeout, nonfinite data, oracle identity failure or excessive cumulative stability error make the pilot **UNRESOLVED**. They never promote a subset to success. Each completed call writes an individual TOML before it is appended to the aggregate, retaining completed work on interruption. Execution errors are retained, and the process exit codes are separate receipts.

Energy uses a paired numerator/phase delete-bin jackknife. Entropy uses the **whole** functional `-sum_k log(sum(numerator_k)/sum(phase_k))`; each deletion removes the same bin label in all four independent node streams, preserving each stream's numerator/phase pairing. Every full and leave-one-out entropy ratio must remain positive. Recompute with 32 adjacent merged bins as a descriptive blocking check. The two independent bundle values and their ordinary sample standard error are descriptive, as are MC-minus-finite-step and MC-minus-continuum differences. No pilot fit, precision margin, equilibrium claim or production confidence interval is evaluated. Native per-node statistics remain in raw receipts but do not replace the whole-functional estimate.

## Resources and handoff

The wrapper requests one CPU, 3500M and 10 minutes. Oracle timeout is 90 seconds, target timeout 480 seconds, leaving nominal 30 seconds for checks/reporting; the scheduler wall cap is authoritative. This is an unmeasured pilot estimate, roughly 2–8 minutes after start, not a runtime guarantee. It uses the existing pinned QuSpin environment and Julia 1.12.1/native project; it installs nothing. Root owns staging, live-queue checks and any submission. No job was submitted when this release was prepared.

Record GNU time resource receipts, all process exits and raw API runtimes. The summary includes only a provisional work-count projection to the separate proposed precision levels 16384/65536/262144, assuming 8 bundles and 2048 burn-in rounds. It scales measured pilot API time by bundle count, total slice count and update count. This also scales cold overhead and omits loading/oracle/reporting; mixing, precision and high-tail interval adequacy remain unknown. A favorable projection authorizes no production job. A pilot eligibility status cannot complete X01; both energy and S2 over the full frozen production grid and all future approved endpoint/adequacy gates would be required.

## Primary references and attribution

- [Versioned JDQMCFramework 1.2.6 propagator source](https://github.com/SmoQySuite/JDQMCFramework.jl/blob/v1.2.6/src/Propagators.jl), Git blob `d7cb41e1c15a018e3c85a2acdaa6f42ac3e05190`; source read before implementation. See also the [official API](https://smoqysuite.github.io/JDQMCFramework.jl/stable/api/).
- [QuSpin spinful basis documentation](https://quspin.github.io/QuSpin/generated/quspin.basis.spinful_fermion_basis_1d.html): full basis, spin-separated operators and spatial partial trace. [Hamiltonian API](https://quspin.github.io/QuSpin/generated/quspin.operators.hamiltonian.html): independent operators/matrices.
- [Grover, Entanglement of Interacting Fermions in Quantum Monte Carlo Calculations](https://arxiv.org/pdf/1307.1486): thermal auxiliary-field Gaussian reduced-state representation. The cut-local similarity and real-matrix even-order derivations above are specific to this release and are not attributed to the paper.

The native implementation and its library dependencies retain their existing authorship. This tranche adds validation/reference/receipt adapters; it does not establish a native continuum-extrapolation feature or universal DQMC reliability.
