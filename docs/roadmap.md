# Project Roadmap

## Milestone 1: Noninteracting Bloch Model

Status: implemented.

- Define immutable hopping parameters and published defaults.
- Construct the two-orbital Bloch Hamiltonian in the `(d_xz, d_yz)` basis.
- Compute its two analytically ordered band energies.
- Validate high-symmetry points, generic momenta, orbital-exchange symmetry,
  and agreement with an independent symmetric eigensolver.

## Milestone 2: Noninteracting Band Analysis

Status: implemented.

- Sample the configurable `Γ-X-M-Γ` band path and its high-symmetry ticks.
- Sample both bands once on a periodic uniform momentum grid.
- Build a normalized per-spin histogram DOS with optional Gaussian smoothing.
- Compute monotone zero-temperature filling including spin degeneracy.
- Locate the two principal internal DOS peaks and map their fillings under all
  three candidate doping conventions.
- Infer `delta = (n-2)/4` as the journal convention by matching both reported
  values, while recording the numerical ambiguity of that inference.
- Compare `Nk = 200`, `400`, and `800`, and generate reproducible CSV,
  metadata, band-structure, DOS, and filling outputs in
  `results/noninteracting/`.
- Test path geometry, Milestone 1 endpoint energies, DOS normalization, filling
  limits and monotonicity, and deterministic fixed-input output.

## Milestone 3: Finite-Size Real-Space Kinetic Matrix

Status: implemented.

- Define checked orbital-fast forward and inverse indexing on an `L x L`
  periodic lattice.
- Assemble the dense one-spin `2L^2 x 2L^2` kinetic matrix from the four
  positive-displacement orbital blocks and their reverse hoppings.
- Accumulate coincident periodic bonds exactly, including the `L=2` aliases.
- Project the real-space matrix onto every allowed finite-lattice momentum and
  reproduce the validated Bloch Hamiltonian element by element.
- Match the complete real-space spectrum to the union of the two Bloch bands
  for `L=2,3,4,5,6`.
- Generate reproducible validation errors, spectral bounds, and convention
  metadata in `results/realspace_kinetic/`.

## Milestone 4: Finite-Temperature Free Fermions

Status: implemented.

- Apply the chemical-potential convention `K_mu = K - mu*I` to the validated
  real-space hopping matrix without changing its hopping or indexing rules.
- Construct dense forward and inverse slice propagators from a symmetric
  eigendecomposition.
- Evaluate the equal-time DQMC Green matrix and complementary density matrix
  with branch-stable spectral Fermi functions rather than a direct inverse.
- Report spin-degenerate orbital fillings and total filling per physical site.
- Evaluate one- and two-spin grand `log Z` with a stable softplus expression
  and record positivity of the noninteracting free weight.
- Match finite-lattice real-space fillings to direct occupations of both Bloch
  bands for `L=2,3,4,5,6`, including spin degeneracy.
- Test Green identities and spectra, propagator inverses and products, filling
  limits and monotonicity, C4-related orbital densities, and finite-temperature
  crossover near both Milestone 2 van Hove chemical potentials.
- Generate the concise `L=2,4,6` validation grid and complete convention
  metadata in `results/free_fermion/`.

## Milestone 5: Correctness-First Two-Orbital DQMC Kernel

Status: implemented.

- Add the particle-hole-centered discrete Ising HS transformation with one
  field per site-orbital and time slice.
- Use dense exact asymmetric propagators with `K_mu = K - mu*I` and separate
  up- and down-spin Green matrices.
- Use JDQMCFramework public routines for stabilized equal-time Green
  calculation, propagation, determinant ratios, rank-1 updates, and periodic
  recomputation.
- Validate the scalar HS identity, the `U=0` free limit, explicit cyclic
  products, full determinant ratios, accepted updates, and seeded sweeps.
- Measure only sign, total and orbital density, orbital double occupancy,
  acceptance, numerical recomputation error, and runtime.
- Generate three tiny deterministic-seed smoke cases in
  `results/dqmc_smoke/`, with runtime explicitly treated as machine dependent.

## Milestone 6: DQMC Sign Feasibility Before Correlation Measurements

### Milestone 6A: Sign And Filling Diagnostics

Status: Milestone 6A tooling implemented. The registered three-seed
fixed-density calibration confirms `mu=-2.5202` at registered precision. The
fixed-mu sign-feasibility campaign and correlation measurements remain pending.

- Reuse the validated Milestone 5 state construction and update kernel without
  changing its HS-update, propagation, acceptance, or sampling behavior.
  Measurement-only thermodynamic additions to `src/two_orbital_dqmc.jl` do
  not alter that kernel.
- Before constructing or sweeping any `U > 0` calibration or sign chain,
  require the deterministic `U = 0`, `L = 4`, `beta*t = 4`,
  `mu/t = -1.897`, `dtau*t = 0.1` density-and-energy benchmark to pass. It
  uses fixed auxiliary fields with no burn-in, sweeps, or Monte Carlo
  sampling and compares DQMC estimators with an independent finite-lattice
  Bloch momentum-grid oracle.
- Validate the consistently named `density`, `doping`,
  `kinetic_energy_per_cell`, `shifted_interaction_energy_per_cell`,
  `internal_energy_per_cell`, and `grand_energy_per_cell`. Internal kinetic
  energy uses the bare `K`, grand energy subtracts `mu*density`, and the
  particle-hole-centered interaction constants are retained. The target
  finite-grid density is `1.820507020579327`, not the dense-momentum-grid
  value near `1.48466`.
- Keep three independently selectable phases: exact `U=0` controls, an
  interacting chemical-potential calibration scan, and fixed-chemical-
  potential sign diagnostics. Phase 1A restricts the configured study to
  `mu/t = -1.897`, near the upper-energy van Hove anchor, with calibration
  candidates centered on that meeting target. It is not claimed to be an
  interacting filling calibration.
- Record the raw per-configuration total-density estimator returned by
  `measure_basic_observables` and form
  `signed_filling_sum = sum(configuration_sign * raw_total_density)` with
  exactly one sign factor. The reweighted filling is a ratio to the measured
  sign sum and is omitted when that denominator is not usable.
- Estimate the average-sign standard error from contiguous block means and
  separately report early-versus-late block stability, the number of blocks,
  and whether the configured precision is reached. Retain and report a final
  short block rather than silently discarding measurements.
- Define the average-sign gate, and only that gate, by the strict condition
  `average_sign - 2*SE > 0.2` under the default configuration.
- Define numerical validity separately: configuration signs and measured
  densities must be finite, acceptance must be finite and in `[0,1]`, derived
  results must be finite when defined, and the maximum Green recomputation
  correction must not exceed `1e-8`.
- Define sampling sufficiency separately: at least eight blocks, early/late
  block drift at most two combined standard errors, and average-sign standard
  error at most `0.02`.
- Define filling resolution separately: the absolute average sign must exceed
  `1e-8` and its configured signal-to-noise ratio must be at least five.
  A zero standard error with nonzero average sign has the exact limiting SNR
  `+Inf`; this is a statistical limit, not a non-finite simulation result.
- Report `correlation_measurement_ready` only as the conjunction of the
  average-sign gate, numerical validity, sampling sufficiency, and filling
  resolution. Passing the average-sign gate alone never authorizes correlation
  measurements.
- Make `scripts/run_dqmc_sign_feasibility.jl` dry-run by default. It prints the
  complete selected case/seed table before any possible chain construction,
  requires `--run`, refuses output replacement without `--overwrite`, and can
  run a selected chain without writing through `--no-write` for benchmarking.

#### Registered Three-Seed Confirmation Rule

This fixed operational protocol applies only to:

- `L=4`
- `beta=4.0`
- `U=2.0`
- `dtau=0.1`
- `Ltau=40`
- `mu=-2.5202`
- target filling `n_target=1.46021875`

The existing strict-pilot calibration evidence is:

- `mu=-2.550`, `seed=6116`: `n=1.4393897420342978`,
  `SE=0.0013312914577584202`, `BELOW_TARGET`.
- `mu=-2.5202`, `seed=6117`: `n=1.4605219571679755`,
  `SE=0.001094369279004781`, `TARGET_COMPATIBLE`.
- `mu=-2.497`, `seed=6111`: `n=1.4764400882679722`,
  `SE=0.0010445205002070845`, `ABOVE_TARGET`.

The endpoints bracketed the target. At registration, the value `mu=-2.5202`
remained provisional; the completed three-seed adjudication below now confirms
it at the registered precision and only for the registered scope.

Seed 6117 is the previously observed pilot. Seeds 6118 and 6119 are the
prospectively registered independent confirmations. All three use:

- `burnin_sweeps=200`
- `measurement_sweeps=800`
- `sweeps_between_measurements=2`
- `block_size=100`
- `stabilization_period=5`

A chain is valid only if its exact registered workload completes, its CSV and
metadata agree and parse successfully, its sign, denominator, uncertainty, and
block-jackknife diagnostics are valid, its block-jackknife standard error
`s_i` is positive and finite, and every strict-pilot gate passes.

A missing, incomplete, inconsistent, or invalid chain may not be discarded,
replaced, or automatically rerun.

For each valid seed `i`, use the unrounded stored filling estimate `n_i` and
standard error `s_i`. Individual compatibility is:

`abs(n_i - n_target) <= 2*s_i`

For any seed set `S`, define:

- `w_i = 1/s_i^2`
- `n_pool(S) = sum(w_i*n_i for i in S) / sum(w_i for i in S)`
- `SE_pool(S) = 1/sqrt(sum(w_i for i in S))`

Require:

`abs(n_pool(S) - n_target) <= 2*SE_pool(S)`

for both:

- `S_confirm = {6118, 6119}`
- `S_all = {6117, 6118, 6119}`

For all three seeds define:

`Q = sum(w_i*(n_i - n_pool(S_all))^2 for i in S_all)`

Require:

`Q <= 5.991464547`

Use unrounded stored values throughout. All `<=` boundaries are inclusive.
Apply classifications in exactly this order:

1. Any missing, incomplete, inconsistent, or invalid chain:
   `NOT_ADJUDICABLE`.
2. All chains valid but `Q > 5.991464547`:
   `INCONCLUSIVE_HETEROGENEITY`.
3. All chains valid, all three individual conditions pass, both pooled
   conditions pass, and `Q <= 5.991464547`:
   `CONFIRMED_AT_REGISTERED_PRECISION`.
4. Otherwise: `NOT_CONFIRMED_AT_REGISTERED_PRECISION`.

Only the 6118+6119 confirmation pool is fully prospective; seed 6117 was
observed before this rule was registered.

This is a fixed operational rule. The inverse-variance pool is a fixed-effect
calculation; passing `Q` does not prove homogeneity, and the combined
`+/-2SE` requirements are not one exact 95% hypothesis test.

Seed 6119 remains required even if seed 6118 is unfavorable.

No classification automatically authorizes another seed, longer sampling,
exclusion or replacement of a chain, a rerun, or a change in `mu`. Do not
extrapolate the result to another lattice size, `beta`, `U`, or `dtau`.

The registered seed-6118 selector was:

```text
--phase=mu-calibration
--case=calibrate_target_mu_m2p5202
--seed=6118
--require-single
--strict-pilot
```

Its identifier is:

`calibrate_target_mu_m2p5202_seed6118_strict_pilot`

Its output directory is:

`results/dqmc_sign_feasibility/calibrate_target_mu_m2p5202_seed6118_strict_pilot`

The registered seed-6119 selector was:

```text
--phase=mu-calibration
--case=calibrate_target_mu_m2p5202
--seed=6119
--require-single
--strict-pilot
```

Its identifier is:

`calibrate_target_mu_m2p5202_seed6119_strict_pilot`

Its output directory is:

`results/dqmc_sign_feasibility/calibrate_target_mu_m2p5202_seed6119_strict_pilot`

Both confirmation chains completed at their registered output directories.

#### Registered Confirmation Outcome

The stored CSV/metadata pairs for seeds 6117, 6118, and 6119 agree, contain the
exact registered workloads, have positive finite block-jackknife standard
errors, and pass every strict-pilot gate. Using the unrounded stored values:

- all three individual compatibility conditions pass;
- the prospective confirmation pool has
  `n_pool=1.4596055452386897` and `SE_pool=0.0008287269978258003`;
- the all-seed pool has `n_pool=1.4599395349665640` and
  `SE_pool=0.0006606707522167680`; and
- `Q=0.4958851219864424 <= 5.991464547`.

The registered classification is therefore
`CONFIRMED_AT_REGISTERED_PRECISION`. Only the seed-6118 and seed-6119 pool is
fully prospective. This confirmation neither changes the fixed-effect caveats
above nor authorizes correlation measurements.

See the [full confirmation provenance](provenance/fixed_density_mu_m2p5202_confirmation.md)
and its [SHA-256 artifact manifest](provenance/fixed_density_mu_m2p5202_artifacts.sha256).

### Milestone 6B: Correlation Measurements

This submilestone remains unimplemented. It may be specified only after actual
Milestone 6A data pass all four diagnostic domains for the intended parameter
region. Its observables, estimators, correctness tests, and data products must
be approved separately.

## Deferred Work

The following capabilities are intentionally outside Milestone 6A and remain
unimplemented:

- spin structure factors and altermagnetic correlation functions;
- altermagnetic susceptibility and transport observables;
- production correlation-estimator and multi-observable jackknife workflows;
- sparse storage, checkerboard decomposition, and performance optimization;
- production-scale checkpointing, restart, and data-management machinery.

Each later milestone should introduce its own correctness tests before adding
the next layer of the simulation.
