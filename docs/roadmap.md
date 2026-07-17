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

Status: IN PROGRESS (Phase 1A). Feasibility chains have not been executed and
no average-sign result is claimed.

- Reuse the validated Milestone 5 state construction and update kernel without
  changing `src/two_orbital_dqmc.jl`.
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
