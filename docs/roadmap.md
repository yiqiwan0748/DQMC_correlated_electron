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

## Deferred Work

The following capabilities are intentionally outside Milestone 4 and remain
unimplemented:

- DQMC evolution, sweeps, and stabilization;
- Hubbard-Stratonovich fields and interactions;
- interacting observables and measurements;
- sparse storage, checkerboard decomposition, and performance optimization.

Each later milestone should introduce its own correctness tests before adding
the next layer of the simulation.
