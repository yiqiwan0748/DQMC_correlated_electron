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

## Deferred Work

The following capabilities are intentionally outside Milestone 2 and remain
unimplemented:

- a real-space kinetic matrix;
- DQMC evolution and stabilization;
- Hubbard-Stratonovich fields and interactions;
- interacting observables and measurements.

Each later milestone should introduce its own correctness tests before adding
the next layer of the simulation.
