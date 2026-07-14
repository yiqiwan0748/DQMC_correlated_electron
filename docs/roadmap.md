# Project Roadmap

## Milestone 1: Noninteracting Bloch Model

Status: implemented by the files covered in this milestone.

- Define immutable hopping parameters and published defaults.
- Construct the two-orbital Bloch Hamiltonian in the `(d_xz, d_yz)` basis.
- Compute its two analytically ordered band energies.
- Validate high-symmetry points, generic momenta, orbital-exchange symmetry,
  and agreement with an independent symmetric eigensolver.

## Deferred Work

The following capabilities are intentionally outside Milestone 1 and remain
unimplemented:

- a real-space kinetic matrix;
- DQMC evolution and stabilization;
- Hubbard-Stratonovich fields and interactions;
- observables, including a density of states;
- plotting and visualization.

Each later milestone should introduce its own correctness tests before adding
the next layer of the simulation.
