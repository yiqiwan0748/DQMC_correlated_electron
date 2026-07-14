# Repository Guidelines

## Scope

- Work milestone by milestone. The current implementation covers only the
  noninteracting two-orbital Bloch Hamiltonian and its band energies.
- Do not add DQMC machinery, Hubbard-Stratonovich fields, a real-space kinetic
  matrix, density-of-states calculations, or plotting until a later milestone
  explicitly requests them.

## Frozen Reference

- `reference/square_hubbard.jl` is a frozen one-orbital reference. Never edit,
  reformat, rename, or delete it.

## Model Conventions

- The orbital basis is ordered as `(d_xz, d_yz)`:
  orbital 1 is `d_xz`, and orbital 2 is `d_yz`.
- Momenta are expressed in radians with lattice spacing set to one.
- Keep the Bloch Hamiltonian real and symmetric for this model.

## Engineering Requirements

- Prefer generic real arithmetic when it does not obscure the implementation;
  the default model parameters and ordinary evaluations use `Float64`.
- Add no package dependency without an explicit milestone requirement.
- Use explicit tolerances for floating-point tests and document why they are
  appropriate.
- Run the milestone test suite with `julia --project=. test/runtests.jl`.
