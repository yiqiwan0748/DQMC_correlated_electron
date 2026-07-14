# Noninteracting Two-Orbital Model

## Basis And Conventions

The Bloch basis is ordered as

1. orbital 1: `d_xz`
2. orbital 2: `d_yz`

Momenta `kx` and `ky` are in radians, and the lattice spacing is one. In this
basis, the noninteracting Hamiltonian is the real symmetric matrix

```math
H_0(k) =
\begin{pmatrix}
\epsilon_x(k) & \epsilon_{xy}(k) \\
\epsilon_{xy}(k) & \epsilon_y(k)
\end{pmatrix}.
```

Its entries are

```math
\begin{aligned}
\epsilon_x(k_x,k_y) &= -2t_1\cos(k_x)-2t_2\cos(k_y)
                       -4t_3\cos(k_x)\cos(k_y), \\
\epsilon_y(k_x,k_y) &= -2t_2\cos(k_x)-2t_1\cos(k_y)
                       -4t_3\cos(k_x)\cos(k_y), \\
\epsilon_{xy}(k_x,k_y) &= -4t_4\sin(k_x)\sin(k_y).
\end{aligned}
```

## Default Parameters

| Parameter | Value |
| --- | ---: |
| `t1` | `-1.0` |
| `t2` | `-1.75` |
| `t3` | `-0.85` |
| `t4` | `-0.65` |

`ModelParameters()` and `DEFAULT_PARAMETERS` both provide these `Float64`
values. `ModelParameters` is immutable and can also hold another common real
type.

## Public API

- `bloch_hamiltonian(kx, ky, parameters=DEFAULT_PARAMETERS)` returns the
  `2 x 2` real symmetric matrix in the orbital convention above.
- `band_energies(kx, ky, parameters=DEFAULT_PARAMETERS)` returns a two-element
  vector of eigenvalues in ascending order. It uses the analytic eigenvalues
  of a real symmetric `2 x 2` matrix.

No chemical-potential shift or interaction term is included in `H_0`.

## Numerical Checks

Float64 comparisons in the milestone tests use `atol = rtol = 1e-12`. This is
comfortably above accumulated roundoff from trigonometric functions and the
independent symmetric eigensolver while remaining small relative to the
order-one hopping and band-energy scales. Exact structural properties such as
matrix dimensions and symmetry are tested exactly.
