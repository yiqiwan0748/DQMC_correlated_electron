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

## Band Path And Momentum Grid

`band_path` samples a configurable piecewise-linear path. Its default is
`Γ = (0,0)`, `X = (π,0)`, `M = (π,π)`, and `Γ = (0,0)`. The cumulative path
coordinate is Euclidean momentum-space distance. Segment junctions appear only
once, and the result includes high-symmetry tick coordinates, labels, and
sample indices.

`sample_band_energies(Nk)` evaluates `band_energies` exactly once at each point
of the periodic grid

```math
k_j=-\pi+\frac{2\pi j}{N_k},\qquad j=0,\ldots,N_k-1.
```

Its two-row energy array contains lower and upper bands, with `kx` varying
fastest. The same `BandGridSamples` object is reused by the DOS, filling, and
van-Hove calculations; those analyses do not diagonalize the Hamiltonian again.

## Density Of States And Filling

`density_of_states(samples; bins, smoothing_width)` uses an `O(Nk^2)` energy
histogram. Optional smoothing is a discrete Gaussian convolution whose width
is its standard deviation in energy units; the kernel is truncated at four
standard deviations. The reported DOS is per spin and is normalized after
smoothing so that its midpoint-rule integral is two states per physical site.

`zero_temperature_filling` counts sampled energies at or below the chemical
potential and multiplies by the spin degeneracy:

```math
n(\mu)=\frac{2}{N_k^2}\sum_{k,b}\Theta(\mu-E_b(k)).
```

Thus `0 <= n <= 4`, half filling is `n = 2`, and the finite-grid result is a
monotone step function. A sorted spectrum should be reused for a filling curve
and van-Hove fillings.

## van Hove Mapping And Doping Convention

`estimate_van_hove_peaks` ranks internal local maxima of the smoothed DOS,
excludes five percent of the energy span at each band edge, requires separated
peaks, and refines each bin maximum with quadratic interpolation. Each peak
energy is interpreted as the noninteracting `mu_vH`. The recorded candidates
are

```math
\delta_a=n-2,\qquad
\delta_b=\frac{n-2}{2},\qquad
\delta_c=\frac{n-2}{4}.
```

For `Nk = 800`, 1600 energy bins, and Gaussian width `0.05`:

| peak | `mu_vH` | `n(mu_vH)` | `delta_c` | journal value |
| ---: | ---: | ---: | ---: | ---: |
| 1 | -3.8369423191 | 0.42288125 | -0.3942796875 | -0.396 |
| 2 | -1.8970820299 | 1.46021875 | -0.1349453125 | -0.145 |

Numerical RMS matching gives `0.00721` for `delta_c`, much less than for
`delta_a` or `delta_b`. The journal convention is therefore inferred to be
`delta = (n-2)/4`. This is an inference from the reported numbers, not an
independent definition from the journal text.

The convergence comparison at the same binning and smoothing is:

| `Nk` | lower `mu_vH` | lower filling | upper `mu_vH` | upper filling |
| ---: | ---: | ---: | ---: | ---: |
| 200 | -3.8373763215 | 0.423700 | -1.8969418111 | 1.461100 |
| 400 | -3.8371649831 | 0.422625 | -1.8972160003 | 1.460175 |
| 800 | -3.8369423191 | 0.42288125 | -1.8970820299 | 1.46021875 |

The absolute `Nk=200` to `400` shifts are `0.0002113` and `0.0002742` in
energy, and `0.001075` and `0.000925` in filling. From `Nk=400` to `800`, they
are `0.0002227` and `0.0001340` in energy, and `0.00025625` and `0.00004375`
in filling. Remaining uncertainty comes from histogram binning, smoothing,
interpolation, and finite-grid occupation steps. The second inferred `delta_c`
differs from the approximate journal value by about `0.010`.

## Reproducible Outputs

`scripts/generate_noninteracting_results.jl` writes the band path, primary-grid
DOS and filling curve, convergence summary, metadata, and two PNG figures to
`results/noninteracting/`. `metadata.txt` records all numerical and physical
conventions used by the outputs.

## Finite-Size Real-Space Kinetic Matrix

Milestone 3 represents one spin sector on an `L x L` periodic square lattice.
Coordinates are one-based, `x,y = 1,...,L`, and the orbital basis remains
`(d_xz,d_yz)`. The orbital-fast linear index is

```math
\mathrm{site}=x+L(y-1),\qquad
p=2(\mathrm{site}-1)+\mathrm{orbital}.
```

Thus the dense real-space matrix has dimension `Nsp = 2L^2`. Integer
coordinates are wrapped independently into `1:L` with periodic boundary
conditions. The public forward and inverse index functions validate their
coordinates, orbital, index, and lattice size before returning an `Int` index
or the tuple `(x,y,orbital)`.

For each source cell, the four stored positive-displacement blocks are

```math
T_x = \begin{pmatrix}-t_1&0\\0&-t_2\end{pmatrix},\quad
T_y = \begin{pmatrix}-t_2&0\\0&-t_1\end{pmatrix},
```

```math
T_+ = \begin{pmatrix}-t_3&+t_4\\+t_4&-t_3\end{pmatrix},\quad
T_- = \begin{pmatrix}-t_3&-t_4\\-t_4&-t_3\end{pmatrix},
```

at displacements `(+1,0)`, `(0,+1)`, `(+1,+1)`, and `(+1,-1)`, respectively.
The block is added from the source to the wrapped destination and its transpose
is added for the reverse hopping. Every entry is accumulated with `+=`. This
is physically necessary when distinct directed bonds alias to one matrix entry
on a small torus; in particular, it retains the doubled nearest-neighbor and
quadrupled diagonal contributions at `L=2`. The opposite signs of `T_+` and
`T_-` make the isolated interorbital `t4` matrix cancel at `L=2`, as required
because every allowed component there is zero or `pi`.

## Finite-Lattice Fourier Convention

The allowed components are

```math
k_n=\frac{2\pi n}{L},\qquad n=0,\ldots,L-1.
```

Using cell phase coordinates `r=(x-1,y-1)`, the projected block is

```math
H_{\alpha\beta}(k)=\frac{1}{L^2}\sum_{r,r'}
e^{-ik\cdot r}K_{r\alpha,r'\beta}e^{+ik\cdot r'}.
```

A common one-cell shift of all phase coordinates cancels. With the hopping
blocks above, this projection reproduces the validated `bloch_hamiltonian`
element by element, including
`epsilon_xy = -4*t4*sin(kx)*sin(ky)`.

The Milestone 3 public API is:

- `periodic_coordinate(coordinate, L)` wraps an integer into `1:L`;
- `single_particle_index(x, y, orbital, L)` and
  `single_particle_coordinates(p, L)` implement the checked index map;
- `realspace_kinetic_matrix(L, parameters=DEFAULT_PARAMETERS)` builds the
  dense matrix;
- `allowed_momenta(L)` returns the finite periodic grid;
- `fourier_block(K, L, kx, ky)` performs the normalized projection; and
- `realspace_band_energies(K)` returns the complete sorted Hermitian spectrum.

No chemical-potential shift is included. Sparse storage, checkerboard
decomposition, matrix exponentials, Green functions, Hubbard-Stratonovich
fields, and all DQMC machinery remain outside this milestone.

`scripts/validate_realspace_kinetic.jl` checks `L=2,3,4,5,6`, writes the
symmetry, Fourier-block, and complete-spectrum errors to
`results/realspace_kinetic/validation_summary.csv`, and records these
conventions in the adjacent `metadata.txt`. The Float64 acceptance threshold
for Fourier and spectral mismatches is `1e-10`; symmetry is checked at a
machine-precision-scale absolute tolerance. The default finite-grid spectrum
must stay inside the validated global interval `[-4.9,8.9]`.

## Finite-Temperature Free Fermions

Milestone 4 shifts the validated one-spin real-space hopping matrix by the
chemical potential,

```math
K_\mu=K-\mu I.
```

Thus positive `mu` lowers every one-particle energy. If
`K_mu = Q*diagm(lambda)*Q^T`, the DQMC equal-time Green matrix is

```math
G_{ij}=\langle c_i c_j^\dagger\rangle
=Q\,\mathrm{diag}\!\left(\frac{1}{1+e^{-\beta\lambda_a}}\right)Q^T.
```

The one-particle density matrix has the deliberately different operator and
index order

```math
C_{ij}=\langle c_j^\dagger c_i\rangle
=Q\,\mathrm{diag}\!\left(\frac{1}{1+e^{\beta\lambda_a}}\right)Q^T.
```

Consequently `C = I-G` for this real-symmetric model. Both matrices describe
one spin sector. Their diagonal entries are one-spin occupations; the public
density summaries multiply their traces by the spin degeneracy two:

```math
n_x=\frac{2}{L^2}\sum_r C_{r x,r x},\qquad
n_y=\frac{2}{L^2}\sum_r C_{r y,r y},\qquad
n=\frac{2\,\mathrm{tr}(C)}{L^2}=n_x+n_y.
```

Here `x` and `y` label orbitals `d_xz` and `d_yz`, not lattice
coordinates. The physical filling range is `0 <= n <= 4`, and `beta = 0`
gives `n = 2` for every finite chemical potential.

The Milestone 4 public API is:

- `shifted_kinetic_matrix(K, mu)` constructs `K_mu` without mutating `K`;
- `dense_slice_propagators(K_mu, dtau)` returns
  `B = exp(-dtau*K_mu)` and `B_inv = exp(+dtau*K_mu)`;
- `free_equal_time_greens(K_mu, beta)` returns `G`;
- `free_density_matrix(K_mu, beta)` returns `C`;
- `orbital_resolved_density(C, L)` returns `(n_x, n_y)`, including spin;
- `total_filling_per_physical_site(C, L)` returns `n`, including spin; and
- `free_grand_log_partition(K_mu, beta; spin_degeneracy=1)` returns the
  one-spin grand `log Z` by default and the physical two-spin value when the
  keyword is two.

All primary matrix functions use a symmetric eigendecomposition. The Green
matrix is never constructed through a direct inverse. Positive- and
negative-argument logistic branches evaluate Fermi factors using only
`exp(-abs(x))`, avoiding overflow and retaining small occupations without
subtractive cancellation. The grand partition function uses

```math
\log Z_\sigma=\sum_a\log(1+e^{-\beta\lambda_a})
```

through the stable softplus identity
`max(x,0) + log1p(exp(-abs(x)))`. Dense slice exponentials are spectral and
raise `OverflowError` if both forward and inverse exponentials cannot be
represented in the working floating-point type.

Every one-spin free determinant factor is
`1+exp(-beta*lambda_a) > 0`; therefore the two-spin free-fermion weight is
positive. This statement is only about the noninteracting model. It does not
imply that the interacting doped model is sign-problem free.

The validation grid uses `L=2,4,6`, `beta=1,2,4,8`,
`mu=-5.5,-3.837,-1.897,0,9.5`, and `dtau=0.1`. Tests additionally cover
all `L=2:6` in real-space versus finite-grid momentum comparisons, limiting
fillings, monotonicity, both Milestone 2 van Hove chemical potentials, and
safe direct-definition checks. Float64 tolerances are separated by purpose:
`2e-12` for symmetry and Green identities, `2e-11` for slice inverses,
`5e-11` for accumulated propagators, and `1e-11` for density comparisons.
The accumulated-propagator error is relative to
`max(1,maxabs(exp(-beta*K_mu)))`.

`scripts/validate_free_fermion.jl` writes the full grid to
`results/free_fermion/validation_summary.csv` and all chemical-potential,
Green-function, density, spin, normalization, stability, partition-function,
and error conventions to the adjacent `metadata.txt`.
