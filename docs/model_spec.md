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
