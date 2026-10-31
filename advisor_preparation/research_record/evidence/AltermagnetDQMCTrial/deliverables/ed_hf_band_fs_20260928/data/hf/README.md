# Matched-density Hartree–Fock comparison for the ISAAC ED update

This is a **separate Hartree–Fock calculation**, at the same Hamiltonian parameters and total density as the eight-cell, twelve-electron ED result. It is not an interpolation or fit of the ED spectrum. The plotted Hartree–Fock states do not have the ED state's fixed spin populations.

The run finished on 28 September 2026 at approximately 12:16:51 p.m. EDT. The scientific calculation and export took 30.817 seconds, with no scheduler queue. Its independently checked quantities passed 109 numerical checks. This is a precision statement about the stated calculation, not proof of the correctness of Hartree–Fock for the many-body model.

## Hamiltonian and conventions

Use orbital basis `(d_xz,d_yz)`, lattice spacing one, momentum in radians, and energy in units of `|t1|=1`. The hopping tuple is `t1=-1, t2=-1.75, t3=-0.85, t4=-0.65`. Define

```
c(k) = -(t1+t2)(cos kx+cos ky)-4 t3 cos kx cos ky
z(k) = -(t1-t2)(cos kx-cos ky)
x(k) = -4 t4 sin kx sin ky
H0(k) = c(k) I + z(k) tau_z + x(k) tau_x.
H = T + U sum_(r,a) (n_(ra,up)-1/2)(n_(ra,down)-1/2), U=20.
```

There is no Hund interaction, interorbital density interaction, spin-orbit coupling, or applied Zeeman term. The Hartree–Fock fields below represent variational order in a chosen ansatz, not a new term added to the physical model.

Hartree–Fock uses an infinite periodic square lattice represented by a Brillouin-zone quadrature. Total density is `n=1.5` electrons per unit cell, equal to `N=12` on eight ED cells. The final Hartree–Fock result uses `T=0.0125`; ED used a zero-temperature initial state. No zero-temperature extrapolation of the Hartree–Fock result is claimed.

## Restricted states and self-consistency

Both collinear uniform ansätze impose equal total orbital densities `n_x=n_y=n/2=0.75`. Write `m_a=n_(a,up)-n_(a,down)`. The AM line has `m_x=-m_y=m`; the FM line has `m_x=m_y=m`. Their source matrices are

```
Hs(k) = H0(k) - s h Gamma,  s=+1 for up and -1 for down
Gamma_AM = tau_z, Gamma_FM = I,  h=U m/2.
```

For each trial `h`, solve the Fermi occupation equation for `mu_source` so the spin-and-band sum of `f((e_source-mu_source)/T)` per cell is `1.5`. The common centered-interaction Hartree shift is `U(n-2)/4=-2.5`, so `mu_centered=mu_source-2.5`. Both the band energies and chemical potential get that same shift; all exported energies are `e_source-mu_source = e_centered-mu_centered`.

The fixed-density variational source functional per cell is

```
F(h) = Omega_source + mu_source*n + 2h^2/U + U(n-2)^2/8
dF/dh = -2m(h) + 4h/U.
```

Each branch is scanned at 128 positive-field intervals, starting at `h=1e-5` and extending to `U*n/4+0.05=7.55`. A self-consistent positive field cannot exceed `U*n/4=7.5`. Sign-changing stationary roots are refined by safeguarded secant/bisection; `h=0` is included explicitly. This sampled scan cannot rigorously exclude arbitrarily narrow pairs of roots between samples. Negative fields are degenerate spin-reversed partners.

The AM and FM states selected for the plots are the positive-field minima along their respective one-dimensional lines. PM is a zero-field stationary reference, unstable along both tested magnetic directions. The calculation does not search mixed FM/AM fields, orbital charge order, finite-momentum ordering, or phase separation; it does not establish a global Hartree–Fock ground state.

## Final values and a numerical checkpoint

At `Nk=384`, `T=0.0125`:

| Branch | h | mu_source | mu_centered | F per cell | up density | down density |
|---|---:|---:|---:|---:|---:|---:|
| AM | 7.378598021947 | -5.694909296117 | -8.194909296117 | -7.676642768677 | 0.75 | 0.75 |
| FM | 7.499999999972 | -5.196133026189 | -7.696133026189 | -7.638249439992 | 1.50 | approximately 0 |
| PM | 0 | -1.866913032920 | -4.366913032920 | -4.276642562362 | 0.75 | 0.75 |

For AM, `(n_xup,n_xdown)=(0.743929901097,0.006070098903)` and the `y` orbital reverses these populations. Thus `m=0.737859802195`, and `U*m/2=7.378598021947` reproduces the self-consistent field. AM has lower free energy than FM by approximately `0.0383933` per cell within the tested lines.

The selected ED state has `(Nup,Ndown)=(10,2)`, or `(1.25,0.25)` per cell, and finite total spin `S=4`. It is neither of these Hartree–Fock states. The comparison matches `H`, `U`, hoppings and total density; it does not match the initial many-body state, net magnetization, system size or temperature exactly.

At `X=(pi,0)`, `c=-3.4`, `z=1.5`, and `x=0`. On the AM branch,

```
e_source(up,lower) = -3.4 - |1.5-h| = -9.278598021947
e_source(up,upper) = -3.4 + |1.5-h| =  2.478598021947
e_source(down,lower) = -3.4 - |1.5+h| = -12.278598021947
e_source(down,upper) = -3.4 + |1.5+h| =  5.478598021947.
```

Subtract `mu_source=-5.694909296117` to get the displayed energies. At `Y` the opposite-spin bands are exchanged by the ansatz symmetry. This property of the separate AM Hartree–Fock solution does not repair or replace a failed spin-rotation relation in a different ED state.

## Numerical controls and what the figures mean

The density quadrature is a midpoint periodic mesh in each momentum direction; it never duplicates the Brillouin-zone boundary. All stationary solutions were retained for `(Nk,T)=(192,0.025),(384,0.025),(384,0.0125)`. At `T=0.025`, increasing `Nk` from 192 to 384 changes AM `h` by `-2.5903e-5` and `mu_source` by `+0.0024900`. At `Nk=384`, lowering `T` from 0.025 to 0.0125 changes AM `h` by `-7.4038e-6` and `mu_source` by `+0.0005608`. These are measured controls, not rigorous error bars or a claim that every last displayed digit is mesh converged. The complete differences for both magnetic branches are in `mesh_temperature_controls.csv`.

The band path is `Gamma-X-M-Y-Gamma-M` with the actual accumulated Euclidean momentum distance in radians. The first four segments have length `pi`; the last diagonal has length `sqrt(2)*pi`. `hf_path_vertices.csv` gives exact tick locations. Each spin has two orbital-mixed eigenvalue bands, ordered lower/upper independently at each momentum. Bands alone have no linewidth and are not the multi-pole correlated ED spectral function. If displayed with a linewidth, that is an additional documented rendering choice.

`hf_fermi_grid.csv` evaluates the selected source bands on a separate 257-by-257 display grid, including both periodic endpoints. A Hartree–Fock Fermi contour is the zero contour of `e_source(k)-mu_source`. It is not inferred from Lorentzian tails or spectral-intensity thresholds. The final chemical potential is density matched at finite `T=0.0125`, so this is a finite-temperature Hartree–Fock band-crossing contour, not a separately solved strict-zero-temperature Fermi surface. Any filled/unfilled shading based on the sign of the energy indicates whether a band lies below/above the chemical potential; exact finite-temperature occupations are Fermi functions.

## Verification and reproduction

`verification.toml` records each check, its absolute error and tolerance. The checks include density and field residuals, the independent kinetic-plus-interaction-minus-entropy energy, a finite-difference derivative of the fixed-density free energy, explicit symmetric 2-by-2 eigensystems at generic and symmetry points, all four spin/orbital populations computed with explicit eigenvectors on the final 384-by-384 mesh, AM rotated opposite-spin band equality, and PM spin degeneracy. Numerical tolerances exceed solver/Float64 roundoff; the looser derivative tolerance includes finite-difference truncation. Mesh and temperature sensitivity are retained separately.

Reproduce from the repository root:

```
OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1 julia --project=. scripts/isaac_band_fs_hf_20260928.jl
```

Only Julia standard libraries are needed. The additive driver reads the definition prefix of `scripts/literature_transfer/hf_check.jl` before the unique `mkpath(OUTPUT_HF)` marker. Its historical execution block is never run and its historical result directory is not modified. Both calculation source SHA-256 hashes are in `hf_summary.toml`. `computed_source.jl` is the exact driver used for the numerical execution, and `helper_source.txt` preserves the historical helper. The current driver and its `driver_source.txt` snapshot have one bookkeeping fix after that execution: exclude `output_hashes.toml` from hashing itself on repeat runs. No scientific logic or data changed. `source_provenance.json` records both hashes and that change. `output_hashes.toml` hashes the saved output files, except itself.

Data mapping: `stationary_points.csv` contains every detected stationary solution; `field_scan.csv` contains the sampled derivatives; `populations.csv` contains final orbital/spin densities; `hf_path.csv` supplies band plots; `hf_path_vertices.csv` supplies path ticks; `hf_fermi_grid.csv` supplies Fermi contours; `hf_summary.toml` supplies the final parameters and explicit limitations.
