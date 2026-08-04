# Finite-`U` Response Campaign Registration

## Registration Basis

This protocol responds to the requested finite-interaction study of the upper
van-Hove region. It is registered from the reviewed Milestone 6 baseline merged
as Git commit `68a5a0cd22774164980e12f5b1e129716df29fe0`.

The operational van-Hove target is the noninteracting upper-peak filling

`n_target=1.46021875`.

This is a fixed-density convention. It is not a claim that the interacting
spectral function has a van-Hove singularity at the Fermi level. That question
requires the separately validated spectral program described below.

Unless a later registration explicitly changes them, all interacting chains in
this program use:

- `U=2.0`;
- `dtau=0.1`;
- `stabilization_period=5`;
- CPU execution with one BLAS thread; and
- ignored, immutable result directories with no overwrite or replacement.

The confirmed value `mu=-2.5202` applies only to `L=4`, `beta=4.0`. Every new
`(L,beta)` requires its own calibration before a fixed-density production
claim.

## Authorized Initial Batch

Only the work listed in this section is authorized by this registration.

### `L=4`, `beta=4` Fixed-Density Sign Chains

| Case | Seed | `mu` | Burn-in | Measurements | Between | Block |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `fixed_density_mu_m2p5202_L4_beta4` | 6311 | -2.5202 | 400 | 1600 | 2 | 100 |
| `fixed_density_mu_m2p5202_L4_beta4` | 6312 | -2.5202 | 400 | 1600 | 2 | 100 |

Both seeds must be run regardless of the first result. Each must complete its
exact workload, satisfy `abs(n-n_target) <= 2*SE_n`, and receive strict-pilot
classification `PASS`. A missing, invalid, failed, or inconclusive chain may
not be replaced or automatically rerun.

#### Recorded Outcome

Both registered chains completed and independently received strict-pilot
classification `PASS`, with average sign `1.0` and sign standard error `0.0`.
Neither production density satisfied the separately registered
`abs(n-n_target) <= 2*SE_n` condition. The exact two-chain gate therefore
returned `NOT_ADJUDICABLE`, `advance=false`, with reasons
`chain_1_density_not_target_compatible` and
`chain_2_density_not_target_compatible`.

This does not establish a fixed-density sign-feasibility point, a lowest
temperature, a sign limit, or an interacting van-Hove location. The chains are
immutable and may not be replaced or rerun. The full calculation is recorded
in the [pair outcome](fixed_density_mu_m2p5202_L4_beta4_pair_outcome.md) and
its [SHA-256 artifact manifest](fixed_density_mu_m2p5202_L4_beta4_pair_artifacts.sha256).

### `L=4`, `beta=6` Calibration Discovery

This scan was registered in advance but is execution-gated on both registered
`L=4`, `beta=4` chains passing every strict and density criterion. That
prerequisite was not met, so the registered `beta=6` discovery scan remains
blocked and must not be executed under this registration.

| Case | Seed | `mu` | Burn-in | Measurements | Between | Block |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `calibrate_L4_beta6_mu_m2p5802` | 6511 | -2.5802 | 200 | 800 | 2 | 100 |
| `calibrate_L4_beta6_mu_m2p5202` | 6512 | -2.5202 | 200 | 800 | 2 | 100 |
| `calibrate_L4_beta6_mu_m2p4602` | 6513 | -2.4602 | 200 | 800 | 2 | 100 |

All three points must complete. Their estimates must be finite, have usable
full and leave-one-out sign denominators, valid jackknife errors, and monotone
density as a function of `mu`. The target must either be compatible with a
sampled point at two standard errors or be bracketed by adjacent points at the
same precision.

If interpolation is required, use the adjacent bracketing means without
pre-rounding,

`mu_star = mu_lo + (n_target-n_lo)*(mu_hi-mu_lo)/(n_hi-n_lo)`,

then round once to four decimal places. The resulting `beta=6` production
chains require a later code-and-documentation registration with fresh seeds.
No bracket extension is automatic.

### Resource-Only Scaling Benchmarks

The cases `resource_L4_beta4`, `resource_L6_beta4`,
`resource_L8_beta4`, `resource_L10_beta4`, and `resource_L12_beta4` use
`beta=4.0`, `U=2.0`, `dtau=0.1`, and `mu=-2.5202`. Each benchmark performs a
warm-up followed by ten complete sweeps. It makes no filling, sign, or physics
claim.

Run benchmarks in increasing `L`. Capture `/usr/bin/time -v` output externally.
A lattice is resource-feasible only if:

- twice the projected time for a 3600-sweep production chain is at most 24
  hours;
- peak resident memory is at most 5.5 GiB;
- there is no swap, out-of-memory event, or numerical failure; and
- the maximum Green-function recomputation correction is at most `1e-8`.

Stop at the first resource failure. `L=14` and `L=16` benchmarks require a later
registration and are permitted only after `L=12` passes.

#### Recorded Resource Outcome

All five registered benchmarks completed in increasing `L` and independently
returned `RESOURCE_FEASIBLE`. The largest resource-feasible size demonstrated
within this registration is `L=12`. The result does not establish a resource
ceiling beyond the tested range: `L=14` and `L=16` were not authorized or run
by this registration. Full metrics and artifact validation are recorded in the
[resource outcome](resource_benchmarks_L4_L12_outcome.md) and its
[SHA-256 manifest](resource_benchmarks_L4_L12_artifacts.sha256).

These benchmarks remain short, unequilibrated timing diagnostics. They do not
measure or qualify density, sign, altermagnetic correlations, spectra, or an
interacting van-Hove location, and they do not override the initial
fixed-density pair's `NOT_ADJUDICABLE` stopping gate.

## Decision Rules For Later Batches

The reserved `L=4` temperature ladder is `beta=4,6,8,10,12,16`, corresponding
to `T=0.25,1/6,0.125,0.1,1/12,0.0625`. The reserved two-seed production pairs
are `(6311,6312)`, `(6313,6314)`, `(6315,6316)`, `(6317,6318)`, `(6319,6320)`,
and `(6321,6322)`. The first pair completed with the recorded
`NOT_ADJUDICABLE` outcome; no later production pair is authorized.

Advance to the next temperature only when both chains at the current
temperature pass all strict, density, denominator, sampling, and numerical
gates. A `DEFINITIVE_SIGN_FAIL` makes the preceding passed temperature the
lowest demonstrated reachable temperature. An `INCONCLUSIVE`, invalid, or
missing result stops the ladder as not adjudicable and does not establish a
sign limit.

The reserved `beta=4` size ladder is `L=6,8,10,12`, followed conditionally by
`L=14,16`, with production seed pairs `(6331,6332)`, `(6351,6352)`,
`(6371,6372)`, `(6391,6392)`, `(6411,6412)`, and `(6431,6432)`. None is
authorized until its resource benchmark and its separate fixed-density
calibration pass. Resource feasibility and sign-qualified scientific
feasibility must be reported separately.

## Correlation And Spectral Boundary

The local operator corresponding to Eq. 4 of
[Giuli, Mejuto-Zaera, and Capone](https://arxiv.org/abs/2410.00909),
Phys. Rev. B 111, L020401 (2025), is

`O_i = ((n_ix_up-n_ix_down) - (n_iy_up-n_iy_down))/2`.

The direct finite-lattice estimator is
`Delta_alm = sum(O_i)/L^2`. Because an unbiased finite system restores the
spin/orbital symmetry, the physics diagnostic is accompanied by

- `S_alm(0) = sum(<O_i O_j> for i,j)/L^2`; and
- `m_alm_squared = S_alm(0)/L^2`.

The two-point quantity must use the full equal-time Wick contraction, including
same-mode contact terms. A squared conditional one-point estimate is not an
acceptable substitute. Implementation and correctness tests are authorized;
production measurement is allowed only for a parameter point whose two strict
sign chains pass.

The implementation entry points are `measure_altermagnetic_observables` and
`summarize_sign_reweighted_altermagnetism`. The latter uses one common sign
denominator and sign-weights every raw field exactly once.

The current code does not yet produce interacting `A(k,omega)` or DOS. The
authorized implementation boundary is stabilized unequal-time
`G_ab,sigma(k,tau)` with per-block sign numerators/denominators, uncertainty,
covariance, and exact `U=0` validation. The arrays use a documented
no-leading-minus convention and remain raw until sign reweighting.
Real-frequency continuation is a distinct post-processing milestone.
[SmoQyDEAC](https://arxiv.org/abs/2407.04568) is the intended covariance-aware
method, but no continuation dependency is added here and it may not be applied
to interacting data until the imaginary-time measurement, sign denominator,
reconstruction, positivity, and spectral-sum rules pass.

The validated imaginary-time entry points are
`measure_unequal_time_greens`, `fourier_transform_unequal_time_greens`, and
`summarize_sign_reweighted_unequal_time_channel`. A valid covariance summary
is necessary but is not sufficient authorization for continuation.

## Immutability And Provenance

The runner remains dry-run by default, requires explicit single-chain
selection for strict pilots, and refuses existing outputs. Results from every
attempt, including unfavorable or inconclusive attempts, remain immutable.
Each completed stage receives a documentation-only provenance checkpoint and
SHA-256 manifest before another dependent batch is registered.

For every finite-`U` response chain, real execution additionally requires the
exact phase, case, and seed selectors, `--require-single`, `--strict-pilot`,
and the canonical output directory
`results/dqmc_sign_feasibility/<case>_seed<seed>_strict_pilot`. The runner
forbids both `--overwrite` and `--no-write` for these cases and refuses the
attempt if that directory already exists.
