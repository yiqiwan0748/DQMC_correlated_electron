# Fixed-Density Calibration Confirmation at `mu=-2.5202`

## Classification

The registered three-seed result is
`CONFIRMED_AT_REGISTERED_PRECISION`.

This record adjudicates the fixed operational rule registered in
`docs/roadmap.md` at Git commit
`6c263caf3b646f2c958240edee9b69ddaea88525`. No chain was rerun, replaced, or
excluded, and no simulation output was changed while creating this checkpoint.

## Registered Scope

| Quantity | Registered value |
| --- | ---: |
| `L` | 4 |
| `beta` | 4.0 |
| `U` | 2.0 |
| `dtau` | 0.1 |
| `Ltau` | 40 |
| `mu` | -2.5202 |
| `n_target` | 1.46021875 |
| `burnin_sweeps` | 200 |
| `measurement_sweeps` | 800 |
| `sweeps_between_measurements` | 2 |
| `block_size` | 100 |
| `stabilization_period` | 5 |

Seed 6117 is the previously observed pilot. Seeds 6118 and 6119 are the two
prospectively registered independent confirmations.

## Artifact Validation

For each seed, the registered output directory contains exactly one
`sign_feasibility.csv` and one `metadata.txt` artifact covered by the
[SHA-256 manifest](fixed_density_mu_m2p5202_artifacts.sha256). Each CSV contains
one completed row. Its selected-chain and completed-diagnostics fields agree
with its metadata, including the physical parameters, workload, seed, target,
filling estimate, uncertainty, update counts, acceptance, sign diagnostics,
Green-function correction, and strict-pilot status.

All three rows report 800 measurements in eight blocks of length 100, a finite
filling estimate, a positive finite block-jackknife standard error, usable sign
denominators, valid numerical diagnostics, no failed strict-pilot gates, and
`strict_pilot_status=PASS`. Thus all three chains meet the registered validity
prerequisite for adjudication.

## Registered Calculation

The source estimates and standard errors below are the unrounded decimal
strings stored in the CSV artifacts. Derived values were evaluated with
high-precision decimal arithmetic without pre-rounding; displayed derived
values were rounded only after each registered condition was evaluated.

Individual compatibility requires
`abs(n_i - n_target) <= 2*s_i`.

| Seed | Role | `n_i` | `s_i` | `abs(n_i - n_target)` | `2*s_i` | Result |
| ---: | --- | ---: | ---: | ---: | ---: | --- |
| 6117 | previously observed | 1.4605219571679755 | 0.001094369279004781 | 0.0003032071679755 | 0.002188738558009562 | pass |
| 6118 | prospective confirmation | 1.4592898649408286 | 0.0016342652723925584 | 0.0009288850591714 | 0.0032685305447851168 | pass |
| 6119 | prospective confirmation | 1.4597148203412997 | 0.000961522441271973 | 0.0005039296587003 | 0.001923044882543946 | pass |

Inverse-variance pooling uses `w_i=1/s_i^2`,
`n_pool=sum(w_i*n_i)/sum(w_i)`, and
`SE_pool=1/sqrt(sum(w_i))`. Both registered pools require
`abs(n_pool - n_target) <= 2*SE_pool`.

| Seed set | `n_pool` | `SE_pool` | `abs(n_pool - n_target)` | `2*SE_pool` | Result |
| --- | ---: | ---: | ---: | ---: | --- |
| `S_confirm={6118,6119}` | 1.4596055452386897 | 0.0008287269978258003 | 0.0006132047613102906 | 0.0016574539956516005 | pass |
| `S_all={6117,6118,6119}` | 1.4599395349665640 | 0.0006606707522167680 | 0.0002792150334360302 | 0.0013213415044335360 | pass |

For all three seeds, the registered heterogeneity statistic is

`Q=0.4958851219864424 <= 5.991464547`.

All three chains are valid, every individual condition passes, both pooled
conditions pass, and the heterogeneity condition passes. The third branch of
the preregistered decision rule therefore applies:

`CONFIRMED_AT_REGISTERED_PRECISION`.

## Interpretation

Only `S_confirm={6118,6119}` is fully prospective; seed 6117 was observed before
the confirmation rule was registered. The pooling calculation is the specified
fixed-effect calculation. Passing `Q` does not prove homogeneity, and the
combined `+/-2SE` rules are not one exact 95% hypothesis test.

The classification applies only to the registered lattice size, temperature,
interaction, time step, chemical potential, workload, and seeds. It does not
authorize another chain, a changed chemical potential, extrapolation to other
parameters, or correlation measurements.
