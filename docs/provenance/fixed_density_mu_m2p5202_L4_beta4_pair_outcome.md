# Fixed-Density `L=4`, `beta=4` Production-Pair Outcome

## Classification

The registered seed-6311 and seed-6312 pair is

`NOT_ADJUDICABLE`.

Both chains completed their exact workloads and independently received
strict-pilot classification `PASS`. Neither chain satisfied the separately
registered fixed-density compatibility condition. The pair therefore does not
advance the temperature ladder.

This record adjudicates the initial pair registered in
`docs/provenance/finite_u_response_registration.md`. No chain was rerun,
replaced, or excluded, and no simulation artifact was changed while creating
this checkpoint. Both attempts record execution from branch
`milestone6-finite-u-execution` at Git commit
`857c7faa59e1529da55c86047a31bbbe6f217428`.

## Registered Scope

| Quantity | Registered value |
| --- | ---: |
| `L` | 4 |
| `beta` | 4.0 |
| `T` | 0.25 |
| `U` | 2.0 |
| `dtau` | 0.1 |
| `Ltau` | 40 |
| `mu` | -2.5202 |
| `n_target` | 1.46021875 |
| `burnin_sweeps` | 400 |
| `measurement_sweeps` | 1600 |
| `sweeps_between_measurements` | 2 |
| `block_size` | 100 |
| `stabilization_period` | 5 |
| seeds | 6311, 6312 |

The target filling is the registered operational proxy for the upper
noninteracting van-Hove filling. It is not a claim about the location of an
interacting spectral singularity.

## Artifact Validation

Each canonical result directory contains `metadata.txt`,
`sign_feasibility.csv`, `runner.log`, `time.txt`, and `attempt_status.txt`, all
covered by the [SHA-256 manifest](fixed_density_mu_m2p5202_L4_beta4_pair_artifacts.sha256).
The CSV and metadata values agree. Both attempts exited with status zero and
record the exact registered selectors and workloads. Each chain contains 1600
measurements in 16 blocks of length 100, usable full and leave-one-out sign
denominators, a valid positive finite density jackknife error, and a maximum
Green-function recomputation correction below `1e-8`.

The two immutable result directories add 12 descendants to the previous 43;
the result tree contains 55 descendants after this stage.

## Registered Calculation

The density estimates and errors below are the unrounded decimal strings
stored in the result artifacts. Derived values were evaluated from those
strings with high-precision decimal arithmetic without pre-rounding.
Individual compatibility requires

`abs(n_i-n_target) <= 2*SE_i`.

| Seed | `n_i` | `SE_i` | `abs(n_i-n_target)` | `2*SE_i` | Density result |
| ---: | ---: | ---: | ---: | ---: | --- |
| 6311 | 1.4585174067362547 | 0.0006899852993155704 | 0.0017013432637453 | 0.0013799705986311408 | fail |
| 6312 | 1.4586825058255015 | 0.0006367611275137972 | 0.0015362441744985 | 0.0012735222550275944 | fail |

The stored strict-pilot diagnostics are:

| Seed | Average sign | Sign `SE` | Sign lower `2SE` | Internal runtime (s) | Wall time | Peak RSS (KiB) | Maximum correction | Strict status |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 6311 | 1.0 | 0.0 | 1.0 | 7.135014961 | 0:15.40 | 635280 | 2.0728307958961523e-11 | `PASS` |
| 6312 | 1.0 | 0.0 | 1.0 | 7.225680869 | 0:15.65 | 634852 | 5.892175636290631e-12 | `PASS` |

Both `/usr/bin/time -v` records report zero swaps and exit status zero. The
registered conservative pair sign diagnostic is

`min(average_sign_i - 2*SE_sign_i) = 1.0`.

With the gate input ordered as seeds 6311 and 6312, the exact pair gate
returned:

```text
status = NOT_ADJUDICABLE
advance = false
reasons = (chain_1_density_not_target_compatible,
           chain_2_density_not_target_compatible)
density_compatible_seeds = ()
conservative_average_sign_lower_2se = 1.0
```

## Interpretation And Stopping Rule

The sign and numerical diagnostics passed for both sampled chains, but those
facts do not override the fixed-density prerequisite. At the longer production
workload, both density confidence intervals miss the registered target. This
also means that the per-chain CSV field `correlation_measurement_ready=true`
does not authorize the campaign-level correlation measurement: the runner does
not apply the separate target-density pair gate. This does not revoke the
earlier `CONFIRMED_AT_REGISTERED_PRECISION` calibration
classification, which applies to its own registered 200/800 workload and seed
set; it means the new 400/1600 production pair does not qualify at its greater
precision.

No fixed-density sign-feasibility point or lowest reachable temperature is
established by this pair. This is not a sign failure and does not establish a
sign-problem limit. The conditionally registered `L=4`, `beta=6` discovery scan
must not run, and no correlation measurement is authorized. Continuing the
scientific ladder would require a separate, prospective `L=4`, `beta=4`
recalibration-extension registration with fresh seeds and immutable outputs;
this checkpoint does not register that extension.

The resource-only size benchmarks remain separately authorized. They make no
filling, sign, fixed-density, or interacting-physics claim.
