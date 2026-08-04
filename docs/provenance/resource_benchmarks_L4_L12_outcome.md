# Resource Benchmarks Through `L=12`

## Classification

Each registered case at `L=4,6,8,10,12` is `RESOURCE_FEASIBLE`.

The largest resource-feasible size demonstrated within the registered ladder
is `L=12`. This is not a measured resource ceiling because larger sizes were
not part of this execution stage.

The benchmarks were executed from branch `milestone6-resource-benchmarks` at
Git commit `077209629c1de6851dc441816eabf1a0eaca33a2`. They were run in
strictly increasing size order. No result was rerun, replaced, or excluded,
and no simulation artifact was changed while creating this checkpoint.

## Registered Scope And Gate

All five cases use `beta=4.0`, `U=2.0`, `dtau=0.1`, `Ltau=40`,
`mu=-2.5202`, and stabilization period 5 on one CPU/BLAS thread with GPU
visibility disabled. The chemical potential is held fixed only to make timing
comparisons reproducible; these runs make no filling or interacting van-Hove
claim.

Each benchmark uses its registered seed, one warm-up sweep, and ten timed
complete sweeps. Its internal ten-sweep elapsed time is projected linearly to
two 3600-sweep production chains. `RESOURCE_FEASIBLE` requires:

- projected two-chain time at most 86400 seconds;
- peak resident memory at most 5905580032 bytes (5.5 GiB);
- zero swap use and no out-of-memory or numerical failure; and
- maximum Green-function recomputation correction at most `1e-8`.

The projection deliberately excludes Julia startup, state construction, and
the warm-up sweep. GNU `time -v` supplies the peak resident-memory and process
exit records, not the projected timing. Its zero `Swaps` count proves no swap
event occurred; no nonzero count is relabeled as a byte measurement.

## Preflight And Artifact Validation

Immediately before every interacting benchmark, the deterministic
`L=4`, `beta=4`, `U=0`, `mu=-1.897` density-and-energy oracle passed. The
largest absolute error in every preflight was
`1.06581410364015e-14`, below its `1e-10` tolerance.

Each immutable result directory contains `u0_preflight.log`, `benchmark.log`,
and `time.txt`, all covered by the
[SHA-256 manifest](resource_benchmarks_L4_L12_artifacts.sha256). Every
benchmark log records the exact registered single-case selector, physical and
workload parameters, finite timing metrics, an observed attempted-update count
consistent with that workload, and `outputs_written=false` from the benchmark
script. Every time record
reports exit status zero and zero swaps. The resource result root, five case
directories, and 15 files add 21 descendants to the previous 55; the complete
result tree contains 76 descendants after this stage.

## Results

Peak RSS bytes are derived exactly as the recorded GNU-time kilobyte value
multiplied by 1024. The gate was applied separately to each row.

| `L` | Seed | Timed 10 sweeps (s) | Projected two-chain time (s) | Projected hours | Wall time | Peak RSS (bytes) | Swaps | Maximum correction | Gate |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 4 | 6704 | 0.020131621 | 14.494767119999999 | 0.0040263242 | 0:07.63 | 567521280 | 0 | 9.241496456979803e-13 | `RESOURCE_FEASIBLE` |
| 6 | 6706 | 0.14875321 | 107.1023112 | 0.029750642 | 0:04.77 | 563392512 | 0 | 1.365546564713327e-12 | `RESOURCE_FEASIBLE` |
| 8 | 6708 | 0.621469531 | 447.45806232 | 0.1242939062 | 0:05.16 | 567803904 | 0 | 3.1188385207769898e-12 | `RESOURCE_FEASIBLE` |
| 10 | 6710 | 2.122549183 | 1528.2354117599998 | 0.42450983659999997 | 0:09.90 | 565915648 | 0 | 6.491166545247173e-12 | `RESOURCE_FEASIBLE` |
| 12 | 6712 | 7.474528187 | 5381.66029464 | 1.4949056374 | 0:15.91 | 619061248 | 0 | 5.881291981202352e-12 | `RESOURCE_FEASIBLE` |

All five attempts have finite metrics, no OOM evidence, no numerical failure,
and peak RSS far below the registered limit. Their projected two-chain times
also remain below 24 hours.

## Interpretation And Boundary

The registered resource-only evidence establishes that `L=12` is feasible for
the projected workload on this laptop under the stated `beta=4`, `U=2`,
single-thread CPU conditions. It does not establish that `L=12` is the largest
possible size. Conditional `L=14` and `L=16` benchmarks require a separate
prospective registration before execution.

This resource classification is distinct from scientific feasibility. The
initial `L=4`, `beta=4` fixed-density production pair is
`NOT_ADJUDICABLE` because both density confidence intervals missed the target.
Consequently no `L>=6` fixed-density calibration or sign chain, no correlation
production, and no spectral production is authorized by this checkpoint.
