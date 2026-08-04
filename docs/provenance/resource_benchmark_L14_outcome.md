# `L=14` Resource Benchmark Outcome

## Classification

The registered `resource_L14_beta4` case is `RESOURCE_FEASIBLE`.

The benchmark was executed from branch `milestone6-resource-l14-execution` at
Git commit `6c77ee781dc56e5a33ec2b4fc126483a32a2b2e6`. No result was rerun,
replaced, or excluded, and no simulation artifact was changed while creating
this checkpoint.

## Registered Scope

| Quantity | Registered value |
| --- | ---: |
| `L` | 14 |
| seed | 6714 |
| `beta` | 4.0 |
| `U` | 2.0 |
| `dtau` | 0.1 |
| `Ltau` | 40 |
| `mu` | -2.5202 |
| stabilization period | 5 |
| warm-up sweeps | 1 |
| timed complete sweeps | 10 |
| projected sweeps per chain | 3600 |
| projected chains | 2 |

Execution used one Julia thread, one BLAS thread, and no visible GPU. The
fixed chemical potential is a reproducible timing input only; this benchmark
makes no density, sign, or interacting van-Hove claim.

## Artifact Validation

The canonical result directory contains `u0_preflight.log`, `benchmark.log`,
and `time.txt`, all covered by the
[SHA-256 manifest](resource_benchmark_L14_artifacts.sha256). Immediately before
the interacting benchmark, the deterministic `U=0` oracle passed with maximum
absolute error `1.06581410364015e-14`, below its `1e-10` tolerance.

The benchmark and GNU-time records agree on the exact registered selector and
successful execution. The time record reports exit status zero and zero
swaps. The completed directory and three files add four descendants to the
previous 76; the complete result tree contains 80 descendants.

## Result And Gate

| Metric | Stored or derived value | Limit |
| --- | ---: | ---: |
| timed ten-sweep time | 20.033858047 s | -- |
| seconds per complete sweep | 2.0033858046999997 | -- |
| projected two-chain time | 14424.377793839998 s | 86400 s |
| projected two-chain time | 4.006771609399999 h | 24 h |
| wall time | 0:30.13 | -- |
| peak RSS | 545540 KiB | -- |
| peak RSS | 558632960 bytes | 5905580032 bytes |
| swap count | 0 | 0 |
| maximum Green-function correction | 1.0183698177269884e-11 | 1e-8 |
| attempted updates | 156800 | -- |
| accepted updates | 141231 | -- |
| acceptance rate | 0.9007079081632653 | -- |

Peak RSS bytes are the recorded GNU-time kilobyte value multiplied by 1024.
Finite benchmark output, exit status zero, and the absence of OOM or numerical
failure evidence support the gate record. Every registered threshold passes,
so the deterministic gate returns `RESOURCE_FEASIBLE` with no reasons.

## Boundary And Next Gate

This result establishes resource feasibility at `L=14` for the projected
workload on this laptop under the stated `beta=4`, `U=2`, CPU-only conditions.
It is not scientific fixed-density feasibility and is not a resource ceiling.

The registered `L=16` case remains code-gated with
`execution_authorized=false`. This checkpoint must be reviewed and merged,
then a separate reviewed activation must set that flag to true before `L=16`
can run. No `L>=18` case, scientific calibration, production sign chain,
correlation measurement, or spectrum is authorized.
