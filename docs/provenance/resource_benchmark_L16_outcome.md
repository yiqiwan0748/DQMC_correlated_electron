# `L=16` Resource Benchmark Outcome

## Classification

The registered `resource_L16_beta4` case is `RESOURCE_FEASIBLE`.

The benchmark was executed from branch `milestone6-resource-l16-execution` at
Git commit `5ca19b43b20c52e8a91a256ef17b6fc90b5e7243`. No result was rerun,
replaced, or excluded, and no simulation artifact was changed while creating
this checkpoint.

## Registered Scope

| Quantity | Registered value |
| --- | ---: |
| `L` | 16 |
| seed | 6716 |
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
[SHA-256 manifest](resource_benchmark_L16_artifacts.sha256). Immediately before
the interacting benchmark, the deterministic `U=0` oracle passed with maximum
absolute error `1.06581410364015e-14`, below its `1e-10` tolerance.

The benchmark and GNU-time records agree on the exact registered selector and
successful execution. The time record reports exit status zero and zero
swaps. The completed directory and three files add four descendants to the
previous 80; the complete result tree contains 84 descendants.

## Result And Gate

| Metric | Stored or derived value | Limit |
| --- | ---: | ---: |
| timed ten-sweep time | 41.617442088 s | -- |
| seconds per complete sweep | 4.1617442088 | -- |
| projected two-chain time | 29964.55830336 s | 86400 s |
| projected two-chain time | 8.3234884176 h | 24 h |
| wall time | 0:54.56 | -- |
| peak RSS | 630672 KiB | -- |
| peak RSS | 645808128 bytes | 5905580032 bytes |
| swap count | 0 | 0 |
| maximum Green-function correction | 1.606504165807543e-11 | 1e-8 |
| attempted updates | 204800 | -- |
| accepted updates | 184808 | -- |
| acceptance rate | 0.9023828125 | -- |

Peak RSS bytes are the recorded GNU-time kilobyte value multiplied by 1024.
Finite benchmark output, exit status zero, and the absence of OOM or numerical
failure evidence support the gate record. Every registered threshold passes,
so the deterministic gate returns `RESOURCE_FEASIBLE` with no reasons.

## Interpretation And Boundary

This result establishes resource feasibility through `L=16`, the largest size
tested in the registered ladder, for the projected workload on this laptop
under the stated `beta=4`, `U=2`, CPU-only conditions. It is not a resource
ceiling and is not scientific fixed-density feasibility.

No `L>=18` case, scientific calibration, production sign chain, correlation
measurement, spectrum, or DOS is authorized. The fixed-density `L=4`,
`beta=4` production pair remains `NOT_ADJUDICABLE`, so no sign-qualified
lattice size or lowest qualified temperature has been established.
