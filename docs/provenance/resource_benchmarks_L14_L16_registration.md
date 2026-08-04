# Conditional `L=14` And `L=16` Resource Registration

## Registration Basis

The registered resource-only cases through `L=12` all returned
`RESOURCE_FEASIBLE`, as recorded in
`docs/provenance/resource_benchmarks_L4_L12_outcome.md` and merged to GitHub
`main` as commit `ddec4ba101d935ae52a8c6626ea158449cf07bc1`.

That result satisfies the prerequisite for this prospective extension. This
record registers only the two resource timing cases below and authorizes real
execution of `L=14` first. `L=16` is prospectively reserved but code-gated
until a later reviewed activation. Registration does not imply that either
case completed or passed.

## Authorized Cases

| Case | `L` | Seed | `beta` | `U` | `dtau` | `Ltau` | `mu` | Warm-up | Timed | Recorded state |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| `resource_L14_beta4` | 14 | 6714 | 4.0 | 2.0 | 0.1 | 40 | -2.5202 | 1 | 10 | completed, feasible |
| `resource_L16_beta4` | 16 | 6716 | 4.0 | 2.0 | 0.1 | 40 | -2.5202 | 1 | 10 | blocked pending activation |

Both cases must be executed with stabilization period 5, one Julia thread, one
BLAS thread, and no visible GPU. The protected execution handoff must set
`JULIA_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, and an empty
`CUDA_VISIBLE_DEVICES` value. The fixed chemical potential is a reproducible
timing input, not a fixed-density or interacting van-Hove claim.

The exact immutable output directories are:

- `results/dqmc_resource_benchmarks/resource_L14_beta4_seed6714`;
- `results/dqmc_resource_benchmarks/resource_L16_beta4_seed6716`.

Each directory must be absent before its attempt and must retain the
deterministic `U=0` preflight log, benchmark log, and GNU `time -v` record.
Existing directories or files may not be overwritten, replaced, or silently
rerun.

## Sequential Gate

Run `L=14` first. The runner registers `L=16` with
`execution_authorized=false` and must reject a real `L=16` invocation. If the
complete `L=14` record returns `RESOURCE_FEASIBLE`, first create and merge its
checksummed outcome checkpoint, then use a separate reviewed code-and-document
commit to set the `L=16` execution flag to true. Stop on an infeasible,
invalid, interrupted, or otherwise not-adjudicable `L=14` result; do not skip
to `L=16`.

Immediately before each interacting attempt, the deterministic
`L=4`, `beta=4`, `U=0`, `mu=-1.897` density-and-energy oracle must pass its
`1e-10` tolerance. Each benchmark must use the exact `--case` selector with
`--run` and `--require-single`, after a dry run confirms one selected case and
zero states, sweeps, and outputs.

The protected handoff, rather than the benchmark script, creates the fresh
canonical output directory and captures the external logs. It must verify the
directory is absent, the tracked worktree and index are clean, no Julia process
or open result handle exists, and all earlier manifests still pass before each
attempt. The benchmark script intentionally writes no files.

The internal ten-sweep time is projected linearly to two 3600-sweep chains.
The resource gate is unchanged:

- projected two-chain time at most 86400 seconds;
- peak resident memory at most 5905580032 bytes (5.5 GiB);
- zero swap use and no out-of-memory or numerical failure; and
- maximum Green-function recomputation correction at most `1e-8`.

GNU-time peak RSS in kilobytes is converted to bytes by multiplying by 1024.
A zero `Swaps` count establishes zero swap use; a nonzero count must not be
misrepresented as a byte value.

## Recorded `L=14` Outcome

The registered `L=14` attempt completed from Git commit
`6c77ee781dc56e5a33ec2b4fc126483a32a2b2e6` and returned
`RESOURCE_FEASIBLE`. Its projected two-chain time is
`14424.377793839998` seconds, peak RSS is `558632960` bytes, the GNU-time
swap count is zero, and the maximum Green-function correction is
`1.0183698177269884e-11`.

The full metrics and interpretation are recorded in the
[L=14 outcome](resource_benchmark_L14_outcome.md) and its
[SHA-256 manifest](resource_benchmark_L14_artifacts.sha256). This completion
satisfies the registered resource prerequisite for a later `L=16` activation,
but the current runner still has `execution_authorized=false`; `L=16` must not
run until the separate reviewed activation commit is merged.

## Boundary

These are short resource diagnostics without equilibration or physics
measurement. They make no density, sign, order-parameter, spectral, DOS, or
interacting van-Hove claim. They do not override the
`NOT_ADJUDICABLE` fixed-density `L=4`, `beta=4` pair and do not authorize any
scientific size chain or correlation production.

No `L>=18` resource case is registered. A completed `L=16` pass would establish
resource feasibility through the largest size in this registered extension,
not a ceiling beyond it.
