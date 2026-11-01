# Matched Rényi node-loop timing

This record compares the same locally lazy two-replica workflow **with and without rank-one Grover caching**. The primary metric is the ratio of median elapsed node-loop times; update-kernel results remain separate.

**Post-primary cut-size exploration.** This campaign was designed after seeing the original half-chain results and does not replace that preregistered grid. Both exploratory cuts must remain visible in the combined public record. A full-system cut (nA=L) describes finite-temperature thermal Rényi purity, not pure-state bipartite entanglement. No primary headline or automatic CV candidate is selected from this exploration.

Completed exploratory case **L=64, nA=64**: **2.862×** node-loop speedup (9.5773 s direct / 3.34693 s optimized).

| L | nA | Status | Paired repeats | Direct median (s) | Optimized median (s) | Node-loop speedup | Timing-jitter 95% interval | Update-kernel speedup |
|---:|---:|---|---:|---:|---:|---:|---|---:|
| 64 | 64 | PASS | 7/7 | 9.5773 | 3.34693 | 2.8615× | [2.7094, 3.0071] | 2.8652× |

## Allocations

These are median cumulative allocated bytes inside each timed node loop, including allocation churn. They do not measure peak resident memory. The reduction ratio is direct median bytes divided by optimized median bytes.

| L | nA | Direct median allocated bytes | Optimized median allocated bytes | Allocation reduction |
|---:|---:|---:|---:|---:|
| 64 | 64 | 2.48971032e+10 | 40805472 | 610.14× |

## Scope and fairness

Fixed physics: nearest-neighbor periodic Hubbard chain, t=1, U=2, centered mu=0, beta=2, dtau=0.05, lambda=0.5, nA as declared in the table and stabilization period 5. Both paths use the corrected 50% independent idle decision before every local proposal. Both start from clones of identical warmed states and RNG streams; execution order alternates. The production protocol uses seven repeats of twelve two-replica rounds. The same forward/backward and endpoint ratio measurements and sweep-boundary snapshots are included.

Initialization, common warmup, compilation, deep copies and explicit GC are outside the primary timer. Allocation-triggered GC inside the measured loop is included. Raw allocation, GC and component timings are preserved in [timing_repeats.csv](timing_repeats.csv). The bootstrap interval measures repeat-to-repeat timing jitter on this host, not performance variation across machines.

The direct baseline rebuilds subsystem Grover matrices and refactors determinants for proposals. The optimized path maintains their determinant and inverse, uses rank-one updates, avoids redundant full candidate-state copying/rebuilding, and reconstructs at registered stabilization boundaries with conditioning safeguards. This is a Grover-update measurement; it does not establish an isolated speedup for ordinary equal/unequal-time Green functions.

## Correctness

Every accepted timing case must reproduce sweep-boundary auxiliary fields, sweep acceptance rates, determinant signs and next RNG draws exactly. Full Green matrices, determinant logs and all ratio observations are checked numerically. Maximum discrepancies across every retained repeat are:

- `maximum_Green_absolute_error`: 0
- `maximum_logdet_absolute_error`: 0
- `maximum_measurement_absolute_error`: 0
- `maximum_measurement_relative_error`: 0
- `maximum_reported_deltaG`: 1.0425029e-12

These checks establish matched sampled trajectories at the recorded sweep boundaries. They are not an independent proof of the final entropy’s physical accuracy, which belongs to the external scientific benchmark campaign. Failed and unresolved sizes remain visible.

## Provenance

- Recorded UTC: 2026-10-22T17:26:40.230 to 2026-10-22T17:28:24.406
- Host / Slurm job: `il1348` / `6457169`
- CPU: `model name	: Intel(R) Xeon(R) Gold 6348 CPU @ 2.60GHz`
- Julia / threads / BLAS threads: `1.12.1` / `1` / `1`
- Seeds: `[940031, 940032]`
- Receipt SHA256: `00145bb4e312426df3fb956c06624ca6a27c3b81cc7eec92ded172b1b22328ec`
- Benchmark runner SHA256: `aefc720f6626f020b455b0eb80d398b96a910b627b17cb1968984ac3bd4c4e8d`
- Analysis script SHA256: `d0e7a6676450404fa02a8fdaa9bb10276acbffebdd3d62d6cd4a61ec7ea422f9`

Source hashes:

- `src/HubbardRenyiDQMC.jl`: `336f4393cd69b04ffc1b8c42f97207d719a4ab88cbbd1d343a8088234bffc721`
- `Project.toml`: `ea3b7eacfcf9a36e4652654a5e7b573ef27529edfbad2e272da2bafca1a7f58f`
- `src/dqmc_kernel.jl`: `b72bbbf158381b0f776c1fa259d532b069094c47e361d3490c19ea0f43ada4be`
- `Manifest.toml`: `5cefaf1acdb4fa23a456c2c58f20fb1fcebd54b3cd27bc448cb65b0fcfecb90e`
- `src/grover_rank1_update.jl`: `79064a0aa868ff7b5c665876bc5dd6eeeee17da2684b7a04807891c81d0d358e`

BLAS configuration: `LBTConfig([ILP64] libopenblas64_.so)`

Full metadata and every row are in [timing_summary.json](timing_summary.json); compact values are in [timing_summary.csv](timing_summary.csv).

Optional figure status: UNAVAILABLE. No module named 'matplotlib'
