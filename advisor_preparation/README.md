# Hubbard numerical research evidence

Numerical results, source code and external checks for research conversations.
Updated 22 October 2026. [Reproduction on ISAAC](isaac/ISAAC_REPRODUCE.md).

## Research results

The [historical research record](research_record/README.md) preserves five
subsystem Rényi benchmarks, finite-size results through 22 sites and related
ED/DMRG evidence. The [evidence ledger](research_record/EVIDENCE_LEDGER.md)
identifies calculation conventions, source hashes and the scope of each claim.

For the periodic half-filled 10-site Hubbard chain with U/t=2, beta*t=40 and a
five-site subsystem, saved DQMC gives S2=1.716442 +/- 0.003040; ground-state ED gives
1.716552. All five subsystem comparisons lie within 0.6 nominal MC standard errors.
Finite-temperature/grand-canonical and ground-state/fixed-sector targets differ;
the ledger retains that distinction and unresolved systematic errors.

## Fresh with-versus-without performance comparison

[Measured timing report](isaac/results/matched_timing_report_20261022/README.md):
the same two-replica Rényi sampling node, with direct determinant recomputation
or optimized cached rank-one Grover updates. On ISAAC, seven matched repeats at
each of five half-chain sizes gave 1.46–2.67x faster node-loop execution. The
128-site case took 23.4104 s direct and 8.75864 s optimized. All 35 paired trajectory
and measurement checks passed. Compilation, initialization and a complete lambda
grid are outside this timing scope. Every size and raw repeat remains visible.

At 128 sites, cumulative allocated bytes fell from 61.58 GB to 40.97 MB per timed loop,
about 1503x less allocation churn. This is not a peak-RAM reduction claim.

The separately labelled [subsystem-size follow-up](isaac/results/matched_timing_cut48_report_20261022/README.md)
measures a 48-site subsystem in a 64-site chain: **3.035x faster**, from 6.44880 s
to 2.12467 s. The [full-system thermal-purity comparison](isaac/results/matched_timing_cut64_report_20261022/README.md)
gives 2.862x. These post-primary workloads retain all seven paired repeats each;
they were chosen after the original grid and use a different CPU model from that
grid. Every individual ratio compares matched variants on the same host.

## External scientific validation

The [initial 100-case batch](validation_design/README.md) compares the user's ED,
DMRG and DQMC against independently constructed QuSpin references, analytical
free-fermion formulas and exact Hubbard atoms. These 100 parameter/observable
comparisons use repeated reference families. They are not 100 different methods.

| Stage | Execution evidence | Status |
|---|---|---|
| Independent QuSpin oracles | Job 6456511; 80 references with analytic guards | Complete |
| ED and thermal ED | Jobs 6456548/6456512; 88 arms; largest error 3.8e-13 | Passed |
| DMRG | Job 6456512; 64 arms; largest observable error 2.4e-7 | Passed |
| Free DQMC | Job 6456512; 12 arms; largest error 3.6e-15 | Passed |
| Corrected interacting DQMC | 24 stochastic comparisons; fixed 4096/16384-measurement stages | Passed |

External atomic tests found a reducible systematic sampling scan. The
[diagnosis and repair](external_oracles/ATOMIC_SAMPLER_DIAGNOSIS.md) retain the
original failure, exact 32-state proof, source patch and corrected execution
receipts. The working copy adds an independent 50% hold before every local
proposal. Original historical source/data remain immutable. The corrected atomic
and finite-hopping stages pass all 24 stochastic comparisons after their
registered precision extensions. Original stochastic passes are not reused in
the corrected summary. The [final summary](isaac/results/summary/summary.json)
records **100 passing cases and 188 passing solver arms**, with no unresolved
or failed arms in this initial batch.

## 100 distinct external benchmark ideas

The [100-idea catalogue](external_benchmark_roadmap/README.md) contains distinct
analytical, journal-result, independent-package and externally calibrated error
questions. It is a broader research programme, with unsupported and unexecuted
questions explicitly identified. The small initial 100-case batch provides
partial coverage; it does not complete this 100-idea programme.

The new [distinct-benchmark campaign](distinct_benchmarks/README.md) has executed
25 questions on ISAAC. Its latest
[independently audited ledger](distinct_benchmarks/reporting/reports/tranche01_journal02_v3/README.md)
records **24 PASS, 1 FAIL and 75 NOT_RUN**. The failure concerns ED spatial
entropy for a noncontiguous union of regions. An independent fermion-permutation
proof identifies missing within-spin regrouping signs; a minimally patched
production-source regression is queued before the failure can be marked resolved.
The original failure and all source versions remain visible. Separate TeNPy/ITensor
DMRG and repeated-campaign DQMC coverage pilots are queued; neither pilot counts
as a completed distinct question. These staged checks keep the method exercised
by each result explicit rather than inferring DQMC/DMRG coverage from ED alone.
