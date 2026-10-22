# 100 distinct external scientific benchmarks

This campaign executes the questions in the
[100-idea catalogue](../external_benchmark_roadmap/README.md), with actual user
solver/source outputs compared against analytical solutions, independent
packages and primary journal references. It extends the earlier 100 passing
parameter comparisons and 188 solver arms. Those earlier counts do not complete
this programme.

One idea receives one status, irrespective of its number of parameters, seeds,
matrix elements or convergence arms. Scientific disagreement is FAIL; missing
execution, an unresolved convention, a failed receipt or inadequate precision
cannot PASS. Source, protocols and acceptance gates are frozen before execution.
Failed attempts remain available beside versioned repairs.

The [latest audited 100-row ledger](reporting/reports/tranche01_repaired_v4/README.md)
records **25 PASS and 75 NOT_RUN**. The original serialization failure
and the original P16 scientific disagreement remain preserved.

## First tranche

The frozen first release contains 25 questions: A01–A14, A22, A23, P01–P06,
P16, J20 and J21. The [release receipt](source_releases/tranche01-freeze.json)
hashes its complete adapter source bundle. ISAAC job 6457635 used one CPU and
3500 MB, and completed in 200 seconds. Scheduler completion alone does not mean
scientific success: each family's process and comparison status is retained.

- [Analytical protocol](analytic/PROTOCOL.md): 16 questions about atoms,
  interacting dimers, free thermal/dynamical observables, complete entanglement
  levels and exact Hubbard theorems. These call actual user ED/thermal-ED
  Hamiltonian and analysis functions with explicitly added measurement adapters.
- [QuSpin protocol](packages/PROTOCOL.md): seven questions about interacting
  entropy indices, spin multiplets, natural occupations, spin/charge structure,
  pairing and separated-region mutual information. Independent target and
  oracle states are constructed separately.
- [Journal protocol](journal/PROTOCOL.md): two questions about the published
  two-orbital band geometry and orbital composition, calling the unchanged user
  production band/real-space source. The reference version is explicitly
  Raghu et al.'s author preprint **0804.1113v1**, with chemical potential 1.45.
  The final APS PDF was inaccessible; this is not a claim to have compared its
  final numerical version or digitized a published curve.

The [original results](results/tranche01_attempt01/) and
[complete log](logs/distinct_tranche01_6457635.log) retain all successes and
failures. P16 exposes a disagreement for the noncontiguous union of two regions;
the individual one-site purities agree. Independent permutation and QuSpin checks identified missing within-spin
regrouping signs. The actual minimal production repair passed all seven package
questions and all 16 spatial subsets on ISAAC job 6457716 (maximum subset
purity disagreement 5.61e-15). The sparse contraction, Hamiltonian and Lanczos
solver remain intact; the working package includes the verified repair. The
[prior ledger](reporting/reports/tranche01_journal02_v3/README.md) preserves
its original FAIL, and the latest ledger retains it as hashed historical evidence.

The first journal attempt completed J20's calculation but failed to serialize
a coordinate represented as Julia's irrational pi type. Its partial output is
preserved. The [versioned repair](journal/v2_repair.json) converts only the
recorded coordinates to Float64, and reruns both journal questions under the
same physical inputs and gates.

## Scope of a pass

The first analytical and QuSpin tranche exercises user ED/thermal-ED source and
newly declared observable adapters. It does not execute corresponding DQMC or
DMRG measurement routines. The journal tranche exercises production free-band
construction and does not certify interacting magnetism. Method-specific
TeNPy/ITensor DMRG and repeated-campaign DQMC checks are separate tranches;
their pilots do not establish final scientific passes.

The reporting tool recomputes comparison errors from raw values, checks frozen
source hashes and produces exactly 100 catalogue rows. Error metrics stay
separate: for example the full entanglement-level distribution uses L1 error,
whereas component comparisons use maximum absolute error. Unrun questions
remain visible rather than being counted from related parameter checks.
