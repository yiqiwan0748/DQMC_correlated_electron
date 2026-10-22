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

The [latest audited 100-row ledger](reporting/reports/full32_public_replay_v7/README.md)
records **32 PASS and 68 NOT_RUN**. The original serialization failure
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

## Thermal and ensemble-reference additions

[Three additional questions](thermal_bridges/README.md) passed on ISAAC job
6457908, with independent provenance/arithmetic audit and the same frozen
acceptance gates: P07 fixed-physical-mu energy response, X02 thermal-to-ground
entropy contamination, and X03 canonical versus matched-mean grand-canonical
finite-size differences. The largest cross-backend discrepancy is 9.48e-12.
The cooling bound certifies beta=64 and 128 only for the four-site fixture; it
is not a thermal-bias certificate for the historical ten-site CV comparison.
Raw values, all tested temperatures, exact source closure and rerun comparator
are in the [updated 100-row ledger](reporting/reports/thermal_bridges_v5/README.md).
The original 25-question tranche and all original failures remain immutable.

## Independent full DMRG additions

ISAAC job 6457944 completed the full frozen TeNPy/ITensor protocol for P11,
P12 and P15: both open eight-site fixtures and both coarse/refined convergence
arms. Symmetry-resolved Schmidt probabilities agree within 8.93e-10,
doped site densities within 8.74e-12, and every registered charge-parity string
within 5.86e-10. Independent same-sector sparse ED energies/residuals/gaps,
MPS norms, particle/spin sectors, energy variances and refinement stability
passed the registered gates. The complete run used about 65 seconds for the
TeNPy oracle and 191 seconds for the user-ITensor adapter, about 1.1 GiB peak RSS.
These are finite-fixture results, with optimizer credit retained for the
independent libraries. They do not certify all DMRG model/size regimes.

The [latest immutable audit](reporting/reports/full32_public_replay_v7/README.md)
replays frozen comparators and preserves all source hashes and earlier reports.
One cross-platform Fourier-error replay difference of 8.85e-18 is recorded as
machine roundoff; scientific acceptance gates remain unchanged. The earlier
strict receipt-equality audit is preserved rather than rewritten.

## Repeated DQMC interval coverage

The frozen X06 study passed for the registered four-site, half-filled Hubbard
fixture (U/t=2, beta*t=1, time step*t=0.1). The actual corrected DQMC kernel
produced 400 independent eight-chain trials; nominal 95% Student intervals
contained the independent finite-time-step energy in **378/400 trials (94.5%)**.
The 99% exact coverage interval is [0.908725, 0.970173], and the 99% mean-energy
bias interval is [-0.000566704, 0.000793669] in hopping units. All 400 trials
passed their frozen convergence and effective-sample gates. No trial was removed
or rerun. Sampling job 6457947 completed in 548 seconds, about 824 MiB peak RSS.

The [raw paired-bin data and frozen audit](results/distinct_X06_final01/) and
[separate standard-library implementation](dqmc_coverage/audits/independent_stdlib_audit.py)
reproduce all 400 trials, 3,200 chains and 204,800 paired bins. The
[latest ledger](reporting/reports/full32_public_replay_v7/README.md) was rebuilt
using analysis code checked out from the exact published commit, with all 53
input hashes verified. This validates this energy-interval procedure for the
specified fixture; it does not certify entropy intervals, other regimes or the
historical ten-site calculation's systematic errors.
