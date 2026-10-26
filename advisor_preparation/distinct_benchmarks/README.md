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

The [latest audited 100-row ledger](reporting/reports/weekend_reconciliation_v17/README.md)
records **49 PASS, 2 UNRESOLVED and 49 NOT_RUN**. The original serialization failure
and the original P16 scientific disagreement remain preserved.

The [26 October runtime prerequisite record](runtime_validation_20261026/README.md)
retains the rejected empty-mapping native PASS, both later Python preflight
successes, the unresolved native check and the independently audited full-width
inode diagnostic. The replacement native check is staged and NOT_RUN at its
recorded resource checkpoint. These checks add zero scientific credit; the
100-row ledger and all prior scientific evidence remain unchanged.
The [later native preflight attempt](runtime_validation_20261026/native_preflight_v3_timeout/README.md)
reached its 180-second limit before producing a native report and remains
UNRESOLVED, with zero scientific credit.

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

The [earlier immutable 32-question audit](reporting/reports/full32_public_replay_v7/README.md)
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
[earlier 32-question ledger](reporting/reports/full32_public_replay_v7/README.md) was rebuilt
using analysis code checked out from the exact published commit, with all 53
input hashes verified. This validates this energy-interval procedure for the
specified fixture; it does not certify entropy intervals, other regimes or the
historical ten-site calculation's systematic errors.

## Independent real-time, spectral and thermal-response additions

P08, P09 and P10 passed on ISAAC job 6458468. The actual user ED/thermal-ED
Hamiltonian and state routines, with declared new measurement adapters, agree
with independent QuSpin calculations for a charge-density quench (25 times),
all eight spin/momentum spectral channels with complete addition/removal sectors,
and 108 imaginary-time spin-response points. Quench observables agree within
8.89e-15, spectral poles/residues within 1.16e-14/1.00e-15, and spin responses
within 2.84e-15. Complete grids, cross-spin fermionic anticommutators, moment
sum rules, conservation, KMS and endpoint checks passed an independent audit;
the frozen comparison replay is byte-identical to its ISAAC receipt.

[Raw data and independent audit](results/spectral_dynamics_6458468/) and the
[updated ledger](reporting/reports/spectral_v8/README.md) retain source hashes
and all previous results. This is a four-site fixture with additive measurement
adapters; it does not test native DQMC dynamics. The job used one CPU3500MB and
completed in3m10s; target process peak RSS was about562MiB.

## Two-orbital analytical limits

A15/A16/A17 passed on job 6458671 after the separately audited resource pilot.
The complete one-electron ED spectra and orbital projectors agree with an
independently constructed one-body Hamiltonian (energy error1.20e-14; projector
error4.06e-15). All held-out2048-mesh DOS curves match the elliptic-integral
reference within8.13e-5 after identical binning/broadening; the final mesh
change is below9.71e-5. All27 decoupled two-orbital thermal fixtures agree
with independent product Gibbs states within2.23e-15, including every connected
spin component and nonzero local moments. This uses unchanged user ED/atomic/
band functions plus disclosed measurement adapters, without DQMC sampling.

The [38-question ledger](reporting/reports/two_orbital_v9/README.md) retains
all35 prior PASS rows and every historical receipt. All796 detailed gates,
including lower bounds and positivity, retain their own inequality directions.
A25 remains NOT_RUN because this release does not provide a native flux-threaded
one-dimensional ring. [Raw data and audit](results/two_orbital_analytic_full_6458671/).


## Ground-manifold comparison and retained precision failure

ISAAC job 6458708 evaluated both original certification questions in 5 min 53 s.
Scheduler completion produced **X12 PASS and X09 UNRESOLVED**. The
[original mixed-status leaf](reporting/reports/ground_initial_v10/README.md)
keeps all 38 preceding passing rows, the other 98 row objects and every prior
receipt. Its independently replayed comparator preserves 15 failed detailed
gates for three native DMRG starts, including signed variance, external energy,
ground-projector fidelity and the resulting certificate. The eight-sweep
attempt and its thresholds are not overwritten. This new certification adapter
used eight sweeps, compared with the production entry point's 14-sweep default.
Full Hamiltonian and MPO matrices agree with independent QuSpin; the observed
shortfall is optimization precision under that adapter budget, without a
demonstrated Hamiltonian or package defect.

X12 compares the complete canonical total-N3 dimer spectrum and ten prescribed
pure/coherent/mixed/Gibbs states, including both subsystem cuts and all physical
spin components. Its limiting state is the full equal ground-manifold mixture,
not an arbitrary single eigenvector. Manufactured excited MPS states in X09
are explicitly diagnostic challenges; they do not claim a natural optimizer
failure. [Raw outputs and independent audit](results/ground_certification_6458708/)
and a [self-contained reproduction command](reporting/GROUND_INITIAL_V10_REPRODUCE.md)
retain source/runtime closure and the complete physical grids.

## User-stop cancellation and portable handoff

The user requested all campaign work stop by 17:00 EDT on 22 October 2026, regardless of outcomes. Automatic follow-up is paused. The three stages below were cancelled before starting at 16:22:12 EDT. L32 pilot 6458711 completed at 16:22:48 EDT after 21m56s, exit 0 and MaxRSS 1,735,276 KiB. Its pilot quality PASS has been [independently audited](reporting/reports/bethe_l32_audit_6458711_portable_v1/README.md) with 177 checks and adds no catalogue completion. All campaign jobs are stopped before the deadline. [The initial stop receipt](submissions/user_stop_deadline_20261022.json) preserves the original submission history; [final accounting](submissions/user_stop_final_state_20261022.json) confirms all campaign jobs are terminal. These cancellations add no scientific failures or passes.

| Stage | Receipt | Scientific status |
|---|---|---|
| X09/X12 precision revision2, same fixtures and gates | [6458825](submissions/ground_certification_precision_v2_6458825.json) | Cancelled unrun for user deadline; original X09 UNRESOLVED retained |
| Native ALF full-matrix observer eligibility revision5 | [6458826](submissions/alf_measurement_proof_v5_6458826.json) | Cancelled unrun; P21 NOT_RUN; proofs and sampling still required |
| Excluded X01 continuum timing pilot | [6458828](submissions/continuum_bias_timing_pilot_6458828.json) | Cancelled unrun; X01 NOT_RUN; zero scientific completions |

The [precision revision](ground_certification_precision_v2/README.md) increases
optimizer precision without changing its original acceptance gates.
ALF revision4's [preserved startup failure](results/alf_adapter_proof_v4_6458709/failure_audit.json)
followed successful compilation/linkage and a pinned-provider check. The
[revision5 seed-input repair](alf_preparation/measurement_preparation_v5/README.md)
changes no Fortran, propagation, observer, proof or tolerance from revision4.
The [X01 pilot](continuum_bias/pilot_v1/README.md) measures timing and numerical
eligibility only; it performs no production continuum fit or equivalence decision.

A [portable handoff](PORTABLE_HANDOFF.md) and [machine-readable state](handoff_state.json) identify exact ledgers, frozen versions, raw receipts and replay entry points for another computer. L32 independent receipt audit passes; its finite-size pilot does not change the ledger. No automatic resume is authorized; a new explicit human instruction and fresh scheduler/resource reconciliation are required.


## October 23 scientific reconciliation

The [new immutable leaf](reporting/reports/pc_reconciliation_v12/README.md) reconciles the apartment-PC continuation against complete original raw files and frozen sources. Only X09, X13, X16, P21 and X01 change; all 95 unrelated row objects and 39 earlier PASS rows remain exact. X09 now passes the unchanged finite ground-certificate gates at sixteen sweeps. X13 passes its complete charged-threshold/visibility grid, and X16 passes all eight finite binding intervals under its declared conditional binary64 budget.

P21's actual native ALF and maintained user buffers remain UNRESOLVED: two component mixing gates and 376 complex precision limits fail across 704 complex entries. X01 remains UNRESOLVED at its complete second reserved stage; quartic entropy halfwidth 0.0051318999469 misses 0.005. The final stage is incomplete and contributes no primary estimate. Pilots, repeated X12, analysis repairs, setup-only PySCF and cancelled attempts add no completions.

[Portable source/data packets and replay commands](pc_continuation_20261023/README.md) preserve physical conventions, complete matrices/bins, all scientific failure history and explicit public selections. No private ZIP, controller policies, CV or contacts are required. The earlier October 22 stop remains fulfilled historical evidence; the user explicitly resumed work on October 23. Automatic follow-up remains paused.

## Independent PySCF spectra and directed transitions

Full production job 6465923 passes P23 and P24 at the frozen gates: seven complete finite fixtures and eight complete low-energy manifolds per fixture. The reference uses actual PySCF 2.14.0 FCI Hamiltonian/spin contractions and transition density matrices; the native target uses the existing two-orbital Hamiltonian with explicitly new dense-spectral and observable adapters. The asymmetric direction witness distinguishes the directed operator from its adjoint. All 3,503 comparator gates and 887 independent checks pass.

The [44-question leaf](reporting/reports/pyscf_v13/README.md) changes only these two IDs from the immutable 42-question leaf. All 98 other rows, every prior PASS and both unresolved histories remain exact. The [compact public scientific packet](pyscf_external_20261023/README.md) retains complete arrays, frozen sources, licences, installed-reference hash metadata and a fresh replay from the actual compressed payload. Public replay passes 847 scientific/integrity checks; private remote-inode/transport checks are explicitly excluded.

Pilot 6465868's native and reference producers succeeded, then comparator JSON serialization failed. Its original failed execution and traceback remain intact. Analysis-only repair 6465910 converts NumPy Boolean values for JSON, preserves all scientific formulas/gates and retains original plus repair costs. Both add zero credit. Production v2 changes that single Boolean conversion only. The interaction is intraorbital Hubbard U without Hund/interorbital/pair-hopping terms; selected finite-cluster agreement does not establish thermodynamic or whole-package reliability.

## Orbital entanglement with unchanged density

X17 passes full production job 6466242: six fixed open-cluster fixtures, all five total-N4 spin sectors per fixture and complete 256-level orbital spectra for both cuts. The actual native Lanczos state remains primary; a basis/mask metadata carrier lets the unchanged repaired native entropy function evaluate orbital partitions without generating a state from its carrier Hamiltonian. Independent QuSpin states/partial traces and explicit fermion/RDM/spin audits agree. All 1116 frozen gates and 703 independent checks pass.

The [45-question leaf](reporting/reports/orbital_entanglement_v14/README.md) changes only X17. The [scientific packet](x17_external_20261023/README.md) retains one unchanged frozen source release, complete pilot/full raw data, original scientific execution/resource receipts, both independent auditors and original reports. Compressed-payload replay passes 672 full public checks with bitwise-identical case metrics. Pilot 6466078 remains zero-credit: 185 frozen gates, 173 original independent checks and 81 selected public mathematical/resource checks.

The registered t3=0 model has uniform spin-resolved densities while t4 produces orbital entropy and nonzero connected density/spin-dot correlations. New dense certificates and full-RDM/observable adapters are attributed; only intraorbital U and selected finite fixtures are covered. Conditional binary64 uniqueness/gap budgets are not directed-rounding proofs. No full-Hund, bulk ALM, native-DQMC or whole-package claim follows; all other 99 row objects and both unresolved histories remain exact.

## L48 DMRG convergence pilot remains unresolved

The [complete original L48 pilot and saved-data replay](bethe_l48_pilot_20261023/README.md) preserve six completed sector/precision solves and all ten original quality failures. Half-filled results are stable; charged-state variance, particle-hole and refinement checks fail. The original strict audit and the complete 215-check audit remain available, with the public transport omissions explicitly documented. The stored final truncation value is a leftward half-sweep maximum, as a separate source-semantics diagnostic explains. No gate or original result changed. This pilot grants zero catalogue credit; A18/A19/A20 remain NOT_RUN. Its publication retained the immutable 45-question checkpoint.


## Prescribed orbital-antialignment trial functional

J23 passes full ISAAC job 6466733: 112 native evaluator rows, 196 independent analytical/quadrature rows and all 3,542 frozen gates. The native returned fields, densities and pure occupations are qualified before comparing the physical kinetic-plus-interaction energy. Independent projectors, Wick contractions and an exact-Fraction moment enclosure reproduce the prescribed small-hybridization energetic mechanism in the journal appendix; the largest leading-coefficient error is 6.16e-9 against the fixed 1e-6 gate. This concerns prescribed density-two Slater trial states, with no self-consistency, gRISB, doped phase or exact-ground-state claim.

The [46-question leaf](reporting/reports/j23_mean_field_v15/README.md) changes only J23; all 99 other row objects, all 45 earlier PASS rows and both unresolved histories remain exact. The [complete scientific packet and portable replay](j23_external_20261023/README.md) retain every scientific array, original source and independent auditor. The original full audit includes 646 source/metadata and 2,842 independent physical checks. Public replay executes the unchanged scientific fragments and verifies the actual pinned compute comparator replay, rather than claiming to repeat private transport checks. The 599-check excluded pilot and two saved-data replay jobs grant zero additional credit.

## Preserved spectral-reference precision failure

The [X14 failure and diagnostic packet](x14_failure_20261023/README.md) retains the original completed native rendering, failed independent reference and trace-only diagnosis. No comparison was produced. A separate exact-rational, saved-context analysis shows that the original derivative padding cannot certify the required bracket width within the retained failed bracket under its declared conditional arithmetic model. This is a bounded diagnostic, with no root-absence or complete-census claim. Original padding, gates and failure remain unchanged; X14 stays NOT_RUN and this pilot grants zero credit. Explicit scientific selections and derived traceback projections exclude private transport and workstation paths.

The X14 replay manifest SHA256 is `1c928e25fa732e599ef41acc123ff8566612eecf7cb04d5f804429059bd9364b`; the metadata and conditional-analysis replay commands are in its README.

## Free-ring flux response

[A25](a25_external_20261023/README.md) adds one analytical comparison for the explicit free one-dimensional ring adapter using the actual native complex-sector Hamiltonian and dense ground-state APIs. The two fixed rings and two gauges produce 196 point records, 980 sector solves and 56 current, curvature and stiffness responses. All 6,072 frozen comparator gates and 13,436 compute audit checks pass; an independent actual-result review checks 9,462 bindings and physical qualifications. The largest propagated conditional error bound is 6.44e-9 against the fixed 1e-8 gate. This establishes one free-ring result; interacting current, bulk stiffness and whole-package reliability require separate comparisons.

The complete sources, arrays, resource receipts, original pinned audits and byte-identical pinned comparator replay are retained. Public saved-data replay calls the unchanged original mathematical helpers; derived receipt details differ by at most 4.44e-16 (two pilot and 50 full floating-point leaves), with every difference recorded under the explicit receipt-only policy. All original physical gates, literal bounds, decisions and report bytes remain unchanged. The excluded pilot and two saved-data analysis jobs add zero further credit. All 99 other catalogue rows, all 46 earlier PASS rows and both unresolved histories remain exact.
