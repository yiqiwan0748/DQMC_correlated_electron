# Hubbard numerical research evidence

Numerical results, source code and external checks for research conversations.
Updated 23 October 2026. [Reproduction on ISAAC](isaac/ISAAC_REPRODUCE.md).

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

The new [distinct-benchmark campaign](distinct_benchmarks/README.md) records 45 passing questions and two unresolved
precision/mixing questions on ISAAC. Its latest
[independently audited ledger](distinct_benchmarks/reporting/reports/orbital_entanglement_v14/README.md)
records **45 PASS, 2 UNRESOLVED and 53 NOT_RUN**. These checks found a missing fermionic
regrouping sign in ED entropy for noncontiguous regions. The minimal production
repair passed all seven original QuSpin questions and all 16 four-site subsets;
only the ED entropy source file changed. The working package includes the
verified repair. Original failures, exact permutation proofs, source versions
and raw receipts remain visible. TeNPy/ITensor DMRG and repeated-campaign DQMC coverage pilots completed
with successful diagnostic/resource audits; neither pilot adds a distinct PASS.
The full independent TeNPy/ITensor DMRG comparison passed all three registered
questions: symmetry-resolved entanglement, doped density profiles and nonlocal
charge-parity strings. Both convergence arms and both eight-site fixtures passed
their energy, variance, sector and observable gates. The fixed 400-trial
DQMC interval-coverage study passed: nominal 95% energy intervals covered the
independent finite-time-step target in 378/400 trials (94.5%). All trial-quality
gates passed, and a separate standard-library audit reproduced the result.
This is one defined energy-error question at the registered four-site fixture;
it does not establish coverage for entropy or other model/size regimes.
Three further thermal-ensemble/response questions passed independent QuSpin
checks: cooling-to-ground entropy, fixed-total-number versus matched-mean grand
ensembles, and fixed-chemical-potential energy response. These four-site fixtures
do not certify the thermal bias in the historical ten-site calculation. Each
result identifies the actual method and functions exercised.

Three further QuSpin comparisons passed for real-time charge-density evolution,
addition/removal spectra and imaginary-time spin response. Complete physical
grids and independent moment/conservation checks passed, with explicit credit
to the existing user ED APIs and newly added measurement adapters.

Three two-orbital analytical limits also passed: complete one-electron orbital
projectors, the square-lattice DOS singularity with a held-out momentum mesh,
and thermal factorization of independent interacting orbitals. Source versions,
normalizations, smoothing and all inequality gates remain explicit.


## Ground-manifold comparison and original precision history

[X12 passed](distinct_benchmarks/reporting/reports/ground_initial_v10/README.md)
for the interacting canonical three-particle dimer: pure and coherent spin states,
equal ground-manifold mixtures, basis-rotated mixtures and finite-temperature
Gibbs states are compared against their own independent QuSpin references.
X09 was **UNRESOLVED** at its original frozen eight-sweep DMRG budget.
The new certification adapter used eight sweeps; the production entry point
has a 14-sweep default. Full Hamiltonian and MPO matrices agree with independent
QuSpin. This attempt shows insufficient optimization precision under the adapter
budget, without demonstrating a Hamiltonian or package defect.
Three of four native starts miss the variance requirement, with their energy,
fidelity and certificate checks retained. Small final truncation does not replace
external ground-state certification. The complete original states and 15 failed
detailed gates remain visible; no acceptance threshold changed.

A [separate precision revision](distinct_benchmarks/ground_certification_precision_v2/README.md)
retains the same questions, physical fixtures, starts and scientific gates.
Its job 6458825 was cancelled before starting for the user deadline; it adds no pass or failure. Native ALF observer
revision4 compiled and passed the pinned-library check, then stopped because
its fresh proof directory lacked the required seed input. The [revision5 input
repair](distinct_benchmarks/alf_preparation/measurement_preparation_v5/README.md)
was cancelled before start as job 6458826; P21 was NOT_RUN at that October 22 checkpoint.
The [X01 timing pilot](distinct_benchmarks/continuum_bias/pilot_v1/README.md),
cancelled before start as job 6458828, is excluded from future production fits and contributes zero
scientific completions. Submission receipts preserve the serial dependency chain.


The user requested a hard stop by 22 October 2026, 17:00 EDT (21:00 UTC), regardless of outcomes. Automatic follow-up is paused. The [stop receipt](distinct_benchmarks/submissions/user_stop_deadline_20261022.json) confirms all three queued stages were cancelled before starting. The L32 pilot completed at 16:22:48 EDT (21m56s, exit 0); its pilot quality PASS has been independently audited with 177 checks and adds no catalogue completion. [Final stop accounting](distinct_benchmarks/submissions/user_stop_final_state_20261022.json) confirms all campaign jobs are stopped before the deadline. A [portable public handoff](distinct_benchmarks/PORTABLE_HANDOFF.md) records source versions, raw receipts and replay commands. The user explicitly resumed work on October 23; the stop receipts remain unchanged historical evidence. Automatic follow-up remains paused.


## Reconciled apartment-PC continuation

The [October 23 scientific leaf](distinct_benchmarks/reporting/reports/pc_reconciliation_v12/README.md) adds three distinct passes while preserving all 39 earlier PASS rows and 95 unrelated objects. X09's frozen sixteen-sweep release passes the original ground-certificate gates; the eight-sweep failure remains intact. X13 passes complete finite charged thresholds and visibility; X16 passes all eight finite binding intervals under an explicitly conditional binary64 error budget. New charged/pair measurement adapters use the unchanged native Hamiltonian. No thermodynamic binding or whole-package certification follows.

P21 remains UNRESOLVED across 704 complex native ALF/user Green entries: two individual-chain ESS gates and 376 precision limits fail. Its frozen production paired-chain/complex-disk protocol explicitly supersedes the older unexecuted sampling candidate. X01 remains UNRESOLVED because stage 2 quartic entropy halfwidth 0.0051318999469 exceeds 0.005. Its incomplete final reserved stage, cancelled parts and failed analyses remain visible; no pooling or fourth look is introduced. The [portable scientific packets](distinct_benchmarks/pc_continuation_20261023/README.md) include full selected text data, source and mathematical replay, with controller/private handoff material excluded.

The pinned [PySCF environment setup](distinct_benchmarks/pyscf_setup_v1/README.md) completed as job 6465769. This is setup only and adds no scientific pass.

## Independent PySCF finite-cluster comparison

The [44-question checkpoint](distinct_benchmarks/reporting/reports/pyscf_v13/README.md) adds exactly P23 and P24: complete low-energy spectra/spin labels and directed orbital-transition strengths, compared against actual PySCF FCI operators on seven finite two-orbital fixtures. All 3,503 frozen comparator checks and 887 independent checks pass. The [public source/data replay](distinct_benchmarks/pyscf_external_20261023/README.md) independently reproduces the complete result from retained operators/eigensystems, without solver reruns.

The original pilot's comparator serialization failure and its analysis-only Boolean repair remain visible with zero catalogue credit. The native Hamiltonian, reference producer, physical fixtures and scientific gates are unchanged. Only intraorbital Hubbard interaction is tested; no whole-package, full-Hund or thermodynamic certification follows. All 98 other catalogue rows, including X01/P21 unresolved histories, remain exact.

## Orbital entanglement beyond density redistribution

The [45-question checkpoint](distinct_benchmarks/reporting/reports/orbital_entanglement_v14/README.md) adds exactly X17. Six frozen finite fixtures compare actual native Lanczos states and the repaired native entropy API with independent QuSpin orbital reduced-density spectra and physical connected correlations. Hybridization changes the quantum correlations while density means remain fixed. All 30 sector certificates, 1116 frozen gates and 703 independent checks pass; a [replay of the compressed public data](distinct_benchmarks/x17_external_20261023/README.md) reproduces identical scientific case metrics.

The excluded pilot adds zero credit. All 99 other catalogue rows, all 44 prior PASS rows and both unresolved histories remain exact. The model uses analytically chosen t3=0 and intraorbital U only; this finite-cluster result does not establish a bulk ALM phase or whole-package reliability. Existing package APIs and new RDM/certificate/observable adapters are distinguished explicitly.
