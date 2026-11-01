# 100 distinct external scientific benchmark ideas

Execution progress is recorded separately in the [distinct scientific benchmark campaign](../distinct_benchmarks/README.md) and its audited 100-row ledgers. The catalogue below preserves its original question definitions and planning statuses.

[CSV catalogue](benchmark_ideas_100.csv) · [Structured JSON](benchmark_ideas_100.json)

This is a longer research programme: 25 analytical questions, 25 independent-package questions, 25 external-reference error/limit questions, and 25 published journal targets. Each row defines a scientific question, model/geometry, state or ensemble, observable normalization, independent oracle, proposed acceptance rule, resource tier, required adapter, execution status and relation to the initial ISAAC batch. Internal metadata, file-format checks and repeated U/beta/size arms are not counted as separate ideas.

The initial ISAAC batch contains **100 quantitative comparisons grouped into a few oracle families**. It supplies early coverage, rather than completing these 100 different ideas. Its 80 QuSpin references are available; the user-package comparisons have separate results. The catalogue's `initial_batch_mapping` column links relevant rows without turning partial coverage into a completed scientific programme. For example, two atomic thermodynamic points cover part of A01; they do not establish the atomic response over its full parameter range. All new roadmap questions remain unexecuted as defined, with one separately identified historical partial journal-reproduction attempt.

The analytical rows use directly derived finite-state/mode sums, published exact solutions or theorems. They are external truth targets, rather than algebra identities copied from the implementation. Independent package rows use [QuSpin](https://arxiv.org/abs/1804.06782), [TeNPy](https://tenpy.readthedocs.io/en/latest/intro/overview.html), [ALF](https://arxiv.org/abs/2012.11914), [PySCF](https://pyscf.org/quickstart.html) and [OpenFermion](https://quantumai.google/openfermion/tutorials/intro_to_openfermion). Their documentation establishes the available building blocks; the proposed fixtures and acceptance rules are this programme's design, not claims that those packages have already produced each reference here.

The journal rows request physical comparisons to primary papers, including [Broecker–Trebst Hubbard Rényi results](https://arxiv.org/abs/1404.3027), [Lieb–Wu exact-chain results](https://arxiv.org/pdf/cond-mat/0207529), and distinct spectral, thermal, magnetic and two-orbital targets. A reproduced picture is insufficient: match the source Hamiltonian, geometry, density/particle sector, temperature, boundary conditions, operators and normalization before evaluating agreement. The one-orbital chain, the intraorbital-U altermagnet model, and the full-Hund pnictide model are separate targets. Unsupported observables/models are explicitly marked.

## Sequence and acceptance

Start with exact small-system targets and reuse certified oracles across observables, while reporting one idea separately from its parameter arms, scalar values and batch jobs. Next add genuinely independent implementations: TeNPy for the ITensor-based DMRG, ALF for finite-temperature QMC, and PySCF/OpenFermion for different fermion/operator constructions. Then pursue the journal and thermodynamic targets that fit the scientific question. Shared JDQMC components make SmoQy comparisons useful partial cross-checks, but they cannot establish independence of the shared kernels.

A convergence trend becomes a benchmark only when it is calibrated against an independent target. For example, X08 measures DMRG truncation bias against external exact values; a larger bond dimension alone is not truth. Similarly, X01 compares a timestep extrapolation to an external continuum thermal trace. Error-control questions concern bias, uncertainty coverage, false ground-state convergence and state/ensemble selection, rather than internal implementation details.

The numeric tolerances in the CSV are **proposed preregistration values**, not published author tolerances. Before a run, freeze the full parameter sheet, reference provenance, equivalence margin, simultaneous comparison count and statistical rule. For stochastic accuracy, require the whole confidence interval of estimate minus oracle to lie inside the scientific margin; propagate oracle uncertainty. A small z score by itself is insufficient. Inadequate precision, unavailable functionality or exhausted convergence budget is unresolved. Preserve failed and unresolved results. Never relax a threshold after seeing a discrepancy without recording a new protocol.

For curves extracted from figures, retain the original extraction coordinates and an independently justified graphics allowance. Match raw numerical quantities before plots or fitted exponents. A finite-size theorem check, an approximate-method comparison or a nonzero finite-cluster magnetic signal does not establish a bulk phase. Renyi combinations need not obey all von-Neumann entropy inequalities.

## Resources and scope

| Tier | Planning scope | Important limits |
|---|---|---|
| T0 | Four-state/dimer/free-mode formulas and tiny dense calculations | Usually inexpensive after startup; not a measured runtime promise. |
| T1 | Small exact clusters, orbital response or modest independent-package calculations | Dense spaces grow exponentially; package installation/compilation and matrix sizes need a pilot. |
| T2 | Medium DMRG/DQMC comparisons, repeated statistical runs, controlled extrapolations | Precision, sign, bond dimension and oracle cost determine the budget. |
| T3 | Large-size, thermodynamic, dynamical or full two-orbital paper campaigns | May require substantial allocations and new implementations; no finish date or scientific completion is implied. |

These are unmeasured resource classes, not runtime forecasts. Give a pilot-based elapsed-time estimate, queue uncertainty and check-back point before submitting a scientific batch. The present catalogue builder launches no calculation and schedules no automatic job checks.

## Related and overlapping targets

A shared observable does not automatically make two questions identical, but changing only a parameter or software label never makes a new idea. The catalogue retains explicit overlap notes. Examples of distinct questions using related data are atomic versus canonical thermal entropy, scalar S2 versus a charge-resolved entanglement spectrum, finite-ring persistent current versus thermodynamic stiffness, and a tiny exact solver comparison versus a published large-system physical crossover. Treat the underlying data as shared evidence rather than independent confirmations.

Two direct duplicates were removed during review: generic orbital-weight analytics did not become a second Raghu orbital-weight reproduction; local spin-charge wavefronts did not become separate ideas merely because one reference was TeNPy and another was a journal. The retained replacements test the square-lattice DOS singularity and global-quench entanglement growth. The published spatial-cut entropy profile remains distinct from the higher-Renyi-index question. Close relationships are recorded in `overlap_note` for reviewer judgment.

The builder preserves [the delegated journal candidates](journal_candidates.json) and [the other candidates](base_candidates.json), then emits exactly 100 rows with unique IDs and titles. Rebuild with `python3 build_roadmap.py`; that checks the count and identifiers, not scientific execution.
