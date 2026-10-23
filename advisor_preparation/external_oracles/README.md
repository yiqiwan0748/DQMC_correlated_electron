# External QuSpin oracle exporter

This exporter constructs independent finite-Hubbard-chain references with QuSpin and SciPy. It imports no user package. The 80 rows in `cases.json` join the 20 analytic comparisons maintained by the campaign runner: 100 observable/model/ensemble cases, grouped into oracle families.

From an environment where QuSpin, NumPy and SciPy are already installed:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python quspin_export.py \
  --manifest cases.json --output-dir /new/campaign/oracles
```

The output directory must be new. Outputs are `oracle_results.jsonl`, `oracle_values.toml` and `source_receipt.json`. The receipt records package versions, source/manifest hashes, threading and Slurm job ID. `ORACLE_READY` only means a reference was computed. The separate user-package comparison must pass its quantitative acceptance gate.

The exporter uses QuSpin's fermionic partial trace for both spin species on each retained spatial site. The finite-dtau density follows the user's V_K measurement cut; its subsystem purity is `Tr(rho_A*rho_A)`. A non-Hermitian asymmetric transfer density is never passed to a Hermitian eigensolver. Ground-state targets with unresolved degeneracy are withheld.

[Scientific sources and conventions](BENCHMARK_SOURCES.md) explain the independent package, free Gaussian and four-state atomic references. The ISAAC execution completed all 80 references in job `6456511` using QuSpin 1.0.1, NumPy 2.4.6 and SciPy 1.17.1. Its five analytic dimer guards agree within `1.12e-16`. The exporter took 0.65 seconds after startup; scheduler elapsed was about 15 seconds. [Execution receipt](../isaac/results/quspin_attempt02/source_receipt.json), [scalar references](../isaac/results/quspin_attempt02/oracle_values.toml), [full rows](../isaac/results/quspin_attempt02/oracle_results.jsonl) and [independent integrity audit](final_oracle_audit.json) retain the evidence.

QuSpin returned mixed-state reduced density matrices with a leading sample axis: `(1,4,4)` or `(1,16,16)`. The exporter now removes only that singleton axis before tracing the physical matrix. The original shapes are saved in the result metadata. Earlier attempt 01 left those four entropy rows unresolved; its 76 other reference values agree with attempt 02. These are 80 available independent references; passing user-package comparisons are reported separately.

The external atomic comparison subsequently exposed a reducible deterministic local-proposal schedule in the user sampler. [Diagnosis and repair review](ATOMIC_SAMPLER_DIAGNOSIS.md) preserve the original ISAAC failure, an independent exact transition-graph proof and the limits of the proposed per-local idle correction. [Retrospective split-Rhat](retrospective_split_rhat.json) rejects three affected original arms; its other diagnostic passes do not certify equilibration and are excluded from corrected-kernel validation summaries. `atomic_diagnosis_manifest.json` records the relevant source and receipt hashes.

The corrected campaign now passes all 24 interacting atomic/thermal DQMC comparisons, after four registered precision extensions. The [final independent raw-bin/count audit](final_campaign_independent_audit.json) verifies the corrected source, all 28 stochastic inspections and the combined **100-case/188-method-arm** result. These are bounded finite-system comparisons; original failures and unresolved stages remain visible.
