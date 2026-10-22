# Reproduce the 38-case evidence roll-up

This is postprocessing only. It runs the frozen comparator on retained results; it does not sample a model, run a solver, or submit a cluster job. Use Python3.11+ with its standard library from the repository root containing `distinct_benchmarks/`. Source and raw-receipt hashes are checked before any completion is added.

```bash
python3 distinct_benchmarks/reporting/merge_two_orbital_tranche.py \
  --independent-audit distinct_benchmarks/results/two_orbital_analytic_full_6458671/independent_audit_attempt01.json \
  --independent-replay distinct_benchmarks/results/two_orbital_analytic_full_6458671/comparison_independent_replay_attempt01.json \
  --independent-auditor distinct_benchmarks/results/two_orbital_analytic_full_6458671/independent_auditor_attempt01.py \
  --output distinct_benchmarks/reporting/reports/two_orbital_v9_reproduction
```

The output directory must be new. From the local `advisor_preparation/` workspace, the same relative paths apply. The immutable published report is `reporting/reports/two_orbital_v9/`; do not overwrite it. A fresh report has its own timestamp and path receipts, while the frozen comparator result must reproduce the original scientific values/statuses within the declared binary64 replay tolerance.

The base is the exact 35-case `spectral_v8` ledger SHA256 `25b6733504af033cdbb6d0a942f86bcc8eacbe33c48451602c1af7b5b1f22c23`. Its 128 files, all 35 passing rows, and the other 97 row objects remain intact. New evidence is the complete A15/A16/A17 run6458671, its 14-member source release, four verbatim production ALM sources verified against the original baseline archive, all raw grids, the independently authored audit and a new replay of the exact frozen comparator. The independent auditor's execution copy is retained alongside its receipt so later edits elsewhere cannot silently replace that audit implementation.

The result is exactly 100 distinct catalogue IDs: **38 PASS, 62 NOT_RUN**. The 236 A15,47 A16 and513 A17 gates are acceptance checks within three questions, not additional benchmark ideas. Every gate retains its own value, inequality, bound and eligibility flag. Positive-density and minimum-moment gates remain lower-bound checks. A25 remains NOT_RUN; no interacting DQMC, Hund coupling or DMRG result is inferred from this analytical tranche. The original failures and corrected-source histories remain in the copied base evidence.
