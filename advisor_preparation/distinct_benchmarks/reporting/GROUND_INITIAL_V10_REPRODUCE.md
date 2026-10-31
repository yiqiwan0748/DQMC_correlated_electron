# Reproduce the original ground-certification record

Run this read-only postprocessing from the directory containing `distinct_benchmarks/`, using Python3.11+ and its standard library. It replays the frozen comparator on existing result matrices and vectors; it does not run DMRG, ED, Monte Carlo or a cluster job.

```bash
python3 distinct_benchmarks/reporting/merge_ground_initial_attempt.py \
  --independent-audit distinct_benchmarks/results/ground_certification_6458708/independent_audit_attempt01.json \
  --independent-replay distinct_benchmarks/results/ground_certification_6458708/comparison_independent_replay_attempt01.json \
  --independent-auditor distinct_benchmarks/results/ground_certification_6458708/independent_auditor_attempt01.py \
  --source-root distinct_benchmarks/reporting/reports/ground_initial_v10/frozen_sources/ground_user \
  --output distinct_benchmarks/reporting/reports/ground_initial_v10_reproduction
```

The explicit source root uses the five published source/environment snapshots in the immutable leaf, so the public reproduction does not need a private development checkout or Julia package installation. Every byte must match the frozen source and actual-execution receipts. Choose a new output directory; do not overwrite `ground_initial_v10`.

The result is **39 PASS,1 UNRESOLVED,60 NOT_RUN**, with exactly100 distinct IDs. X12 passes its complete prescribed-state dimer comparison. X09 remains unresolved at the original frozen eight-sweep budget: three native starts miss the precision requirements. Its original matrices, complex state coefficients, eligibility reasons and15 failed detailed variance/fidelity/energy/certificate gates remain visible. The thresholds are unchanged; low final truncation alone is not accepted as ground-state certification. Manufactured excited MPS states are labeled diagnostic challenges, not natural optimizer failures.

All38 prior passing rows and the other98 row objects are preserved. The source archive has eight members; actual target provenance closes the ED, thermal, DMRG, Project and Manifest files. Original and independently replayed comparator output agree exactly. The receipt audit passes as an audit of faithful mixed scientific outcomes; it reports one scientific completion and does not turn X09 into a pass. A subsequent precision attempt requires its own versioned source, results and leaf record, retaining this initial result.

The original ledger SHA256 is `abee1b7af4263979880d24afc679d1553567e2141caf16c41888c7c46cf2e228`. A new leaf has its own timestamp/path receipts; its statuses, preserved rows, raw hashes and frozen comparator must reproduce the original evidence.
