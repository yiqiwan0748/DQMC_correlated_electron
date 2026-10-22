# Reproduce the 35-case spectral roll-up

This is read-only receipt analysis, with outputs written to a new directory. It preserves the independently audited 32-case report at its exact SHA and replays only the frozen P08/P09/P10 comparator. It does not rerun physics or submit cluster jobs. Python 3.11 or newer is sufficient; NumPy and SciPy are not required for this leaf report.

From the `distinct_benchmarks` directory in the complete published record, run:

```bash
python3 reporting/merge_spectral_tranche.py \
  --independent-audit results/spectral_dynamics_6458468/independent_audit_attempt01.json \
  --independent-replay results/spectral_dynamics_6458468/comparison_independent_replay_attempt01.json \
  --output /tmp/spectral_v8_replay_new
```

The output directory must be unused. The tool checks the immutable base ledger and all its evidence, the exact eight-file spectral archive and source freeze, the actual user ED/thermal source, runtime/adapter receipts, complete fixture grids, and the independent audit and its exact authored source. It then executes the frozen stdlib comparator on copied raw receipts. Only P08/P09/P10 may change; all 97 other row objects, including the previous 32 PASS results, must remain semantically identical. Expected counts are 35 PASS and 65 NOT_RUN among exactly 100 unique catalogue IDs.

The archive, raw results, previous report, independent auditor, complete comparison and new reporting source are retained in the generated evidence manifest. Any missing source, hash mismatch, incomplete grid, failed scientific gate or disagreement stops the merge and writes `ROLLUP_FAILURE.json`; no new successful ledger is written. The additional results concern the declared four-site Hubbard fixtures and new observable adapters around actual user ED/thermal APIs. They do not establish native DQMC dynamics or larger-system accuracy.
