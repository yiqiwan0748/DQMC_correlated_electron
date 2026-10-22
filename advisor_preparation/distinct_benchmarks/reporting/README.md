# Read-only 100-idea status ledger

`merge_status_ledger.py` uses Python 3.11+ standard libraries and runs no physics. It audits the frozen source archive, receipt hashes, exact target/oracle values and metric-specific tolerances, then emits exactly one row per catalogue ID. Parameters, matrix elements and convergence arms do not increase the distinct-idea count.

From the workspace root, audit the original attempt:

```sh
python3 advisor_preparation/distinct_benchmarks/reporting/merge_status_ledger.py \
  --result-dir advisor_preparation/distinct_benchmarks/results/tranche01_attempt01 \
  --output-dir /new/report/revision
```

Select the valid journal serialization rerun explicitly:

```sh
python3 advisor_preparation/distinct_benchmarks/reporting/merge_status_ledger.py \
  --result-dir advisor_preparation/distinct_benchmarks/results/tranche01_attempt01 \
  --journal advisor_preparation/distinct_benchmarks/results/journal_attempt02/journal.toml \
  --journal-driver advisor_preparation/distinct_benchmarks/journal/run_raghu_v2.jl \
  --output-dir /another/new/report/revision
```

The original malformed journal remains copied and hashed as a predecessor. The tool verifies the exact `Float64.(center)` serialization-only source transformation. Malformed or missing execution receipts produce `UNRESOLVED`; a logged PASS alone cannot establish completion.

For a later **complete seven-ID package-suite rerun**, add `--packages NEW_ATTEMPT/packages.json`, `--package-ed-source .../compute_s2_ed_v2.jl` and `--package-ed-patch .../compute_s2_ed_v2.patch`. The target must be `NEW_ATTEMPT/target.toml`; the independent oracle defaults to the original attempt's `oracle/` directory, or can be supplied with `--package-oracle-dir`. The patch is applied only in memory to the original audited source, checked against the selected source bytes and the independently reviewed patch hash. Original package results, including P16's scientific FAIL, remain predecessor evidence. An isolated candidate result cannot substitute for the complete suite. The source catalogue and raw receipts remain unchanged.

Every output revision requires a new directory. `inputs/` contains exact copied bytes, including the reporting tool itself; `input_manifest.json` records relative/logical paths, hashes and sizes. No CV files are collected. `ledger_100.json` includes model/ensemble scope, existing functions versus new adapters, metric-specific errors, provenance amendments and historical statuses. `ledger_100.csv` and `README.md` provide compact views. A11 uses its actual L1 spectrum error rather than an empty generic maximum-absolute field. No ED result implies a DQMC sampling or DMRG pass.

The first audit with valid journal results is `reports/tranche01_journal02_v3`: **24 PASS, 1 FAIL, 75 NOT_RUN**. P16 remains FAIL until the actual source repair passes an independently matched full-suite rerun. The internal initial `tranche01_attempt01_v1` report predates logical-path handling and is excluded from Git; publication should use later revisions.
