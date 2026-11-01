This frozen tranche implements catalogue **X06**, empirical uncertainty-coverage calibration against the independently computed QuSpin finite-dtau energy. Read [PROTOCOL.md](PROTOCOL.md) and [manifest.toml](manifest.toml) for the exact fixture and acceptance rules. The pilot is four separate repetitions; the primary experiment is 400 repetitions. No run is represented as completed until its actual ISAAC receipt exists.

After staging the whole directory under `$VALIDATION_ROOT/source/distinct_benchmarks/dqmc_coverage/`, the parent operator submits the pilot:

```bash
sbatch --export=ALL,VALIDATION_ROOT="$VALIDATION_ROOT" \
  "$VALIDATION_ROOT/source/distinct_benchmarks/dqmc_coverage/pilot.sbatch"
```

It writes a new `results/distinct_X06_pilot01/` and independently audits its raw bins. Inspect `summary.json` for `PILOT_ONLY`, raw-audit status, quality failures, and measured cost. After that measurement, fixed final index ranges can run with `final_shard.sbatch`, for example a range 1–50 (seven further disjoint ranges are required for all400):

```bash
sbatch --export=ALL,VALIDATION_ROOT="$VALIDATION_ROOT",COVERAGE_FIRST=1,COVERAGE_LAST=50,COVERAGE_PILOT_SUMMARY="$VALIDATION_ROOT/results/distinct_X06_pilot01/summary.json" \
  "$VALIDATION_ROOT/source/distinct_benchmarks/dqmc_coverage/final_shard.sbatch"
```

The final scientific summary is generated only after all400 fixed IDs exist:

```bash
"$VALIDATION_ROOT/quspin-venv/bin/python" \
  "$VALIDATION_ROOT/source/distinct_benchmarks/dqmc_coverage/summarize_coverage.py" \
  --mode final --inputs "$VALIDATION_ROOT"/results/distinct_X06_final01/shard_*.toml \
  --output "$VALIDATION_ROOT/results/distinct_X06_final01/summary.json"
```

Source, external-oracle, seed and protocol hashes must match. Output files are never overwritten. Pilot and final seeds never overlap. No inference about the earlier Bonferroni-tail coverage or full package reliability follows from this single-fixture nominal95% calibration.
