# Additive later-tranche audit

`merge_later_tranches.py` runs postprocessing only. It never submits a job, reruns physics, edits the catalogue, or changes the immutable v4 report. It starts from exactly 25 PASS and 75 NOT_RUN in `reports/tranche01_repaired_v4/ledger_100.json`, pinned to SHA256 `d32231a48f54a46d5233af6df9d7c0d275adec4f0d0134ee8cf4f2f21d696de0`.

The latest completed selection is `reports/thermal_tenpy_v6_attempt02`: **31 PASS, 0 FAIL, 0 UNRESOLVED, 69 NOT_RUN**. It adds the three thermal and three full TeNPy questions to the original 25. All 49 copied input identities and every analysis-receipt hash were verified; the original 25 rows, including A11, remain byte-for-byte equivalent as JSON objects. The first v6 audit's strict replay-equality rejection remains preserved: one platform-dependent Fourier-error difference of `8.851379265772331e-18` required the explicit roundoff comparison below. Both raw values remain available; no scientific acceptance margin changed. The immutable v5 thermal-only report and both pilot histories are also retained. A later X06 report must still start from the original 25-row-pass v4 base and select thermal, full TeNPy and final X06 explicitly, avoiding duplicate additions.

Every report uses a **new output directory**. It copies the exact tool, v4 ledger and all v4 evidence; verifies the three selected frozen source archives against fixed SHA256 anchors and their complete member maps; checks actual source, adapter and manifest hashes; then reruns the frozen comparators on copied raw data. Raw comparison receipts, when present, must agree with the independently recomputed comparison: IDs, statuses, text, integer values and structure match exactly; floating-point replay allows only `64*binary64_epsilon*max(1,abs(a),abs(b))` and records every difference. The frozen scientific tolerances remain unchanged. Completeness checks require the observables declared for each fixture: P11 charge spectra in the doped fixture and P15 parity intervals in the half-filled fixture. They reject empty required spectra, missing charge sectors, absent required parity intervals, missing convergence arms and pilot receipts. The ITensor runtime versions are checked against the copied, pinned Project/Manifest. Source/archive failures, malformed files and missing dependencies produce UNRESOLVED for the affected ideas while retaining the other 100 rows.

Supported additions are P07/X02/X03 (thermal ED and ensemble/response adapters versus QuSpin), P11/P12/P15 (full two-fixture, two-arm ITensorMPS versus TeNPy), and X06 (exactly 400 final independent DQMC replicates). Parameters, seeds and convergence arms never add distinct ideas. Pilot TeNPy and X06 runs cannot PASS. Old FAIL and malformed receipts remain in the immutable base evidence; supply any later pilot/predecessor receipts with `--history` to preserve their exact bytes too.

From the workspace root, after full results have been downloaded, use explicit actual result directories:

```sh
python3 advisor_preparation/distinct_benchmarks/reporting/merge_later_tranches.py \
  --thermal-dir advisor_preparation/distinct_benchmarks/results/thermal_bridges_6457908 \
  --tenpy-dir advisor_preparation/distinct_benchmarks/results/tenpy_full_6457944 \
  --coverage-final advisor_preparation/distinct_benchmarks/results/distinct_X06_final01/shard_1_400.toml \
  --history advisor_preparation/distinct_benchmarks/results/tenpy_pilot_6457724/comparison.json \
            advisor_preparation/distinct_benchmarks/results/dqmc_X06_pilot01/summary.json \
  --output-dir advisor_preparation/distinct_benchmarks/reporting/reports/later_full_attempt01
```

The full result names above follow the submitted job names; select the actual downloaded paths. Thermal and TeNPy directories must contain `oracle.json` and `target.toml`; `comparison.json` is also copied and checked when present. For split X06 jobs, list **every** immutable final shard after `--coverage-final`; their IDs must cover exactly 1 through 400 with no duplicates or replacement. To audit only pilot evidence, replace `--coverage-final ...` with `--coverage-pilot .../pilot.toml`; this cannot increase the pass count. Omit a tranche flag to retain its existing NOT_RUN rows. An attempted but incomplete supplied tranche becomes UNRESOLVED.

Python 3.11+ standard libraries suffice for thermal/TeNPy analysis. X06 requires SciPy for the exact binomial, Clopper–Pearson and Student-t gates. Use `--coverage-python /absolute/path/to/quspin-venv/bin/python` if the main interpreter lacks SciPy; the interpreter is an executable path, not a shell command. Missing SciPy yields UNRESOLVED, never a trusted-summary fallback.

On ISAAC, run with the **existing QuSpin virtual environment** after copying only public scientific evidence. For example, with `VALIDATION_ROOT` pointing to the staged `advisor_preparation` directory:

```sh
cd "$VALIDATION_ROOT"
"$VALIDATION_ROOT/quspin-venv/bin/python" distinct_benchmarks/reporting/merge_later_tranches.py \
  --thermal-dir distinct_benchmarks/results/thermal_bridges_6457908 \
  --tenpy-dir distinct_benchmarks/results/tenpy_full_6457944 \
  --coverage-final distinct_benchmarks/results/distinct_X06_final01/shard_1_400.toml \
  --output-dir distinct_benchmarks/reporting/reports/later_full_attempt01
```

Required staging is limited to this reporting tool, the complete v4 report and its `inputs/`, the selected raw result directories, the three matching `source_releases/*-freeze.json` and `*-source.tar.gz` pairs, `isaac/results/quspin_attempt02/oracle_values.toml`, and the corrected package source closure: `analysis/compute_s2_ed.jl`, `analysis/audit_u4_l4_thermal_ed.jl`, `analysis/compute_s2_itensor_dmrg.jl`, `src/HubbardRenyiDQMC.jl`, `src/dqmc_kernel.jl`, `src/grover_rank1_update.jl`, `Project.toml`, and `Manifest.toml` under `validated_package/HubbardRenyiDQMC/`. No CV, private working files, accounts, home directory or complete personal workspace is needed. `--renyi-source-root` and `--coverage-oracle` can select differently located **identical hashed** scientific inputs.

Output contains `ledger_100.json`, `ledger_100.csv`, `README.md`, `input_manifest.json`, copied evidence in `inputs/`, verified source code in `frozen_sources/`, and hashed recomputed receipts/process logs in `analysis/`. Inputs use public relative/logical source paths. Existing A11 L1 error and all original per-metric units remain intact; no ED pass certifies DQMC or DMRG. All input identities and the exact reporting source remain available for independent replay.

Local postprocessing-only checks: a base-only audit preserved 25 PASS / 75 NOT_RUN; supplying the two actual pilot receipts preserved 25 PASS and classified the four attempted ideas UNRESOLVED. Local Python lacks SciPy, which is correctly an explicit X06 analysis dependency rather than a scientific pass. No new physics was executed to test this reporting tool.
