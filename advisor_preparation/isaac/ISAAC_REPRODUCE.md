# Reproduce the ISAAC validation and timing

The later [distinct-question campaign](../distinct_benchmarks/README.md) has a
separate frozen source/receipt history. Its verified noncontiguous ED-entropy
repair changes only `analysis/compute_s2_ed.jl`; the three archives below retain
their original bytes. To reconstruct the new working source, first restore the
original corrected archive below, copy `source/renyi_corrected` to a new
`source/renyi_disjoint_corrected`, and replace that new tree's
`analysis/compute_s2_ed.jl` with
`distinct_benchmarks/packages/p16_diagnosis/compute_s2_ed_v2.jl`. Verify SHA256
`f8f94cfc13e11bb8b5e64e8ac8f4ade5ed3a13b3976c83b5ab5db89d86ba5f28` and the
single-file source change recorded in the raw source-release receipt. Preserve
both trees: the original and repaired stages use their own frozen hashes. See
the [reporting instructions](../distinct_benchmarks/reporting/README.md) for
explicit selection of repaired-source results and preservation of failures.

This guide restores the exact original package, its original external-comparison consumer, and the separately corrected package. It also stages the independent QuSpin exporter and the later timing script. It does not claim a run passed: use each result's status, source hashes and retained diagnostics. Historical advisor-facing values have a separate [reproduction note](../research_record/REPRODUCE.md).

## Published layout and frozen inputs

New execution receipts and summaries live under **`advisor_preparation/isaac/results/`**, matching the local working tree. The older `validation_results/` ED/QuSpin links are compatibility copies. Scientific reports should cite the canonical `isaac/results/` path. Slurm wrappers, runtime setup, protocols and diagnostics live under `isaac/`; downloadable source archives live under `validation_source/`.

| Frozen input | SHA-256 | Contents |
|---|---|---|
| [baseline-source.tar.gz](../validation_source/baseline-source.tar.gz) | `11c3c33adaa70a78b8ab10437fba04651e3951c040bb289c4a3ce7ec293420aa` | Original `renyi/`, `alm/`, early `isaac/`, `SHA256SUMS`, source provenance |
| [baseline-consumer-20261022.tar.gz](../validation_source/baseline-consumer-20261022.tar.gz) | `d25375c07acd2aea687d38da84dd506490135e2bfc4d63388a13f7fd1928888f` | Original `validation_design/`, retrieved unchanged from ISAAC |
| [corrected-source.tar.gz](../validation_source/corrected-source.tar.gz) | `ae174444854849faf5ddb604a2b4a3e8dc88ccdb70f4fb0d9b9ff065dc9636e6` | `renyi_corrected/`, `validation_design_corrected/`, selected `isaac/` support |

[Original archive receipt](../validation_source/baseline-package-receipt.json) and [corrected archive receipt](../validation_source/corrected_package_receipt.json) retain byte counts and hashes. The public baseline filename differs from its original `validation-source.tar.gz` filename; its bytes are identical. The original consumer uses `alpha=0.01/100`; the corrected consumer adds split-Rhat and uses `alpha=0.01/400`. Do not substitute the current consumer when claiming exact reproduction of the original statistical protocol.

The original source closure includes working-tree files identified by its hashes. It is distinct from the exact historical Git tree under `research_record/source_snapshot/`. The latter preserves the August research record; it is not silently replaced by the October corrected package.

## Restore a new scratch campaign

Run from a checkout containing `advisor_preparation/` on ISAAC. Use a new scratch directory; these commands never replace an accepted campaign. The published batch wrappers retain the original account `acf-utk0011` and ISAAC partitions/QoS. If your allocation differs, adjust a new submission copy or Slurm command-line overrides.

```bash
set -euo pipefail
export RECORD="$PWD/advisor_preparation"
export VALIDATION_ROOT="${SCRATCHDIR:?}/advisor-reproduce-$(date -u +%Y%m%dT%H%M%SZ)"
test ! -e "$VALIDATION_ROOT"
mkdir -p "$VALIDATION_ROOT"/{source,runtime,logs,results,jobs}
sha256sum "$RECORD/validation_source/"*.tar.gz
# Compare all three hashes with the table above before extraction.
tar -xzf "$RECORD/validation_source/baseline-source.tar.gz" -C "$VALIDATION_ROOT/source"
(cd "$VALIDATION_ROOT/source" && sha256sum -c SHA256SUMS)
tar -xzf "$RECORD/validation_source/baseline-consumer-20261022.tar.gz" -C "$VALIDATION_ROOT/source"
tar -xzf "$RECORD/validation_source/corrected-source.tar.gz" -C "$VALIDATION_ROOT/source"
```

The baseline checksum check precedes later support overlays. Neither source archive contains the full later exporter/timing closure. Stage those published files explicitly:

```bash
mkdir -p "$VALIDATION_ROOT/source/external_oracles" "$VALIDATION_ROOT/source/isaac/diagnostics"
cp "$RECORD/external_oracles/quspin_export.py" "$RECORD/external_oracles/cases.json" \
   "$VALIDATION_ROOT/source/external_oracles/"
cp "$RECORD/isaac/"*.sbatch "$VALIDATION_ROOT/source/isaac/"
cp "$RECORD/isaac/diagnostics/"*.jl "$VALIDATION_ROOT/source/isaac/diagnostics/"
cp "$RECORD/isaac/logs/quspin_requirements.txt" "$VALIDATION_ROOT/source/oracle_requirements.txt"
```

The expected layout is now `source/{renyi,alm,validation_design,renyi_corrected,validation_design_corrected,external_oracles,isaac}`. Preserve both package trees. All actual-run source hashes are recorded in results, including the later timing script's own hash.

## Prepare the runtimes

The historical [bootstrap script](bootstrap_runtime.sh) downloads Julia 1.12.1 and checks the official checksum, but hardcodes its original campaign directory. For the new directory above, perform the same versioned steps with the new root:

```bash
cd "$VALIDATION_ROOT/runtime"
curl --fail --location --retry 2 --output julia-1.12.1-linux-x86_64.tar.gz \
  https://julialang-s3.julialang.org/bin/linux/x64/1.12/julia-1.12.1-linux-x86_64.tar.gz
curl --fail --location --retry 2 --output julia-1.12.1.sha256 \
  https://julialang-s3.julialang.org/bin/checksums/julia-1.12.1.sha256
awk '$2=="julia-1.12.1-linux-x86_64.tar.gz" {print}' julia-1.12.1.sha256 > selected.sha256
test -s selected.sha256
sha256sum -c selected.sha256
tar -xzf julia-1.12.1-linux-x86_64.tar.gz
cd "$VALIDATION_ROOT/source"
sbatch --wait --export=ALL --output="$VALIDATION_ROOT/logs/setup-%j.out" isaac/setup.sbatch
```

Package installation/precompilation is separate from numerical validation. [setup.sbatch](setup.sbatch) instantiates the pinned original Julia environments in this campaign's isolated depot; the corrected package retains those dependency pins. Run dependency-backed work on compute nodes.

For QuSpin, use the published Python distribution list rather than allowing all transitive versions to float. The recorded environment was Python 3.11.11, QuSpin 1.0.1, NumPy 2.4.6 and SciPy 1.17.1. Make a new setup-job copy; the frozen original remains unchanged:

```bash
cp isaac/setup_quspin.sbatch "$VALIDATION_ROOT/jobs/setup_quspin_locked.sbatch"
python3 - <<'PY'
import os
from pathlib import Path
p=Path(os.environ['VALIDATION_ROOT'])/'jobs/setup_quspin_locked.sbatch'
s=p.read_text()
old='quspin==1.0.1'
assert s.count(old)==1
p.write_text(s.replace(old,'-r "$VALIDATION_ROOT/source/oracle_requirements.txt"'))
PY
sbatch --wait --export=ALL --output="$VALIDATION_ROOT/logs/quspin-setup-%j.out" \
  "$VALIDATION_ROOT/jobs/setup_quspin_locked.sbatch"
```

The historical module is `anaconda3/2024.06`. A different Python/platform may require compatible wheels; retain the actual distribution receipt instead of claiming binary identity across platforms.

## Generate independent references

The frozen [oracle wrapper](oracle.sbatch) writes `quspin_attempt01`. The consumer expects `quspin_attempt02`, so make a new wrapper with that one output-directory change inside the fresh campaign:

```bash
sed 's@results/quspin_attempt01@results/quspin_attempt02@g' \
  isaac/oracle.sbatch > "$VALIDATION_ROOT/jobs/oracle_reproduce.sbatch"
sbatch --wait --export=ALL --output="$VALIDATION_ROOT/logs/oracle-%j.out" \
  "$VALIDATION_ROOT/jobs/oracle_reproduce.sbatch"
```

The exporter refuses an existing output directory. Its 80 references must be `ORACLE_READY`, with analytic dimer guards passing. Reference availability does not certify any user-package solver. The published [source receipt](results/quspin_attempt02/source_receipt.json) and [independent audit](../external_oracles/final_oracle_audit.json) document the completed historical reference run; newly generated receipts record their own versions, hashes and job ID.

## Reproduce the original baseline and failure

[baseline.sbatch](baseline.sbatch) runs the original package regressions and the original direct/optimized equivalence checks. Its log files are supporting tests, excluded from the 100 external comparison count. The external deterministic stage uses the original consumer restored above:

```bash
sbatch --wait --export=ALL --output="$VALIDATION_ROOT/logs/baseline-%j.out" isaac/baseline.sbatch
VALIDATION_RUN_DQMC=0 VALIDATION_ATTEMPT=external_deterministic_reproduction \
  sbatch --wait --export=ALL --output="$VALIDATION_ROOT/logs/deterministic-%j.out" isaac/external.sbatch
```

To reproduce the diagnosed atomic failure, use the original kernel and original statistical consumer in a separate receipt. A numerical/statistical result is not required to reproduce every saved decimal:

```bash
VALIDATION_RUN_DQMC=1 VALIDATION_ID_PREFIX=ATOM- \
VALIDATION_MC_MEASUREMENTS=4096 VALIDATION_ATTEMPT=external_atomic_original_reproduction \
  sbatch --wait --export=ALL --output="$VALIDATION_ROOT/logs/atomic-original-%j.out" isaac/external.sbatch
```

The original [atomic receipt](results/external_atomic_mc_4096.toml), [exact diagnosis](../external_oracles/ATOMIC_SAMPLER_DIAGNOSIS.md) and [retrospective mixing audit](../external_oracles/retrospective_split_rhat.json) remain failure history. Do not include original-kernel stochastic PASS rows in a corrected summary.

## Validate the corrected kernel

Read the [registered corrected protocol](corrected_monte_carlo_protocol.md). The original [atomic wrapper](atomic_corrected.sbatch) passed **382 supporting checks** in job `6456921`, then stopped on a Julia world-age error in the diagnostic harness before running the eight external atomic comparisons. The [failed log](logs/atomic_corrected_6456921.log) remains preserved. This harness failure did not change the physical sampling kernel or the frozen corrected archive.

For a new reproduction, run the corrected package tests through a new tests-only submission copy, then use [atomic_corrected_resume.sbatch](atomic_corrected_resume.sbatch), which invokes the versioned [v2 diagnostic](diagnostics/diagnose_atomic_odd_sweep_v2.jl) and all eight external atomic comparisons. The second scientific job runs the 16 thermal comparisons:

```bash
python3 - <<'PY'
import os
from pathlib import Path
r=Path(os.environ['VALIDATION_ROOT'])
s=(r/'source/isaac/atomic_corrected.sbatch').read_text()
marker='export ATOMIC_DIAGNOSTIC_OUTPUT='
assert s.count(marker)==1
(r/'jobs/atomic_tests_only.sbatch').write_text(s.split(marker,1)[0])
PY
sbatch --wait --export=ALL --output="$VALIDATION_ROOT/logs/atomic-tests-%j.out" \
  "$VALIDATION_ROOT/jobs/atomic_tests_only.sbatch"
sbatch --wait --export=ALL --output="$VALIDATION_ROOT/logs/atomic-corrected-%j.out" \
  isaac/atomic_corrected_resume.sbatch
VALIDATION_ID_PREFIX=QT- VALIDATION_MC_MEASUREMENTS=4096 \
VALIDATION_ATTEMPT=external_thermal_corrected_4096 \
  sbatch --wait --export=ALL --output="$VALIDATION_ROOT/logs/thermal-corrected-%j.out" isaac/external_corrected.sbatch
```

Both load `renyi_corrected` with `validation_design_corrected`, use eight registered seeds, and write new receipts. Original snapshots stay unchanged. The exact atomic transition proof establishes one communicating class for strictly positive atomic field support; fresh solver comparisons establish the observable accuracy actually tested. Neither proves global ergodicity for every determinant-zero/sign-changing problem.

Keep unresolved cases visible. Only the registered extensions to 16,384 and 65,536 measurements per chain are allowed in this campaign. Set `VALIDATION_ID_LIST` to the unresolved IDs, clear any earlier `VALIDATION_ID_PREFIX`, set a new `VALIDATION_ATTEMPT`, and keep all original seeds and margins. Failed mixing is investigated before extending. Confidence statements are nominal under the stated equilibrated independent-chain approximation.

Build a corrected summary from the deterministic receipt and the corrected stochastic receipts only, in chronological order. The original pilot and original stochastic receipts are excluded:

```bash
python3 "$RECORD/isaac/summarize_campaign.py" \
  --manifest validation_design_corrected/benchmark_manifest.toml \
  --output-dir "$VALIDATION_ROOT/results/summary_corrected_reproduction" \
  "$VALIDATION_ROOT/results/external_deterministic_reproduction.toml" \
  "$VALIDATION_ROOT/results/external_atomic_corrected_4096.toml" \
  "$VALIDATION_ROOT/results/external_thermal_corrected_4096.toml"
```

Append registered extension receipts if used. Read both case and method-arm counts. The executable 100 comparisons repeat observable/model/ensemble families; the [100 distinct scientific ideas](../external_benchmark_roadmap/README.md) are a separate, largely unexecuted research programme.

In the completed recorded campaign, the [atomic](results/external_atomic_corrected_16384.toml) and [thermal](results/external_thermal_corrected_16384.toml) 16,384-measurement extensions resolve the four initially unresolved cases. The [final summary](results/summary/summary.json) has 100 passing cases and 188 passing method arms. The [independent final audit](../external_oracles/final_campaign_independent_audit.json), reproducible with `python3 external_oracles/audit_final_campaign.py` from `advisor_preparation`, reconstructs the saved-bin statistics and count composition without starting simulations. Keep the first-stage and original-kernel receipts alongside these final results.

## Measure the optimization fairly

The later [matched timing script](diagnostics/benchmark_matched_lazy_grover.jl) is staged explicitly above because it postdates the corrected archive. [timing_matched.sbatch](timing_matched.sbatch) runs a separate pilot and the registered size grid:

```bash
sbatch --wait --export=ALL --output="$VALIDATION_ROOT/logs/timing-%j.out" isaac/timing_matched.sbatch
```

The [timing protocol](diagnostics/MATCHED_TIMING_PROTOCOL.md) compares actual direct fused production updates with rank-one Grover caching, with the same local-idle policy, warmed fields, RNG state, measurements and one-thread controls. Seven repetitions alternate order for `L=8,16,32,64,128`, half-chain subsystems, `U=2,β=2,Δτ=0.05,λ=0.5`. The primary statistic is the ratio of median node-loop times; all sizes and failures remain reported. A 600-second soft budget may leave later sizes unresolved, and the batch exit status alone does not replace inspecting the receipt.

The correctness gates compare boundary fields, RNG draws, acceptance, signs, Green matrices and measurements. Initialization, compilation, deep copying and shared warm-up are excluded from the primary timing. The result measures this matched two-replica node loop, including measurements and shared snapshots; it does not measure a full λ-grid entropy campaign or isolate ordinary Green-function Sherman–Morrison updates. Report the measured scope, size, hardware and date with any CV speedup number.

The completed primary execution is [job 6456993's raw receipt](results/matched_timing_full_20261022.toml), with [all five size results](results/matched_timing_report_20261022/README.md) and an [independent median/source/replay audit](../external_oracles/independent_primary_timing_audit.json). All 35 paired repetitions passed the recorded checks. The largest registered size, `L=128,nA=64`, gives `23.410449674/8.758641974 = 2.672840×` matched node-loop speedup. Later cut-size explorations belong to separate receipts and must retain their post-primary disclosure.

The separate [larger-cut wrapper](timing_large_cuts.sbatch) uses the explicitly exploratory [runner](diagnostics/benchmark_matched_large_cuts.jl), already included by the support-overlay commands above:

```bash
sbatch --wait --export=ALL --output="$VALIDATION_ROOT/logs/timing-cuts-%j.out" isaac/timing_large_cuts.sbatch
```

It retains both `L=64,nA=48` and `L=64,nA=64` outcomes. Their recorded ratios are respectively [3.035193×](results/matched_timing_cut48_report_20261022/README.md) and [2.861522×](results/matched_timing_cut64_report_20261022/README.md), with an [independent audit](../external_oracles/independent_large_cut_timing_audit.json). These cuts were selected after observing the primary campaign. They ran on the same Xeon 6348 host, while the primary used Xeon 6248R; every speedup uses its own matched same-host baseline. The full-system cut measures a thermal-purity node. Preserve this scope and both results when quoting the 48-site-subsystem timing.
