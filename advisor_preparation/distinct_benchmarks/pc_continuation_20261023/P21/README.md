# P21: audited evidence, scientific result unresolved

ISAAC production job 6461010 completed the fixed four-site periodic Hubbard comparison (t=1, U=2, beta=1, centered mu=0, dtau=.1). Native ALF 2.6 continuous-Gaussian Mz sampling and the user's discrete-Hirsch DQMC use the same finite-step V_K density operator. Both actual maintained GT0/G0T channels retain 11 times, 2 spins and every 4×4 entry: 704 complex/1408 real components. No overall Matsubara minus, averaging, filtering or dense replacement of primary buffers is introduced. The frozen correspondence/source proof is retained.

The independent raw audit passes 2,634 checks but P21 remains **UNRESOLVED**. Two user components miss the minimum per-chain ESS 16 (observed 15.590770315); 376 complex intervals miss the .01 half-radius precision requirement and 170 rectangles are not contained in the .02 equivalence disk. No rectangle is wholly outside the disk. These are inadequate mixing/precision findings, not a demonstrated implementation disagreement. Eight fixed independent chain-pair estimates, 64 retained bins each, 4,096 burn cycles and 65,536 retained cycles per backend chain remain unchanged.

The actual frozen production protocol supersedes the earlier unimplemented candidate. In particular it uses fixed-pair Student t with 7 degrees of freedom intervals instead of candidate Welch intervals and explicit complex-disk/structural-zero gates. The production source/freeze/archive preceded sampling. Audit PASS is distinct from scientific PASS; there is zero catalogue credit. The v5 fixed-configuration proof and excluded short pilot also add zero credit. Pilot data never enter production estimators.

## Reproduce without sampling

Use Python 3.12 with NumPy 2.4.6/SciPy 1.17.1 in an isolated environment, then:

```bash
python replay.py --output-dir /tmp/new-p21-analysis
```

The helper checks archive/member SHA256 values, rejects absolute/parent/link/duplicate paths, extracts to a fresh directory, and runs unchanged independent mathematics plus the frozen arithmetic comparator. It never invokes ALF, Julia or Monte Carlo. Every scientific raw/source byte is preserved. The derived transfer manifest excludes three controller submission/queue records and records the original transfer hash; administrative-selection checks are not statistical gates. Source-root aliases and existing 64-epsilon arithmetic replay differences are disclosed in audit_summary.json. Native binaries/libraries remain on ISAAC and are represented by execution/source/provider hashes; none is attached. ALF GPL3-or-later and origin notices are inside the retained source evidence.

Raw production, proof and excluded-pilot archives have complete readable member manifests. independent_statistics.json.gz contains only the independently reconstructed numerical grid. Earlier public ALF setup failures, cancelled job 6458826 and the old unexecuted protocol remain immutable in their original locations.

The [excluded pilot audit history](pilot_auditor_history/history.json) retains both initial arithmetic-conditioning rejections and the separately versioned auditor repairs. The final pilot audit is contained in the pilot data archive and the frozen production planning source; no samples were added or pooled.
