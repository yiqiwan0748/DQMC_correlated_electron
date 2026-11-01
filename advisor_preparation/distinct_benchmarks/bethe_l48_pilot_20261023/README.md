# L48 DMRG pilot: convergence remains unresolved

ISAAC job **6466266** completed all six planned calculations for the open
48-site Hubbard chain, U/t=4: half filling and one electron removed/added,
each at the fixed coarse/refined settings. The original quality summary is
**PILOT_QUALITY_UNRESOLVED**. This pilot earns **zero** catalogue completions;
A18, A19 and A20 remain NOT_RUN. No larger chain is admitted by this record.

The half-filled calculations are stable. The charged calculations fail the
original energy-variance, stored-observer truncation, particle-hole and
refinement checks. In particular, the charge gap changes from 1.3620264593 to
1.3525508511 between the two fixed arms, a change of 0.0094756081 against the
unchanged 0.00001 stability margin. Neither a thermodynamic gap nor agreement
with the Lieb–Wu bulk solution is established. Six finite solves do not replace
the full predeclared size series, extrapolation or held-out checks.

The original strict audit stopped at a genuine charged-sector variance failure.
Its source is retained; `strict_failure_projection.json` retains all original
receipt fields except the private workstation traceback, names that omission
and binds the unchanged private original by hash. The subsequent complete audit checks all
receipts and retains **10 quality failures across 215 checks**, with integrity
PASS and scientific quality UNRESOLVED. Original source, values, thresholds,
logs, source freeze and failed history are unchanged.

The target calls the user's existing Hubbard MPO; ITensorMPS supplies DMRG
optimization. Six sector/precision and local-observable adapters are explicitly
added benchmark code. The job took 1,908 seconds; native process peak RSS was
1,684,840 KiB. This is a resource/quality pilot and no speedup claim follows.

An additional diagnostic identifies a reporting distinction in ITensorMPS
0.4.1: `truncerrors(DMRGObserver())` records the maximum over the leftward half
sweep, while the console reports the maximum over both passes. The original
`final_truncerr` is the stored leftward value. `observer_semantics.json` retains
all six comparisons, exact installed-source hashes and line references. The
original field and gate are preserved. The independent variance, particle-hole
and refinement failures remain. Installed dependency files stay in the private
provenance record; their full source is not bundled here.

Reproduce the saved-data public audit using Python 3.11 or newer:

```sh
python3 -B replay_public.py --output /tmp/l48-public-replay
```

This command checks the published bytes, safely unpacks text inputs, applies
every original mathematical/quality statement, and reruns only the frozen
standard-library summarizer. It starts no DMRG, eigensolver, quadrature oracle or
sampler. It requires exact scientific fields and byte-identical summary output
relative to the original complete audit; the result must remain UNRESOLVED.

`public_input_inventory.json` explicitly selects scientific source/raw/logs and
own-job accounting, with hashes. Private transport/inode receipts are omitted;
their hashes and the original complete audit remain available for attribution.
The public auditor therefore reports its own check count and does not claim
to replay those private checks. Postexecution runtime hashes are disclosed as
postexecution observations, not preexecution binary pins. The original source
archive has 12 regular text members; the separate scientific-input archive
contains only regular UTF-8 files. No credentials, CV/contact content, whole-user
queues, private controller bodies or runtime binaries are included.
