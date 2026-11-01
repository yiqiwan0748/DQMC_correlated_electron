# Execution history

[Slurm parent-job history](slurm_job_history.csv) retains completed, failed and
cancelled submissions. Start/submit/end fields use the cluster's local clock;
scientific TOML receipts use their explicitly recorded UTC timestamps. Parent
TotalCPU fields may be zero when accounting is stored on job steps; they are not
used to calculate the reported speedups.

The initial oracle exporter failed on four mixed-state entropy outputs. Its
versioned repair produced all 80 references; the first attempt remains preserved.
Original-kernel atomic sampling then exposed a genuine reducibility defect.
Corrected-kernel results retain all earlier failed or unresolved evidence.

Job 6456921 passed 382 package checks and stopped in a diagnostic world-age loading
error. Its exact failed helper is `../diagnostics/diagnose_atomic_odd_sweep_failed_v1.jl`.
The versioned v2 helper completed job 6457036's 16 actual-kernel diagnostics, followed
by the external atomic cases. This repair changed helper loading, not physics.

Job 6457271 was cancelled after 21 seconds because Slurm's comma-separated export
syntax truncated the intended case selector. It produced no accepted comparison
receipt. The replacement wrapper records the selector inside the script;
job 6457272 completed the three intended precision cases. No seeds or acceptance
margins were changed.
