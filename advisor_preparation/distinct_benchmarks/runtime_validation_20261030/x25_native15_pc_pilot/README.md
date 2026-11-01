# Native fifteen-case timing pilot

This record contains one completed native calculation attempt and saved-data audits for a small two-orbital fermion model. All 15 requested cases produced complete records. The available original native criteria and independent saved-array consistency checks passed. This is a pilot record with **zero scientific benchmark credit**; X25, the independent QuSpin comparison, and the external ground-state and response certificates remain **NOT_RUN**.

The model has an open 2×2 arrangement of four cells, two orbitals per cell, eight spatial modes, four particles, and centered interaction U=6. Its hopping parameters are t1=−1, t2=−1.75, t3=0, and t4=−0.65. Each of h=0, −0.02, and +0.02 has all five sectors N_up=0…4, N_down=4−N_up. Their dimensions are 70, 448, 784, 448, and 70. The field convention is H(h)=H(0)−hO, with O=Σ_cells(Sz_x−Sz_y); M=O/4 is the value per cell. This field conserves the two spin particle counts.

The unchanged native API used Lanczos with atol=rtol=1e−13, maxiter equal to the sector dimension, and required convergence. The largest recorded residual was 1.92330347614013e−12, below the original 1e−11 native criterion. These residual and saved-array checks do not independently establish the lowest eigenstate, uniqueness, a gap, or external accuracy. Unbalanced zero-field sectors may be degenerate: their single-vector observables can depend on the chosen vector. No ground-subspace averaging was performed. Balanced-sector symmetry criteria were applied only to the balanced sector.

The whole external attempt took 41.063506921 seconds (rounded here; the JSON preserves its exact recorded scalar). The native Julia entry took 33.876808775 seconds, including startup, definitions, guards, metadata and output work. Producer elapsed time was 9.136321249 seconds; its scope includes checkpoints and later qualification work, and excludes startup and initial qualification. Native API timings include the observable API and compilation where it occurs; they are not isolated Lanczos timings. Compilation was not measured separately. There is no reference timing or speedup claim.

The attempt used a one-CPU limit and a memory limit of 3,499,999,232 bytes. The aggregate unit peak was 758,894,592 bytes (723.73828125 MiB), with no observed OOM event. This is an aggregate peak, not per-process or per-case RSS.

`scientific_cases.csv` and `component_timings.csv` preserve the audited reported scalar values. `results.json` adds the saved spatial observables, model matrix, conventions, scoped timings and audit summary. `native_criteria.json` distinguishes original native criteria from saved-array diagnostics. `sources.json` records content hashes of logical frozen source members; this scalar projection does not bundle their implementation or execution environment. It is not a complete reproduction package.

Run `python3 -B check_public_packet.py` to check this packet's file hashes, exact case coverage, and CSV/JSON scalar consistency. That check performs no model construction, eigensolve, or scientific comparison. See [AUDIT_SCOPE.md](AUDIT_SCOPE.md) for the remaining scientific work and timing definitions.

The 100-idea campaign ledger remains **50 PASS / 0 FAIL / 2 UNRESOLVED / 48 NOT_RUN**. This pilot adds zero ideas and changes no ledger result.
