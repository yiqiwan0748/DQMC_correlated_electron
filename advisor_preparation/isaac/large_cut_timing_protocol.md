# Follow-up subsystem-size timing study

Registered before this follow-up runs, after seeing the first half-chain timings
on 22 October 2026. This is a transparent follow-up to the original size grid,
not an originally preregistered result. All original timing rows remain visible.

Hold the chain at 64 sites and compare subsystem sizes48 and64, each in a new
receipt. Keep U/t=2, beta*t=2, dtau*t=0.05, PBC, lambda=0.5, stabilization5,
the same two seeds, identical common warmed states and matched direct/optimized
paths, seven alternating-order repeats and12 two-replica rounds. Both variants
use the same corrected50% local idle proposals. Repeat and correctness rules
are unchanged. Report each cut, elapsed times and allocations separately.

The larger subsystem increases the Grover-matrix dimension at fixed propagation
cost. This asks where the optimization is most useful. A48-site cut remains a
nontrivial spatial subsystem; a64-site cut is the **full-system finite-temperature
second Renyi purity**, not an entanglement entropy of a pure-state bipartition.
Both are legitimate Renyi computations, but their scopes must remain explicit.

These are node-loop throughput timings, including replica updates and matching
measurements. They exclude initialization, common warmup, compilation, external
I/O and the complete lambda grid. Matching trajectories checks computational
equivalence and does not establish that the short timing run is equilibrated.
The source's existing direct determinant path is the baseline; no deliberately
inefficient replacement is introduced.

Estimated sampling runtime from the completed half-chain64-site measurement is
roughly2–8minutes for both cuts together, plus loading/queue delay. The batch has
a30-minute safety limit, and the benchmark retains partial repeats at its fixed
soft budget. No cut is repeated with new seeds or discarded for a small speedup.
