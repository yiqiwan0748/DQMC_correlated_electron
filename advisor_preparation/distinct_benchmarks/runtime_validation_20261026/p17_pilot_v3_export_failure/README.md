# P17 excluded pilot: strict TeNPy JSON export failure

The excluded six-site Hubbard boundary-spin pilot stopped at the independent TeNPy producer's final JSON export. The scheduler recorded FAILED (exit 2:0). The native producer completed six samples; those samples have no completed independent comparison and establish no scientific PASS. TeNPy rejected a NaN under `allow_nan=False` before writing its final report. The QuSpin reference and comparator were not run.

The pilot used six particles, open boundaries, hopping t=1, onsite U=4, and edge fields h=0 and ±0.01 coupled as −h Sz, with Sz=(Nup−Ndn)/2. It used two prescribed starting states and the refined64 arm. The frozen source and failure receipts are preserved unchanged. Native and TeNPy processes took 263.1501754075289 s and 155.2780849058181 s; final workflow elapsed time was 583.1100067906082 s. These are failed-pilot timings, not a speed comparison.

The pinned TeNPy 1.1.1 source initializes the first previous-sweep energy to NaN and appends the resulting undefined first energy change to `Delta_E`. That source behavior motivates a narrow export repair. The original in-memory report was not persisted, so the exact failing NaN leaf and completed TeNPy sample count remain unknown. This record does not claim to have recovered them.

The strict export error must remain distinct from a physical-data verdict. Any subsequent attempt needs its own reviewed source and preserved failure history. A proposed encoding for an undefined initial diagnostic must retain its nonfinite identity and must keep nonfinite physical results rejected; no tolerance change or successful recovery is claimed here.

`result.json` contains the public projection and hashes binding the preserved private original archive, execution receipt, native producer output and frozen source. Private operational records are not included. The original terminal evidence is not rewritten by this projection.

The catalogue remains **49 PASS / 0 FAIL / 2 UNRESOLVED / 49 NOT_RUN**. This excluded pilot adds **zero** distinct completions, and P17's full seven-field comparison remains **NOT_RUN**. The ledger is unchanged.

An independent 210-check source/hash/failure-boundary audit passed (SHA256 `393cfd930e856d798ce512e26de43d601c76fdf4e9c913313f267ae106246d09`). Its scope is receipt integrity and failure boundaries, with zero scientific credit. Runtime status is preserved as reported and hash-bound; this projection does not claim a new full runtime-provider audit or a successful final workflow source postcheck.
