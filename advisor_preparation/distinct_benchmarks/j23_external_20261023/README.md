# Prescribed mean-field trial functional

J23 compares the user's actual `HFCompetition.evaluate` / `solve_fields` and
Bloch APIs with independently derived projectors, the physical Wick functional
and Gaussian quadrature. It tests the orbital-antialignment mechanism described
in [Appendix D of the primary paper](https://arxiv.org/html/2410.00909v1#A4).
This is one specified variational-functional benchmark at density n=2. It does
not establish a self-consistent solution, a gRISB or doped phase diagram, an
exact ground state, or reliability of the whole package.

All 112 native rows and 196 reference rows, including 84 Gaussian-integration
rows, are retained. The fixed meshes are 16/32/64/128, quadrature orders
64/96/128, fields m=12/16 with U=2m, seven registered hybridizations and both
parallel/orbital-antiparallel trial inputs. Every actual returned field,
chemical potential, density and pure occupation is qualified before the
physical kinetic-plus-interaction energy is interpreted. The 3,542 original
gates include all primary and held-out Richardson, mesh and quadrature tests.

The limiting coefficient of the energy difference divided by t4² is:

| m | Independent coefficient | Native primary Richardson |
|---:|---:|---:|
| 12 | -0.3339875758251529 | -0.33398757803689466 |
| 16 | -0.2502754154885741 | -0.2502754216404668 |

The maximum coefficient discrepancy is about 6.16e-9, below the registered
1e-6 gate. These are physical trial-functional differences, rather than a sum
of auxiliary single-particle eigenvalues. The complete arrays and every gate
are available in the [original comparison](production/original_comparison.json).

The full ISAAC run is 6466733. The independent original saved-data audit
preserves 646 source/metadata checks and 2,842 independent physical checks,
including closed 2x2 projectors, the Wick functional, a separate exact-Fraction
cosine-moment coefficient enclosure, and reconstruction of all original gates.
Its [unchanged report](production/original_independent_audit.json) is bound to
the [original auditor](independent/audit_production_original.py). The actual
pinned-runtime comparator replay 6466878 is byte-identical to the original
comparison. The excluded cold/warm pilot 6466377, its 599 checks and its
planning assessment remain [separate](pilot/original_independent_audit.json),
with zero catalogue credit. Its analysis-only replay is 6466655.

`sources/` contains the exact 19-member scientific release and both original
saved-analysis source releases. `pilot/data.tar.gz` and
`production/data.tar.gz` retain complete scientific text arrays, own execution,
resource logs and own-job accounting with per-member SHA256/size inventories.
No success subset is selected. The full decoded payload is approximately
365 MB and compresses to approximately 45.5 MB. Extraction rejects absolute or
parent paths, links, duplicate names, unexpected members, changed bytes and
non-UTF8/runtime binary data. Original source/freezes are not rewritten to
describe later execution; their initial preparation statuses remain historical.

The public selection excludes private transport/controller/admission bodies,
unrelated queues, CV/contact material and runtime binaries. Their omitted
original hashes remain provenance anchors. The two replay execution documents
are explicitly named **derived scientific projections**: their retained fields
are unchanged, their original receipt hashes are recorded, and their before/
after transport bodies are omitted. The original full/pilot audit reports stay
byte-exact and retain authority for those original private checks.

Replay the complete public scientific math with Python 3.11 or later and only
the standard library, from this directory:

```sh
python3 replay.py --role pilot --output-dir /tmp/j23-pilot-public-replay
python3 replay.py --role production --output-dir /tmp/j23-full-public-replay
```

Use fresh output directories. This safely reconstructs every selected source/
data byte and executes the **unchanged complete scientific main fragments** of
the original independent auditors. It rederives the full projectors, Wick
energies, exact-moment coefficient, all fixed gates and pilot resource math.
It verifies the retained actual pinned compute replay; it does not pretend to
rerun that frozen comparator or the remote transport checks on the local PC.
The full raw arrays can require several GB of memory during parsing.

To rerun the original frozen comparator itself, extract the packet, use the
actual pinned Python/39-library environment in a separately approved ISAAC
compute allocation, and execute `source/compare.py` on the extracted original
native/reference files to a fresh output. The original saved-analysis source
archives retain the exact bounded runner. `--frozen-replay NEW_COMPARISON` can
check such a new comparison, but only bytes matching the original historical
replay receive the original provenance flag. No renderer, evaluator,
quadrature oracle, eigensolver, sampler, job submission or full admission is
performed by this public replay.
