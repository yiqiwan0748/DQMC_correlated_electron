# Scope of the saved audits

The saved native audit passed 158 original native-only criterion checks and 285 saved-array consistency checks. A separate audit passed 2,405 runtime, record-integrity, resource and lifecycle checks. These are check instances within one pilot, not completed benchmark ideas. The frozen model, solver, observer and acceptance criteria were unchanged.

The native checks cover the requested field/sector coverage, convergence and method, normalization, reported and action residuals, extensive O and O² closure, field-energy sign and normalization, spatial sums, and the registered balanced-sector zero-field and spin-flip controls. The original residual ceiling is 1e−11, norm-error ceiling 1e−12, O/field/spatial closure ceiling 1e−10, balanced zero-field spin ceiling 1e−10, and balanced spin-flip ceiling 1e−9. The complete per-criterion summary is in `native_criteria.json`.

Independent read-only saved-array checks covered the physical basis and particle counts, recorded kinetic coefficients and fermionic signs, centered U and −hO diagonal, trace and Hermiticity, matrix-vector residual and Rayleigh consistency, spatial observers and energy decomposition. No new target or reference model, eigensolve, or native API call was made by those saved audits. Their floating-point diagnostics are not outward-rounded uncertainty certificates.

The runtime and integrity audit joined the saved source/runtime identities, initial output-file identities, all 15 completed records and actual API-dispatch witnesses, resource limits, and terminal cleanup. This supports a source-conditioned interpretation of the saved results; it is not an instruction-level attestation.

The independent QuSpin construction and comparison remain NOT_RUN. So do the global ground-state/gap/fidelity certificates, extended-precision uncertainty work, complete spectral and independent-resolvent response, and full nine-field study. Full X25 needs all nine fixed fields and all five sectors, the registered central and Richardson stencils, held-out and stress controls, energy curvature, and original error bounds. This three-field native record cannot replace them. No scientific completion, full reproduction, external comparison, or speedup follows from the pilot audit status.

## Timing definitions

All cases ran sequentially in one attempt. `hamiltonian_seconds` measures sector construction; `dense_materialization_seconds` measures dense conversion. `native_ground_api_seconds` includes the full native ground-state-observable API and compilation where it occurs. `spatial_observer_seconds` measures the copied spatial observer and can also include compilation.

`case_elapsed_excluding_final_serialization_seconds` ends before completion-record encoding, flushing, closing and readback. That work is recorded separately as `serialize_flush_close_readback_seconds`. Case elapsed includes the started record and COO/basis/record work; the subsequent per-case scope check lies outside it. The components are not a complete disjoint partition of wall time.

Producer elapsed begins after the initial runtime-before hook, before model construction. It includes model creation, the cases, checkpoints, per-case hooks and the runtime-after hook. The whole native Julia entry includes startup, definitions and qualification. Syntax-only parsing is separately recorded. Whole external time includes transport and surrounding execution work. Compilation has no separately measured duration, so first-call and later-call values must retain their distinct contexts.

## Data and validation limits

The CSV and JSON values are an explicit scientific projection of saved records. No matrix/vector eigensolve is performed by the public packet checker. `MANIFEST.json` hashes every other packet member; its own hash is intentionally external to that manifest. The checker accepts only the nine listed regular files, with a 1 MiB limit per file and 4 MiB limit for the whole packet, and parses the same bounded bytes that it hashes. It verifies bytes and scalar consistency, not the private raw evidence or the original scientific audits. Unknown or unsupported future work remains NOT_RUN.
