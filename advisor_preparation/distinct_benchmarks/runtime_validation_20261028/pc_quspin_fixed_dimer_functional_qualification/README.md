# Fixed Hubbard dimer: QuSpin functional check

A single fixed-case QuSpin reference-process check completed and passed its saved-record audit. This packet preserves the scientific routine, frozen protocol, complete finite output, and recorded gate values.

The model is the open two-site Hubbard dimer,

\[
H=-\sum_{\sigma}(c^\dagger_{0\sigma}c_{1\sigma}+c^\dagger_{1\sigma}c_{0\sigma})
  +4\sum_{i=0}^{1}n_{i\uparrow}n_{i\downarrow},
\]

with zero chemical potential and energy offset, uncentered onsite interaction, and fixed particle numbers \(N_\uparrow=N_\downarrow=1\). The sector has dimension four. Its registered analytic levels are \((4-\sqrt{32})/2,0,4,(4+\sqrt{32})/2\).

| Sorted level | Saved energy |
| --- | ---: |
| 0 | -0.8284271247461863 |
| 1 | 0.0 |
| 2 | 4.000000000000001 |
| 3 | 4.82842712474619 |

All four recorded conservative absolute-energy bounds passed the fixed limit of \(10^{-11}\). The recorded stored-matrix Hermitian defect is zero, within its \(10^{-12}\) limit. Exact rational bounds and the complete saved arrays appear in `scientific_data.json`.

This checks QuSpin's constructor, number-operator, dense-matrix, and complete-eigensolver path for this case. It does not test the user's ED, DQMC, or DMRG implementation, complete full X08, or qualify another model or system size. New ideas: **0**. Catalogue count: **unchanged**.

`dimer_science.py` and `DIMER_PROTOCOL.json` are unchanged scientific-registration sources used for the recorded check. Their generic admission annotations and future-tense comments describe the original pre-data registration; they are not the outcome of the completed check. The routine accepts qualified package exports and bookkeeping callbacks. It is not a standalone rerun launcher.

`saved_scientific_record.json` preserves the audited scientific checkpoint byte for byte. `scientific_data.json` presents its physical fields with a public schema and explicit scope. `PROJECTION_AUDIT.json` records the field mapping and scientific-file hashes. `MANIFEST.json` inventories the packet.

The packet has no dependency lock, machine-qualification code, or physical execution driver. It documents one audited result. Its byte checker does not repeat the QuSpin calculation or recreate the process qualification.

Run `python3 verify_packet.py` to check file hashes, the projection, finite array shapes, and recorded rational bounds. This checker uses the Python standard library and performs no package imports, reference generation, diagonalization, or residual/Gram calculation.

See `AUDIT_SCOPE.md` for the qualification boundary and the two defects corrected before numerical execution.
