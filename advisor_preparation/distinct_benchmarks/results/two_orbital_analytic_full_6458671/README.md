# ISAAC A15/A16/A17 full attempt6458671

Independent audit: **3 PASS**, A25 **NOT_RUN**. These are three distinct catalogue questions; the parameter points and796 checks do not increase that count. No new model, mesh, eigensolver or simulation was run during the audit.

The frozen comparator replay is byte-for-byte identical to `comparison.json`. `independent_audit_attempt01.json` additionally verifies28 input hashes, all14 release-archive members, four verbatim user-source files against the original baseline archive, and nine remote/local text transfers. Its exact auditor source is preserved as `independent_auditor_attempt01.py`.

| Question | Complete coverage | Metric-specific result |
|---|---|---|
| A15 one-electron limit | All six U/spin-sector fixtures, eight eigenstates and whole degeneracy projectors | Maximum physical-energy error1.1990e−14; projector error4.0523e−15; independently reconstructed H-projector residual6.6613e−16 |
| A16 continuum DOS | All four meshes256/512/1024/held-out2048; both widths0.12/0.24;512 bins each | Final maximum curve errors8.1215e−5 /2.1096e−5; preceding-mesh changes9.7064e−5 /2.4295e−5; edge-mass error1.6844e−6; edge-mass mesh change2.3842e−7 |
| A17 independent atoms | All27 U/mu/beta fixtures; nine spin components; complete16-state/nine-sector trace and native atom helper | Every independently recalculated analytic/native observable error at most2.2205e−15 |

The A16 density is per site and spin for one orbital, obtained by dividing the two identical production bands by2. The comparison uses the declared nonzero Gaussian broadening, zero extension and final area renormalization. It compares that observable with continuum elliptic-integral bin masses. The independent audit reconstructs the entire broadened curve from the raw masses with maximum error1.6653e−16 and recomputes peak/area metadata. It does not equate a finite sampled delta peak with the unsmoothed continuum singularity.

A15 subtracts the centered interaction shift1.5U at one electron. A17 uses centered mu, with physical mu=mu_centered+U/2, and the full grand-canonical trace of two independent orbitals. The fixed source has intraorbital U; no Hund/interorbital-interaction result is inferred. A25 remains unsupported in this tranche because its requested arbitrary-flux one-dimensional ring is not the native two-orbital torus geometry.

ISAAC completed with exit0 in125seconds. GNU time recorded target79.81seconds /587952KiB maximum RSS and oracle36.74seconds /70424KiB. Keep these command-level peaks distinct from the scheduler batch-step RSS. Julia1.12.1, NumPy2.4.6, SciPy1.17.1 and single Julia/BLAS threads match the frozen manifest. Source archive SHA256 is `3ac8111fcbf2b04e8201cd4602da80f55b94b672ff335562942671cb3c06d983`.

The original source archive has logical input name `baseline-source.tar.gz` in the audit receipt; its public workspace location is `advisor_preparation/validation_source/baseline-source.tar.gz`, SHA256 `11c3c33adaa70a78b8ab10437fba04651e3951c040bb289c4a3ce7ec293420aa`. The frozen `source_provenance.json` records its relative location and exact four member paths. Every failed/pilot predecessor remains separate evidence. This audit changes no canonical catalogue status or publication itself.
