# Journal-anchored tranche 01

Frozen before execution on 22 October 2026. Two distinct scientific questions,
J20 and J21, use actual verbatim production band/real-space source from the
user's AltermagnetDQMCTrial package. Scalar momenta and mesh refinements are
components of these two questions, not additional benchmark ideas.

The primary reference is Raghu et al., *Physical Review B* 77, 220503 (2008),
DOI 10.1103/PhysRevB.77.220503, with the explicit accessible reference version
[arXiv 0804.1113v1](https://arxiv.org/pdf/0804.1113v1). Its Eqs. 2–4 and Fig. 3 use
t1=-1,t2=1.3,t3=t4=-0.85 and mu=1.45. The final APS PDF was inaccessible during
preparation. This run tests the stated author-preprint version, and does not
silently substitute the commonly used mu=1.54 or the user's different ALM
default hoppings. No digitized numerical curve is claimed.

J20 checks two hole pockets in the lower band at Gamma and (pi,pi), and two
electron pockets in the upper band at (pi,0) and (0,pi), in the unfolded one-Fe
Brillouin zone. Periodic connected components are measured on 128,256,512,1024
square meshes. Every mesh must have the published connectivity and centers
within one mesh spacing. Target band energies/filling must match a separately
written evaluation of the published equations within 1e-10/1e-12. The last two
fillings must differ by at most 1e-3 per cell including both spins. Sixty-four
radial contour crossings must agree within 1e-9. A separate Fourier projection
of the user's real-space kinetic matrix on an 8x8 torus must reproduce the
published Bloch matrix within 1e-10. These tests do not run interacting DQMC.

J21 checks orbital-resolved band residues along two registered paths crossing
avoided-crossing regions. Numerical eigenvectors of the user's production
Bloch matrix are compared with the independently evaluated 2x2 spectral
projectors from the paper equations. Every residue must agree within 1e-8;
normalization error must be at most 1e-12. Degenerate points are excluded by a
fixed 1e-8 splitting guard because individual degenerate-band residues are
gauge-dependent. Both sides of each crossing must exchange dominant orbital
character. This is an equation-derived scientific benchmark, not a claim that
the paper publishes these individual residue samples.

Each question receives PASS only if all its registered components pass. Missing
crossings/execution errors remain UNRESOLVED; numerical disagreement is FAIL.
No thresholds or parameters change after seeing results. All source hashes,
measurements and failed outcomes are retained. Estimated runtime is 1–3 minutes
on one ISAAC CPU after loading, separate from scheduler delays; the first
execution supplies the measured runtime for later planning.

The source folder `../frozen_alm` contains unchanged user files recovered from
the already preserved baseline archive. The driver includes those actual files
inside an isolated module. Its independent analytic expressions import none of
the target helpers. This tests production Hamiltonian/band construction, not an
interacting solver or a many-body phase claim.
