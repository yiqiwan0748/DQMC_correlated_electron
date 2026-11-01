# Analytical two-orbital tranche: A15, A16 and A17

This is a new, not-yet-executed release for three distinct questions in the exact100 catalogue. A25 remains NOT_RUN. Parameter points, spin sectors, momentum meshes and scalar components do not increase the number of questions. Earlier releases, source repairs and receipts are unchanged. The CV is outside this release.

The four files in `source_snapshot/` are byte-for-byte existing user production source recovered from the frozen baseline archive. `source_provenance.json` records archive/member hashes. The target calls those functions directly. The separately written Python oracle imports none of them. New measurement adapters, rather than new native package capabilities, provide full spectral projectors and the cross-orbital thermal observables. These tests do not run a DQMC Markov chain or establish a bulk ordered phase.

## A15: interaction-free one-electron states inside interacting ED

The fixture is the open2x2 two-orbital cluster with the existing ALM hoppings(-1,-1.75,-0.85,-0.65), lattice spacing1, and no chemical potential or magnetic field. The primary model's intraorbital interaction and hopping convention are specified in [Giuli et al., author preprintv1, Eqs1–2](https://arxiv.org/pdf/2410.00909v1), associated with [PRB111,L020401](https://doi.org/10.1103/PhysRevB.111.L020401). The tests below are independently derived limits, not reproduced figures from that paper.

The actual `open_cluster_kinetic_matrix`, `sector_hamiltonian`, `dense_sector_hamiltonian` and `ground_state_observables` are exercised at U=0,4,12 in sectors(1,0) and(0,1). Each has8 states. The source uses

`Hc=K+U*sum_a(naup−1/2)(nadown−1/2)`.

Since no electron can be doubly occupied, its one-electron spectrum equals the independent one-body spectrum after subtracting the constant `U*(8/4−1/2)=1.5U`. The oracle constructs the open one-body matrix with Kronecker products, rather than copying the target's bond-enumeration loop. In x-fast spatial order define Ax=I⊗X, Ay=X⊗I, J=[[0,1],[-1,0]]. Then

`K=−Ax⊗diag(t1,t2)−Ay⊗diag(t2,t1)−t3*AxAy⊗I+t4*(J⊗J)⊗X`.

The complete sorted physical spectrum and every complete degenerate-group projector are compared. Grouping tolerance1e−8 and comparison tolerance1e−9 are frozen. Site probabilities and total orbital weights are diagonal traces of those projectors; arbitrary eigenvector phases or rotations inside a degenerate space do not matter. Native ground energy, residual, convergence and vanishing double occupancy provide additional gates. This is a canonical interacting-ED limit; J21 instead tested Bloch-band orbital residues at individual momenta.

## A16: the square-lattice van-Hove DOS and band edges

Set t1=t2=1, t3=t4=0 in the production band/DOS API. The two orbitals are identical copies of `epsilon(k)=−2(cos kx+cos ky)`. Divide the returned two-orbital per-spin DOS and counts by2, giving one state per site and spin. The continuum reference is

`rho(E)=K(1−(E/4)^2)/(2*pi^2)`, for |E|<4.

Here K uses the parameter m convention. This is the square limit of the rectangular-lattice derivation in [Piasecki, Eqs1–12](https://arxiv.org/pdf/0804.1037); [DLMF19.2](https://dlmf.nist.gov/19.2) defines the elliptic integral. The implementation uses [SciPy ellipkm1](https://docs.scipy.org/doc/scipy-1.17.0/reference/generated/scipy.special.ellipkm1.html) with p=(E/4)^2 to retain accuracy near the logarithmic singularity. The singularity sits on a bin boundary and is integrated, never sampled as a finite point value.

There are512 fixed bins on[-4,4], sigma=0.12/0.24, and meshes256,512,1024,2048. The last mesh is the preregistered held-out size. The observable matches the native code exactly: first obtain bin masses, convolve their bin densities with the sampled Gaussian kernel truncated at ceil(4sigma/binwidth), zero-extend at the band edges, then renormalize the retained[-4,4] area. An untruncated Gaussian on the whole real line is a different observable and is not silently compared here.

The oracle integrates the elliptic DOS over every bin. A second independent representation integrates the conditional one-dimensional CDF: `F2(E)=int_0^pi F1(E+2cos k)dk/pi`, where F1(s)=1/2+asin(s/2)/pi on[-2,2], extended by0/1 outside. Known endpoint changes are supplied to quadrature. The two bin-mass representations and reported quadrature errors must agree within1e−9; these are numerical checks, not rigorous quadrature bounds.

Every required curve and mesh must be present. Final full curves must agree within0.002 absolute DOS units; the last mesh change must be at most0.001. Each bare outer-band mass (|E|≥3.75, each side separately) must agree within0.001 and stabilize within0.001. Area must be1 within1e−10. The broadened peak must lie within2bins of zero and grow by at least0.01 when sigma narrows. No pointwise equality to the unbroadened singularity is asserted. This question concerns a continuum singularity and edge weight, rather than the previously tested two-orbital Fermi-pocket topology or orbital mixing.

## A17: independent interacting atoms must factorize

One cell contains two local Hubbard orbitals with no hopping, hybridization, Hund coupling or interorbital interaction. The full grand-canonical16-state space includes all9(Nup,Ndown) sectors; it is not a single fixed-spin sector. The grid is U=0,4,9; centered mu=−1.2,0,0.7; beta=0.3,2,6, giving27 fixtures.

Actual production sector Hamiltonians at zero2x2 kinetic provide all blocks. A new adapter subtracts mu_centered*N from each block, forms the full generator and computes its normalized thermal density matrix. For two orbitals, `G=Hphysical−mu_physical*N+U/2`, with `mu_physical=mu_centered+U/2`; the constant cancels in normalized thermal expectations. Target mode order is x-up,y-up,x-down,y-down. Native atomic annihilation operators build the even density and physical spin operators; `Sz=(nup−ndown)/2`, `Sx=(Splus+Sminus)/2`, `Sy=(Splus−Sminus)/(2i)`.

The oracle uses two independent local four-state Gibbs distributions in order(empty,up,down,double), with energies(U/4,−U/4−mu,−U/4−mu,U/4−2mu), and explicit local spin matrices. Their tensor product gives the physical probabilities after a documented bit-order permutation. Every connected interorbital density correlation and all9 connected spin components must vanish within1e−10. Full probabilities, density, double occupancy, local density variance and all local spin-squared components must match within1e−10. Nonzero local moments are required, preventing a zero-operator construction from passing the cross-correlation checks. The native `solve_two_orbital_atom` density, double occupancy, ALM first moment and `S_alm_q0=<(Sxz−Syz)^2>` are independently compared as well. This is a decoupling/factorization question, not another single-orbital atomic equation-of-state test.

## Pilot, full run and status rules

No numerical physics has been run locally or submitted by this preparation. Python AST/TOML parsing, Julia syntax parsing, shell syntax and verbatim-source hashes are checked without evaluating the scientific routines. Exact runtime pins are Julia1.12.1, NumPy2.4.6 and SciPy1.17.1; the Python interpreter is at least3.11. Only standard libraries are needed on the Julia side. Source hashes and Julia version are checked before source inclusion. The wrapper verifies the complete release hash list before launching.

`run.sbatch` defaults to **pilot**: A15/A17 small fixtures and only256/512 DOS meshes, no oracle or comparator, and explicitly zero scientific completions. Resources are1CPU,3500MiB,15min wall time, with a10min target timeout. The unmeasured planning estimate is1–3min for the pilot and1–10min for a full job, excluding queue delay. Actual per-mesh timings and `/usr/bin/time -v` receipts are retained. The2048 mesh has4,194,304 momenta and8,388,608 band values; one energies array occupies64MiB, while Julia startup/temporary allocations remain unmeasured.

After the pilot finishes, the parent may run the standard-library assessment:

```
python pilot_summary.py --target PILOT/target.toml --max-rss-kib MEASURED_RSS --output NEW_ASSESSMENT.json
```

Its conservative empirical target-time projection is `2*t512*sum(Nk^2)/512^2+120sec`. The projection must not exceed600sec, and measured pilot RSS must not exceed1,500,000KiB, before a separate full-stage review. These are planning gates, not guaranteed full-size resource bounds or scientific acceptance. The wrapper never submits another job. A full stage requires a separate explicit parent launch with `ANALYTIC_RUN_MODE=full`, unchanged grids and tolerances, and a new job-specific directory. Oracle timeout is180sec. No missing2048 result or failed gate is dropped or replaced by a more favorable mesh.

The comparator is standard-library postprocessing. It rejects wrong source/runtime/manifest identities, duplicates, missing grids and pilot receipts. Missing or malformed evidence, unconverged eigensystems, invalid density matrices, and failed last-mesh stability gates are UNRESOLVED; eligible numerical disagreement is FAIL. A question passes only when all its frozen components pass. Any resource limit or scientific failure is retained, and a revised protocol would require a new release. The full job can at most add3 distinct completions after independent audit; pilot results add0.

## A25 remains unsupported here

The frozen complex-twist constructor provides an LxL two-orbital torus, rather than a native one-dimensional ring. At specially decoupled hoppings a principal-block reduction could isolate a ring, but that would be a new geometry adapter requiring separate assessment. It is not implemented or counted here. The one-orbital real PBC/APBC seam also cannot supply arbitrary flux. A25 stays NOT_RUN without substituting a torus stiffness or a newly constructed target.
