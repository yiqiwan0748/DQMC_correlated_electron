# Independent benchmark references and conventions

The executable campaign contains 100 observable/model/ensemble comparisons, grouped by their independent oracle families. Repeating a parameter point or measuring several observables does not create a new independent method. An exported oracle value is a reference, not evidence that the user package passed.

## QuSpin: 80 finite interacting-chain comparisons

The exporter builds the Hamiltonian and fermionic reduced density matrix with the published QuSpin package, without importing user code. Its basis and operator conventions follow the [official spinful-fermion basis documentation](https://quspin.github.io/QuSpin/generated/quspin.basis.spinful_fermion_basis_1d.html) and [Hubbard example](https://quspin.github.io/QuSpin/examples/example6.html); the package paper is [Weinberg and Bukov, SciPost Physics 7, 020 (2019)](https://arxiv.org/abs/1804.06782).

Eight fixed-sector ground-state models provide 64 comparisons of energy, double occupancy, two longitudinal correlations, two one-body correlations and two subsystem Rényi entropies. Two grand-canonical finite-dtau models provide another 16 comparisons. Model parameters, zero-based site indices and observable normalization are explicit in `cases.json`.

Ground-state ED and DMRG must use the same fixed Nup/Ndown, boundary seam, hopping sign and interaction. Pure-state observables require a unique ground state; the exporter withholds a degenerate target. DQMC comparisons instead use the full grand-canonical trace at the same finite beta and dtau. Their external density is proportional to `[exp(-dtau V) exp(-dtau K)]^M`. This asymmetric transfer density must not be replaced by `exp(-beta H)` or silently symmetrized before taking its subsystem trace.

## Analytic free Gaussian references: 12 comparisons

For U=0, diagonalize the one-body hopping matrix independently. With eigenpairs `(epsilon_a,v_a)`, the thermal occupation matrix is `C_ij=sum_a v_ia conjugate(v_ja) f(beta*(epsilon_a-mu))`. The physical kinetic energy sums `epsilon_a f_a` over both spins. If `nu_a` are subsystem eigenvalues of C for one spin, then `S2=-2 sum_a log(nu_a^2+(1-nu_a)^2)` in a spin-balanced state. These formulas are an application of [Peschel's correlation-matrix construction](https://arxiv.org/abs/cond-mat/0212631). Use both spin species, grand-canonical thermal occupations and natural logarithms.

## Analytic Hubbard atom: 8 comparisons

This oracle is a direct four-state calculation, derived here rather than read from user implementation. For `Hc=U(nup-1/2)(ndown-1/2)-mu_c N`, the empty, two singly occupied and doubly occupied states have energies `U/4`, `-U/4-mu_c` (twice), and `U/4-2mu_c`.

Write their Boltzmann weights as `w0,w1,w1,w2`, and `Z=w0+2w1+w2`. Then `N=(2w1+2w2)/Z`, `D=w2/Z`, and the centered physical energy excluding the chemical-potential term is `Ecenter=U*(w0+w2-2w1)/(4Z)`. The atom Rényi entropy is `-log((w0^2+2w1^2+w2^2)/Z^2)`. For independent atoms, subsystem S2 adds across retained sites. A canonical fixed-N atom has a different target and must not be substituted.

## Analytic dimer guard for the exporter

Before exporting any value, the Python program checks its open, two-site, Nup=Ndown=1 Hamiltonian against `E0=(U-sqrt(U^2+16t^2))/2` and total double occupancy `(1-U/sqrt(U^2+16t^2))/2`, using uncentered physical energy. The energy formula is given in [Carrascal et al., Eq. 36](https://arxiv.org/html/1502.02194); the double occupancy follows by differentiating that energy with respect to U. This plumbing guard is not counted among the 100 user-package comparisons.

## Published curves for a later paper-reproduction stage

The current external package/analytic campaign can be extended with these directly relevant journal targets. They are proposed comparisons, not completed results:

- [Broecker and Trebst (2014), Fig. 9](https://arxiv.org/pdf/1404.3027): periodic half-filled L10 Hubbard chain, beta=40, U/t=1,2,4; reproduce S2 versus all contiguous block sizes. The plotted ED references use beta=infinity, so finite-temperature and Trotter effects need separate control.
- [Broecker and Trebst, Figs. 5–7](https://arxiv.org/pdf/1404.3027): ground-state L64 U/t=2 subsystem curve, finite-temperature L32 U/t=2 crossover, and the high-temperature half-chain limit S2/L=log(2). Use a figure-specific digitization allowance where source ordinates are unavailable.
- [Lieb and Wu (1968)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.20.1445), with their [detailed derivation](https://arxiv.org/pdf/cond-mat/0207529): thermodynamic half-filled 1D energy and charge gap; finite chains require a size extrapolation. The uncentered energy per site is `-4t integral_0^infinity J0(x)J1(x)/(x*(1+exp(U*x/(2t)))) dx`. Double occupancy is its U derivative.
- [LeBlanc et al. (2015), Tables II–III](https://arxiv.org/pdf/1505.02290): square-lattice thermodynamic-limit half-filled energy and double occupancy for U/t=2,4,6,8,12. These require a two-dimensional, zero-temperature, t'=0 calculation and a controlled size extrapolation. They must not be used as finite-cluster pass/fail thresholds.
- [TeNPy's FermiHubbardModel](https://tenpy.readthedocs.io/en/latest/reference/tenpy.models.hubbard.FermiHubbardModel.html): an additional independently implemented tensor-network comparison for matched finite chains. If the user's solver is ITensor-based, a new ITensor call alone does not supply implementation independence for the DMRG algorithm.
- [ALF](https://arxiv.org/abs/2012.11914) supplies a separately implemented finite-temperature AFQMC comparison. [SmoQyDQMC](https://arxiv.org/abs/2311.09395) is useful too, but shares JDQMC components with this project and should not be described as wholly implementation-independent for those shared kernels.
