# Altermagnetism in the two orbital Hubbard model

Pedagogical project note. Updated October 21, 2026. Live calculation status was checked at 1:46 p.m. EDT.

The question is whether intraorbital Hubbard repulsion produces altermagnetism (ALM) in this two-orbital model. The candidate pattern has opposite orbital magnetizations and zero total magnetization, together with the appropriate spatial and spin symmetry. Our restricted mean-field calculation supports such a solution. Our finite-cluster ED calculations have verified several eigenstates and observables, but have not established an ALM phase.

## Model and parameters

Each square-lattice cell contains two orbitals, ordered as `(d_xz, d_yz)`. The physical Hamiltonian is

```text
H = Σ(k,σ) c†(k,σ) h(k) c(k,σ)
    + U Σ(i,a) [n(i,a,↑) − 1/2][n(i,a,↓) − 1/2]

h11(k) = −2t1 cos(kx) − 2t2 cos(ky) − 4t3 cos(kx)cos(ky)
h22(k) = −2t2 cos(kx) − 2t1 cos(ky) − 4t3 cos(kx)cos(ky)
h12(k) = h21(k) = −4t4 sin(kx)sin(ky)
```

Here `c(k,σ)` is the two-component orbital annihilation operator. The hopping matrix is real and symmetric. There is no applied field, spin-orbit coupling, Hund interaction or interorbital density interaction.

| Setting | Value |
|---|---|
| Hoppings `(t1,t2,t3,t4)` | `(−1, −1.75, −0.85, −0.65)` |
| Units | `abs(t1)=1`, lattice spacing `1`, spin unit `hbar=1` |
| Main interaction | `U=20`; requested continuation to `U=100,1000` |
| Additional controls | `U=0` and `U=12` |
| ED particle numbers | `N↑=N↓=6`, total `N=12` |
| Density | `12/8=1.5` electrons per cell, counting both orbitals and spins |
| ED state prescription | A normalized pure eigenstate; no thermal average |
| Current chosen branch | The verified excited X singlet, followed as U increases |

The centered interaction convention matters for energies. With 16 spatial orbitals and N=12,

```text
E_centered = E_uncentered − U*N/2 + U*16/4
           = E_uncentered − 2U.
```

At U=20, add 40 to a reported centered cluster energy to obtain the energy for `T + U Σ n↑n↓`.

## Lattice and the meaning of momentum

ED uses an eight-cell tilted periodic cluster, with boundary vectors `a1=(2,2)` and `a2=(−2,2)`. It has 16 spatial orbitals and 32 spin modes. Its eight allowed momenta include Γ=(0,0), X=(π,0), Y=(0,π), M=(π,π), and the four combinations `(±π/2,±π/2)`, modulo 2π.

The cell ordering used in our configuration pictures is

```text
Cell   C1     C2       C3      C4       C5      C6     C7     C8
(x,y) (0,0) (−1,−1) (0,−1) (1,−1) (−1,0) (1,0) (0,1) (0,−2)
```

Periodic translations can move these representatives to other equivalent positions in a drawing.

**Total many-body momentum K** labels how the whole 12-electron wavefunction transforms under translation. **Electron momentum k** labels the electron added or removed in a spectral measurement. A many-body state with K=X contains many electron arrangements; it does not place all electrons at X. A correlation momentum q is a third label: it describes the spatial pattern of a measured correlation.

The two eigenvalues of the 2×2 matrix `h(k)` are bare single-electron bands. At finite density and U=20, the interacting electron spectrum requires the many-body calculation or a stated approximation such as mean field.

Three bases also need to stay distinct. The occupation basis lists electron arrangements. The momentum basis combines translated arrangements with momentum phases. An energy eigenvector is a superposition of those momentum-basis vectors. Our configuration tables report literal occupation amplitudes, including the conversion from the momentum basis. In its own energy eigenbasis, a chosen eigenstate would simply have one unit coefficient.

## Mean field and the reference paper

The reference is Giuli, Mejuto-Zaera and Capone, *Altermagnetism from interaction-driven itinerant magnetism*, PRB 111, L020401 (2025). It uses the same hopping parameters and intraorbital interaction. Its Figure 2 phase diagram uses **gRISB**, a correlated slave-boson method with auxiliary ghost orbitals. Figure 3 compares HF, RISB and gRISB ALM spectra at U=20 and the paper's δ=−0.15.

We report our density explicitly as 1.5 electrons per cell. The supplied main paper does not explicitly define the normalization of δ, so we do not assume that its δ=−0.15 panel is exactly our filling.

Our own completed calculation is restricted Hartree–Fock (HF) on an infinite periodic lattice, using a 384×384 momentum mesh and temperature T=0.0125. It compares uniform collinear AM and FM solutions with equal orbital densities, plus a paramagnetic (PM) reference. Its self-consistent magnetic field represents the chosen mean-field order; it is not an applied field in the physical Hamiltonian.

Define the number polarization per orbital and cell by

```text
m_a = n(a,↑) − n(a,↓)
Δ_ALM = (m_xz − m_yz)/2.
```

The physical orbital spin per cell is `m_a/2`. Our AM solution has

```text
n_xz = n_yz = 0.75
m_xz = −m_yz = 0.737859802195
n(xz,↑) = 0.743929901097; n(xz,↓) = 0.006070098903
```

The yz populations are reversed. Thus the net magnetization cancels while Δ_ALM is nonzero.

| HF solution | Centered free energy per cell |
|---|---:|
| AM | −7.676642768677 |
| FM | −7.638249439992 |
| PM | −4.276642562362 |

AM is lower than FM by 0.038393328685 per cell **within the tested uniform families**. The calculation passed 109 numerical checks. It does not establish a global HF minimum: mixed orders, finite-q orders and phase separation were not searched. There is also no strict T→0 extrapolation. Its infinite-lattice finite-temperature free energies should not be directly equated with finite-cluster ED energies.

The HF AM bands have opposite-spin branches exchanged between X and Y. That result belongs to a mean-field state with a selected magnetic orientation. It is not automatically a property of an exact finite-cluster singlet.

## Verified ED results

ED here uses matrix-free Krylov solvers rather than a dense diagonalization of the whole Hamiltonian. Each accepted state is checked with fresh applications of the original H and full S². A small residual verifies an eigenpair; it does not prove that no lower eigenstate exists.

All interacting states in the table have N↑=N↓=6. Energies are centered cluster energies, rounded for readability.

| U | Total K | S | Energy | State |
|---:|---|---:|---:|---|
| 12 | Γ | 0 | −52.635500032 | Selected singlet candidate |
| 12 | M | 4 | −52.659001054 | Lowest among seven calculated records |
| 20 | M | 4 | −66.689705180 | Verified selected M state |
| 20 | X | 3 | −66.438426619 | Earlier verified X state |
| 20 | X | 0 | −65.598249718 | Current chosen singlet |
| 20 | Y | 0 | −65.598249718 | Constructed by rotating that X singlet |

The lower S=3 X state proves that our selected U20 singlet is excited. Its lowest-singlet status and the global ground state remain unestablished. A finite-S excited state alone does not establish a ferromagnetic phase.

The current X singlet has original-H residual about `1.14×10⁻¹⁰` and full-S² residual about `3.99×10⁻¹⁴`, within the existing `10⁻⁹` and `10⁻⁶` gates. Norm and relevant Gram gates are `10⁻¹⁰` and `10⁻⁹`. The rotated Y passed corresponding fresh checks; it was not independently diagonalized.

Older magnetic correlations are also accepted. For physical spin operators, define

```text
D = Σ_i (S(i,xz) − S(i,yz))
S_L(0) = <D·D>/Nc,                  Nc=8
S_M(0) = <S_total²>/Nc = S(S+1)/8
C12(r=0) = (1/Nc) Σ_i <S(i,xz)·S(i,yz)>.
```

Here r=0 means the same cell, whereas the zero in S_L(0) and S_M(0) means momentum q=0. Positive C12 means onsite parallel spin correlation; negative C12 means onsite antialignment.

| Selected static state | S_M(0) | S_L(0) | C12(r=0) |
|---|---:|---:|---:|
| U0 exact closed-shell singlet | 0 | 0.750000 | Not listed here |
| U12 Γ singlet | 0 | 0.538110 | −0.083288 |
| U12 M, S=4 | 2.5 | 0.596907 | +0.061273 |
| U20 M, S=4 | 2.5 | 0.581035 | +0.058648 |

These are different selected states, not a controlled scan along the current singlet branch. Local antialignment occurs in the U12 Γ candidate, but that does not establish long-range ALM. The older U20 M state instead has positive onsite interorbital correlation.

A useful counterexample is U0: **S_M(0)=0 and S_L(0)>0 already occur without interactions.** Those two finite-cluster conditions alone cannot prove ALM.

## What S=0 means for coefficients and spectra

Equal up/down particle numbers fix `Sz=(N↑−N↓)/2=0`. They do not fix S. An S=3 multiplet has seven members, including Sz=0. This explains why the earlier balanced calculation could produce S=3.

An exact singlet is invariant under every global spin rotation. Let F flip all spins through a unitary π rotation; it does not move cells or orbitals and does not complex-conjugate coefficients. In our occupation basis, all up creators precede all down creators. For six of each, the six rotation minus signs and 36 fermion swaps give total sign +1. Therefore

```text
ψ_K = Σ_C a_K(C)|C>
Fψ_K = ψ_K  implies  a_K(FC) = a_K(C).
```

This holds separately in X and Y, for ground or excited singlets. With different basis phases the signed relation can change, but equal probabilities remain. For definite S and Sz=0, our convention gives `a_K(FC)=(−1)^S a_K(C)`: S=3 has opposite signs but equal probabilities. Coefficient equality by itself does not identify a singlet either; even S can give the same sign.

An exact singlet has zero mean spin polarization in each orbital, yet can have strong magnetic correlations. Its spin-up and spin-down electron spectra at the same k are equal. The older S=4, Sz=0 spin-flip eigenstate also has equal diagonal spin spectra. Equal spectra alone do not determine S or establish the absence of magnetic correlations.

## Why X and Y coefficients can differ

A spatial 90° rotation R sends `(x,y)` to `(−y,x)` and transforms `dxz→dyz`, `dyz→−dxz`; it leaves spin labels unchanged. It preserves the torus and Hamiltonian, and sends total momentum X to Y. Thus

```text
ψ_Y = Rψ_X
R|C> = s_C |RC>
a_Y(RC) = s_C a_X(C).
```

The sign s_C includes orbital and fermion reordering signs. In the reduced momentum basis, translation phases are included too. This relation compares a configuration with its rotated arrangement, rather than imposing equality for the unchanged arrangement in X and Y.

Moreo's examples, in C1–C8 order, are

```text
A: dxz  ↑ ↑ ↑ ↑ 0 ↓ ↓ 0      B: dxz  ↑ ↑ 0 ↑ ↑ 0 ↓ ↓
   dyz  ↓ ↓ ↓ ↓ 0 ↑ ↑ 0         dyz  0 ↓ ↓ 0 ↓ ↓ ↑ ↑
```

Each has opposite orbital number polarizations `(+0.25,−0.25)` per cell. They are individual arrangements, not uniformly polarized orbitals or the complete wavefunction. Their occurrence does not imply nonzero orbital moments in a singlet.

| Coefficient ×10⁻³ | X original | Y unchanged | Y rotated |
|---|---:|---:|---:|
| A | −0.431544 | +0.614352 | −0.431544 |
| B | −2.763094 | −0.079250 | −2.763094 |

These are amplitudes; the configuration probability is the squared magnitude. Cross-state signed coefficients use the fixed phase convention ψ_Y=Rψ_X. A separately solved Y state would need phase alignment before comparing signs.

The Y arrangement carrying B's X coefficient is

```text
Rotated B: dxz  0 0 ↓ ↓ ↑ ↓ ↓ ↑
           dyz  ↑ ↑ ↑ ↑ ↓ 0 0 ↓
```

Its sign s_B is +1. In each singlet, the spin-flipped counterpart has the same coefficient. The proposed exchange of two unequal spin-flip coefficients between X and Y cannot occur in these exact singlets.

Our Y was constructed from X, so the rotated-coefficient match follows by construction. Fresh original-H, S² and input-vector rotation-covariance checks passed as additional tests. The covariance check is on this vector, not a bound on the full operator. C4 coefficient matching is a general symmetry check, not an ALM-specific signature.

## What the spectral calculation measures

`A_aσ(k,ω)` gives the energies and weights for adding or removing an electron of momentum k, orbital a and spin σ from the chosen many-body state. Addition poles are `E(N+1)−E0`; removal poles are `E0−E(N−1)`. This is a many-body electron spectrum, not a plot of all N=12 energy levels.

The Fourier operator uses `1/sqrt(8)` normalization. Each orbital/spin add-plus-remove pair has full-frequency weight one. We use Lorentzian broadening η=0.15 and μ=0 for the present balanced-state spectra. Both actual ω and the requested ω−U display are retained; the latter is only an axis translation.

The older U20 M, S=4, Sz=0 spectral package is complete and accepted at 0.1% normalized L1 and maximum-error targets. Its static and spectral vectors have separate provenance; matching energy and spin does not establish vector identity. Neither is the current X singlet.

The current request is to follow the excited X singlet to U=100 and 1000 and calculate electron spectra at X and Y from the same initial vector at each U. Initial total K=X gives final N±1 total K=Γ for electron k=X and M for k=Y. At the 1:46 p.m. check, the original singlet remained accepted. The comparison solve was still unconverged. A short test verified a repair for small spin leakage in its unfinished vectors, and the longer continuation was running. All 12 U20 spectral channels and a final audit are in one serial queue; each next channel requires earlier accuracy checks to pass. No new-U state or full singlet spectral package was accepted.

Increasing U cannot remove singlet spin symmetry. These spectra can reveal changes in orbital weights and excitation peaks, but cannot show spin splitting of this exact singlet. C4 also rotates the initial state: X/Y spectral differences within the same initial X vector alone are not ALM evidence.

## Main advances and what remains

The main advances so far are:

1. A matched-density restricted HF calculation supports a compensated AM solution at our actual parameters.
2. ED has exposed high-spin competitors and established a verified excited X singlet for Moreo's configuration test.
3. We separated Sz from S, total momentum K from electron momentum k, and occupation configurations from energy eigenstates.
4. The X/Y coefficient puzzle is resolved by rotating cells and orbitals. The apparent mismatch was not an inconsistent coefficient in the constructed pair.
5. The U0 correlation control showed why a nonzero orbital-difference factor on one cluster is insufficient evidence of order.
6. A selected M-state spectral package is validated; the requested singlet and higher-U spectra are a separate calculation.

Moreo's approach is to inspect candidate compensated configurations, check their symmetry relations, and then examine interaction dependence and electron spectra. Her statement that X need not be the ground state is useful for this exploratory work. An equilibrium phase claim still needs a justified connection to low-energy states.

The next magnetic test is to measure the full same-orbital and cross-orbital spin correlations of the chosen singlet. For size scaling, use `m_orb²=S_L(0)/Nc`, retaining the onsite contribution and its normalization explicitly. Keep the same filling and a consistent, stated state or ensemble prescription and cluster-geometry sequence across sizes. A positive limit as Nc grows, together with the appropriate long-range pattern and symmetry, would support order. A nonzero value at Nc=8 does not. In an exact singlet, vanishing total-spin factor is automatic and supplies no second independent ALM test.

The selected excited branch can be studied on its own terms. Establishing a zero-field equilibrium ALM phase remains open. Automatic job checks are off; this note is updated when relevant verified results or interpretations change during project work.

## Evidence and calculation details

The local links below refer to this project's workspace. They preserve the data, exact source, state prescription, normalization, reproduction commands and verification evidence.

- [Reference paper](https://doi.org/10.1103/PhysRevB.111.L020401): model Eqs. 1–4, phase diagram Fig. 2, spectra Fig. 3.
- [Our HF calculation and reproduction](/home/yiqiwan/code/AltermagnetDQMCTrial/deliverables/ed_hf_band_fs_20260928/data/hf/README.md), [numerical summary](/home/yiqiwan/code/AltermagnetDQMCTrial/deliverables/ed_hf_band_fs_20260928/data/hf/hf_summary.toml), and [109 checks](/home/yiqiwan/code/AltermagnetDQMCTrial/deliverables/ed_hf_band_fs_20260928/data/hf/verification.toml).
- [U0 and selected U20 static correlations](/home/yiqiwan/code/AltermagnetDQMCTrial/deliverables/u_only_magnetic_correlations_20260930/CALCULATION.md).
- [U12 independent-start results and seven records](/home/yiqiwan/code/AltermagnetDQMCTrial/deliverables/u12_magnetic_evidence_20260930/CALCULATION.md).
- [Accepted U20 X singlet](/home/yiqiwan/code/AltermagnetDQMCTrial/results/isaac_20261014/moreo_configurations/solve003_terminal_validation_20261015.json) and [earlier S=3 comparison](/home/yiqiwan/code/AltermagnetDQMCTrial/results/isaac_20261014/moreo_configurations/independent_validation.json).
- [Accepted coefficient and rotated-Y audit](/home/yiqiwan/code/AltermagnetDQMCTrial/results/isaac_20261014/moreo_configurations/singlet_adapter_independent_validation.json), [A/B data](/home/yiqiwan/code/AltermagnetDQMCTrial/results/isaac_20261021/moreo_followup/rotation_audit/rotation_coefficients.csv), and [geometry, phases and reproduction](/home/yiqiwan/code/AltermagnetDQMCTrial/results/isaac_20261021/moreo_followup/rotation_audit/calculation_notes.txt).
- [One-page singlet coefficient proof](/home/yiqiwan/code/AltermagnetDQMCTrial/deliverables/spin_flip_coefficient_proof_20261020/equal_spin_flip_coefficients.pdf).
- [Accepted older balanced-M spectra](/home/yiqiwan/code/AltermagnetDQMCTrial/deliverables/u20_balanced_xy_spectra_20260930/README.md).
- [Current spectral definition and numerical method](/home/yiqiwan/code/AltermagnetDQMCTrial/results/isaac_20261021/moreo_followup/SPECTRAL_METHOD.txt) and [current continuation record](/home/yiqiwan/code/AltermagnetDQMCTrial/results/isaac_20261021/moreo_followup/CONTINUATION.txt).

## Update record

October 21, 2026: First consolidated note. Added the model and eight-cell geometry, distinguished reference gRISB from our restricted HF result, collected accepted ED and magnetic results, explained spin versus spatial rotation, and recorded the chosen excited-singlet continuation. No ALM phase claim is made.

October 21, 1:46 p.m. EDT: A seven-minute repair test passed its audit and saved an unqualified checkpoint. The singlet continuation resumed with the same physical checks. The first U20 spectrum and the remaining serial sequence continue; U100 and U1000 await accepted states. Details are in the continuation record.
