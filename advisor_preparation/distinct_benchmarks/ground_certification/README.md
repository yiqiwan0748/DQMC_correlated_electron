# X09 / X12: external ground certification and degenerate ensembles

Prepared protocol; **not executed**. This release can establish at most two distinct catalogue ideas. No state, seed, sweep, cut or temperature is counted as an extra idea. X08 and X10 remain outside this tranche because their catalogue definitions require two-orbital coverage.

The targets import the existing, hash-pinned user ED, thermal-ED and Hubbard-MPO source. They add explicit state initialization, matrix/amplitude export, projector certification and prescribed-ensemble adapters. Native DMRG is provided by ITensorMPS. No production kernel or historical result is changed. QuSpin independently constructs the fermion operators and Hamiltonians; analytic dimer formulas supply another reference for X12. This is not a native DQMC benchmark.

## X09: a ground state needs more than small variance

The fixture is the interacting four-site periodic Hubbard ring, hopping1, U4, physical chemical potential0, in exactly Nup=Ndown2. QuSpin diagonalizes the complete36-dimensional sector and supplies its unique ground-state projector and gap. The actual user `build_space` / `apply_hamiltonian!` Hamiltonian and the actual user `hubbard_mpo`, with both QN and dense physical indices, must match that whole matrix after the declared fermionic basis transformation.

The native optimizer runs all eight declared sweeps from four seeded starts: Néel, doublon/hole charge order, spin domains and shuffled spins. Every resulting state must have the right norm and sector, variance in [−1e−9,1e−8], absolute energy error at most1e−8 and fidelity at least0.99999999 to the independent ground projector. Source eigensystem residuals must be at most1e−10. All settings, states and seeds are fixed in `manifest.toml` before execution.

Two additional states are manufactured from the first excited and highest eigenvectors of the user ED Hamiltonian and converted to dense MPS without truncation. Both must have small variance and be rejected by the external ground certificate. They test the certification logic; they do **not** demonstrate that the native optimizer naturally became trapped. Successful stress cases alone cannot pass X09: all four actual DMRG starts must also pass ground certification.

For any normalized exact eigenstate, Var(H)=0 independently of its energy rank. That elementary spectral identity is the reason for the projector and energy requirements. Variance remains useful as a residual diagnostic; its use in MPS error estimates is discussed by [Hubig, Haegeman and Schollwöck](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.97.045125).

The user ED ordering lists all up modes before down modes. ITensor Electron indices order the local states Emp, Up, Dn, UpDn, with the first tensor index varying fastest in the dense array. The adapter applies the parity of moving each occupied down mode past occupied up modes to its right. Both full Hamiltonian comparisons and the exported raw state amplitudes make this convention auditable. The MPO imaginary maximum is retained and must also be below1e−10; taking a real part cannot conceal a discrepancy.

## X12: match the ensemble within a degenerate manifold

The fixture is the interacting dimer at U4, effective hopping1 and canonical **total N3**, including both spin sectors. The source L2 periodic bond loop counts the same physical edge twice, so its input hopping is0.5. Uncentered physical energies are3,3,5,5. The full-Fock thermal source uses the centered interaction; the adapter restores Hphysical=Hcentered+(U/2)N−UL/4 before constructing the canonical thermal density. At N3,L2 the restoration adds4. The fixed-sector ED solver already uses the uncentered Hamiltonian.

Each backend separately constructs and checks these prescriptions:

- The actual fixed-sector Lanczos ground state with Nup2,Ndown1.
- Its coherent spin rotation cos(theta)|up>+exp(i phi)sin(theta)|down>, with theta=pi/7, phi=pi/5 and |down> fixed by normalized total Sminus acting on |up>.
- The equal mixture of the complete two-dimensional ground manifold.
- The same mixture after a fixed unitary rotation of an orthonormal ground basis.
- The canonical Gibbs density at each beta0.5,1,2,4,8,16, formed from all four canonical states.

The target calls the existing `ground_state_lanczos`, `exact_thermal_density` and `renyi2` / `reduced_operator` functions. Both one-site and full-system RDM spectra, purity, S2, particle number, energy, all local Sx/Sy/Sz means, Czz, connected Czz and trace norm to the equal-ground mixture are compared with QuSpin and the independent formulas below. Local spin is S=(1/2)c-dagger sigma c. The transverse spin means make the complex rotation phase measurable. Complete RDM eigenvalue spectra, trace and positivity are required; no arbitrary eigenvector gauge is compared.

Writing epsilon=1/(1+exp(2 beta)), the thermal excited-manifold weight is epsilon, E=3+2epsilon and ||rho_beta−Pground/2||_1=2epsilon. The norm here is the full trace norm, without a factor1/2. A pure ground vector has one-site spectrum(1/2,1/2,0,0), purity1/2 and S2=log2. The equal mixture and this fixture's canonical Gibbs states have one-site spectrum(1/2,1/4,1/4,0), purity3/8 and S2=log(8/3). Their full spectra distinguish finite-temperature excited-state admixture even though their one-site spectra agree.

For the spin-up component, both sites have <Sz>=1/4 and connected <Sz1 Sz2>=−1/16; for the coherent rotation the local spin vector is(sin(2theta)cos(phi),sin(2theta)sin(phi),cos(2theta))/4 and connected Czz=−cos(2theta)^2/16. Mixtures and Gibbs states have zero local spin and connected Czz. Every prescription has unconnected Czz0. These expressions follow by diagonalizing the two-state hopping block in each N3 spin sector and tracing its one-site probabilities. They are analytic references, not fitted outputs or a claim of reproducing a journal plot.

The last beta must approach the equal-ground mixture within trace norm1e−10. The fixture must also distinguish a pure component from the mixture by at least0.1 in S2 and0.01 in connected Czz. The comparator preserves this failure sensitivity rather than substituting a single ground vector for a mixture.

## Execution and evidence

Root owns staging, resource checks and ISAAC submission. `run.sbatch` specifies one CPU,3500MB and a15-minute unmeasured cap, one BLAS/Julia thread, the existing isolated QuSpin environment and pinned Julia project. No nested jobs are submitted. Preserve every failed or incomplete attempt; use a new result path for each run. An execution, API, version, source-pin, convergence, RDM-positivity or missing-grid problem is UNRESOLVED. Eligible scientific disagreement is FAIL. A pilot, if requested, cannot pass a catalogue idea.

Run with the release copied to an absolute `CERTIFICATION_TRANCHE` path and campaign root in `VALIDATION_ROOT`. The batch script writes `oracle.json`, `target.toml`, resource logs and `comparison.json` under a new result directory. Postprocessing can be replayed without physics:

```bash
export RENYI_ROOT=/absolute/pinned/user/source
python3 compare.py --manifest manifest.toml --oracle /absolute/result/oracle.json \
  --target /absolute/result/target.toml --output /absolute/new/comparison_replay.json
```

The comparator uses Python's standard library. On ISAAC, the existing `quspin-venv/bin/python` can run it. It checks the immutable raw-file, manifest, adapter, user-source and runtime identities and all required arms. Local syntax checks do not execute a model or count as benchmark completion. Publication must preserve the source closure, full raw receipts, resource logs and exact comparator version, then independently replay the comparator before updating the100-row ledger.

Primary implementation references: [QuSpin Part II](https://arxiv.org/abs/1804.06782), [QuSpin fermion-basis documentation](https://quspin.github.io/QuSpin/basis.html), [native DMRG API](https://itensor.github.io/ITensorMPS.jl/stable/DMRG.html), and the pinned [ITensorMPSv0.4.1 dense-tensor MPS constructor](https://raw.githubusercontent.com/ITensor/ITensorMPS.jl/v0.4.1/src/abstractmps.jl). QuSpin supplies the independent package oracle; the target optimizer and its MPO belong to the existing user workflow and are not counted as an external DMRG implementation.
