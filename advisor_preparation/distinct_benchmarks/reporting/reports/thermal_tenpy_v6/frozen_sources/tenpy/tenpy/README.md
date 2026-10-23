# Independent TeNPy/ITensor finite-chain tranche

Execution status: **NOT RUN**. Only syntax, source interfaces, and the mathematical conventions have been reviewed. The root agent stages and launches these scripts; none submit jobs themselves.

The frozen manifest defines three different scientific comparisons. Convergence arms and eigenpair checks do not count as additional ideas.

| Catalogue ID | Fixture and question | Observable |
|---|---|---|
| P11 | Open 8-site Hubbard chain, 6 electrons, zero spin projection: are Schmidt levels assigned to the correct physical charge sectors? | Prefix of 4 sites; probabilities and natural-log entanglement levels grouped by left electron number and twice spin projection |
| P12 | Same doped chain: does boundary-induced charge modulation agree between independent MPS implementations? | All 8 local densities and their complex Fourier amplitudes on the fixed grid q=pi*m/9, m=1..8 |
| P15 | Open 8-site chain, 8 electrons, zero spin projection: are nonlocal charge-fluctuation strings correct? | All 36 contiguous intervals, exp[i*pi*sum(n-1)] |

Both models use hopping 1, repulsion 4, no boundary potentials, and the uncentered interaction U*n_up*n_down. The chemical-potential coefficient is zero. TeNPy's `FermiHubbardChain` and the user's existing `hubbard_mpo(..., boundary="obc")` represent the same Hamiltonian. Finite MPS boundary conditions and physical open boundaries are both explicit. The new target adapter supplies a doped product state; it does not change the user's Hamiltonian or DMRG kernel. ITensorMPS supplies the optimizer, while the user's implementation supplies the Hubbard MPO and established workflow.

P11 obtains physical sectors independently of virtual charge gauges. TeNPy contracts its prefix RDM and groups basis states using the diagonal local Ntot and 2*Sz operators in the actual sorted onsite basis. It also checks `entanglement_spectrum(by_charge=True)`. Tensor conservation gives qtotal_i=qphysical_i+Q_i-Q_(i+1), so the native cut charge becomes physical prefix charge by adding sum(prefix qtotal_i)-Q_left. This fixed offset and raw charge names are recorded; it is never fitted to measured spectra. The target instead contracts its left-canonical prefix amplitude and computes singular values within physical occupation sectors. Julia's first physical index runs fastest, so its row decoding uses little-endian base-4 digits. No disjoint-region fermion permutation is involved.

P15 uses the identity exp[i*pi*sum(n-1)]=(-1)^length*product((-1)^n). TeNPy contracts local `JW` operators; ITensor builds the product MPO from local `F` operators. These are even diagonal operators and require no added fermion swap strings.

The independent sparse-ED gate uses a separately built interleaved spin-orbital Fock Hamiltonian and two lowest SciPy eigenpairs. The largest fixed sector has 4,900 states. This supplies a numerically resolved same-sector gap and eigenpair residuals; it is not a mathematical ground-state certificate. Both MPS energies must match that ground energy, remain separated from the next level, and pass norm, particle-number, spin, and energy-variance gates. Coarse/refined bond caps are 128/256, with different seeds. Probability, density, and string tolerances are fixed in the manifest. Tiny Schmidt levels are compared as probabilities; log-level errors are assessed only above the declared probability floor. Nonconvergence, execution failure, an unresolved native charge convention, or missing refinement gives UNRESOLVED. An observable discrepancy after eligibility gates gives FAIL.

The bounded pilot runs only doped/coarse. It measures runtime and peak resident memory, and cannot complete any catalogue idea. Setup uses a new `$VALIDATION_ROOT/tenpy-venv`, never the QuSpin environment. TeNPy, NumPy, and SciPy are pinned; `pip_freeze.txt` records all resolved transitive versions. The pilot has separate 240-second driver limits and a 15-minute allocation including installation. The full run has a 20-minute allocation at one CPU and 3500 MB. These are unmeasured hard caps, not performance predictions; inspect the pilot before launching full.

Root submission after staging:

```bash
export VALIDATION_ROOT=/lustre/isaac24/scratch/ywang287/advisor-validation-20261022
export TENPY_TRANCHE=/absolute/path/to/deployed/tenpy
sbatch --dependency=afterany:6457686 setup_pilot.sbatch
# Inspect the pilot's JSON/TOML and *.resources.txt before any full run.
sbatch full.sbatch
```

Primary API and version sources:

- [TeNPy 1.1.1 distribution](https://pypi.org/project/physics-tenpy/1.1.1/)
- [FermiHubbardModel conventions](https://tenpy.readthedocs.io/en/v1.1.1/reference/tenpy.models.hubbard.FermiHubbardModel.html)
- [MPS RDM, charge, spectrum, and expectation APIs](https://tenpy.readthedocs.io/en/v1.1.1/reference/tenpy.networks.mps.MPS.html)
- [Tagged MPS source, including the boundary-charge convention](https://raw.githubusercontent.com/tenpy/tenpy/v1.1.1/tenpy/networks/mps.py)
- [Tagged onsite charge/operator source](https://raw.githubusercontent.com/tenpy/tenpy/v1.1.1/tenpy/networks/site.py)
- [MPO variance](https://tenpy.readthedocs.io/en/v1.1.1/reference/tenpy.networks.mpo.MPO.html)

All conclusions are confined to these frozen finite-size fixtures. A 12-site extension or additional fillings would extend coverage, without becoming extra scientific ideas.
