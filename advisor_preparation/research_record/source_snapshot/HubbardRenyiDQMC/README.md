# Hubbard Renyi DQMC

A compact Julia workflow for the second Renyi entanglement entropy of the
one-dimensional Hubbard model.

## Current result

- The lambda-grid DQMC estimator agrees with ED and DMRG benchmarks.
- The accepted rank-one Grover update is 3.8236x faster on the matched L=10
  production workload.
- Finite-size data are validated through L=22.
- Strong mod-4 corrections prevent a reliable central-charge extraction, so no
  central-charge value is claimed.
- Closed-shell APBC and exact U=4 references support a finite-size Mott
  crossover: the L=8,12 diagnostic moves from `c_eff=1.989` at U=2 to
  `c_eff=1.108` at U=4. The corresponding DQMC runs are screening pilots.

## Start here

- [`docs/dqmc_renyi_ee.md`](docs/dqmc_renyi_ee.md): the authoritative
  scientific result and method.
- [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md): the one safe, copy-paste workflow.

The tracked Julia environment is self-contained; no machine-specific project
is required.

## Validation

Run the fast simulation-free suite:

```powershell
julia --project=. --startup-file=no --history-file=no .\test\runtests.jl
```

## Repository map

- `src/`: active DQMC and Grover implementation.
- `run/`: the node runner and rolling local launcher.
- `analysis/summarize_grid.jl`: the active raw-data validator and summarizer.
- `parameters/`: tracked run configurations.
- `test/`: simulation-free regression tests.
- `docs/dqmc_renyi_ee.md`: authoritative study note.
- `results/processed/` and `results/figures/`: accepted artifacts.
- `archive/validated_old_pipeline/`: immutable historical reference.

Other analysis files and detailed notes are retained as provenance for the
accepted benchmarks; they are not additional user-facing workflows.
