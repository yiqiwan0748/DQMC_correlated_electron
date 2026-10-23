# AGENTS.md

Guidance for agents working in this repository.

## Project Context

This is a Julia DQMC research project for Hubbard Renyi entropy benchmarks.
The active workflow studies lambda-grid estimators for the 1D Hubbard model,
with raw per-node outputs summarized into processed benchmark tables.

Important active paths:

- `run/run_lambda_node.jl`: active single-lambda-node runner.
- `run/run_grid_local.jl`: local parallel lambda-grid launcher.
- `analysis/summarize_grid.jl`: summarizes raw bin CSVs into benchmark tables.
- `parameters/`: benchmark TOML files.
- `results/raw/`: raw per-node simulation outputs.
- `results/processed/`: processed benchmark summaries, logs, and snapshots.
- `archive/validated_old_pipeline/`: validated reference implementation.

## Hard Constraints

- Do not modify `archive/validated_old_pipeline/` unless explicitly asked.
- Do not modify `results/raw/`.
- Do not overwrite processed benchmark snapshots unless explicitly asked.
- Do not change estimator logic without explaining the mathematical identity being preserved.
- Before running Julia simulations, ask first.
- It is okay to run lightweight read-only PowerShell inspection commands.
- Preserve PowerShell-compatible commands in docs, scripts, and examples.
- Prefer small, reviewable patches.

## Physics Conventions

- `ratio` means `Z(lambda_R) / Z(lambda_L)`.
- `dS2 = -log(ratio)`.
- The product of adjacent ratios gives `Z1/Z0`.
- `S2 = -log(Z1/Z0)`.
- ED exact benchmark value: `S2 = 1.716551661433`.

## Current Workflow

Typical PowerShell-compatible commands:

```powershell
julia .\run\run_lambda_node.jl .\parameters\L10_U2_beta40_Ninc8.toml --lambda-index 4 --root-seed 1234
julia .\run\run_grid_local.jl .\parameters\L10_U2_beta40_Ninc8.toml --root-seed 1234 --max-jobs 4
julia .\analysis\summarize_grid.jl .\results\raw\L10_U2_beta40_Ninc8 --parameter-file .\parameters\L10_U2_beta40_Ninc8.toml
```

Do not run the simulation commands without first asking the user.
Summary commands may still be expensive if they process many files; use judgment and
ask first if there is any risk of rewriting benchmark outputs.

## Editing Guidance

- Keep patches narrowly scoped to the requested files.
- Treat `archive/validated_old_pipeline/` as a reference, not the active editing target.
- Treat `results/raw/` as immutable experimental data.
- Treat named processed CSVs such as `FRIDAY_BENCHMARK_*`, `BETA_CONVERGENCE_*`,
  and `lambda_grid_*_highstat.csv` as benchmark snapshots.
- If a processed output must be regenerated, write to a new explicit path unless the
  user explicitly asks to overwrite an existing snapshot.
- Maintain the one user-facing `root_seed` convention; replica RNG streams are derived
  internally by the runner.
