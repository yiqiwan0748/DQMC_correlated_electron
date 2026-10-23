# Reproduction

This is the one supported workflow. It uses the tracked Julia environment and
writes new simulation data to a unique temporary directory. Accepted raw data
and processed benchmark snapshots are never overwritten.

## 1. Set up and test

Use Julia 1.12.1 from the repository root in PowerShell:

```powershell
julia --project=. --startup-file=no --history-file=no -e "using Pkg; Pkg.instantiate()"
julia --project=. --startup-file=no --history-file=no .\test\runtests.jl
```

`Project.toml` and `Manifest.toml` pin the DQMC, DMRG, and plotting
dependencies. No external Julia project is needed.

## 2. Run one complete check

The following command runs the L=10, Ninc=8, four-seed grid with at most four
Julia workers:

```powershell
$raw = Join-Path $env:TEMP ('HubbardRenyiDQMC_' + [guid]::NewGuid().ToString('N'))
New-Item -ItemType Directory -Path $raw | Out-Null

julia --project=. --startup-file=no --history-file=no .\run\run_grid_local.jl .\parameters\L10_U2_beta40_Ninc8.toml --all-seeds --max-jobs 4 --output-dir $raw
```

Rerunning the same command skips every valid complete node. If a node has a
partial or invalid output pair, the launcher stops for inspection instead of
guessing or deleting data. Prefer a new temporary directory over `--force`.

The Julia processes run on the machine where this command is executed. The
launcher does not detach them: over SSH, use `tmux` or `screen` before
disconnecting, and keep the remote machine awake.

## 3. Summarize

Write the processed summary beside the temporary raw files:

```powershell
$summary = Join-Path $raw 'summary.csv'
julia --project=. --startup-file=no --history-file=no .\analysis\summarize_grid.jl $raw --parameter-file .\parameters\L10_U2_beta40_Ninc8.toml --output $summary
Get-Content $summary
```

The summarizer rejects missing nodes or bins, duplicates, wrong metadata,
nonfinite fields, and invalid endpoint fields before calculating the result.

## Optional exact APBC Mott-crossover check

This fast ED-only check reproduces the U=4 values used in the finite-size
crossover diagnostic without starting Monte Carlo production:

```powershell
$diag = Join-Path $env:TEMP ('HubbardRenyiDQMC_U4_ED_' + [guid]::NewGuid().ToString('N'))
New-Item -ItemType Directory -Path $diag | Out-Null

julia --project=. --startup-file=no --history-file=no .\analysis\compute_s2_ed.jl --L 8 --U 4 --nA 4 --boundary apbc --output (Join-Path $diag 'L8.csv')
julia --project=. --startup-file=no --history-file=no .\analysis\compute_s2_ed.jl --L 12 --U 4 --nA 6 --boundary apbc --output (Join-Path $diag 'L12.csv')

Get-Content (Join-Path $diag 'L8.csv')
Get-Content (Join-Path $diag 'L12.csv')

julia --project=. --startup-file=no --history-file=no .\analysis\plot_u2_u4_apbc_mott_crossover.jl --output (Join-Path $diag 'u2_u4_apbc_mott_crossover.svg')
```

The full APBC screening commands and their limitations are recorded in
[`docs/apbc_shell_pilot.md`](docs/apbc_shell_pilot.md) and
[`docs/u4_mott_crossover.md`](docs/u4_mott_crossover.md).

## Detached U=4 finite-size campaign

The recoverable PBC campaign through L=22 has its fixed matrix, Windows task
launcher, status command, and completion gates in
[`docs/u4_fss_L4_to_L22_campaign.md`](docs/u4_fss_L4_to_L22_campaign.md).
It writes only to a persistent external campaign directory; it does not modify
the repository's accepted raw data or processed snapshots. One scheduled-task
launch owns the run, verifies a live heartbeat, and automatically resumes
validated work after interruption.

## Optional kernel benchmark

This checks the optimized rank-one update against the direct implementation
and reports timings without writing accepted results:

```powershell
julia --project=. --startup-file=no --history-file=no .\analysis\validate_grover_rank1_update.jl --benchmark
```

The scientific result, validation numbers, finite-size conclusion, and artifact
provenance are in [`docs/dqmc_renyi_ee.md`](docs/dqmc_renyi_ee.md).
