# Reproduction and verification

## Check this record without a simulation

From this directory:

```bash
python3 verify_evidence.py
```

This verifies every preserved file against `provenance.json`, recalculates the five ED-comparison z scores, checks the L20/L22 DMRG comparisons and finite-size fit counts, and compares saved final DMRG extrema between the independent audit and returned summary. It writes a new `verification.json` for this assembled record. A pass means file integrity and these arithmetic checks passed; it does not certify every scientific calculation.

## Restore the original calculation source

`source_snapshot/HubbardRenyiDQMC/` is the exact committed tree at `b0adedf042304a72bc1aea618abdcda1962242e9` from the complete 28 August backup. All 183 committed files were restored with `git archive HEAD`; uncommitted U4 thermal-audit work was excluded. The separate git bundle advertises `c103a8b226c0f2d441b45a60c7736357b933c42c`; it is a different saved history and was not substituted for the archived working tree.

The original SHA-256 hashes, archive paths, Git blob IDs and per-file hashes are in `provenance.json`. Git normalizes checkout line endings: the canonical five-subsystem CSV is 460 bytes, while the archived worktree copy is 461 bytes with one CRLF. Both are preserved with distinct hashes and identical numerical rows. The June Downloads directory is an older empty-module scaffold and should not replace this source.

## Rényi package regression suite

The source snapshot targets Julia 1.12.1, pins its package environment, and contains the project's own simulation-free regression suite. From the snapshot directory, in PowerShell:

```powershell
julia --project=. --startup-file=no --history-file=no .\test\runtests.jl
```

On Linux/ISAAC the equivalent is:

```bash
julia --project=. --startup-file=no --history-file=no test/runtests.jl
```

Use the already installed environment when available. If dependencies are missing, the source's `REPRODUCIBILITY.md` records the explicit environment-instantiation step. A package installation is distinct from test execution.

## Fresh Monte Carlo reproduction

The complete copy-paste workflow is the source snapshot's [REPRODUCIBILITY.md](source_snapshot/HubbardRenyiDQMC/REPRODUCIBILITY.md). The example below starts the L10/Ninc8 four-root grid, which is a separate lower-statistics checkpoint from the high-statistics Ninc16 table quoted in this record. It uses a new output directory and preserves accepted data.

```powershell
$raw = Join-Path $env:TEMP ('HubbardRenyiDQMC_' + [guid]::NewGuid().ToString('N'))
New-Item -ItemType Directory -Path $raw | Out-Null
julia --project=. --startup-file=no --history-file=no .\run\run_grid_local.jl .\parameters\L10_U2_beta40_Ninc8.toml --all-seeds --max-jobs 4 --output-dir $raw
julia --project=. --startup-file=no --history-file=no .\analysis\summarize_grid.jl $raw --parameter-file .\parameters\L10_U2_beta40_Ninc8.toml --output (Join-Path $raw 'summary.csv')
```

For the quoted high-statistics half-chain result use `parameters/L10_U2_beta40_Ninc16_highstat.toml` in both commands. For nA1…4 use the correspondingly named parameter files. For new stochastic runs, agreement should be assessed with predeclared tolerances and uncertainty; exact reproduction of a stored decimal is not the statistical acceptance criterion. Give a measured runtime estimate and account for the ISAAC queue before starting a production grid.

## ED and DMRG source records

The ALM evidence files retain the exact returned DMRG runner/environment, 20 stages, four final arms, all returned spin matrices and local observables, and the independent CSV audit. Full reproduction additionally needs the original ED serialized reference `results/literature_transfer/alm_exact_locking_20260911/dmrg_reference_target.jls` in `/home/yiqiwan/code/AltermagnetDQMCTrial`; its SHA-256 is recorded in the metadata and ledger. It is about 109 MB and was not copied into this bounded record.

The larger eight-cell selected-singlet record is a small ISAAC receipt with pinned source-manifest hashes and reported fresh original-H and S² actions. Repeating that solve requires its registered source release, native Julia version and large remote checkpoint. It cannot be recreated from the receipt alone. These limits are explicit so the record remains useful without overstating portability.
