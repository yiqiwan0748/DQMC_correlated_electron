Set-Location "C:\Users\YQW\RemoteCode\HubbardRenyiDQMC"
$ErrorActionPreference = "Stop"

$expectedCount = 68
$nAValues = 1, 2, 3, 4, 5, 6, 7, 8

foreach ($nA in $nAValues) {
    $param = ".\parameters\L16_U2_beta64_nA$($nA)_Ninc16_highstat.toml"
    $rawDir = ".\results\raw\L16_U2_beta64_nA$($nA)_Ninc16_highstat"
    $output = ".\results\processed\L16_U2_beta64_subsystem_nA$($nA)_Ninc16.csv"

    Write-Host "=== nA=$nA: running lambda grid ==="
    julia .\run\run_grid_local.jl $param --all-seeds --max-jobs 4

    $summaryCount = 0
    $binsCount = 0

    if (Test-Path $rawDir) {
        $summaryCount = @(Get-ChildItem -Path $rawDir -Filter "*summary*.csv" -File).Count
        $binsCount = @(Get-ChildItem -Path $rawDir -Filter "*bins*.csv" -File).Count
    }

    Write-Host "nA=$nA summary CSVs: $summaryCount / $expectedCount"
    Write-Host "nA=$nA bins CSVs:    $binsCount / $expectedCount"

    if (($summaryCount -eq $expectedCount) -and ($binsCount -eq $expectedCount)) {
        Write-Host "=== nA=$nA: summarizing complete raw output ==="
        julia .\analysis\summarize_grid.jl $rawDir --parameter-file $param --output $output
    } else {
        Write-Warning "nA=$nA incomplete; skipping summarization."
    }
}

