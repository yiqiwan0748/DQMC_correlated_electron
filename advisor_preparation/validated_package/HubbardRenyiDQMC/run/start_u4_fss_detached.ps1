param(
    [string]$RunDirectory = 'C:\Users\YQW\HubbardRenyiDQMC_runs\u4_fss_L4_L22_20260822',
    [string]$TaskName = 'HubbardRenyiDQMC-U4-FSS',
    [ValidateRange(1, 16)]
    [int]$MaxJobs = 4
)

$ErrorActionPreference = 'Stop'

$allowedRoot = [IO.Path]::GetFullPath('C:\Users\YQW\HubbardRenyiDQMC_runs')
$fullRunDirectory = [IO.Path]::GetFullPath($RunDirectory)
$rootPrefix = $allowedRoot.TrimEnd('\') + '\'
if (-not $fullRunDirectory.StartsWith($rootPrefix, [StringComparison]::OrdinalIgnoreCase)) {
    throw "RunDirectory must be below $allowedRoot"
}

if ($fullRunDirectory -notmatch '^([A-Za-z]):\\(.*)$') {
    throw 'RunDirectory must be a local Windows drive path'
}
$drive = $Matches[1].ToLowerInvariant()
$tail = $Matches[2].Replace('\', '/')
$wslRunDirectory = "/mnt/$drive/$tail"

New-Item -ItemType Directory -Path $fullRunDirectory -Force | Out-Null

$wsl = Join-Path $env:SystemRoot 'System32\wsl.exe'
$arguments = @(
    '-d Ubuntu'
    '-u yqw'
    '--cd /home/yqw/code/HubbardRenyiDQMC'
    '-- /bin/bash run/run_u4_fss_campaign.sh'
    ('"{0}"' -f $wslRunDirectory)
    $MaxJobs
    ('"{0}"' -f $TaskName)
) -join ' '

$action = New-ScheduledTaskAction -Execute $wsl -Argument $arguments
$trigger = New-ScheduledTaskTrigger `
    -Once `
    -At (Get-Date).AddMinutes(1) `
    -RepetitionInterval (New-TimeSpan -Minutes 5) `
    -RepetitionDuration (New-TimeSpan -Days 14)
$settings = New-ScheduledTaskSettingsSet `
    -MultipleInstances IgnoreNew `
    -StartWhenAvailable `
    -AllowStartIfOnBatteries `
    -DontStopIfGoingOnBatteries `
    -WakeToRun `
    -ExecutionTimeLimit ([TimeSpan]::Zero) `
    -Priority 7
$principal = New-ScheduledTaskPrincipal `
    -UserId "$env:USERDOMAIN\$env:USERNAME" `
    -LogonType Interactive `
    -RunLevel Limited

$existing = Get-ScheduledTask -TaskName $TaskName -ErrorAction SilentlyContinue
if ($null -ne $existing) {
    $existingAction = $existing.Actions | Select-Object -First 1
    if ($existingAction.Execute -ne $wsl -or $existingAction.Arguments -ne $arguments) {
        throw "Task $TaskName already exists with a different action; refusing to replace it"
    }
    Set-ScheduledTask `
        -TaskName $TaskName `
        -Action $action `
        -Trigger $trigger `
        -Settings $settings `
        -Principal $principal | Out-Null
} else {
    $task = New-ScheduledTask `
        -Action $action `
        -Trigger $trigger `
        -Settings $settings `
        -Principal $principal `
        -Description 'Recoverable detached U=4 Hubbard Renyi FSS campaign through L=22'
    Register-ScheduledTask -TaskName $TaskName -InputObject $task | Out-Null
}

Enable-ScheduledTask -TaskName $TaskName | Out-Null

$heartbeatPath = Join-Path $fullRunDirectory 'state\heartbeat.txt'
$launchStartedUtc = [DateTime]::UtcNow
Start-ScheduledTask -TaskName $TaskName

$deadlineUtc = $launchStartedUtc.AddSeconds(45)
$launchVerified = $false
$heartbeat = $null
do {
    $taskState = Get-ScheduledTask -TaskName $TaskName
    $heartbeat = Get-Item -LiteralPath $heartbeatPath -ErrorAction SilentlyContinue
    if (
        $taskState.State -eq 'Running' -and
        $null -ne $heartbeat -and
        $heartbeat.LastWriteTimeUtc -ge $launchStartedUtc
    ) {
        $launchVerified = $true
        break
    }
    Start-Sleep -Seconds 2
} while ([DateTime]::UtcNow -lt $deadlineUtc)

$taskInfo = Get-ScheduledTaskInfo -TaskName $TaskName
if (-not $launchVerified) {
    $heartbeatStatus = if ($null -eq $heartbeat) {
        'missing'
    } else {
        "last written $($heartbeat.LastWriteTimeUtc.ToString('o')) UTC"
    }
    throw "Task $TaskName did not prove a detached launch within 45 seconds: state=$($taskState.State), heartbeat=$heartbeatStatus, last result=$($taskInfo.LastTaskResult)"
}

[pscustomobject]@{
    TaskName = $TaskName
    State = $taskState.State
    LastTaskResult = $taskInfo.LastTaskResult
    Heartbeat = $heartbeat.FullName
    HeartbeatUtc = $heartbeat.LastWriteTimeUtc.ToString('o')
    RunDirectory = $fullRunDirectory
    WslRunDirectory = $wslRunDirectory
    Repetition = '5 minutes for 14 days'
    Recovery = 'completed nodes and sizes are validated and skipped'
}
