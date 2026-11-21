param([Parameter(Mandatory=$true)][string]$WorkDir)
$ErrorActionPreference = 'Stop'
if ($env:OS -ne 'Windows_NT') { throw 'Native Windows is required; WSL/macOS/Linux do not count.' }
if (Test-Path -LiteralPath $WorkDir) { throw 'Use a new evidence directory; existing output will not be changed.' }
$Repo = (Resolve-Path (Join-Path $PSScriptRoot '../..')).Path
$WorkDir = [IO.Path]::GetFullPath($WorkDir)
if ($Repo -match '[^\x20-\x7e]' -or $Repo.Length -gt 150) { throw 'Use a short ASCII source path, such as C:\saf-source.' }
if ($WorkDir -match '[^\x20-\x7e]' -or $WorkDir.Length -gt 120) { throw 'Use a short ASCII evidence path, such as C:\saf-check.' }
New-Item -ItemType Directory -Path $WorkDir | Out-Null
$Evidence = (Resolve-Path -LiteralPath $WorkDir).Path
$OldTemp = $env:TEMP; $OldTmp = $env:TMP; $OldTmpdir = $env:TMPDIR
$Scratch = Join-Path $Evidence 'temp'
New-Item -ItemType Directory -Path $Scratch | Out-Null
$env:TEMP = $Scratch; $env:TMP = $Scratch; $env:TMPDIR = $Scratch
$Checks = @()
$Failure = $null
$Pushed = $false
$Result = [ordered]@{
    status = 'RUNNING'; native_windows = $true
    os = [System.Environment]::OSVersion.VersionString
    architecture = $env:PROCESSOR_ARCHITECTURE
    scope = 'bounded native Perl platform, kernel, spool, whole/direct matrix, CRLF and audit checks'
    exclusions = @('full panel', 'network paths', 'non-ASCII filesystem paths', 'performance')
}
function Invoke-Logged([string]$Name, [string]$Executable, [string[]]$Arguments) {
    # Windows PowerShell 5.1 can turn a native stderr record into a terminating
    # error under Stop. Capture the command's actual exit status instead.
    $PreviousPreference = $ErrorActionPreference
    # -CommandType Application may return multiple PATH matches. Invoke the
    # first (PATH-priority) executable, not a string containing all matches.
    $ResolvedExecutable = (Get-Command $Executable -CommandType Application -ErrorAction Stop | Select-Object -First 1).Source
    try {
        $ErrorActionPreference = 'Continue'
        $global:LASTEXITCODE = $null
        & $ResolvedExecutable @Arguments 2>&1 | Tee-Object -FilePath (Join-Path $Evidence ($Name + '.log')) -ErrorAction Stop
        $Exit = $global:LASTEXITCODE
    } finally { $ErrorActionPreference = $PreviousPreference }
    $Log = Get-Content -LiteralPath (Join-Path $Evidence ($Name + '.log')) -Raw
    $SkipCount = ([regex]::Matches($Log, '(?im)^ok[^\r\n]*#\s*skip\b')).Count
    $script:Checks += @{case = $Name; exit = $Exit; TAP_skip_count = $SkipCount}
    if ($null -eq $Exit -or $Exit -ne 0) { throw "$Name failed with exit code $Exit." }
}
try {
    Push-Location $Repo
    $Pushed = $true
    $Native = & perl -e 'print $^O'
    if ($LASTEXITCODE -ne 0 -or $Native -ne 'MSWin32') { throw 'Use native Windows Perl, not Cygwin/WSL Perl.' }
    Invoke-Logged 'perl-version' 'perl' @('-V')
    Invoke-Logged 'platform-tests' 'perl' @('tests/windows_platform/run.pl')
    Invoke-Logged 'kernel-tests' 'perl' @('tests/compact_kernel/run.pl')
    Invoke-Logged 'spool-tests' 'perl' @('tests/materialized_spool/run.pl')
    Invoke-Logged 'spool-utf8-tests' 'perl' @('tests/materialized_spool/utf8_bytes.pl')
    Invoke-Logged 'pipeline-smoke' 'python' @('-I', 'tests/windows_platform/smoke.py', '--require-windows', '--work-dir', (Join-Path $Evidence 'smoke space'))
    $Smoke = Get-Content -LiteralPath (Join-Path $Evidence 'smoke space/RESULTS.json') -Raw | ConvertFrom-Json
    if ($Smoke.status -ne 'PASS' -or -not $Smoke.native_windows -or $Smoke.cases.Count -ne 18) { throw 'Native smoke receipt is missing or incomplete.' }
    $Residue = @(Get-ChildItem -LiteralPath $Scratch -Force)
    if ($Residue.Count -ne 0) { throw 'Private temporary files were not fully cleaned; see the temp evidence folder.' }
    $Result.status = 'PASS'
} catch {
    $Failure = $_
    $Result.status = 'FAIL'
    $Result.error = $_.Exception.Message
} finally {
    $env:TEMP = $OldTemp; $env:TMP = $OldTmp; $env:TMPDIR = $OldTmpdir
    if ($Pushed) { Pop-Location }
    $Result.checks = $Checks
    $Result | ConvertTo-Json -Depth 6 | Out-File -Encoding utf8 (Join-Path $Evidence 'RESULTS.json')
    ($Result.status + ': ' + $Result.scope) | Out-File -Encoding utf8 (Join-Path $Evidence 'RESULT.txt')
}
if ($null -ne $Failure) { throw $Failure }
