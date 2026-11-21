param(
    [Parameter(Mandatory = $true)][string]$Results,
    [string]$OpenSSLRoot = "",
    [string]$ToolchainFile = "",
    [string]$Python = "python",
    [string]$Generator = "Visual Studio 17 2022"
)
$ErrorActionPreference = "Stop"
if ($env:OS -ne "Windows_NT") { throw "Native Windows is required; WSL is a separate target." }
if (-not [Environment]::Is64BitOperatingSystem) { throw "The first Windows port targets 64-bit systems." }
if (Test-Path -LiteralPath $Results) { throw "Choose a new results directory; existing paths are not reused." }
$source = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
$resultsPath = [IO.Path]::GetFullPath($Results)
if ($resultsPath.StartsWith($source + [IO.Path]::DirectorySeparatorChar, [StringComparison]::OrdinalIgnoreCase)) {
    throw "Put results outside the source tree."
}
New-Item -ItemType Directory -Path $resultsPath | Out-Null
$build = Join-Path $resultsPath "build"
$status = [ordered]@{status="STARTED"; platform="native Windows"; started_utc=[DateTime]::UtcNow.ToString("o");
    os=[Environment]::OSVersion.VersionString; architecture=$env:PROCESSOR_ARCHITECTURE;
    source=$source; generator=$Generator; tests="pending"; windows_validation=$true}
$oldPath = $env:PATH
function Run-Logged([string]$Program, [string[]]$Arguments, [string]$Log) {
    # Windows PowerShell 5.1 wraps native stderr as ErrorRecord objects.
    # Record diagnostics without interpreting a compiler warning as process failure.
    $savedPreference = $ErrorActionPreference
    try {
        $ErrorActionPreference = "Continue"
        $global:LASTEXITCODE = $null
        & $Program @Arguments 2>&1 | Tee-Object -FilePath (Join-Path $resultsPath $Log) -ErrorAction Stop | Out-Host
        $code = $LASTEXITCODE
    } finally { $ErrorActionPreference = $savedPreference }
    if ($null -eq $code -or $code -ne 0) { throw "$Program did not finish successfully (exit $code); see $Log" }
}
try {
    $cmake = (Get-Command cmake -ErrorAction Stop).Source
    $ctest = (Get-Command ctest -ErrorAction Stop).Source
    $pythonCommand = (Get-Command $Python -ErrorAction Stop).Source
    $pythonPath = (& $pythonCommand -I -B -c "import sys, struct; print(sys.executable); sys.exit(0 if sys.version_info >= (3,8) and struct.calcsize('P') == 8 else 1)")
    if ($LASTEXITCODE -ne 0) { throw "Python 3.8 or later is required for the complete test suite." }
    $pythonPath = ($pythonPath | Out-String).Trim()
    Run-Logged -Program $cmake -Arguments @("--version") -Log "cmake-version.log"
    Run-Logged -Program $ctest -Arguments @("--version") -Log "ctest-version.log"
    foreach ($versionLog in @("cmake-version.log","ctest-version.log")) {
        $first = Get-Content -LiteralPath (Join-Path $resultsPath $versionLog) -TotalCount 1
        if ($first -notmatch 'version (\d+\.\d+(?:\.\d+)?)' -or [version]$Matches[1] -lt [version]"3.26") {
            throw "This verification script requires CMake and CTest 3.26 or later."
        }
    }
    Run-Logged -Program $pythonPath -Arguments @("--version") -Log "python-version.log"
    $configure = @("-S",$source,"-B",$build,"-G",$Generator,"-DCMAKE_BUILD_TYPE=Release",
        "-DSAF_BUILD_COMMIT=fast-v0.2.0-rc.3","-DSAF_REQUIRE_PYTHON_TESTS=ON","-DPython3_EXECUTABLE=$pythonPath")
    if ($Generator -like "Visual Studio*") { $configure += @("-A","x64") }
    if ($OpenSSLRoot) {
        $configure += "-DOPENSSL_ROOT_DIR=$OpenSSLRoot"
        $env:PATH = (Join-Path $OpenSSLRoot "bin") + ";" + $oldPath
    }
    if ($ToolchainFile) { $configure += "-DCMAKE_TOOLCHAIN_FILE=$ToolchainFile" }
    $hashes = [ordered]@{}
    foreach ($file in Get-ChildItem -LiteralPath $source -Recurse -File) {
        $relative = $file.FullName.Substring($source.Length + 1)
        if ($relative -match '^(src|include|adapters|tests|examples)[\\/]' -or $relative -eq "CMakeLists.txt") {
            if ($relative -notmatch '[\\/]__pycache__[\\/]') {
                $hashes[$relative] = (Get-FileHash -LiteralPath $file.FullName -Algorithm SHA256).Hash.ToLowerInvariant()
            }
        }
    }
    $hashes | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $resultsPath "SOURCE_SHA256.json") -Encoding UTF8
    $configure | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $resultsPath "CONFIGURE_ARGUMENTS.json") -Encoding UTF8
    Run-Logged -Program $cmake -Arguments $configure -Log "configure.log"
    Run-Logged -Program $cmake -Arguments @("--build",$build,"--config","Release","--parallel","2") -Log "build.log"
    $testList = & $ctest --test-dir $build -C Release --show-only=json-v1
    if ($LASTEXITCODE -ne 0) { throw "Cannot enumerate CTest tests." }
    $testList | Set-Content -LiteralPath (Join-Path $resultsPath "TESTS.json") -Encoding UTF8
    $names = (($testList | Out-String | ConvertFrom-Json).tests | ForEach-Object { $_.name })
    foreach ($required in @("native","collision","c_abi","platform","python_binding","cli_layered_groups",
        "cli_inputs","cli_support_tree","cli_fusion_matrix","cli_quartet","cli_windows_port")) {
        if ($names -notcontains $required) { throw "Required test is missing: $required" }
    }
    Run-Logged -Program $ctest -Arguments @("--test-dir",$build,"-C","Release","--no-tests=error","--verbose","--output-on-failure") -Log "ctest.log"
    $install = Join-Path $resultsPath "install"
    Run-Logged -Program $cmake -Arguments @("--install",$build,"--config","Release","--prefix",$install) -Log "install.log"
    $installedCLI = Join-Path $install "bin/saf.exe"
    Run-Logged -Program $installedCLI -Arguments @("--version") -Log "installed-version.log"
    Run-Logged -Program $pythonPath -Arguments @("-I","-B",(Join-Path $source "tests/test_windows_port.py"),$installedCLI) -Log "installed-cli.log"
    $installedDLL = Join-Path $install "bin/splitaligner_fast.dll"
    if (-not (Test-Path -LiteralPath $installedDLL)) { $installedDLL = Join-Path $install "bin/libsplitaligner_fast.dll" }
    $adapter = Join-Path $install "share/SplitAlignerFast/python"
    Run-Logged -Program $pythonPath -Arguments @("-I","-B","-c",
        "import sys; sys.path.insert(0,sys.argv[1]); from splitaligner_fast import Library; Library(sys.argv[2]); print('installed DLL and Python adapter loaded')",
        $adapter,$installedDLL) -Log "installed-library.log"
    foreach ($relative in $hashes.Keys) {
        $now = (Get-FileHash -LiteralPath (Join-Path $source $relative) -Algorithm SHA256).Hash.ToLowerInvariant()
        if ($now -ne $hashes[$relative]) { throw "Source changed during validation: $relative" }
    }
    $status.tests = "CTest and Python suites completed; inspect ctest.log for individual skips"
    $status.status = "PASS_NATIVE_WINDOWS_SUITE"
} catch {
    $status.status = "FAILED"
    $status.error = $_.Exception.Message
    throw
} finally {
    $env:PATH = $oldPath
    $status.finished_utc = [DateTime]::UtcNow.ToString("o")
    $status | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $resultsPath "RESULTS.json") -Encoding UTF8
}
