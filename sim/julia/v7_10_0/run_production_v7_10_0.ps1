param(
    [string]$SmokeOutput = (Join-Path $PSScriptRoot "outputs\smoke"),
    [string]$Output = (Join-Path $PSScriptRoot "outputs\production")
)

$ErrorActionPreference = "Stop"

$repo = (Resolve-Path (Join-Path $PSScriptRoot "..\..\..")).Path
. (Join-Path $repo "scripts\runtime.ps1")
$python = Resolve-PythonExecutable
$runner = Join-Path $PSScriptRoot "run_raw_v7_10_0.ps1"
$packager = Join-Path $PSScriptRoot "package_production_v7_10_0.py"
$validator = Join-Path $PSScriptRoot "validate_production_v7_10_0.py"
$smoke = [System.IO.Path]::GetFullPath($SmokeOutput)
$root = [System.IO.Path]::GetFullPath($Output)
$raw = Join-Path $root "raw"
$media = Join-Path $root "media"
$concentrations = "10,20,30,40,50,60,70,80"

$smokeStatus = Join-Path $smoke "workflow_status.txt"
if (-not (Test-Path -LiteralPath $smokeStatus) -or
        (Get-Content -LiteralPath $smokeStatus -Raw).Trim() -ne "smoke_passed") {
    throw "The V7.10.0 smoke workflow has not passed."
}
if (Test-Path -LiteralPath $root) {
    $existing = Get-ChildItem -LiteralPath $root -Force
    if ($existing.Count -gt 0) {
        throw "Refusing to overwrite nonempty V7.10.0 production directory: $root"
    }
}

& $runner `
    -Output $raw `
    -Concentrations $concentrations `
    -Nodes 641 `
    -MaxIterations 16000 `
    -CallbackStride 10 `
    -MaxCallbackFrames 12
if ($LASTEXITCODE -ne 0) {
    throw "V7.10.0 production minimization failed."
}

& $python $packager --run $raw --output $media
if ($LASTEXITCODE -ne 0) {
    throw "V7.10.0 production packaging failed."
}

& $python $validator `
    --run $raw `
    --media $media `
    --nodes 641 `
    --concentrations $concentrations `
    --require-restoration
if ($LASTEXITCODE -ne 0) {
    throw "V7.10.0 production validation failed."
}

Set-Content -LiteralPath (
    Join-Path $root "workflow_status.txt") -Value "production_passed"
