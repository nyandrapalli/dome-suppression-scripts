param(
    [string]$Output = (Join-Path $PSScriptRoot "outputs\smoke")
)

$ErrorActionPreference = "Stop"

$repo = (Resolve-Path (Join-Path $PSScriptRoot "..\..\..")).Path
. (Join-Path $repo "scripts\runtime.ps1")
$python = Resolve-PythonExecutable
$runner = Join-Path $PSScriptRoot "run_raw_v7_10_0.ps1"
$packager = Join-Path $PSScriptRoot "package_production_v7_10_0.py"
$validator = Join-Path $PSScriptRoot "validate_production_v7_10_0.py"
$root = [System.IO.Path]::GetFullPath($Output)
$raw = Join-Path $root "raw"
$media = Join-Path $root "media"

if (Test-Path -LiteralPath $root) {
    $existing = Get-ChildItem -LiteralPath $root -Force
    if ($existing.Count -gt 0) {
        throw "Refusing to overwrite nonempty V7.10.0 smoke directory: $root"
    }
}

& $runner `
    -Output $raw `
    -Concentrations "10,80" `
    -Nodes 321 `
    -MaxIterations 12000 `
    -CallbackStride 10 `
    -MaxCallbackFrames 4
if ($LASTEXITCODE -ne 0) {
    throw "V7.10.0 smoke minimization failed."
}

& $python $packager --run $raw --output $media
if ($LASTEXITCODE -ne 0) {
    throw "V7.10.0 smoke packaging failed."
}

& $python $validator `
    --run $raw `
    --media $media `
    --nodes 321 `
    --concentrations "10,80"
if ($LASTEXITCODE -ne 0) {
    throw "V7.10.0 smoke validation failed."
}

Set-Content -LiteralPath (
    Join-Path $root "workflow_status.txt") -Value "smoke_passed"
