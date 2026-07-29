param(
    [string]$Output = (Join-Path $PSScriptRoot "outputs\audit")
)

$ErrorActionPreference = "Stop"

$runner = Join-Path $PSScriptRoot "run_raw_v7_10_0.ps1"
$output = [System.IO.Path]::GetFullPath($Output)

& $runner `
    -Output $output `
    -Mode audit `
    -Concentrations "10,20,30,40,50,60,70,80" `
    -Nodes 641 `
    -MaxIterations 16000 `
    -CallbackStride 10 `
    -MaxCallbackFrames 12
if ($LASTEXITCODE -ne 0) {
    throw "V7.10.0 audit-only initialization failed."
}

if ((Get-Content -LiteralPath (
        Join-Path $output "status.txt") -Raw).Trim() -ne
        "audit_only: no minimization was run.") {
    throw "V7.10.0 audit unexpectedly entered minimization mode."
}

Set-Content -LiteralPath (
    Join-Path $output "audit_gate_status.txt") -Value "passed_no_minimization"
