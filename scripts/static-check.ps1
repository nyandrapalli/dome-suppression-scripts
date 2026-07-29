$ErrorActionPreference = "Stop"

$root = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
. (Join-Path $PSScriptRoot "runtime.ps1")
$julia = Resolve-JuliaExecutable
$python = Resolve-PythonExecutable

& $julia (Join-Path $PSScriptRoot "static-check.jl")
if ($LASTEXITCODE -ne 0) {
    throw "Julia static checks failed."
}

$pythonFiles = @(
    Join-Path $root "sim\julia\v7_10_0\package_production_v7_10_0.py"
    Join-Path $root "sim\julia\v7_10_0\validate_production_v7_10_0.py"
)
& $python -m py_compile $pythonFiles
if ($LASTEXITCODE -ne 0) {
    throw "Python syntax checks failed."
}

$parseErrors = @()
Get-ChildItem -LiteralPath $root -Recurse -Filter *.ps1 | ForEach-Object {
    $tokens = $null
    $errors = $null
    [System.Management.Automation.Language.Parser]::ParseFile(
        $_.FullName, [ref]$tokens, [ref]$errors) | Out-Null
    $parseErrors += $errors
}
if ($parseErrors.Count -gt 0) {
    $parseErrors | Format-List
    throw "PowerShell syntax checks failed."
}

Write-Host "Repository static checks passed."
