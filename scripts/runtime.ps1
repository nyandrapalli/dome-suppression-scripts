function Resolve-JuliaExecutable {
    $candidates = @()
    if ($env:JULIA_EXE) {
        $candidates += $env:JULIA_EXE
    }
    $candidates += Get-Command julia -All -ErrorAction SilentlyContinue |
        Select-Object -ExpandProperty Source

    foreach ($candidate in $candidates | Select-Object -Unique) {
        if (-not (Test-Path -LiteralPath $candidate)) {
            continue
        }
        & $candidate -e "exit(VERSION.major > 1 || (VERSION.major == 1 && VERSION.minor >= 12) ? 0 : 1)" *> $null
        if ($LASTEXITCODE -eq 0) {
            return $candidate
        }
    }
    throw "Julia 1.12 or newer was not found. Set JULIA_EXE to julia.exe."
}

function Resolve-PythonExecutable {
    $candidates = @()
    if ($env:PYTHON_EXE) {
        $candidates += $env:PYTHON_EXE
    }
    $candidates += Get-Command python -All -ErrorAction SilentlyContinue |
        Select-Object -ExpandProperty Source
    $candidates += Get-Command python3 -All -ErrorAction SilentlyContinue |
        Select-Object -ExpandProperty Source

    foreach ($candidate in $candidates | Select-Object -Unique) {
        if (-not (Test-Path -LiteralPath $candidate)) {
            continue
        }
        & $candidate -c "import sys; raise SystemExit(0 if sys.version_info >= (3, 10) else 1)" *> $null
        if ($LASTEXITCODE -eq 0) {
            return $candidate
        }
    }
    throw "Python 3.10 or newer was not found. Set PYTHON_EXE to python.exe."
}
