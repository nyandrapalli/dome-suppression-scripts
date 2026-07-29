param(
    [Parameter(Mandatory = $true)]
    [string]$Output,

    [ValidateSet("audit", "smoke")]
    [string]$Mode = "smoke",
    [string]$Concentrations = "10,80",
    [int]$Nodes = 321,
    [int]$MaxIterations = 12000,
    [int]$CallbackStride = 10,
    [int]$MaxCallbackFrames = 8
)

$ErrorActionPreference = "Stop"

$repo = (Resolve-Path (Join-Path $PSScriptRoot "..\..\..")).Path
. (Join-Path $repo "scripts\runtime.ps1")
$julia = Resolve-JuliaExecutable
$project = Join-Path $repo "sim\julia"
$solver = Join-Path $PSScriptRoot "run_fig7_production_v7_10_0.jl"
$resolvedOutput = [System.IO.Path]::GetFullPath($Output)

if ($Nodes -notin @(321, 641)) {
    throw "V7.10.0 accepts only 321-node smoke or 641-node production meshes."
}
if ($CallbackStride -lt 1 -or $MaxCallbackFrames -lt 1) {
    throw "CallbackStride and MaxCallbackFrames must both be positive."
}
if (-not (Test-Path -LiteralPath $solver)) {
    throw "V7.10.0 solver not found: $solver"
}
if (Test-Path -LiteralPath $resolvedOutput) {
    $existing = Get-ChildItem -LiteralPath $resolvedOutput -Force
    if ($existing.Count -gt 0) {
        throw "Refusing to overwrite nonempty raw output directory: $resolvedOutput"
    }
}

$env:FIG7_V5_N_NODES = "$Nodes"
$env:FIG7_V5_KAPPA_LO_KBT = "80.0"
$env:FIG7_V5_KAPPA_LD_KBT = "25.0"
$env:FIG7_V5_KA_LO_J_M2 = "0.20"
$env:FIG7_V5_KA_LD_J_M2 = "0.20"
$env:FIG7_V5_KDELTA_LO_J_M2 = "0.20"
$env:FIG7_V5_KDELTA_LD_J_M2 = "0.20"
$env:FIG7_V5_NEUTRAL_OFFSET_UM = "0.002"
$env:FIG7_V5_KBT_J = "4.114e-21"
$env:FIG7_V5_PRE_SPHERE_TOL_UM = "0.25"

$env:FIG7_V770_MODE = $Mode
$env:FIG7_V770_RUN_TAG = "v7_10_0_gamma4"
$env:FIG7_V770_OUT = $resolvedOutput
$env:FIG7_V770_CONCENTRATIONS_NM = $Concentrations
$env:FIG7_V770_LINE_TENSION_PN = "0.3"
$env:FIG7_V770_GAUSSIAN_RATIO = "0.9"
$env:FIG7_V770_EPS_CTB_OUTER_LO = "0.0916"
$env:FIG7_V770_INNER_LD_BIOTIN_MOLE_FRACTION = "0.002"
$env:FIG7_V770_BOUND_SITE_AREA_COUPLING = "4.0"
$env:FIG7_V770_TARGET_C0_LD_COMBINED10_UM_INV = "-0.015"
$env:FIG7_V770_TARGET_C0_LO_POST_UM_INV = "0.3333333333333333"
$env:FIG7_V770_DOME_INCREASE_TARGET_UM = "2.95"
$env:FIG7_V770_DOME_INCREASE_TOL_UM = "0.50"
$env:FIG7_V770_STREP_SATURATION_NM = "80.0"
$env:FIG7_V770_MONOTONIC_TOL_UM = "0.02"
$env:FIG7_V770_INTERFACE_MOMENT_REL_TOL = "0.10"
$env:FIG7_V770_INTERFACE_TRANSVERSE_REL_TOL = "0.50"
$env:FIG7_V770_IPOPT_TOL = "1e-8"
$env:FIG7_V770_IPOPT_MAXITER = "$MaxIterations"

$env:FIG7_V710_CAPTURE_CALLBACK = "true"
$env:FIG7_V710_CALLBACK_STRIDE = "$CallbackStride"
$env:FIG7_V710_MAX_CALLBACK_FRAMES = "$MaxCallbackFrames"

& $julia "--project=$project" $solver
if ($LASTEXITCODE -ne 0) {
    throw "V7.10.0 raw solve failed with exit code $LASTEXITCODE."
}

@(
    [pscustomobject]@{ parameter = "workflow"; value = "V7.10.0" }
    [pscustomobject]@{ parameter = "mode"; value = $Mode }
    [pscustomobject]@{ parameter = "nodes"; value = "$Nodes" }
    [pscustomobject]@{ parameter = "gamma_A"; value = "4.0" }
    [pscustomobject]@{ parameter = "concentrations_nm"; value = $Concentrations }
    [pscustomobject]@{ parameter = "callback_stride"; value = "$CallbackStride" }
    [pscustomobject]@{ parameter = "maximum_callback_frames"; value = "$MaxCallbackFrames" }
    [pscustomobject]@{ parameter = "trajectory_interpretation"; value = "quasi-static minimized equilibria; optimizer iterates are not physical time" }
) | Export-Csv -LiteralPath (
    Join-Path $resolvedOutput "v7_10_wrapper_config.csv") -NoTypeInformation

$driverStatus = if ($Mode -eq "audit") {
    "audit_completed"
} else {
    "solver_completed"
}
Set-Content -LiteralPath (
    Join-Path $resolvedOutput "v7_10_driver_status.txt") `
    -Value $driverStatus
