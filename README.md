# Dome Suppression Scripts

Reproducible Julia and Python workflow for the validated V7.10.0
concentration-dependent dome-suppression model in a phase-separated vesicle.

The model represents the membrane as an axisymmetric, two-domain meridian with
liquid-ordered (Lo) and liquid-disordered (Ld) regions. It minimizes a
two-phase membrane energy at each streptavidin concentration. The concentration
states are independent quasi-static equilibria; neither the concentration
sequence nor the saved optimizer iterates represent physical time.

## Validated release

- Release: V7.10.0
- Smoke mesh: 321 nodes
- Production mesh: 641 nodes
- Production concentrations: 10, 20, 30, 40, 50, 60, 70, and 80 nM
- Validated local production run: `production_passed`
- Color convention: Lo red (`#d9362b`), Ld green (`#168447`)

Small validation tables from the successful local production run are retained
under `reference/`. Raw meshes, videos, and optimizer frames are excluded
because they are generated artifacts.

## Repository layout

```text
sim/julia/
  Project.toml
  Manifest.toml
  fig7_geometry_contract.toml
  v5/run_fig7_fluid_domain_area_mismatch_v5.jl
  v7_7_0/run_fig7_local_leaflet_interface_v7_7_0.jl
  v7_10_0/
reference/
docs/model-audit.md
scripts/static-check.ps1
```

V7.10.0 loads V7.7.0, which loads V5. These files and their relative paths are
therefore part of the executable model and must remain together.

## Requirements

- Windows PowerShell 5.1 or PowerShell 7
- Julia 1.12
- Python 3.10 or newer
- FFmpeg on `PATH` for MP4 generation

Install dependencies from the repository root:

```powershell
julia --project=sim/julia -e 'using Pkg; Pkg.instantiate()'
python -m pip install -r requirements.txt
```

If Julia or Python is not on `PATH`, set `JULIA_EXE` or `PYTHON_EXE` to the
full executable path.

## Required execution order

First verify the repository without minimizing a membrane:

```powershell
& .\scripts\static-check.ps1
& .\sim\julia\v7_10_0\run_audit_v7_10_0.ps1
```

Inspect these audit outputs before a numerical run:

```text
sim/julia/v7_10_0/outputs/audit/condition_audit.csv
sim/julia/v7_10_0/outputs/audit/run_config.csv
sim/julia/v7_10_0/outputs/audit/state_table.csv
```

Run the 321-node smoke workflow:

```powershell
& .\sim\julia\v7_10_0\run_smoke_v7_10_0.ps1
```

Only after `outputs/smoke/workflow_status.txt` contains `smoke_passed`, run the
641-node production workflow:

```powershell
& .\sim\julia\v7_10_0\run_production_v7_10_0.ps1
```

Every output directory is overwrite-protected. Pass a new `-Output` path for a
new run rather than replacing an earlier result.

## Interpretation boundary

This repository establishes mechanistic, concentration-dependent suppression
of an Lo protrusion through a streptavidin-dependent preferred inner-Ld
leaflet-area shift. It does not contain the later exploratory dimple-forming
branches, and it must not be cited as a physical-time dynamics model.

The exact active and inactive terms are listed in
[`docs/model-audit.md`](docs/model-audit.md).
