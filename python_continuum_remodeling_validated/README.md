# Python Continuum Vesicle Remodeling Sweeps

This directory contains the streamlined Python continuum mechanics solvers, multi-parameter sweeps, and compact summary data tables for two-phase ($L_o / L_d$) vesicle remodeling driven by asymmetric protein binding.

## Directory Structure

```text
python_continuum_remodeling/
├── inner_leaflet_Ld/                # Protein binding on the INNER leaflet of the Ld domain (PRL standard)
│   ├── run_all_sweeps.py            # Master runner for single- and two-domain expanded sweeps
│   ├── audit_all_traces.py          # Equilibrium & geometric constraint audit across sweep traces
│   ├── generate_comparison_markdown.py
│   ├── singledomain/                # Single Lo cap + Ld bulk solver, parameter sweeps, and plotting
│   ├── twodomain/                   # Two antipodal Lo caps + equatorial Ld belt solver and sweeps
│   ├── scaling_and_3d/              # Vesicle size scaling analysis, phase diagrams, and 3D mesh exporters
│   ├── exploratory_sweeps/          # Exploratory sweeps (large domains, plateau scaling, two-patch, size sweeps)
│   └── results/                     # Compact summary trace CSVs and publication figures
│       ├── param_sweeps/            # Single-domain gamma, D_Lo, and reduced-volume (v) sweep traces
│       ├── two_domains_sweeps/      # Two-domain gamma, D_Lo, v, and titration sweep traces
│       ├── scaling_analysis/        # Size-scaling trace CSV and figures
│       └── exploratory/             # Summary traces and figures from exploratory sweeps
└── outer_leaflet_Lo/                # Control/comparison: Protein binding on the OUTER leaflet of the Lo domain
    ├── run_all_sweeps.py            # Master runner for outer-leaflet Lo sweeps
    ├── smoke_test.py                # Fast verification test
    ├── plot_comparison_single_two.py
    ├── plot_local_saddle_bending.py
    ├── plot_comparison_2d_vs_3d.py  # 2D axisymmetric continuum vs. 3D triangulated mesh comparison
    ├── singledomain/                # Single-domain outer-leaflet Lo mechanics and solver
    ├── twodomain/                   # Two-domain outer-leaflet Lo mechanics and solver
    └── results/                     # Summary trace CSVs and comparative figures
        ├── single/
        └── double/
```

## Physical Models

1. **Inner-Leaflet $L_d$ Binding (`inner_leaflet_Ld/`)**:
   - **Canham–Helfrich Bending Elasticity**: Stiff $L_o$ cap ($\kappa_{Lo} = 80\,k_B T$) and protein-stiffened $L_d$ domain ($\kappa_{Ld}(\theta) = 20 + 70\,\theta\,k_B T$, spontaneous curvature $c_0(\theta) = -\alpha_{c0}\theta$).
   - **Non-Local Area-Difference Elasticity (ADE)**: Coupled to inner-leaflet expansion $\eta_{\text{inner}} = 0.08\,\theta$.
   - **Jülicher–Lipowsky Line Tension & Gauss–Bonnet Jump**: Boundary energy at the $L_o / L_d$ interface with regularized pole curvatures.
2. **Outer-Leaflet $L_o$ Binding (`outer_leaflet_Lo/`)**:
   - Protein binds to the outer leaflet of the $L_o$ domain ($\kappa_{Lo}(\theta) = 80 + 70\,\theta\,k_B T$, positive spontaneous curvature $c_0(\theta) = +\alpha_{c0}\theta$, outer-leaflet expansion $\eta_{\text{outer}} = 0.08\,\theta$) while the $L_d$ domain remains bare ($\kappa_{Ld} = 20\,k_B T$).

## Running the Sweeps

From `python_continuum_remodeling/inner_leaflet_Ld`:
```powershell
python run_all_sweeps.py
python audit_all_traces.py
```

From `python_continuum_remodeling/outer_leaflet_Lo`:
```powershell
python smoke_test.py
python run_all_sweeps.py
```

> **Note on Excluded Artifacts**: Raw per-node meridian spatial profiles (`profiles/prof_*.csv`), 3D `.obj`/`.mtl` meshes, and presentation `.mp4`/`.wav` media are excluded from version control to keep the repository compact. All summary `trace_*.csv` tables and publication `.png` figures are retained under `results/`.
