"""
Production Titration Sweep for Two-Phase Vesicle Remodeling (PRL Standard)
Sweeps streptavidin coverage theta in [0.0, 1.0] across competing morphological branches
(sphere, inward dimple, outward bud).
Saves detailed trace CSV and per-condition (r, z, psi, H) coordinate profiles.
"""

import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))

try:
    from mechanics import MembraneParameters
    from solver import PRLVesicleSolver
except ImportError:
    from sim_prl.mechanics import MembraneParameters
    from sim_prl.solver import PRLVesicleSolver

def run_production_sweep(output_dir="simulation_results_prl"):
    os.makedirs(output_dir, exist_ok=True)
    profiles_dir = os.path.join(output_dir, "profiles")
    os.makedirs(profiles_dir, exist_ok=True)
    
    params = MembraneParameters()
    solver = PRLVesicleSolver(params=params, n_points=160, n_ctrl_Lo=8, n_ctrl_Ld=13)
    
    theta_values = np.linspace(0.0, 1.0, 11)  # 0.0, 0.1, 0.2, ..., 1.0
    branches = ["sphere", "inward", "outward"]
    
    trace_rows = []
    
    print("=" * 85)
    print("STARTING PRL PRODUCTION TITRATION SWEEP")
    print(f"Coverage Grid (11 points): {list(np.round(theta_values, 2))}")
    print(f"Competing Branches: {branches}")
    print("=" * 85)
    
    for theta in theta_values:
        for branch in branches:
            print(f"Solving theta = {theta:.2f} | branch = {branch:<8} ...", end=" ", flush=True)
            res = solver.optimize_shape(theta=theta, branch=branch, maxiter=700)
            
            # Save profile coordinate CSV
            prof_df = pd.DataFrame({
                "s_um": solver.s_grid,
                "r_um": res["r"],
                "z_um": res["z"],
                "psi_rad": res["psi"],
                "H_um_inv": res["H"],
                "domain": np.where(solver.s_grid <= solver.s_split, "Lo", "Ld")
            })
            prof_filename = f"profile_theta_{theta:.2f}_{branch}.csv"
            prof_path = os.path.join(profiles_dir, prof_filename)
            prof_df.to_csv(prof_path, index=False)
            
            # Collect trace row
            row = {
                "theta": round(float(theta), 3),
                "branch": branch,
                "converged": res["converged"],
                "opt_success": res["opt_success"],
                "geom_converged": res["geom_converged"],
                "max_geom_error": res["max_geom_error"],
                "res_V": res["res_V"],
                "res_A_Lo": res["res_A_Lo"],
                "res_A_tot": res["res_A_tot"],
                "total_energy_kBT": res["total_energy"],
                "E_bend_Lo": res["E_bend_Lo"],
                "E_bend_Ld": res["E_bend_Ld"],
                "E_bend_tot": res["E_bend_tot"],
                "E_line": res["E_line"],
                "E_gauss": res["E_gauss"],
                "E_ade": res["E_ade"],
                "r_b_um": res["r_b"],
                "z_b_um": res["z_b"],
                "psi_b_rad": res["psi_b"],
                "dimple_depth_um": res["dimple_depth"],
                "moment_residual": res["moment_residual"],
                "iterations": res["iterations"],
                "profile_file": prof_filename
            }
            trace_rows.append(row)
            print(f"Done. E_tot = {res['total_energy']:8.2f} kBT | Max err = {res['max_geom_error']:.2e} | Depth = {res['dimple_depth']:+.3f} um")
            
    trace_df = pd.DataFrame(trace_rows)
    trace_path = os.path.join(output_dir, "trace_remodeling_sweep_prl.csv")
    trace_df.to_csv(trace_path, index=False)
    print("=" * 85)
    print(f"Production sweep finished. Saved trace to {trace_path}")
    print("=" * 85)
    return trace_df

if __name__ == "__main__":
    run_production_sweep()
