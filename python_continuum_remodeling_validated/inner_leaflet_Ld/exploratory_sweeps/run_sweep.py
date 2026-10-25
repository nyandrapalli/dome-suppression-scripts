"""
Production Sweep Script: Remodeling of Two-Phase GUV with Streptavidin Inner-Leaflet Binding.
Evaluates energy landscapes across coverage theta in [0.0, 1.0] for:
  - Branch 1: Sphere / Quasi-flat Lo domain
  - Branch 2: Inward Dimple (Invaginated Lo domain)
  - Branch 3: Outward Bud (Bulged Lo domain)
"""

import os
import csv
import numpy as np
from sim.mechanics import MembraneParameters
from sim.solver import VesicleShapeSolver

def run_production_sweep():
    os.makedirs("simulation_results", exist_ok=True)
    trace_path = os.path.join("simulation_results", "trace_remodeling_sweep.csv")
    
    params = MembraneParameters()
    solver = VesicleShapeSolver(params=params, n_points=120)
    
    theta_list = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    branches = ["sphere", "inward", "outward"]
    
    results_all = []
    
    print("==========================================================================")
    print("STARTING FULL BIOPHYSICAL PARAMETRIC SWEEP (theta: 0.0 -> 1.0)")
    print("==========================================================================")
    print(f"{'Theta':<6} | {'Branch':<8} | {'kappa_Ld':<8} | {'E_total (kBT)':<14} | {'E_line':<10} | {'E_bend_Lo':<10} | {'E_ade':<8} | {'r_b (um)':<8} | {'Dimple (um)':<10}")
    print("-" * 90)
    
    for theta in theta_list:
        for branch in branches:
            res = solver.optimize_shape(theta=theta, branch=branch, maxiter=400)
            results_all.append(res)
            
            # Save profile coordinate data
            profile_fn = f"profile_theta_{theta:.1f}_{branch}.csv"
            profile_path = os.path.join("simulation_results", profile_fn)
            with open(profile_path, "w", newline="") as pf:
                p_writer = csv.writer(pf)
                p_writer.writerow(["s", "psi", "r", "z", "H"])
                for s_val, psi_val, r_val, z_val, H_val in zip(
                    res['s_grid'], res['psi_grid'], res['r_profile'], res['z_profile'], res['H_profile']
                ):
                    p_writer.writerow([s_val, psi_val, r_val, z_val, H_val])
            
            print(
                f"{theta:<6.1f} | {branch:<8} | {res['props']['kappa_Ld']:<8.1f} | "
                f"{res['total_energy']:<14.2f} | {res['E_line']:<10.2f} | "
                f"{res['E_bend_Lo']:<10.2f} | {res['E_ade']:<8.2f} | "
                f"{res['r_b']:<8.3f} | {res['dimple_depth']:<10.3f}"
            )
            
    # Write summary trace CSV
    headers = [
        "theta", "branch", "kappa_Lo", "kappa_Ld", "c0_Ld", "gamma",
        "E_total", "E_bend_Lo", "E_bend_Ld", "E_line", "E_ade",
        "r_b", "z_b", "dimple_depth", "A_Lo", "A_Ld", "A_tot", "V_tot",
        "res_A_tot", "res_V", "opt_success"
    ]
    
    with open(trace_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        for r in results_all:
            writer.writerow([
                r['theta'],
                r['branch'],
                r['props']['kappa_Lo'],
                r['props']['kappa_Ld'],
                r['props']['c0_Ld'],
                r['props']['gamma'],
                r['total_energy'],
                r['E_bend_Lo'],
                r['E_bend_Ld'],
                r['E_line'],
                r['E_ade'],
                r['r_b'],
                r['z_b'],
                r['dimple_depth'],
                r['A_Lo'],
                r['A_Ld'],
                r['A_tot'],
                r['V_tot'],
                r['res_A_tot'],
                r['res_V'],
                r['opt_success']
            ])
            
    print("==========================================================================")
    print(f"Production sweep finished. Summary saved to: {trace_path}")
    print("==========================================================================")

if __name__ == "__main__":
    run_production_sweep()
