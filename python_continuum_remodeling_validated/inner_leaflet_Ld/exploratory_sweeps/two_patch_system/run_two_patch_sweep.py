"""
Production Runner for Two-Patch Antipodal Vesicle Sweep (D_Lo = 2.0 um to 7.0 um).
Evaluates dual North/South inward dimpling and biconcave discocyte remodeling.
"""

import os
import csv
import numpy as np
from sim.two_patch_system.two_patch_solver import TwoPatchVesicleSolver

def run_two_patch_sweep():
    results_dir = os.path.join("simulation_results", "two_patch_system")
    os.makedirs(results_dir, exist_ok=True)
    trace_path = os.path.join(results_dir, "trace_two_patch_sweep.csv")
    
    diameters = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
    thetas = [0.0, 0.5, 1.0]
    branches = ["sphere", "inward", "outward"]
    
    all_results = []
    
    print("=================================================================================================")
    print("STARTING TWO-PATCH ANTIPODAL SWEEP (D_Lo: 2.0 um to 7.0 um per patch)")
    print("=================================================================================================")
    print(f"{'D_Lo (um)':<10} | {'Theta':<6} | {'Branch':<8} | {'E_total (kBT)':<14} | {'E_line':<10} | {'E_ade':<8} | {'r_b1 (um)':<9} | {'Depth Mean':<12}")
    print("-" * 95)
    
    for D in diameters:
        solver = TwoPatchVesicleSolver(D_Lo=D, R_ves=5.0, n_points=140)
        for theta in thetas:
            for branch in branches:
                res = solver.optimize_shape(theta=theta, branch=branch, maxiter=450)
                all_results.append(res)
                
                p_fn = f"profile_twopatch_D{D:.1f}um_theta{theta:.1f}_{branch}.csv"
                p_path = os.path.join(results_dir, p_fn)
                with open(p_path, "w", newline="") as pf:
                    p_writer = csv.writer(pf)
                    p_writer.writerow(["s", "psi", "r", "z", "H", "n_Lo1", "n_Lo2"])
                    for s_v, psi_v, r_v, z_v, H_v in zip(
                        res['s_grid'], res['psi_grid'], res['r_profile'], res['z_profile'], res['H_profile']
                    ):
                        p_writer.writerow([s_v, psi_v, r_v, z_v, H_v, res['n_Lo1'], res['n_Lo2']])
                
                depth_m = res['dimple_depth_mean']
                print(
                    f"{D:<10.1f} | {theta:<6.1f} | {branch:<8} | "
                    f"{res['total_energy']:<14.2f} | {res['E_line']:<10.2f} | "
                    f"{res['E_ade']:<8.2f} | {res['r_b1']:<9.3f} | {depth_m:<12.3f}"
                )
                
    headers = [
        "D_Lo", "theta", "branch", "kappa_Lo", "kappa_Ld", "E_total",
        "E_bend_Lo", "E_bend_Ld", "E_line", "E_ade", "r_b1", "r_b2", "z_b1", "z_b2",
        "dimple_depth_north", "dimple_depth_south", "dimple_depth_mean",
        "A_Lo1", "A_Lo2", "A_Ld", "A_tot", "V_tot", "res_A_tot", "res_V", "opt_success",
        "n_Lo1", "n_Lo2"
    ]
    
    with open(trace_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        for r in all_results:
            writer.writerow([
                r['D_Lo'],
                r['theta'],
                r['branch'],
                r['props']['kappa_Lo'],
                r['props']['kappa_Ld'],
                r['total_energy'],
                r['E_bend_Lo'],
                r['E_bend_Ld'],
                r['E_line'],
                r['E_ade'],
                r['r_b1'],
                r['r_b2'],
                r['z_b1'],
                r['z_b2'],
                r['dimple_depth_north'],
                r['dimple_depth_south'],
                r['dimple_depth_mean'],
                r['A_Lo1'],
                r['A_Lo2'],
                r['A_Ld'],
                r['A_tot'],
                r['V_tot'],
                r['res_A_tot'],
                r['res_V'],
                r['opt_success'],
                r['n_Lo1'],
                r['n_Lo2']
            ])
            
    print("=================================================================================================")
    print(f"Two-patch sweep finished. Results saved to: {trace_path}")
    print("=================================================================================================")

if __name__ == "__main__":
    run_two_patch_sweep()
