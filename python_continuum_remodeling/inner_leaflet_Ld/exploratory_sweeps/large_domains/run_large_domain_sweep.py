"""
Production Runner for Large Lo Domain Sweep (D_Lo = 6.0 um, 7.0 um, 8.0 um).
Generates detailed trace CSV and profile files for large domain phase separation.
"""

import os
import csv
import numpy as np
from sim.large_domains.large_domain_solver import LargeDomainVesicleSolver

def run_large_domain_sweep():
    results_dir = os.path.join("simulation_results", "large_domains")
    os.makedirs(results_dir, exist_ok=True)
    trace_path = os.path.join(results_dir, "trace_large_domain_sweep.csv")
    
    diameters = [6.0, 7.0, 8.0]
    thetas = [0.0, 0.5, 1.0]
    branches = ["sphere", "inward", "outward"]
    
    all_results = []
    
    print("=================================================================================================")
    print("STARTING LARGE DOMAIN SWEEP (D_Lo: 6.0 um, 7.0 um, 8.0 um)")
    print("=================================================================================================")
    print(f"{'D_Lo (um)':<10} | {'Theta':<6} | {'Branch':<8} | {'E_total (kBT)':<14} | {'E_line':<10} | {'E_ade':<8} | {'r_b (um)':<8} | {'Dimple Depth':<12}")
    print("-" * 95)
    
    for D in diameters:
        solver = LargeDomainVesicleSolver(D_Lo=D, R_ves=5.0, n_points=140)
        for theta in thetas:
            for branch in branches:
                res = solver.optimize_shape(theta=theta, branch=branch, maxiter=450)
                all_results.append(res)
                
                # Save profile CSV
                p_fn = f"profile_large_D{D:.1f}um_theta{theta:.1f}_{branch}.csv"
                p_path = os.path.join(results_dir, p_fn)
                with open(p_path, "w", newline="") as pf:
                    p_writer = csv.writer(pf)
                    p_writer.writerow(["s", "psi", "r", "z", "H", "n_Lo"])
                    for s_v, psi_v, r_v, z_v, H_v in zip(
                        res['s_grid'], res['psi_grid'], res['r_profile'], res['z_profile'], res['H_profile']
                    ):
                        p_writer.writerow([s_v, psi_v, r_v, z_v, H_v, res['n_Lo']])
                
                depth = res['measured_dimple_depth']
                print(
                    f"{D:<10.1f} | {theta:<6.1f} | {branch:<8} | "
                    f"{res['total_energy']:<14.2f} | {res['E_line']:<10.2f} | "
                    f"{res['E_ade']:<8.2f} | {res['r_b']:<8.3f} | {depth:<12.3f}"
                )
                
    headers = [
        "D_Lo", "theta", "branch", "kappa_Lo", "kappa_Ld", "E_total",
        "E_bend_Lo", "E_bend_Ld", "E_line", "E_ade", "r_b", "z_b",
        "dimple_depth", "measured_dimple_depth", "A_Lo", "A_Ld", "A_tot", "V_tot",
        "res_A_tot", "res_V", "opt_success", "n_Lo"
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
                r['r_b'],
                r['z_b'],
                r['dimple_depth'],
                r['measured_dimple_depth'],
                r['A_Lo'],
                r['A_Ld'],
                r['A_tot'],
                r['V_tot'],
                r['res_A_tot'],
                r['res_V'],
                r['opt_success'],
                r['n_Lo']
            ])
            
    print("=================================================================================================")
    print(f"Large domain sweep finished. Results saved to: {trace_path}")
    print("=================================================================================================")

if __name__ == "__main__":
    run_large_domain_sweep()
