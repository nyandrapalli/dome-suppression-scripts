"""
Smoke test for Two-Patch Antipodal System (D_Lo = 3.0 um, theta = 0.0 & 1.0).
"""

import os
import csv
from sim.two_patch_system.two_patch_solver import TwoPatchVesicleSolver

def run_smoke_test_two_patch():
    results_dir = os.path.join("simulation_results", "two_patch_system")
    os.makedirs(results_dir, exist_ok=True)
    trace_path = os.path.join(results_dir, "trace_smoke_test_two_patch.csv")
    
    print("===========================================================================")
    print("RUNNING TWO-PATCH SMOKE TEST: D_Lo = 3.0 um (North & South Lo, theta=0 & 1)")
    print("===========================================================================")
    
    solver = TwoPatchVesicleSolver(D_Lo=3.0, R_ves=5.0, n_points=140)
    
    res_bare = solver.optimize_shape(theta=0.0, branch="sphere", maxiter=300)
    res_inward = solver.optimize_shape(theta=1.0, branch="inward", maxiter=300)
    
    headers = [
        "D_Lo", "theta", "branch", "kappa_Ld", "E_total", "E_line", "E_bend_Lo", "E_ade",
        "r_b1", "r_b2", "depth_north", "depth_south", "depth_mean", "res_A_tot", "res_V", "opt_success"
    ]
    
    with open(trace_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        for r in [res_bare, res_inward]:
            writer.writerow([
                r['D_Lo'],
                r['theta'],
                r['branch'],
                r['props']['kappa_Ld'],
                r['total_energy'],
                r['E_line'],
                r['E_bend_Lo'],
                r['E_ade'],
                r['r_b1'],
                r['r_b2'],
                r['dimple_depth_north'],
                r['dimple_depth_south'],
                r['dimple_depth_mean'],
                r['res_A_tot'],
                r['res_V'],
                r['opt_success']
            ])
            
    print(f"Smoke test completed. Trace saved to: {trace_path}")
    print(f"Bare (theta=0) Total Energy: {res_bare['total_energy']:.2f} kBT, r_b1: {res_bare['r_b1']:.3f} um, r_b2: {res_bare['r_b2']:.3f} um")
    print(f"Inward (theta=1) Total Energy: {res_inward['total_energy']:.2f} kBT, North Depth: {abs(res_inward['dimple_depth_north']):.3f} um, South Depth: {abs(res_inward['dimple_depth_south']):.3f} um")
    print(f"Optimization Success: {res_bare['opt_success']} and {res_inward['opt_success']}")
    print("===========================================================================")

if __name__ == "__main__":
    run_smoke_test_two_patch()
