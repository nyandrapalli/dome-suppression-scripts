"""
Smoke test for Robust Size Solver across D = 2.0 um.
"""

import os
import csv
from sim.size_sweep.size_solver_v2 import RobustSizeVesicleSolver

def run_smoke_test_size_v2():
    results_dir = os.path.join("simulation_results", "size_sweep")
    os.makedirs(results_dir, exist_ok=True)
    trace_path = os.path.join(results_dir, "trace_smoke_test_size_v2.csv")
    
    print("====================================================================")
    print("RUNNING ROBUST SIZE SWEEP SMOKE TEST: D_Lo = 2.0 um (theta=0 & 1.0)")
    print("====================================================================")
    
    solver = RobustSizeVesicleSolver(D_Lo=2.0, R_ves=5.0, n_points=120)
    
    res_bare = solver.optimize_shape(theta=0.0, branch="sphere", maxiter=300)
    res_inward = solver.optimize_shape(theta=1.0, branch="inward", maxiter=300)
    
    headers = [
        "D_Lo", "theta", "branch", "kappa_Ld", "E_total", "E_line", "E_bend_Lo", "E_ade",
        "r_b", "dimple_depth", "measured_depth", "res_A_tot", "res_V", "opt_success"
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
                r['r_b'],
                r['dimple_depth'],
                r['measured_dimple_depth'],
                r['res_A_tot'],
                r['res_V'],
                r['opt_success']
            ])
            
    print(f"Smoke test completed. Trace saved to: {trace_path}")
    print(f"Bare (theta=0) Total Energy: {res_bare['total_energy']:.2f} kBT, r_b: {res_bare['r_b']:.3f} um")
    print(f"Inward (theta=1) Total Energy: {res_inward['total_energy']:.2f} kBT, Dimple Depth: {res_inward['measured_dimple_depth']:.3f} um")
    print(f"Optimization Success: {res_bare['opt_success']} and {res_inward['opt_success']}")
    print("====================================================================")

if __name__ == "__main__":
    run_smoke_test_size_v2()
