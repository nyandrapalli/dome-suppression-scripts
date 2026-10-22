"""
One-case Smoke Test for PRL Variational Shape Solver.
Audits the physical parameters, exact constraint residuals, and energy components
for baseline (theta = 0.0) and partially covered (theta = 0.5) conditions.
"""

import os
import sys
import pandas as pd
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))

try:
    from mechanics import MembraneParameters
    from solver import PRLVesicleSolver
except ImportError:
    from sim_prl.mechanics import MembraneParameters
    from sim_prl.solver import PRLVesicleSolver

def run_smoke_test(output_dir="simulation_results_prl"):
    os.makedirs(output_dir, exist_ok=True)
    params = MembraneParameters()
    solver = PRLVesicleSolver(params=params, n_points=140, n_ctrl_Lo=7, n_ctrl_Ld=11)
    
    test_conditions = [
        {"theta": 0.0, "branch": "sphere"},
        {"theta": 0.0, "branch": "inward"},
        {"theta": 0.5, "branch": "sphere"},
        {"theta": 0.5, "branch": "inward"},
    ]
    
    results = []
    print("=" * 80)
    print("RUNNING PRL PYTHON SOLVER SMOKE TEST")
    print("=" * 80)
    
    for cond in test_conditions:
        theta = cond["theta"]
        branch = cond["branch"]
        print(f"Solving condition: theta = {theta:.2f}, branch = {branch} ...")
        
        res = solver.optimize_shape(theta=theta, branch=branch, maxiter=500)
        
        row = {
            "theta": theta,
            "branch": branch,
            "converged": res["converged"],
            "max_geom_error": res["max_geom_error"],
            "res_V": res["res_V"],
            "res_A_Lo": res["res_A_Lo"],
            "res_A_tot": res["res_A_tot"],
            "total_energy_kBT": res["total_energy"],
            "E_bend_Lo": res["E_bend_Lo"],
            "E_bend_Ld": res["E_bend_Ld"],
            "E_line": res["E_line"],
            "E_gauss": res["E_gauss"],
            "E_ade": res["E_ade"],
            "r_b_um": res["r_b"],
            "z_b_um": res["z_b"],
            "psi_b_rad": res["psi_b"],
            "dimple_depth_um": res["dimple_depth"],
            "moment_residual": res["moment_residual"],
            "iterations": res["iterations"]
        }
        results.append(row)
        print(f"  -> Total E: {res['total_energy']:.2f} kBT | Max Geom Error: {res['max_geom_error']:.2e} | Dimple Depth: {res['dimple_depth']:.3f} um")
        
    df = pd.DataFrame(results)
    out_csv = os.path.join(output_dir, "smoke_test_prl.csv")
    df.to_csv(out_csv, index=False)
    print("=" * 80)
    print(f"Smoke test completed successfully. Output saved to {out_csv}")
    print("=" * 80)
    return df

if __name__ == "__main__":
    run_smoke_test()
