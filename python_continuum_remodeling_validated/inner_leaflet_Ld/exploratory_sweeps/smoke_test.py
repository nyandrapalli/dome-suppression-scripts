"""
Smoke Test for Two-Phase Vesicle Remodeling Simulation Suite
Validates bare vesicle (theta = 0.0) equilibrium and energy terms.
"""

import os
import csv
import numpy as np
from sim.mechanics import MembraneParameters
from sim.solver import VesicleShapeSolver

def run_smoke_test():
    os.makedirs("simulation_results", exist_ok=True)
    trace_path = os.path.join("simulation_results", "trace_smoke_test.csv")
    
    print("==========================================================")
    print("RUNNING BIOPHYSICAL SMOKE TEST: BARE VESICLE (theta = 0.0)")
    print("==========================================================")
    
    params = MembraneParameters()
    solver = VesicleShapeSolver(params=params, n_points=100)
    
    # Solve unperturbed sphere branch at theta = 0.0
    res_sphere = solver.optimize_shape(theta=0.0, branch="sphere", maxiter=200)
    
    # Write trace CSV
    headers = [
        "theta", "branch", "kappa_Lo", "kappa_Ld", "c0_Ld", "gamma",
        "E_total", "E_bend_Lo", "E_bend_Ld", "E_line", "E_ade",
        "r_b", "dimple_depth", "res_A_tot", "res_V", "opt_success"
    ]
    
    row = [
        res_sphere['theta'],
        res_sphere['branch'],
        res_sphere['props']['kappa_Lo'],
        res_sphere['props']['kappa_Ld'],
        res_sphere['props']['c0_Ld'],
        res_sphere['props']['gamma'],
        res_sphere['total_energy'],
        res_sphere['E_bend_Lo'],
        res_sphere['E_bend_Ld'],
        res_sphere['E_line'],
        res_sphere['E_ade'],
        res_sphere['r_b'],
        res_sphere['dimple_depth'],
        res_sphere['res_A_tot'],
        res_sphere['res_V'],
        res_sphere['opt_success']
    ]
    
    with open(trace_path, mode="w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerow(row)
        
    print(f"Smoke test completed successfully. Trace saved to: {trace_path}")
    print(f"Total Energy: {res_sphere['total_energy']:.3f} k_B T")
    print(f"E_bend_Lo: {res_sphere['E_bend_Lo']:.3f} k_B T, E_bend_Ld: {res_sphere['E_bend_Ld']:.3f} k_B T")
    print(f"E_line: {res_sphere['E_line']:.3f} k_B T, E_ade: {res_sphere['E_ade']:.3f} k_B T")
    print(f"Boundary radius r_b: {res_sphere['r_b']:.3f} um")
    print(f"Dimple depth: {res_sphere['dimple_depth']:.3f} um")
    print(f"Area residual: {res_sphere['res_A_tot']:.2e}, Volume residual: {res_sphere['res_V']:.2e}")
    print(f"Optimization Success: {res_sphere['opt_success']}")
    print("==========================================================")

if __name__ == "__main__":
    run_smoke_test()
