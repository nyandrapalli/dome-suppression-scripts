"""
One-case Smoke Test for Multi-Size Vesicle Scaling Analysis (PRL Standard).
Tests a 20 um diameter vesicle with a 6.0 um Lo domain at theta = 0.60, v = 0.98.
"""

import os
import sys
import pandas as pd
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))

from mechanics import MembraneParameters
from solver import PRLVesicleSolver

def run_scaling_smoke_test(out_dir="simulation_results_prl/scaling_analysis"):
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(os.path.join(out_dir, "profiles"), exist_ok=True)
    
    # Test case: 20 um vesicle (R = 10 um), 6.0 um Lo domain (R_Lo = 3.0 um)
    # Matching ratio D_Lo / D_ves = 0.30
    params = MembraneParameters(
        R_ves=10.0,
        D_Lo=6.0,
        gamma_0=1.5e-12,
        reduced_volume=0.98
    )
    
    solver = PRLVesicleSolver(params=params, n_points=140, n_ctrl_Lo=7, n_ctrl_Ld=11)
    
    print("=" * 80)
    print("RUNNING SCALING SMOKE TEST: D_ves = 20 um, D_Lo = 6.0 um, theta = 0.60")
    print(f"R_ves: {params.R_ves:.1f} um | D_Lo: {params.D_Lo:.1f} um | s_split: {params.s_split:.3f} um | S_tot: {params.S_tot:.3f} um")
    print(f"Target Area: {params.A_tot_target:.2f} um^2 | Target Lo Area: {params.A_Lo_target:.2f} um^2 | Target Volume: {params.V_tot_target:.2f} um^3")
    print("=" * 80)
    
    theta = 0.60
    branch = "inward"
    
    res = solver.optimize_shape(theta=theta, branch=branch, maxiter=500, ftol=1e-9)
    
    prof_file = f"prof_smoke_Dves_20.0_DLo_6.0_theta_{theta:.2f}.csv"
    prof_path = os.path.join(out_dir, "profiles", prof_file)
    df_prof = pd.DataFrame({
        "s_um": solver.s_grid,
        "psi_rad": res["psi"],
        "r_um": res["r"],
        "z_um": res["z"],
        "H_um_inv": res["H"],
        "K_um_inv2": res["K"]
    })
    df_prof["domain"] = np.where(solver.s_grid <= solver.s_split, "Lo", "Ld")
    df_prof.to_csv(prof_path, index=False)
    
    row = {
        "D_ves_um": 20.0,
        "D_Lo_um": 6.0,
        "theta": theta,
        "branch": branch,
        "converged": res["converged"],
        "opt_success": res["opt_success"],
        "geom_converged": res["geom_converged"],
        "max_geom_error": res["max_geom_error"],
        "res_V": res["res_V"],
        "res_A_Lo": res["res_A_Lo"],
        "res_A_tot": res["res_A_tot"],
        "dimple_depth_um": res["dimple_depth"],
        "norm_dimple_depth": res["dimple_depth"] / 6.0,
        "r_b_um": res["r_b"],
        "neck_aperture_2rb_um": 2.0 * res["r_b"],
        "norm_neck_aperture": (2.0 * res["r_b"]) / 6.0,
        "total_energy_kBT": res["total_energy"],
        "E_bend_tot": res["E_bend_tot"],
        "E_bend_Lo": res["E_bend_Lo"],
        "E_bend_Ld": res["E_bend_Ld"],
        "E_line": res["E_line"],
        "E_gauss": res["E_gauss"],
        "E_ade": res["E_ade"],
        "iterations": res["iterations"],
        "profile_file": prof_file
    }
    
    df = pd.DataFrame([row])
    out_csv = os.path.join(out_dir, "smoke_test_scaling.csv")
    df.to_csv(out_csv, index=False)
    
    print(f"Optimization Converged: {res['converged']} (Status: {res['solver_status']}, Geom Error: {res['max_geom_error']:.2e})")
    print(f"Dimple Depth: {res['dimple_depth']:+.3f} um (Normalized delta_z / D_Lo: {res['dimple_depth'] / 6.0:.3f})")
    print(f"Boundary Radius r_b: {res['r_b']:.3f} um (Neck Aperture 2*r_b / D_Lo: {(2.0 * res['r_b']) / 6.0:.3f})")
    print(f"Energy Breakdown: Bend_Lo={res['E_bend_Lo']:.2f}, Bend_Ld={res['E_bend_Ld']:.2f}, Line={res['E_line']:.2f}, ADE={res['E_ade']:.2f}, Gauss={res['E_gauss']:.2f}")
    print(f"Total Energy: {res['total_energy']:.2f} k_BT")
    print("=" * 80)
    print(f"Saved smoke test trace to: {out_csv}")
    print(f"Saved profile to: {prof_path}")
    print("=" * 80)
    return df

if __name__ == "__main__":
    run_scaling_smoke_test()
