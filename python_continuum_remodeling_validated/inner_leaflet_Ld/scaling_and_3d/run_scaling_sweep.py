"""
Comprehensive Scaling Analysis Sweep for Multi-Size Two-Phase Vesicles (PRL Standard).
Investigates whether geometric scaling (D_Lo / D_ves = 0.30) behaves linearly across:
  - 5 um vesicle with 1.5 um Lo domain (R_ves = 2.5 um, R_Lo = 0.75 um)
  - 10 um vesicle with 3.0 um Lo domain (R_ves = 5.0 um, R_Lo = 1.50 um)
  - 20 um vesicle with 6.0 um Lo domain (R_ves = 10.0 um, R_Lo = 3.00 um)
Across coverage theta in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0].
"""

import os
import sys
import time
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(__file__))

from mechanics import MembraneParameters
from solver import PRLVesicleSolver

VESICLE_SIZES = [
    {"name": "5um", "D_ves": 5.0, "R_ves": 2.5, "D_Lo": 1.5, "R_Lo": 0.75},
    {"name": "10um", "D_ves": 10.0, "R_ves": 5.0, "D_Lo": 3.0, "R_Lo": 1.50},
    {"name": "20um", "D_ves": 20.0, "R_ves": 10.0, "D_Lo": 6.0, "R_Lo": 3.00},
]

THETA_GRID = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

def _solve_scaling_case(task):
    size_info, theta, v_val, out_dir = task
    
    params = MembraneParameters(
        R_ves=size_info["R_ves"],
        D_Lo=size_info["D_Lo"],
        gamma_0=1.5e-12,
        reduced_volume=v_val
    )
    
    solver = PRLVesicleSolver(params=params, n_points=140, n_ctrl_Lo=7, n_ctrl_Ld=11)
    
    branch = "sphere" if theta == 0.0 else "inward"
    res = solver.optimize_shape(theta=theta, branch=branch, maxiter=500, ftol=1e-9)
    
    prof_file = f"prof_scale_Dves_{size_info['D_ves']:.1f}_DLo_{size_info['D_Lo']:.1f}_theta_{theta:.2f}.csv"
    prof_path = os.path.join(out_dir, "profiles", prof_file)
    os.makedirs(os.path.dirname(prof_path), exist_ok=True)
    
    df_prof = pd.DataFrame({
        "s_um": solver.s_grid,
        "s_norm": solver.s_grid / params.R_ves,
        "psi_rad": res["psi"],
        "r_um": res["r"],
        "r_norm": res["r"] / params.R_ves,
        "z_um": res["z"],
        "z_norm": res["z"] / params.R_ves,
        "H_um_inv": res["H"],
        "H_norm": res["H"] * params.R_ves,
        "K_um_inv2": res["K"]
    })
    df_prof["domain"] = np.where(solver.s_grid <= solver.s_split, "Lo", "Ld")
    df_prof.to_csv(prof_path, index=False)
    
    return {
        "size_name": size_info["name"],
        "D_ves_um": size_info["D_ves"],
        "R_ves_um": size_info["R_ves"],
        "D_Lo_um": size_info["D_Lo"],
        "R_Lo_um": size_info["R_Lo"],
        "theta": theta,
        "branch": branch,
        "reduced_volume_v": v_val,
        "dimple_depth_um": res["dimple_depth"],
        "norm_dimple_depth": res["dimple_depth"] / size_info["D_Lo"],
        "r_b_um": res["r_b"],
        "neck_aperture_2rb_um": 2.0 * res["r_b"],
        "norm_neck_aperture": (2.0 * res["r_b"]) / size_info["D_Lo"],
        "total_energy_kBT": res["total_energy"],
        "E_bend_tot": res["E_bend_tot"],
        "E_bend_Lo": res["E_bend_Lo"],
        "E_bend_Ld": res["E_bend_Ld"],
        "E_line": res["E_line"],
        "E_gauss": res["E_gauss"],
        "E_ade": res["E_ade"],
        "A_tot_um2": res["A_tot"],
        "V_tot_um3": res["V_tot"],
        "converged": res["converged"],
        "opt_success": res["opt_success"],
        "geom_converged": res["geom_converged"],
        "max_geom_error": res["max_geom_error"],
        "iterations": res["iterations"],
        "profile_file": prof_file
    }

def run_all_scaling_sweeps(out_dir="simulation_results_prl/scaling_analysis", max_workers=None):
    if max_workers is None:
        max_workers = os.cpu_count() or 4
        
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(os.path.join(out_dir, "profiles"), exist_ok=True)
    
    print("=" * 80)
    print("STARTING MULTI-SIZE VESICLE SCALING PRODUCTION SWEEPS")
    print(f"Parallel Workers: {max_workers} CPU cores")
    print(f"Sizes: 5 um (Lo 1.5 um), 10 um (Lo 3.0 um), 20 um (Lo 6.0 um)")
    print(f"Theta Grid: {len(THETA_GRID)} values in [0.0, 1.0]")
    print(f"Total production runs: {len(VESICLE_SIZES) * len(THETA_GRID)}")
    print("=" * 80)
    
    tasks = []
    for s_info in VESICLE_SIZES:
        for th in THETA_GRID:
            tasks.append((s_info, th, 0.98, out_dir))
            
    t0 = time.time()
    results = []
    
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(_solve_scaling_case, t): t for t in tasks}
        completed = 0
        total = len(tasks)
        for future in as_completed(futures):
            res = future.result()
            results.append(res)
            completed += 1
            elapsed = time.time() - t0
            print(f"Progress: [{completed:2d}/{total:2d}] ({completed/total*100:5.1f}%) | {res['size_name']:4s} th={res['theta']:.2f} -> Depth={res['dimple_depth_um']:+.3f} um (Norm: {res['norm_dimple_depth']:+.3f}) | Neck: {res['norm_neck_aperture']:.3f} | Converged={res['converged']}")

    df = pd.DataFrame(results).sort_values(["D_ves_um", "theta"]).reset_index(drop=True)
    
    trace_path = os.path.join(out_dir, "trace_scaling_sweep.csv")
    df.to_csv(trace_path, index=False)
    
    t_tot = time.time() - t0
    print("=" * 80)
    print(f"SCALING PRODUCTION SWEEPS COMPLETED in {t_tot:.1f} seconds!")
    print(f"Saved consolidated trace: {trace_path}")
    print("=" * 80)
    return df

if __name__ == "__main__":
    run_all_scaling_sweeps()
