"""
Production Parameter Sweeps for Single-Domain HermesLoouter (Protein on outer leaflet of Lo).
"""

import os
import sys
import time
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(__file__))
from mechanics import MembraneParameters, get_effective_properties
from solver import HermesSingleSolver

GAMMA_GRID = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
D_LO_GRID  = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
V_GRID     = [0.90, 0.95, 0.98, 1.00]
THETA_GRID = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]

def _solve_single_case(task):
    gamma_pN, D_Lo, v_val, theta, sweep_name, prof_file, out_dir = task
    
    params = MembraneParameters(
        R_ves=5.0,
        D_Lo=D_Lo,
        gamma_0=gamma_pN * 1e-12,
        reduced_volume=v_val
    )
    solver = HermesSingleSolver(params=params, n_points=120)
    
    branch = "sphere" if theta == 0.0 else "auto"
    res = solver.optimize_shape(theta=theta, branch=branch, maxiter=350, ftol=1e-9)
    
    # Save spatial profile CSV
    if prof_file:
        prof_path = os.path.join(out_dir, "profiles", prof_file)
        os.makedirs(os.path.dirname(prof_path), exist_ok=True)
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
        
    props = get_effective_properties(params, theta)
    
    return {
        "sweep": sweep_name,
        "gamma_pN": gamma_pN,
        "D_Lo_um": D_Lo,
        "reduced_volume_v": v_val,
        "theta": theta,
        "c0_Lo_um_inv": props["c0_Lo"],
        "c0_Ld_um_inv": props["c0_Ld"],
        "kappa_Lo_kBT": props["kappa_Lo"],
        "kappa_Ld_kBT": props["kappa_Ld"],
        "branch": res.get("branch", branch),
        "dimple_depth_um": res["dimple_depth"],
        "cap_elevation_um": res["cap_elevation"],
        "r_b_um": res["r_b"],
        "neck_aperture_2rb_um": 2.0 * res["r_b"],
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

def run_all_single_domain_sweeps(out_dir="HermesLoouter/single", max_workers=None):
    if max_workers is None:
        max_workers = os.cpu_count() or 4
        
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(os.path.join(out_dir, "profiles"), exist_ok=True)
    
    print("=" * 80)
    print("STARTING HERMES SINGLE-DOMAIN PRODUCTION PARAMETER SWEEPS (Lo Outer Leaflet)")
    print(f"Parallel Workers: {max_workers} CPU cores")
    print(f"1. Gamma sweep:       {len(GAMMA_GRID)} values in [0.5, 5.0] pN x {len(THETA_GRID)} thetas = {len(GAMMA_GRID)*len(THETA_GRID)} runs")
    print(f"2. Domain Size sweep: {len(D_LO_GRID)} values in [1.0, 7.0] um x {len(THETA_GRID)} thetas = {len(D_LO_GRID)*len(THETA_GRID)} runs")
    print(f"3. Deflation sweep:   {len(V_GRID)} values in [0.90, 1.00]   x {len(THETA_GRID)} thetas = {len(V_GRID)*len(THETA_GRID)} runs")
    print(f"Total production runs: {len(GAMMA_GRID)*len(THETA_GRID) + len(D_LO_GRID)*len(THETA_GRID) + len(V_GRID)*len(THETA_GRID)}")
    print("=" * 80)
    
    tasks = []
    
    # 1. Line Tension Sweep (D_Lo = 3.0 um, v = 0.98)
    for g in GAMMA_GRID:
        for th in THETA_GRID:
            p_file = f"prof_gamma_{g:.1f}_theta_{th:.2f}.csv"
            tasks.append((g, 3.0, 0.98, th, "gamma", p_file, out_dir))
            
    # 2. Domain Size Sweep (gamma = 1.5 pN, v = 0.98)
    for d in D_LO_GRID:
        for th in THETA_GRID:
            p_file = f"prof_DLo_{d:.1f}_theta_{th:.2f}.csv"
            tasks.append((1.5, d, 0.98, th, "D_Lo", p_file, out_dir))
            
    # 3. Deflation Sweep (gamma = 1.5 pN, D_Lo = 3.0 um)
    for v_val in V_GRID:
        for th in THETA_GRID:
            p_file = f"prof_v_{v_val:.2f}_theta_{th:.2f}.csv"
            tasks.append((1.5, 3.0, v_val, th, "deflation", p_file, out_dir))
            
    t0 = time.time()
    results = []
    
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(_solve_single_case, t): t for t in tasks}
        completed = 0
        total = len(tasks)
        for future in as_completed(futures):
            res = future.result()
            results.append(res)
            completed += 1
            if completed % 10 == 0 or completed == total:
                elapsed = time.time() - t0
                print(f"Progress: [{completed:3d}/{total:3d}] ({completed/total*100:5.1f}%) | Elapsed: {elapsed:5.1f}s | Last: {res['sweep']} g={res['gamma_pN']:.1f} D={res['D_Lo_um']:.1f} v={res['reduced_volume_v']:.2f} th={res['theta']:.2f} -> Elev={res['cap_elevation_um']:+.3f} um")

    df = pd.DataFrame(results)
    
    # Save partitioned CSVs
    df_gamma = df[df["sweep"] == "gamma"].sort_values(["gamma_pN", "theta"])
    df_dlo   = df[df["sweep"] == "D_Lo"].sort_values(["D_Lo_um", "theta"])
    df_v     = df[df["sweep"] == "deflation"].sort_values(["reduced_volume_v", "theta"])
    
    df_gamma.to_csv(os.path.join(out_dir, "trace_gamma_sweep.csv"), index=False)
    df_dlo.to_csv(os.path.join(out_dir, "trace_domain_size_sweep.csv"), index=False)
    df_v.to_csv(os.path.join(out_dir, "trace_deflation_sweep.csv"), index=False)
    df.to_csv(os.path.join(out_dir, "trace_all_single_domain.csv"), index=False)
    
    t_tot = time.time() - t0
    print("=" * 80)
    print(f"HERMES SINGLE-DOMAIN SWEEPS COMPLETED in {t_tot:.1f} seconds!")
    print(f"Saved: {os.path.join(out_dir, 'trace_gamma_sweep.csv')}")
    print(f"Saved: {os.path.join(out_dir, 'trace_domain_size_sweep.csv')}")
    print(f"Saved: {os.path.join(out_dir, 'trace_deflation_sweep.csv')}")
    print(f"Saved: {os.path.join(out_dir, 'trace_all_single_domain.csv')}")
    print("=" * 80)
    return df

if __name__ == "__main__":
    run_all_single_domain_sweeps()
