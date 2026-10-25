"""
Two-Domain Fine Baseline Titration Sweep (PRL Standard).
Computes 33 simulations across 11 coverage steps theta in [0.0, 1.0] (step 0.1)
for 3 competing morphological branches ('sphere', 'inward', 'outward')
at baseline vesicle conditions (R_ves = 5.0 um, D_Lo = 3.0 um, gamma = 1.5 pN, v = 0.98).

Also includes comparison with deflated v = 0.90 state to analyze bifurcation energetics.
"""

import os
import sys
import time
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(__file__))

try:
    from mechanics_two_domains import TwoDomainParameters, TwoDomainEnergyCalculator
    from solver_two_domains import TwoDomainVesicleSolver
except ImportError:
    from sim_prl.mechanics_two_domains import TwoDomainParameters, TwoDomainEnergyCalculator
    from sim_prl.solver_two_domains import TwoDomainVesicleSolver

def run_single_titration_case(task):
    theta, branch, v_val, D_Lo, gamma, out_dir = task
    
    params = TwoDomainParameters(
        R_ves=5.0,
        D_Lo=D_Lo,
        gamma_0=gamma * 1e-12,
        reduced_volume=v_val
    )
    
    solver = TwoDomainVesicleSolver(params=params, n_points=101)
    
    t0 = time.time()
    res = solver.optimize_shape(theta=theta, branch=branch, maxiter=400)
    elapsed = time.time() - t0
    
    prof_name = f"prof_titr_v_{v_val:.2f}_theta_{theta:.2f}_{branch}.csv"
    prof_path = os.path.join(out_dir, "profiles", prof_name)
    
    # Save profile
    df_prof = pd.DataFrame({
        "s_um": solver.s_grid,
        "psi_rad": res["psi_profile"],
        "r_um": res["r_profile"],
        "z_um": res["z_profile"],
        "H_um_inv": res["H_profile"],
        "K_um_inv2": res["K_profile"]
    })
    domain_tags = []
    for s_val in solver.s_grid:
        if s_val <= solver.s_split1:
            domain_tags.append("Lo_north")
        elif s_val >= solver.s_split2:
            domain_tags.append("Lo_south")
        else:
            domain_tags.append("Ld")
    df_prof["domain"] = domain_tags
    df_prof.to_csv(prof_path, index=False)
    
    return {
        "theta": theta,
        "branch": branch,
        "reduced_volume_v": v_val,
        "D_Lo_um": D_Lo,
        "gamma_pN": gamma,
        "converged": res["converged"],
        "opt_success": res["opt_success"],
        "geom_converged": res["geom_converged"],
        "max_geom_error": res["max_geom_error"],
        "res_V": res["res_V"],
        "res_A_tot": res["res_A_tot"],
        "res_A_Lo": res["res_A_Lo1"],
        "total_energy_kBT": res["total_energy"],
        "E_bend_tot": res["E_bend_tot"],
        "E_bend_Lo_north": res["E_bend_Lo_north"],
        "E_bend_Ld": res["E_bend_Ld"],
        "E_bend_Lo_south": res["E_bend_Lo_south"],
        "E_line": res["E_line"],
        "E_gauss": res["E_gauss"],
        "E_ade": res["E_ade"],
        "depth_north_um": res["dimple_depth_north"],
        "depth_south_um": res["dimple_depth_south"],
        "r_b1_um": res["r_b1"],
        "r_b2_um": res["r_b2"],
        "iterations": res["iterations"],
        "elapsed_sec": elapsed,
        "profile_file": prof_name
    }

def run_titration_suite(smoke_only=False, max_workers=6):
    out_dir = "simulation_results_prl/two_domains_sweeps"
    prof_dir = os.path.join(out_dir, "profiles")
    os.makedirs(prof_dir, exist_ok=True)
    
    if smoke_only:
        print("=== RUNNING TWO-DOMAIN TITRATION SMOKE TEST ===")
        tasks = [
            (0.0, "sphere", 0.98, 3.0, 1.5, out_dir),
            (0.5, "inward", 0.98, 3.0, 1.5, out_dir),
            (1.0, "outward", 0.98, 3.0, 1.5, out_dir)
        ]
        results = [run_single_titration_case(t) for t in tasks]
        df_smoke = pd.DataFrame(results)
        smoke_file = os.path.join(out_dir, "trace_two_domains_titration_smoke.csv")
        df_smoke.to_csv(smoke_file, index=False)
        print(f"Smoke test results saved to: {smoke_file}")
        print(df_smoke[["theta", "branch", "converged", "max_geom_error", "total_energy_kBT", "depth_north_um"]].to_string())
        return df_smoke
        
    print("=== RUNNING PRODUCTION 33-POINT TWO-DOMAIN BASE TITRATION ===")
    thetas = np.round(np.linspace(0.0, 1.0, 11), 2)
    branches = ["sphere", "inward", "outward"]
    
    tasks = []
    for th in thetas:
        for br in branches:
            tasks.append((th, br, 0.98, 3.0, 1.5, out_dir))
            
    print(f"Total tasks in queue: {len(tasks)}")
    
    results = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(run_single_titration_case, t): t for t in tasks}
        for future in as_completed(futures):
            res = future.result()
            results.append(res)
            print(f"  [Done] theta={res['theta']:.2f}, branch={res['branch']:<7} | E={res['total_energy_kBT']:7.1f} kBT | d_N={res['depth_north_um']:+.3f} um | err={res['max_geom_error']:.1e}")
            
    df_prod = pd.DataFrame(results)
    df_prod.sort_values(by=["theta", "branch"], inplace=True)
    out_csv = os.path.join(out_dir, "trace_two_domains_titration.csv")
    df_prod.to_csv(out_csv, index=False)
    print(f"\nProduction titration trace saved to: {out_csv}")
    return df_prod

if __name__ == "__main__":
    is_smoke = "--smoke" in sys.argv
    run_titration_suite(smoke_only=is_smoke)
