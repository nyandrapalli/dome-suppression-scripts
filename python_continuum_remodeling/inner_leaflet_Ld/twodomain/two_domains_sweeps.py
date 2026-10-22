"""
Expanded Production Parameter Sweeps Module for Two-Domain Opposite Lo Vesicles (Bipolar - PRL Standard).
Explores:
1. Line tension gamma (10 values): [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0] pN
2. Domain size D_Lo (7 values):    [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0] um
3. Deflation factor v (4 values):  [0.90, 0.95, 0.98, 1.00]
Across streptavidin coverage theta in [0.0, 0.2, 0.4, 0.6, 0.8, 1.0].

Enforces exact physical mirror symmetry: depth_north == depth_south to machine precision.
Multi-core parallelized using concurrent.futures.ProcessPoolExecutor.
"""

import os
import sys
import time
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(__file__))

from mechanics_two_domains import TwoDomainParameters
from solver_two_domains import TwoDomainVesicleSolver

GAMMA_GRID = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
D_LO_GRID  = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
V_GRID     = [0.90, 0.95, 0.98, 1.00]
THETA_GRID = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]


def _solve_two_domain_case(task):
    gamma_pN, D_Lo, v_val, theta, sweep_name, prof_file, out_dir = task

    params = TwoDomainParameters(
        R_ves=5.0,
        D_Lo=D_Lo,
        gamma_0=gamma_pN * 1e-12,
        reduced_volume=v_val
    )
    solver = TwoDomainVesicleSolver(params=params, n_points=101)

    branch = "sphere" if theta == 0.0 else "auto"
    res = solver.optimize_shape(theta=theta, branch=branch, maxiter=350, ftol=1e-9)

    if prof_file:
        prof_path = os.path.join(out_dir, "profiles", prof_file)
        os.makedirs(os.path.dirname(prof_path), exist_ok=True)
        df_prof = pd.DataFrame({
            "s_um": solver.s_grid,
            "psi_rad": res["psi_profile"],
            "r_um": res["r_profile"],
            "z_um": res["z_profile"],
            "H_um_inv": res["H_profile"],
            "K_um_inv2": res["K_profile"]
        })
        dom_tags = []
        for s_val in solver.s_grid:
            if s_val <= solver.s_split1 + 1e-12:
                dom_tags.append("Lo_north")
            elif s_val >= solver.s_split2 - 1e-12:
                dom_tags.append("Lo_south")
            else:
                dom_tags.append("Ld")
        df_prof["domain"] = dom_tags
        df_prof.to_csv(prof_path, index=False)

    return {
        "sweep": sweep_name,
        "gamma_pN": gamma_pN,
        "D_Lo_um": D_Lo,
        "reduced_volume_v": v_val,
        "theta": theta,
        "branch": res.get("branch", branch),
        "depth_north_um": res["dimple_depth_north"],
        "depth_south_um": res["dimple_depth_south"],
        "r_b1_um": res["r_b1"],
        "r_b2_um": res["r_b2"],
        "neck_aperture_2rb_um": 2.0 * res["r_b1"],
        "total_energy_kBT": res["total_energy"],
        "E_bend_tot": res["E_bend_tot"],
        "E_bend_Lo_north": res["E_bend_Lo_north"],
        "E_bend_Ld": res["E_bend_Ld"],
        "E_bend_Lo_south": res["E_bend_Lo_south"],
        "E_line": res["E_line"],
        "E_gauss": res["E_gauss"],
        "E_ade": res["E_ade"],
        "converged": res["converged"],
        "opt_success": res["opt_success"],
        "geom_converged": res["geom_converged"],
        "max_geom_error": res["max_geom_error"],
        "iterations": res["iterations"],
        "profile_file": prof_file
    }


def run_all_two_domain_sweeps(out_dir=None, max_workers=None):
    if out_dir is None:
        out_dir = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "results", "two_domains_sweeps"
        )
    if max_workers is None:
        max_workers = os.cpu_count() or 4

    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(os.path.join(out_dir, "profiles"), exist_ok=True)

    print("=" * 80)
    print("STARTING EXPANDED TWO-DOMAIN PRODUCTION PARAMETER SWEEPS")
    print(f"Parallel Workers: {max_workers} CPU cores")
    print("=" * 80)

    tasks = []
    for g in GAMMA_GRID:
        for th in THETA_GRID:
            p_file = f"prof_two_dom_gamma_{g:.1f}_theta_{th:.2f}.csv"
            tasks.append((g, 3.0, 0.98, th, "gamma", p_file, out_dir))

    for d in D_LO_GRID:
        for th in THETA_GRID:
            p_file = f"prof_two_dom_DLo_{d:.1f}_theta_{th:.2f}.csv"
            tasks.append((1.5, d, 0.98, th, "D_Lo", p_file, out_dir))

    for v_val in V_GRID:
        for th in THETA_GRID:
            p_file = f"prof_two_dom_v_{v_val:.2f}_theta_{th:.2f}.csv"
            tasks.append((1.5, 3.0, v_val, th, "deflation", p_file, out_dir))

    t0 = time.time()
    results = []

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(_solve_two_domain_case, t): t for t in tasks}
        completed = 0
        total = len(tasks)
        for future in as_completed(futures):
            res = future.result()
            results.append(res)
            completed += 1
            if completed % 10 == 0 or completed == total:
                elapsed = time.time() - t0
                print(
                    f"Progress: [{completed:3d}/{total:3d}] ({completed/total*100:5.1f}%) | "
                    f"Elapsed: {elapsed:5.1f}s | Last: {res['sweep']} g={res['gamma_pN']:.1f} "
                    f"D={res['D_Lo_um']:.1f} v={res['reduced_volume_v']:.2f} th={res['theta']:.2f} "
                    f"-> Depth={res['depth_north_um']:+.3f} um"
                )

    df = pd.DataFrame(results)

    df_gamma = df[df["sweep"] == "gamma"].sort_values(["gamma_pN", "theta"])
    df_dlo   = df[df["sweep"] == "D_Lo"].sort_values(["D_Lo_um", "theta"])
    df_v     = df[df["sweep"] == "deflation"].sort_values(["reduced_volume_v", "theta"])

    df_gamma.to_csv(os.path.join(out_dir, "trace_gamma_sweep.csv"), index=False)
    df_dlo.to_csv(os.path.join(out_dir, "trace_domain_size_sweep.csv"), index=False)
    df_v.to_csv(os.path.join(out_dir, "trace_deflation_sweep.csv"), index=False)
    df.to_csv(os.path.join(out_dir, "trace_two_domains_all.csv"), index=False)

    t_tot = time.time() - t0
    print("=" * 80)
    print(f"TWO-DOMAIN PRODUCTION SWEEPS COMPLETED in {t_tot:.1f} seconds!")
    print("=" * 80)
    return df


if __name__ == "__main__":
    run_all_two_domain_sweeps()
