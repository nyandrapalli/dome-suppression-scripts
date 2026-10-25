"""
One-case Smoke Test for PRL Variational Shape Solver (Inner-Leaflet Ld Binding).
Audits the physical parameters, exact constraint residuals (including South-pole closure r(S_tot)=0),
and energy components for baseline (theta = 0.0) and partially covered (theta = 0.5) conditions.
"""

import os
import sys
import pandas as pd
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))

from mechanics import MembraneParameters, get_effective_properties
from solver import PRLVesicleSolver


def run_smoke_test(output_dir=None):
    if output_dir is None:
        output_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
    os.makedirs(output_dir, exist_ok=True)
    params = MembraneParameters()
    solver = PRLVesicleSolver(params=params, n_points=140, n_ctrl_Lo=6, n_ctrl_Ld=10)

    props_05 = get_effective_properties(params, theta=0.5)
    print("=" * 80)
    print("RUNNING PRL PYTHON SOLVER SMOKE TEST (Inner-Leaflet Ld Binding)")
    print("=" * 80)
    print(f"[Parameter Audit: theta = 0.5]")
    print(f"  - c0_Ld (inner leaflet): {props_05['c0_Ld']:+.4f} um^-1  (Expected > 0)")
    print(f"  - c0_Lo (bare):          {props_05['c0_Lo']:+.4f} um^-1  (Expected == 0)")
    print(f"  - kappa_Ld:              {props_05['kappa_Ld']:.1f} k_B T (Expected 55.0 k_B T)")
    print(f"  - kappa_Lo:              {props_05['kappa_Lo']:.1f} k_B T (Expected 80.0 k_B T)")
    print(f"  - delta_A0_shift:        {props_05['delta_A0_shift']:+.5f} um^2 (Expected > 0)")
    print(f"  - k_ade_eff:             {props_05['k_ade_eff']:.2f} k_B T")

    test_conditions = [
        {"theta": 0.0, "branch": "sphere"},
        {"theta": 0.0, "branch": "inward"},
        {"theta": 0.5, "branch": "sphere"},
        {"theta": 0.5, "branch": "inward"},
        {"theta": 0.5, "branch": "auto"},
    ]

    results = []
    for cond in test_conditions:
        theta = cond["theta"]
        branch = cond["branch"]
        print(f"\nSolving condition: theta = {theta:.2f}, branch = {branch} ...")

        res = solver.optimize_shape(theta=theta, branch=branch, maxiter=500)

        row = {
            "theta": theta,
            "requested_branch": branch,
            "selected_branch": res["branch"],
            "converged": res["converged"],
            "max_geom_error": res["max_geom_error"],
            "res_V": res["res_V"],
            "res_A_Lo": res["res_A_Lo"],
            "res_A_tot": res["res_A_tot"],
            "res_r_pole": res["res_r_pole"],
            "r_south_pole_um": float(res["r"][-1]),
            "min_r_interior_um": res["min_r_interior"],
            "s_split_um": res["s_split"],
            "S_tot_um": res["S_tot"],
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
            "dome_peak_to_dimple_um": res["dome_peak_to_dimple"],
            "peak_to_dimple_crest_um": res["peak_to_dimple_crest"],
            "inward_d_raw_um": res.get("inward_branch_result", res)["dimple_depth"],
            "inward_dome_peak_to_dimple_um": res.get("inward_branch_result", res)["dome_peak_to_dimple"],
            "moment_residual": res["moment_residual"],
            "iterations": res["iterations"]
        }
        results.append(row)
        print(
            f"  -> Converged: {res['converged']} | Total E: {res['total_energy']:.2f} kBT | "
            f"Max Geom Err: {res['max_geom_error']:.2e} | r(S_tot): {res['r'][-1]:+.2e} um | "
            f"Raw Depth: {res['dimple_depth']:+.3f} um | Peak-to-Dimple: {res['dome_peak_to_dimple']:+.3f} um | Moment Res: {res['moment_residual']:+.2f}"
        )

    df = pd.DataFrame(results)
    d0 = float(df[(df["theta"] == 0.0) & (df["requested_branch"] == "sphere")]["dimple_depth_um"].iloc[0])
    df["stroke_from_theta0_dome_um"] = df["dimple_depth_um"] - d0
    df["inward_stroke_from_theta0_dome_um"] = df["inward_d_raw_um"] - d0
    out_csv = os.path.join(output_dir, "smoke_test_prl.csv")
    df.to_csv(out_csv, index=False)
    print("=" * 80)
    print(f"Smoke test completed successfully. Output saved to {out_csv}")
    print("=" * 80)
    return df


if __name__ == "__main__":
    run_smoke_test()
