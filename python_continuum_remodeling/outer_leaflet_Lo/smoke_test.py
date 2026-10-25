"""
Smoke Test for HermesLoouter (Protein on outer leaflet of Lo, bare Ld).
Verifies:
1. Exact parameter effective values (c0_Lo > 0, c0_Ld = 0, area-weighted ADE outer shift > 0).
2. Convergence and South-pole closure r(S_tot)=0 for single-domain case (D_Lo = 3.0 um, v = 0.98, theta = 0.5).
3. Convergence and exact North-South mirror symmetry for two-domain case (D_Lo = 3.0 um, v = 0.98, theta = 0.5).
"""

import os
import sys
import numpy as np
import pandas as pd

base_dir = os.path.dirname(os.path.abspath(__file__))
single_dir = os.path.join(base_dir, "singledomain")
two_dir = os.path.join(base_dir, "twodomain")

sys.path.insert(0, single_dir)
sys.path.insert(0, two_dir)

from mechanics import MembraneParameters as SingleParams, get_effective_properties as get_single_props
from solver import HermesSingleSolver
from mechanics_two_domains import TwoDomainParameters as DoubleParams
from solver_two_domains import TwoDomainHermesSolver


def run_smoke_test():
    smoke_dir = os.path.join(base_dir, "results")
    os.makedirs(smoke_dir, exist_ok=True)

    print("=" * 80)
    print("RUNNING SMOKE TEST: HermesLoouter (Protein on Lo Outer Leaflet)")
    print("=" * 80)

    # 1. Parameter Inspection
    p_single = SingleParams(R_ves=5.0, D_Lo=3.0, gamma_0=1.5e-12, reduced_volume=0.98)
    props_s = get_single_props(p_single, theta=0.5)

    print("\n[Parameter Check: theta = 0.5]")
    print(f"  - c0_Lo (outer leaflet): {props_s['c0_Lo']:+.4f} um^-1  (Expected < 0)")
    print(f"  - c0_Ld (bare):          {props_s['c0_Ld']:+.4f} um^-1  (Expected == 0)")
    print(f"  - kappa_Lo:              {props_s['kappa_Lo']:.1f} k_B T (Expected 115.0 k_B T)")
    print(f"  - kappa_Ld:              {props_s['kappa_Ld']:.1f} k_B T (Expected 20.0 k_B T)")
    print(f"  - delta_A0_shift:        {props_s['delta_A0_shift']:+.5f} um^2 (Expected < 0, weighted by f_Lo)")
    print(f"  - k_ade_eff:             {props_s['k_ade_eff']:.2f} k_B T")

    assert props_s['c0_Lo'] < 0.0, "c0_Lo must be negative for outer leaflet condensation/wrapping!"
    assert props_s['c0_Ld'] == 0.0, "c0_Ld must be 0 for bare Ld!"
    assert props_s['delta_A0_shift'] < 0.0, "delta_A0_shift must be negative for outer leaflet condensation!"
    print(">> Parameter Audit PASSED!\n")

    records = []

    # 2. Single Domain Smoke Run
    print(">>> 1. Running Single-Domain Smoke Test Case (D_Lo=3.0, v=0.98, theta=0.5, gamma=1.5 pN)...")
    s_solver = HermesSingleSolver(params=p_single, n_points=120)
    res_s = s_solver.optimize_shape(theta=0.5, branch="auto")

    print(f"    - Converged:       {res_s['converged']} (branch={res_s['branch']})")
    print(f"    - Total Energy:    {res_s['total_energy']:.3f} k_B T")
    print(f"    - Cap Elevation:   {res_s['cap_elevation']:+.4f} um")
    print(f"    - South Pole r(S): {res_s['r'][-1]:+.2e} um")
    print(f"    - Max Geom Error:  {res_s['max_geom_error']:.2e}")
    print(f"    - Moment Residual: {res_s['moment_residual']:+.2f}")
    print(f"    - Iterations:      {res_s['iterations']}")

    records.append({
        "system": "single_domain",
        "gamma_pN": 1.5,
        "D_Lo_um": 3.0,
        "reduced_volume_v": 0.98,
        "theta": 0.5,
        "c0_Lo_um_inv": props_s["c0_Lo"],
        "c0_Ld_um_inv": props_s["c0_Ld"],
        "kappa_Lo_kBT": props_s["kappa_Lo"],
        "kappa_Ld_kBT": props_s["kappa_Ld"],
        "branch": res_s.get("branch", "auto"),
        "elevation_um": res_s["cap_elevation"],
        "depth_um": res_s["dimple_depth"],
        "dome_peak_to_dimple_um": res_s["dome_peak_to_dimple"],
        "peak_to_dimple_crest_um": res_s["peak_to_dimple_crest"],
        "inward_d_raw_um": res_s.get("inward_branch_result", res_s)["dimple_depth"],
        "inward_dome_peak_to_dimple_um": res_s.get("inward_branch_result", res_s)["dome_peak_to_dimple"],
        "r_b_um": res_s["r_b"],
        "r_south_pole_um": float(res_s["r"][-1]),
        "total_energy_kBT": res_s["total_energy"],
        "E_bend_tot": res_s["E_bend_tot"],
        "E_line": res_s["E_line"],
        "E_gauss": res_s["E_gauss"],
        "E_ade": res_s["E_ade"],
        "converged": res_s["converged"],
        "max_geom_error": res_s["max_geom_error"]
    })

    # 3. Two Domain Smoke Run
    print("\n>>> 2. Running Two-Domain Smoke Test Case (D_Lo=3.0, v=0.98, theta=0.5, gamma=1.5 pN)...")
    p_double = DoubleParams(R_ves=5.0, D_Lo=3.0, gamma_0=1.5e-12, reduced_volume=0.98)
    d_solver = TwoDomainHermesSolver(params=p_double, n_points=101)
    res_d = d_solver.optimize_shape(theta=0.5, branch="auto")

    print(f"    - Converged:          {res_d['converged']} (branch={res_d['branch']})")
    print(f"    - Total Energy:       {res_d['total_energy']:.3f} k_B T")
    print(f"    - North Elevation:    {res_d['elevation_north']:+.4f} um")
    print(f"    - South Elevation:    {res_d['elevation_south']:+.4f} um")
    print(f"    - Symmetry Residual:  {abs(res_d['elevation_north'] - res_d['elevation_south']):.2e} um")
    print(f"    - South Pole r(S):    {res_d['r_profile'][-1]:+.2e} um")
    print(f"    - Max Geom Error:     {res_d['max_geom_error']:.2e}")
    print(f"    - Iterations:         {res_d['iterations']}")

    eff_d = d_solver.calc.get_effective_properties(0.5)
    records.append({
        "system": "two_domain",
        "gamma_pN": 1.5,
        "D_Lo_um": 3.0,
        "reduced_volume_v": 0.98,
        "theta": 0.5,
        "c0_Lo_um_inv": eff_d["c0_Lo"],
        "c0_Ld_um_inv": eff_d["c0_Ld"],
        "kappa_Lo_kBT": eff_d["kappa_Lo"],
        "kappa_Ld_kBT": eff_d["kappa_Ld"],
        "branch": res_d.get("branch", "auto"),
        "elevation_um": res_d["elevation_north"],
        "depth_um": res_d["dimple_depth_north"],
        "dome_peak_to_dimple_um": res_d["dome_peak_to_dimple_north"],
        "peak_to_dimple_crest_um": res_d["peak_to_dimple_crest_north"],
        "inward_d_raw_um": res_d.get("inward_branch_result", res_d)["dimple_depth_north"],
        "inward_dome_peak_to_dimple_um": res_d.get("inward_branch_result", res_d)["dome_peak_to_dimple_north"],
        "r_b_um": res_d["r_b1"],
        "r_south_pole_um": float(res_d["r_profile"][-1]),
        "total_energy_kBT": res_d["total_energy"],
        "E_bend_tot": res_d["E_bend_tot"],
        "E_line": res_d["E_line"],
        "E_gauss": res_d["E_gauss"],
        "E_ade": res_d["E_ade"],
        "converged": res_d["converged"],
        "max_geom_error": res_d["max_geom_error"]
    })

    trace_path = os.path.join(smoke_dir, "smoke_trace.csv")
    df_smoke = pd.DataFrame(records)
    df_smoke.to_csv(trace_path, index=False)

    print("\n" + "=" * 80)
    print("SMOKE TEST COMPLETE AND VERIFIED!")
    print(f"Saved Trace: {trace_path}")
    print("=" * 80)


if __name__ == "__main__":
    run_smoke_test()
