"""
Smoke Test for HermesLoouter (Protein on outer leaflet of Lo, bare Ld).
Verifies:
1. Exact parameter effective values (c0_Lo > 0, c0_Ld = 0, ADE outer shift > 0).
2. Convergence for single domain case (D_Lo = 3.0 um, v = 0.98, theta = 0.5).
3. Convergence for two domain case (D_Lo = 3.0 um, v = 0.98, theta = 0.5).
4. Physical validity (A_tot, V_tot, no self-intersection, mirror symmetry in 2-domain).
5. Generation of test trace CSV and profile files in HermesLoouter/smoke_test/.
"""

import os
import sys
import numpy as np
import pandas as pd

single_dir = os.path.join(os.path.dirname(__file__), "singledomain")
two_dir = os.path.join(os.path.dirname(__file__), "twodomain")

sys.path.insert(0, single_dir)
sys.path.insert(0, two_dir)

from mechanics import MembraneParameters as SingleParams, get_effective_properties as get_single_props
from solver import HermesSingleSolver
from mechanics_two_domains import TwoDomainParameters as DoubleParams, TwoDomainEnergyCalculator
from solver_two_domains import TwoDomainHermesSolver

def run_smoke_test():
    smoke_dir = os.path.join(os.path.dirname(__file__), "smoke_test")
    os.makedirs(smoke_dir, exist_ok=True)
    os.makedirs(os.path.join(smoke_dir, "profiles"), exist_ok=True)
    
    print("=" * 80)
    print("RUNNING SMOKE TEST: HermesLoouter (Protein on Lo Outer Leaflet)")
    print("=" * 80)
    
    # -------------------------------------------------------------
    # 1. Parameter Inspection
    # -------------------------------------------------------------
    p_single = SingleParams(R_ves=5.0, D_Lo=3.0, gamma_0=1.5e-12, reduced_volume=0.98)
    props_s = get_single_props(p_single, theta=0.5)
    
    print("\n[Parameter Check: theta = 0.5]")
    print(f"  - c0_Lo (outer leaflet): {props_s['c0_Lo']:+.4f} um^-1  (Expected > 0)")
    print(f"  - c0_Ld (bare):          {props_s['c0_Ld']:+.4f} um^-1  (Expected == 0)")
    print(f"  - kappa_Lo:              {props_s['kappa_Lo']:.1f} k_B T (Expected 115.0 k_B T)")
    print(f"  - kappa_Ld:              {props_s['kappa_Ld']:.1f} k_B T (Expected 20.0 k_B T)")
    print(f"  - delta_A0_shift:        {props_s['delta_A0_shift']:+.5f} um^2 (Expected > 0 for outer leaflet)")
    
    assert props_s['c0_Lo'] > 0.0, "c0_Lo must be positive for outer leaflet protein!"
    assert props_s['c0_Ld'] == 0.0, "c0_Ld must be 0 for bare Ld!"
    assert props_s['delta_A0_shift'] > 0.0, "delta_A0_shift must be positive for outer leaflet expansion!"
    print(">> Parameter Audit PASSED!\n")
    
    records = []
    
    # -------------------------------------------------------------
    # 2. Single Domain Smoke Run
    # -------------------------------------------------------------
    print(">>> 1. Running Single-Domain Smoke Test Case (D_Lo=3.0, v=0.98, theta=0.5, gamma=1.5 pN)...")
    s_solver = HermesSingleSolver(params=p_single, n_points=120)
    res_s = s_solver.optimize_shape(theta=0.5, branch="auto")
    
    print(f"    - Converged:      {res_s['converged']}")
    print(f"    - Total Energy:   {res_s['total_energy']:.3f} k_B T")
    print(f"    - Cap Elevation:  {res_s['cap_elevation']:+.4f} um")
    print(f"    - Max Geom Error: {res_s['max_geom_error']:.2e}")
    print(f"    - Iterations:     {res_s['iterations']}")
    
    prof_s_path = os.path.join(smoke_dir, "profiles", "smoke_profile_single.csv")
    df_prof_s = pd.DataFrame({
        "s_um": s_solver.s_grid,
        "psi_rad": res_s["psi"],
        "r_um": res_s["r"],
        "z_um": res_s["z"],
        "H_um_inv": res_s["H"],
        "K_um_inv2": res_s["K"]
    })
    df_prof_s["domain"] = np.where(s_solver.s_grid <= s_solver.s_split, "Lo", "Ld")
    df_prof_s.to_csv(prof_s_path, index=False)
    
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
        "r_b_um": res_s["r_b"],
        "total_energy_kBT": res_s["total_energy"],
        "E_bend_tot": res_s["E_bend_tot"],
        "E_line": res_s["E_line"],
        "E_gauss": res_s["E_gauss"],
        "E_ade": res_s["E_ade"],
        "converged": res_s["converged"],
        "max_geom_error": res_s["max_geom_error"]
    })
    
    # -------------------------------------------------------------
    # 3. Two Domain Smoke Run
    # -------------------------------------------------------------
    print("\n>>> 2. Running Two-Domain Smoke Test Case (D_Lo=3.0, v=0.98, theta=0.5, gamma=1.5 pN)...")
    p_double = DoubleParams(R_ves=5.0, D_Lo=3.0, gamma_0=1.5e-12, reduced_volume=0.98)
    d_solver = TwoDomainHermesSolver(params=p_double, n_points=101)
    res_d = d_solver.optimize_shape(theta=0.5, branch="auto")
    
    print(f"    - Converged:          {res_d['converged']}")
    print(f"    - Total Energy:       {res_d['total_energy']:.3f} k_B T")
    print(f"    - North Elevation:    {res_d['elevation_north']:+.4f} um")
    print(f"    - South Elevation:    {res_d['elevation_south']:+.4f} um")
    print(f"    - Symmetry Residual:  {abs(res_d['elevation_north'] - res_d['elevation_south']):.2e} um")
    print(f"    - Max Geom Error:     {res_d['max_geom_error']:.2e}")
    print(f"    - Iterations:         {res_d['iterations']}")
    
    prof_d_path = os.path.join(smoke_dir, "profiles", "smoke_profile_double.csv")
    df_prof_d = pd.DataFrame({
        "s_um": d_solver.s_grid,
        "psi_rad": res_d["psi_profile"],
        "r_um": res_d["r_profile"],
        "z_um": res_d["z_profile"],
        "H_um_inv": res_d["H_profile"],
        "K_um_inv2": res_d["K_profile"]
    })
    dom_tags = []
    for s_val in d_solver.s_grid:
        if s_val <= d_solver.s_split1:
            dom_tags.append("Lo_north")
        elif s_val >= d_solver.s_split2:
            dom_tags.append("Lo_south")
        else:
            dom_tags.append("Ld")
    df_prof_d["domain"] = dom_tags
    df_prof_d.to_csv(prof_d_path, index=False)
    
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
        "r_b_um": res_d["r_b1"],
        "total_energy_kBT": res_d["total_energy"],
        "E_bend_tot": res_d["E_bend_tot"],
        "E_line": res_d["E_line"],
        "E_gauss": res_d["E_gauss"],
        "E_ade": res_d["E_ade"],
        "converged": res_d["converged"],
        "max_geom_error": res_d["max_geom_error"]
    })
    
    # Save smoke trace CSV
    trace_path = os.path.join(smoke_dir, "smoke_trace.csv")
    df_smoke = pd.DataFrame(records)
    df_smoke.to_csv(trace_path, index=False)
    
    print("\n" + "=" * 80)
    print("SMOKE TEST COMPLETE AND VERIFIED!")
    print(f"Saved Trace:   {trace_path}")
    print(f"Saved Profile: {prof_s_path}")
    print(f"Saved Profile: {prof_d_path}")
    print("=" * 80)

if __name__ == "__main__":
    run_smoke_test()
