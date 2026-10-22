"""
Universal Markdown Benchmark Generator for Single-Domain vs. Two-Domain Systems (PRL Standard).
Dynamically parses trace_all_single_domain.csv and trace_two_domains_all.csv
and generates publication-quality Markdown comparative tables.
"""

import os
import sys
import numpy as np
import pandas as pd

def generate_markdown(artifact_path=r"C:\Users\nares\.gemini\antigravity\brain\847e2b35-470f-472b-a973-538b97b8e7f9\single_vs_double_domain_comprehensive_tables.md"):
    f_single = "simulation_results_prl/param_sweeps/trace_all_single_domain.csv"
    f_two = "simulation_results_prl/two_domains_sweeps/trace_two_domains_all.csv"
    
    if not os.path.exists(f_single) or not os.path.exists(f_two):
        print(f"Waiting for simulation traces: {f_single} and {f_two}")
        return
        
    df1 = pd.read_csv(f_single)
    df2 = pd.read_csv(f_two)
    
    lines = []
    lines.append("# Comprehensive Simulation Benchmark: Single-Domain vs. Two-Domain Systems")
    lines.append("")
    lines.append("> **Dataset Summary**: Complete audit of Single-Domain and Two-Domain production runs across the expanded parameter space:")
    lines.append("> - **Line Tension**: $\\gamma \\in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]\\,\\text{pN}$ (10 values)")
    lines.append("> - **Domain Size**: $D_{Lo} \\in [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]\\,\\mu\\text{m}$ (7 values)")
    lines.append("> - **Deflation Factor**: $v \\in [0.90, 0.95, 0.98, 1.00]$ (4 values)")
    lines.append("> - **Protein Coverage**: $\\theta \\in [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]$")
    lines.append("> ")
    lines.append("> All runs performed under exact NLP hard constraints ($|\\Delta V/V| < 10^{-12}$, $|\\Delta A/A| < 10^{-12}$) with canonical ADE energy and regularized South-pole curvature.")
    lines.append("")
    lines.append("---")
    lines.append("")
    
    # ---------------------------------------------------------
    # Helper to generate comparative table
    # ---------------------------------------------------------
    def format_table(df_s, df_t, title, param_name, param_unit):
        tbl = []
        tbl.append(f"### {title}")
        tbl.append("")
        tbl.append(f"| {param_name} ({param_unit}) | $\\theta$ | 1-Dom Depth ($\\mu\\mathrm{{m}}$) | 1-Dom Morph | 1-Dom $E_{{\\mathrm{{tot}}}}$ ($k_BT$) | 2-Dom Depth ($\\mu\\mathrm{{m}}$) | 2-Dom Morph | 2-Dom $E_{{\\mathrm{{tot}}}}$ ($k_BT$) | $\\Delta E$ ($k_BT$) |")
        tbl.append("| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |")
        
        merged = pd.merge(
            df_s, df_t,
            on=[col for col in ["gamma_pN", "D_Lo_um", "reduced_volume_v", "theta"] if col in df_s.columns and col in df_t.columns],
            suffixes=("_1", "_2")
        )
        
        for _, row in merged.iterrows():
            pval = row["gamma_pN"] if "gamma" in param_name.lower() else (row["D_Lo_um"] if "size" in param_name.lower() or "d_lo" in param_name.lower() else row["reduced_volume_v"])
            th = row["theta"]
            
            d1 = row["dimple_depth_um"]
            m1 = "Inward Dimple" if d1 < 0 else "Outward Dome"
            e1 = row["total_energy_kBT_1"]
            
            d2 = row["depth_north_um"]
            m2 = "Dual Inward" if d2 < 0 else "Dual Outward"
            e2 = row["total_energy_kBT_2"]
            
            de = e2 - e1
            tbl.append(f"| **{pval:.2f}** | {th:.2f} | `{d1:+.3f}` | {m1} | {e1:,.1f} | `{d2:+.3f}` | {m2} | {e2:,.1f} | `{de:+,.1f}` |")
            
        tbl.append("")
        return "\n".join(tbl)

    # 1. Line Tension Sweep
    df_g1 = df1[df1["sweep"] == "gamma"].sort_values(["gamma_pN", "theta"])
    df_g2 = df2[df2["sweep"] == "gamma"].sort_values(["gamma_pN", "theta"])
    lines.append("## 1. Line Tension ($\\gamma$) Sweep ($D_{Lo} = 3.0\\,\\mu\\text{m}, v = 0.98$)")
    lines.append(format_table(df_g1, df_g2, "Comparison Across Line Tensions 0.5 to 5.0 pN", "$\\gamma$", "pN"))
    lines.append("---")
    
    # 2. Domain Size Sweep
    df_d1 = df1[df1["sweep"] == "D_Lo"].sort_values(["D_Lo_um", "theta"])
    df_d2 = df2[df2["sweep"] == "D_Lo"].sort_values(["D_Lo_um", "theta"])
    lines.append("## 2. Domain Size ($D_{Lo}$) Sweep ($\\gamma = 1.5\\,\\text{pN}, v = 0.98$)")
    lines.append(format_table(df_d1, df_d2, "Comparison Across Domain Diameters 1.0 to 7.0 um", "$D_{Lo}$", "$\\mu$m"))
    lines.append("---")
    
    # 3. Deflation Sweep
    df_v1 = df1[df1["sweep"] == "deflation"].sort_values(["reduced_volume_v", "theta"])
    df_v2 = df2[df2["sweep"] == "deflation"].sort_values(["reduced_volume_v", "theta"])
    lines.append("## 3. Deflation Factor ($v$) Sweep ($\\gamma = 1.5\\,\\text{pN}, D_{Lo} = 3.0\\,\\mu\\text{m}$)")
    lines.append(format_table(df_v1, df_v2, "Comparison Across Reduced Volumes 0.90 to 1.00", "$v$", "rel"))
    lines.append("---")
    
    content = "\n".join(lines)
    os.makedirs(os.path.dirname(artifact_path), exist_ok=True)
    with open(artifact_path, "w", encoding="utf-8") as f:
        f.write(content)
        
    print(f"Comprehensive markdown table generated at: {artifact_path}")

if __name__ == "__main__":
    generate_markdown()
