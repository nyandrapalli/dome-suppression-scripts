"""
Comparison Plotting Script: Single Domain vs Two Domains for HermesLoouter.
Plots direct morphological and energetic comparisons across:
- Protein coverage theta
- Domain size D_Lo
- Deflation factor v
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def plot_comparison(root_dir="HermesLoouter"):
    single_dir = os.path.join(root_dir, "single")
    double_dir = os.path.join(root_dir, "double")
    out_fig = os.path.join(root_dir, "fig_comparison_single_vs_two_domains.png")
    
    single_trace_file = os.path.join(single_dir, "trace_domain_size_sweep.csv")
    double_trace_file = os.path.join(double_dir, "trace_domain_size_sweep.csv")
    
    if not (os.path.exists(single_trace_file) and os.path.exists(double_trace_file)):
        print("Traces not yet generated. Run sweeps first.")
        return
        
    df_s = pd.read_csv(single_trace_file)
    df_d = pd.read_csv(double_trace_file)
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), dpi=300)
    
    # 1. Total Energy vs Domain Size at theta = 1.0
    df_s_th1 = df_s[df_s["theta"] == 1.0].sort_values("D_Lo_um")
    df_d_th1 = df_d[df_d["theta"] == 1.0].sort_values("D_Lo_um")
    
    axes[0].plot(df_s_th1["D_Lo_um"], df_s_th1["total_energy_kBT"], "o-", color="#d73027", linewidth=2.2, label=r"Single $L_o$ Domain (Outer Prot)")
    axes[0].plot(df_d_th1["D_Lo_um"], df_d_th1["total_energy_kBT"], "s--", color="#4575b4", linewidth=2.2, label=r"Two $L_o$ Domains (Outer Prot)")
    axes[0].set_xlabel(r"Domain Diameter $D_{Lo}$ ($\mu\mathrm{m}$)", fontsize=11)
    axes[0].set_ylabel(r"Total Energy $E_{tot}$ ($k_B T$)", fontsize=11)
    axes[0].set_title(r"Total Free Energy vs Domain Size ($\theta = 1.0$, $v = 0.98$)", fontsize=12, fontweight="bold")
    axes[0].grid(True, linestyle=":", alpha=0.6)
    axes[0].legend(fontsize=10)
    
    # 2. Elevation / Budding vs Domain Size
    axes[1].plot(df_s_th1["D_Lo_um"], df_s_th1["cap_elevation_um"], "o-", color="#d73027", linewidth=2.2, label="Single Cap Elevation")
    axes[1].plot(df_d_th1["D_Lo_um"], df_d_th1["elevation_north_um"], "s--", color="#4575b4", linewidth=2.2, label="Two Caps North Elevation")
    axes[1].set_xlabel(r"Domain Diameter $D_{Lo}$ ($\mu\mathrm{m}$)", fontsize=11)
    axes[1].set_ylabel(r"Cap Elevation above Boundary ($\mu\mathrm{m}$)", fontsize=11)
    axes[1].set_title(r"Cap Outward Protrusion vs Domain Size ($\theta = 1.0$)", fontsize=12, fontweight="bold")
    axes[1].grid(True, linestyle=":", alpha=0.6)
    axes[1].legend(fontsize=10)
    
    plt.tight_layout()
    fig.savefig(out_fig)
    plt.close(fig)
    print(f"Saved comparison figure: {out_fig}")

if __name__ == "__main__":
    plot_comparison()
