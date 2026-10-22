"""
Phase Diagram and Morphological Landscape for Reduced Volume (v = V / V_sphere)
Comparing Single-Domain and Two-Domain Vesicles (D_Lo = 3.0 um, R_ves = 5.0 um, gamma = 1.5 pN).
Ultra-clean publication style: Only panel labels (a, b, c, d), minimal standard axes, no text boxes or explanatory phrases.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.interpolate import griddata

def generate_reduced_volume_phase_diagram(
    out_path="simulation_results/phase_diagram_reduced_volume.png"
):
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    
    # -------------------------------------------------------------------------
    # 1. Load Data
    # -------------------------------------------------------------------------
    f_single = "simulation_results_prl/param_sweeps/trace_deflation_sweep.csv"
    f_two = "simulation_results_prl/two_domains_sweeps/trace_deflation_sweep.csv"
    
    if not os.path.exists(f_single):
        f_single = "C:/Users/nares/Documents/Antigravity/Hermesall/simulation_results_prl/param_sweeps/trace_deflation_sweep.csv"
    if not os.path.exists(f_two):
        f_two = "C:/Users/nares/Documents/Antigravity/Hermesall/simulation_results_prl/two_domains_sweeps/trace_deflation_sweep.csv"

    df_s_all = pd.read_csv(f_single)
    df_d_all = pd.read_csv(f_two)
    
    df_s = df_s_all[np.isclose(df_s_all["D_Lo_um"], 3.0)].copy()
    df_d = df_d_all[np.isclose(df_d_all["D_Lo_um"], 3.0)].copy()

    # -------------------------------------------------------------------------
    # 2. Compute Depth Calculated from Peak Baseline (theta = 0.0)
    # -------------------------------------------------------------------------
    stroke_s_list = []
    for _, row in df_s.iterrows():
        v_val = row["reduced_volume_v"]
        sub = df_s[np.isclose(df_s["reduced_volume_v"], v_val)]
        d0 = sub[np.isclose(sub["theta"], 0.0)]["dimple_depth_um"].values[0]
        stroke_s_list.append(row["dimple_depth_um"] - d0)
    df_s["depth_from_peak"] = stroke_s_list

    stroke_d_list = []
    for _, row in df_d.iterrows():
        v_val = row["reduced_volume_v"]
        sub = df_d[np.isclose(df_d["reduced_volume_v"], v_val)]
        d0 = sub[np.isclose(sub["theta"], 0.0)]["depth_north_um"].values[0]
        stroke_d_list.append(row["depth_north_um"] - d0)
    df_d["depth_from_peak"] = stroke_d_list

    # -------------------------------------------------------------------------
    # 3. Clean Styling with Prominent Bold Fonts for A4 Publication
    # -------------------------------------------------------------------------
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "axes.labelweight": "bold",
        "axes.titleweight": "bold",
        "figure.titleweight": "bold",
        "axes.labelsize": 17,
        "axes.titlesize": 20,
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
        "legend.fontsize": 12.5,
        "mathtext.fontset": "cm",
        "axes.linewidth": 1.6,
        "xtick.major.width": 1.5,
        "ytick.major.width": 1.5,
        "xtick.major.size": 6.0,
        "ytick.major.size": 6.0,
        "figure.dpi": 300
    })

    fig = plt.figure(figsize=(12.2, 9.2))
    gs = GridSpec(2, 2, figure=fig, hspace=0.32, wspace=0.36)

    vi = np.linspace(0.895, 1.005, 150)
    thi = np.linspace(0.0, 1.0, 150)
    VI, THI = np.meshgrid(vi, thi)
    cmap = matplotlib.colormaps["Spectral"]
    xticks_common = np.arange(0.90, 1.01, 0.02)
    xticklabels_common = [f"{x:.2f}" for x in xticks_common]

    # =========================================================================
    # Panel (a): 2D Phase Diagram - Single Domain
    # =========================================================================
    ax_a = fig.add_subplot(gs[0, 0])
    
    v_s = df_s["reduced_volume_v"].values
    th_s = df_s["theta"].values
    depth_peak_s = df_s["depth_from_peak"].values
    
    grid_s = griddata((v_s, th_s), depth_peak_s, (VI, THI), method="cubic")
    grid_s_near = griddata((v_s, th_s), depth_peak_s, (VI, THI), method="nearest")
    grid_s = np.where(np.isnan(grid_s), grid_s_near, grid_s)
    
    cs1 = ax_a.contourf(VI, THI, grid_s, levels=30, cmap=cmap, alpha=0.85)
    cbar1 = plt.colorbar(cs1, ax=ax_a, pad=0.03)
    cbar1.set_label(r"$\Delta d\ (\mu\mathrm{m})$", fontsize=15.5, fontweight="bold", labelpad=8)
    cbar1.ax.tick_params(labelsize=12.5, width=1.5)
    plt.setp(cbar1.ax.get_yticklabels(), weight="bold")
    
    # Boundary contour: zero-crossing
    ax_a.contour(VI, THI, grid_s, levels=[-0.05], colors=["#1a1a1a"], linewidths=2.6, linestyles=["--"])
    ax_a.scatter(v_s, th_s, c=depth_peak_s, cmap=cmap, edgecolors="black", s=65, linewidth=1.2, zorder=5)
    
    ax_a.set_title("(a)", loc="left", fontweight="bold", fontsize=20, pad=9)
    ax_a.set_xlabel(r"$v$", fontsize=17, fontweight="bold")
    ax_a.set_ylabel(r"$\theta$", fontsize=17, fontweight="bold")
    ax_a.set_xlim(0.895, 1.005)
    ax_a.set_ylim(-0.02, 1.02)
    ax_a.set_xticks(xticks_common)
    ax_a.set_xticklabels(xticklabels_common)
    ax_a.grid(True, linestyle=":", alpha=0.55)

    # =========================================================================
    # Panel (b): 2D Phase Diagram - Two Domains
    # =========================================================================
    ax_b = fig.add_subplot(gs[0, 1])
    
    v_d = df_d["reduced_volume_v"].values
    th_d = df_d["theta"].values
    depth_peak_d = df_d["depth_from_peak"].values
    
    grid_d = griddata((v_d, th_d), depth_peak_d, (VI, THI), method="cubic")
    grid_d_near = griddata((v_d, th_d), depth_peak_d, (VI, THI), method="nearest")
    grid_d = np.where(np.isnan(grid_d), grid_d_near, grid_d)
    
    cs2 = ax_b.contourf(VI, THI, grid_d, levels=30, cmap=cmap, alpha=0.85)
    cbar2 = plt.colorbar(cs2, ax=ax_b, pad=0.03)
    cbar2.set_label(r"$\Delta d\ (\mu\mathrm{m})$", fontsize=15.5, fontweight="bold", labelpad=8)
    cbar2.ax.tick_params(labelsize=12.5, width=1.5)
    plt.setp(cbar2.ax.get_yticklabels(), weight="bold")
    
    ax_b.contour(VI, THI, grid_d, levels=[-0.05], colors=["#1a1a1a"], linewidths=2.6, linestyles=["--"])
    ax_b.scatter(v_d, th_d, c=depth_peak_d, cmap=cmap, edgecolors="black", s=60, linewidth=1.1, zorder=5)
    
    ax_b.set_title("(b)", loc="left", fontweight="bold", fontsize=20, pad=9)
    ax_b.set_xlabel(r"$v$", fontsize=17, fontweight="bold")
    ax_b.set_ylabel(r"$\theta$", fontsize=17, fontweight="bold")
    ax_b.set_xlim(0.895, 1.005)
    ax_b.set_ylim(-0.02, 1.02)
    ax_b.set_xticks(xticks_common)
    ax_b.set_xticklabels(xticklabels_common)
    ax_b.grid(True, linestyle=":", alpha=0.55)

    # =========================================================================
    # Panel (c): Depth vs v
    # =========================================================================
    ax_c = fig.add_subplot(gs[1, 0])
    
    s_th04 = df_s[np.isclose(df_s["theta"], 0.4)].sort_values("reduced_volume_v")
    s_th02 = df_s[np.isclose(df_s["theta"], 0.2)].sort_values("reduced_volume_v")
    d_th04 = df_d[np.isclose(df_d["theta"], 0.4)].sort_values("reduced_volume_v")
    d_th02 = df_d[np.isclose(df_d["theta"], 0.2)].sort_values("reduced_volume_v")

    ax_c.plot(s_th04["reduced_volume_v"], s_th04["depth_from_peak"], "o-", color="#d73027", linewidth=2.6, markersize=7,
              label=r"1-Domain ($\theta = 0.4$)")
    ax_c.plot(s_th02["reduced_volume_v"], s_th02["depth_from_peak"], "o--", color="#f46d43", linewidth=2.2, markersize=6,
              label=r"1-Domain ($\theta = 0.2$)")
              
    ax_c.plot(d_th04["reduced_volume_v"], d_th04["depth_from_peak"], "s-", color="#2171b5", linewidth=2.6, markersize=7,
              label=r"2-Domain ($\theta = 0.4$)")
    ax_c.plot(d_th02["reduced_volume_v"], d_th02["depth_from_peak"], "s--", color="#6baed6", linewidth=2.2, markersize=6,
              label=r"2-Domain ($\theta = 0.2$)")

    ax_c.plot(d_th04["reduced_volume_v"], 2.0 * d_th04["depth_from_peak"], "^:", color="#08306b", linewidth=2.4, markersize=7,
              label=r"2-Domain (Total)")

    ax_c.axhline(0, color="gray", linestyle="-", linewidth=1.0, alpha=0.7)

    ax_c.set_title("(c)", loc="left", fontweight="bold", fontsize=20, pad=9)
    ax_c.set_xlabel(r"$v$", fontsize=17, fontweight="bold")
    ax_c.set_ylabel(r"$\Delta d\ (\mu\mathrm{m})$", fontsize=17, fontweight="bold")
    ax_c.set_xlim(0.895, 1.005)
    ax_c.set_xticks(xticks_common)
    ax_c.set_xticklabels(xticklabels_common)
    ax_c.grid(True, linestyle=":", alpha=0.55)
    ax_c.legend(frameon=True, facecolor="white", edgecolor="#cccccc", loc="lower right",
                prop={"weight": "bold", "size": 11.5})

    # =========================================================================
    # Panel (d): Energy F(v)
    # =========================================================================
    ax_d = fig.add_subplot(gs[1, 1])
    
    s_th04_e = df_s[np.isclose(df_s["theta"], 0.4)].sort_values("reduced_volume_v")
    d_th04_e = df_d[np.isclose(df_d["theta"], 0.4)].sort_values("reduced_volume_v")
    s_th02_e = df_s[np.isclose(df_s["theta"], 0.2)].sort_values("reduced_volume_v")
    d_th02_e = df_d[np.isclose(df_d["theta"], 0.2)].sort_values("reduced_volume_v")

    ax_d.plot(s_th04_e["reduced_volume_v"], s_th04_e["total_energy_kBT"], "o-", color="#d73027", linewidth=2.6, markersize=7,
              label=r"1-Domain ($\theta = 0.4$)")
    ax_d.plot(s_th02_e["reduced_volume_v"], s_th02_e["total_energy_kBT"], "o--", color="#f46d43", linewidth=2.2, markersize=6,
              label=r"1-Domain ($\theta = 0.2$)")
              
    ax_d.plot(d_th04_e["reduced_volume_v"], d_th04_e["total_energy_kBT"], "s-", color="#2171b5", linewidth=2.6, markersize=7,
              label=r"2-Domain ($\theta = 0.4$)")
    ax_d.plot(d_th02_e["reduced_volume_v"], d_th02_e["total_energy_kBT"], "s--", color="#6baed6", linewidth=2.2, markersize=6,
              label=r"2-Domain ($\theta = 0.2$)")
              
    ax_d.set_title("(d)", loc="left", fontweight="bold", fontsize=20, pad=9)
    ax_d.set_xlabel(r"$v$", fontsize=17, fontweight="bold")
    ax_d.set_ylabel(r"$\mathcal{F}\ (k_\mathrm{B} T)$", fontsize=17, fontweight="bold")
    ax_d.set_xlim(0.895, 1.005)
    ax_d.set_ylim(4400, 10600)
    ax_d.set_xticks(xticks_common)
    ax_d.set_xticklabels(xticklabels_common)
    ax_d.grid(True, linestyle=":", alpha=0.55)
    ax_d.legend(frameon=True, facecolor="white", edgecolor="#cccccc", loc="upper right",
                prop={"weight": "bold", "size": 11.5})

    # Ensure bold ticks on all axes
    for ax in [ax_a, ax_b, ax_c, ax_d]:
        plt.setp(ax.get_xticklabels(), weight="bold")
        plt.setp(ax.get_yticklabels(), weight="bold")

    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Phase diagram successfully saved: {out_path}")

if __name__ == "__main__":
    generate_reduced_volume_phase_diagram()
