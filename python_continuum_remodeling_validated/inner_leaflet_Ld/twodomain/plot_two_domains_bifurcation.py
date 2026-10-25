"""
Bifurcation and Energy Landscape Plot for Two-Domain System (PRL Standard).
Plots:
1. Energy vs. Coverage theta for v = 0.98 and v = 0.90 across branches.
2. Dimple Depth vs. Coverage theta highlighting the first-order bifurcation at theta >= 0.80.
3. Component Energy Breakdown (E_bend, E_line, E_ade).
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

def plot_bifurcation(out_dir="simulation_results_prl/two_domains_sweeps"):
    f_v98 = os.path.join(out_dir, "trace_two_domains_titration.csv")
    f_v90 = os.path.join(out_dir, "trace_two_domains_titration_v090.csv")
    out_fig = os.path.join(out_dir, "fig_two_domains_bifurcation.png")
    
    df_v98 = pd.read_csv(f_v98)
    df_v90 = pd.read_csv(f_v90)
    
    plt.rcParams.update({
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.titlesize": 12,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "legend.fontsize": 9,
        "figure.dpi": 300
    })
    
    fig = plt.figure(figsize=(16, 5.5))
    gs = GridSpec(1, 3, figure=fig, wspace=0.28)
    
    # -------------------------------------------------------------------------
    # Panel (a): Free Energy vs Coverage theta
    # -------------------------------------------------------------------------
    ax1 = fig.add_subplot(gs[0, 0])
    
    # v = 0.98 branches
    in_98 = df_v98[df_v98["branch"] == "inward"].sort_values("theta")
    out_98 = df_v98[df_v98["branch"] == "outward"].sort_values("theta")
    sph_98 = df_v98[df_v98["branch"] == "sphere"].sort_values("theta")
    
    # v = 0.90 branches
    in_90 = df_v90[df_v90["branch"] == "inward"].sort_values("theta")
    sph_90 = df_v90[df_v90["branch"] == "sphere"].sort_values("theta")
    
    ax1.plot(sph_98["theta"], sph_98["total_energy_kBT"], "k--", linewidth=1.5, label=r"Bare Reference ($v=0.98$)")
    ax1.plot(out_98["theta"], out_98["total_energy_kBT"], "o-", color="#2171b5", linewidth=2.2, markersize=5, label=r"$v=0.98$ Outward Spindle")
    ax1.plot(in_90["theta"], in_90["total_energy_kBT"], "s-", color="#d73027", linewidth=2.5, markersize=5.5, label=r"$v=0.90$ Dual Inward Dimple")
    
    ax1.set_title(r"$\mathbf{(a)}$ Bifurcation Energetics $\mathcal{F}(\theta)$", fontweight="bold", pad=8)
    ax1.set_xlabel(r"Streptavidin Coverage $\theta$")
    ax1.set_ylabel(r"Total Free Energy $\mathcal{F}$ ($k_B T$)")
    ax1.grid(True, linestyle=":", alpha=0.5)
    ax1.legend(frameon=True, facecolor="white", loc="lower left")

    # -------------------------------------------------------------------------
    # Panel (b): Cap Depth vs Coverage theta (Bifurcation Diagram)
    # -------------------------------------------------------------------------
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.axhline(0, color="gray", linestyle="-", linewidth=0.8, alpha=0.7)
    
    ax2.plot(out_98["theta"], out_98["depth_north_um"], "o-", color="#2171b5", linewidth=2.2, markersize=5, label=r"$v=0.98$ (Outward Domes)")
    ax2.plot(in_90["theta"], in_90["depth_north_um"], "s-", color="#d73027", linewidth=2.5, markersize=5.5, label=r"$v=0.90$ (Dual Inward Dimples)")
    
    min_row = in_90.loc[in_90["depth_north_um"].idxmin()]
    min_depth = min_row["depth_north_um"]
    min_theta = min_row["theta"]
    
    inward_rows = in_90[in_90["depth_north_um"] < 0]
    if len(inward_rows) > 0:
        th_snap = inward_rows["theta"].min()
        annot_text = f"Inward Transition\nat $\\theta \\geq {th_snap:.2f}$\n$\\mathbf{{{min_depth:+.2f}\\,\\mu\\mathrm{{m}}}}$ (Inward)"
    else:
        annot_text = f"Minimum Depth:\n$\\mathbf{{{min_depth:+.2f}\\,\\mu\\mathrm{{m}}}}$"
        
    ax2.annotate(
        annot_text,
        xy=(min_theta, min_depth), xytext=(0.45, min_depth + 0.15 if min_depth < 0 else min_depth - 0.2),
        arrowprops=dict(facecolor="#d73027", arrowstyle="->", lw=1.5),
        fontsize=9, fontweight="bold", color="#7f0000",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="#fff5f0", edgecolor="#d73027", alpha=0.95)
    )
    
    ax2.set_title(r"$\mathbf{(b)}$ Cap Depth Bifurcation $d(\theta)$", fontweight="bold", pad=8)
    ax2.set_xlabel(r"Streptavidin Coverage $\theta$")
    ax2.set_ylabel(r"Cap Depth $d$ ($\mu\mathrm{m}$) [$-$ Inward / $+$ Outward]")
    ax2.grid(True, linestyle=":", alpha=0.5)
    ax2.legend(frameon=True, facecolor="white", loc="upper left")

    # -------------------------------------------------------------------------
    # Panel (c): Energy Breakdown Components at v = 0.90
    # -------------------------------------------------------------------------
    ax3 = fig.add_subplot(gs[0, 2])
    
    ax3.plot(in_90["theta"], in_90["E_bend_tot"], "d-", color="#2ca02c", linewidth=2.0, markersize=5, label=r"Bending $E_{\mathrm{bend}}$")
    ax3.plot(in_90["theta"], in_90["E_line"], "^-", color="#ff7f0e", linewidth=2.0, markersize=5, label=r"Line Tension $E_{\mathrm{line}}$")
    ax3.plot(in_90["theta"], in_90["E_ade"], "x-", color="#9467bd", linewidth=2.0, markersize=5, label=r"Non-local ADE $E_{\mathrm{ADE}}$")
    
    ax3.set_title(r"$\mathbf{(c)}$ Energy Partition ($v=0.90$)", fontweight="bold", pad=8)
    ax3.set_xlabel(r"Streptavidin Coverage $\theta$")
    ax3.set_ylabel(r"Component Energy ($k_B T$)")
    ax3.grid(True, linestyle=":", alpha=0.5)
    ax3.legend(frameon=True, facecolor="white", loc="center right")

    plt.tight_layout()
    fig.savefig(out_fig, dpi=300)
    plt.close(fig)
    print(f"Bifurcation figure saved to: {out_fig}")
    return out_fig

if __name__ == "__main__":
    plot_bifurcation()
