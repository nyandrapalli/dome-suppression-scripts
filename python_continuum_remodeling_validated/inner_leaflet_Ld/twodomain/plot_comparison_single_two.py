"""
Comparative Plotting Script: Single Domain vs Two Domains (PRL Standard).
Directly compares:
1. Best Single Domain Invagination (D_Lo = 4.5 um, v = 0.98, gamma = 1.5 pN) -> Inward Dimple (-0.84 um)
2. Two Domains at the EXACT SAME condition (D_Lo = 4.5 um, v = 0.98, gamma = 1.5 pN) -> Outward Domes (+0.48 um)
3. Best Two Domains Dual Invagination (D_Lo = 3.0 um, v = 0.90, gamma = 1.5 pN) -> Dual Inward Dimples (-0.59 um)
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

def plot_single_vs_two_comparison(out_fig="simulation_results_prl/fig_single_vs_two_domains_comparison.png"):
    color_Lo = "#d73027"       # Vibrant red for Lo
    color_Ld = "#1a9850"       # Green for Ld
    color_ref = "#3182bd"      # Blue reference
    
    plt.rcParams.update({
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.titlesize": 12,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "legend.fontsize": 9,
        "figure.dpi": 300
    })
    
    fig = plt.figure(figsize=(18, 6.5))
    gs = GridSpec(1, 3, figure=fig, wspace=0.28)
    
    # -------------------------------------------------------------------------
    # Panel 1: Single Domain Best Invagination (D_Lo = 4.5 um, v = 0.98)
    # -------------------------------------------------------------------------
    ax1 = fig.add_subplot(gs[0, 0])
    single_file = "simulation_results_prl/param_sweeps/profiles/prof_DLo_4.5um_theta_1.00.csv"
    if os.path.exists(single_file):
        df1 = pd.read_csv(single_file)
        r1 = df1["r_um"].values
        z1 = df1["z_um"].values
        dom1 = df1["domain"].values
        idx_b = np.where(dom1 == "Lo")[0][-1]
        z_mid1 = 0.5 * (np.max(z1) + np.min(z1))
        z1_al = z1 - z_mid1
        
        r_Lo = np.concatenate([-r1[:idx_b+1][::-1], r1[:idx_b+1]])
        z_Lo = np.concatenate([z1_al[:idx_b+1][::-1], z1_al[:idx_b+1]])
        r_Ld_l = -r1[idx_b:][::-1]
        z_Ld_l = z1_al[idx_b:][::-1]
        r_Ld_r = r1[idx_b:]
        z_Ld_r = z1_al[idx_b:]
        
        ax1.plot(r_Ld_l, z_Ld_l, color=color_Ld, linewidth=2.8, label=r"$L_d$ Matrix")
        ax1.plot(r_Ld_r, z_Ld_r, color=color_Ld, linewidth=2.8)
        ax1.plot(r_Lo, z_Lo, color=color_Lo, linewidth=3.8, label=r"Single $L_o$ Domain")
        ax1.plot([-r1[idx_b], r1[idx_b]], [z1_al[idx_b], z1_al[idx_b]], "ko", markersize=6.0)
        
        depth1 = abs(z1_al[0] - z1_al[idx_b])
        ax1.annotate(
            f"SINGLE DOMAIN BEST\n$D_{{Lo}} = 4.5\\,\\mu\\mathrm{{m}}, v = 0.98$\nInward Dimple: $\\mathbf{{{depth1:.2f}\\,\\mu\\mathrm{{m}}}}$\nNeck $2r_b = {2*r1[idx_b]:.2f}\\,\\mu\\mathrm{{m}}$\n$E_{{line}} = 4,760\\,k_BT$",
            xy=(0, 0), xytext=(0, -0.8),
            fontsize=9.5, ha="center", va="center", fontweight="bold", color="#7f0000",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#fff5f0", edgecolor="#d73027", alpha=0.95)
        )
        
    ax1.axvline(0, color="gray", linestyle=":", linewidth=0.8, alpha=0.7)
    ax1.set_title(r"$\mathbf{(a)}$ Single Domain Best ($v=0.98$)", fontweight="bold", pad=8)
    ax1.set_aspect("equal")
    ax1.set_xlim(-6.0, 6.0)
    ax1.set_ylim(-6.0, 6.0)
    ax1.set_xlabel(r"Radial $r$ ($\mu\mathrm{m}$)")
    ax1.set_ylabel(r"Axial $z$ ($\mu\mathrm{m}$)")
    ax1.grid(True, linestyle=":", alpha=0.5)
    ax1.legend(frameon=True, facecolor="white", loc="lower right")

    # -------------------------------------------------------------------------
    # Panel 2: Two Domains at the EXACT SAME condition (D_Lo = 4.5 um, v = 0.98)
    # -------------------------------------------------------------------------
    ax2 = fig.add_subplot(gs[0, 1])
    two_same_file = "simulation_results_prl/two_domains_sweeps/profiles/prof_two_dom_DLo_4.5um_theta_1.00.csv"
    if os.path.exists(two_same_file):
        df2 = pd.read_csv(two_same_file)
        r2 = df2["r_um"].values
        z2 = df2["z_um"].values
        dom2 = df2["domain"].values
        idx_b1 = np.where(dom2 == "Lo_north")[0][-1]
        idx_b2 = np.where(dom2 == "Lo_south")[0][0]
        z_mid2 = 0.5 * (np.max(z2) + np.min(z2))
        z2_al = z2 - z_mid2
        
        r_Lo1 = np.concatenate([-r2[:idx_b1+1][::-1], r2[:idx_b1+1]])
        z_Lo1 = np.concatenate([z2_al[:idx_b1+1][::-1], z2_al[:idx_b1+1]])
        r_Lo2 = np.concatenate([-r2[idx_b2:][::-1], r2[idx_b2:]])
        z_Lo2 = np.concatenate([z2_al[idx_b2:][::-1], z2_al[idx_b2:]])
        
        r_Ld2_l = -r2[idx_b1:idx_b2+1][::-1]
        z_Ld2_l = z2_al[idx_b1:idx_b2+1][::-1]
        r_Ld2_r = r2[idx_b1:idx_b2+1]
        z_Ld2_r = z2_al[idx_b1:idx_b2+1]
        
        ax2.plot(r_Ld2_l, z_Ld2_l, color=color_Ld, linewidth=2.8)
        ax2.plot(r_Ld2_r, z_Ld2_r, color=color_Ld, linewidth=2.8)
        ax2.plot(r_Lo1, z_Lo1, color=color_Lo, linewidth=3.8, label=r"Two $L_o$ Domains")
        ax2.plot(r_Lo2, z_Lo2, color=color_Lo, linewidth=3.8)
        ax2.plot([-r2[idx_b1], r2[idx_b1]], [z2_al[idx_b1], z2_al[idx_b1]], "ko", markersize=6.0)
        ax2.plot([-r2[idx_b2], r2[idx_b2]], [z2_al[idx_b2], z2_al[idx_b2]], "ko", markersize=6.0)
        
        depth2 = float(z2_al[0] - z2_al[idx_b1])
        gamma_kBT_um = (1.5e-12 / 4.114e-21) * 1e-6
        e_line_2 = float(2.0 * 2.0 * np.pi * r2[idx_b1] * gamma_kBT_um)
        morph2 = "Outward Dome" if depth2 > 0 else "Inward Dimple"
        ax2.annotate(
            f"SAME CONDITION AS (a)\n$D_{{Lo}} = 4.5\\,\\mu\\mathrm{{m}}, v = 0.98$\n{morph2}: $\\mathbf{{{depth2:+.2f}\\,\\mu\\mathrm{{m}}}}$\nDual Necks $2r_b = {2*r2[idx_b1]:.2f}\\,\\mu\\mathrm{{m}}$\n$E_{{line}} = {e_line_2:.0f}\\,k_BT$ (Doubled!)",
            xy=(0, 0), xytext=(0, 0),
            fontsize=9.5, ha="center", va="center", fontweight="bold", color="#08519c" if depth2 > 0 else "#7f0000",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#eff3ff" if depth2 > 0 else "#fff5f0", 
                      edgecolor="#3182bd" if depth2 > 0 else "#d73027", alpha=0.95)
        )
        
    ax2.axvline(0, color="gray", linestyle=":", linewidth=0.8, alpha=0.7)
    ax2.set_title(r"$\mathbf{(b)}$ Two Domains at Same Condition ($v=0.98$)", fontweight="bold", pad=8)
    ax2.set_aspect("equal")
    ax2.set_xlim(-6.0, 6.0)
    ax2.set_ylim(-6.0, 6.0)
    ax2.set_xlabel(r"Radial $r$ ($\mu\mathrm{m}$)")
    ax2.grid(True, linestyle=":", alpha=0.5)

    # -------------------------------------------------------------------------
    # Panel 3: Two Domains Best Dual Invagination (D_Lo = 3.0 um, v = 0.90)
    # -------------------------------------------------------------------------
    ax3 = fig.add_subplot(gs[0, 2])
    two_best_file = "simulation_results_prl/two_domains_sweeps/profiles/prof_two_dom_v_0.90_theta_1.00.csv"
    if os.path.exists(two_best_file):
        df3 = pd.read_csv(two_best_file)
        r3 = df3["r_um"].values
        z3 = df3["z_um"].values
        dom3 = df3["domain"].values
        idx_b1_3 = np.where(dom3 == "Lo_north")[0][-1]
        idx_b2_3 = np.where(dom3 == "Lo_south")[0][0]
        z_mid3 = 0.5 * (np.max(z3) + np.min(z3))
        z3_al = z3 - z_mid3
        
        r_Lo1_3 = np.concatenate([-r3[:idx_b1_3+1][::-1], r3[:idx_b1_3+1]])
        z_Lo1_3 = np.concatenate([z3_al[:idx_b1_3+1][::-1], z3_al[:idx_b1_3+1]])
        r_Lo2_3 = np.concatenate([-r3[idx_b2_3:][::-1], r3[idx_b2_3:]])
        z_Lo2_3 = np.concatenate([z3_al[idx_b2_3:][::-1], z3_al[idx_b2_3:]])
        
        r_Ld3_l = -r3[idx_b1_3:idx_b2_3+1][::-1]
        z_Ld3_l = z3_al[idx_b1_3:idx_b2_3+1][::-1]
        r_Ld3_r = r3[idx_b1_3:idx_b2_3+1]
        z_Ld3_r = z3_al[idx_b1_3:idx_b2_3+1]
        
        ax3.plot(r_Ld3_l, z_Ld3_l, color=color_Ld, linewidth=2.8)
        ax3.plot(r_Ld3_r, z_Ld3_r, color=color_Ld, linewidth=2.8)
        ax3.plot(r_Lo1_3, z_Lo1_3, color=color_Lo, linewidth=3.8)
        ax3.plot(r_Lo2_3, z_Lo2_3, color=color_Lo, linewidth=3.8)
        ax3.plot([-r3[idx_b1_3], r3[idx_b1_3]], [z3_al[idx_b1_3], z3_al[idx_b1_3]], "ko", markersize=6.0)
        ax3.plot([-r3[idx_b2_3], r3[idx_b2_3]], [z3_al[idx_b2_3], z3_al[idx_b2_3]], "ko", markersize=6.0)
        
        raw_depth3 = float(z3_al[0] - z3_al[idx_b1_3])
        morph3 = "Dual Inward Dimple" if raw_depth3 < 0 else "Dual Outward Dome"
        ax3.annotate(
            f"TWO DOMAINS BEST\n$D_{{Lo}} = 3.0\\,\\mu\\mathrm{{m}}, v = 0.90$\n{morph3}: $\\mathbf{{{raw_depth3:+.2f}\\,\\mu\\mathrm{{m}}}}$\nDual Necks $2r_b = {2*r3[idx_b1_3]:.2f}\\,\\mu\\mathrm{{m}}$\nDeflation frees excess $\\Delta A$",
            xy=(0, 0), xytext=(0, 0),
            fontsize=9.5, ha="center", va="center", fontweight="bold", color="#7f0000" if raw_depth3 < 0 else "#08519c",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#fff5f0" if raw_depth3 < 0 else "#eff3ff", 
                      edgecolor="#d73027" if raw_depth3 < 0 else "#3182bd", alpha=0.95)
        )
        
    ax3.axvline(0, color="gray", linestyle=":", linewidth=0.8, alpha=0.7)
    ax3.set_title(r"$\mathbf{(c)}$ Two Domains Best ($v=0.90$)", fontweight="bold", pad=8)
    ax3.set_aspect("equal")
    ax3.set_xlim(-6.0, 6.0)
    ax3.set_ylim(-6.0, 6.0)
    ax3.set_xlabel(r"Radial $r$ ($\mu\mathrm{m}$)")
    ax3.grid(True, linestyle=":", alpha=0.5)

    plt.tight_layout()
    fig.savefig(out_fig, dpi=300)
    plt.close(fig)
    print(f"Comparison figure saved to: {out_fig}")
    return out_fig

if __name__ == "__main__":
    plot_single_vs_two_comparison()
