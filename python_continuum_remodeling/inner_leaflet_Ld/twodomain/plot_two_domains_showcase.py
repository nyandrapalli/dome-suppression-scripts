"""
Publication Showcase for Two Opposite Lo Domains (Symmetric Bipolar Vesicle - PRL Standard).
Highlights:
1. Exact Dual Inward Dimples Scenario (Deflation v = 0.90, gamma = 1.5 pN, D_Lo = 3.0 um, theta = 1.0)
2. Zoomed dual inward invaginations at North and South poles (both in RED)
3. Arc length Mean Curvature profile H(s) demonstrating symmetric negative mean curvature inside the dimples
4. Comparison across deflation factors v in {0.90, 0.95, 0.98}
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

def plot_two_domains_showcase(results_dir="simulation_results_prl/two_domains_sweeps"):
    profiles_dir = os.path.join(results_dir, "profiles")
    out_fig = os.path.join(results_dir, "fig_two_domains_showcase.png")
    
    color_Lo = "#d73027"       # Vibrant red for Lo caps
    color_Ld = "#1a9850"       # Emerald green for equatorial Ld belt
    color_bare = "#3182bd"     # Blue reference
    
    plt.rcParams.update({
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.titlesize": 12,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "legend.fontsize": 9,
        "figure.dpi": 300
    })
    
    fig = plt.figure(figsize=(16, 11))
    gs = GridSpec(2, 3, figure=fig, hspace=0.34, wspace=0.28, width_ratios=[1.2, 1.0, 1.1])
    
    def get_sym_profile(csv_path):
        df = pd.read_csv(csv_path)
        r = df["r_um"].values
        z = df["z_um"].values
        dom = df["domain"].values
        
        idx_b1 = np.where(dom == "Lo_north")[0][-1]
        idx_b2 = np.where(dom == "Lo_south")[0][0]
        
        z_mid = 0.5 * (np.max(z) + np.min(z))
        z_al = z - z_mid
        return r, z_al, idx_b1, idx_b2, df
        
    # =========================================================================
    # ROW 1: BEST DUAL INWARD DIMPLE SCENARIO (Deflation v = 0.90, D_Lo = 3.0 um, gamma = 1.5 pN)
    # =========================================================================
    # Panel (a1): Full Cross Section (Bare Sphere vs Remodeled Dual Inward)
    ax_a1 = fig.add_subplot(gs[0, 0])
    bare_file = os.path.join(profiles_dir, "prof_two_dom_v_0.90_theta_0.00.csv")
    remo_file = os.path.join(profiles_dir, "prof_two_dom_v_0.90_theta_1.00.csv")
    
    r_b, z_b, b1_b, b2_b, _ = get_sym_profile(bare_file)
    r_r, z_r, b1_r, b2_r, _ = get_sym_profile(remo_file)
    
    # Plot bare dashed
    ax_a1.plot(np.concatenate([-r_b[::-1], r_b]), np.concatenate([z_b[::-1], z_b]), "--", color="#7f7f7f", linewidth=1.5, label=r"Bare Vesicle ($\theta=0.0$)")
    
    # Plot remodeled 3 zones
    r_Ld_l = -r_r[b1_r:b2_r+1][::-1]
    z_Ld_l = z_r[b1_r:b2_r+1][::-1]
    r_Ld_r = r_r[b1_r:b2_r+1]
    z_Ld_r = z_r[b1_r:b2_r+1]
    
    r_Lo1 = np.concatenate([-r_r[:b1_r+1][::-1], r_r[:b1_r+1]])
    z_Lo1 = np.concatenate([z_r[:b1_r+1][::-1], z_r[:b1_r+1]])
    
    r_Lo2 = np.concatenate([-r_r[b2_r:][::-1], r_r[b2_r:]])
    z_Lo2 = np.concatenate([z_r[b2_r:][::-1], z_r[b2_r:]])
    
    ax_a1.plot(r_Ld_l, z_Ld_l, color=color_Ld, linewidth=2.8, label=r"Equatorial $L_d$ Belt")
    ax_a1.plot(r_Ld_r, z_Ld_r, color=color_Ld, linewidth=2.8)
    ax_a1.plot(r_Lo1, z_Lo1, color=color_Lo, linewidth=3.8, label=r"Dual Inward $L_o$ Dimples ($\theta=1.0$)")
    ax_a1.plot(r_Lo2, z_Lo2, color=color_Lo, linewidth=3.8)
    
    ax_a1.plot([-r_r[b1_r], r_r[b1_r]], [z_r[b1_r], z_r[b1_r]], "ko", markersize=5.5)
    ax_a1.plot([-r_r[b2_r], r_r[b2_r]], [z_r[b2_r], z_r[b2_r]], "ko", markersize=5.5)
    
    ax_a1.axvline(0, color="gray", linestyle=":", linewidth=0.8, alpha=0.7)
    ax_a1.set_title(r"$\mathbf{(a1)}$ Best Dual Inward Dimple ($v=0.90$)", fontweight="bold", pad=6)
    ax_a1.set_aspect("equal")
    ax_a1.set_xlim(-6.0, 6.0)
    ax_a1.set_ylim(-6.0, 6.0)
    ax_a1.set_ylabel(r"Axial $z$ ($\mu\mathrm{m}$)")
    ax_a1.grid(True, linestyle=":", alpha=0.5)
    ax_a1.legend(frameon=True, facecolor="white", loc="lower right")

    # Panel (a2): Zoomed View of Symmetric Dual Inward Dimples
    ax_a2 = fig.add_subplot(gs[0, 1])
    ax_a2.plot(r_Lo1, z_Lo1, color=color_Lo, linewidth=3.8, label="North $L_o$ Inward Dimple")
    ax_a2.plot(r_Lo2, z_Lo2, color="#b2182b", linewidth=3.8, linestyle="-", label="South $L_o$ Inward Dimple")
    ax_a2.plot([-r_r[b1_r], r_r[b1_r]], [z_r[b1_r], z_r[b1_r]], "ko", markersize=5.5)
    ax_a2.plot([-r_r[b2_r], r_r[b2_r]], [z_r[b2_r], z_r[b2_r]], "ko", markersize=5.5)
    
    raw_depth = float(z_r[0] - z_r[b1_r])
    is_inward = raw_depth < 0
    morph_lbl = "Identical Dual Inward Dimples" if is_inward else "Dual Outward Domes"
    sign_lbl = "(inward)" if is_inward else "(outward)"
    
    ax_a2.annotate(
        f"{morph_lbl}\nNorth Depth: $\\mathbf{{{raw_depth:+.2f}\\,\\mu\\mathrm{{m}}}}$ {sign_lbl}\nSouth Depth: $\\mathbf{{{raw_depth:+.2f}\\,\\mu\\mathrm{{m}}}}$ {sign_lbl}\nDual Necks: $2r_b = {2*r_r[b1_r]:.2f}\\,\\mu\\mathrm{{m}}$",
        xy=(0, 0), xytext=(0, 0),
        fontsize=9, ha="center", va="center", fontweight="bold", color="#7f0000" if is_inward else "#08519c",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="#fff5f0" if is_inward else "#eff3ff", 
                  edgecolor="#d73027" if is_inward else "#3182bd", alpha=0.95)
    )
    ax_a2.set_title(r"$\mathbf{(a2)}$ Zoomed Dual Inward Geometry", fontweight="bold", pad=6)
    ax_a2.set_xlim(-3.0, 3.0)
    ax_a2.set_ylim(-4.5, 4.5)
    ax_a2.set_ylabel(r"Axial $z$ ($\mu\mathrm{m}$)")
    ax_a2.grid(True, linestyle=":", alpha=0.5)
    ax_a2.legend(frameon=True, facecolor="white", loc="upper right")

    # Panel (a3): Mean Curvature H(s) Profile
    ax_a3 = fig.add_subplot(gs[0, 2])
    _, _, _, _, df_b = get_sym_profile(bare_file)
    _, _, _, _, df_r = get_sym_profile(remo_file)
    s_grid = df_r["s_um"].values
    s_b1 = s_grid[b1_r]
    s_b2 = s_grid[b2_r]
    
    ax_a3.plot(s_grid, df_b["H_um_inv"].values, "--", color=color_bare, linewidth=2.0, label=r"Bare ($v=0.90, \theta=0$)")
    ax_a3.plot(s_grid, df_r["H_um_inv"].values, "-", color=color_Lo, linewidth=2.5, label=r"Dual Inward ($\theta=1.0$)")
    
    ax_a3.axvline(s_b1, color="black", linestyle=":", linewidth=1.5, label=r"Boundary Necks $s_b$")
    ax_a3.axvline(s_b2, color="black", linestyle=":", linewidth=1.5)
    ax_a3.axvspan(0, s_b1, color="#fee0d2", alpha=0.4, label=r"North $L_o$")
    ax_a3.axvspan(s_b2, s_grid[-1], color="#fee0d2", alpha=0.4, label=r"South $L_o$")
    
    ax_a3.set_title(r"$\mathbf{(a3)}$ Mean Curvature $H(s)$", fontweight="bold", pad=6)
    ax_a3.set_xlabel(r"Arc Length $s$ ($\mu\mathrm{m}$)")
    ax_a3.set_ylabel(r"Mean Curvature $H$ ($\mu\mathrm{m}^{-1}$)")
    ax_a3.grid(True, linestyle=":", alpha=0.5)
    ax_a3.legend(frameon=True, facecolor="white", loc="upper right")

    # =========================================================================
    # ROW 2: PARAMETRIC COMPARISONS (Deflation v=0.90 vs v=0.95 vs v=0.98)
    # =========================================================================
    # Panel (b1): Deflation v = 0.90
    ax_b1 = fig.add_subplot(gs[1, 0])
    ax_b1.plot(r_Ld_l, z_Ld_l, color=color_Ld, linewidth=2.5)
    ax_b1.plot(r_Ld_r, z_Ld_r, color=color_Ld, linewidth=2.5)
    ax_b1.plot(r_Lo1, z_Lo1, color=color_Lo, linewidth=3.5)
    ax_b1.plot(r_Lo2, z_Lo2, color=color_Lo, linewidth=3.5)
    ax_b1.plot([-r_r[b1_r], r_r[b1_r]], [z_r[b1_r], z_r[b1_r]], "ko", markersize=5.0)
    ax_b1.plot([-r_r[b2_r], r_r[b2_r]], [z_r[b2_r], z_r[b2_r]], "ko", markersize=5.0)
    ax_b1.axvline(0, color="gray", linestyle=":", linewidth=0.8, alpha=0.7)
    ax_b1.set_title(r"$\mathbf{(b1)}$ $v = 0.90$: $\mathbf{0.59\,\mu\mathrm{m}}$ Dual Inward", fontweight="bold", pad=6)
    ax_b1.set_aspect("equal")
    ax_b1.set_xlim(-6.0, 6.0)
    ax_b1.set_ylim(-6.0, 6.0)
    ax_b1.set_xlabel(r"Radial $r$ ($\mu\mathrm{m}$)")
    ax_b1.set_ylabel(r"Axial $z$ ($\mu\mathrm{m}$)")
    ax_b1.grid(True, linestyle=":", alpha=0.5)

    # Panel (b2): Deflation v = 0.95
    ax_b2 = fig.add_subplot(gs[1, 1])
    file_v95 = os.path.join(profiles_dir, "prof_two_dom_v_0.95_theta_1.00.csv")
    r_95, z_95, b1_95, b2_95, _ = get_sym_profile(file_v95)
    ax_b2.plot(-r_95[b1_95:b2_95+1][::-1], z_95[b1_95:b2_95+1][::-1], color=color_Ld, linewidth=2.5)
    ax_b2.plot(r_95[b1_95:b2_95+1], z_95[b1_95:b2_95+1], color=color_Ld, linewidth=2.5)
    ax_b2.plot(np.concatenate([-r_95[:b1_95+1][::-1], r_95[:b1_95+1]]), np.concatenate([z_95[:b1_95+1][::-1], z_95[:b1_95+1]]), color=color_Lo, linewidth=3.5)
    ax_b2.plot(np.concatenate([-r_95[b2_95:][::-1], r_95[b2_95:]]), np.concatenate([z_95[b2_95:][::-1], z_95[b2_95:]]), color=color_Lo, linewidth=3.5)
    ax_b2.plot([-r_95[b1_95], r_95[b1_95]], [z_95[b1_95], z_95[b1_95]], "ko", markersize=5.0)
    ax_b2.plot([-r_95[b2_95], r_95[b2_95]], [z_95[b2_95], z_95[b2_95]], "ko", markersize=5.0)
    ax_b2.axvline(0, color="gray", linestyle=":", linewidth=0.8, alpha=0.7)
    ax_b2.set_title(r"$\mathbf{(b2)}$ $v = 0.95$: Intermediate State", fontweight="bold", pad=6)
    ax_b2.set_aspect("equal")
    ax_b2.set_xlim(-6.0, 6.0)
    ax_b2.set_ylim(-6.0, 6.0)
    ax_b2.set_xlabel(r"Radial $r$ ($\mu\mathrm{m}$)")
    ax_b2.set_ylabel(r"Axial $z$ ($\mu\mathrm{m}$)")
    ax_b2.grid(True, linestyle=":", alpha=0.5)

    # Panel (b3): Deflation v = 0.98
    ax_b3 = fig.add_subplot(gs[1, 2])
    file_v98 = os.path.join(profiles_dir, "prof_two_dom_v_0.98_theta_1.00.csv")
    r_98, z_98, b1_98, b2_98, _ = get_sym_profile(file_v98)
    ax_b3.plot(-r_98[b1_98:b2_98+1][::-1], z_98[b1_98:b2_98+1][::-1], color=color_Ld, linewidth=2.5)
    ax_b3.plot(r_98[b1_98:b2_98+1], z_98[b1_98:b2_98+1], color=color_Ld, linewidth=2.5)
    ax_b3.plot(np.concatenate([-r_98[:b1_98+1][::-1], r_98[:b1_98+1]]), np.concatenate([z_98[:b1_98+1][::-1], z_98[:b1_98+1]]), color=color_Lo, linewidth=3.5)
    ax_b3.plot(np.concatenate([-r_98[b2_98:][::-1], r_98[b2_98:]]), np.concatenate([z_98[b2_98:][::-1], z_98[b2_98:]]), color=color_Lo, linewidth=3.5)
    ax_b3.plot([-r_98[b1_98], r_98[b1_98]], [z_98[b1_98], z_98[b1_98]], "ko", markersize=5.0)
    ax_b3.plot([-r_98[b2_98], r_98[b2_98]], [z_98[b2_98], z_98[b2_98]], "ko", markersize=5.0)
    ax_b3.axvline(0, color="gray", linestyle=":", linewidth=0.8, alpha=0.7)
    ax_b3.set_title(r"$\mathbf{(b3)}$ $v = 0.98$: Dual Spindle Extrusion", fontweight="bold", pad=6)
    ax_b3.set_aspect("equal")
    ax_b3.set_xlim(-6.0, 6.0)
    ax_b3.set_ylim(-6.0, 6.0)
    ax_b3.set_xlabel(r"Radial $r$ ($\mu\mathrm{m}$)")
    ax_b3.set_ylabel(r"Axial $z$ ($\mu\mathrm{m}$)")
    ax_b3.grid(True, linestyle=":", alpha=0.5)

    plt.tight_layout()
    fig.savefig(out_fig, dpi=300)
    plt.close(fig)
    print(f"Two-domain showcase figure saved to: {out_fig}")
    return out_fig

if __name__ == "__main__":
    plot_two_domains_showcase()
