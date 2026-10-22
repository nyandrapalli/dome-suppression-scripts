"""
Publication Figure Showcase for Optimal Invagination Regimes (PRL Standard).
Presents:
1. Regime 1: Domain Size Effect (D_Lo = 4.5 um, gamma = 1.5 pN, v = 0.98)
2. Regime 2: Constricted Neck Pouch (D_Lo = 3.0 um, gamma = 1.5 pN, v = 0.95)
Includes full-vesicle cross-sections, zoomed domain views, and mean curvature profiles H(s).
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

def _load_profile(profiles_dir, candidate_names):
    for name in candidate_names:
        p = os.path.join(profiles_dir, name)
        if os.path.exists(p):
            return pd.read_csv(p)
    return None

def plot_best_invaginations(out_path="simulation_results_prl/fig_best_invaginations_showcase.png"):
    profiles_dir = "simulation_results_prl/param_sweeps/profiles"
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    
    color_Lo = "#d73027"       # Vibrant crimson for Lo domain
    color_Ld = "#1a9850"       # Vibrant emerald green for Ld matrix
    
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
    gs = GridSpec(2, 3, figure=fig, hspace=0.32, wspace=0.28, width_ratios=[1.2, 1.0, 1.1])
    
    # =========================================================================
    # ROW 1: REGIME 1 — (D_Lo = 4.5 um or 5.0 um, v = 0.98)
    # =========================================================================
    df_bare_1 = _load_profile(profiles_dir, ["prof_DLo_4.5_theta_0.00.csv", "prof_DLo_4.5um_theta_0.00.csv", "prof_DLo_5.0_theta_0.00.csv"])
    df_remo_1 = _load_profile(profiles_dir, ["prof_DLo_4.5_theta_1.00.csv", "prof_DLo_4.5um_theta_1.00.csv", "prof_DLo_5.0_theta_1.00.csv"])
    
    if df_bare_1 is not None and df_remo_1 is not None:
        ax_a1 = fig.add_subplot(gs[0, 0])
        z_mid_bare_1 = 0.5 * (df_bare_1["z_um"].max() + df_bare_1["z_um"].min())
        z_bare_1 = df_bare_1["z_um"].values - z_mid_bare_1
        r_bare_1 = df_bare_1["r_um"].values
        
        z_mid_remo_1 = 0.5 * (df_remo_1["z_um"].max() + df_remo_1["z_um"].min())
        z_remo_1 = df_remo_1["z_um"].values - z_mid_remo_1
        r_remo_1 = df_remo_1["r_um"].values
        dom_remo_1 = df_remo_1["domain"].values
        n_Lo_1 = np.sum(dom_remo_1 == "Lo")
        
        r_bare_sym_1 = np.concatenate([-r_bare_1[::-1], r_bare_1])
        z_bare_sym_1 = np.concatenate([z_bare_1[::-1], z_bare_1])
        ax_a1.plot(r_bare_sym_1, z_bare_sym_1, "--", color="#7f7f7f", linewidth=1.5, label=r"Bare State ($\theta=0.0$)")
        
        r_Lo_1 = np.concatenate([-r_remo_1[:n_Lo_1][::-1], r_remo_1[:n_Lo_1]])
        z_Lo_1 = np.concatenate([z_remo_1[:n_Lo_1][::-1], z_remo_1[:n_Lo_1]])
        r_Ld_1_l = -r_remo_1[n_Lo_1-1:][::-1]
        z_Ld_1_l = z_remo_1[n_Lo_1-1:][::-1]
        r_Ld_1_r = r_remo_1[n_Lo_1-1:]
        z_Ld_1_r = z_remo_1[n_Lo_1-1:]
        
        ax_a1.plot(r_Ld_1_l, z_Ld_1_l, color=color_Ld, linewidth=2.8, label=r"Remodeled $L_d$ Matrix")
        ax_a1.plot(r_Ld_1_r, z_Ld_1_r, color=color_Ld, linewidth=2.8)
        ax_a1.plot(r_Lo_1, z_Lo_1, color=color_Lo, linewidth=3.8, label=r"Invaginated $L_o$ Domain")
        ax_a1.plot([-r_remo_1[n_Lo_1-1], r_remo_1[n_Lo_1-1]], [z_remo_1[n_Lo_1-1], z_remo_1[n_Lo_1-1]], "ko", markersize=5.5)
        
        depth_1 = z_remo_1[0] - z_remo_1[n_Lo_1-1]
        ax_a1.axvline(0, color="gray", linestyle=":", linewidth=0.8, alpha=0.7)
        ax_a1.set_title(r"$\mathbf{(a1)}$ Single-Domain Full Morphology", fontweight="bold", pad=8)
        ax_a1.set_aspect("equal")
        ax_a1.set_xlim(-6.5, 6.5)
        ax_a1.set_ylim(-6.5, 6.5)
        ax_a1.set_ylabel(r"Axial $z$ ($\mu\mathrm{m}$)")
        ax_a1.grid(True, linestyle=":", alpha=0.5)
        ax_a1.legend(frameon=True, facecolor="white", loc="lower right")
        
        # Panel (a2): Zoom
        ax_a2 = fig.add_subplot(gs[0, 1])
        ax_a2.plot(r_Lo_1, z_Lo_1, color=color_Lo, linewidth=3.8, label=r"$L_o$ Invagination")
        ax_a2.plot(r_Ld_1_l, z_Ld_1_l, color=color_Ld, linewidth=2.5)
        ax_a2.plot(r_Ld_1_r, z_Ld_1_r, color=color_Ld, linewidth=2.5)
        ax_a2.plot([-r_remo_1[n_Lo_1-1], r_remo_1[n_Lo_1-1]], [z_remo_1[n_Lo_1-1], z_remo_1[n_Lo_1-1]], "ko", markersize=5.5)
        
        ax_a2.annotate(
            f"Invagination Depth: $\\mathbf{{{depth_1:+.2f}\\,\\mu\\mathrm{{m}}}}$\nNeck Diameter: $2r_b = {2*r_remo_1[n_Lo_1-1]:.2f}\\,\\mu\\mathrm{{m}}$",
            xy=(0, z_remo_1[0]), xytext=(0, z_remo_1[0] - 1.2),
            fontsize=9.5, ha="center", fontweight="bold", color="#7f0000" if depth_1 <= 0 else "#08519c",
            bbox=dict(boxstyle="round,pad=0.35", facecolor="#fff5f0" if depth_1 <= 0 else "#eff3ff", edgecolor="#d73027" if depth_1 <= 0 else "#3182bd", alpha=0.95)
        )
        ax_a2.set_title(r"$\mathbf{(a2)}$ Zoomed Cap Profile", fontweight="bold", pad=8)
        ax_a2.set_xlim(-3.5, 3.5)
        ax_a2.set_ylim(-5.5, 5.5)
        ax_a2.grid(True, linestyle=":", alpha=0.5)

    # =========================================================================
    # ROW 2: REGIME 2 — CONSTRICTED NECK POUCH (v = 0.95, D_Lo = 3.0 um)
    # =========================================================================
    df_bare_2 = _load_profile(profiles_dir, ["prof_v_0.95_theta_0.00.csv", "prof_v_0.98_theta_0.00.csv"])
    df_remo_2 = _load_profile(profiles_dir, ["prof_v_0.95_theta_1.00.csv", "prof_v_0.98_theta_1.00.csv"])
    
    if df_bare_2 is not None and df_remo_2 is not None:
        ax_b1 = fig.add_subplot(gs[1, 0])
        z_mid_bare_2 = 0.5 * (df_bare_2["z_um"].max() + df_bare_2["z_um"].min())
        z_bare_2 = df_bare_2["z_um"].values - z_mid_bare_2
        r_bare_2 = df_bare_2["r_um"].values
        
        z_mid_remo_2 = 0.5 * (df_remo_2["z_um"].max() + df_remo_2["z_um"].min())
        z_remo_2 = df_remo_2["z_um"].values - z_mid_remo_2
        r_remo_2 = df_remo_2["r_um"].values
        dom_remo_2 = df_remo_2["domain"].values
        n_Lo_2 = np.sum(dom_remo_2 == "Lo")
        
        r_bare_sym_2 = np.concatenate([-r_bare_2[::-1], r_bare_2])
        z_bare_sym_2 = np.concatenate([z_bare_2[::-1], z_bare_2])
        ax_b1.plot(r_bare_sym_2, z_bare_sym_2, "--", color="#7f7f7f", linewidth=1.5, label=r"Bare State ($\theta=0.0$)")
        
        r_Lo_2 = np.concatenate([-r_remo_2[:n_Lo_2][::-1], r_remo_2[:n_Lo_2]])
        z_Lo_2 = np.concatenate([z_remo_2[:n_Lo_2][::-1], z_remo_2[:n_Lo_2]])
        r_Ld_2_l = -r_remo_2[n_Lo_2-1:][::-1]
        z_Ld_2_l = z_remo_2[n_Lo_2-1:][::-1]
        r_Ld_2_r = r_remo_2[n_Lo_2-1:]
        z_Ld_2_r = z_remo_2[n_Lo_2-1:]
        
        ax_b1.plot(r_Ld_2_l, z_Ld_2_l, color=color_Ld, linewidth=2.8)
        ax_b1.plot(r_Ld_2_r, z_Ld_2_r, color=color_Ld, linewidth=2.8)
        ax_b1.plot(r_Lo_2, z_Lo_2, color=color_Lo, linewidth=3.8)
        ax_b1.plot([-r_remo_2[n_Lo_2-1], r_remo_2[n_Lo_2-1]], [z_remo_2[n_Lo_2-1], z_remo_2[n_Lo_2-1]], "ko", markersize=5.5)
        
        depth_2 = z_remo_2[0] - z_remo_2[n_Lo_2-1]
        ax_b1.axvline(0, color="gray", linestyle=":", linewidth=0.8, alpha=0.7)
        ax_b1.set_title(r"$\mathbf{(b1)}$ Constricted Pouch Morphology", fontweight="bold", pad=8)
        ax_b1.set_aspect("equal")
        ax_b1.set_xlim(-6.5, 6.5)
        ax_b1.set_ylim(-6.5, 6.5)
        ax_b1.set_xlabel(r"Radial $r$ ($\mu\mathrm{m}$)")
        ax_b1.set_ylabel(r"Axial $z$ ($\mu\mathrm{m}$)")
        ax_b1.grid(True, linestyle=":", alpha=0.5)
        
        # Panel (b2): Zoom
        ax_b2 = fig.add_subplot(gs[1, 1])
        ax_b2.plot(r_Lo_2, z_Lo_2, color=color_Lo, linewidth=3.8)
        ax_b2.plot(r_Ld_2_l, z_Ld_2_l, color=color_Ld, linewidth=2.5)
        ax_b2.plot(r_Ld_2_r, z_Ld_2_r, color=color_Ld, linewidth=2.5)
        ax_b2.plot([-r_remo_2[n_Lo_2-1], r_remo_2[n_Lo_2-1]], [z_remo_2[n_Lo_2-1], z_remo_2[n_Lo_2-1]], "ko", markersize=5.5)
        
        ax_b2.annotate(
            f"Pouch Depth: $\\mathbf{{{depth_2:+.2f}\\,\\mu\\mathrm{{m}}}}$\nNeck Diameter: $2r_b = {2*r_remo_2[n_Lo_2-1]:.2f}\\,\\mu\\mathrm{{m}}$",
            xy=(0, z_remo_2[0]), xytext=(0, z_remo_2[0] - 1.2),
            fontsize=9.5, ha="center", fontweight="bold", color="#7f0000" if depth_2 <= 0 else "#08519c",
            bbox=dict(boxstyle="round,pad=0.35", facecolor="#fff5f0" if depth_2 <= 0 else "#eff3ff", edgecolor="#d73027" if depth_2 <= 0 else "#3182bd", alpha=0.95)
        )
        ax_b2.set_title(r"$\mathbf{(b2)}$ Zoomed Constricted Neck", fontweight="bold", pad=8)
        ax_b2.set_xlim(-3.5, 3.5)
        ax_b2.set_ylim(-5.5, 5.5)
        ax_b2.set_xlabel(r"Radial $r$ ($\mu\mathrm{m}$)")
        ax_b2.grid(True, linestyle=":", alpha=0.5)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Showcase figure saved to: {out_path}")

if __name__ == "__main__":
    plot_best_invaginations()
