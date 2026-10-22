"""
Plotting Module for Two Opposite Lo Domains (Bipolar Vesicle - PRL Standard).
Renders publication figures showing:
- Panel A: Line Tension Sweep (gamma in {0.5, 1.5, 3.0, 5.0} pN)
- Panel B: Domain Size Sweep (D_Lo in {1.0, 3.0, 5.0, 7.0} um)
- Panel C: Deflation Factor Sweep (v in {0.90, 0.95, 0.98, 1.00})
North and South Lo caps are highlighted in bold RED, Equatorial Ld matrix in GREEN,
with boundary contact points marked by black dots.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

def _find_profile(profiles_dir, names):
    for name in names:
        p = os.path.join(profiles_dir, name)
        if os.path.exists(p):
            return p
    return None

def plot_two_domains_sweeps(results_dir="simulation_results_prl/two_domains_sweeps"):
    profiles_dir = os.path.join(results_dir, "profiles")
    out_fig = os.path.join(results_dir, "fig_two_domains_comprehensive.png")
    
    color_Lo = "#d73027"  # Deep red for Lo caps
    color_Ld = "#1a9850"  # Forest green for Ld matrix
    
    plt.rcParams.update({
        "font.size": 9.5,
        "axes.labelsize": 10.5,
        "axes.titlesize": 11.5,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 8.5,
        "figure.dpi": 300
    })
    
    fig = plt.figure(figsize=(18, 14))
    gs = GridSpec(3, 4, figure=fig, hspace=0.32, wspace=0.25)
    
    def render_bipolar_profile(ax, df, title_str, param_str, is_first=False):
        r = df["r_um"].values
        z = df["z_um"].values
        dom = df["domain"].values
        
        mask_lo_north = (dom == "Lo_north")
        mask_ld = (dom == "Ld")
        mask_lo_south = (dom == "Lo_south")
        
        idx_b1 = np.where(mask_lo_north)[0][-1]
        idx_b2 = np.where(mask_lo_south)[0][0]
        
        z_mid = 0.5 * (np.max(z) + np.min(z))
        z_al = z - z_mid
        
        # Symmetrize 3 zones
        # 1. North Lo cap
        r_Lo1 = np.concatenate([-r[:idx_b1+1][::-1], r[:idx_b1+1]])
        z_Lo1 = np.concatenate([z_al[:idx_b1+1][::-1], z_al[:idx_b1+1]])
        
        # 2. Ld Matrix Belt
        r_Ld_l = -r[idx_b1:idx_b2+1][::-1]
        z_Ld_l = z_al[idx_b1:idx_b2+1][::-1]
        r_Ld_r = r[idx_b1:idx_b2+1]
        z_Ld_r = z_al[idx_b1:idx_b2+1]
        
        # 3. South Lo cap
        r_Lo2 = np.concatenate([-r[idx_b2:][::-1], r[idx_b2:]])
        z_Lo2 = np.concatenate([z_al[idx_b2:][::-1], z_al[idx_b2:]])
        
        ax.plot(r_Ld_l, z_Ld_l, color=color_Ld, linewidth=2.5, label=r"Equatorial $L_d$ Belt" if is_first else "")
        ax.plot(r_Ld_r, z_Ld_r, color=color_Ld, linewidth=2.5)
        ax.plot(r_Lo1, z_Lo1, color=color_Lo, linewidth=3.5, label=r"Bipolar $L_o$ Caps" if is_first else "")
        ax.plot(r_Lo2, z_Lo2, color=color_Lo, linewidth=3.5)
        ax.plot([-r[idx_b1], r[idx_b1]], [z_al[idx_b1], z_al[idx_b1]], "ko", markersize=4.5)
        ax.plot([-r[idx_b2], r[idx_b2]], [z_al[idx_b2], z_al[idx_b2]], "ko", markersize=4.5)
        
        depth_n = z_al[0] - z_al[idx_b1]
        depth_s = z_al[idx_b2] - z_al[-1]
        
        txt_n = f"North: {'Inward ' + str(abs(depth_n))[:4] if depth_n <= 0 else 'Dome +' + str(depth_n)[:4]} $\\mu\\mathrm{{m}}$"
        txt_s = f"South: {'Inward ' + str(abs(depth_s))[:4] if depth_s <= 0 else 'Dome +' + str(depth_s)[:4]} $\\mu\\mathrm{{m}}$"
        
        ax.annotate(
            f"{param_str}\n{txt_n}\n{txt_s}",
            xy=(0, 0), xytext=(0, 0),
            fontsize=8, ha="center", va="center", fontweight="bold", color="#7f0000" if depth_n <= 0 else "#08519c",
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.9, edgecolor="#cccccc")
        )
        
        ax.axvline(0, color="gray", linestyle=":", linewidth=0.8, alpha=0.7)
        ax.set_title(title_str, fontweight="bold", pad=4)
        ax.set_aspect("equal")
        ax.set_xlim(-6.0, 6.0)
        ax.set_ylim(-6.0, 6.0)
        ax.set_ylabel(r"Axial $z$ ($\mu\mathrm{m}$)")
        ax.grid(True, linestyle=":", alpha=0.5)
        if is_first:
            ax.legend(frameon=True, facecolor="white", loc="lower right")

    # -------------------------------------------------------------
    # Row 1: Line Tension Sweep (gamma = 0.5, 1.5, 3.0, 5.0 pN)
    # -------------------------------------------------------------
    gamma_vals = [0.5, 1.5, 3.0, 5.0]
    for col_idx, gamma in enumerate(gamma_vals):
        ax = fig.add_subplot(gs[0, col_idx])
        prof_path = _find_profile(profiles_dir, [
            f"prof_two_dom_gamma_{gamma:.1f}_theta_1.00.csv",
            f"prof_two_dom_gamma_{gamma:.1f}pN_theta_1.00.csv"
        ])
        if prof_path:
            df = pd.read_csv(prof_path)
            render_bipolar_profile(ax, df, f"(a{col_idx+1}) $\\gamma = {gamma:.1f}\\,\\mathrm{{pN}}$ ($\\theta=1.0$)", f"$\\gamma = {gamma:.1f}\\,\\mathrm{{pN}}$", is_first=(col_idx==0))

    # -------------------------------------------------------------
    # Row 2: Domain Size Sweep (D_Lo = 1.0, 3.0, 5.0, 7.0 um)
    # -------------------------------------------------------------
    dlo_vals = [1.0, 3.0, 5.0, 7.0]
    for col_idx, D_Lo in enumerate(dlo_vals):
        ax = fig.add_subplot(gs[1, col_idx])
        prof_path = _find_profile(profiles_dir, [
            f"prof_two_dom_DLo_{D_Lo:.1f}_theta_1.00.csv",
            f"prof_two_dom_DLo_{D_Lo:.1f}um_theta_1.00.csv"
        ])
        if prof_path:
            df = pd.read_csv(prof_path)
            render_bipolar_profile(ax, df, f"(b{col_idx+1}) $D_{{Lo}} = {D_Lo:.1f}\\,\\mu\\mathrm{{m}}$ ($\\theta=1.0$)", f"$D_{{Lo}} = {D_Lo:.1f}\\,\\mu\\mathrm{{m}}$")

    # -------------------------------------------------------------
    # Row 3: Deflation Factor Sweep (v = 0.90, 0.95, 0.98, 1.00)
    # -------------------------------------------------------------
    v_vals = [0.90, 0.95, 0.98, 1.00]
    for col_idx, v in enumerate(v_vals):
        ax = fig.add_subplot(gs[2, col_idx])
        prof_path = _find_profile(profiles_dir, [
            f"prof_two_dom_v_{v:.2f}_theta_1.00.csv"
        ])
        if prof_path:
            df = pd.read_csv(prof_path)
            render_bipolar_profile(ax, df, f"(c{col_idx+1}) Deflation $v = {v:.2f}$ ($\\theta=1.0$)", f"$v = {v:.2f}$")
        ax.set_xlabel(r"Radial $r$ ($\mu\mathrm{m}$)")
        
    plt.tight_layout()
    fig.savefig(out_fig, dpi=300)
    plt.close(fig)
    print(f"Two-domain comprehensive figure saved to: {out_fig}")
    return out_fig

if __name__ == "__main__":
    plot_two_domains_sweeps()
