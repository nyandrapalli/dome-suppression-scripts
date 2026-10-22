"""
Plotting Module for Two-Domain HermesLoouter (Protein on outer leaflet of Lo caps).
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

def plot_two_domains(results_dir="HermesLoouter/double"):
    profiles_dir = os.path.join(results_dir, "profiles")
    out_fig = os.path.join(results_dir, "fig_two_domains_sweeps.png")
    
    color_Lo = "#d73027"  # Red for Lo
    color_Ld = "#1a9850"  # Green for Ld
    
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
    
    # 1. Line Tension Sweep
    gamma_vals = [0.5, 1.5, 3.0, 5.0]
    for col_idx, gamma in enumerate(gamma_vals):
        ax = fig.add_subplot(gs[0, col_idx])
        prof_path = _find_profile(profiles_dir, [f"prof_two_dom_gamma_{gamma:.1f}_theta_1.00.csv"])
        if prof_path:
            df = pd.read_csv(prof_path)
            r, z, dom = df["r_um"].values, df["z_um"].values, df["domain"].values
            z_aligned = z - 0.5 * (np.max(z) + np.min(z))
            
            mask_n = (dom == "Lo_north")
            mask_s = (dom == "Lo_south")
            mask_ld = (dom == "Ld")
            
            # Plot Ld equatorial belt
            ax.plot(-r[mask_ld][::-1], z_aligned[mask_ld][::-1], color=color_Ld, linewidth=2.5, label=r"$L_d$ Belt (Bare)" if col_idx == 0 else "")
            ax.plot(r[mask_ld], z_aligned[mask_ld], color=color_Ld, linewidth=2.5)
            
            # Plot Lo caps
            ax.plot(np.concatenate([-r[mask_n][::-1], r[mask_n]]), np.concatenate([z_aligned[mask_n][::-1], z_aligned[mask_n]]), color=color_Lo, linewidth=3.5, label=r"$L_o$ Caps (Outer Prot)" if col_idx == 0 else "")
            ax.plot(np.concatenate([-r[mask_s][::-1], r[mask_s]]), np.concatenate([z_aligned[mask_s][::-1], z_aligned[mask_s]]), color=color_Lo, linewidth=3.5)
            
        ax.axvline(0, color="gray", linestyle=":", linewidth=0.8, alpha=0.7)
        ax.set_title(f"(a{col_idx+1}) $\\gamma = {gamma:.1f}\\,\\mathrm{{pN}}$ ($\\theta=1.0$)", fontweight="bold", pad=4)
        ax.set_aspect("equal")
        ax.set_xlim(-6.0, 6.0)
        ax.set_ylim(-6.0, 6.0)
        ax.set_ylabel(r"Axial $z$ ($\mu\mathrm{m}$)")
        ax.grid(True, linestyle=":", alpha=0.5)
        if col_idx == 0:
            ax.legend(loc="lower left", framealpha=0.9)

    # 2. Domain Size Sweep
    d_lo_vals = [1.0, 3.0, 5.0, 7.0]
    for col_idx, d_lo in enumerate(d_lo_vals):
        ax = fig.add_subplot(gs[1, col_idx])
        prof_path = _find_profile(profiles_dir, [f"prof_two_dom_DLo_{d_lo:.1f}_theta_1.00.csv"])
        if prof_path:
            df = pd.read_csv(prof_path)
            r, z, dom = df["r_um"].values, df["z_um"].values, df["domain"].values
            z_aligned = z - 0.5 * (np.max(z) + np.min(z))
            
            mask_n = (dom == "Lo_north")
            mask_s = (dom == "Lo_south")
            mask_ld = (dom == "Ld")
            
            ax.plot(-r[mask_ld][::-1], z_aligned[mask_ld][::-1], color=color_Ld, linewidth=2.5)
            ax.plot(r[mask_ld], z_aligned[mask_ld], color=color_Ld, linewidth=2.5)
            ax.plot(np.concatenate([-r[mask_n][::-1], r[mask_n]]), np.concatenate([z_aligned[mask_n][::-1], z_aligned[mask_n]]), color=color_Lo, linewidth=3.5)
            ax.plot(np.concatenate([-r[mask_s][::-1], r[mask_s]]), np.concatenate([z_aligned[mask_s][::-1], z_aligned[mask_s]]), color=color_Lo, linewidth=3.5)
            
        ax.axvline(0, color="gray", linestyle=":", linewidth=0.8, alpha=0.7)
        ax.set_title(f"(b{col_idx+1}) $D_{{Lo}} = {d_lo:.1f}\\,\\mu\\mathrm{{m}}$ ($\\theta=1.0$)", fontweight="bold", pad=4)
        ax.set_aspect("equal")
        ax.set_xlim(-6.0, 6.0)
        ax.set_ylim(-6.0, 6.0)
        ax.set_ylabel(r"Axial $z$ ($\mu\mathrm{m}$)")
        ax.grid(True, linestyle=":", alpha=0.5)

    # 3. Deflation Sweep
    v_vals = [1.00, 0.98, 0.95, 0.90]
    for col_idx, v_val in enumerate(v_vals):
        ax = fig.add_subplot(gs[2, col_idx])
        prof_path = _find_profile(profiles_dir, [f"prof_two_dom_v_{v_val:.2f}_theta_1.00.csv"])
        if prof_path:
            df = pd.read_csv(prof_path)
            r, z, dom = df["r_um"].values, df["z_um"].values, df["domain"].values
            z_aligned = z - 0.5 * (np.max(z) + np.min(z))
            
            mask_n = (dom == "Lo_north")
            mask_s = (dom == "Lo_south")
            mask_ld = (dom == "Ld")
            
            ax.plot(-r[mask_ld][::-1], z_aligned[mask_ld][::-1], color=color_Ld, linewidth=2.5)
            ax.plot(r[mask_ld], z_aligned[mask_ld], color=color_Ld, linewidth=2.5)
            ax.plot(np.concatenate([-r[mask_n][::-1], r[mask_n]]), np.concatenate([z_aligned[mask_n][::-1], z_aligned[mask_n]]), color=color_Lo, linewidth=3.5)
            ax.plot(np.concatenate([-r[mask_s][::-1], r[mask_s]]), np.concatenate([z_aligned[mask_s][::-1], z_aligned[mask_s]]), color=color_Lo, linewidth=3.5)
            
        ax.axvline(0, color="gray", linestyle=":", linewidth=0.8, alpha=0.7)
        ax.set_title(f"(c{col_idx+1}) $v = {v_val:.2f}$ ($\\theta=1.0$)", fontweight="bold", pad=4)
        ax.set_aspect("equal")
        ax.set_xlim(-6.0, 6.0)
        ax.set_ylim(-6.0, 6.0)
        ax.set_xlabel(r"Radial $r$ ($\mu\mathrm{m}$)")
        ax.set_ylabel(r"Axial $z$ ($\mu\mathrm{m}$)")
        ax.grid(True, linestyle=":", alpha=0.5)

    plt.suptitle("HermesLoouter: Two-Domain Parameter Sweeps (Protein on Outer Leaflet of Both $L_o$ Caps)", fontsize=14, fontweight="bold", y=0.98)
    fig.savefig(out_fig, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved figure: {out_fig}")

if __name__ == "__main__":
    plot_two_domains()
