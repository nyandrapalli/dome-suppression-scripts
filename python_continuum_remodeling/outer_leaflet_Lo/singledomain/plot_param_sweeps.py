"""
Plotting Module for Single-Domain HermesLoouter (Protein on outer leaflet of Lo).
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

def plot_param_sweeps(results_dir=None):
    if results_dir is None:
        results_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "single")
    profiles_dir = os.path.join(results_dir, "profiles")
    out_fig = os.path.join(results_dir, "fig_single_domain_sweeps.png")
    
    color_Lo = "#d73027"  # Red for Lo (protein on outer leaflet)
    color_Ld = "#1a9850"  # Green for bare Ld
    
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
        prof_path = _find_profile(profiles_dir, [f"prof_gamma_{gamma:.1f}_theta_1.00.csv"])
        if prof_path:
            df = pd.read_csv(prof_path)
            r, z, dom = df["r_um"].values, df["z_um"].values, df["domain"].values
            mask_lo = (dom == "Lo")
            n_Lo = np.sum(mask_lo)
            z_aligned = z - 0.5 * (np.max(z) + np.min(z))
            
            z_apex, z_rim, r_b_val = z_aligned[0], z_aligned[n_Lo - 1], r[n_Lo - 1]
            elev_val = z_rim - z_apex  # or apex vs rim
            
            r_Lo = np.concatenate([-r[:n_Lo][::-1], r[:n_Lo]])
            z_Lo = np.concatenate([z_aligned[:n_Lo][::-1], z_aligned[:n_Lo]])
            r_Ld_l, z_Ld_l = -r[n_Lo-1:][::-1], z_aligned[n_Lo-1:][::-1]
            r_Ld_r, z_Ld_r = r[n_Lo-1:], z_aligned[n_Lo-1:]
            
            ax.plot(r_Ld_l, z_Ld_l, color=color_Ld, linewidth=2.5, label=r"$L_d$ Bare" if col_idx == 0 else "")
            ax.plot(r_Ld_r, z_Ld_r, color=color_Ld, linewidth=2.5)
            ax.plot(r_Lo, z_Lo, color=color_Lo, linewidth=3.5, label=r"$L_o$ (Outer Prot)" if col_idx == 0 else "")
            ax.plot([-r_b_val, r_b_val], [z_rim, z_rim], "ko", markersize=4.5)
            
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
        prof_path = _find_profile(profiles_dir, [f"prof_DLo_{d_lo:.1f}_theta_1.00.csv"])
        if prof_path:
            df = pd.read_csv(prof_path)
            r, z, dom = df["r_um"].values, df["z_um"].values, df["domain"].values
            mask_lo = (dom == "Lo")
            n_Lo = np.sum(mask_lo)
            z_aligned = z - 0.5 * (np.max(z) + np.min(z))
            
            z_rim, r_b_val = z_aligned[n_Lo - 1], r[n_Lo - 1]
            r_Lo = np.concatenate([-r[:n_Lo][::-1], r[:n_Lo]])
            z_Lo = np.concatenate([z_aligned[:n_Lo][::-1], z_aligned[:n_Lo]])
            
            ax.plot(-r[n_Lo-1:][::-1], z_aligned[n_Lo-1:][::-1], color=color_Ld, linewidth=2.5)
            ax.plot(r[n_Lo-1:], z_aligned[n_Lo-1:], color=color_Ld, linewidth=2.5)
            ax.plot(r_Lo, z_Lo, color=color_Lo, linewidth=3.5)
            ax.plot([-r_b_val, r_b_val], [z_rim, z_rim], "ko", markersize=4.5)
            
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
        prof_path = _find_profile(profiles_dir, [f"prof_v_{v_val:.2f}_theta_1.00.csv"])
        if prof_path:
            df = pd.read_csv(prof_path)
            r, z, dom = df["r_um"].values, df["z_um"].values, df["domain"].values
            mask_lo = (dom == "Lo")
            n_Lo = np.sum(mask_lo)
            z_aligned = z - 0.5 * (np.max(z) + np.min(z))
            
            z_rim, r_b_val = z_aligned[n_Lo - 1], r[n_Lo - 1]
            r_Lo = np.concatenate([-r[:n_Lo][::-1], r[:n_Lo]])
            z_Lo = np.concatenate([z_aligned[:n_Lo][::-1], z_aligned[:n_Lo]])
            
            ax.plot(-r[n_Lo-1:][::-1], z_aligned[n_Lo-1:][::-1], color=color_Ld, linewidth=2.5)
            ax.plot(r[n_Lo-1:], z_aligned[n_Lo-1:], color=color_Ld, linewidth=2.5)
            ax.plot(r_Lo, z_Lo, color=color_Lo, linewidth=3.5)
            ax.plot([-r_b_val, r_b_val], [z_rim, z_rim], "ko", markersize=4.5)
            
        ax.axvline(0, color="gray", linestyle=":", linewidth=0.8, alpha=0.7)
        ax.set_title(f"(c{col_idx+1}) $v = {v_val:.2f}$ ($\\theta=1.0$)", fontweight="bold", pad=4)
        ax.set_aspect("equal")
        ax.set_xlim(-6.0, 6.0)
        ax.set_ylim(-6.0, 6.0)
        ax.set_xlabel(r"Radial $r$ ($\mu\mathrm{m}$)")
        ax.set_ylabel(r"Axial $z$ ($\mu\mathrm{m}$)")
        ax.grid(True, linestyle=":", alpha=0.5)

    plt.suptitle("HermesLoouter: Single-Domain Parameter Sweeps (Protein on Outer Leaflet of $L_o$)", fontsize=14, fontweight="bold", y=0.98)
    fig.savefig(out_fig, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved figure: {out_fig}")

if __name__ == "__main__":
    plot_param_sweeps()
