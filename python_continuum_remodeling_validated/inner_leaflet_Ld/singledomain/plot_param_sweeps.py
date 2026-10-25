"""
Plotting Module for Multi-Parameter Sweeps (PRL Standard).
Generates 3-panel publication figures comparing:
- Panel A: Line Tension Sweep (gamma in {0.5, 1.5, 3.0, 5.0} pN)
- Panel B: Domain Size Sweep (D_Lo in {1.0, 3.0, 5.0, 7.0} um)
- Panel C: Deflation Factor Sweep (v in {0.90, 0.95, 0.98, 1.00})
Each panel plots symmetric vertical cross-sections (Lo in RED, Ld in GREEN).
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
        results_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "param_sweeps")
    profiles_dir = os.path.join(results_dir, "profiles")
    out_fig = os.path.join(results_dir, "fig_param_sweeps_comprehensive.png")
    
    color_Lo = "#d73027"  # Deep red for Lo
    color_Ld = "#1a9850"  # Forest green for Ld
    
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
    
    # -------------------------------------------------------------
    # Row 1: Line Tension Sweep (gamma = 0.5, 1.5, 3.0, 5.0 pN)
    # -------------------------------------------------------------
    gamma_vals = [0.5, 1.5, 3.0, 5.0]
    for col_idx, gamma in enumerate(gamma_vals):
        ax = fig.add_subplot(gs[0, col_idx])
        prof_path = _find_profile(profiles_dir, [
            f"prof_gamma_{gamma:.1f}_theta_1.00.csv",
            f"prof_gamma_{gamma:.1f}pN_theta_1.00.csv"
        ])
        
        if prof_path:
            df = pd.read_csv(prof_path)
            r = df["r_um"].values
            z = df["z_um"].values
            dom = df["domain"].values
            
            mask_lo = (dom == "Lo")
            n_Lo = np.sum(mask_lo)
            z_mid = 0.5 * (np.max(z) + np.min(z))
            z_aligned = z - z_mid
            
            z_apex = z_aligned[0]
            z_rim = z_aligned[n_Lo - 1]
            r_b_val = r[n_Lo - 1]
            depth_val = z_apex - z_rim
            
            r_Lo = np.concatenate([-r[:n_Lo][::-1], r[:n_Lo]])
            z_Lo = np.concatenate([z_aligned[:n_Lo][::-1], z_aligned[:n_Lo]])
            
            r_Ld_l = -r[n_Lo-1:][::-1]
            z_Ld_l = z_aligned[n_Lo-1:][::-1]
            r_Ld_r = r[n_Lo-1:]
            z_Ld_r = z_aligned[n_Lo-1:]
            
            ax.plot(r_Ld_l, z_Ld_l, color=color_Ld, linewidth=2.5, label=r"$L_d$ Matrix" if col_idx == 0 else "")
            ax.plot(r_Ld_r, z_Ld_r, color=color_Ld, linewidth=2.5)
            ax.plot(r_Lo, z_Lo, color=color_Lo, linewidth=3.5, label=r"$L_o$ Domain" if col_idx == 0 else "")
            ax.plot([-r_b_val, r_b_val], [z_rim, z_rim], "ko", markersize=4.5)
            
            depth_text = f"Inward: {abs(depth_val):.2f} $\\mu\\mathrm{{m}}$" if depth_val <= 0 else f"Outward: +{depth_val:.2f} $\\mu\\mathrm{{m}}$"
            ax.annotate(
                f"$\\gamma = {gamma:.1f}\\,\\mathrm{{pN}}$\n{depth_text}\n$2r_b = {2*r_b_val:.2f}\\,\\mu\\mathrm{{m}}$",
                xy=(0, z_apex), xytext=(0, z_apex - 1.5),
                fontsize=8, ha="center", fontweight="bold", color="#7f0000" if depth_val <= 0 else "#08519c",
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.9, edgecolor="#cccccc")
            )
            
        ax.axvline(0, color="gray", linestyle=":", linewidth=0.8, alpha=0.7)
        ax.set_title(f"(a{col_idx+1}) $\\gamma = {gamma:.1f}\\,\\mathrm{{pN}}$ ($\\theta=1.0$)", fontweight="bold", pad=4)
        ax.set_aspect("equal")
        ax.set_xlim(-6.0, 6.0)
        ax.set_ylim(-6.0, 6.0)
        ax.set_ylabel(r"Axial $z$ ($\mu\mathrm{m}$)")
        ax.grid(True, linestyle=":", alpha=0.5)
        if col_idx == 0:
            ax.legend(frameon=True, facecolor="white", loc="lower right")

    # -------------------------------------------------------------
    # Row 2: Domain Size Sweep (D_Lo = 1.0, 3.0, 5.0, 7.0 um)
    # -------------------------------------------------------------
    dlo_vals = [1.0, 3.0, 5.0, 7.0]
    for col_idx, D_Lo in enumerate(dlo_vals):
        ax = fig.add_subplot(gs[1, col_idx])
        prof_path = _find_profile(profiles_dir, [
            f"prof_DLo_{D_Lo:.1f}_theta_1.00.csv",
            f"prof_DLo_{D_Lo:.1f}um_theta_1.00.csv"
        ])
        
        if prof_path:
            df = pd.read_csv(prof_path)
            r = df["r_um"].values
            z = df["z_um"].values
            dom = df["domain"].values
            
            mask_lo = (dom == "Lo")
            n_Lo = np.sum(mask_lo)
            z_mid = 0.5 * (np.max(z) + np.min(z))
            z_aligned = z - z_mid
            
            z_apex = z_aligned[0]
            z_rim = z_aligned[n_Lo - 1]
            r_b_val = r[n_Lo - 1]
            depth_val = z_apex - z_rim
            
            r_Lo = np.concatenate([-r[:n_Lo][::-1], r[:n_Lo]])
            z_Lo = np.concatenate([z_aligned[:n_Lo][::-1], z_aligned[:n_Lo]])
            
            r_Ld_l = -r[n_Lo-1:][::-1]
            z_Ld_l = z_aligned[n_Lo-1:][::-1]
            r_Ld_r = r[n_Lo-1:]
            z_Ld_r = z_aligned[n_Lo-1:]
            
            ax.plot(r_Ld_l, z_Ld_l, color=color_Ld, linewidth=2.5)
            ax.plot(r_Ld_r, z_Ld_r, color=color_Ld, linewidth=2.5)
            ax.plot(r_Lo, z_Lo, color=color_Lo, linewidth=3.5)
            ax.plot([-r_b_val, r_b_val], [z_rim, z_rim], "ko", markersize=4.5)
            
            depth_text = f"Inward: {abs(depth_val):.2f} $\\mu\\mathrm{{m}}$" if depth_val <= 0 else f"Outward: +{depth_val:.2f} $\\mu\\mathrm{{m}}$"
            ax.annotate(
                f"$D_{{Lo}} = {D_Lo:.1f}\\,\\mu\\mathrm{{m}}$\n{depth_text}\n$2r_b = {2*r_b_val:.2f}\\,\\mu\\mathrm{{m}}$",
                xy=(0, z_apex), xytext=(0, z_apex - 1.5),
                fontsize=8, ha="center", fontweight="bold", color="#7f0000" if depth_val <= 0 else "#08519c",
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.9, edgecolor="#cccccc")
            )
            
        ax.axvline(0, color="gray", linestyle=":", linewidth=0.8, alpha=0.7)
        ax.set_title(f"(b{col_idx+1}) $D_{{Lo}} = {D_Lo:.1f}\\,\\mu\\mathrm{{m}}$ ($\\theta=1.0$)", fontweight="bold", pad=4)
        ax.set_aspect("equal")
        ax.set_xlim(-6.0, 6.0)
        ax.set_ylim(-6.0, 6.0)
        ax.set_ylabel(r"Axial $z$ ($\mu\mathrm{m}$)")
        ax.grid(True, linestyle=":", alpha=0.5)

    # -------------------------------------------------------------
    # Row 3: Deflation Sweep (v = 0.90, 0.95, 0.98, 1.00)
    # -------------------------------------------------------------
    v_vals = [0.90, 0.95, 0.98, 1.00]
    for col_idx, v in enumerate(v_vals):
        ax = fig.add_subplot(gs[2, col_idx])
        prof_path = _find_profile(profiles_dir, [
            f"prof_v_{v:.2f}_theta_1.00.csv"
        ])
        
        if prof_path:
            df = pd.read_csv(prof_path)
            r = df["r_um"].values
            z = df["z_um"].values
            dom = df["domain"].values
            
            mask_lo = (dom == "Lo")
            n_Lo = np.sum(mask_lo)
            z_mid = 0.5 * (np.max(z) + np.min(z))
            z_aligned = z - z_mid
            
            z_apex = z_aligned[0]
            z_rim = z_aligned[n_Lo - 1]
            r_b_val = r[n_Lo - 1]
            depth_val = z_apex - z_rim
            
            r_Lo = np.concatenate([-r[:n_Lo][::-1], r[:n_Lo]])
            z_Lo = np.concatenate([z_aligned[:n_Lo][::-1], z_aligned[:n_Lo]])
            
            r_Ld_l = -r[n_Lo-1:][::-1]
            z_Ld_l = z_aligned[n_Lo-1:][::-1]
            r_Ld_r = r[n_Lo-1:]
            z_Ld_r = z_aligned[n_Lo-1:]
            
            ax.plot(r_Ld_l, z_Ld_l, color=color_Ld, linewidth=2.5)
            ax.plot(r_Ld_r, z_Ld_r, color=color_Ld, linewidth=2.5)
            ax.plot(r_Lo, z_Lo, color=color_Lo, linewidth=3.5)
            ax.plot([-r_b_val, r_b_val], [z_rim, z_rim], "ko", markersize=4.5)
            
            depth_text = f"Inward: {abs(depth_val):.2f} $\\mu\\mathrm{{m}}$" if depth_val <= 0 else f"Outward: +{depth_val:.2f} $\\mu\\mathrm{{m}}$"
            ax.annotate(
                f"$v = {v:.2f}$\n{depth_text}\n$2r_b = {2*r_b_val:.2f}\\,\\mu\\mathrm{{m}}$",
                xy=(0, z_apex), xytext=(0, z_apex - 1.5),
                fontsize=8, ha="center", fontweight="bold", color="#7f0000" if depth_val <= 0 else "#08519c",
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.9, edgecolor="#cccccc")
            )
            
        ax.axvline(0, color="gray", linestyle=":", linewidth=0.8, alpha=0.7)
        ax.set_title(f"(c{col_idx+1}) $v = {v:.2f}$ ($\\theta=1.0$)", fontweight="bold", pad=4)
        ax.set_aspect("equal")
        ax.set_xlim(-6.0, 6.0)
        ax.set_ylim(-6.0, 6.0)
        ax.set_xlabel(r"Radial $r$ ($\mu\mathrm{m}$)")
        ax.set_ylabel(r"Axial $z$ ($\mu\mathrm{m}$)")
        ax.grid(True, linestyle=":", alpha=0.5)

    plt.tight_layout()
    plt.savefig(out_fig, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Comprehensive parameter sweep figure saved to: {out_fig}")

if __name__ == "__main__":
    plot_param_sweeps()
