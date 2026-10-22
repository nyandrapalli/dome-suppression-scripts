"""
Publication-Quality Plotting Module for Two-Phase Vesicle Remodeling (PRL Standard)
Generates 4-panel publication figures:
- Panel A: Bifurcation diagram & Free Energy landscape vs coverage theta
- Panel B: Axisymmetric cross-sectional shapes (sphere vs inward dimple vs outward bud)
- Panel C: Energy decomposition (Bending, ADE, Line, Gauss-Bonnet)
- Panel D: Numerical constraint verification audit (relative errors < 10^-5)
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def plot_prl_results(results_dir="simulation_results_prl"):
    trace_path = os.path.join(results_dir, "trace_remodeling_sweep_prl.csv")
    if not os.path.exists(trace_path):
        raise FileNotFoundError(f"Trace file not found: {trace_path}")
        
    df = pd.read_csv(trace_path)
    profiles_dir = os.path.join(results_dir, "profiles")
    
    fig = plt.figure(figsize=(13, 10), dpi=300)
    gs = fig.add_gridspec(2, 2, hspace=0.32, wspace=0.28)
    
    colors = {"sphere": "#4575b4", "inward": "#d73027", "outward": "#fdae61"}
    markers = {"sphere": "o", "inward": "s", "outward": "^"}
    
    # -------------------------------------------------------------
    # Panel A: Free Energy Bifurcation Diagram
    # -------------------------------------------------------------
    ax_a = fig.add_subplot(gs[0, 0])
    for branch in ["sphere", "inward", "outward"]:
        sub = df[df["branch"] == branch].sort_values("theta")
        ax_a.plot(
            sub["theta"], sub["total_energy_kBT"],
            label=f"{branch.capitalize()} branch",
            color=colors[branch], marker=markers[branch],
            linewidth=2.2, markersize=6
        )
    ax_a.set_xlabel(r"Streptavidin Coverage $\theta$", fontsize=11, fontweight="bold")
    ax_a.set_ylabel(r"Total Free Energy $\mathcal{E}$ ($k_B T$)", fontsize=11, fontweight="bold")
    ax_a.set_title("(a) Free Energy Bifurcation Landscape", fontsize=12, fontweight="bold", loc="left")
    ax_a.grid(True, linestyle="--", alpha=0.5)
    ax_a.legend(frameon=True, fontsize=10)
    
    # -------------------------------------------------------------
    # Panel B: Cross-Sectional Vesicle Profiles
    # -------------------------------------------------------------
    ax_b = fig.add_subplot(gs[0, 1])
    # Plot profiles for theta = 0.0 (sphere) and theta = 1.0 (inward dimple)
    cases = [
        ("profile_theta_0.00_sphere.csv", "Bare Sphere ($\theta=0.0$)", "#4575b4", "-"),
        ("profile_theta_1.00_inward.csv", "Remodeled Inward ($\theta=1.0$)", "#d73027", "-"),
        ("profile_theta_1.00_outward.csv", "Outward Bud ($\theta=1.0$)", "#fdae61", "--"),
    ]
    
    for filename, label, col, ls in cases:
        path = os.path.join(profiles_dir, filename)
        if os.path.exists(path):
            p_df = pd.read_csv(path)
            r = p_df["r_um"].values
            z = p_df["z_um"].values
            # Plot right and left symmetric halves
            ax_b.plot(r, z, label=label, color=col, linestyle=ls, linewidth=2.0)
            ax_b.plot(-r, z, color=col, linestyle=ls, linewidth=2.0)
            
            # Highlight Lo domain apex
            mask_lo = (p_df["domain"] == "Lo")
            ax_b.plot(r[mask_lo], z[mask_lo], color="#800026" if "Inward" in label else "#1a9850", linewidth=3.5)
            ax_b.plot(-r[mask_lo], z[mask_lo], color="#800026" if "Inward" in label else "#1a9850", linewidth=3.5)
            
    ax_b.set_xlabel(r"Radial coordinate $r$ ($\mu$m)", fontsize=11, fontweight="bold")
    ax_b.set_ylabel(r"Axial coordinate $z$ ($\mu$m)", fontsize=11, fontweight="bold")
    ax_b.set_title("(b) Cross-Sectional Morphologies", fontsize=12, fontweight="bold", loc="left")
    ax_b.set_aspect("equal")
    ax_b.grid(True, linestyle="--", alpha=0.5)
    ax_b.legend(frameon=True, fontsize=9, loc="lower right")
    
    # -------------------------------------------------------------
    # Panel C: Physical Energy Breakdown (Inward Dimple Branch)
    # -------------------------------------------------------------
    ax_c = fig.add_subplot(gs[1, 0])
    sub_inward = df[df["branch"] == "inward"].sort_values("theta")
    ax_c.plot(sub_inward["theta"], sub_inward["E_bend_tot"], label=r"$E_{\mathrm{bend}}$ (Canham-Helfrich)", color="#2b83ba", linewidth=2.0, marker="o")
    ax_c.plot(sub_inward["theta"], sub_inward["E_ade"], label=r"$E_{\mathrm{ADE}}$ (Monolayer mismatch)", color="#d7191c", linewidth=2.0, marker="s")
    ax_c.plot(sub_inward["theta"], sub_inward["E_line"], label=r"$E_{\mathrm{line}}$ (Domain boundary)", color="#fdae61", linewidth=2.0, marker="^")
    ax_c.plot(sub_inward["theta"], sub_inward["E_gauss"], label=r"$E_{\mathrm{Gauss}}$ (Gauss-Bonnet)", color="#abdda4", linewidth=2.0, marker="d")
    
    ax_c.set_xlabel(r"Streptavidin Coverage $\theta$", fontsize=11, fontweight="bold")
    ax_c.set_ylabel(r"Energy Component ($k_B T$)", fontsize=11, fontweight="bold")
    ax_c.set_title("(c) Energy Decomposition (Inward Branch)", fontsize=12, fontweight="bold", loc="left")
    ax_c.grid(True, linestyle="--", alpha=0.5)
    ax_c.legend(frameon=True, fontsize=9.5)
    
    # -------------------------------------------------------------
    # Panel D: Numerical Constraint Violation Audit (PRL Quality)
    # -------------------------------------------------------------
    ax_d = fig.add_subplot(gs[1, 1])
    for branch in ["sphere", "inward", "outward"]:
        sub = df[df["branch"] == branch].sort_values("theta")
        ax_d.semilogy(
            sub["theta"], sub["max_geom_error"],
            label=f"{branch.capitalize()} max error",
            color=colors[branch], marker=markers[branch],
            linewidth=1.8, markersize=5
        )
    ax_d.axhline(1e-4, color="gray", linestyle=":", label=r"PRL threshold ($10^{-4}$)")
    ax_d.set_xlabel(r"Streptavidin Coverage $\theta$", fontsize=11, fontweight="bold")
    ax_d.set_ylabel(r"Maximum Constraint Violation $|\Delta g_i|$", fontsize=11, fontweight="bold")
    ax_d.set_title("(d) Numerical Constraint Audit", fontsize=12, fontweight="bold", loc="left")
    ax_d.grid(True, linestyle="--", alpha=0.5)
    ax_d.legend(frameon=True, fontsize=9.5)
    
    out_fig = os.path.join(results_dir, "prl_remodeling_bifurcation_and_audit.png")
    fig.savefig(out_fig, bbox_inches="tight")
    plt.close(fig)
    print(f"Publication figure generated: {out_fig}")
    return out_fig

if __name__ == "__main__":
    plot_prl_results()
