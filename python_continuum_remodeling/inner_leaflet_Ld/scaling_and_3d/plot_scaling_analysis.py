"""
Plotting and Visualization Module for Multi-Size Vesicle Scaling Analysis.
Generates publication-quality comparison figures demonstrating whether domain bending
scales linearly with vesicle and domain size.
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# Set publication-style formatting
plt.rcParams.update({
    'font.size': 11,
    'font.family': 'sans-serif',
    'axes.labelsize': 12,
    'axes.titlesize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.titlesize': 14,
    'lines.linewidth': 2.0,
    'axes.grid': True,
    'grid.alpha': 0.3,
    'grid.linestyle': '--'
})

def plot_scaling_summary(csv_path="simulation_results_prl/scaling_analysis/trace_scaling_sweep.csv",
                         out_dir="simulation_results_prl/scaling_analysis"):
    os.makedirs(out_dir, exist_ok=True)
    df = pd.read_csv(csv_path)
    
    colors = {
        "5um": "#2ca02c",    # Green
        "10um": "#1f77b4",   # Blue (Reference)
        "20um": "#d62728"    # Red
    }
    markers = {"5um": "o", "10um": "s", "20um": "^"}
    labels = {
        "5um": r"$D_{\rm ves}=5\,\mu{\rm m},\ D_{\rm Lo}=1.5\,\mu{\rm m}$",
        "10um": r"$D_{\rm ves}=10\,\mu{\rm m},\ D_{\rm Lo}=3.0\,\mu{\rm m}$ (Ref)",
        "20um": r"$D_{\rm ves}=20\,\mu{\rm m},\ D_{\rm Lo}=6.0\,\mu{\rm m}$"
    }

    fig = plt.figure(figsize=(16, 12))
    gs = GridSpec(2, 2, figure=fig, hspace=0.3, wspace=0.25)
    
    # -------------------------------------------------------------
    # Panel 1: Absolute Dimple Depth vs Coverage theta
    # -------------------------------------------------------------
    ax1 = fig.add_subplot(gs[0, 0])
    for s_name, group in df.groupby("size_name", sort=False):
        ax1.plot(group["theta"], -group["dimple_depth_um"], 
                 color=colors[s_name], marker=markers[s_name], 
                 label=labels[s_name], markersize=6)
    
    ax1.set_xlabel(r"Streptavidin Coverage $\theta$ on Ld Inner Leaflet")
    ax1.set_ylabel(r"Absolute Inward Dimple Depth $|\delta z|$ ($\mu$m)")
    ax1.set_title(r"{\bf (A)} Absolute Invagination Depth vs. Coverage", fontweight='bold')
    ax1.legend(loc="upper left", frameon=True)
    ax1.set_xlim(-0.02, 1.02)
    
    # -------------------------------------------------------------
    # Panel 2: Normalized Dimple Depth (delta_z / D_Lo) vs Coverage theta
    # -------------------------------------------------------------
    ax2 = fig.add_subplot(gs[0, 1])
    for s_name, group in df.groupby("size_name", sort=False):
        ax2.plot(group["theta"], -group["norm_dimple_depth"], 
                 color=colors[s_name], marker=markers[s_name], 
                 label=labels[s_name], markersize=6)
        
    ax2.set_xlabel(r"Streptavidin Coverage $\theta$ on Ld Inner Leaflet")
    ax2.set_ylabel(r"Normalized Depth $|\delta z| / D_{\rm Lo}$")
    ax2.set_title(r"{\bf (B)} Geometric Collapse Test (Linear Scaling Check)", fontweight='bold')
    ax2.annotate("If scaling were linear,\nall 3 curves would collapse!", 
                 xy=(0.5, 0.12), xytext=(0.15, 0.20),
                 arrowprops=dict(arrowstyle="->", color="black", lw=1.2),
                 bbox=dict(boxstyle="round,pad=0.3", fc="yellow", alpha=0.3, ec="orange"))
    ax2.legend(loc="upper left", frameon=True)
    ax2.set_xlim(-0.02, 1.02)

    # -------------------------------------------------------------
    # Panel 3: Normalized Domain Neck Aperture (2*r_b / D_Lo)
    # -------------------------------------------------------------
    ax3 = fig.add_subplot(gs[1, 0])
    for s_name, group in df.groupby("size_name", sort=False):
        ax3.plot(group["theta"], group["norm_neck_aperture"], 
                 color=colors[s_name], marker=markers[s_name], 
                 label=labels[s_name], markersize=6)
        
    ax3.set_xlabel(r"Streptavidin Coverage $\theta$ on Ld Inner Leaflet")
    ax3.set_ylabel(r"Normalized Neck Aperture $2 r_b / D_{\rm Lo}$")
    ax3.set_title(r"{\bf (C)} Domain Boundary Constriction ($2r_b / D_{\rm Lo}$)", fontweight='bold')
    ax3.legend(loc="lower left", frameon=True)
    ax3.set_xlim(-0.02, 1.02)
    
    # -------------------------------------------------------------
    # Panel 4: Energy Scaling (Total & Line Tension relative)
    # -------------------------------------------------------------
    ax4 = fig.add_subplot(gs[1, 1])
    for s_name, group in df.groupby("size_name", sort=False):
        # Line tension energy fraction
        e_frac_line = group["E_line"] / group["total_energy_kBT"] * 100.0
        ax4.plot(group["theta"], e_frac_line, 
                 color=colors[s_name], marker=markers[s_name], 
                 label=f"{s_name} (Line E %)", markersize=6)
                 
    ax4.set_xlabel(r"Streptavidin Coverage $\theta$ on Ld Inner Leaflet")
    ax4.set_ylabel(r"Line Tension Energy Contribution (\%)")
    ax4.set_title(r"{\bf (D)} Line Tension Fraction $E_{\rm line} / E_{\rm tot}$", fontweight='bold')
    ax4.legend(loc="upper left", frameon=True)
    ax4.set_xlim(-0.02, 1.02)
    
    plt.tight_layout()
    out_fig = os.path.join(out_dir, "fig_scaling_comprehensive_analysis.png")
    plt.savefig(out_fig, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved comprehensive scaling figure to: {out_fig}")


def plot_profile_overlays(out_dir="simulation_results_prl/scaling_analysis"):
    """
    Overlays normalized cross-sectional shapes (r/R_ves, z/R_ves) for theta=0.0, 0.6, 1.0
    """
    thetas_to_plot = [0.00, 0.60, 1.00]
    sizes = [
        {"name": "5um", "D_ves": 5.0, "D_Lo": 1.5, "color": "#2ca02c", "style": "--"},
        {"name": "10um", "D_ves": 10.0, "D_Lo": 3.0, "color": "#1f77b4", "style": "-"},
        {"name": "20um", "D_ves": 20.0, "D_Lo": 6.0, "color": "#d62728", "style": "-."}
    ]
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)
    
    for ax_idx, th in enumerate(thetas_to_plot):
        ax = axes[ax_idx]
        ax.set_aspect('equal')
        
        for s in sizes:
            p_file = f"prof_scale_Dves_{s['D_ves']:.1f}_DLo_{s['D_Lo']:.1f}_theta_{th:.2f}.csv"
            p_path = os.path.join(out_dir, "profiles", p_file)
            if not os.path.exists(p_path):
                continue
            df_p = pd.read_csv(p_path)
            
            # Plot right and left axisymmetric halves in normalized coordinates
            r_norm = df_p["r_norm"].values
            z_norm = df_p["z_norm"].values
            
            # Lo domain
            mask_Lo = (df_p["domain"] == "Lo")
            # Ld domain
            mask_Ld = (df_p["domain"] == "Ld")
            
            # Plot Ld
            ax.plot(r_norm[mask_Ld], z_norm[mask_Ld], color=s["color"], linestyle=s["style"], lw=2.0)
            ax.plot(-r_norm[mask_Ld], z_norm[mask_Ld], color=s["color"], linestyle=s["style"], lw=2.0)
            
            # Plot Lo with thicker highlight
            ax.plot(r_norm[mask_Lo], z_norm[mask_Lo], color=s["color"], lw=4.0, 
                    label=f"{s['name']} (Lo)" if ax_idx==0 else None)
            ax.plot(-r_norm[mask_Lo], z_norm[mask_Lo], color=s["color"], lw=4.0)

        ax.set_title(rf"$\theta = {th:.2f}$", fontsize=13, fontweight='bold')
        ax.set_xlabel(r"Normalized Radius $r / R_{\rm ves}$")
        if ax_idx == 0:
            ax.set_ylabel(r"Normalized Axial Position $z / R_{\rm ves}$")
            ax.legend(loc="lower center", frameon=True)
            
        ax.set_xlim(-1.15, 1.15)
        ax.set_ylim(-0.1, 2.1)
        ax.grid(True, linestyle="--", alpha=0.3)
        
    plt.suptitle("Normalized Vesicle Shape Overlays Across Sizes ($D_{\\rm Lo}/D_{\\rm ves} = 0.30$)", 
                 fontsize=15, fontweight='bold', y=1.02)
    plt.tight_layout()
    out_overlay = os.path.join(out_dir, "fig_scaling_profiles_overlay.png")
    plt.savefig(out_overlay, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved shape overlay figure to: {out_overlay}")

if __name__ == "__main__":
    plot_scaling_summary()
    plot_profile_overlays()
