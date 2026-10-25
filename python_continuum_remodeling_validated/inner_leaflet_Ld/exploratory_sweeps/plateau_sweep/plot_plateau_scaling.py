"""
Plotting Module for Complete Plateau Scaling (D_Lo = 2.0 um to 14.0 um).
Renders 3D cutaway sequence and the complete asymptotic plateau curve up to the R_ves limit.
"""

import os
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from mpl_toolkits.mplot3d import Axes3D

def generate_plateau_scaling_figure():
    plateau_dir = os.path.join("simulation_results", "plateau_sweep")
    large_dir = os.path.join("simulation_results", "large_domains")
    size_dir = os.path.join("simulation_results", "size_sweep")
    
    # All tested diameters
    all_diameters = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0]
    
    # Read depths from trace files
    depth_dict = {}
    
    # 2, 3, 4, 5 um
    trace_size = os.path.join(size_dir, "trace_size_sweep.csv")
    if os.path.exists(trace_size):
        with open(trace_size, "r") as f:
            for row in csv.DictReader(f):
                if float(row["theta"]) == 1.0 and row["branch"] == "inward":
                    depth_dict[float(row["D_Lo"])] = abs(float(row["measured_dimple_depth"]))
                    
    # 6, 7, 8 um
    trace_large = os.path.join(large_dir, "trace_large_domain_sweep.csv")
    if os.path.exists(trace_large):
        with open(trace_large, "r") as f:
            for row in csv.DictReader(f):
                if float(row["theta"]) == 1.0 and row["branch"] == "inward":
                    depth_dict[float(row["D_Lo"])] = abs(float(row["measured_dimple_depth"]))
                    
    # 8 to 14 um
    trace_plat = os.path.join(plateau_dir, "trace_plateau_sweep.csv")
    if os.path.exists(trace_plat):
        with open(trace_plat, "r") as f:
            for row in csv.DictReader(f):
                if float(row["theta"]) == 1.0 and row["branch"] == "inward":
                    depth_dict[float(row["D_Lo"])] = abs(float(row["measured_dimple_depth"]))
                    
    depths_ordered = [depth_dict[D] for D in all_diameters]
    
    plt.rcParams.update({
        'font.size': 10,
        'axes.labelsize': 11,
        'axes.titlesize': 12,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 9,
        'figure.dpi': 300
    })
    
    fig = plt.figure(figsize=(18, 12))
    gs = GridSpec(2, 4, figure=fig, height_ratios=[1.2, 1.0], wspace=0.25, hspace=0.28)
    
    phi_cut = np.linspace(0, 1.5*np.pi, 50)
    
    # -------------------------------------------------------------
    # Top Row: 3D Cutaway Progression (2um, 6um, 10um, 14um)
    # -------------------------------------------------------------
    rep_diameters = [2.0, 6.0, 10.0, 14.0]
    rep_colors = {2.0: '#33a02c', 6.0: '#8c564b', 10.0: '#d95f02', 14.0: '#7570b3'}
    
    for i, D in enumerate(rep_diameters):
        ax_3d = fig.add_subplot(gs[0, i], projection='3d')
        
        if D <= 5.0:
            p_path = os.path.join(size_dir, f"profile_D{D:.1f}um_theta1.0_inward.csv")
        elif D <= 8.0:
            p_path = os.path.join(large_dir, f"profile_large_D{D:.1f}um_theta1.0_inward.csv")
        else:
            p_path = os.path.join(plateau_dir, f"profile_plateau_D{D:.1f}um_theta1.0_inward.csv")
            
        data = np.genfromtxt(p_path, delimiter=",", skip_header=1)
        s, psi, r, z, H, n_Lo_val = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5]
        n_Lo = int(n_Lo_val[0])
        
        z_max = np.max(z)
        z_realigned = z_max - z - (z_max / 2.0)
        r_realigned = r
        
        r_Lo = r_realigned[:n_Lo]
        z_Lo = z_realigned[:n_Lo]
        R_Lo_cut, Phi_Lo_cut = np.meshgrid(r_Lo, phi_cut)
        Z_Lo_cut = np.tile(z_Lo, (len(phi_cut), 1))
        X_Lo_cut = R_Lo_cut * np.cos(Phi_Lo_cut)
        Y_Lo_cut = R_Lo_cut * np.sin(Phi_Lo_cut)
        
        r_Ld = r_realigned[n_Lo-1:]
        z_Ld = z_realigned[n_Lo-1:]
        R_Ld_cut, Phi_Ld_cut = np.meshgrid(r_Ld, phi_cut)
        Z_Ld_cut = np.tile(z_Ld, (len(phi_cut), 1))
        X_Ld_cut = R_Ld_cut * np.cos(Phi_Ld_cut)
        Y_Ld_cut = R_Ld_cut * np.sin(Phi_Ld_cut)
        
        ax_3d.plot_surface(X_Lo_cut, Y_Lo_cut, Z_Lo_cut, color=rep_colors[D], alpha=0.95, edgecolor='none')
        ax_3d.plot_surface(X_Ld_cut, Y_Ld_cut, Z_Ld_cut, color='#4575b4', alpha=0.40, edgecolor='none')
        
        ax_3d.plot(r_realigned, np.zeros_like(r_realigned), z_realigned, color='black', linewidth=1.5)
        ax_3d.plot(np.zeros_like(r_realigned), -r_realigned, z_realigned, color='black', linewidth=1.5)
        
        regime_label = "RBC Dimple" if D <= 3 else ("Stomatocyte" if D <= 8 else "Vesosome Limit")
        depth_val = depth_dict[D]
        ax_3d.set_title(f"\\textbf{{{chr(97+i)}}} $D_{{L_o}} = {D:.0f}\\,\\mu\\mathrm{{m}}$ ({regime_label})\nDepth: ${depth_val:.2f}\\,\\mu\\mathrm{{m}}$", pad=5)
        ax_3d.view_init(elev=25, azim=40)
        ax_3d.set_box_aspect([1, 1, 1])
        ax_3d.set_xlabel('$X$')
        ax_3d.set_ylabel('$Y$')
        ax_3d.set_zlabel('$Z$')
        
    # -------------------------------------------------------------
    # Bottom Row Left: 2D Cross-Sections Comparison (Spans 2 cols)
    # -------------------------------------------------------------
    ax_prof = fig.add_subplot(gs[1, 0:2])
    
    colormap = plt.cm.viridis(np.linspace(0, 0.9, len(all_diameters)))
    
    for D, color in zip(all_diameters, colormap):
        if D <= 5.0:
            p_path = os.path.join(size_dir, f"profile_D{D:.1f}um_theta1.0_inward.csv")
        elif D <= 8.0:
            p_path = os.path.join(large_dir, f"profile_large_D{D:.1f}um_theta1.0_inward.csv")
        else:
            p_path = os.path.join(plateau_dir, f"profile_plateau_D{D:.1f}um_theta1.0_inward.csv")
            
        data = np.genfromtxt(p_path, delimiter=",", skip_header=1)
        s, psi, r, z, H, n_Lo_val = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5]
        n_Lo = int(n_Lo_val[0])
        
        z_max = np.max(z)
        z_realigned = z_max - z - (z_max / 2.0)
        r_realigned = r
        
        r_Lo_symm = np.concatenate([-r_realigned[:n_Lo][::-1], r_realigned[:n_Lo]])
        z_Lo_symm = np.concatenate([z_realigned[:n_Lo][::-1], z_realigned[:n_Lo]])
        
        ax_prof.plot(r_Lo_symm, z_Lo_symm, color=color, linewidth=2.0,
                     label=f"${D:.0f}\\,\\mu\\mathrm{{m}}$ (${depth_dict[D]:.2f}\\,\\mu\\mathrm{{m}}$)")
        
    ax_prof.set_xlabel(r'Radial Coordinate $r$ ($\mu\mathrm{m}$)')
    ax_prof.set_ylabel(r'Axial Coordinate $z$ ($\mu\mathrm{m}$)')
    ax_prof.set_title(r'\textbf{(e)} Morphological Sequence ($D_{L_o} = 2\,\mu\mathrm{m} \to 14\,\mu\mathrm{m}$)')
    ax_prof.axis('equal')
    ax_prof.grid(True, linestyle=':', alpha=0.6)
    ax_prof.legend(frameon=True, facecolor='white', framealpha=0.9, loc='lower right', ncol=3, fontsize=8)
    
    # -------------------------------------------------------------
    # Bottom Row Right: Asymptotic Plateau Scaling Curve (Spans 2 cols)
    # -------------------------------------------------------------
    ax_scale = fig.add_subplot(gs[1, 2:4])
    
    ax_scale.plot(all_diameters, depths_ordered, 'o-', color='#d73027', linewidth=2.8, markersize=8, label='Simulated Invagination Depth')
    
    # Geometric plateau limit: R_ves = 5.0 um
    ax_scale.axhline(5.0, color='#000000', linestyle='--', linewidth=2.0, label=r'Geometric Asymptotic Limit ($R_{\mathrm{ves}} = 5.0\,\mu\mathrm{m}$)')
    ax_scale.axhline(0.50, color='#2b5c8f', linestyle=':', linewidth=1.5, label=r'Dimple Criterion ($0.50\,\mu\mathrm{m}$)')
    
    # Shaded plateau saturation band
    ax_scale.axhspan(4.5, 5.0, color='#fee08b', alpha=0.35, label='Vesosome Plateau Saturation Regime')
    
    # Annotate plateau points
    for D, depth in zip(all_diameters[::2], depths_ordered[::2]):
        ax_scale.annotate(
            f"${depth:.2f}\\,\\mu\\mathrm{{m}}$",
            xy=(D, depth), xytext=(D - 0.35, depth + 0.18),
            fontweight='bold', fontsize=8.5, color='#a50026'
        )
        
    ax_scale.set_xlabel(r'$L_o$ Domain Diameter $D_{L_o}$ ($\mu\mathrm{m}$)')
    ax_scale.set_ylabel(r'Inward Invagination Depth ($\mu\mathrm{m}$)')
    ax_scale.set_title(r'\textbf{(f)} Asymptotic Invagination Plateau ($D_{L_o} = 2 \to 14\,\mu\mathrm{m}$)')
    ax_scale.set_ylim(0.2, 5.5)
    ax_scale.set_xlim(1.5, 14.5)
    ax_scale.grid(True, linestyle=':', alpha=0.6)
    ax_scale.legend(frameon=True, facecolor='white', framealpha=0.9, loc='lower right', fontsize=8.5)
    
    output_fig = os.path.join(plateau_dir, "3d_plateau_invagination_scaling.png")
    plt.tight_layout()
    plt.savefig(output_fig, dpi=300)
    plt.close()
    print(f"Plateau scaling figure saved to: {output_fig}")

if __name__ == "__main__":
    generate_plateau_scaling_figure()
