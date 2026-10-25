"""
Comparative Plotting Module for Large Lo Domain Sweep (D = 6.0um, 7.0um, 8.0um)
and Unified Scaling across all sizes (D = 2.0um to 8.0um).
"""

import os
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from mpl_toolkits.mplot3d import Axes3D

def generate_large_domain_comparison():
    large_results_dir = os.path.join("simulation_results", "large_domains")
    size_results_dir = os.path.join("simulation_results", "size_sweep")
    
    diameters_large = [6.0, 7.0, 8.0]
    colors_large = {6.0: '#8c564b', 7.0: '#e377c2', 8.0: '#17becf'}
    
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
    gs = GridSpec(2, 3, figure=fig, height_ratios=[1.2, 1.0], wspace=0.25, hspace=0.28)
    
    phi_cut = np.linspace(0, 1.5*np.pi, 50)
    
    depths_large = []
    
    # -------------------------------------------------------------
    # Top Row: 3D Cutaway Renderings for 6um, 7um, 8um at theta = 1.0
    # -------------------------------------------------------------
    for i, D in enumerate(diameters_large):
        ax_3d = fig.add_subplot(gs[0, i], projection='3d')
        p_fn = f"profile_large_D{D:.1f}um_theta1.0_inward.csv"
        p_path = os.path.join(large_results_dir, p_fn)
        
        data = np.genfromtxt(p_path, delimiter=",", skip_header=1)
        s, psi, r, z, H, n_Lo_val = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5]
        n_Lo = int(n_Lo_val[0])
        
        z_max = np.max(z)
        z_realigned = z_max - z - (z_max / 2.0)
        r_realigned = r
        
        z_apex = z_realigned[0]
        z_rim = z_realigned[n_Lo - 1]
        depth_val = z_rim - z_apex  # Positive depth into vesicle
        depths_large.append(depth_val)
        
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
        
        ax_3d.plot_surface(X_Lo_cut, Y_Lo_cut, Z_Lo_cut, color=colors_large[D], alpha=0.95, edgecolor='none')
        ax_3d.plot_surface(X_Ld_cut, Y_Ld_cut, Z_Ld_cut, color='#4575b4', alpha=0.40, edgecolor='none')
        
        # Cut edges
        ax_3d.plot(r_realigned, np.zeros_like(r_realigned), z_realigned, color='black', linewidth=1.5)
        ax_3d.plot(np.zeros_like(r_realigned), -r_realigned, z_realigned, color='black', linewidth=1.5)
        
        ax_3d.set_title(f"\\textbf{{{chr(97+i)}}} $D_{{L_o}} = {D:.0f}\\,\\mu\\mathrm{{m}}$ (Stomatocyte)\nDepth: ${depth_val:.3f}\\,\\mu\\mathrm{{m}}$", pad=5)
        ax_3d.view_init(elev=25, azim=40)
        ax_3d.set_box_aspect([1, 1, 1])
        ax_3d.set_xlabel('$X$')
        ax_3d.set_ylabel('$Y$')
        ax_3d.set_zlabel('$Z$')
        
    # -------------------------------------------------------------
    # Bottom Row Left: 2D Cross-Section Profiles Comparison (Spans 1.5 cols)
    # -------------------------------------------------------------
    ax_prof = fig.add_subplot(gs[1, 0:2])
    
    # Plot all sizes from 2 to 8 um
    all_diameters = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    all_colors = {
        2.0: '#33a02c', 3.0: '#e31a1c', 4.0: '#ff7f00', 5.0: '#6a3d9a',
        6.0: '#8c564b', 7.0: '#e377c2', 8.0: '#17becf'
    }
    
    for D in all_diameters:
        if D <= 5.0:
            p_path = os.path.join(size_results_dir, f"profile_D{D:.1f}um_theta1.0_inward.csv")
        else:
            p_path = os.path.join(large_results_dir, f"profile_large_D{D:.1f}um_theta1.0_inward.csv")
            
        data = np.genfromtxt(p_path, delimiter=",", skip_header=1)
        s, psi, r, z, H, n_Lo_val = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5]
        n_Lo = int(n_Lo_val[0])
        
        z_max = np.max(z)
        z_realigned = z_max - z - (z_max / 2.0)
        r_realigned = r
        
        r_Lo_symm = np.concatenate([-r_realigned[:n_Lo][::-1], r_realigned[:n_Lo]])
        z_Lo_symm = np.concatenate([z_realigned[:n_Lo][::-1], z_realigned[:n_Lo]])
        
        ax_prof.plot(r_Lo_symm, z_Lo_symm, color=all_colors[D], linewidth=2.5,
                     label=f"$D_{{L_o}} = {D:.0f}\\,\\mu\\mathrm{{m}}$ (Depth: ${z_realigned[n_Lo-1] - z_realigned[0]:.2f}\\,\\mu\\mathrm{{m}}$)")
        
    ax_prof.set_xlabel(r'Radial Coordinate $r$ ($\mu\mathrm{m}$)')
    ax_prof.set_ylabel(r'Axial Coordinate $z$ ($\mu\mathrm{m}$)')
    ax_prof.set_title(r'\textbf{(d)} Remodeled Profiles from Local Dimples ($2\,\mu\mathrm{m}$) to Deep Cavities ($8\,\mu\mathrm{m}$)')
    ax_prof.axis('equal')
    ax_prof.grid(True, linestyle=':', alpha=0.6)
    ax_prof.legend(frameon=True, facecolor='white', framealpha=0.9, loc='lower right', ncol=2)
    
    # -------------------------------------------------------------
    # Bottom Row Right: Full Spectrum Scaling Curve (2um to 8um)
    # -------------------------------------------------------------
    ax_scale = fig.add_subplot(gs[1, 2])
    
    # Collate all depths
    all_depths = [0.561, 0.837, 1.161, 1.493] + depths_large
    
    ax_scale.plot(all_diameters, all_depths, 'o-', color='#d73027', linewidth=2.5, markersize=8, label='Simulated Depth')
    ax_scale.axhline(0.50, color='#2b5c8f', linestyle='--', linewidth=2.0, label='Dimple Criterion ($0.50\\,\\mu\\mathrm{m}$)')
    
    for D, depth in zip(all_diameters, all_depths):
        ax_scale.annotate(
            f"${depth:.2f}\\,\\mu\\mathrm{{m}}$",
            xy=(D, depth), xytext=(D - 0.25, depth + 0.09),
            fontweight='bold', fontsize=8.5, color='#a50026'
        )
        
    ax_scale.set_xlabel(r'$L_o$ Domain Diameter $D_{L_o}$ ($\mu\mathrm{m}$)')
    ax_scale.set_ylabel(r'Inward Dimple / Invagination Depth ($\mu\mathrm{m}$)')
    ax_scale.set_title(r'\textbf{(e)} Comprehensive Inward Depth Scaling ($2 \to 8\,\mu\mathrm{m}$)')
    ax_scale.set_ylim(0.3, 2.7)
    ax_scale.grid(True, linestyle=':', alpha=0.6)
    ax_scale.legend(frameon=True, facecolor='white', framealpha=0.9, loc='upper left')
    
    output_fig = os.path.join(large_results_dir, "3d_large_domain_dimple_comparison.png")
    plt.tight_layout()
    plt.savefig(output_fig, dpi=300)
    plt.close()
    print(f"Large domain comparison figure saved to: {output_fig}")

if __name__ == "__main__":
    generate_large_domain_comparison()
