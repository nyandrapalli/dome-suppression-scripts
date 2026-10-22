"""
Plotting Module for Two-Patch Antipodal System.
Generates vertical central cross-sections (Lo in Red, Ld in Green) for D = 2um to 7um,
and renders 3D biconcave discocyte models.
"""

import os
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from mpl_toolkits.mplot3d import Axes3D

def generate_two_patch_figures():
    out_dir = os.path.join("simulation_results", "two_patch_system")
    os.makedirs(out_dir, exist_ok=True)
    
    diameters = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
    
    color_Lo = '#d73027'  # Red for Lo
    color_Ld = '#1a9850'  # Green for Ld
    
    plt.rcParams.update({
        'font.size': 10,
        'axes.labelsize': 11,
        'axes.titlesize': 12,
        'xtick.labelsize': 9.5,
        'ytick.labelsize': 9.5,
        'legend.fontsize': 9.5,
        'figure.dpi': 300
    })
    
    # -------------------------------------------------------------
    # Figure 1: Grid of 2D Vertical Central Cross-Sections
    # -------------------------------------------------------------
    fig, axes = plt.subplots(2, 3, figsize=(18, 12), sharex=True, sharey=True)
    axes_flat = axes.flatten()
    
    profiles_data = []
    
    for i, D in enumerate(diameters):
        ax = axes_flat[i]
        p_path = os.path.join(out_dir, f"profile_twopatch_D{D:.1f}um_theta1.0_inward.csv")
        
        data = np.genfromtxt(p_path, delimiter=",", skip_header=1)
        s, psi, r, z, H, n_Lo1_val, n_Lo2_val = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5], data[:,6]
        n_Lo1 = int(n_Lo1_val[0])
        n_Lo2 = int(n_Lo2_val[0])
        idx_b1 = n_Lo1 - 1
        idx_b2 = len(s) - n_Lo2
        
        # Center vertical coordinate
        z_max = np.max(z)
        z_realigned = z_max - z - (z_max / 2.0)
        r_realigned = r
        
        z_north_apex = z_realigned[0]
        z_north_rim = z_realigned[idx_b1]
        depth_north = abs(z_north_rim - z_north_apex)
        
        z_south_apex = z_realigned[-1]
        z_south_rim = z_realigned[idx_b2]
        depth_south = abs(z_south_rim - z_south_apex)
        
        r_b1 = r_realigned[idx_b1]
        r_b2 = r_realigned[idx_b2]
        
        profiles_data.append({
            'D': D,
            'r': r_realigned,
            'z': z_realigned,
            'idx_b1': idx_b1,
            'idx_b2': idx_b2,
            'depth_north': depth_north,
            'depth_south': depth_south
        })
        
        # North Lo symmetric
        r_Lo1_symm = np.concatenate([-r_realigned[:idx_b1+1][::-1], r_realigned[:idx_b1+1]])
        z_Lo1_symm = np.concatenate([z_realigned[:idx_b1+1][::-1], z_realigned[:idx_b1+1]])
        
        # South Lo symmetric
        r_Lo2_symm = np.concatenate([-r_realigned[idx_b2:][::-1], r_realigned[idx_b2:]])
        z_Lo2_symm = np.concatenate([z_realigned[idx_b2:][::-1], z_realigned[idx_b2:]])
        
        # Equatorial Ld belt symmetric (left & right sides)
        r_Ld_left = -r_realigned[idx_b1:idx_b2+1][::-1]
        z_Ld_left = z_realigned[idx_b1:idx_b2+1][::-1]
        r_Ld_right = r_realigned[idx_b1:idx_b2+1]
        z_Ld_right = z_realigned[idx_b1:idx_b2+1]
        
        # Plot Ld Belt in Green
        ax.plot(r_Ld_left, z_Ld_left, color=color_Ld, linewidth=3.2, label=r'Equatorial $L_d$ Belt (Green)')
        ax.plot(r_Ld_right, z_Ld_right, color=color_Ld, linewidth=3.2)
        
        # Plot North & South Lo in Red
        ax.plot(r_Lo1_symm, z_Lo1_symm, color=color_Lo, linewidth=3.8, label=r'Antipodal $L_o$ Caps (Red)')
        ax.plot(r_Lo2_symm, z_Lo2_symm, color=color_Lo, linewidth=3.8)
        
        # Boundary dots
        ax.plot([-r_b1, r_b1], [z_north_rim, z_north_rim], 'ko', markersize=5.0, zorder=5)
        ax.plot([-r_b2, r_b2], [z_south_rim, z_south_rim], 'ko', markersize=5.0, zorder=5)
        
        # Symmetry axis
        ax.axvline(0, color='gray', linestyle=':', linewidth=0.8, alpha=0.7)
        
        # Annotate depths
        ax.text(
            0, z_north_apex + 0.6,
            f"North Depth: {depth_north:.2f} $\\mu\\mathrm{{m}}$",
            fontsize=8.5, ha='center', fontweight='bold', color='#7f0000',
            bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.9, edgecolor='#cccccc')
        )
        ax.text(
            0, z_south_apex - 0.6,
            f"South Depth: {depth_south:.2f} $\\mu\\mathrm{{m}}$",
            fontsize=8.5, ha='center', fontweight='bold', color='#7f0000',
            bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.9, edgecolor='#cccccc')
        )
        
        ax.set_title(f"\\textbf{{({chr(97+i)})}} $D_{{L_o}} = {D:.0f}\\,\\mu\\mathrm{{m}}$ (Dual Dimples: {0.5*(depth_north+depth_south):.2f} $\\mu\\mathrm{{m}}$ avg)", pad=4)
        ax.set_aspect('equal')
        ax.set_xlim(-6.0, 6.0)
        ax.set_ylim(-6.0, 6.0)
        ax.grid(True, linestyle=':', alpha=0.5)
        
        if i == 0:
            ax.legend(frameon=True, facecolor='white', framealpha=0.9, loc='upper right', fontsize=8)
            
    for r_idx in range(2):
        axes[r_idx, 0].set_ylabel(r'Axial Coordinate $z$ ($\mu\mathrm{m}$)')
    for c_idx in range(3):
        axes[1, c_idx].set_xlabel(r'Radial Coordinate $r$ ($\mu\mathrm{m}$)')
        
    out_fig1 = os.path.join(out_dir, "two_patch_vertical_slices.png")
    plt.tight_layout()
    plt.savefig(out_fig1, dpi=300)
    plt.close()
    print(f"Two-patch vertical slices saved to: {out_fig1}")
    
    # -------------------------------------------------------------
    # Figure 2: 3D Biconcave Discocyte Renderings & Cutaways
    # -------------------------------------------------------------
    fig3d = plt.figure(figsize=(16, 8))
    gs3d = GridSpec(1, 2, figure=fig3d, wspace=0.18)
    
    # Select representative D = 4.0 um
    p_rep_path = os.path.join(out_dir, "profile_twopatch_D4.0um_theta1.0_inward.csv")
    data_rep = np.genfromtxt(p_rep_path, delimiter=",", skip_header=1)
    s_r, r_r, z_r, n1_r, n2_r = data_rep[:,0], data_rep[:,2], data_rep[:,3], int(data_rep[0,5]), int(data_rep[0,6])
    idx1_r = n1_r - 1
    idx2_r = len(s_r) - n2_r
    
    z_max_r = np.max(z_r)
    z_real_r = z_max_r - z_r - (z_max_r / 2.0)
    
    phi = np.linspace(0, 2*np.pi, 60)
    phi_cut = np.linspace(0, 1.5*np.pi, 50)
    
    # Panel A: 3D Full Surface
    ax_full = fig3d.add_subplot(gs3d[0, 0], projection='3d')
    
    # North Lo
    R1, P1 = np.meshgrid(r_r[:idx1_r+1], phi)
    Z1 = np.tile(z_real_r[:idx1_r+1], (len(phi), 1))
    ax_full.plot_surface(R1*np.cos(P1), R1*np.sin(P1), Z1, color=color_Lo, alpha=0.95, edgecolor='none')
    
    # Equatorial Ld
    R_Ld, P_Ld = np.meshgrid(r_r[idx1_r:idx2_r+1], phi)
    Z_Ld = np.tile(z_real_r[idx1_r:idx2_r+1], (len(phi), 1))
    ax_full.plot_surface(R_Ld*np.cos(P_Ld), R_Ld*np.sin(P_Ld), Z_Ld, color=color_Ld, alpha=0.45, edgecolor='none')
    
    # South Lo
    R2, P2 = np.meshgrid(r_r[idx2_r:], phi)
    Z2 = np.tile(z_real_r[idx2_r:], (len(phi), 1))
    ax_full.plot_surface(R2*np.cos(P2), R2*np.sin(P2), Z2, color=color_Lo, alpha=0.95, edgecolor='none')
    
    ax_full.set_title(r"\textbf{(a)} 3D Remodeled Biconcave Discocyte ($D_{L_o} = 4\,\mu\mathrm{m}$)", pad=8)
    ax_full.view_init(elev=20, azim=45)
    ax_full.set_box_aspect([1, 1, 1])
    ax_full.set_xlabel('$X$')
    ax_full.set_ylabel('$Y$')
    ax_full.set_zlabel('$Z$')
    
    # Panel B: 3D Cutaway View
    ax_cut = fig3d.add_subplot(gs3d[0, 1], projection='3d')
    
    R1_c, P1_c = np.meshgrid(r_r[:idx1_r+1], phi_cut)
    Z1_c = np.tile(z_real_r[:idx1_r+1], (len(phi_cut), 1))
    ax_cut.plot_surface(R1_c*np.cos(P1_c), R1_c*np.sin(P1_c), Z1_c, color=color_Lo, alpha=0.98, edgecolor='none')
    
    R_Ld_c, P_Ld_c = np.meshgrid(r_r[idx1_r:idx2_r+1], phi_cut)
    Z_Ld_c = np.tile(z_real_r[idx1_r:idx2_r+1], (len(phi_cut), 1))
    ax_cut.plot_surface(R_Ld_c*np.cos(P_Ld_c), R_Ld_c*np.sin(P_Ld_c), Z_Ld_c, color=color_Ld, alpha=0.45, edgecolor='none')
    
    R2_c, P2_c = np.meshgrid(r_r[idx2_r:], phi_cut)
    Z2_c = np.tile(z_real_r[idx2_r:], (len(phi_cut), 1))
    ax_cut.plot_surface(R2_c*np.cos(P2_c), R2_c*np.sin(P2_c), Z2_c, color=color_Lo, alpha=0.98, edgecolor='none')
    
    ax_cut.plot(r_r, np.zeros_like(r_r), z_real_r, color='black', linewidth=1.8)
    ax_cut.plot(np.zeros_like(r_r), -r_r, z_real_r, color='black', linewidth=1.8)
    
    ax_cut.set_title(r"\textbf{(b)} 3D Cutaway: Opposing Dual RBC Inward Cavities", pad=8)
    ax_cut.view_init(elev=20, azim=45)
    ax_cut.set_box_aspect([1, 1, 1])
    ax_cut.set_xlabel('$X$')
    ax_cut.set_ylabel('$Y$')
    ax_cut.set_zlabel('$Z$')
    
    out_fig2 = os.path.join(out_dir, "3d_two_patch_biconcave_proof.png")
    plt.tight_layout()
    plt.savefig(out_fig2, dpi=300)
    plt.close()
    print(f"Two-patch 3D figure saved to: {out_fig2}")

if __name__ == "__main__":
    generate_two_patch_figures()
