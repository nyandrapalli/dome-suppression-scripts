"""
Vertical Central Slice Generator for Remodeled Vesicles (D_Lo = 2.0 um to 14.0 um).
Renders the 2D central vertical cross-section for each vesicle with:
  - Lo domain in RED
  - Ld domain in GREEN
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

def generate_vertical_slices():
    out_dir = os.path.join("simulation_results", "vertical_slices")
    os.makedirs(out_dir, exist_ok=True)
    
    size_dir = os.path.join("simulation_results", "size_sweep")
    large_dir = os.path.join("simulation_results", "large_domains")
    plateau_dir = os.path.join("simulation_results", "plateau_sweep")
    
    diameters = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0]
    
    plt.rcParams.update({
        'font.size': 9,
        'axes.labelsize': 10,
        'axes.titlesize': 11,
        'xtick.labelsize': 8.5,
        'ytick.labelsize': 8.5,
        'legend.fontsize': 8.5,
        'figure.dpi': 300
    })
    
    # 4 rows x 4 cols grid (13 individual vesicles + 1 composite overlay panel)
    fig, axes = plt.subplots(4, 4, figsize=(20, 18), sharex=True, sharey=True)
    axes_flat = axes.flatten()
    
    color_Lo = '#d73027'  # Red for Lo
    color_Ld = '#1a9850'  # Green for Ld
    
    profiles_data = []
    
    for i, D in enumerate(diameters):
        ax = axes_flat[i]
        
        # Determine path
        if D <= 5.0:
            p_path = os.path.join(size_dir, f"profile_D{D:.1f}um_theta1.0_inward.csv")
        elif D <= 8.0:
            p_path = os.path.join(large_dir, f"profile_large_D{D:.1f}um_theta1.0_inward.csv")
        else:
            p_path = os.path.join(plateau_dir, f"profile_plateau_D{D:.1f}um_theta1.0_inward.csv")
            
        data = np.genfromtxt(p_path, delimiter=",", skip_header=1)
        s, psi, r, z, H, n_Lo_val = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5]
        n_Lo = int(n_Lo_val[0])
        
        # Center vertical alignment
        z_max = np.max(z)
        z_realigned = z_max - z - (z_max / 2.0)
        r_realigned = r
        
        z_apex = z_realigned[0]
        z_rim = z_realigned[n_Lo - 1]
        depth_val = abs(z_rim - z_apex)
        r_b_val = r_realigned[n_Lo - 1]
        
        profiles_data.append({
            'D': D,
            'r': r_realigned,
            'z': z_realigned,
            'n_Lo': n_Lo,
            'depth': depth_val
        })
        
        # Symmetric vertical cross-section profiles
        r_Lo_symm = np.concatenate([-r_realigned[:n_Lo][::-1], r_realigned[:n_Lo]])
        z_Lo_symm = np.concatenate([z_realigned[:n_Lo][::-1], z_realigned[:n_Lo]])
        
        r_Ld_left = -r_realigned[n_Lo-1:][::-1]
        z_Ld_left = z_realigned[n_Lo-1:][::-1]
        r_Ld_right = r_realigned[n_Lo-1:]
        z_Ld_right = z_realigned[n_Lo-1:]
        
        # Plot Ld in Green
        ax.plot(r_Ld_left, z_Ld_left, color=color_Ld, linewidth=2.8, label=r'$L_d$ Matrix (Green)')
        ax.plot(r_Ld_right, z_Ld_right, color=color_Ld, linewidth=2.8)
        
        # Plot Lo in Red
        ax.plot(r_Lo_symm, z_Lo_symm, color=color_Lo, linewidth=3.5, label=r'$L_o$ Domain (Red)')
        
        # Domain boundary contact points (Black dots)
        ax.plot([-r_b_val, r_b_val], [z_rim, z_rim], 'ko', markersize=4.5, zorder=5)
        
        # Central axis of symmetry
        ax.axvline(0, color='gray', linestyle=':', linewidth=0.8, alpha=0.7)
        
        # Depth callout
        ax.annotate(
            f"Depth: {depth_val:.2f} $\\mu\\mathrm{{m}}$\n$2r_b = {2*r_b_val:.2f} \\mu\\mathrm{{m}}$",
            xy=(0, z_apex), xytext=(0, z_apex + (0.9 if depth_val > 2.0 else -1.5)),
            fontsize=8.5, ha='center', fontweight='bold', color='#7f0000',
            bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.85, edgecolor='#cccccc')
        )
        
        ax.set_title(f"\\textbf{{({chr(97+i)})}} $D_{{L_o}} = {D:.0f}\\,\\mu\\mathrm{{m}}$ ({depth_val:.2f} $\\mu\\mathrm{{m}}$ depth)", pad=4)
        ax.set_aspect('equal')
        ax.set_xlim(-6.2, 6.2)
        ax.set_ylim(-6.2, 6.2)
        ax.grid(True, linestyle=':', alpha=0.5)
        
        if i == 0:
            ax.legend(frameon=True, facecolor='white', framealpha=0.9, loc='upper right', fontsize=7.5)
            
    # Subplot 14: Combined Overlay of all vertical slices
    ax_comp = axes_flat[13]
    for p in profiles_data:
        r_v = p['r']
        z_v = p['z']
        n_L = p['n_Lo']
        r_Lo_s = np.concatenate([-r_v[:n_L][::-1], r_v[:n_L]])
        z_Lo_s = np.concatenate([z_v[:n_L][::-1], z_v[:n_L]])
        ax_comp.plot(r_Lo_s, z_Lo_s, color=color_Lo, alpha=0.6, linewidth=1.5)
        
    # Plot outer Ld shell of 14um for reference
    r_last = profiles_data[-1]['r']
    z_last = profiles_data[-1]['z']
    n_last = profiles_data[-1]['n_Lo']
    ax_comp.plot(-r_last[n_last-1:], z_last[n_last-1:], color=color_Ld, linewidth=2.0)
    ax_comp.plot(r_last[n_last-1:], z_last[n_last-1:], color=color_Ld, linewidth=2.0)
    ax_comp.set_title(r"\textbf{(n)} Composite Overlay ($2 \to 14\,\mu\mathrm{m}$)", pad=4)
    ax_comp.set_aspect('equal')
    ax_comp.set_xlim(-6.2, 6.2)
    ax_comp.set_ylim(-6.2, 6.2)
    ax_comp.grid(True, linestyle=':', alpha=0.5)
    
    # Hide unused subplots (indices 14 and 15)
    axes_flat[14].axis('off')
    axes_flat[15].axis('off')
    
    # Add common labels
    for r_idx in range(4):
        axes[r_idx, 0].set_ylabel(r'Axial Coordinate $z$ ($\mu\mathrm{m}$)')
    for c_idx in range(4):
        axes[3, c_idx].set_xlabel(r'Radial Coordinate $r$ ($\mu\mathrm{m}$)')
    # If bottom cells are empty, put label on cell above
    axes[2, 2].set_xlabel(r'Radial Coordinate $r$ ($\mu\mathrm{m}$)')
    axes[2, 3].set_xlabel(r'Radial Coordinate $r$ ($\mu\mathrm{m}$)')
    
    output_path = os.path.join(out_dir, "vertical_slices_2um_to_14um.png")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Vertical slices figure saved successfully to: {output_path}")

if __name__ == "__main__":
    generate_vertical_slices()
