"""
Plotting and Visualization Module for Two-Phase Vesicle Remodeling.
Generates comprehensive multi-panel publication-grade figures.
"""

import os
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

def generate_figures():
    results_dir = "simulation_results"
    trace_path = os.path.join(results_dir, "trace_remodeling_sweep.csv")
    
    # Read trace data
    data = {}
    with open(trace_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            branch = row["branch"]
            if branch not in data:
                data[branch] = {
                    "theta": [],
                    "kappa_Ld": [],
                    "E_total": [],
                    "E_bend_Lo": [],
                    "E_bend_Ld": [],
                    "E_line": [],
                    "E_ade": [],
                    "r_b": [],
                    "dimple_depth": []
                }
            data[branch]["theta"].append(float(row["theta"]))
            data[branch]["kappa_Ld"].append(float(row["kappa_Ld"]))
            data[branch]["E_total"].append(float(row["E_total"]))
            data[branch]["E_bend_Lo"].append(float(row["E_bend_Lo"]))
            data[branch]["E_bend_Ld"].append(float(row["E_bend_Ld"]))
            data[branch]["E_line"].append(float(row["E_line"]))
            data[branch]["E_ade"].append(float(row["E_ade"]))
            data[branch]["r_b"].append(float(row["r_b"]))
            data[branch]["dimple_depth"].append(float(row["dimple_depth"]))
            
    # Set high-quality styling
    plt.rcParams.update({
        'font.size': 11,
        'axes.labelsize': 12,
        'axes.titlesize': 13,
        'xtick.labelsize': 11,
        'ytick.labelsize': 11,
        'legend.fontsize': 10,
        'lines.linewidth': 2.0,
        'figure.dpi': 300
    })
    
    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(2, 2, figure=fig, wspace=0.28, hspace=0.32)
    
    # Panel A: Free Energy Landscape vs Streptavidin Coverage
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.plot(data["sphere"]["theta"], data["sphere"]["E_total"], 'o--', color='#2b5c8f', label=r'Quasi-Spherical $L_o$')
    ax1.plot(data["outward"]["theta"], data["outward"]["E_total"], 's--', color='#d95f02', label=r'Outward Budded $L_o$')
    ax1.plot(data["inward"]["theta"], data["inward"]["E_total"], '^-', color='#1b9e77', linewidth=2.8, label=r'Inward Dimpled $L_o$ (Ground State)')
    
    ax1.set_xlabel(r'Streptavidin Coverage $\theta$ on Inner Leaflet of $L_d$')
    ax1.set_ylabel(r'Total Free Energy $\mathcal{F}_{\mathrm{tot}}$ ($k_B T$)')
    ax1.set_title(r'\textbf{(a)} Thermodynamic Bifurcation: Ground State Selection')
    ax1.grid(True, linestyle=':', alpha=0.6)
    ax1.legend(frameon=True, facecolor='white', framealpha=0.9)
    
    # Panel B: Ld Bending Modulus & Boundary Neck Radius
    ax2 = fig.add_subplot(gs[0, 1])
    color_k = '#7570b3'
    color_rb = '#e7298a'
    
    ax2.plot(data["inward"]["theta"], data["inward"]["kappa_Ld"], 'd-', color=color_k, label=r'$\kappa_{L_d}(\theta)$')
    ax2.set_xlabel(r'Streptavidin Coverage $\theta$')
    ax2.set_ylabel(r'$L_d$ Bending Rigidity $\kappa_{L_d}$ ($k_B T$)', color=color_k)
    ax2.tick_params(axis='y', labelcolor=color_k)
    ax2.set_ylim(15, 100)
    
    ax2_twin = ax2.twinx()
    ax2_twin.plot(data["inward"]["theta"], data["inward"]["r_b"], 'v--', color=color_rb, label=r'Neck Radius $r_b$')
    ax2_twin.set_ylabel(r'Boundary Neck Radius $r_b$ ($\mu\mathrm{m}$)', color=color_rb)
    ax2_twin.tick_params(axis='y', labelcolor=color_rb)
    ax2.set_title(r'\textbf{(b)} $L_d$ Stiffening \& Line-Tension Neck Constriction')
    ax2.grid(True, linestyle=':', alpha=0.6)
    
    # Panel C: Cross-Sectional Axisymmetric Profiles
    ax3 = fig.add_subplot(gs[1, 0])
    
    # Load theta = 0.0 sphere and theta = 1.0 inward profiles
    p0_path = os.path.join(results_dir, "profile_theta_0.0_sphere.csv")
    p1_in_path = os.path.join(results_dir, "profile_theta_1.0_inward.csv")
    p1_out_path = os.path.join(results_dir, "profile_theta_1.0_outward.csv")
    
    p0 = np.genfromtxt(p0_path, delimiter=",", skip_header=1)
    p1_in = np.genfromtxt(p1_in_path, delimiter=",", skip_header=1)
    p1_out = np.genfromtxt(p1_out_path, delimiter=",", skip_header=1)
    
    # Split index (first 40 points are Lo)
    n_Lo = 40
    
    # Bare vesicle profile
    r0, z0 = p0[:, 2], p0[:, 3]
    # Saturated inward profile
    r1_in, z1_in = p1_in[:, 2], p1_in[:, 3]
    # Saturated outward profile
    r1_out, z1_out = p1_out[:, 2], p1_out[:, 3]
    
    # Plot symmetric profiles (r and -r)
    # Theta = 0.0 (Bare)
    ax3.plot(np.concatenate([-r0[:n_Lo][::-1], r0[:n_Lo]]),
             np.concatenate([z0[:n_Lo][::-1], z0[:n_Lo]]),
             color='#e41a1c', linewidth=2.5, linestyle=':', label=r'Bare $L_o$ ($\theta=0$)')
    ax3.plot(np.concatenate([-r0[n_Lo:][::-1], r0[n_Lo:]]),
             np.concatenate([z0[n_Lo:][::-1], z0[n_Lo:]]),
             color='#377eb8', linewidth=1.5, linestyle=':', label=r'Bare $L_d$ Matrix ($\theta=0$)')
             
    # Theta = 1.0 (Inward Dimple)
    ax3.plot(np.concatenate([-r1_in[:n_Lo][::-1], r1_in[:n_Lo]]),
             np.concatenate([z1_in[:n_Lo][::-1], z1_in[:n_Lo]]),
             color='#e41a1c', linewidth=3.5, label=r'Remodeled $L_o$ Dimple ($\theta=1$)')
    ax3.plot(np.concatenate([-r1_in[n_Lo:][::-1], r1_in[n_Lo:]]),
             np.concatenate([z1_in[n_Lo:][::-1], z1_in[n_Lo:]]),
             color='#377eb8', linewidth=2.5, label=r'Stiffened $L_d$ Sphere ($\theta=1$)')
             
    ax3.set_xlabel(r'Radial Coordinate $r$ ($\mu\mathrm{m}$)')
    ax3.set_ylabel(r'Axial Coordinate $z$ ($\mu\mathrm{m}$)')
    ax3.set_title(r'\textbf{(c)} Remodeled Vesicle Cross-Section: $L_o$ Inward Dimple')
    ax3.axis('equal')
    ax3.grid(True, linestyle=':', alpha=0.6)
    ax3.legend(frameon=True, facecolor='white', framealpha=0.9, loc='lower right')
    
    # Panel D: 3D Surface Reconstruction of the Dimpled Vesicle
    ax4 = fig.add_subplot(gs[1, 1], projection='3d')
    
    # Create 3D mesh by revolving r(s), z(s) around z-axis
    phi = np.linspace(0, 2*np.pi, 60)
    
    # Lo domain surface
    r_Lo = r1_in[:n_Lo]
    z_Lo = z1_in[:n_Lo]
    R_Lo_grid, Phi_Lo_grid = np.meshgrid(r_Lo, phi)
    Z_Lo_grid = np.tile(z_Lo, (len(phi), 1))
    X_Lo = R_Lo_grid * np.cos(Phi_Lo_grid)
    Y_Lo = R_Lo_grid * np.sin(Phi_Lo_grid)
    
    # Ld domain surface
    r_Ld = r1_in[n_Lo:]
    z_Ld = z1_in[n_Lo:]
    R_Ld_grid, Phi_Ld_grid = np.meshgrid(r_Ld, phi)
    Z_Ld_grid = np.tile(z_Ld, (len(phi), 1))
    X_Ld = R_Ld_grid * np.cos(Phi_Ld_grid)
    Y_Ld = R_Ld_grid * np.sin(Phi_Ld_grid)
    
    # Plot surfaces
    ax4.plot_surface(X_Lo, Y_Lo, Z_Lo_grid, color='#e41a1c', alpha=0.95, edgecolor='none')
    ax4.plot_surface(X_Ld, Y_Ld, Z_Ld_grid, color='#377eb8', alpha=0.45, edgecolor='none')
    
    ax4.set_xlabel(r'$x$ ($\mu\mathrm{m}$)', labelpad=2)
    ax4.set_ylabel(r'$y$ ($\mu\mathrm{m}$)', labelpad=2)
    ax4.set_zlabel(r'$z$ ($\mu\mathrm{m}$)', labelpad=2)
    ax4.set_title(r'\textbf{(d)} 3D Morphology: $L_o$ (Red) Dimple in $L_d$ (Blue)')
    ax4.view_init(elev=20, azim=45)
    
    output_fig = os.path.join(results_dir, "bifurcation_and_dimple_profiles.png")
    plt.tight_layout()
    plt.savefig(output_fig, dpi=300)
    plt.close()
    print(f"Publication-grade multi-panel figure saved to: {output_fig}")

if __name__ == "__main__":
    generate_figures()
