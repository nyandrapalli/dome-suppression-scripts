"""
High-Resolution 3D Visualization of Maximum Dimpled Vesicles:
Comparing Single-Domain vs Two-Domain Remodeling (D_Lo = 3.0 um, R_ves = 5.0 um).
Non-transparent (opaque alpha=1.0) 3D looking mesh with realistic directional lighting,
coloring Lo in Red and Ld in Green, completely clean with only panel labels.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LightSource

def render_3d_max_dimpled_comparison(
    out_path="simulation_results/3d_max_dimpled_vesicle_comparison.png"
):
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    
    # -------------------------------------------------------------------------
    # 1. Load Single Domain Numerical Profiles
    # -------------------------------------------------------------------------
    f_s_0 = "simulation_results_prl/param_sweeps/profiles/prof_DLo_3.0um_baseline_theta_0.00.csv"
    f_s_dimple = "simulation_results_prl/param_sweeps/profiles/prof_DLo_3.0um_gamma_1.5pN_theta_0.20.csv"
    
    p1_0 = pd.read_csv(f_s_0)
    p1_d = pd.read_csv(f_s_dimple)
    
    z_mid_s0 = (p1_0["z_um"].max() + p1_0["z_um"].min()) / 2.0
    z_mid_sd = (p1_d["z_um"].max() + p1_d["z_um"].min()) / 2.0
    
    z1_0 = p1_0["z_um"].values - z_mid_s0
    r1_0 = p1_0["r_um"].values
    
    z1_d = p1_d["z_um"].values - z_mid_sd
    r1_d = p1_d["r_um"].values
    
    idx_b1 = int(np.argmin(np.abs(p1_d["s_um"].values - 1.50)))
    r_rim_s = r1_d[idx_b1]
    z_rim_s = z1_d[idx_b1]

    # -------------------------------------------------------------------------
    # 2. Load Two Domain Numerical Profiles
    # -------------------------------------------------------------------------
    f_t_0 = "simulation_results_prl/two_domains_sweeps/profiles/prof_two_dom_DLo_3.0um_baseline_theta_0.00.csv"
    f_t_dimple = "simulation_results_prl/two_domains_sweeps/profiles/prof_two_dom_DLo_3.0um_gamma_1.5pN_theta_0.20.csv"
    
    p2_0 = pd.read_csv(f_t_0)
    p2_d = pd.read_csv(f_t_dimple)
    
    z_mid_t0 = (p2_0["z_um"].iloc[0] + p2_0["z_um"].iloc[-1]) / 2.0
    z_mid_td = (p2_d["z_um"].iloc[0] + p2_d["z_um"].iloc[-1]) / 2.0
    
    z2_0 = p2_0["z_um"].values - z_mid_t0
    r2_0 = p2_0["r_um"].values
    
    z2_d = p2_d["z_um"].values - z_mid_td
    r2_d = p2_d["r_um"].values
    
    idx_b_north = int(np.argmin(np.abs(p2_d["s_um"].values - 1.50)))
    s_tot = p2_d["s_um"].values[-1]
    idx_b_south = int(np.argmin(np.abs(p2_d["s_um"].values - (s_tot - 1.50))))
    
    r_rim_n = r2_d[idx_b_north]
    z_rim_n = z2_d[idx_b_north]
    r_rim_so = r2_d[idx_b_south]
    z_rim_so = z2_d[idx_b_south]

    # -------------------------------------------------------------------------
    # 3. Styling & Color Palette (Lo in Red, Ld in Green, Opaque alpha=1.0)
    # -------------------------------------------------------------------------
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 11,
        "axes.labelsize": 11.5,
        "axes.titlesize": 13.5,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "figure.dpi": 300
    })

    fig = plt.figure(figsize=(17, 11))
    gs = GridSpec(2, 3, figure=fig, hspace=0.20, wspace=0.18, width_ratios=[1.2, 1.2, 1.0])

    color_Lo = "#e41a1c"       # Vibrant red for Lo domains
    edge_Lo = "#7f0000"        # Dark red wireframe lines
    color_Ld = "#41ab5d"       # Vibrant emerald green for Ld membrane
    edge_Ld = "#1b4f24"        # Dark green wireframe lines
    
    # Directional lighting source for realistic 3D appearance
    ls = LightSource(azdeg=315, altdeg=50)

    phi_full = np.linspace(0, 2 * np.pi, 90)
    phi_cut = np.linspace(-np.pi / 2, np.pi, 70)  # 270 degree cutaway

    # =========================================================================
    # ROW 1: SINGLE DOMAIN (a1, a2, a3)
    # =========================================================================
    
    # Panel (a1): 3D Surface Mesh (Non-transparent)
    ax_a1 = fig.add_subplot(gs[0, 0], projection="3d")
    
    r_Lo_s = r1_d[:idx_b1+1]
    z_Lo_s = z1_d[:idx_b1+1]
    R_Lo_s, P_Lo_s = np.meshgrid(r_Lo_s, phi_full)
    Z_Lo_s = np.tile(z_Lo_s, (len(phi_full), 1))
    X_Lo_s = R_Lo_s * np.cos(P_Lo_s)
    Y_Lo_s = R_Lo_s * np.sin(P_Lo_s)
    
    r_Ld_s = r1_d[idx_b1:]
    z_Ld_s = z1_d[idx_b1:]
    R_Ld_s, P_Ld_s = np.meshgrid(r_Ld_s, phi_full)
    Z_Ld_s = np.tile(z_Ld_s, (len(phi_full), 1))
    X_Ld_s = R_Ld_s * np.cos(P_Ld_s)
    Y_Ld_s = R_Ld_s * np.sin(P_Ld_s)
    
    # Opaque 3D mesh with directional lighting and wireframe edges
    ax_a1.plot_surface(X_Lo_s, Y_Lo_s, Z_Lo_s, color=color_Lo, alpha=1.0, 
                       shade=True, lightsource=ls, rstride=1, cstride=1,
                       edgecolor=edge_Lo, linewidth=0.25, antialiased=True)
    ax_a1.plot_surface(X_Ld_s, Y_Ld_s, Z_Ld_s, color=color_Ld, alpha=1.0, 
                       shade=True, lightsource=ls, rstride=1, cstride=1,
                       edgecolor=edge_Ld, linewidth=0.20, antialiased=True)
    
    ax_a1.set_title(r"$\mathbf{(a1)}$", loc="left", fontweight="bold", pad=8)
    ax_a1.set_xlabel(r"$X$ ($\mu\mathrm{m}$)", labelpad=3)
    ax_a1.set_ylabel(r"$Y$ ($\mu\mathrm{m}$)", labelpad=3)
    ax_a1.set_zlabel(r"$Z$ ($\mu\mathrm{m}$)", labelpad=3)
    ax_a1.view_init(elev=36, azim=45)
    ax_a1.set_box_aspect([1, 1, 1])

    # Panel (a2): 3D Cutaway View (Non-transparent)
    ax_a2 = fig.add_subplot(gs[0, 1], projection="3d")
    
    R_Lo_cut, P_Lo_cut = np.meshgrid(r_Lo_s, phi_cut)
    Z_Lo_cut = np.tile(z_Lo_s, (len(phi_cut), 1))
    X_Lo_cut = R_Lo_cut * np.cos(P_Lo_cut)
    Y_Lo_cut = R_Lo_cut * np.sin(P_Lo_cut)
    
    R_Ld_cut, P_Ld_cut = np.meshgrid(r_Ld_s, phi_cut)
    Z_Ld_cut = np.tile(z_Ld_s, (len(phi_cut), 1))
    X_Ld_cut = R_Ld_cut * np.cos(P_Ld_cut)
    Y_Ld_cut = R_Ld_cut * np.sin(P_Ld_cut)
    
    ax_a2.plot_surface(X_Lo_cut, Y_Lo_cut, Z_Lo_cut, color=color_Lo, alpha=1.0, 
                       shade=True, lightsource=ls, rstride=1, cstride=1,
                       edgecolor=edge_Lo, linewidth=0.25, antialiased=True)
    ax_a2.plot_surface(X_Ld_cut, Y_Ld_cut, Z_Ld_cut, color=color_Ld, alpha=1.0, 
                       shade=True, lightsource=ls, rstride=1, cstride=1,
                       edgecolor=edge_Ld, linewidth=0.20, antialiased=True)
    
    # Boundary cut lines
    ax_a2.plot(r1_d, np.zeros_like(r1_d), z1_d, color="#111111", linewidth=2.0)
    ax_a2.plot(np.zeros_like(r1_d), -r1_d, z1_d, color="#111111", linewidth=2.0)
    
    ax_a2.set_title(r"$\mathbf{(a2)}$", loc="left", fontweight="bold", pad=8)
    ax_a2.set_xlabel(r"$X$ ($\mu\mathrm{m}$)", labelpad=3)
    ax_a2.set_ylabel(r"$Y$ ($\mu\mathrm{m}$)", labelpad=3)
    ax_a2.set_zlabel(r"$Z$ ($\mu\mathrm{m}$)", labelpad=3)
    ax_a2.view_init(elev=20, azim=40)
    ax_a2.set_box_aspect([1, 1, 1])

    # Panel (a3): Cross-Section (Lo in Red, Ld in Green, Clean)
    ax_a3 = fig.add_subplot(gs[0, 2])
    
    ax_a3.plot(r1_0, z1_0, "--", color="#7f8c8d", lw=1.8)
    ax_a3.plot(-r1_0, z1_0, "--", color="#7f8c8d", lw=1.8)
    
    ax_a3.plot(r_Ld_s, z_Ld_s, color=color_Ld, linewidth=2.8)
    ax_a3.plot(-r_Ld_s, z_Ld_s, color=color_Ld, linewidth=2.8)
    
    r_full_Lo_s = np.concatenate([-r_Lo_s[::-1], r_Lo_s])
    z_full_Lo_s = np.concatenate([z_Lo_s[::-1], z_Lo_s])
    ax_a3.plot(r_full_Lo_s, z_full_Lo_s, color=color_Lo, linewidth=3.5)
    
    ax_a3.plot([-r_rim_s, r_rim_s], [z_rim_s, z_rim_s], "ko", markersize=5.5)
    
    ax_a3.set_title(r"$\mathbf{(a3)}$", loc="left", fontweight="bold", pad=8)
    ax_a3.set_xlabel(r"$r$ ($\mu\mathrm{m}$)")
    ax_a3.set_ylabel(r"$z$ ($\mu\mathrm{m}$)")
    ax_a3.set_aspect("equal")
    ax_a3.set_xlim(-6.2, 6.2)
    ax_a3.set_ylim(-6.2, 6.2)
    ax_a3.grid(True, linestyle=":", alpha=0.5)

    # =========================================================================
    # ROW 2: TWO DOMAINS (b1, b2, b3)
    # =========================================================================
    
    # Panel (b1): 3D Surface View (Non-transparent)
    ax_b1 = fig.add_subplot(gs[1, 0], projection="3d")
    
    r_Lo1 = r2_d[:idx_b_north+1]
    z_Lo1 = z2_d[:idx_b_north+1]
    R_Lo1, P_Lo1 = np.meshgrid(r_Lo1, phi_full)
    Z_Lo1 = np.tile(z_Lo1, (len(phi_full), 1))
    X_Lo1 = R_Lo1 * np.cos(P_Lo1)
    Y_Lo1 = R_Lo1 * np.sin(P_Lo1)
    
    r_Ld_t = r2_d[idx_b_north:idx_b_south+1]
    z_Ld_t = z2_d[idx_b_north:idx_b_south+1]
    R_Ld_t, P_Ld_t = np.meshgrid(r_Ld_t, phi_full)
    Z_Ld_t = np.tile(z_Ld_t, (len(phi_full), 1))
    X_Ld_t = R_Ld_t * np.cos(P_Ld_t)
    Y_Ld_t = R_Ld_t * np.sin(P_Ld_t)
    
    r_Lo2 = r2_d[idx_b_south:]
    z_Lo2 = z2_d[idx_b_south:]
    R_Lo2, P_Lo2 = np.meshgrid(r_Lo2, phi_full)
    Z_Lo2 = np.tile(z_Lo2, (len(phi_full), 1))
    X_Lo2 = R_Lo2 * np.cos(P_Lo2)
    Y_Lo2 = R_Lo2 * np.sin(P_Lo2)
    
    ax_b1.plot_surface(X_Lo1, Y_Lo1, Z_Lo1, color=color_Lo, alpha=1.0, 
                       shade=True, lightsource=ls, rstride=1, cstride=1,
                       edgecolor=edge_Lo, linewidth=0.25, antialiased=True)
    ax_b1.plot_surface(X_Ld_t, Y_Ld_t, Z_Ld_t, color=color_Ld, alpha=1.0, 
                       shade=True, lightsource=ls, rstride=1, cstride=1,
                       edgecolor=edge_Ld, linewidth=0.20, antialiased=True)
    ax_b1.plot_surface(X_Lo2, Y_Lo2, Z_Lo2, color=color_Lo, alpha=1.0, 
                       shade=True, lightsource=ls, rstride=1, cstride=1,
                       edgecolor=edge_Lo, linewidth=0.25, antialiased=True)
    
    ax_b1.set_title(r"$\mathbf{(b1)}$", loc="left", fontweight="bold", pad=8)
    ax_b1.set_xlabel(r"$X$ ($\mu\mathrm{m}$)", labelpad=3)
    ax_b1.set_ylabel(r"$Y$ ($\mu\mathrm{m}$)", labelpad=3)
    ax_b1.set_zlabel(r"$Z$ ($\mu\mathrm{m}$)", labelpad=3)
    ax_b1.view_init(elev=36, azim=45)
    ax_b1.set_box_aspect([1, 1, 1])

    # Panel (b2): 3D Cutaway View (Non-transparent)
    ax_b2 = fig.add_subplot(gs[1, 1], projection="3d")
    
    R1_cut, P1_cut = np.meshgrid(r_Lo1, phi_cut)
    Z1_cut = np.tile(z_Lo1, (len(phi_cut), 1))
    X1_cut = R1_cut * np.cos(P1_cut)
    Y1_cut = R1_cut * np.sin(P1_cut)
    
    R_Ld_t_cut, P_Ld_t_cut = np.meshgrid(r_Ld_t, phi_cut)
    Z_Ld_t_cut = np.tile(z_Ld_t, (len(phi_cut), 1))
    X_Ld_t_cut = R_Ld_t_cut * np.cos(P_Ld_t_cut)
    Y_Ld_t_cut = R_Ld_t_cut * np.sin(P_Ld_t_cut)
    
    R2_cut, P2_cut = np.meshgrid(r_Lo2, phi_cut)
    Z2_cut = np.tile(z_Lo2, (len(phi_cut), 1))
    X2_cut = R2_cut * np.cos(P2_cut)
    Y2_cut = R2_cut * np.sin(P2_cut)
    
    ax_b2.plot_surface(X1_cut, Y1_cut, Z1_cut, color=color_Lo, alpha=1.0, 
                       shade=True, lightsource=ls, rstride=1, cstride=1,
                       edgecolor=edge_Lo, linewidth=0.25, antialiased=True)
    ax_b2.plot_surface(X_Ld_t_cut, Y_Ld_t_cut, Z_Ld_t_cut, color=color_Ld, alpha=1.0, 
                       shade=True, lightsource=ls, rstride=1, cstride=1,
                       edgecolor=edge_Ld, linewidth=0.20, antialiased=True)
    ax_b2.plot_surface(X2_cut, Y2_cut, Z2_cut, color=color_Lo, alpha=1.0, 
                       shade=True, lightsource=ls, rstride=1, cstride=1,
                       edgecolor=edge_Lo, linewidth=0.25, antialiased=True)
    
    # Boundary cut lines
    ax_b2.plot(r2_d, np.zeros_like(r2_d), z2_d, color="#111111", linewidth=2.0)
    ax_b2.plot(np.zeros_like(r2_d), -r2_d, z2_d, color="#111111", linewidth=2.0)
    
    ax_b2.set_title(r"$\mathbf{(b2)}$", loc="left", fontweight="bold", pad=8)
    ax_b2.set_xlabel(r"$X$ ($\mu\mathrm{m}$)", labelpad=3)
    ax_b2.set_ylabel(r"$Y$ ($\mu\mathrm{m}$)", labelpad=3)
    ax_b2.set_zlabel(r"$Z$ ($\mu\mathrm{m}$)", labelpad=3)
    ax_b2.view_init(elev=18, azim=40)
    ax_b2.set_box_aspect([1, 1, 1])

    # Panel (b3): Cross-Section (Lo in Red, Ld in Green, Clean)
    ax_b3 = fig.add_subplot(gs[1, 2])
    
    ax_b3.plot(r2_0, z2_0, "--", color="#7f8c8d", lw=1.8)
    ax_b3.plot(-r2_0, z2_0, "--", color="#7f8c8d", lw=1.8)
    
    ax_b3.plot(r_Ld_t, z_Ld_t, color=color_Ld, linewidth=2.8)
    ax_b3.plot(-r_Ld_t, z_Ld_t, color=color_Ld, linewidth=2.8)
    
    r_full_Lo1 = np.concatenate([-r_Lo1[::-1], r_Lo1])
    z_full_Lo1 = np.concatenate([z_Lo1[::-1], z_Lo1])
    r_full_Lo2 = np.concatenate([-r_Lo2[::-1], r_Lo2])
    z_full_Lo2 = np.concatenate([z_Lo2[::-1], z_Lo2])
    
    ax_b3.plot(r_full_Lo1, z_full_Lo1, color=color_Lo, linewidth=3.5)
    ax_b3.plot(r_full_Lo2, z_full_Lo2, color=color_Lo, linewidth=3.5)
    
    ax_b3.plot([-r_rim_n, r_rim_n], [z_rim_n, z_rim_n], "ko", markersize=5.5)
    ax_b3.plot([-r_rim_so, r_rim_so], [z_rim_so, z_rim_so], "ko", markersize=5.5)
    
    ax_b3.set_title(r"$\mathbf{(b3)}$", loc="left", fontweight="bold", pad=8)
    ax_b3.set_xlabel(r"$r$ ($\mu\mathrm{m}$)")
    ax_b3.set_ylabel(r"$z$ ($\mu\mathrm{m}$)")
    ax_b3.set_aspect("equal")
    ax_b3.set_xlim(-6.2, 6.2)
    ax_b3.set_ylim(-6.2, 6.2)
    ax_b3.grid(True, linestyle=":", alpha=0.5)

    plt.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    print(f"3D comparison successfully saved: {out_path}")

if __name__ == "__main__":
    render_3d_max_dimpled_comparison()
