"""
Comparative Visualization: 2D Axisymmetric Continuum Simulation vs 3D Triangulated Vertex Simulation.
Demonstrates:
1. Left Panel: 2D Continuum Profile from 1D Euler-Lagrange BVP solver (r(s), z(s))
2. Middle Panel: Direct Cross-Sectional Overlay (1D ODE curve vs 3D mesh planar slice through y=0)
3. Right Panel: Full 3D Triangulated Mesh in Perspective View
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def slice_mesh_plane_y(vertices, faces):
    """
    Computes the exact cross-section intersection of the 3D triangular mesh
    with the vertical cutting plane y = 0.
    Returns line segments (x, z) lying on the cutting plane.
    """
    segments = []
    for tri in faces:
        v = vertices[tri] # (3, 3)
        y = v[:, 1]
        
        # Check if triangle crosses y = 0
        if np.all(y > 0) or np.all(y < 0):
            continue
            
        # Find edge intersections with y = 0
        pts = []
        edges = [(0, 1), (1, 2), (2, 0)]
        for i1, i2 in edges:
            y1, y2 = y[i1], y[i2]
            if (y1 <= 0 and y2 >= 0) or (y1 >= 0 and y2 <= 0):
                if abs(y2 - y1) > 1e-12:
                    t = -y1 / (y2 - y1)
                    pt = v[i1] + t * (v[i2] - v[i1])
                    pts.append(pt)
                    
        if len(pts) >= 2:
            segments.append((pts[0][[0, 2]], pts[1][[0, 2]]))
            
    return segments

def plot_2d_vs_3d_comparison():
    root_dir = "HermesLoouter"
    out_fig = os.path.join(root_dir, "fig_comparative_2d_vs_3d.png")
    
    # 1. Load 2D profile
    prof_2d_path = os.path.join(root_dir, "double", "profiles", "prof_two_dom_DLo_3.0_theta_1.00.csv")
    if not os.path.exists(prof_2d_path):
        print(f"File not found: {prof_2d_path}")
        return
        
    df_2d = pd.read_csv(prof_2d_path)
    r_2d = df_2d["r_um"].values
    z_2d = df_2d["z_um"].values
    dom_2d = df_2d["domain"].values
    
    # Center z at equator
    z_mid_2d = 0.5 * (np.max(z_2d) + np.min(z_2d))
    z_aligned_2d = z_2d - z_mid_2d
    
    # 2. Load 3D mesh
    mesh_3d_obj = os.path.join(root_dir, "triangulated_3d", "mesh_equilibrated_sub3_v0.98_DLo3.0.obj")
    if not os.path.exists(mesh_3d_obj):
        print(f"File not found: {mesh_3d_obj}")
        return
        
    verts_3d = []
    faces_3d = []
    with open(mesh_3d_obj, "r") as f:
        for line in f:
            if line.startswith("v "):
                verts_3d.append([float(x) for x in line.strip().split()[1:4]])
            elif line.startswith("f "):
                # 1-indexed to 0-indexed
                faces_3d.append([int(x.split("/")[0]) - 1 for x in line.strip().split()[1:4]])
    verts_3d = np.array(verts_3d)
    faces_3d = np.array(faces_3d)
    
    # Domain assignments for 3D faces
    sys.path.insert(0, os.path.join(root_dir, "triangulated_3d"))
    from mesh_icosphere import assign_domains
    domains_3d = assign_domains(verts_3d, R_ves=5.0, D_Lo=3.0)
    
    # 3. Setup Figure
    fig = plt.figure(figsize=(20, 6.5), dpi=300)
    gs = GridSpec(1, 3, figure=fig, width_ratios=[1.0, 1.0, 1.15], wspace=0.25)
    
    color_Lo = "#d73027" # Red
    color_Ld = "#1a9850" # Green
    
    # -------------------------------------------------------------
    # PANEL A: 2D Axisymmetric Continuum Simulation
    # -------------------------------------------------------------
    ax_a = fig.add_subplot(gs[0, 0])
    mask_lo_n = (dom_2d == "Lo_north")
    mask_lo_s = (dom_2d == "Lo_south")
    mask_ld   = (dom_2d == "Ld")
    
    # North cap
    ax_a.plot(np.concatenate([-r_2d[mask_lo_n][::-1], r_2d[mask_lo_n]]),
              np.concatenate([z_aligned_2d[mask_lo_n][::-1], z_aligned_2d[mask_lo_n]]),
              color=color_Lo, linewidth=3.2, label=r"North $L_o$ Cap ($c_0 > 0$)")
    # Ld belt
    ax_a.plot(r_2d[mask_ld], z_aligned_2d[mask_ld], color=color_Ld, linewidth=2.8, label=r"Equatorial $L_d$ Belt (Bare)")
    ax_a.plot(-r_2d[mask_ld][::-1], z_aligned_2d[mask_ld][::-1], color=color_Ld, linewidth=2.8)
    # South cap
    ax_a.plot(np.concatenate([-r_2d[mask_lo_s][::-1], r_2d[mask_lo_s]]),
              np.concatenate([z_aligned_2d[mask_lo_s][::-1], z_aligned_2d[mask_lo_s]]),
              color=color_Lo, linewidth=3.2, label=r"South $L_o$ Cap ($c_0 > 0$)")
    
    ax_a.axvline(0, color="gray", linestyle=":", linewidth=0.8, alpha=0.7)
    ax_a.set_aspect("equal")
    ax_a.set_xlim(-6.5, 6.5)
    ax_a.set_ylim(-6.5, 6.5)
    ax_a.set_xlabel(r"Radial coordinate $r$ ($\mu\mathrm{m}$)", fontsize=11)
    ax_a.set_ylabel(r"Axial coordinate $z$ ($\mu\mathrm{m}$)", fontsize=11)
    ax_a.set_title(r"(a) 2D Mathematical Model" + "\n(Axisymmetric Continuum ODE)", fontsize=12, fontweight="bold")
    ax_a.grid(True, linestyle=":", alpha=0.5)
    ax_a.legend(loc="lower left", fontsize=8.5, framealpha=0.95)
    
    # -------------------------------------------------------------
    # PANEL B: Direct Overlay (2D Continuum vs 3D Mesh Slice)
    # -------------------------------------------------------------
    ax_b = fig.add_subplot(gs[0, 1])
    
    # Plot 2D curve as faint reference
    ax_b.plot(np.concatenate([-r_2d[::-1], r_2d]), 
              np.concatenate([z_aligned_2d[::-1], z_aligned_2d]), 
              color="black", linestyle="--", linewidth=1.5, alpha=0.7, label="2D ODE Continuum")
    
    # Compute 3D mesh cross-section at y = 0
    slice_segs = slice_mesh_plane_y(verts_3d, faces_3d)
    for p1, p2 in slice_segs:
        # Determine dominant z for color
        z_avg = 0.5 * (p1[1] + p2[1])
        c = color_Lo if abs(z_avg) > 4.2 else color_Ld
        ax_b.plot([p1[0], p2[0]], [p1[1], p2[1]], color=c, linewidth=2.8)
        
    # Dummy lines for legend
    ax_b.plot([], [], color=color_Lo, linewidth=2.8, label=r"3D Slice: $L_o$ Caps")
    ax_b.plot([], [], color=color_Ld, linewidth=2.8, label=r"3D Slice: $L_d$ Belt")
    
    ax_b.axvline(0, color="gray", linestyle=":", linewidth=0.8, alpha=0.7)
    ax_b.set_aspect("equal")
    ax_b.set_xlim(-6.5, 6.5)
    ax_b.set_ylim(-6.5, 6.5)
    ax_b.set_xlabel(r"Cartesian $x$ ($\mu\mathrm{m}$)", fontsize=11)
    ax_b.set_ylabel(r"Cartesian $z$ ($\mu\mathrm{m}$)", fontsize=11)
    ax_b.set_title(r"(b) Direct Overlay Comparison" + "\n(2D Continuum Line vs 3D Mesh Slice)", fontsize=12, fontweight="bold")
    ax_b.grid(True, linestyle=":", alpha=0.5)
    ax_b.legend(loc="lower left", fontsize=8.5, framealpha=0.95)
    
    # -------------------------------------------------------------
    # PANEL C: 3D Triangulated Vertex Simulation
    # -------------------------------------------------------------
    ax_c = fig.add_subplot(gs[0, 2], projection="3d")
    face_domains = np.mean(domains_3d[faces_3d], axis=1)
    face_colors = []
    for fd in face_domains:
        if (fd >= 0.8 and fd <= 1.2) or fd >= 1.8:
            face_colors.append(color_Lo)
        elif fd < 0.2:
            face_colors.append(color_Ld)
        else:
            face_colors.append("#fee08b")
            
    poly3d = [[verts_3d[idx] for idx in face] for face in faces_3d]
    collection = Poly3DCollection(poly3d, facecolors=face_colors, edgecolors="k", linewidths=0.12, alpha=0.92)
    ax_c.add_collection3d(collection)
    
    lim = 6.0
    ax_c.set_xlim(-lim, lim)
    ax_c.set_ylim(-lim, lim)
    ax_c.set_zlim(-lim, lim)
    ax_c.set_xlabel("X (um)", fontsize=9)
    ax_c.set_ylabel("Y (um)", fontsize=9)
    ax_c.set_zlabel("Z (um)", fontsize=9)
    ax_c.set_title(r"(c) 3D Triangulated Simulation" + "\n(Discrete Helfrich on 2-Manifold Mesh)", fontsize=12, fontweight="bold")
    ax_c.view_init(elev=18, azim=-42)
    ax_c.grid(False)
    
    plt.tight_layout()
    fig.savefig(out_fig)
    plt.close(fig)
    print(f">> Saved comparative figure: {out_fig}")

if __name__ == "__main__":
    plot_2d_vs_3d_comparison()
