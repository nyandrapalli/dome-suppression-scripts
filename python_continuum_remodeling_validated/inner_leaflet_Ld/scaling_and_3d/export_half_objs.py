"""
Generate 3D Sliced Wavefront OBJ Models for the Two-Domain Vesicle:
1. Equatorial Half (Cut along Z = 0 equator into two identical top/bottom halves):
   - Red Lo dimple at the pole, Green Ld membrane down to the equator.
2. Vertical / Meridian Half (Cut along Y = 0 meridian plane into two identical halves):
   - Shows both North and South Red Lo dimples and the interior lumen cross-section.

Both models feature:
- Colored Lo (Red) and Ld (Green) with Wavefront .mtl materials and embedded vertex colors.
- Exact outward normals for smooth 3D shading.
"""

import os
import shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LightSource

def export_equatorial_half_obj(
    profile_csv,
    out_obj_path,
    out_mtl_path,
    n_phi=96
):
    """Exports top half (Z >= 0) sliced along the equator."""
    df = pd.read_csv(profile_csv)
    r_full = df["r_um"].values.copy()
    z_full = df["z_um"].values.copy()
    s_full = df["s_um"].values.copy()
    
    z_mid = (z_full[0] + z_full[-1]) / 2.0
    z_full -= z_mid
    
    # Equatorial midpoint index (z = 0)
    mid_idx = len(df) // 2
    r = r_full[:mid_idx + 1]
    z = z_full[:mid_idx + 1]
    s = s_full[:mid_idx + 1]
    
    r[0] = 0.0
    r = np.abs(r)
    n_pts = len(r)
    
    is_Lo = (s <= 1.50)
    
    dr = np.gradient(r, s)
    dz = np.gradient(z, s)
    psi = np.arctan2(dz, dr)
    
    # MTL file
    mtl_filename = os.path.basename(out_mtl_path)
    with open(out_mtl_path, "w", encoding="utf-8") as fmtl:
        fmtl.write("# Material definitions for Equatorial Half Vesicle\n")
        fmtl.write("newmtl Lo_Red\nKa 0.15 0.02 0.02\nKd 0.894 0.102 0.110\nKs 0.40 0.40 0.40\nNs 40.0\nd 1.0\nillum 2\n\n")
        fmtl.write("newmtl Ld_Green\nKa 0.02 0.15 0.05\nKd 0.255 0.671 0.365\nKs 0.35 0.35 0.35\nNs 35.0\nd 1.0\nillum 2\n")

    phi = np.linspace(0, 2 * np.pi, n_phi, endpoint=False)
    cos_phi = np.cos(phi)
    sin_phi = np.sin(phi)
    
    color_Lo_rgb = (0.894, 0.102, 0.110)
    color_Ld_rgb = (0.255, 0.671, 0.365)
    
    vertices = []
    normals = []
    ring_start_idx = []
    
    # Ring 0: Pole
    c0 = color_Lo_rgb if is_Lo[0] else color_Ld_rgb
    vertices.append((0.0, 0.0, z[0], *c0))
    normals.append((0.0, 0.0, 1.0))
    ring_start_idx.append(0)
    
    for i in range(1, n_pts):
        ring_start_idx.append(len(vertices))
        ri = r[i]
        zi = z[i]
        ci = color_Lo_rgb if is_Lo[i] else color_Ld_rgb
        sin_psi = np.sin(psi[i])
        cos_psi = np.cos(psi[i])
        
        for j in range(n_phi):
            x = ri * cos_phi[j]
            y = ri * sin_phi[j]
            vertices.append((x, y, zi, *ci))
            nx = -sin_psi * cos_phi[j]
            ny = -sin_psi * sin_phi[j]
            nz = cos_psi
            norm = np.sqrt(nx*nx + ny*ny + nz*nz) + 1e-12
            normals.append((nx/norm, ny/norm, nz/norm))
            
    faces_Lo = []
    faces_Ld = []
    
    # Cap fan
    r1_start = ring_start_idx[1]
    is_cap_Lo = is_Lo[0] or is_Lo[1]
    for j in range(n_phi):
        j_next = (j + 1) % n_phi
        v_pole = 1
        v1 = r1_start + j + 1
        v2 = r1_start + j_next + 1
        face = (v_pole, v1, v2)
        if is_cap_Lo:
            faces_Lo.append(face)
        else:
            faces_Ld.append(face)
            
    # Quad bands down to equator ring (n_pts - 1)
    for i in range(1, n_pts - 1):
        curr_start = ring_start_idx[i]
        next_start = ring_start_idx[i + 1]
        band_Lo = is_Lo[i] and is_Lo[i + 1]
        target_list = faces_Lo if band_Lo else faces_Ld
        
        for j in range(n_phi):
            j_next = (j + 1) % n_phi
            v_curr_j = curr_start + j + 1
            v_curr_next = curr_start + j_next + 1
            v_next_j = next_start + j + 1
            v_next_next = next_start + j_next + 1
            
            target_list.append((v_curr_j, v_next_j, v_next_next))
            target_list.append((v_curr_j, v_next_next, v_curr_next))
            
    with open(out_obj_path, "w", encoding="utf-8") as fobj:
        fobj.write(f"# Wavefront OBJ: Two-Domain Vesicle (Equatorial Top Half, Z >= 0)\n")
        fobj.write(f"# Vertices: {len(vertices)}, Faces: {len(faces_Lo) + len(faces_Ld)}\n")
        fobj.write(f"mtllib {mtl_filename}\n\n")
        for v in vertices:
            fobj.write(f"v {v[0]:.6f} {v[1]:.6f} {v[2]:.6f} {v[3]:.4f} {v[4]:.4f} {v[5]:.4f}\n")
        fobj.write("\n")
        for n in normals:
            fobj.write(f"vn {n[0]:.6f} {n[1]:.6f} {n[2]:.6f}\n")
        fobj.write("\n")
        if len(faces_Lo) > 0:
            fobj.write("g Lo_Domain\nusemtl Lo_Red\n")
            for f in faces_Lo:
                fobj.write(f"f {f[0]}//{f[0]} {f[1]}//{f[1]} {f[2]}//{f[2]}\n")
            fobj.write("\n")
        if len(faces_Ld) > 0:
            fobj.write("g Ld_Membrane\nusemtl Ld_Green\n")
            for f in faces_Ld:
                fobj.write(f"f {f[0]}//{f[0]} {f[1]}//{f[1]} {f[2]}//{f[2]}\n")
            fobj.write("\n")
            
    print(f"Exported {out_obj_path} ({len(vertices)} verts, {len(faces_Lo)+len(faces_Ld)} faces)")


def export_vertical_half_obj(
    profile_csv,
    out_obj_path,
    out_mtl_path,
    n_phi=64
):
    """Exports vertical half (Y >= 0, phi in [0, pi]) showing cross-section of both dimples."""
    df = pd.read_csv(profile_csv)
    r = df["r_um"].values.copy()
    z = df["z_um"].values.copy()
    s = df["s_um"].values.copy()
    
    z_mid = (z[0] + z[-1]) / 2.0
    z -= z_mid
    r[0] = 0.0
    r[-1] = 0.0
    r = np.abs(r)
    n_pts = len(r)
    
    s_tot = s[-1]
    is_Lo = (s <= 1.50) | (s >= s_tot - 1.50)
    
    dr = np.gradient(r, s)
    dz = np.gradient(z, s)
    psi = np.arctan2(dz, dr)
    
    mtl_filename = os.path.basename(out_mtl_path)
    with open(out_mtl_path, "w", encoding="utf-8") as fmtl:
        fmtl.write("# Material definitions for Vertical Half Vesicle\n")
        fmtl.write("newmtl Lo_Red\nKa 0.15 0.02 0.02\nKd 0.894 0.102 0.110\nKs 0.40 0.40 0.40\nNs 40.0\nd 1.0\nillum 2\n\n")
        fmtl.write("newmtl Ld_Green\nKa 0.02 0.15 0.05\nKd 0.255 0.671 0.365\nKs 0.35 0.35 0.35\nNs 35.0\nd 1.0\nillum 2\n")

    # phi from 0 to pi (Y >= 0)
    phi = np.linspace(0, np.pi, n_phi)
    cos_phi = np.cos(phi)
    sin_phi = np.sin(phi)
    
    color_Lo_rgb = (0.894, 0.102, 0.110)
    color_Ld_rgb = (0.255, 0.671, 0.365)
    
    vertices = []
    normals = []
    ring_start_idx = []
    
    # Ring 0: North pole (single vertex)
    c0 = color_Lo_rgb if is_Lo[0] else color_Ld_rgb
    vertices.append((0.0, 0.0, z[0], *c0))
    normals.append((0.0, 0.0, 1.0))
    ring_start_idx.append(0)
    
    for i in range(1, n_pts - 1):
        ring_start_idx.append(len(vertices))
        ri = r[i]
        zi = z[i]
        ci = color_Lo_rgb if is_Lo[i] else color_Ld_rgb
        sin_psi = np.sin(psi[i])
        cos_psi = np.cos(psi[i])
        
        for j in range(n_phi):
            x = ri * cos_phi[j]
            y = ri * sin_phi[j]
            vertices.append((x, y, zi, *ci))
            nx = -sin_psi * cos_phi[j]
            ny = -sin_psi * sin_phi[j]
            nz = cos_psi
            norm = np.sqrt(nx*nx + ny*ny + nz*nz) + 1e-12
            normals.append((nx/norm, ny/norm, nz/norm))
            
    # South pole (single vertex)
    ring_start_idx.append(len(vertices))
    cN = color_Lo_rgb if is_Lo[-1] else color_Ld_rgb
    vertices.append((0.0, 0.0, z[-1], *cN))
    normals.append((0.0, 0.0, -1.0))
    
    faces_Lo = []
    faces_Ld = []
    
    # North cap fan: j from 0 to n_phi - 2 (non-wrapping)
    r1_start = ring_start_idx[1]
    is_cap_Lo = is_Lo[0] or is_Lo[1]
    for j in range(n_phi - 1):
        v_pole = 1
        v1 = r1_start + j + 1
        v2 = r1_start + j + 2
        face = (v_pole, v1, v2)
        if is_cap_Lo:
            faces_Lo.append(face)
        else:
            faces_Ld.append(face)
            
    # Intermediate quad bands
    for i in range(1, n_pts - 2):
        curr_start = ring_start_idx[i]
        next_start = ring_start_idx[i + 1]
        band_Lo = is_Lo[i] and is_Lo[i + 1]
        target_list = faces_Lo if band_Lo else faces_Ld
        
        for j in range(n_phi - 1):
            v_curr_j = curr_start + j + 1
            v_curr_next = curr_start + j + 2
            v_next_j = next_start + j + 1
            v_next_next = next_start + j + 2
            
            target_list.append((v_curr_j, v_next_j, v_next_next))
            target_list.append((v_curr_j, v_next_next, v_curr_next))
            
    # South cap fan
    last_ring_start = ring_start_idx[n_pts - 2]
    south_pole_idx = len(vertices)
    is_south_Lo = is_Lo[-1] or is_Lo[-2]
    for j in range(n_phi - 1):
        v1 = last_ring_start + j + 1
        v2 = last_ring_start + j + 2
        face = (v1, south_pole_idx, v2)
        if is_south_Lo:
            faces_Lo.append(face)
        else:
            faces_Ld.append(face)
            
    with open(out_obj_path, "w", encoding="utf-8") as fobj:
        fobj.write(f"# Wavefront OBJ: Two-Domain Vesicle (Vertical Half, Y >= 0)\n")
        fobj.write(f"# Vertices: {len(vertices)}, Faces: {len(faces_Lo) + len(faces_Ld)}\n")
        fobj.write(f"mtllib {mtl_filename}\n\n")
        for v in vertices:
            fobj.write(f"v {v[0]:.6f} {v[1]:.6f} {v[2]:.6f} {v[3]:.4f} {v[4]:.4f} {v[5]:.4f}\n")
        fobj.write("\n")
        for n in normals:
            fobj.write(f"vn {n[0]:.6f} {n[1]:.6f} {n[2]:.6f}\n")
        fobj.write("\n")
        if len(faces_Lo) > 0:
            fobj.write("g Lo_Domain\nusemtl Lo_Red\n")
            for f in faces_Lo:
                fobj.write(f"f {f[0]}//{f[0]} {f[1]}//{f[1]} {f[2]}//{f[2]}\n")
            fobj.write("\n")
        if len(faces_Ld) > 0:
            fobj.write("g Ld_Membrane\nusemtl Ld_Green\n")
            for f in faces_Ld:
                fobj.write(f"f {f[0]}//{f[0]} {f[1]}//{f[1]} {f[2]}//{f[2]}\n")
            fobj.write("\n")
            
    print(f"Exported {out_obj_path} ({len(vertices)} verts, {len(faces_Lo)+len(faces_Ld)} faces)")


def render_preview_figure():
    """Render a clean 3D preview of both halved pieces."""
    fig = plt.figure(figsize=(12, 6))
    ls = LightSource(azdeg=315, altdeg=50)
    
    # 1. Equatorial Half Preview
    ax1 = fig.add_subplot(121, projection="3d")
    df = pd.read_csv("simulation_results_prl/two_domains_sweeps/profiles/prof_two_dom_DLo_3.0um_gamma_1.5pN_theta_0.20.csv")
    r_full = df["r_um"].values
    z_full = df["z_um"].values - (df["z_um"].iloc[0] + df["z_um"].iloc[-1])/2.0
    s_full = df["s_um"].values
    
    mid_idx = len(df) // 2
    r_eq = r_full[:mid_idx + 1]
    z_eq = z_full[:mid_idx + 1]
    
    idx_b = int(np.argmin(np.abs(s_full[:mid_idx+1] - 1.50)))
    phi = np.linspace(0, 2*np.pi, 60)
    
    # Lo
    R_lo, P_lo = np.meshgrid(r_eq[:idx_b+1], phi)
    Z_lo = np.tile(z_eq[:idx_b+1], (len(phi), 1))
    ax1.plot_surface(R_lo*np.cos(P_lo), R_lo*np.sin(P_lo), Z_lo, color="#e41a1c",
                     shade=True, lightsource=ls, edgecolor="#7f0000", linewidth=0.2)
    # Ld
    R_ld, P_ld = np.meshgrid(r_eq[idx_b:], phi)
    Z_ld = np.tile(z_eq[idx_b:], (len(phi), 1))
    ax1.plot_surface(R_ld*np.cos(P_ld), R_ld*np.sin(P_ld), Z_ld, color="#41ab5d",
                     shade=True, lightsource=ls, edgecolor="#1b4f24", linewidth=0.2)
    
    ax1.set_title(r"$\mathbf{(a)}$", loc="left", fontweight="bold")
    ax1.set_xlabel(r"$X$ ($\mu\mathrm{m}$)")
    ax1.set_ylabel(r"$Y$ ($\mu\mathrm{m}$)")
    ax1.set_zlabel(r"$Z$ ($\mu\mathrm{m}$)")
    ax1.view_init(elev=25, azim=45)
    ax1.set_box_aspect([1, 1, 0.8])
    
    # 2. Vertical Half Preview
    ax2 = fig.add_subplot(122, projection="3d")
    idx_bn = int(np.argmin(np.abs(s_full - 1.50)))
    idx_bs = int(np.argmin(np.abs(s_full - (s_full[-1] - 1.50))))
    phi_v = np.linspace(0, np.pi, 45)
    
    # North Lo
    R_lo1, P_lo1 = np.meshgrid(r_full[:idx_bn+1], phi_v)
    Z_lo1 = np.tile(z_full[:idx_bn+1], (len(phi_v), 1))
    ax2.plot_surface(R_lo1*np.cos(P_lo1), R_lo1*np.sin(P_lo1), Z_lo1, color="#e41a1c",
                     shade=True, lightsource=ls, edgecolor="#7f0000", linewidth=0.2)
    # Middle Ld
    R_ldm, P_ldm = np.meshgrid(r_full[idx_bn:idx_bs+1], phi_v)
    Z_ldm = np.tile(z_full[idx_bn:idx_bs+1], (len(phi_v), 1))
    ax2.plot_surface(R_ldm*np.cos(P_ldm), R_ldm*np.sin(P_ldm), Z_ldm, color="#41ab5d",
                     shade=True, lightsource=ls, edgecolor="#1b4f24", linewidth=0.2)
    # South Lo
    R_lo2, P_lo2 = np.meshgrid(r_full[idx_bs:], phi_v)
    Z_lo2 = np.tile(z_full[idx_bs:], (len(phi_v), 1))
    ax2.plot_surface(R_lo2*np.cos(P_lo2), R_lo2*np.sin(P_lo2), Z_lo2, color="#e41a1c",
                     shade=True, lightsource=ls, edgecolor="#7f0000", linewidth=0.2)
    
    ax2.set_title(r"$\mathbf{(b)}$", loc="left", fontweight="bold")
    ax2.set_xlabel(r"$X$ ($\mu\mathrm{m}$)")
    ax2.set_ylabel(r"$Y$ ($\mu\mathrm{m}$)")
    ax2.set_zlabel(r"$Z$ ($\mu\mathrm{m}$)")
    ax2.view_init(elev=20, azim=40)
    ax2.set_box_aspect([1, 1, 1])
    
    plt.tight_layout()
    preview_path = "simulation_results/half_vesicle_preview.png"
    fig.savefig(preview_path, dpi=300)
    plt.close(fig)
    print(f"Saved preview figure to {preview_path}")


def main():
    profile_csv = "simulation_results_prl/two_domains_sweeps/profiles/prof_two_dom_DLo_3.0um_gamma_1.5pN_theta_0.20.csv"
    
    # 1. Equatorial Half (Cut along Z = 0)
    eq_obj = "simulation_results/vesicle_two_domain_half_equatorial.obj"
    eq_mtl = "simulation_results/vesicle_two_domain_half_equatorial.mtl"
    export_equatorial_half_obj(profile_csv, eq_obj, eq_mtl)
    
    # 2. Vertical Half (Cut along Y = 0 meridian plane)
    v_obj = "simulation_results/vesicle_two_domain_half_vertical.obj"
    v_mtl = "simulation_results/vesicle_two_domain_half_vertical.mtl"
    export_vertical_half_obj(profile_csv, v_obj, v_mtl)
    
    # 3. Preview render
    render_preview_figure()
    
    # Copy to Brain Artifacts
    brain_dir = r"C:\Users\nares\.gemini\antigravity\brain\f9a4e5d0-ff7a-425c-87a5-41ce15a21b6e"
    for f in [eq_obj, eq_mtl, v_obj, v_mtl, "simulation_results/half_vesicle_preview.png"]:
        dst = os.path.join(brain_dir, os.path.basename(f))
        shutil.copyfile(f, dst)
        print(f"Copied {os.path.basename(f)} to brain artifacts")

if __name__ == "__main__":
    main()
