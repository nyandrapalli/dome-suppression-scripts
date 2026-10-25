"""
Export Watertight 3D Wavefront OBJ and MTL Files for Dimpled Vesicles:
1. Single-Domain Dimpled Vesicle (Lo in Red, Ld in Green)
2. Two-Domain Biconcave Discocyte Vesicle (Dual Lo in Red, Equatorial Ld in Green)

Includes both standard MTL material references (usemtl) and vertex colors (v x y z r g b)
for maximum compatibility across Blender, MeshLab, 3D Viewer, ParaView, and CAD software.
"""

import os
import shutil
import numpy as np
import pandas as pd

def generate_vesicle_obj(
    profile_csv,
    domain_type, # 'single' or 'two'
    out_obj_path,
    out_mtl_path,
    n_phi=96
):
    os.makedirs(os.path.dirname(out_obj_path), exist_ok=True)
    
    # -------------------------------------------------------------------------
    # 1. Load and Clean Profile
    # -------------------------------------------------------------------------
    df = pd.read_csv(profile_csv)
    
    r = df["r_um"].values.copy()
    z = df["z_um"].values.copy()
    s = df["s_um"].values.copy()
    
    # Center equator at z = 0
    z_mid = (z.max() + z.min()) / 2.0
    z -= z_mid
    
    # Ensure exact apex pole closure
    r[0] = 0.0
    r[-1] = 0.0
    r = np.abs(r) # radial coordinate strictly non-negative
    
    n_pts = len(r)
    
    # Domain boundary segmentation
    # Single domain: s <= 1.50 um is Lo (Red), remainder is Ld (Green)
    # Two domain: s <= 1.50 um is Lo (Red), s >= s_tot - 1.50 um is Lo (Red), middle is Ld (Green)
    s_tot = s[-1]
    is_Lo = np.zeros(n_pts, dtype=bool)
    if domain_type == "single":
        is_Lo[s <= 1.50] = True
    else: # two domain
        is_Lo[(s <= 1.50) | (s >= s_tot - 1.50)] = True

    # Compute tangent angles psi and outward normals
    # ds = diff(s), dr = diff(r), dz = diff(z)
    dr = np.gradient(r, s)
    dz = np.gradient(z, s)
    psi = np.arctan2(dz, dr) # profile slope
    
    # -------------------------------------------------------------------------
    # 2. Write Material File (.mtl)
    # -------------------------------------------------------------------------
    mtl_filename = os.path.basename(out_mtl_path)
    with open(out_mtl_path, "w", encoding="utf-8") as fmtl:
        fmtl.write("# Material definitions for Dimpled Giant Unilamellar Vesicle\n\n")
        
        # Lo Material: Vibrant Red
        fmtl.write("newmtl Lo_Red\n")
        fmtl.write("Ka 0.15 0.02 0.02\n")
        fmtl.write("Kd 0.894 0.102 0.110\n")  # #e41a1c
        fmtl.write("Ks 0.40 0.40 0.40\n")
        fmtl.write("Ns 40.0\n")
        fmtl.write("d 1.0\n")
        fmtl.write("illum 2\n\n")
        
        # Ld Material: Vibrant Emerald Green
        fmtl.write("newmtl Ld_Green\n")
        fmtl.write("Ka 0.02 0.15 0.05\n")
        fmtl.write("Kd 0.255 0.671 0.365\n")  # #41ab5d
        fmtl.write("Ks 0.35 0.35 0.35\n")
        fmtl.write("Ns 35.0\n")
        fmtl.write("d 1.0\n")
        fmtl.write("illum 2\n")

    # -------------------------------------------------------------------------
    # 3. Generate 3D Vertex Grid & Normals
    # -------------------------------------------------------------------------
    phi = np.linspace(0, 2 * np.pi, n_phi, endpoint=False)
    cos_phi = np.cos(phi)
    sin_phi = np.sin(phi)
    
    # Color palette
    color_Lo_rgb = (0.894, 0.102, 0.110)
    color_Ld_rgb = (0.255, 0.671, 0.365)
    
    vertices = []       # [(x, y, z, r, g, b), ...]
    normals = []        # [(nx, ny, nz), ...]
    ring_start_idx = [] # index in vertex list of the first vertex of ring i
    
    # Ring 0: North pole (single vertex)
    c0 = color_Lo_rgb if is_Lo[0] else color_Ld_rgb
    vertices.append((0.0, 0.0, z[0], *c0))
    normals.append((0.0, 0.0, 1.0))
    ring_start_idx.append(0)
    
    # Intermediate rings: i = 1 to n_pts - 2
    for i in range(1, n_pts - 1):
        ring_start_idx.append(len(vertices))
        ri = r[i]
        zi = z[i]
        ci = color_Lo_rgb if is_Lo[i] else color_Ld_rgb
        
        # Normal components: n = (-sin(psi)*cos(phi), -sin(psi)*sin(phi), cos(psi))
        # Note: dz/ds = sin(psi), dr/ds = cos(psi). Outward normal for clockwise profile:
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
            
    # Final Ring: South pole (single vertex)
    ring_start_idx.append(len(vertices))
    cN = color_Lo_rgb if is_Lo[-1] else color_Ld_rgb
    vertices.append((0.0, 0.0, z[-1], *cN))
    normals.append((0.0, 0.0, -1.0))
    
    # -------------------------------------------------------------------------
    # 4. Generate Triangular Faces Grouped by Material
    # -------------------------------------------------------------------------
    faces_Lo = []
    faces_Ld = []
    
    # North cap fan: connect North Pole (vertex 1, 1-based indexing) to ring 1
    # Ring 1 has n_phi vertices starting at index ring_start_idx[1]
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
            
    # Intermediate quad bands between ring i and ring i + 1
    for i in range(1, n_pts - 2):
        curr_start = ring_start_idx[i]
        next_start = ring_start_idx[i + 1]
        
        # Band material: Lo if either ring is Lo
        band_Lo = is_Lo[i] and is_Lo[i + 1]
        target_list = faces_Lo if band_Lo else faces_Ld
        
        for j in range(n_phi):
            j_next = (j + 1) % n_phi
            
            # 1-based vertex indices
            v_curr_j = curr_start + j + 1
            v_curr_next = curr_start + j_next + 1
            v_next_j = next_start + j + 1
            v_next_next = next_start + j_next + 1
            
            # Two triangles forming quad with consistent outward winding
            target_list.append((v_curr_j, v_next_j, v_next_next))
            target_list.append((v_curr_j, v_next_next, v_curr_next))
            
    # South cap fan: connect last ring (n_pts - 2) to South Pole (last vertex)
    last_ring_start = ring_start_idx[n_pts - 2]
    south_pole_idx = len(vertices) # 1-based index
    is_south_Lo = is_Lo[-1] or is_Lo[-2]
    for j in range(n_phi):
        j_next = (j + 1) % n_phi
        v1 = last_ring_start + j + 1
        v2 = last_ring_start + j_next + 1
        face = (v1, south_pole_idx, v2)
        if is_south_Lo:
            faces_Lo.append(face)
        else:
            faces_Ld.append(face)
            
    # -------------------------------------------------------------------------
    # 5. Write OBJ File
    # -------------------------------------------------------------------------
    with open(out_obj_path, "w", encoding="utf-8") as fobj:
        fobj.write(f"# Wavefront OBJ: Dimpled Vesicle ({domain_type} domain)\n")
        fobj.write(f"# Vertices: {len(vertices)}, Faces: {len(faces_Lo) + len(faces_Ld)}\n")
        fobj.write(f"mtllib {mtl_filename}\n\n")
        
        # Vertices with vertex colors
        for v in vertices:
            fobj.write(f"v {v[0]:.6f} {v[1]:.6f} {v[2]:.6f} {v[3]:.4f} {v[4]:.4f} {v[5]:.4f}\n")
            
        fobj.write("\n")
        # Normals
        for n in normals:
            fobj.write(f"vn {n[0]:.6f} {n[1]:.6f} {n[2]:.6f}\n")
            
        fobj.write("\n")
        
        # Lo Domain Group
        if len(faces_Lo) > 0:
            fobj.write("g Lo_Domain\n")
            fobj.write("usemtl Lo_Red\n")
            for f in faces_Lo:
                fobj.write(f"f {f[0]}//{f[0]} {f[1]}//{f[1]} {f[2]}//{f[2]}\n")
            fobj.write("\n")
            
        # Ld Membrane Group
        if len(faces_Ld) > 0:
            fobj.write("g Ld_Membrane\n")
            fobj.write("usemtl Ld_Green\n")
            for f in faces_Ld:
                fobj.write(f"f {f[0]}//{f[0]} {f[1]}//{f[1]} {f[2]}//{f[2]}\n")
            fobj.write("\n")
            
    print(f"Successfully generated {out_obj_path} ({len(vertices)} vertices, {len(faces_Lo) + len(faces_Ld)} faces)")

def main():
    f_s = "simulation_results_prl/param_sweeps/profiles/prof_DLo_3.0um_gamma_1.5pN_theta_0.20.csv"
    f_t = "simulation_results_prl/two_domains_sweeps/profiles/prof_two_dom_DLo_3.0um_gamma_1.5pN_theta_0.20.csv"
    
    # 1. Single Domain
    obj_s = "simulation_results/vesicle_single_domain_dimple.obj"
    mtl_s = "simulation_results/vesicle_single_domain_dimple.mtl"
    generate_vesicle_obj(f_s, "single", obj_s, mtl_s)
    
    # 2. Two Domain
    obj_t = "simulation_results/vesicle_two_domain_dimple.obj"
    mtl_t = "simulation_results/vesicle_two_domain_dimple.mtl"
    generate_vesicle_obj(f_t, "two", obj_t, mtl_t)
    
    # Copy to Brain Artifacts folder for easy access
    brain_dir = r"C:\Users\nares\.gemini\antigravity\brain\f9a4e5d0-ff7a-425c-87a5-41ce15a21b6e"
    for f in [obj_s, mtl_s, obj_t, mtl_t]:
        dst = os.path.join(brain_dir, os.path.basename(f))
        shutil.copyfile(f, dst)
        print(f"Copied {os.path.basename(f)} to brain artifacts directory")

if __name__ == "__main__":
    main()
