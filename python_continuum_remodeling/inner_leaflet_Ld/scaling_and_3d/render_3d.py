"""
3D Surface Rendering and Wavefront OBJ Exporter for Two-Phase Vesicles (PRL Standard).
Generates publication-quality 3D renders and exported 3D meshes (.obj / .mtl)
for bare vs remodeled dimpled vesicle geometries.
"""

import os
import sys
import math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource, to_rgba
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

LD_COLOR = "#2b83ba"  # Clean blue for Ld domain
LO_COLOR = "#d7191c"  # Vibrant red/coral for Lo domain

def mesh_from_profile(r_prof, z_prof, domains, n_azimuth=64):
    """
    Constructs a 3D surface mesh of vertices and faces from axisymmetric profile (r, z).
    """
    phi = np.linspace(0.0, 2.0 * np.pi, n_azimuth, endpoint=False)
    vertices = [[0.0, 0.0, float(z_prof[0])]]
    
    for i in range(1, len(r_prof) - 1):
        for angle in phi:
            vertices.append([
                float(r_prof[i] * math.cos(angle)),
                float(r_prof[i] * math.sin(angle)),
                float(z_prof[i])
            ])
    south_idx = len(vertices)
    vertices.append([0.0, 0.0, float(z_prof[-1])])
    
    def ring_start(idx):
        return 1 + (idx - 1) * n_azimuth
        
    faces = []
    face_domains = []
    
    # North pole cap
    first_ring = ring_start(1)
    for a in range(n_azimuth):
        nxt = (a + 1) % n_azimuth
        faces.append([0, first_ring + a, first_ring + nxt])
        face_domains.append(domains[1])
        
    # Body quads
    for i in range(1, len(r_prof) - 2):
        upper = ring_start(i)
        lower = ring_start(i + 1)
        for a in range(n_azimuth):
            nxt = (a + 1) % n_azimuth
            faces.append([upper + a, upper + nxt, lower + nxt, lower + a])
            face_domains.append(domains[i + 1])
            
    # South pole cap
    last_ring = ring_start(len(r_prof) - 2)
    for a in range(n_azimuth):
        nxt = (a + 1) % n_azimuth
        faces.append([last_ring + a, last_ring + nxt, south_idx])
        face_domains.append(domains[-1])
        
    return np.array(vertices), faces, face_domains

def export_obj_mesh(vertices, faces, face_domains, obj_path, mtl_path):
    """
    Exports a 3D Wavefront OBJ and MTL mesh with separate Lo and Ld material groups.
    """
    with open(mtl_path, "w", encoding="ascii") as mtl:
        mtl.write("newmtl Lo\nKd 0.843 0.098 0.110\nKa 0.1 0.02 0.02\nKs 0.15 0.15 0.15\nNs 40\n\n")
        mtl.write("newmtl Ld\nKd 0.169 0.514 0.729\nKa 0.02 0.08 0.12\nKs 0.15 0.15 0.15\nNs 40\n")
        
    with open(obj_path, "w", encoding="ascii") as obj:
        obj.write(f"# PRL Axisymmetric Vesicle Mesh\nmtllib {os.path.basename(mtl_path)}\n")
        for v in vertices:
            obj.write(f"v {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n")
        current_dom = None
        for face, dom in zip(faces, face_domains):
            if dom != current_dom:
                obj.write(f"g {dom}\nusemtl {dom}\n")
                current_dom = dom
            indices = " ".join(str(idx + 1) for idx in face)
            obj.write(f"f {indices}\n")

def render_3d_vesicle(profile_csv, output_png, title="3D Vesicle Morphology"):
    df = pd.read_csv(profile_csv)
    r = df["r_um"].values
    z = df["z_um"].values
    domains = df["domain"].values
    
    sub_idx = np.linspace(0, len(r) - 1, 80, dtype=int)
    r_sub = r[sub_idx]
    z_sub = z[sub_idx]
    dom_sub = domains[sub_idx]
    
    vertices, faces, face_domains = mesh_from_profile(r_sub, z_sub, dom_sub, n_azimuth=48)
    
    polygons = [[vertices[idx] for idx in face] for face in faces]
    colors = [to_rgba(LO_COLOR if dom == "Lo" else LD_COLOR, alpha=0.95) for dom in face_domains]
    
    fig = plt.figure(figsize=(7, 7), dpi=250)
    ax = fig.add_subplot(111, projection="3d")
    
    light = LightSource(azdeg=315, altdeg=45)
    poly_col = Poly3DCollection(
        polygons,
        facecolors=colors,
        edgecolors=colors,
        linewidths=0.0,
        shade=True,
        lightsource=light
    )
    ax.add_collection3d(poly_col)
    
    limit = 5.8
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_zlim(-limit, limit)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=22, azim=-60)
    ax.set_axis_off()
    ax.set_title(title, fontsize=12, fontweight="bold", pad=-15)
    
    fig.savefig(output_png, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Rendered 3D image saved to {output_png}")
    
    # Export OBJ
    obj_path = output_png.replace(".png", ".obj")
    mtl_path = output_png.replace(".png", ".mtl")
    export_obj_mesh(vertices, faces, face_domains, obj_path, mtl_path)
    print(f"Exported OBJ mesh: {obj_path}")

if __name__ == "__main__":
    prof_inward = "simulation_results_prl/profiles/profile_theta_0.90_inward.csv"
    prof_sphere = "simulation_results_prl/profiles/profile_theta_0.00_sphere.csv"
    
    if os.path.exists(prof_inward):
        render_3d_vesicle(prof_inward, "simulation_results_prl/3d_dimple_prl.png", "Remodeled Inward Dimple (theta = 0.9)")
    if os.path.exists(prof_sphere):
        render_3d_vesicle(prof_sphere, "simulation_results_prl/3d_bare_sphere_prl.png", "Bare Spherical Vesicle (theta = 0.0)")
