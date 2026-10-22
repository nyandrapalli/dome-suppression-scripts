"""
Multi-Panel Cross-Section Figure for Two Domains at High Deflation (v = 0.90, theta = 1.0, gamma = 1.5 pN).
Compares Domain Diameters:
  Panel (a): D_Lo = 1.5 um (Shallow Stiff Caps, d = +0.31 um)
  Panel (b): D_Lo = 3.0 um (Optimal Deep Dual Inward Dimples, d = -0.59 um)
  Panel (c): D_Lo = 4.5 um (Broad Flat Domes, d = +0.22 um)
  Panel (d): D_Lo = 6.0 um (Large Spherical Domes, d = +0.52 um)
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

def plot_two_domains_v090_sweep(out_fig="simulation_results_prl/two_domains_sweeps/fig_two_domains_v090_DLo_sweep.png"):
    color_Lo = "#d73027"       # Vibrant red for Lo domains
    color_Ld = "#1a9850"       # Emerald green for equatorial Ld belt
    
    plt.rcParams.update({
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.titlesize": 12,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "legend.fontsize": 9,
        "figure.dpi": 300
    })
    
    fig = plt.figure(figsize=(18, 5.0))
    gs = GridSpec(1, 4, figure=fig, wspace=0.22)
    
    diameters = [1.5, 3.0, 4.5, 6.0]
    panel_labels = ["(a)", "(b)", "(c)", "(d)"]
    
    prof_dir = "simulation_results_prl/two_domains_sweeps/profiles"
    
    for i, d in enumerate(diameters):
        ax = fig.add_subplot(gs[0, i])
        candidates = [
            os.path.join(prof_dir, f"prof_two_dom_DLo_{d:.1f}_theta_1.00.csv"),
            os.path.join(prof_dir, f"prof_two_dom_v090_DLo_{d:.1f}um.csv"),
            os.path.join(prof_dir, f"prof_two_dom_DLo_{d:.1f}um_theta_1.00.csv"),
            os.path.join(prof_dir, "prof_two_dom_v_0.90_theta_1.00.csv") if d == 3.0 else ""
        ]
        prof_path = next((p for p in candidates if p and os.path.exists(p)), None)
            
        if os.path.exists(prof_path):
            df = pd.read_csv(prof_path)
            r = df["r_um"].values
            z = df["z_um"].values
            dom = df["domain"].values
            
            idx_b1 = np.where(dom == "Lo_north")[0][-1]
            idx_b2 = np.where(dom == "Lo_south")[0][0]
            
            z_mid = 0.5 * (np.max(z) + np.min(z))
            z_al = z - z_mid
            
            r_Ld_l = -r[idx_b1:idx_b2+1][::-1]
            z_Ld_l = z_al[idx_b1:idx_b2+1][::-1]
            r_Ld_r = r[idx_b1:idx_b2+1]
            z_Ld_r = z_al[idx_b1:idx_b2+1]
            
            r_Lo1 = np.concatenate([-r[:idx_b1+1][::-1], r[:idx_b1+1]])
            z_Lo1 = np.concatenate([z_al[:idx_b1+1][::-1], z_al[:idx_b1+1]])
            
            r_Lo2 = np.concatenate([-r[idx_b2:][::-1], r[idx_b2:]])
            z_Lo2 = np.concatenate([z_al[idx_b2:][::-1], z_al[idx_b2:]])
            
            ax.plot(r_Ld_l, z_Ld_l, color=color_Ld, linewidth=2.6, label=r"Equatorial $L_d$" if i == 0 else "")
            ax.plot(r_Ld_r, z_Ld_r, color=color_Ld, linewidth=2.6)
            ax.plot(r_Lo1, z_Lo1, color=color_Lo, linewidth=3.6, label=r"Polar $L_o$ Caps" if i == 0 else "")
            ax.plot(r_Lo2, z_Lo2, color=color_Lo, linewidth=3.6)
            
            ax.plot([-r[idx_b1], r[idx_b1]], [z_al[idx_b1], z_al[idx_b1]], "ko", markersize=5.0)
            ax.plot([-r[idx_b2], r[idx_b2]], [z_al[idx_b2], z_al[idx_b2]], "ko", markersize=5.0)
            
            depth = z_al[0] - z_al[idx_b1]
            rb = r[idx_b1]
            
            if d == 3.0:
                ann_color = "#7f0000"
                bg_color = "#fee0d2"
                edge_col = "#d73027"
                status_txt = r"$\mathbf{DEEP\ INWARD\ DIMPLE}$"
            else:
                ann_color = "#08519c"
                bg_color = "#f7fbff"
                edge_col = "#6baed6"
                status_txt = "Outward Dome" if depth > 0 else "Inward Dimple"
                
            ax.annotate(
                f"{status_txt}\nDepth: $\\mathbf{{{depth:+.2f}\\,\\mu\\mathrm{{m}}}}$\nNeck: $2r_b = {2*rb:.2f}\\,\\mu\\mathrm{{m}}$",
                xy=(0, 0), xytext=(0, 0),
                fontsize=9, ha="center", va="center", fontweight="bold", color=ann_color,
                bbox=dict(boxstyle="round,pad=0.35", facecolor=bg_color, edgecolor=edge_col, alpha=0.95)
            )
            
        ax.axvline(0, color="gray", linestyle=":", linewidth=0.8, alpha=0.7)
        ax.set_title(f"{panel_labels[i]} $D_{{Lo}} = {d:.1f}\\,\\mu\\mathrm{{m}}$ ($v=0.90$)", fontweight="bold", pad=8)
        ax.set_aspect("equal")
        ax.set_xlim(-6.0, 6.0)
        ax.set_ylim(-6.0, 6.0)
        ax.set_xlabel(r"Radial $r$ ($\mu\mathrm{m}$)")
        if i == 0:
            ax.set_ylabel(r"Axial $z$ ($\mu\mathrm{m}$)")
            ax.legend(frameon=True, facecolor="white", loc="lower right")
        ax.grid(True, linestyle=":", alpha=0.5)

    plt.tight_layout()
    fig.savefig(out_fig, dpi=300)
    plt.close(fig)
    print(f"Two-domain v=0.90 size sweep figure saved to: {out_fig}")
    return out_fig

if __name__ == "__main__":
    # First ensure profile files exist for all 4 diameters
    from sim_prl.mechanics_two_domains import TwoDomainParameters
    from sim_prl.solver_two_domains import TwoDomainVesicleSolver
    prof_dir = "simulation_results_prl/two_domains_sweeps/profiles"
    for d in [1.5, 3.0, 4.5, 6.0]:
        pfile = os.path.join(prof_dir, f"prof_two_dom_v090_DLo_{d:.1f}um.csv")
        params = TwoDomainParameters(R_ves=5.0, D_Lo=d, gamma_0=1.5e-12, reduced_volume=0.90)
        solver = TwoDomainVesicleSolver(params=params, n_points=101)
        res = solver.optimize_shape(theta=1.0, branch="inward", maxiter=400)
        
        df_p = pd.DataFrame({
            "s_um": solver.s_grid,
            "psi_rad": res["psi_profile"],
            "r_um": res["r_profile"],
            "z_um": res["z_profile"],
            "H_um_inv": res["H_profile"],
            "K_um_inv2": res["K_profile"]
        })
        dom_tags = []
        for s_val in solver.s_grid:
            if s_val <= solver.s_split1:
                dom_tags.append("Lo_north")
            elif s_val >= solver.s_split2:
                dom_tags.append("Lo_south")
            else:
                dom_tags.append("Ld")
        df_p["domain"] = dom_tags
        df_p.to_csv(pfile, index=False)
        
    plot_two_domains_v090_sweep()
