"""
Vertical Cross-Section Visualizer for PRL Vesicle Remodeling.
Generates publication-quality 2D central vertical cross-sections showing:
- Lo domain in RED (#d73027)
- Ld matrix in GREEN (#1a9850)
- Domain boundary neck points (black dots)
- Central vertical symmetry axis
- Step-by-step morphological inversion across streptavidin coverage theta (0.0 to 1.0)
- Detailed zoomed comparison of bare vs remodeled states
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

def generate_vertical_crossections(results_dir="simulation_results_prl"):
    profiles_dir = os.path.join(results_dir, "profiles")
    out_path = os.path.join(results_dir, "vertical_cross_section_remodeling.png")
    
    color_Lo = "#d73027"  # Deep red for Lo domain
    color_Ld = "#1a9850"  # Forest green for Ld matrix
    
    plt.rcParams.update({
        "font.size": 9.5,
        "axes.labelsize": 10.5,
        "axes.titlesize": 11.5,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "figure.dpi": 300
    })
    
    thetas = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    
    # 2 rows x 4 cols layout:
    # 6 individual titration states + 1 Zoomed Apex comparison + 1 Full Overlay panel
    fig = plt.figure(figsize=(19, 10))
    gs = GridSpec(2, 4, figure=fig, hspace=0.28, wspace=0.25)
    
    profiles_history = []
    
    for i, theta in enumerate(thetas):
        row_idx = i // 3
        col_idx = i % 3
        ax = fig.add_subplot(gs[row_idx, col_idx])
        
        # Load profile
        branch = "sphere" if theta == 0.0 else "inward"
        prof_file = f"profile_theta_{theta:.2f}_{branch}.csv"
        prof_path = os.path.join(profiles_dir, prof_file)
        
        df = pd.read_csv(prof_path)
        r = df["r_um"].values
        z = df["z_um"].values
        dom = df["domain"].values
        
        mask_lo = (dom == "Lo")
        n_Lo = np.sum(mask_lo)
        
        # Realign z such that the vesicle center is at z = 0
        z_mid = 0.5 * (np.max(z) + np.min(z))
        z_aligned = z - z_mid
        
        z_apex = z_aligned[0]
        z_rim = z_aligned[n_Lo - 1]
        r_b_val = r[n_Lo - 1]
        depth_val = z_apex - z_rim  # Positive = outward dome, Negative = inward dimple
        
        profiles_history.append({
            "theta": theta,
            "branch": branch,
            "r": r,
            "z": z_aligned,
            "n_Lo": n_Lo,
            "r_b": r_b_val,
            "z_b": z_rim,
            "depth": depth_val
        })
        
        # Symmetric left and right halves
        r_Lo_symm = np.concatenate([-r[:n_Lo][::-1], r[:n_Lo]])
        z_Lo_symm = np.concatenate([z_aligned[:n_Lo][::-1], z_aligned[:n_Lo]])
        
        r_Ld_left = -r[n_Lo-1:][::-1]
        z_Ld_left = z_aligned[n_Lo-1:][::-1]
        r_Ld_right = r[n_Lo-1:]
        z_Ld_right = z_aligned[n_Lo-1:]
        
        # Plot Ld in green
        ax.plot(r_Ld_left, z_Ld_left, color=color_Ld, linewidth=2.6, label=r"$L_d$ Matrix (Green)")
        ax.plot(r_Ld_right, z_Ld_right, color=color_Ld, linewidth=2.6)
        
        # Plot Lo in red
        ax.plot(r_Lo_symm, z_Lo_symm, color=color_Lo, linewidth=3.6, label=r"$L_o$ Domain (Red)")
        
        # Mark domain boundary contact neck (black dots)
        ax.plot([-r_b_val, r_b_val], [z_rim, z_rim], "ko", markersize=4.5, zorder=5)
        
        # Symmetry axis
        ax.axvline(0, color="gray", linestyle=":", linewidth=0.9, alpha=0.7)
        
        # State annotation
        state_label = "Bare Prolate Dome" if theta == 0.0 else "Remodeled Dimple"
        depth_text = f"Outward Dome: +{depth_val:.2f} $\\mu\\mathrm{{m}}$" if depth_val > 0 else f"Inward Dimple: {abs(depth_val):.2f} $\\mu\\mathrm{{m}}$"
        callout_y = (z_apex - 1.7) if depth_val > 0 else (z_apex - 1.5)
        ax.annotate(
            f"$\\theta = {theta:.1f}$\n{depth_text}\n$2r_b = {2*r_b_val:.2f}\\,\\mu\\mathrm{{m}}$",
            xy=(0, z_apex), xytext=(0, callout_y),
            fontsize=8.5, ha="center", fontweight="bold",
            color="#990000" if depth_val <= 0 else "#004488",
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.9, edgecolor="#cccccc")
        )
        
        ax.set_title(f"({chr(97+i)}) $\\theta = {theta:.1f}$ ({state_label})", fontweight="bold", pad=6)
        ax.set_aspect("equal")
        ax.set_xlim(-6.0, 6.0)
        ax.set_ylim(-6.0, 6.0)
        ax.set_xlabel(r"Radial coordinate $r$ ($\mu\mathrm{m}$)")
        ax.set_ylabel(r"Axial coordinate $z$ ($\mu\mathrm{m}$)")
        ax.grid(True, linestyle=":", alpha=0.5)
        
        if i == 0:
            ax.legend(frameon=True, facecolor="white", framealpha=0.95, loc="lower right", fontsize=8)
            
    # -------------------------------------------------------------
    # Panel 7: Zoomed Lo Cap Inversion Comparison (Bare vs Remodeled)
    # -------------------------------------------------------------
    ax_zoom = fig.add_subplot(gs[0, 3])
    
    # Plot theta = 0.0 (bare) vs theta = 0.8 (remodeled)
    p_bare = profiles_history[0]
    p_remo = profiles_history[4]  # theta = 0.8
    
    # Bare profile
    n_b = p_bare["n_Lo"]
    r_b_symm = np.concatenate([-p_bare["r"][:n_b][::-1], p_bare["r"][:n_b]])
    z_b_symm = np.concatenate([p_bare["z"][:n_b][::-1], p_bare["z"][:n_b]])
    ax_zoom.plot(r_b_symm, z_b_symm, color="#4575b4", linestyle="--", linewidth=2.8, label=r"Bare Dome ($\theta=0.0$)")
    
    # Remodeled profile
    n_r = p_remo["n_Lo"]
    r_r_symm = np.concatenate([-p_remo["r"][:n_r][::-1], p_remo["r"][:n_r]])
    z_r_symm = np.concatenate([p_remo["z"][:n_r][::-1], p_remo["z"][:n_r]])
    ax_zoom.plot(r_r_symm, z_r_symm, color="#d73027", linestyle="-", linewidth=3.5, label=r"Remodeled Dimple ($\theta=0.8$)")
    
    # Boundary markers
    ax_zoom.plot([-p_bare["r_b"], p_bare["r_b"]], [p_bare["z_b"], p_bare["z_b"]], "o", color="#4575b4", markersize=5)
    ax_zoom.plot([-p_remo["r_b"], p_remo["r_b"]], [p_remo["z_b"], p_remo["z_b"]], "s", color="#d73027", markersize=5)
    
    ax_zoom.axvline(0, color="gray", linestyle=":", linewidth=0.9, alpha=0.7)
    ax_zoom.set_title(r"(g) Zoomed $L_o$ Apex Inversion", fontweight="bold", pad=6)
    ax_zoom.set_aspect("equal")
    ax_zoom.set_xlim(-2.5, 2.5)
    ax_zoom.set_ylim(2.8, 6.2)
    ax_zoom.set_xlabel(r"Radial coordinate $r$ ($\mu\mathrm{m}$)")
    ax_zoom.set_ylabel(r"Axial coordinate $z$ ($\mu\mathrm{m}$)")
    ax_zoom.grid(True, linestyle=":", alpha=0.5)
    ax_zoom.legend(frameon=True, facecolor="white", framealpha=0.95, loc="lower center", fontsize=8)
    
    # -------------------------------------------------------------
    # Panel 8: Composite Full-Vesicle Overlay
    # -------------------------------------------------------------
    ax_overlay = fig.add_subplot(gs[1, 3])
    cmap = matplotlib.colormaps["viridis"].resampled(len(profiles_history))
    
    for idx, p in enumerate(profiles_history):
        n_L = p["n_Lo"]
        r_L = np.concatenate([-p["r"][:n_L][::-1], p["r"][:n_L]])
        z_L = np.concatenate([p["z"][:n_L][::-1], p["z"][:n_L]])
        c = cmap(idx)
        ax_overlay.plot(r_L, z_L, color=c, linewidth=2.0, label=f"$\\theta = {p['theta']:.1f}$")
        
    # Plot outer Ld shell of final state
    p_last = profiles_history[4]
    n_last = p_last["n_Lo"]
    ax_overlay.plot(-p_last["r"][n_last-1:], p_last["z"][n_last-1:], color="#1a9850", linewidth=2.0, alpha=0.6)
    ax_overlay.plot(p_last["r"][n_last-1:], p_last["z"][n_last-1:], color="#1a9850", linewidth=2.0, alpha=0.6)
    
    ax_overlay.axvline(0, color="gray", linestyle=":", linewidth=0.9, alpha=0.7)
    ax_overlay.set_title(r"(h) Composite Trajectory Overlay", fontweight="bold", pad=6)
    ax_overlay.set_aspect("equal")
    ax_overlay.set_xlim(-6.2, 6.2)
    ax_overlay.set_ylim(-6.2, 6.2)
    ax_overlay.set_xlabel(r"Radial coordinate $r$ ($\mu\mathrm{m}$)")
    ax_overlay.set_ylabel(r"Axial coordinate $z$ ($\mu\mathrm{m}$)")
    ax_overlay.grid(True, linestyle=":", alpha=0.5)
    ax_overlay.legend(frameon=True, facecolor="white", framealpha=0.9, loc="lower right", fontsize=7.5, ncol=2)
    
    plt.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"Vertical cross-sections saved successfully to: {out_path}")
    return out_path

if __name__ == "__main__":
    generate_vertical_crossections()
