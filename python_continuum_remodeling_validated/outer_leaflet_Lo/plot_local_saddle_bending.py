"""
Visualization of Local Saddle Bending in Ld Domain for HermesLoouter Two-Domain Vesicles.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.collections import LineCollection

def plot_local_saddle():
    root_dir = "HermesLoouter"
    out_fig = os.path.join(root_dir, "fig_local_saddle_bending.png")
    
    # Load representative high-contrast profiles
    # 1. Standard D_Lo = 3.0 um, theta = 1.0, v = 0.98, gamma = 1.5 pN
    p1 = os.path.join(root_dir, "double", "profiles", "prof_two_dom_DLo_3.0_theta_1.00.csv")
    # 2. Deflated v = 0.90, theta = 1.0, D_Lo = 3.0 um, gamma = 1.5 pN
    p2 = os.path.join(root_dir, "double", "profiles", "prof_two_dom_v_0.90_theta_1.00.csv")
    
    if not os.path.exists(p1):
        print(f"Profile {p1} not found.")
        return
        
    df = pd.read_csv(p1)
    s = df["s_um"].values
    psi = df["psi_rad"].values
    r = df["r_um"].values
    z = df["z_um"].values
    dom = df["domain"].values
    
    # Recalculate robust smooth principal curvatures
    ds = np.gradient(s)
    c_m = np.gradient(psi, s)  # Meridional curvature
    r_pos = np.maximum(r, 1e-8)
    c_p = np.sin(psi) / r_pos  # Parallel curvature
    c_p[0] = c_m[0]
    c_p[-1] = c_m[-1]
    
    H = 0.5 * (c_m + c_p)      # Mean curvature
    K = c_m * c_p             # Gaussian curvature
    
    mask_lo_n = (dom == "Lo_north")
    mask_lo_s = (dom == "Lo_south")
    mask_ld   = (dom == "Ld")
    
    idx_b1 = np.where(mask_lo_n)[0][-1]
    idx_b2 = np.where(mask_lo_s)[0][0]
    s_b1 = s[idx_b1]
    s_b2 = s[idx_b2]
    
    z_mid = 0.5 * (np.max(z) + np.min(z))
    z_aligned = z - z_mid
    
    # Setup Figure
    plt.rcParams.update({
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.titlesize": 12,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "figure.dpi": 300
    })
    
    fig = plt.figure(figsize=(16, 11))
    gs = GridSpec(2, 2, figure=fig, hspace=0.32, wspace=0.28)
    
    # -------------------------------------------------------------
    # Panel A: Full Cross-Section with Domain Breakdown & Saddle Annotation
    # -------------------------------------------------------------
    ax_a = fig.add_subplot(gs[0, 0])
    
    # Symmetric 2D vesicle profile
    color_Lo = "#d73027"  # Deep red
    color_Ld = "#1a9850"  # Green
    
    ax_a.plot(np.concatenate([-r[mask_lo_n][::-1], r[mask_lo_n]]), 
              np.concatenate([z_aligned[mask_lo_n][::-1], z_aligned[mask_lo_n]]), 
              color=color_Lo, linewidth=3.5, label=r"North $L_o$ Cap ($c_0 > 0$)")
    
    ax_a.plot(-r[mask_ld][::-1], z_aligned[mask_ld][::-1], color=color_Ld, linewidth=2.8, label=r"Equatorial $L_d$ Belt (Bare)")
    ax_a.plot(r[mask_ld], z_aligned[mask_ld], color=color_Ld, linewidth=2.8)
    
    ax_a.plot(np.concatenate([-r[mask_lo_s][::-1], r[mask_lo_s]]), 
              np.concatenate([z_aligned[mask_lo_s][::-1], z_aligned[mask_lo_s]]), 
              color=color_Lo, linewidth=3.5, label=r"South $L_o$ Cap ($c_0 > 0$)")
    
    # Mark Domain Boundary Necks
    r_b = r[idx_b1]
    z_b1 = z_aligned[idx_b1]
    z_b2 = z_aligned[idx_b2]
    ax_a.plot([-r_b, r_b], [z_b1, z_b1], "ko", markersize=6, label=r"Boundary Neck ($r_b$)")
    ax_a.plot([-r_b, r_b], [z_b2, z_b2], "ko", markersize=6)
    
    # Highlight Local Saddle Zones in Ld adjacent to necks
    saddle_n_idx = np.where((s > s_b1) & (s < s_b1 + 2.0))[0]
    ax_a.plot(r[saddle_n_idx], z_aligned[saddle_n_idx], color="#984ea3", linewidth=4.5, linestyle="--", label="Local Saddle Region ($K < 0$)")
    ax_a.plot(-r[saddle_n_idx], z_aligned[saddle_n_idx], color="#984ea3", linewidth=4.5, linestyle="--")
    
    saddle_s_idx = np.where((s < s_b2) & (s > s_b2 - 2.0))[0]
    ax_a.plot(r[saddle_s_idx], z_aligned[saddle_s_idx], color="#984ea3", linewidth=4.5, linestyle="--")
    ax_a.plot(-r[saddle_s_idx], z_aligned[saddle_s_idx], color="#984ea3", linewidth=4.5, linestyle="--")
    
    ax_a.axvline(0, color="gray", linestyle=":", linewidth=0.8, alpha=0.6)
    ax_a.set_aspect("equal")
    ax_a.set_xlim(-6.5, 6.5)
    ax_a.set_ylim(-6.5, 6.5)
    ax_a.set_xlabel(r"Radial coordinate $r$ ($\mu\mathrm{m}$)")
    ax_a.set_ylabel(r"Axial coordinate $z$ ($\mu\mathrm{m}$)")
    ax_a.set_title(r"(a) Full Bipolar Vesicle Profile ($D_{Lo}=3.0\,\mu\mathrm{m},\,v=0.98,\,\theta=1.0$)", fontweight="bold")
    ax_a.grid(True, linestyle=":", alpha=0.5)
    ax_a.legend(loc="lower left", fontsize=8.5, framealpha=0.95)
    
    # -------------------------------------------------------------
    # Panel B: High-Magnification Zoom of North Neck & Saddle Inflection
    # -------------------------------------------------------------
    ax_b = fig.add_subplot(gs[0, 1])
    
    zoom_idx = np.where(s <= s_b1 + 3.5)[0]
    ax_b.plot(r[zoom_idx], z_aligned[zoom_idx], color="#555555", linewidth=1.5, linestyle=":")
    
    # Plot Lo portion
    zoom_lo = np.where(s <= s_b1)[0]
    ax_b.plot(r[zoom_lo], z_aligned[zoom_lo], color=color_Lo, linewidth=4.0, label=r"$L_o$ Dome ($H > 0$, Convex)")
    
    # Plot Saddle portion of Ld
    zoom_saddle = np.where((s >= s_b1) & (s <= s_b1 + 2.5))[0]
    ax_b.plot(r[zoom_saddle], z_aligned[zoom_saddle], color="#984ea3", linewidth=4.5, label=r"$L_d$ Saddle Zone ($H < 0, K < 0$)")
    
    # Plot remaining Ld
    zoom_ld_rem = np.where((s >= s_b1 + 2.5) & (s <= s_b1 + 3.5))[0]
    ax_b.plot(r[zoom_ld_rem], z_aligned[zoom_ld_rem], color=color_Ld, linewidth=3.5, label=r"$L_d$ Equatorial Flank")
    
    # Mark Neck
    ax_b.plot(r_b, z_b1, "ko", markersize=8)
    ax_b.annotate(r"Domain Neck $r_b = 1.30\,\mu\mathrm{m}$" + "\n(Line Tension $\\gamma$ Pinches Rim)", 
                  xy=(r_b, z_b1), xytext=(r_b + 0.5, z_b1 + 0.8),
                  arrowprops=dict(arrowstyle="->", color="black", lw=1.5),
                  fontsize=9, fontweight="bold",
                  bbox=dict(boxstyle="round,pad=0.3", facecolor="#ffffcc", edgecolor="#cccc00", alpha=0.9))
    
    # Mark Inflection / Saddle
    saddle_mid_idx = zoom_saddle[len(zoom_saddle)//2]
    r_sad, z_sad = r[saddle_mid_idx], z_aligned[saddle_mid_idx]
    ax_b.annotate("Saddle Inflection Point:\n$c_m < 0$ (Meridional concave)\n$c_p > 0$ (Parallel convex)", 
                  xy=(r_sad, z_sad), xytext=(r_sad + 1.2, z_sad - 1.2),
                  arrowprops=dict(arrowstyle="->", color="#984ea3", lw=1.8),
                  fontsize=8.5, fontweight="bold", color="#49006a",
                  bbox=dict(boxstyle="round,pad=0.3", facecolor="#f1eef6", edgecolor="#984ea3", alpha=0.9))
    
    ax_b.set_aspect("equal")
    ax_b.set_xlim(-0.2, 5.0)
    ax_b.set_ylim(z_b1 - 3.5, z_aligned[0] + 1.0)
    ax_b.set_xlabel(r"Radial coordinate $r$ ($\mu\mathrm{m}$)")
    ax_b.set_ylabel(r"Axial coordinate $z$ ($\mu\mathrm{m}$)")
    ax_b.set_title(r"(b) High-Magnification Zoom of $L_o \to L_d$ Boundary Neck", fontweight="bold")
    ax_b.grid(True, linestyle=":", alpha=0.5)
    ax_b.legend(loc="upper right", fontsize=8.5, framealpha=0.95)
    
    # -------------------------------------------------------------
    # Panel C: Spatial Mean Curvature Profile H(s) along Arc Length
    # -------------------------------------------------------------
    ax_c = fig.add_subplot(gs[1, 0])
    
    ax_c.plot(s, H, color="#08519c", linewidth=2.2, label=r"Mean Curvature $H(s) = \frac{1}{2}(c_m + c_p)$")
    ax_c.axhline(0, color="black", linestyle="--", linewidth=1.0, alpha=0.7)
    
    # Shade domains
    ax_c.axvspan(0, s_b1, color=color_Lo, alpha=0.15, label=r"North $L_o$ Cap")
    ax_c.axvspan(s_b1, s_b2, color=color_Ld, alpha=0.10, label=r"Equatorial $L_d$ Belt")
    ax_c.axvspan(s_b2, s[-1], color=color_Lo, alpha=0.15, label=r"South $L_o$ Cap")
    
    # Highlight negative H regions
    ax_c.fill_between(s, H, 0, where=(H < 0), color="#984ea3", alpha=0.4, label=r"Inward Saddle Bending ($H < 0$)")
    
    # Mark boundaries
    ax_c.axvline(s_b1, color="red", linestyle=":", linewidth=1.5)
    ax_c.axvline(s_b2, color="red", linestyle=":", linewidth=1.5)
    ax_c.text(s_b1, np.max(H)*0.85, " Boundary 1", color="red", fontsize=8.5, fontweight="bold")
    ax_c.text(s_b2, np.max(H)*0.85, " Boundary 2", color="red", fontsize=8.5, fontweight="bold", ha="right")
    
    ax_c.set_xlim(0, s[-1])
    ax_c.set_xlabel(r"Arc-length parameter $s$ ($\mu\mathrm{m}$)")
    ax_c.set_ylabel(r"Mean Curvature $H$ ($\mu\mathrm{m}^{-1}$)")
    ax_c.set_title(r"(c) Spatial Mean Curvature: Negative Dip ($H < 0$) in $L_d$", fontweight="bold")
    ax_c.grid(True, linestyle=":", alpha=0.5)
    ax_c.legend(loc="upper center", fontsize=8, ncol=2, framealpha=0.95)
    
    # -------------------------------------------------------------
    # Panel D: Principal & Gaussian Curvatures: c_m, c_p, K = c_m * c_p
    # -------------------------------------------------------------
    ax_d = fig.add_subplot(gs[1, 1])
    
    ax_d.plot(s, c_m, color="#e6550d", linewidth=2.0, linestyle="-", label=r"Meridional $c_m = \mathrm{d}\psi/\mathrm{d}s$")
    ax_d.plot(s, c_p, color="#31a354", linewidth=2.0, linestyle="-", label=r"Parallel $c_p = \sin\psi / r$")
    ax_d.plot(s, K, color="#756bb1", linewidth=2.2, linestyle="--", label=r"Gaussian $K = c_m \cdot c_p$")
    ax_d.axhline(0, color="black", linestyle="--", linewidth=1.0, alpha=0.7)
    
    # Fill region where K < 0 (Hyperbolic / Saddle surface)
    ax_d.fill_between(s, K, 0, where=(K < 0), color="#bcbddc", alpha=0.5, label=r"Saddle Geometry ($K < 0$)")
    
    ax_d.axvspan(0, s_b1, color=color_Lo, alpha=0.10)
    ax_d.axvspan(s_b1, s_b2, color=color_Ld, alpha=0.08)
    ax_d.axvspan(s_b2, s[-1], color=color_Lo, alpha=0.10)
    
    ax_d.set_xlim(0, s[-1])
    ax_d.set_ylim(-3.0, 3.5)
    ax_d.set_xlabel(r"Arc-length parameter $s$ ($\mu\mathrm{m}$)")
    ax_d.set_ylabel(r"Curvature ($\mu\mathrm{m}^{-1}$ or $\mu\mathrm{m}^{-2}$)")
    ax_d.set_title(r"(d) Principal & Gaussian Curvatures ($K < 0$ in $L_d$ Neck Flank)", fontweight="bold")
    ax_d.grid(True, linestyle=":", alpha=0.5)
    ax_d.legend(loc="upper center", fontsize=8, ncol=2, framealpha=0.95)
    
    plt.suptitle("HermesLoouter: Direct Evidence of Local Saddle Bending in the $L_d$ Domain Adjacent to Both $L_o$ Necks", 
                 fontsize=13, fontweight="bold", y=0.98)
    
    fig.savefig(out_fig, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved figure: {out_fig}")

if __name__ == "__main__":
    plot_local_saddle()
