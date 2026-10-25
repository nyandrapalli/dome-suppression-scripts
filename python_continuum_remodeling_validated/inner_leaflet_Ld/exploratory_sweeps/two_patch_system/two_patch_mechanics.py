"""
Biophysical Mechanics Module for Two-Patch Antipodal Vesicles (North & South Lo caps, Equatorial Ld belt).
Implements Canham-Helfrich elasticity for 3 regions, dual boundary line tension, and ADE.
"""

import numpy as np

K_B_T = 4.114e-21  # J at 298.15 K

class TwoPatchMembraneParameters:
    def __init__(
        self,
        D_Lo=3.0,             # Diameter of each Lo patch in um
        R_ves=5.0,            # Vesicle radius in um
        kappa_Lo=80.0,        # in k_B T (stiff Lo)
        d_Lo=5.0e-9,          # 5.0 nm thick
        c0_Lo=0.0,            # spontaneous curvature of Lo
        kappa_Ld_0=20.0,      # bare soft Ld in k_B T
        kappa_prot_sat=70.0,  # additional stiffness from streptavidin (k_B T)
        d_Ld=4.0e-9,          # 4.0 nm thick
        c0_Ld_0=0.0,          # bare Ld spontaneous curvature
        alpha_c0=0.20,        # spontaneous curvature coupling in Ld belt (um^-1)
        gamma_0=1.5e-12,      # Line tension in N (1.5 pN) at each boundary
        k_ade_ratio=2.0,      # k_ade / kappa ratio
        eta_inner_expansion=0.08,
        reduced_volume=0.98
    ):
        self.D_Lo = D_Lo
        self.R_ves = R_ves
        self.R_Lo = D_Lo / 2.0
        self.kappa_Lo = kappa_Lo
        self.d_Lo = d_Lo
        self.c0_Lo = c0_Lo
        self.kappa_Ld_0 = kappa_Ld_0
        self.kappa_prot_sat = kappa_prot_sat
        self.d_Ld = d_Ld
        self.c0_Ld_0 = c0_Ld_0
        self.alpha_c0 = alpha_c0
        self.gamma_0 = gamma_0
        self.k_ade_ratio = k_ade_ratio
        self.eta_inner_expansion = eta_inner_expansion
        self.reduced_volume = reduced_volume
        
        # Target areas
        self.A_Lo_single = np.pi * (self.R_Lo)**2  # Area of each Lo patch
        self.A_tot_target = 4.0 * np.pi * (self.R_ves)**2  # ~ 314.16 um^2
        self.A_Lo_tot_target = 2.0 * self.A_Lo_single
        self.A_Ld_target = self.A_tot_target - self.A_Lo_tot_target
        self.V_tot_target = self.reduced_volume * (4.0 / 3.0) * np.pi * (self.R_ves)**3
        
    def get_effective_properties(self, theta):
        kappa_Ld = self.kappa_Ld_0 + self.kappa_prot_sat * theta
        c0_Ld = self.c0_Ld_0 + self.alpha_c0 * theta
        gamma_kBT_um = (self.gamma_0 / K_B_T) * 1e-6
        
        d_m_um = self.d_Ld * 1e6  # 0.004 um
        delta_A0_bare = 8.0 * np.pi * d_m_um * self.R_ves
        delta_A0_shift = - theta * self.eta_inner_expansion * (self.A_Ld_target / self.A_tot_target) * delta_A0_bare
        
        return {
            'kappa_Lo': self.kappa_Lo,
            'kappa_Ld': kappa_Ld,
            'c0_Lo': self.c0_Lo,
            'c0_Ld': c0_Ld,
            'gamma': gamma_kBT_um,
            'delta_A0_bare': delta_A0_bare,
            'delta_A0_shift': delta_A0_shift
        }


def reconstruct_profile_two_patch(s_grid, psi_grid, s_b1, s_b2):
    ds = np.diff(s_grid)
    cos_psi = np.cos(psi_grid)
    sin_psi = np.sin(psi_grid)
    
    r = np.zeros_like(s_grid)
    z = np.zeros_like(s_grid)
    
    r[1:] = np.cumsum(0.5 * (cos_psi[:-1] + cos_psi[1:]) * ds)
    z[1:] = np.cumsum(0.5 * (sin_psi[:-1] + sin_psi[1:]) * ds)
    
    idx_b1 = np.searchsorted(s_grid, s_b1)
    idx_b2 = np.searchsorted(s_grid, s_b2)
    return r, z, idx_b1, idx_b2


def compute_curvatures_two_patch(s_grid, psi_grid, r_grid):
    dpsi_ds = np.gradient(psi_grid, s_grid)
    c2 = np.zeros_like(psi_grid)
    small_r_mask = (r_grid < 1e-5)
    
    c2[~small_r_mask] = np.sin(psi_grid[~small_r_mask]) / r_grid[~small_r_mask]
    c2[small_r_mask] = dpsi_ds[small_r_mask]
    
    c1 = dpsi_ds
    H = 0.5 * (c1 + c2)
    K = c1 * c2
    return c1, c2, H, K


def compute_two_patch_energy(s_grid, psi_grid, s_b1, s_b2, params, theta):
    props = params.get_effective_properties(theta)
    r, z, idx_b1, idx_b2 = reconstruct_profile_two_patch(s_grid, psi_grid, s_b1, s_b2)
    c1, c2, H, K = compute_curvatures_two_patch(s_grid, psi_grid, r)
    
    ds = np.gradient(s_grid)
    dA = 2.0 * np.pi * r * ds
    dV = np.pi * (r**2) * np.sin(psi_grid) * ds
    
    # 3 Regions:
    # 1. North Lo: [0, idx_b1]
    # 2. Equatorial Ld: [idx_b1, idx_b2]
    # 3. South Lo: [idx_b2, end]
    mask_Lo1 = np.arange(len(s_grid)) <= idx_b1
    mask_Lo2 = np.arange(len(s_grid)) >= idx_b2
    mask_Ld = (~mask_Lo1) & (~mask_Lo2)
    
    A_Lo1 = np.sum(dA[mask_Lo1])
    A_Lo2 = np.sum(dA[mask_Lo2])
    A_Ld = np.sum(dA[mask_Ld])
    A_tot = A_Lo1 + A_Lo2 + A_Ld
    V_tot = np.sum(dV)
    
    # Bending energy densities
    f_bend_Lo1 = 0.5 * props['kappa_Lo'] * ((2.0 * H[mask_Lo1] - props['c0_Lo'])**2)
    f_bend_Lo2 = 0.5 * props['kappa_Lo'] * ((2.0 * H[mask_Lo2] - props['c0_Lo'])**2)
    f_bend_Ld = 0.5 * props['kappa_Ld'] * ((2.0 * H[mask_Ld] - props['c0_Ld'])**2)
    
    E_bend_Lo = np.sum(f_bend_Lo1 * dA[mask_Lo1]) + np.sum(f_bend_Lo2 * dA[mask_Lo2])
    E_bend_Ld = np.sum(f_bend_Ld * dA[mask_Ld])
    E_bend_tot = E_bend_Lo + E_bend_Ld
    
    # Dual boundary line tension: 2 * pi * r_b1 * gamma + 2 * pi * r_b2 * gamma
    r_b1 = r[idx_b1]
    r_b2 = r[idx_b2]
    E_line = (2.0 * np.pi * r_b1 + 2.0 * np.pi * r_b2) * props['gamma']
    
    # ADE Energy
    d_m_um = params.d_Ld * 1e6
    delta_A_geom = 2.0 * d_m_um * np.sum(H * dA)
    delta_A0_eff = props['delta_A0_bare'] + props['delta_A0_shift']
    
    k_ade = params.k_ade_ratio * props['kappa_Ld']
    E_ade = (np.pi * k_ade / (2.0 * A_tot * (d_m_um**2))) * ((delta_A_geom - delta_A0_eff)**2)
    
    total_physical_energy = E_bend_tot + E_line + E_ade
    
    # Residuals
    res_A_Lo1 = (A_Lo1 - params.A_Lo_single) / params.A_Lo_single
    res_A_Lo2 = (A_Lo2 - params.A_Lo_single) / params.A_Lo_single
    res_A_tot = (A_tot - params.A_tot_target) / params.A_tot_target
    res_V = (V_tot - params.V_tot_target) / params.V_tot_target
    res_pole = psi_grid[-1] - np.pi
    res_r_pole = r[-1]
    
    # Depths at both poles
    depth_north = z[idx_b1] - z[0]
    depth_south = z[-1] - z[idx_b2]
    
    return {
        'total_energy': total_physical_energy,
        'E_bend_Lo': E_bend_Lo,
        'E_bend_Ld': E_bend_Ld,
        'E_line': E_line,
        'E_ade': E_ade,
        'r_b1': r_b1,
        'r_b2': r_b2,
        'z_b1': z[idx_b1],
        'z_b2': z[idx_b2],
        'depth_north': depth_north,
        'depth_south': depth_south,
        'A_Lo1': A_Lo1,
        'A_Lo2': A_Lo2,
        'A_Ld': A_Ld,
        'A_tot': A_tot,
        'V_tot': V_tot,
        'res_A_Lo1': res_A_Lo1,
        'res_A_Lo2': res_A_Lo2,
        'res_A_tot': res_A_tot,
        'res_V': res_V,
        'res_pole': res_pole,
        'res_r_pole': res_r_pole,
        'r_profile': r,
        'z_profile': z,
        'H_profile': H,
        'idx_b1': idx_b1,
        'idx_b2': idx_b2,
        'props': props
    }
