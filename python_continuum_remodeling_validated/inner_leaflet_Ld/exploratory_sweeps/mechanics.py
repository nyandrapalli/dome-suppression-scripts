"""
Biophysical Mechanics Module for Two-Phase Vesicle Remodeling
Implements Canham-Helfrich-Jülicher-Lipowsky elasticity, Area-Difference Elasticity (ADE),
and domain boundary line tension.
"""

import numpy as np

# Physical constants in SI / biophysical units
K_B_T = 4.114e-21  # J at 298.15 K (1 k_B T in Joules)

class MembraneParameters:
    def __init__(
        self,
        # Lo domain properties
        kappa_Lo=80.0,        # in k_B T (stiff Lo)
        d_Lo=5.0e-9,          # 5.0 nm thick
        c0_Lo=0.0,            # spontaneous curvature of bare Lo in um^-1
        # Ld domain properties (bare)
        kappa_Ld_0=20.0,      # in k_B T (soft Ld)
        kappa_prot_sat=70.0,  # additional stiffness from saturated streptavidin layer (k_B T -> up to 90 k_B T)
        d_Ld=4.0e-9,          # 4.0 nm thick
        c0_Ld_0=0.0,          # bare Ld spontaneous curvature
        alpha_c0=0.20,        # spontaneous curvature coupling (um^-1) matching 1/R_ves
        # Line tension and ADE
        gamma_0=1.5e-12,      # Line tension in N (1.5 pN)
        k_ade_ratio=2.0,      # k_ade / kappa ratio
        eta_inner_expansion=0.08, # 8% inner leaflet expansion at saturation
        # Geometric scale
        R_ves=5.0,            # Vesicle radius in um
        D_Lo=3.0,             # Lo domain diameter in um (radius = 1.5 um)
        reduced_volume=0.98   # Deflation / reduced volume v = V / V_sphere
    ):
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
        self.R_ves = R_ves
        self.D_Lo = D_Lo
        self.R_Lo = D_Lo / 2.0
        self.reduced_volume = reduced_volume
        
        # Computed target areas
        self.A_Lo_target = np.pi * (self.R_Lo)**2  # ~ 7.068 um^2
        self.A_tot_target = 4.0 * np.pi * (self.R_ves)**2  # ~ 314.16 um^2
        self.A_Ld_target = self.A_tot_target - self.A_Lo_target
        self.V_tot_target = self.reduced_volume * (4.0 / 3.0) * np.pi * (self.R_ves)**3
        
    def get_effective_properties(self, theta):
        """
        Returns the effective mechanical parameters for a given streptavidin coverage theta in [0, 1].
        """
        kappa_Ld = self.kappa_Ld_0 + self.kappa_prot_sat * theta
        c0_Ld = self.c0_Ld_0 + self.alpha_c0 * theta
        # Line tension in reduced units (k_B T / um)
        # 1.5 pN = 1.5e-12 N = 1.5e-18 J/um => ~ 364.6 k_B T / um
        gamma_kBT_um = (self.gamma_0 / K_B_T) * 1e-6
        
        # Relaxed area difference shift from inner leaflet protein binding:
        # Delta A_0_bare = 8 * pi * d_m * R_ves
        d_m_um = self.d_Ld * 1e6  # 0.004 um
        delta_A0_bare = 8.0 * np.pi * d_m_um * self.R_ves  # ~ 0.50265 um^2
        
        # Inner leaflet expansion reduces Delta A_0 (outer - inner):
        delta_A0_shift = - theta * self.eta_inner_expansion * delta_A0_bare
        
        return {
            'kappa_Lo': self.kappa_Lo,
            'kappa_Ld': kappa_Ld,
            'c0_Lo': self.c0_Lo,
            'c0_Ld': c0_Ld,
            'gamma': gamma_kBT_um,
            'delta_A0_bare': delta_A0_bare,
            'delta_A0_shift': delta_A0_shift
        }


def reconstruct_profile(s_grid, psi_grid, s_split):
    """
    Reconstructs r(s) and z(s) coordinates from the tangent angle profile psi(s).
    Boundary conditions: r(0) = 0, z(0) = 0.
    """
    ds = np.diff(s_grid)
    cos_psi = np.cos(psi_grid)
    sin_psi = np.sin(psi_grid)
    
    r = np.zeros_like(s_grid)
    z = np.zeros_like(s_grid)
    
    r[1:] = np.cumsum(0.5 * (cos_psi[:-1] + cos_psi[1:]) * ds)
    z[1:] = np.cumsum(0.5 * (sin_psi[:-1] + sin_psi[1:]) * ds)
    
    idx_b = np.searchsorted(s_grid, s_split)
    return r, z, idx_b


def compute_curvatures(s_grid, psi_grid, r_grid):
    """
    Computes principal curvatures c1 = dpsi/ds, c2 = sin(psi)/r,
    mean curvature H = 0.5 * (c1 + c2), and Gaussian curvature K = c1 * c2.
    """
    dpsi_ds = np.gradient(psi_grid, s_grid)
    c2 = np.zeros_like(psi_grid)
    small_r_mask = (r_grid < 1e-5)
    
    c2[~small_r_mask] = np.sin(psi_grid[~small_r_mask]) / r_grid[~small_r_mask]
    c2[small_r_mask] = dpsi_ds[small_r_mask]
    
    c1 = dpsi_ds
    H = 0.5 * (c1 + c2)
    K = c1 * c2
    return c1, c2, H, K


def compute_energy_and_constraints(s_grid, psi_grid, s_split, params, theta):
    """
    Evaluates the complete Jülicher-Lipowsky free energy functional,
    line energy, ADE term, and geometric constraints.
    """
    props = params.get_effective_properties(theta)
    r, z, idx_b = reconstruct_profile(s_grid, psi_grid, s_split)
    c1, c2, H, K = compute_curvatures(s_grid, psi_grid, r)
    
    # Surface area element dA = 2 * pi * r * ds
    ds = np.gradient(s_grid)
    dA = 2.0 * np.pi * r * ds
    
    # Volume element dV = pi * r^2 * sin(psi) * ds
    dV = np.pi * (r**2) * np.sin(psi_grid) * ds
    
    # Split between Lo domain [0, idx_b] and Ld domain [idx_b, end]
    mask_Lo = np.arange(len(s_grid)) <= idx_b
    mask_Ld = ~mask_Lo
    
    A_Lo = np.sum(dA[mask_Lo])
    A_Ld = np.sum(dA[mask_Ld])
    A_tot = A_Lo + A_Ld
    V_tot = np.sum(dV)
    
    # Bending energy densities: 0.5 * kappa * (2H - c0)^2
    f_bend_Lo = 0.5 * props['kappa_Lo'] * ((2.0 * H[mask_Lo] - props['c0_Lo'])**2)
    f_bend_Ld = 0.5 * props['kappa_Ld'] * ((2.0 * H[mask_Ld] - props['c0_Ld'])**2)
    
    E_bend_Lo = np.sum(f_bend_Lo * dA[mask_Lo])
    E_bend_Ld = np.sum(f_bend_Ld * dA[mask_Ld])
    E_bend_tot = E_bend_Lo + E_bend_Ld
    
    # Line energy at domain boundary: 2 * pi * r_b * gamma
    r_b = r[idx_b]
    E_line = 2.0 * np.pi * r_b * props['gamma']
    
    # ADE energy:
    d_m_um = params.d_Ld * 1e6  # 0.004 um
    delta_A_geom = 2.0 * d_m_um * np.sum(H * dA)
    delta_A0_eff = props['delta_A0_bare'] + props['delta_A0_shift']
    
    k_ade = params.k_ade_ratio * props['kappa_Ld']
    E_ade = (np.pi * k_ade / (2.0 * A_tot * (d_m_um**2))) * ((delta_A_geom - delta_A0_eff)**2)
    
    # Total Free Energy
    total_physical_energy = E_bend_tot + E_line + E_ade
    
    # Constraints residuals
    res_A_Lo = (A_Lo - params.A_Lo_target) / params.A_Lo_target
    res_A_tot = (A_tot - params.A_tot_target) / params.A_tot_target
    res_V = (V_tot - params.V_tot_target) / params.V_tot_target
    res_pole = psi_grid[-1] - np.pi
    res_r_pole = r[-1]
    
    return {
        'total_energy': total_physical_energy,
        'E_bend_Lo': E_bend_Lo,
        'E_bend_Ld': E_bend_Ld,
        'E_line': E_line,
        'E_ade': E_ade,
        'r_b': r_b,
        'z_b': z[idx_b],
        'dimple_depth': z[idx_b] - z[0],
        'A_Lo': A_Lo,
        'A_Ld': A_Ld,
        'A_tot': A_tot,
        'V_tot': V_tot,
        'res_A_Lo': res_A_Lo,
        'res_A_tot': res_A_tot,
        'res_V': res_V,
        'res_pole': res_pole,
        'res_r_pole': res_r_pole,
        'r_profile': r,
        'z_profile': z,
        'H_profile': H,
        'props': props
    }
