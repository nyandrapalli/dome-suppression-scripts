"""
Biophysical Mechanics Module for HermesLoouter (PRL Standard).
Case: Protein located on the OUTER LEAFLET of the Liquid-Ordered (Lo) domain.
      Liquid-Disordered (Ld) domain is bare (no protein).

Physical characteristics:
1. Lo domain (with outer-leaflet protein):
   - kappa_Lo = kappa_Lo_0 + theta * kappa_prot_sat (80 to 150 k_B T)
   - c0_Lo = + alpha_c0 * theta (positive spontaneous curvature towards exterior, +0.20 * theta um^-1)
   - Outer-leaflet area expansion: delta_A0_eff = delta_A0_bare * (1 + eta_outer * theta)
2. Ld domain (bare):
   - kappa_Ld = 20 k_B T
   - c0_Ld = 0.0 um^-1
3. Jülicher-Lipowsky line tension gamma at domain boundary.
4. Gauss-Bonnet boundary modulus difference jump.
5. Area-Difference Elasticity (ADE).
"""

import numpy as np

# Physical constants in SI / biophysical units
K_B_T = 4.114e-21  # J at 298.15 K (1 k_B T in Joules)

class MembraneParameters:
    def __init__(
        self,
        # Lo domain properties (with outer-leaflet protein binding)
        kappa_Lo_0=80.0,          # bare Lo bending rigidity in k_B T
        kappa_prot_sat=70.0,      # protein stiffness saturation bonus in k_B T (up to 150 k_B T)
        d_Lo=5.0e-9,              # 5.0 nm thick
        c0_Lo_0=0.0,              # bare Lo spontaneous curvature
        alpha_c0=0.20,            # outer-leaflet spontaneous curvature coupling (um^-1)
        
        # Ld domain properties (bare, no protein)
        kappa_Ld=20.0,            # bare Ld bending rigidity in k_B T
        d_Ld=4.0e-9,              # 4.0 nm thick
        c0_Ld=0.0,                # bare Ld spontaneous curvature (no protein)
        
        # Line tension and Gauss-Bonnet
        gamma_0=1.5e-12,          # Line tension in N (1.5 pN)
        gaussian_ratio=0.9,       # Delta kappa_G / (kappa_Lo - kappa_Ld)
        
        # Area-Difference Elasticity (ADE)
        k_ade_ratio=2.0,          # k_ade / kappa_Ld ratio
        eta_outer_expansion=0.08, # +8% outer leaflet expansion at saturation
        
        # Geometric scale
        R_ves=5.0,                # Vesicle radius in um
        D_Lo=3.0,                 # Lo domain diameter in um (radius = 1.5 um)
        reduced_volume=0.98       # Deflation / reduced volume v = V / V_sphere
    ):
        self.kappa_Lo_0 = float(kappa_Lo_0)
        self.kappa_prot_sat = float(kappa_prot_sat)
        self.d_Lo = float(d_Lo)
        self.c0_Lo_0 = float(c0_Lo_0)
        self.alpha_c0 = float(alpha_c0)
        self.kappa_Ld = float(kappa_Ld)
        self.d_Ld = float(d_Ld)
        self.c0_Ld = float(c0_Ld)
        self.gamma_0 = float(gamma_0)
        self.gaussian_ratio = float(gaussian_ratio)
        self.k_ade_ratio = float(k_ade_ratio)
        self.eta_outer_expansion = float(eta_outer_expansion)
        self.R_ves = float(R_ves)
        self.D_Lo = float(D_Lo)
        self.R_Lo = self.D_Lo / 2.0
        self.reduced_volume = float(reduced_volume)
        
        # Geometric target references
        self.A_tot_target = 4.0 * np.pi * (self.R_ves ** 2)               # ~314.159 um^2
        self.V_sphere = (4.0 / 3.0) * np.pi * (self.R_ves ** 3)          # ~523.599 um^3
        self.V_tot_target = self.reduced_volume * self.V_sphere           # Target volume under deflation
        self.A_Lo_target = np.pi * (self.R_Lo ** 2)                       # Cap area
        self.A_Ld_target = self.A_tot_target - self.A_Lo_target
        
        # Arc length targets
        self.S_tot = np.pi * self.R_ves                                   # ~15.708 um
        self.phi_Lo_sphere = 2.0 * np.arcsin(self.R_Lo / (2.0 * self.R_ves))
        self.s_split = self.R_ves * self.phi_Lo_sphere                   # Arc-length boundary position


def get_effective_properties(params, theta):
    """
    Computes effective physical properties for outer-leaflet protein on Lo domain.
    theta: protein coverage on Lo in [0, 1].
    """
    theta = float(np.clip(theta, 0.0, 1.0))
    kBT_to_J = K_B_T
    
    # 1. Bending rigidities: protein is on Lo domain, Ld is bare
    kappa_Lo = params.kappa_Lo_0 + theta * params.kappa_prot_sat
    kappa_Ld = params.kappa_Ld
    
    # 2. Spontaneous curvatures (um^-1): outer leaflet Lo is positive, Ld is 0
    c0_Lo = params.c0_Lo_0 + params.alpha_c0 * theta
    c0_Ld = params.c0_Ld  # 0.0
    
    # 3. Line tension in k_B T / um
    gamma_SI = params.gamma_0
    gamma_kBT_per_um = (gamma_SI / kBT_to_J) * 1e-6
    
    # 4. Gauss-Bonnet bending modulus difference
    delta_kappa_G = params.gaussian_ratio * (kappa_Lo - kappa_Ld)
    
    # 5. ADE Parameters: Outer leaflet expansion increases Delta A_0 (outer - inner)
    d_m_nm = (params.d_Lo + params.d_Ld) / 2.0 * 1e9
    d_m_um = d_m_nm * 1e-3
    delta_A0_bare = 8.0 * np.pi * params.R_ves * d_m_um
    delta_A0_shift = + params.eta_outer_expansion * theta * delta_A0_bare  # Positive for outer leaflet
    
    return {
        'kappa_Lo': kappa_Lo,
        'kappa_Ld': kappa_Ld,
        'c0_Lo': c0_Lo,
        'c0_Ld': c0_Ld,
        'gamma': gamma_kBT_per_um,
        'delta_kappa_G': delta_kappa_G,
        'd_m_um': d_m_um,
        'delta_A0_bare': delta_A0_bare,
        'delta_A0_shift': delta_A0_shift
    }


def reconstruct_profile(s_grid, psi_grid):
    """
    Reconstructs axisymmetric coordinates r(s), z(s), area elements dA(s),
    and cumulative volume from tangent angle profile psi(s).
    """
    N = len(s_grid)
    ds = np.gradient(s_grid)
    cos_psi = np.cos(psi_grid)
    sin_psi = np.sin(psi_grid)
    
    dr = cos_psi * ds
    r = np.zeros(N)
    r[1:] = np.cumsum(0.5 * (dr[:-1] + dr[1:]))
    r[0] = 0.0
    r_pos = np.maximum(r, 1e-8)
    
    dz = -sin_psi * ds
    z = np.zeros(N)
    z[1:] = np.cumsum(0.5 * (dz[:-1] + dz[1:]))
    z[0] = 0.0
    
    dA = 2.0 * np.pi * r_pos * ds
    dV = np.pi * (r_pos ** 2) * sin_psi * ds
    
    A_tot = float(np.sum(dA))
    V_tot = float(np.sum(dV))
    
    return r, z, dA, dV, A_tot, V_tot


def compute_curvatures(s_grid, psi_grid, r_pos):
    """
    Computes principal curvatures with L'Hopital regularization at both poles.
    """
    dpsi_ds = np.gradient(psi_grid, s_grid)
    c1 = dpsi_ds
    
    r_safe = np.maximum(r_pos, 1e-8)
    c2 = np.sin(psi_grid) / r_safe
    c2[0] = dpsi_ds[0]      # L'Hopital limit at North pole
    c2[-1] = dpsi_ds[-1]    # L'Hopital limit at South pole
    
    H = 0.5 * (c1 + c2)
    K = c1 * c2
    return c1, c2, H, K


def compute_energy_breakdown(s_grid, psi_grid, s_split, params, theta):
    """
    Computes the physical energy breakdown for HermesLoouter:
    E_total = E_bend_Lo + E_bend_Ld + E_line + E_gauss + E_ade
    """
    props = get_effective_properties(params, theta)
    r, z, dA, dV, A_tot, V_tot = reconstruct_profile(s_grid, psi_grid)
    r_pos = np.maximum(r, 1e-8)
    c1, c2, H, K = compute_curvatures(s_grid, psi_grid, r_pos)
    
    mask_Lo = (s_grid <= s_split)
    mask_Ld = (s_grid > s_split)
    idx_b = np.where(mask_Lo)[0][-1]
    
    A_Lo = float(np.sum(dA[mask_Lo]))
    A_Ld = float(np.sum(dA[mask_Ld]))
    
    # 1. Helfrich Bending Energy
    f_bend_Lo = 0.5 * props['kappa_Lo'] * ((2.0 * H[mask_Lo] - props['c0_Lo'])**2)
    f_bend_Ld = 0.5 * props['kappa_Ld'] * ((2.0 * H[mask_Ld] - props['c0_Ld'])**2)
    
    E_bend_Lo = float(np.sum(f_bend_Lo * dA[mask_Lo]))
    E_bend_Ld = float(np.sum(f_bend_Ld * dA[mask_Ld]))
    E_bend_tot = E_bend_Lo + E_bend_Ld
    
    # 2. Line Tension Energy at domain boundary
    r_b = float(r_pos[idx_b])
    psi_b = float(psi_grid[idx_b])
    E_line = float(2.0 * np.pi * r_b * props['gamma'])
    
    # 3. Gauss-Bonnet Boundary Jump Energy
    E_gauss = float(2.0 * np.pi * props['delta_kappa_G'] * np.cos(psi_b))
    
    # 4. Non-local Area-Difference Elasticity (ADE)
    d_m_um = props['d_m_um']
    delta_A_geom = 2.0 * d_m_um * float(np.sum(H * dA))
    delta_A0_eff = props['delta_A0_bare'] + props['delta_A0_shift']
    k_ade = params.k_ade_ratio * props['kappa_Ld']
    E_ade = float((np.pi * k_ade / (2.0 * A_tot * (d_m_um**2))) * ((delta_A_geom - delta_A0_eff)**2))
    
    # Total exact physical free energy
    total_physical_energy = E_bend_tot + E_line + E_gauss + E_ade
    
    # Exact geometric constraint residuals
    res_V = (V_tot - params.V_tot_target) / params.V_tot_target
    res_A_Lo = (A_Lo - params.A_Lo_target) / params.A_Lo_target
    res_A_tot = (A_tot - params.A_tot_target) / params.A_tot_target
    res_r_pole = float(r[-1])
    res_psi_pole = float(psi_grid[-1] - np.pi)
    
    # Cap depth and elevation:
    # cap_elevation: relative axial displacement of Lo apex above boundary rim
    # dimple_depth: relative axial displacement of Lo apex below boundary rim
    dimple_depth = float(z[0] - z[idx_b])
    cap_elevation = float(z[idx_b] - z[0])
    z_b = float(z[idx_b])
    
    return {
        'total_energy': total_physical_energy,
        'E_bend_Lo': E_bend_Lo,
        'E_bend_Ld': E_bend_Ld,
        'E_bend_tot': E_bend_tot,
        'E_line': E_line,
        'E_gauss': E_gauss,
        'E_ade': E_ade,
        'A_Lo': A_Lo,
        'A_Ld': A_Ld,
        'A_tot': A_tot,
        'V_tot': V_tot,
        'r_b': r_b,
        'z_b': z_b,
        'psi_b': psi_b,
        'dimple_depth': dimple_depth,
        'cap_elevation': cap_elevation,
        'res_V': res_V,
        'res_A_Lo': res_A_Lo,
        'res_A_tot': res_A_tot,
        'res_r_pole': res_r_pole,
        'res_psi_pole': res_psi_pole,
        'r': r,
        'z': z,
        'H': H,
        'K': K,
        'c1': c1,
        'c2': c2
    }
