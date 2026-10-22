"""
Biophysical Mechanics Module for HermesLoouter (PRL Standard).
Case: Protein located on the OUTER LEAFLET of the Liquid-Ordered (Lo) domain.
      Liquid-Disordered (Ld) domain is bare (no protein).

Physical characteristics:
1. Lo domain (with outer-leaflet protein):
   - kappa_Lo = kappa_Lo_0 + theta * kappa_prot_sat (80 to 150 k_B T)
   - c0_Lo = + alpha_c0 * theta (positive spontaneous curvature towards exterior, +0.20 * theta um^-1)
   - Outer-leaflet area expansion weighted by Lo area fraction f_Lo = A_Lo_target / A_tot_target
2. Ld domain (bare):
   - kappa_Ld = 20 k_B T
   - c0_Ld = 0.0 um^-1
3. Jülicher-Lipowsky line tension gamma at domain boundary.
4. Gauss-Bonnet boundary modulus difference jump with one-sided curvature jump diagnostics.
5. Area-Difference Elasticity (ADE) with area-weighted harmonic-mean modulus k_ade_eff.
6. Exact South-pole closure residual r(S_tot) / R_ves.
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
        k_ade_ratio=2.0,          # k_ade / kappa ratio
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

        # Reference spherical arc lengths
        self.S_tot = np.pi * self.R_ves                                   # ~15.708 um
        self.phi_Lo_sphere = 2.0 * np.arcsin(np.clip(self.R_Lo / (2.0 * self.R_ves), 0.0, 1.0))
        self.s_split = self.R_ves * self.phi_Lo_sphere                    # Reference arc-length boundary


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

    # 5. ADE Parameters weighted by protein-bound Lo domain area fraction
    f_Lo = params.A_Lo_target / params.A_tot_target
    f_Ld = params.A_Ld_target / params.A_tot_target
    d_m_nm = (params.d_Lo + params.d_Ld) / 2.0 * 1e9
    d_m_um = d_m_nm * 1e-3
    delta_A0_bare = 8.0 * np.pi * params.R_ves * d_m_um
    delta_A0_shift = + params.eta_outer_expansion * theta * f_Lo * delta_A0_bare

    k_ade_Lo = params.k_ade_ratio * kappa_Lo
    k_ade_Ld = params.k_ade_ratio * kappa_Ld
    k_ade_eff = 1.0 / (f_Lo / k_ade_Lo + f_Ld / k_ade_Ld)

    return {
        'kappa_Lo': kappa_Lo,
        'kappa_Ld': kappa_Ld,
        'c0_Lo': c0_Lo,
        'c0_Ld': c0_Ld,
        'gamma': gamma_kBT_per_um,
        'delta_kappa_G': delta_kappa_G,
        'd_m_um': d_m_um,
        'delta_A0_bare': delta_A0_bare,
        'delta_A0_shift': delta_A0_shift,
        'k_ade_eff': k_ade_eff
    }


def reconstruct_profile(s_grid, psi_grid):
    """
    Reconstructs axisymmetric coordinates r(s), z(s), area elements dA(s),
    and cumulative volume from tangent angle profile psi(s) using exact trapezoidal quadrature.
    """
    N = len(s_grid)
    ds_seg = np.diff(s_grid)
    cos_psi = np.cos(psi_grid)
    sin_psi = np.sin(psi_grid)

    r = np.zeros(N)
    r[1:] = np.cumsum(0.5 * (cos_psi[:-1] + cos_psi[1:]) * ds_seg)
    r[0] = 0.0

    z = np.zeros(N)
    z[1:] = np.cumsum(-0.5 * (sin_psi[:-1] + sin_psi[1:]) * ds_seg)
    z[0] = 0.0

    r_pos = np.maximum(r, 1e-8)
    ds_node = np.gradient(s_grid)
    dA = 2.0 * np.pi * r_pos * ds_node
    dV = np.pi * (r_pos ** 2) * sin_psi * ds_node

    A_tot = float(np.trapezoid(2.0 * np.pi * r_pos, s_grid))
    V_tot = float(np.trapezoid(np.pi * (r_pos ** 2) * sin_psi, s_grid))

    return r, z, dA, dV, A_tot, V_tot


def compute_curvatures(s_grid, psi_grid, r_pos, c1_override=None):
    """
    Computes principal curvatures with L'Hopital regularization at both poles.
    """
    if c1_override is not None:
        c1 = np.asarray(c1_override, dtype=float)
    else:
        c1 = np.gradient(psi_grid, s_grid)

    r_safe = np.maximum(r_pos, 1e-8)
    c2 = np.sin(psi_grid) / r_safe
    c2[0] = c1[0]      # L'Hopital limit at North pole
    c2[-1] = c1[-1]    # L'Hopital limit at South pole

    H = 0.5 * (c1 + c2)
    K = c1 * c2
    return c1, c2, H, K


def compute_energy_breakdown_piecewise(s_Lo, psi_Lo, c1_Lo, s_Ld, psi_Ld, c1_Ld, params, theta):
    """
    Exact piecewise-domain Hamiltonian energy and constraint evaluation for HermesLoouter
    using separate C^2 splines on Lo [0, S_Lo] and Ld [S_Lo, S_tot] joined C^0 at S_Lo.
    """
    props = get_effective_properties(params, theta)

    # Coordinate integration on Lo [0, S_Lo]
    ds_Lo = np.diff(s_Lo)
    cos_Lo, sin_Lo = np.cos(psi_Lo), np.sin(psi_Lo)
    r_Lo = np.zeros_like(s_Lo)
    r_Lo[1:] = np.cumsum(0.5 * (cos_Lo[:-1] + cos_Lo[1:]) * ds_Lo)
    z_Lo = np.zeros_like(s_Lo)
    z_Lo[1:] = np.cumsum(-0.5 * (sin_Lo[:-1] + sin_Lo[1:]) * ds_Lo)

    # Coordinate integration on Ld [S_Lo, S_tot]
    ds_Ld = np.diff(s_Ld)
    cos_Ld, sin_Ld = np.cos(psi_Ld), np.sin(psi_Ld)
    r_Ld = np.zeros_like(s_Ld)
    r_Ld[0] = r_Lo[-1]
    r_Ld[1:] = r_Lo[-1] + np.cumsum(0.5 * (cos_Ld[:-1] + cos_Ld[1:]) * ds_Ld)
    z_Ld = np.zeros_like(s_Ld)
    z_Ld[0] = z_Lo[-1]
    z_Ld[1:] = z_Lo[-1] + np.cumsum(-0.5 * (sin_Ld[:-1] + sin_Ld[1:]) * ds_Ld)

    r_Lo_pos = np.maximum(r_Lo, 1e-8)
    r_Ld_pos = np.maximum(r_Ld, 1e-8)

    # Curvatures on Lo (North pole limit)
    c2_Lo = np.sin(psi_Lo) / r_Lo_pos
    c2_Lo[0] = c1_Lo[0]
    H_Lo = 0.5 * (c1_Lo + c2_Lo)
    K_Lo = c1_Lo * c2_Lo

    # Curvatures on Ld (South pole limit)
    c2_Ld = np.sin(psi_Ld) / r_Ld_pos
    c2_Ld[-1] = c1_Ld[-1]
    H_Ld = 0.5 * (c1_Ld + c2_Ld)
    K_Ld = c1_Ld * c2_Ld

    # Exact piecewise trapezoidal areas and volumes
    A_Lo = float(np.trapezoid(2.0 * np.pi * r_Lo_pos, s_Lo))
    A_Ld = float(np.trapezoid(2.0 * np.pi * r_Ld_pos, s_Ld))
    A_tot = A_Lo + A_Ld

    V_Lo = float(np.trapezoid(np.pi * (r_Lo_pos ** 2) * sin_Lo, s_Lo))
    V_Ld = float(np.trapezoid(np.pi * (r_Ld_pos ** 2) * sin_Ld, s_Ld))
    V_tot = V_Lo + V_Ld

    # 1. Helfrich Bending Energy (piecewise trapezoidal)
    f_bend_Lo = 0.5 * props['kappa_Lo'] * ((2.0 * H_Lo - props['c0_Lo']) ** 2) * (2.0 * np.pi * r_Lo_pos)
    f_bend_Ld = 0.5 * props['kappa_Ld'] * ((2.0 * H_Ld - props['c0_Ld']) ** 2) * (2.0 * np.pi * r_Ld_pos)
    E_bend_Lo = float(np.trapezoid(f_bend_Lo, s_Lo))
    E_bend_Ld = float(np.trapezoid(f_bend_Ld, s_Ld))
    E_bend_tot = E_bend_Lo + E_bend_Ld

    # 2. Line Tension Energy at domain boundary
    r_b = float(r_Lo_pos[-1])
    z_b = float(z_Lo[-1])
    psi_b = float(psi_Lo[-1])
    E_line = float(2.0 * np.pi * r_b * props['gamma'])

    # 3. Gauss-Bonnet Boundary Jump Energy
    E_gauss = float(2.0 * np.pi * props['delta_kappa_G'] * np.cos(psi_b))

    # 4. Non-local Area-Difference Elasticity (ADE)
    d_m_um = props['d_m_um']
    int_H_dA = float(
        np.trapezoid(H_Lo * 2.0 * np.pi * r_Lo_pos, s_Lo) +
        np.trapezoid(H_Ld * 2.0 * np.pi * r_Ld_pos, s_Ld)
    )
    delta_A_geom = 2.0 * d_m_um * int_H_dA
    delta_A0_eff = props['delta_A0_bare'] + props['delta_A0_shift']
    k_ade = props['k_ade_eff']
    E_ade = float((np.pi * k_ade / (2.0 * A_tot * (d_m_um ** 2))) * ((delta_A_geom - delta_A0_eff) ** 2))

    # Soft barrier against interior axis crossing r(s) < 0
    r_interior = np.concatenate([r_Lo[1:], r_Ld[:-1]])
    min_r_interior = float(np.min(r_interior))
    neg_r_violation = np.minimum(0.0, r_interior)
    E_axis_barrier = float(1.0e6 * np.sum(neg_r_violation ** 2))

    total_physical_energy = E_bend_tot + E_line + E_gauss + E_ade + E_axis_barrier

    # Exact geometric constraint residuals
    res_V = (V_tot - params.V_tot_target) / params.V_tot_target
    res_A_Lo = (A_Lo - params.A_Lo_target) / params.A_Lo_target
    res_A_tot = (A_tot - params.A_tot_target) / params.A_tot_target
    res_r_pole = float(r_Ld[-1]) / params.R_ves
    res_psi_pole = float(psi_Ld[-1] - np.pi)

    # One-sided Jülicher-Lipowsky boundary moment jump residual
    M_Lo = props['kappa_Lo'] * (2.0 * H_Lo[-1] - props['c0_Lo'])
    M_Ld = props['kappa_Ld'] * (2.0 * H_Ld[0] - props['c0_Ld'])
    expected_moment_jump = props['delta_kappa_G'] * (np.cos(psi_b) / (r_b + 1e-8))
    moment_residual = float((M_Lo - M_Ld) - expected_moment_jump)

    # Combined meridian arrays
    s_grid = np.concatenate([s_Lo, s_Ld[1:]])
    psi_grid = np.concatenate([psi_Lo, psi_Ld[1:]])
    r = np.concatenate([r_Lo, r_Ld[1:]])
    z = np.concatenate([z_Lo, z_Ld[1:]])
    H = np.concatenate([H_Lo, H_Ld[1:]])
    K = np.concatenate([K_Lo, K_Ld[1:]])
    c1 = np.concatenate([c1_Lo, c1_Ld[1:]])
    c2 = np.concatenate([c2_Lo, c2_Ld[1:]])

    dimple_depth = float(z_Lo[0] - z_b)
    cap_elevation = float(z_b - z_Lo[0])

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
        's_split': float(s_Lo[-1]),
        'S_tot': float(s_Ld[-1]),
        'dimple_depth': dimple_depth,
        'cap_elevation': cap_elevation,
        'res_V': res_V,
        'res_A_Lo': res_A_Lo,
        'res_A_tot': res_A_tot,
        'res_r_pole': res_r_pole,
        'res_psi_pole': res_psi_pole,
        'min_r_interior': min_r_interior,
        'moment_residual': moment_residual,
        's': s_grid,
        'psi': psi_grid,
        'r': r,
        'z': z,
        'H': H,
        'K': K,
        'c1': c1,
        'c2': c2
    }


def compute_energy_breakdown(s_grid, psi_grid, s_split, params, theta):
    """
    Backward-compatible energy breakdown entry point for a single combined s_grid.
    """
    mask_Lo = (s_grid <= s_split + 1e-12)
    idx_b = int(np.where(mask_Lo)[0][-1])

    s_Lo = s_grid[:idx_b + 1]
    psi_Lo = psi_grid[:idx_b + 1]
    c1_Lo = np.gradient(psi_Lo, s_Lo)

    s_Ld = s_grid[idx_b:]
    psi_Ld = psi_grid[idx_b:]
    c1_Ld = np.gradient(psi_Ld, s_Ld)

    return compute_energy_breakdown_piecewise(s_Lo, psi_Lo, c1_Lo, s_Ld, psi_Ld, c1_Ld, params, theta)
