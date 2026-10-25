"""
Membrane Mechanics Module for Two Opposite Lo Domains (Bipolar Vesicle - HermesLoouter).
Case: Protein located on the OUTER LEAFLET of BOTH Lo domain caps (North & South).
      Equatorial Ld belt is bare (no protein).

Defines energy functionals for a three-zone symmetric vesicle:
- Zone 1 (North Lo Cap): s in [0, S_Lo] -> outer-leaflet protein
- Zone 2 (Equatorial Ld Belt): s in [S_Lo, S_Lo + S_Ld] -> bare
- Zone 3 (South Lo Cap): s in [S_Lo + S_Ld, 2*S_Lo + S_Ld] -> outer-leaflet protein
"""

import numpy as np

# Physical constants in SI / biophysical units
K_B_T = 4.114e-21  # J at 298.15 K (1 k_B T in Joules)


class TwoDomainParameters:
    def __init__(
        self,
        # Lo domain properties (with outer-leaflet protein binding)
        kappa_Lo_0=80.0,          # bare Lo bending rigidity in k_B T
        kappa_prot_sat=70.0,      # protein stiffness bonus in k_B T (up to 150 k_B T)
        d_Lo=5.0e-9,              # 5.0 nm thick
        c0_Lo_0=0.0,              # bare Lo spontaneous curvature
        alpha_c0=0.20,            # outer-leaflet spontaneous curvature coupling (um^-1)

        # Ld domain properties (bare, no protein)
        kappa_Ld=20.0,            # bare Ld bending rigidity in k_B T
        d_Ld=4.0e-9,              # 4.0 nm thick
        c0_Ld=0.0,                # bare Ld spontaneous curvature

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
        self.A_Lo_target = np.pi * (self.R_Lo ** 2)                       # Single cap area
        self.A_Ld_target = self.A_tot_target - 2.0 * self.A_Lo_target     # Two caps

        # Reference spherical arc lengths
        self.S_tot = np.pi * self.R_ves                                   # ~15.708 um
        self.phi_Lo_sphere = 2.0 * np.arcsin(np.clip(self.R_Lo / (2.0 * self.R_ves), 0.0, 1.0))
        self.s_split1 = self.R_ves * self.phi_Lo_sphere                   # North boundary position
        self.s_split2 = self.S_tot - self.s_split1                        # South boundary position


class TwoDomainEnergyCalculator:
    def __init__(self, params: TwoDomainParameters = None):
        self.p = params if params is not None else TwoDomainParameters()

    def get_effective_properties(self, theta: float):
        theta = float(np.clip(theta, 0.0, 1.0))

        kappa_Lo = self.p.kappa_Lo_0 + theta * self.p.kappa_prot_sat
        kappa_Ld = self.p.kappa_Ld

        c0_Lo = self.p.c0_Lo_0 - self.p.alpha_c0 * theta
        c0_Ld = self.p.c0_Ld  # 0.0

        delta_kappa_G = self.p.gaussian_ratio * (kappa_Lo - kappa_Ld)

        f_Lo_tot = (2.0 * self.p.A_Lo_target) / self.p.A_tot_target
        f_Ld = self.p.A_Ld_target / self.p.A_tot_target

        d_m_nm = (self.p.d_Lo + self.p.d_Ld) / 2.0 * 1e9
        d_m_um = d_m_nm * 1e-3
        delta_A0_bare = 8.0 * np.pi * self.p.R_ves * d_m_um
        delta_A0_shift = - self.p.eta_outer_expansion * theta * f_Lo_tot * delta_A0_bare

        k_ade_Lo = self.p.k_ade_ratio * kappa_Lo
        k_ade_Ld = self.p.k_ade_ratio * kappa_Ld
        k_ade_eff = 1.0 / (f_Lo_tot / k_ade_Lo + f_Ld / k_ade_Ld)

        return {
            "kappa_Lo": kappa_Lo,
            "kappa_Ld": kappa_Ld,
            "c0_Lo": c0_Lo,
            "c0_Ld": c0_Ld,
            "delta_kappa_G": delta_kappa_G,
            "d_m_um": d_m_um,
            "delta_A0_bare": delta_A0_bare,
            "delta_A0_shift": delta_A0_shift,
            "k_ade_eff": k_ade_eff
        }

    def compute_energy_piecewise(
        self,
        s_Lo1: np.ndarray, psi_Lo1: np.ndarray, c1_Lo1: np.ndarray,
        s_Ld: np.ndarray, psi_Ld: np.ndarray, c1_Ld: np.ndarray,
        s_Lo2: np.ndarray, psi_Lo2: np.ndarray, c1_Lo2: np.ndarray,
        theta: float
    ):
        """
        Exact piecewise 3-zone Hamiltonian energy and geometry calculation for HermesLoouter.
        """
        eff = self.get_effective_properties(theta)
        kBT_to_J = K_B_T

        # 1. Zone 1 (North Lo)
        ds1 = np.diff(s_Lo1)
        cos1, sin1 = np.cos(psi_Lo1), np.sin(psi_Lo1)
        r1 = np.zeros_like(s_Lo1)
        r1[1:] = np.cumsum(0.5 * (cos1[:-1] + cos1[1:]) * ds1)
        z1 = np.zeros_like(s_Lo1)
        z1[1:] = np.cumsum(-0.5 * (sin1[:-1] + sin1[1:]) * ds1)

        # 2. Zone 2 (Equatorial Ld)
        ds2 = np.diff(s_Ld)
        cos2, sin2 = np.cos(psi_Ld), np.sin(psi_Ld)
        r2 = np.zeros_like(s_Ld)
        r2[0] = r1[-1]
        r2[1:] = r1[-1] + np.cumsum(0.5 * (cos2[:-1] + cos2[1:]) * ds2)
        z2 = np.zeros_like(s_Ld)
        z2[0] = z1[-1]
        z2[1:] = z1[-1] + np.cumsum(-0.5 * (sin2[:-1] + sin2[1:]) * ds2)

        # 3. Zone 3 (South Lo)
        ds3 = np.diff(s_Lo2)
        cos3, sin3 = np.cos(psi_Lo2), np.sin(psi_Lo2)
        r3 = np.zeros_like(s_Lo2)
        r3[0] = r2[-1]
        r3[1:] = r2[-1] + np.cumsum(0.5 * (cos3[:-1] + cos3[1:]) * ds3)
        z3 = np.zeros_like(s_Lo2)
        z3[0] = z2[-1]
        z3[1:] = z2[-1] + np.cumsum(-0.5 * (sin3[:-1] + sin3[1:]) * ds3)

        r1_pos = np.maximum(r1, 1e-8)
        r2_pos = np.maximum(r2, 1e-8)
        r3_pos = np.maximum(r3, 1e-8)

        # Curvatures
        c2_1 = np.sin(psi_Lo1) / r1_pos
        c2_1[0] = c1_Lo1[0]
        H1 = 0.5 * (c1_Lo1 + c2_1)
        K1 = c1_Lo1 * c2_1

        c2_2 = np.sin(psi_Ld) / r2_pos
        H2 = 0.5 * (c1_Ld + c2_2)
        K2 = c1_Ld * c2_2

        c2_3 = np.sin(psi_Lo2) / r3_pos
        c2_3[-1] = c1_Lo2[-1]
        H3 = 0.5 * (c1_Lo2 + c2_3)
        K3 = c1_Lo2 * c2_3

        # Piecewise trapezoidal areas and volumes
        A_Lo_north = float(np.trapezoid(2.0 * np.pi * r1_pos, s_Lo1))
        A_Ld = float(np.trapezoid(2.0 * np.pi * r2_pos, s_Ld))
        A_Lo_south = float(np.trapezoid(2.0 * np.pi * r3_pos, s_Lo2))
        A_total = A_Lo_north + A_Ld + A_Lo_south

        V_total = float(
            np.trapezoid(np.pi * (r1_pos ** 2) * sin1, s_Lo1) +
            np.trapezoid(np.pi * (r2_pos ** 2) * sin2, s_Ld) +
            np.trapezoid(np.pi * (r3_pos ** 2) * sin3, s_Lo2)
        )

        r_b1, z_b1, psi_b1 = float(r1_pos[-1]), float(z1[-1]), float(psi_Lo1[-1])
        r_b2, z_b2, psi_b2 = float(r3_pos[0]), float(z3[0]), float(psi_Lo2[0])

        dimple_depth_north = float(z1[0] - z_b1)
        dimple_depth_south = float(z_b2 - z3[-1])
        elevation_north = float(z_b1 - z1[0])
        elevation_south = float(z3[-1] - z_b2)

        z_north_hem = np.concatenate([z1, z2[1:len(z2) // 2 + 1]])
        z_peak_Lo_north = float(np.max(z1))
        z_peak_north = float(np.max(z_north_hem))
        peak_to_dimple_Lo_north = float(z1[0] - z_peak_Lo_north)
        peak_to_dimple_crest_north = float(z1[0] - z_peak_north)
        has_dimple_crest = (z_peak_north > z1[0] + 1e-3) and (z_peak_north > z_b1 + 1e-3)
        dome_peak_to_dimple_north = peak_to_dimple_crest_north if has_dimple_crest else dimple_depth_north
        dome_peak_to_dimple_south = dome_peak_to_dimple_north

        # 1. Helfrich Bending Energy
        kappa_Lo, kappa_Ld = eff["kappa_Lo"], eff["kappa_Ld"]
        c0_Lo, c0_Ld = eff["c0_Lo"], eff["c0_Ld"]

        fbend_1 = 0.5 * kappa_Lo * ((2.0 * H1 - c0_Lo) ** 2) * (2.0 * np.pi * r1_pos)
        fbend_2 = 0.5 * kappa_Ld * ((2.0 * H2 - c0_Ld) ** 2) * (2.0 * np.pi * r2_pos)
        fbend_3 = 0.5 * kappa_Lo * ((2.0 * H3 - c0_Lo) ** 2) * (2.0 * np.pi * r3_pos)

        E_bend_Lo_north = float(np.trapezoid(fbend_1, s_Lo1))
        E_bend_Ld = float(np.trapezoid(fbend_2, s_Ld))
        E_bend_Lo_south = float(np.trapezoid(fbend_3, s_Lo2))
        E_bend_tot = E_bend_Lo_north + E_bend_Ld + E_bend_Lo_south

        # 2. Line Tension Energy
        gamma_kBT_per_um = (self.p.gamma_0 / kBT_to_J) * 1e-6
        E_line = float(2.0 * np.pi * (r_b1 + r_b2) * gamma_kBT_per_um)

        # 3. Gauss-Bonnet boundary jump
        delta_kappa_G = eff["delta_kappa_G"]
        E_gauss = float(2.0 * np.pi * delta_kappa_G * (np.cos(psi_b1) - np.cos(psi_b2)))

        # 4. Non-local ADE Energy
        d_m_um = eff["d_m_um"]
        k_ade = eff["k_ade_eff"]
        delta_A0_target = eff["delta_A0_bare"] + eff["delta_A0_shift"]
        int_H_dA = float(
            np.trapezoid(H1 * 2.0 * np.pi * r1_pos, s_Lo1) +
            np.trapezoid(H2 * 2.0 * np.pi * r2_pos, s_Ld) +
            np.trapezoid(H3 * 2.0 * np.pi * r3_pos, s_Lo2)
        )
        delta_A_geom = 2.0 * d_m_um * int_H_dA
        E_ade = float((np.pi * k_ade / (2.0 * A_total * (d_m_um ** 2))) * ((delta_A_geom - delta_A0_target) ** 2))

        # Soft barrier against interior axis crossing r(s) < 0
        r_interior = np.concatenate([r1[1:], r2[1:], r3[1:-1]])
        min_r_interior = float(np.min(r_interior))
        neg_r_violation = np.minimum(0.0, r_interior)
        E_axis_barrier = float(1.0e6 * np.sum(neg_r_violation ** 2))

        E_total = E_bend_tot + E_line + E_gauss + E_ade + E_axis_barrier

        # One-sided Jülicher-Lipowsky boundary moment jump residual at North interface
        M_Lo1 = kappa_Lo * (2.0 * H1[-1] - c0_Lo)
        M_Ld1 = kappa_Ld * (2.0 * H2[0] - c0_Ld)
        expected_moment_jump = delta_kappa_G * (np.cos(psi_b1) / (r_b1 + 1e-8))
        moment_residual = float((M_Lo1 - M_Ld1) - expected_moment_jump)

        # Combined arrays
        s_grid = np.concatenate([s_Lo1, s_Ld[1:], s_Lo2[1:]])
        psi_grid = np.concatenate([psi_Lo1, psi_Ld[1:], psi_Lo2[1:]])
        r = np.concatenate([r1, r2[1:], r3[1:]])
        z = np.concatenate([z1, z2[1:], z3[1:]])
        H = np.concatenate([H1, H2[1:], H3[1:]])
        K = np.concatenate([K1, K2[1:], K3[1:]])
        c_m = np.concatenate([c1_Lo1, c1_Ld[1:], c1_Lo2[1:]])
        c_p = np.concatenate([c2_1, c2_2[1:], c2_3[1:]])

        r_pos_all = np.maximum(r, 1e-8)
        ds_all = np.gradient(s_grid)
        dA = 2.0 * np.pi * r_pos_all * ds_all
        dV = np.pi * (r_pos_all ** 2) * np.sin(psi_grid) * ds_all

        metrics = {
            "s": s_grid,
            "psi": psi_grid,
            "r": r,
            "z": z,
            "H": H,
            "K": K,
            "c_m": c_m,
            "c_p": c_p,
            "dA": dA,
            "dV": dV,
            "A_total": A_total,
            "A_Lo_north": A_Lo_north,
            "A_Ld": A_Ld,
            "A_Lo_south": A_Lo_south,
            "V_total": V_total,
            "r_b1": r_b1,
            "z_b1": z_b1,
            "psi_b1": psi_b1,
            "r_b2": r_b2,
            "z_b2": z_b2,
            "psi_b2": psi_b2,
            "s_split1": float(s_Lo1[-1]),
            "s_split2": float(s_Lo2[0]),
            "S_tot": float(s_Lo2[-1]),
            "min_r_interior": min_r_interior,
            "moment_residual": moment_residual,
            "dimple_depth_north": dimple_depth_north,
            "dimple_depth_south": dimple_depth_south,
            "elevation_north": elevation_north,
            "elevation_south": elevation_south,
            "peak_to_dimple_Lo_north": peak_to_dimple_Lo_north,
            "peak_to_dimple_crest_north": peak_to_dimple_crest_north,
            "dome_peak_to_dimple_north": dome_peak_to_dimple_north,
            "dome_peak_to_dimple_south": dome_peak_to_dimple_south
        }

        return {
            "E_total": E_total,
            "E_bend_tot": E_bend_tot,
            "E_bend_Lo_north": E_bend_Lo_north,
            "E_bend_Ld": E_bend_Ld,
            "E_bend_Lo_south": E_bend_Lo_south,
            "E_line": E_line,
            "E_gauss": E_gauss,
            "E_ade": E_ade,
            "metrics": metrics
        }

    def compute_shape_metrics(self, s_grid: np.ndarray, psi_grid: np.ndarray, idx_b1: int, idx_b2: int):
        res = self.compute_energy(s_grid, psi_grid, idx_b1, idx_b2, theta=0.0)
        return res["metrics"]

    def compute_energy(self, s_grid: np.ndarray, psi_grid: np.ndarray, idx_b1: int, idx_b2: int, theta: float):
        s_Lo1 = s_grid[:idx_b1 + 1]
        psi_Lo1 = psi_grid[:idx_b1 + 1]
        c1_Lo1 = np.gradient(psi_Lo1, s_Lo1)

        s_Ld = s_grid[idx_b1:idx_b2 + 1]
        psi_Ld = psi_grid[idx_b1:idx_b2 + 1]
        c1_Ld = np.gradient(psi_Ld, s_Ld)

        s_Lo2 = s_grid[idx_b2:]
        psi_Lo2 = psi_grid[idx_b2:]
        c1_Lo2 = np.gradient(psi_Lo2, s_Lo2)

        return self.compute_energy_piecewise(
            s_Lo1, psi_Lo1, c1_Lo1,
            s_Ld, psi_Ld, c1_Ld,
            s_Lo2, psi_Lo2, c1_Lo2,
            theta
        )
