"""
Membrane Mechanics Module for Two Opposite Lo Domains (Bipolar Vesicle - HermesLoouter).
Case: Protein located on the OUTER LEAFLET of BOTH Lo domain caps (North & South).
      Equatorial Ld belt is bare (no protein).

Defines energy functionals for a three-zone symmetric vesicle:
- Zone 1 (North Lo Cap): s in [0, s_b1] -> outer-leaflet protein
- Zone 2 (Equatorial Ld Belt): s in [s_b1, s_b2] -> bare
- Zone 3 (South Lo Cap): s in [s_b2, S_tot] -> outer-leaflet protein
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
        self.A_Ld_target = self.A_tot_target - 2.0 * self.A_Lo_target     # Two caps
        
        # Arc length targets
        self.S_tot = np.pi * self.R_ves                                   # ~15.708 um
        self.phi_Lo_sphere = 2.0 * np.arcsin(self.R_Lo / (2.0 * self.R_ves))
        self.s_split1 = self.R_ves * self.phi_Lo_sphere                  # North boundary position
        self.s_split2 = self.S_tot - self.s_split1                       # South boundary position

class TwoDomainEnergyCalculator:
    def __init__(self, params: TwoDomainParameters = None):
        self.p = params if params is not None else TwoDomainParameters()

    def compute_shape_metrics(self, s_grid: np.ndarray, psi_grid: np.ndarray, idx_b1: int, idx_b2: int):
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
        
        dpsi_ds = np.gradient(psi_grid, s_grid)
        c_m = dpsi_ds
        
        c_p = np.zeros(N)
        c_p[1:-1] = sin_psi[1:-1] / r_pos[1:-1]
        c_p[0] = dpsi_ds[0]    # North pole limit
        c_p[-1] = dpsi_ds[-1]  # South pole limit
        
        H = 0.5 * (c_m + c_p)
        K = c_m * c_p
        
        dA = 2.0 * np.pi * r_pos * ds
        A_total = float(np.sum(dA))
        A_Lo_north = float(np.sum(dA[:idx_b1 + 1]))
        A_Ld = float(np.sum(dA[idx_b1 + 1:idx_b2]))
        A_Lo_south = float(np.sum(dA[idx_b2:]))
        
        dV = np.pi * (r_pos ** 2) * sin_psi * ds
        V_total = float(np.sum(dV))
        
        r_b1 = float(r_pos[idx_b1])
        z_b1 = float(z[idx_b1])
        psi_b1 = float(psi_grid[idx_b1])
        
        r_b2 = float(r_pos[idx_b2])
        z_b2 = float(z[idx_b2])
        psi_b2 = float(psi_grid[idx_b2])
        
        # North and South cap elevations / dimple depths
        dimple_depth_north = float(z[0] - z_b1)
        dimple_depth_south = float(z_b2 - z[-1])
        elevation_north = float(z_b1 - z[0])
        elevation_south = float(z[-1] - z_b2)
        
        return {
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
            "dimple_depth_north": dimple_depth_north,
            "dimple_depth_south": dimple_depth_south,
            "elevation_north": elevation_north,
            "elevation_south": elevation_south
        }

    def get_effective_properties(self, theta: float):
        theta = float(np.clip(theta, 0.0, 1.0))
        kBT_to_J = K_B_T
        
        kappa_Lo = self.p.kappa_Lo_0 + theta * self.p.kappa_prot_sat
        kappa_Ld = self.p.kappa_Ld
        
        c0_Lo = self.p.c0_Lo_0 + self.p.alpha_c0 * theta
        c0_Ld = self.p.c0_Ld  # 0.0
        
        delta_kappa_G = self.p.gaussian_ratio * (kappa_Lo - kappa_Ld)
        
        d_m_nm = (self.p.d_Lo + self.p.d_Ld) / 2.0 * 1e9
        d_m_um = d_m_nm * 1e-3
        delta_A0_bare = 8.0 * np.pi * self.p.R_ves * d_m_um
        delta_A0_shift = + self.p.eta_outer_expansion * theta * delta_A0_bare  # Positive for outer leaflet
        
        return {
            "kappa_Lo": kappa_Lo,
            "kappa_Ld": kappa_Ld,
            "c0_Lo": c0_Lo,
            "c0_Ld": c0_Ld,
            "delta_kappa_G": delta_kappa_G,
            "d_m_um": d_m_um,
            "delta_A0_bare": delta_A0_bare,
            "delta_A0_shift": delta_A0_shift
        }

    def compute_energy(self, s_grid: np.ndarray, psi_grid: np.ndarray, idx_b1: int, idx_b2: int, theta: float):
        metrics = self.compute_shape_metrics(s_grid, psi_grid, idx_b1, idx_b2)
        eff = self.get_effective_properties(theta)
        kBT_to_J = K_B_T
        
        H = metrics["H"]
        dA = metrics["dA"]
        A_total = metrics["A_total"]
        
        kappa_Lo = eff["kappa_Lo"]
        kappa_Ld = eff["kappa_Ld"]
        c0_Lo = eff["c0_Lo"]
        c0_Ld = eff["c0_Ld"]
        
        # 1. Helfrich Bending Energy
        fbend_Lo1 = 0.5 * kappa_Lo * ((2.0 * H[:idx_b1 + 1] - c0_Lo) ** 2)
        fbend_Ld  = 0.5 * kappa_Ld * ((2.0 * H[idx_b1 + 1:idx_b2] - c0_Ld) ** 2)
        fbend_Lo2 = 0.5 * kappa_Lo * ((2.0 * H[idx_b2:] - c0_Lo) ** 2)
        
        E_bend_Lo_north = float(np.sum(fbend_Lo1 * dA[:idx_b1 + 1]))
        E_bend_Ld       = float(np.sum(fbend_Ld  * dA[idx_b1 + 1:idx_b2]))
        E_bend_Lo_south = float(np.sum(fbend_Lo2 * dA[idx_b2:]))
        E_bend_tot = E_bend_Lo_north + E_bend_Ld + E_bend_Lo_south
        
        # 2. Line Tension Energy at both boundaries: 2*pi*(r_b1 + r_b2)*gamma
        gamma_kBT_per_um = (self.p.gamma_0 / kBT_to_J) * 1e-6
        E_line = float(2.0 * np.pi * (metrics["r_b1"] + metrics["r_b2"]) * gamma_kBT_per_um)
        
        # 3. Gauss-Bonnet boundary energy jump at both interfaces
        delta_kappa_G = eff["delta_kappa_G"]
        E_gauss = float(2.0 * np.pi * delta_kappa_G * (np.cos(metrics["psi_b1"]) - np.cos(metrics["psi_b2"])))
        
        # 4. Non-local Area Difference Elasticity (ADE)
        d_m_um = eff["d_m_um"]
        k_ade = self.p.k_ade_ratio * kappa_Ld
        delta_A0_target = eff["delta_A0_bare"] + eff["delta_A0_shift"]
        delta_A_geom = 2.0 * d_m_um * float(np.sum(H * dA))
        E_ade = float((np.pi * k_ade / (2.0 * A_total * (d_m_um**2))) * ((delta_A_geom - delta_A0_target)**2))
        
        # Total free energy in k_B T
        E_total = E_bend_tot + E_line + E_gauss + E_ade
        
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
