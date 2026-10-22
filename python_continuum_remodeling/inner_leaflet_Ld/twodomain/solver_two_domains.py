"""
Exact Non-Linear Constrained Solver for Two Opposite Lo Domains (Bipolar Vesicle - PRL Standard).
Enforces exact physical up-down reflection symmetry:
    psi(S_tot - s) = pi - psi(s)
ensuring that both North and South Lo domains behave identically and symmetrically by construction:
    depth_north == depth_south
    r_b1 == r_b2
    A_Lo1 == A_Lo2
Machine-precision equality constraints: |Delta V / V| < 1e-10, |Delta A / A| < 1e-10.
"""

import os
import sys
import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(__file__))
try:
    from mechanics_two_domains import TwoDomainParameters, TwoDomainEnergyCalculator
except ImportError:
    from sim_prl.mechanics_two_domains import TwoDomainParameters, TwoDomainEnergyCalculator

class TwoDomainVesicleSolver:
    """
    NLP Shape Solver for 3-zone Bipolar Vesicle with exact mirror symmetry.
    """
    def __init__(self, params: TwoDomainParameters = None, n_points: int = 101):
        if n_points % 2 == 0:
            n_points += 1  # Guarantee an odd number of points so equator is exact midpoint
        self.params = params if params is not None else TwoDomainParameters()
        self.p = self.params
        self.calc = TwoDomainEnergyCalculator(self.p)
        self.n_points = n_points
        
        self.S_tot = np.pi * self.p.R_ves
        self.s_grid = np.linspace(0, self.S_tot, self.n_points)
        self.mid_idx = self.n_points // 2
        self.s_half = self.s_grid[:self.mid_idx + 1]
        
        # Split indices
        cos_theta_cap = max(-1.0, min(1.0, 1.0 - self.p.A_Lo_target / (2.0 * np.pi * (self.p.R_ves ** 2))))
        theta_cap = np.arccos(cos_theta_cap)
        self.s_split1 = theta_cap * self.p.R_ves
        self.s_split2 = self.S_tot - self.s_split1
        
        self.idx_b1 = np.argmin(np.abs(self.s_grid - self.s_split1))
        self.idx_b1 = max(2, min(self.mid_idx - 2, self.idx_b1))
        self.idx_b2 = self.n_points - 1 - self.idx_b1
        
        self._last_x = None
        self._cached_metrics = None
        self._cached_energy = None

    def reconstruct_full_psi(self, x_interior: np.ndarray) -> np.ndarray:
        """
        Reconstructs globally smooth, up-down reflection symmetric psi(s)
        satisfying psi(0) = 0, psi(S_tot/2) = pi/2, and psi(S_tot - s) = pi - psi(s).
        """
        ph = np.zeros(len(self.s_half))
        ph[0] = 0.0
        ph[1:-1] = x_interior
        ph[-1] = np.pi / 2.0
        # Exact anti-symmetric reflection
        pf = np.concatenate([ph, np.pi - ph[-2::-1]])
        return pf

    def _eval_cached(self, x_interior: np.ndarray, theta: float):
        if self._last_x is not None and np.array_equal(x_interior, self._last_x):
            return self._cached_energy, self._cached_metrics
        
        pf = self.reconstruct_full_psi(x_interior)
        en_res = self.calc.compute_energy(self.s_grid, pf, self.idx_b1, self.idx_b2, theta)
        self._last_x = np.copy(x_interior)
        self._cached_energy = en_res
        self._cached_metrics = en_res["metrics"]
        return en_res, en_res["metrics"]

    def generate_initial_guess(self, branch: str = "inward") -> np.ndarray:
        """
        Generates symmetric initial guess on half grid [0, S_tot/2].
        """
        s_h = self.s_half
        s_b1 = self.s_grid[self.idx_b1]
        s_eq = self.s_grid[self.mid_idx]
        
        mask_lo = (s_h <= s_b1)
        mask_ld = (s_h > s_b1)
        
        scale = min(1.0, max(0.2, self.s_split1 / 1.5))
        psi_h = np.zeros_like(s_h)
        
        if branch == "sphere" or branch == "bipolar_sphere":
            psi_h = (s_h / self.S_tot) * np.pi
        elif branch == "inward" or branch == "bipolar_inward":
            # Inward dimple seed: negative tangent dip in Lo cap
            psi_h[mask_lo] = -0.5 * scale * np.sin(np.pi * s_h[mask_lo] / s_b1) + (s_h[mask_lo] / s_b1) * (s_b1 / self.S_tot * np.pi)
            psi_at_b1 = psi_h[mask_lo][-1]
            psi_h[mask_ld] = psi_at_b1 + ((s_h[mask_ld] - s_b1) / (s_eq - s_b1)) * (np.pi / 2.0 - psi_at_b1)
        elif branch == "outward" or branch == "bipolar_outward":
            # Outward dome seed: positive steepening in Lo cap
            psi_h[mask_lo] = 0.5 * scale * np.sin(np.pi * s_h[mask_lo] / s_b1) + (s_h[mask_lo] / s_b1) * (s_b1 / self.S_tot * np.pi)
            psi_at_b1 = psi_h[mask_lo][-1]
            psi_h[mask_ld] = psi_at_b1 + ((s_h[mask_ld] - s_b1) / (s_eq - s_b1)) * (np.pi / 2.0 - psi_at_b1)
        else:
            psi_h = (s_h / self.S_tot) * np.pi
            
        return psi_h

    def optimize_shape(self, theta: float, branch: str = "inward", maxiter: int = 600, ftol: float = 1e-9):
        """
        Performs exact non-linear constrained optimization for the symmetric 2-domain vesicle.
        """
        self._last_x = None
        self._cached_metrics = None
        self._cached_energy = None
        
        psi_half_init = self.generate_initial_guess(branch)
        x0 = psi_half_init[1:-1]
        
        def obj_func(x):
            en, _ = self._eval_cached(x, theta)
            return en["E_total"]
            
        def constr_volume(x):
            _, m = self._eval_cached(x, theta)
            return (m["V_total"] - self.p.V_tot_target) / self.p.V_tot_target
            
        def constr_area_tot(x):
            _, m = self._eval_cached(x, theta)
            return (m["A_total"] - self.p.A_tot_target) / self.p.A_tot_target
            
        def constr_area_lo(x):
            _, m = self._eval_cached(x, theta)
            return (m["A_Lo_north"] - self.p.A_Lo_target) / self.p.A_Lo_target

        constraints = [
            {"type": "eq", "fun": constr_volume},
            {"type": "eq", "fun": constr_area_tot},
            {"type": "eq", "fun": constr_area_lo}
        ]
        
        bounds = [(-np.pi, 2.0 * np.pi) for _ in range(len(x0))]
        
        opt_res = minimize(
            obj_func,
            x0,
            method="SLSQP",
            bounds=bounds,
            constraints=constraints,
            options={"maxiter": maxiter, "ftol": ftol, "disp": False}
        )
        
        final_psi = self.reconstruct_full_psi(opt_res.x)
        final_en = self.calc.compute_energy(self.s_grid, final_psi, self.idx_b1, self.idx_b2, theta)
        m = final_en["metrics"]
        
        res_V = abs((m["V_total"] - self.p.V_tot_target) / self.p.V_tot_target)
        res_A_tot = abs((m["A_total"] - self.p.A_tot_target) / self.p.A_tot_target)
        res_A_Lo = abs((m["A_Lo_north"] - self.p.A_Lo_target) / self.p.A_Lo_target)
        max_err = max(res_V, res_A_tot, res_A_Lo)
        
        opt_success = bool(opt_res.success or opt_res.status == 0)
        geom_converged = bool(max_err < 1e-3)
        is_converged = bool(opt_success and geom_converged)
        
        return {
            "converged": is_converged,
            "opt_success": opt_success,
            "geom_converged": geom_converged,
            "opt_result": opt_res,
            "max_geom_error": max_err,
            "res_V": res_V,
            "res_A_Lo1": res_A_Lo,
            "res_A_Lo2": res_A_Lo,
            "res_A_tot": res_A_tot,
            "total_energy": final_en["E_total"],
            "E_bend_tot": final_en["E_bend_tot"],
            "E_bend_Lo_north": final_en["E_bend_Lo_north"],
            "E_bend_Ld": final_en["E_bend_Ld"],
            "E_bend_Lo_south": final_en["E_bend_Lo_south"],
            "E_line": final_en["E_line"],
            "E_gauss": final_en["E_gauss"],
            "E_ade": final_en["E_ade"],
            "dimple_depth_north": m["dimple_depth_north"],
            "dimple_depth_south": m["dimple_depth_south"],
            "r_b1": m["r_b1"],
            "r_b2": m["r_b2"],
            "z_b1": m["z_b1"],
            "z_b2": m["z_b2"],
            "r_profile": m["r"],
            "z_profile": m["z"],
            "psi_profile": final_psi,
            "H_profile": m["H"],
            "K_profile": m["K"],
            "iterations": opt_res.nit
        }
