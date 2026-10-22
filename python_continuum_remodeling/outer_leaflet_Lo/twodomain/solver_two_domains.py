"""
Exact Non-Linear Constrained Solver for Two Opposite Lo Domains (HermesLoouter).
Enforces exact physical up-down reflection symmetry:
    psi(S_tot - s) = pi - psi(s)
ensuring that both North and South Lo domains behave identically:
    elevation_north == elevation_south
    depth_north == depth_south
    r_b1 == r_b2
"""

import os
import sys
import numpy as np
from scipy.optimize import minimize
from scipy.interpolate import CubicSpline

sys.path.insert(0, os.path.dirname(__file__))
from mechanics_two_domains import TwoDomainParameters, TwoDomainEnergyCalculator

class TwoDomainHermesSolver:
    def __init__(self, params: TwoDomainParameters = None, n_points: int = 101):
        if n_points % 2 == 0:
            n_points += 1
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

    def generate_initial_guess(self, branch: str = "outward") -> np.ndarray:
        s_h = self.s_half
        s_b1 = self.s_grid[self.idx_b1]
        s_eq = self.s_grid[self.mid_idx]
        
        mask_lo = (s_h <= s_b1)
        mask_ld = (s_h > s_b1)
        
        scale = min(1.0, max(0.2, self.s_split1 / 1.5))
        psi_h = np.zeros_like(s_h)
        
        if branch in ["sphere", "bipolar_sphere"]:
            psi_h = (s_h / self.S_tot) * np.pi
        elif branch in ["outward", "bipolar_outward"]:
            # Outward dome seed
            psi_h[mask_lo] = 0.5 * scale * np.sin(np.pi * s_h[mask_lo] / s_b1) + (s_h[mask_lo] / s_b1) * (s_b1 / self.S_tot * np.pi)
            psi_at_b1 = psi_h[mask_lo][-1]
            psi_h[mask_ld] = psi_at_b1 + ((s_h[mask_ld] - s_b1) / (s_eq - s_b1)) * (np.pi / 2.0 - psi_at_b1)
        elif branch in ["inward", "bipolar_inward"]:
            # Inward dimple seed
            psi_h[mask_lo] = -0.5 * scale * np.sin(np.pi * s_h[mask_lo] / s_b1) + (s_h[mask_lo] / s_b1) * (s_b1 / self.S_tot * np.pi)
            psi_at_b1 = psi_h[mask_lo][-1]
            psi_h[mask_ld] = psi_at_b1 + ((s_h[mask_ld] - s_b1) / (s_eq - s_b1)) * (np.pi / 2.0 - psi_at_b1)
        else:
            psi_h = (s_h / self.S_tot) * np.pi
            
        return psi_h[1:-1]

    def optimize_shape(self, theta: float = 0.0, branch: str = "auto", maxiter: int = 350, ftol: float = 1e-9) -> dict:
        if branch == "auto":
            branches_to_try = ["outward", "sphere", "inward"] if theta > 0 else ["sphere"]
            best_res = None
            for b in branches_to_try:
                res = self._optimize_single_branch(theta, branch=b, maxiter=maxiter, ftol=ftol)
                if best_res is None or (res['converged'] and (not best_res['converged'] or res['total_energy'] < best_res['total_energy'])):
                    best_res = res
                elif not best_res['converged'] and res['opt_success'] and res['total_energy'] < best_res['total_energy']:
                    best_res = res
            return best_res
        else:
            return self._optimize_single_branch(theta, branch=branch, maxiter=maxiter, ftol=ftol)

    def _optimize_single_branch(self, theta: float = 0.0, branch: str = "outward", maxiter: int = 350, ftol: float = 1e-9) -> dict:
        x0 = self.generate_initial_guess(branch=branch)
        
        def obj_func(x):
            en, _ = self._eval_cached(x, theta)
            return en["E_total"]
            
        def constr_area_tot(x):
            _, m = self._eval_cached(x, theta)
            return (m["A_total"] - self.p.A_tot_target) / self.p.A_tot_target
            
        def constr_vol(x):
            _, m = self._eval_cached(x, theta)
            return (m["V_total"] - self.p.V_tot_target) / self.p.V_tot_target
            
        def constr_area_lo_north(x):
            _, m = self._eval_cached(x, theta)
            return (m["A_Lo_north"] - self.p.A_Lo_target) / self.p.A_Lo_target

        constraints = [
            {"type": "eq", "fun": constr_area_tot},
            {"type": "eq", "fun": constr_vol},
            {"type": "eq", "fun": constr_area_lo_north}
        ]
        
        bounds = [(-np.pi * 0.5, np.pi * 1.5)] * len(x0)
        options = {"maxiter": int(maxiter), "ftol": float(ftol), "disp": False}
        
        res = minimize(obj_func, x0, method="SLSQP", bounds=bounds, constraints=constraints, options=options)
        
        psi_opt = self.reconstruct_full_psi(res.x)
        en_final = self.calc.compute_energy(self.s_grid, psi_opt, self.idx_b1, self.idx_b2, theta)
        m_final = en_final["metrics"]
        
        res_A = abs((m_final["A_total"] - self.p.A_tot_target) / self.p.A_tot_target)
        res_V = abs((m_final["V_total"] - self.p.V_tot_target) / self.p.V_tot_target)
        res_Lo = abs((m_final["A_Lo_north"] - self.p.A_Lo_target) / self.p.A_Lo_target)
        max_geom_err = max(res_A, res_V, res_Lo)
        
        opt_success = bool(res.success or res.status == 0)
        geom_converged = bool(max_geom_err < 1e-3)
        converged = bool(opt_success and geom_converged)
        
        return {
            "converged": converged,
            "opt_success": opt_success,
            "geom_converged": geom_converged,
            "max_geom_error": max_geom_err,
            "total_energy": en_final["E_total"],
            "E_bend_tot": en_final["E_bend_tot"],
            "E_bend_Lo_north": en_final["E_bend_Lo_north"],
            "E_bend_Ld": en_final["E_bend_Ld"],
            "E_bend_Lo_south": en_final["E_bend_Lo_south"],
            "E_line": en_final["E_line"],
            "E_gauss": en_final["E_gauss"],
            "E_ade": en_final["E_ade"],
            "dimple_depth_north": m_final["dimple_depth_north"],
            "dimple_depth_south": m_final["dimple_depth_south"],
            "elevation_north": m_final["elevation_north"],
            "elevation_south": m_final["elevation_south"],
            "r_b1": m_final["r_b1"],
            "r_b2": m_final["r_b2"],
            "psi_profile": psi_opt,
            "r_profile": m_final["r"],
            "z_profile": m_final["z"],
            "H_profile": m_final["H"],
            "K_profile": m_final["K"],
            "branch": branch,
            "iterations": int(res.nit) if hasattr(res, 'nit') else 0
        }
