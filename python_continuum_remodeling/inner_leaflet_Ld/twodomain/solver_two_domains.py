"""
Exact Non-Linear Constrained Solver for Two Opposite Lo Domains (Bipolar Vesicle - PRL Standard).
Uses piecewise C^0-joined cubic splines on North Lo [0, S_Lo] and northern half-Ld [S_Lo, S_eq]
with free domain arc lengths (S_Lo, S_Ld_half) and exact up-down reflection symmetry:
    psi(S_tot - s) = pi - psi(s),   c1(S_tot - s) = c1(s)
ensuring that both North and South Lo domains behave identically and symmetrically by construction:
    depth_north == depth_south
    r_b1 == r_b2
    A_Lo1 == A_Lo2
    r(S_tot) == 0 (Exact South-pole closure by antisymmetry of cos(psi))
"""

import os
import sys
import numpy as np
from scipy.optimize import minimize
from scipy.interpolate import CubicSpline

sys.path.insert(0, os.path.dirname(__file__))
try:
    from mechanics_two_domains import TwoDomainParameters, TwoDomainEnergyCalculator
except ImportError:
    from sim_prl.mechanics_two_domains import TwoDomainParameters, TwoDomainEnergyCalculator


class TwoDomainVesicleSolver:
    """
    NLP Shape Solver for 3-zone Bipolar Vesicle with piecewise C^0-joined cubic splines
    and exact mirror symmetry.
    """
    def __init__(self, params: TwoDomainParameters = None, n_points: int = 101, n_ctrl_Lo: int = 5, n_ctrl_Ld_half: int = 6):
        if n_points % 2 == 0:
            n_points += 1
        self.params = params if params is not None else TwoDomainParameters()
        self.p = self.params
        self.calc = TwoDomainEnergyCalculator(self.p)
        self.n_points = int(n_points)
        self.n_ctrl_Lo = int(n_ctrl_Lo)
        self.n_ctrl_Ld_half = int(n_ctrl_Ld_half)

        # Reference spherical arc lengths
        self.S_tot_0 = np.pi * self.p.R_ves
        cos_theta_cap = max(-1.0, min(1.0, 1.0 - self.p.A_Lo_target / (2.0 * np.pi * (self.p.R_ves ** 2))))
        theta_cap = np.arccos(cos_theta_cap)
        self.s_split1_0 = theta_cap * self.p.R_ves
        self.S_Ld_half_0 = 0.5 * self.S_tot_0 - self.s_split1_0

        # Sub-grid node counts
        self.n_pts_Lo = max(25, int(round(self.n_points * 0.25)))
        self.n_pts_Ld_half = max(30, (self.n_points // 2) - self.n_pts_Lo + 2)

        self.u_grid_Lo = np.linspace(0.0, 1.0, self.n_pts_Lo)
        self.u_grid_Ld_half = np.linspace(0.0, 1.0, self.n_pts_Ld_half)

        # Normalized spline control nodes on [0, 1]
        # Lo: interior nodes + shared boundary u=1.0
        self.u_ctrl_Lo = np.linspace(0.0, 1.0, self.n_ctrl_Lo + 2)[1:]  # length n_ctrl_Lo + 1
        # Ld_half: interior nodes on (0, 1) (u=0 is psi_b1, u=1 is equator psi=pi/2)
        self.u_ctrl_Ld_half = np.linspace(0.0, 1.0, self.n_ctrl_Ld_half + 2)[1:-1]  # length n_ctrl_Ld_half

        self.n_angle_params = (self.n_ctrl_Lo + 1) + self.n_ctrl_Ld_half
        self.n_params = self.n_angle_params

        # Public attributes
        self.s_split1 = self.s_split1_0
        self.s_split2 = self.S_tot_0 - self.s_split1_0
        self.S_tot = self.S_tot_0
        self.idx_b1 = self.n_pts_Lo - 1
        self.idx_b2 = self.n_pts_Lo + 2 * (self.n_pts_Ld_half - 1) - 1
        self.mid_idx = self.n_pts_Lo + self.n_pts_Ld_half - 2

        # Initialize fixed piecewise grids
        self.s_Lo1 = self.u_grid_Lo * self.s_split1_0
        self.s_Ld_n = self.s_split1_0 + self.u_grid_Ld_half * self.S_Ld_half_0
        s_eq = self.s_Ld_n[-1]
        self.s_Ld_s = 2.0 * s_eq - self.s_Ld_n[-2::-1]
        self.s_Ld = np.concatenate([self.s_Ld_n, self.s_Ld_s])
        self.s_Lo2 = 2.0 * s_eq - self.s_Lo1[::-1]
        self.s_grid = np.concatenate([self.s_Lo1, self.s_Ld[1:], self.s_Lo2[1:]])
        self.s_half = np.concatenate([self.s_Lo1, self.s_Ld_n[1:]])

        self._last_x = None
        self._cached_theta = None
        self._cached_metrics = None
        self._cached_energy = None

    def _unpack_and_build_zones(self, x: np.ndarray):
        """
        Builds piecewise C^2 splines on North Lo [0, S_Lo] and northern half-Ld [S_Lo, S_eq]
        joined C^0 at psi_b1, then reflects with exact up-down antisymmetry to South Lo.
        """
        psi_Lo_ctrl = x[:self.n_ctrl_Lo + 1]
        psi_b1 = psi_Lo_ctrl[-1]
        psi_Ld_ctrl = x[self.n_ctrl_Lo + 1:self.n_angle_params]
        S_Lo = self.s_split1_0
        S_Ld_half = self.S_Ld_half_0

        # Zone 1: North Lo [0, S_Lo] with psi''(0)=0 and free boundary slope/curvature at S_Lo
        s_nodes_Lo = np.concatenate([[0.0], self.u_ctrl_Lo * S_Lo])
        psi_nodes_Lo = np.concatenate([[0.0], psi_Lo_ctrl])
        cs_Lo = CubicSpline(s_nodes_Lo, psi_nodes_Lo, bc_type=((2, 0.0), 'not-a-knot'))

        s_Lo1 = self.s_Lo1
        psi_Lo1 = cs_Lo(s_Lo1)
        c1_Lo1 = cs_Lo(s_Lo1, 1)

        # Zone 2 (North half): [S_Lo, S_Lo + S_Ld_half] with free slope at S_Lo and psi''(S_eq)=0
        s_nodes_Ld_n = S_Lo + np.concatenate([[0.0], self.u_ctrl_Ld_half * S_Ld_half, [S_Ld_half]])
        psi_nodes_Ld_n = np.concatenate([[psi_b1], psi_Ld_ctrl, [0.5 * np.pi]])
        cs_Ld_n = CubicSpline(s_nodes_Ld_n, psi_nodes_Ld_n, bc_type=('not-a-knot', (2, 0.0)))

        s_Ld_n = self.s_Ld_n
        psi_Ld_n = cs_Ld_n(s_Ld_n)
        c1_Ld_n = cs_Ld_n(s_Ld_n, 1)

        # Exact reflection across equator s_eq = S_Lo + S_Ld_half
        s_Ld_s = self.s_Ld_s
        psi_Ld_s = np.pi - psi_Ld_n[-2::-1]
        c1_Ld_s = c1_Ld_n[-2::-1]

        s_Ld = self.s_Ld
        psi_Ld = np.concatenate([psi_Ld_n, psi_Ld_s])
        c1_Ld = np.concatenate([c1_Ld_n, c1_Ld_s])

        # Zone 3: South Lo [S_tot - S_Lo, S_tot]
        s_Lo2 = self.s_Lo2
        psi_Lo2 = np.pi - psi_Lo1[::-1]
        c1_Lo2 = c1_Lo1[::-1]

        return s_Lo1, psi_Lo1, c1_Lo1, s_Ld, psi_Ld, c1_Ld, s_Lo2, psi_Lo2, c1_Lo2

    def reconstruct_full_psi(self, x_interior: np.ndarray) -> np.ndarray:
        """
        Reconstructs full meridian psi(s) from parameter vector x_interior.
        Supports both spline parameter vectors (len == n_params) and legacy half-grid vectors.
        """
        if len(x_interior) == self.n_params:
            _, psi_Lo1, _, _, psi_Ld, _, _, psi_Lo2, _ = self._unpack_and_build_zones(x_interior)
            return np.concatenate([psi_Lo1, psi_Ld[1:], psi_Lo2[1:]])
        ph = np.zeros(len(x_interior) + 2)
        ph[0] = 0.0
        ph[1:-1] = x_interior
        ph[-1] = np.pi / 2.0
        return np.concatenate([ph, np.pi - ph[-2::-1]])

    def _eval_cached(self, x: np.ndarray, theta: float):
        if (self._last_x is not None and np.array_equal(x, self._last_x) and
                self._cached_theta == theta):
            return self._cached_energy, self._cached_metrics

        zones = self._unpack_and_build_zones(x)
        en_res = self.calc.compute_energy_piecewise(*zones, theta=theta)
        self._last_x = np.copy(x)
        self._cached_theta = theta
        self._cached_energy = en_res
        self._cached_metrics = en_res["metrics"]
        return en_res, en_res["metrics"]

    def generate_initial_guess(self, branch: str = "inward") -> np.ndarray:
        """
        Generates initial spline parameter vector for the symmetric 2-domain vesicle.
        """
        S_Lo_init = self.s_split1_0
        S_tot_init = self.S_tot_0

        psi_b_sphere = (S_Lo_init / S_tot_init) * np.pi
        scale = min(1.0, max(0.2, S_Lo_init / 1.5))

        if branch in ["sphere", "bipolar_sphere"]:
            psi_Lo_ctrl = self.u_ctrl_Lo * psi_b_sphere
            psi_b1 = psi_Lo_ctrl[-1]
            psi_Ld_ctrl = psi_b1 + self.u_ctrl_Ld_half * (0.5 * np.pi - psi_b1)
        elif branch in ["inward", "bipolar_inward"]:
            psi_Lo_ctrl = -0.85 * self.u_ctrl_Lo * psi_b_sphere - 0.35 * scale * np.sin(np.pi * self.u_ctrl_Lo)
            psi_b1 = psi_Lo_ctrl[-1]
            psi_Ld_ctrl = psi_b1 + self.u_ctrl_Ld_half * (0.5 * np.pi - psi_b1)
        elif branch in ["outward", "bipolar_outward"]:
            psi_Lo_ctrl = 1.15 * self.u_ctrl_Lo * psi_b_sphere + 0.35 * scale * np.sin(np.pi * self.u_ctrl_Lo)
            psi_b1 = psi_Lo_ctrl[-1]
            psi_Ld_ctrl = psi_b1 + self.u_ctrl_Ld_half * (0.5 * np.pi - psi_b1)
        else:
            psi_Lo_ctrl = self.u_ctrl_Lo * psi_b_sphere
            psi_b1 = psi_Lo_ctrl[-1]
            psi_Ld_ctrl = psi_b1 + self.u_ctrl_Ld_half * (0.5 * np.pi - psi_b1)

        return np.concatenate([psi_Lo_ctrl, psi_Ld_ctrl])

    def optimize_shape(self, theta: float, branch: str = "inward", maxiter: int = 400, ftol: float = 1e-9):
        """
        Performs exact non-linear constrained optimization for the symmetric 2-domain vesicle.
        Supports branch='auto' to select the global minimum energy across inward, sphere, and outward seeds.
        """
        if branch == "auto":
            branches_to_try = ["inward", "sphere", "outward"] if theta > 0 else ["sphere"]
            best_res = None
            inward_res = None
            for b in branches_to_try:
                res = self._optimize_single_branch(theta=theta, branch=b, maxiter=maxiter, ftol=ftol)
                if b == "inward":
                    inward_res = res
                if best_res is None or (res['converged'] and (not best_res['converged'] or res['total_energy'] < best_res['total_energy'])):
                    best_res = res
                elif not best_res['converged'] and res['opt_success'] and res['total_energy'] < best_res['total_energy']:
                    best_res = res
            best_res = dict(best_res)
            best_res["inward_branch_result"] = inward_res if inward_res is not None else best_res
            return best_res
        else:
            res = self._optimize_single_branch(theta=theta, branch=branch, maxiter=maxiter, ftol=ftol)
            res["inward_branch_result"] = res
            return res

    def _optimize_single_branch(self, theta: float, branch: str = "inward", maxiter: int = 400, ftol: float = 1e-9):
        self._last_x = None
        self._cached_theta = None
        self._cached_metrics = None
        self._cached_energy = None

        x0 = self.generate_initial_guess(branch)

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

        bounds = [(-0.5 * np.pi, 1.5 * np.pi)] * self.n_angle_params

        opt_res = minimize(
            obj_func,
            x0,
            method="SLSQP",
            bounds=bounds,
            constraints=constraints,
            options={"maxiter": int(maxiter), "ftol": float(ftol), "disp": False}
        )

        final_en, m = self._eval_cached(opt_res.x, theta)

        res_V = abs((m["V_total"] - self.p.V_tot_target) / self.p.V_tot_target)
        res_A_tot = abs((m["A_total"] - self.p.A_tot_target) / self.p.A_tot_target)
        res_A_Lo = abs((m["A_Lo_north"] - self.p.A_Lo_target) / self.p.A_Lo_target)
        res_r_pole = abs(float(m["r"][-1]) / self.p.R_ves)
        max_err = max(res_V, res_A_tot, res_A_Lo, res_r_pole)

        opt_success = bool(opt_res.success or opt_res.status == 0)
        geom_converged = bool(max_err < 1e-3 and m["min_r_interior"] > 0.0)
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
            "res_r_pole": res_r_pole,
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
            "peak_to_dimple_Lo_north": m["peak_to_dimple_Lo_north"],
            "peak_to_dimple_crest_north": m["peak_to_dimple_crest_north"],
            "dome_peak_to_dimple_north": m["dome_peak_to_dimple_north"],
            "dome_peak_to_dimple_south": m["dome_peak_to_dimple_south"],
            "r_b1": m["r_b1"],
            "r_b2": m["r_b2"],
            "z_b1": m["z_b1"],
            "z_b2": m["z_b2"],
            "r_profile": m["r"],
            "z_profile": m["z"],
            "psi_profile": m["psi"],
            "H_profile": m["H"],
            "K_profile": m["K"],
            "branch": branch,
            "iterations": int(opt_res.nit) if hasattr(opt_res, "nit") else 0
        }
