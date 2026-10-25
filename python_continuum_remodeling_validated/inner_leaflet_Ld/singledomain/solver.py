"""
Variational Shape Solver for Two-Phase Vesicles (PRL Exact Constrained NLP Standard)
Solves the exact constrained optimization problem over piecewise C^0-joined cubic splines:
  Minimize E_total(psi)
  subject to:
    V(psi) = V_target
    A_Lo(psi) = A_Lo_target
    A_tot(psi) = A_tot_target
    r(S_tot) = 0  (Exact South-pole axis closure)
    psi(0) = 0, psi(s_split^-) = psi(s_split^+) = psi_b, psi(S_tot) = pi
"""

import os
import sys
import numpy as np
from scipy.optimize import minimize
from scipy.interpolate import CubicSpline

sys.path.insert(0, os.path.dirname(__file__))

try:
    from mechanics import (
        MembraneParameters,
        reconstruct_profile,
        compute_curvatures,
        compute_energy_breakdown,
        compute_energy_breakdown_piecewise
    )
except ImportError:
    from sim_prl.mechanics import (
        MembraneParameters,
        reconstruct_profile,
        compute_curvatures,
        compute_energy_breakdown,
        compute_energy_breakdown_piecewise
    )


class PRLVesicleSolver:
    def __init__(self, params=None, n_points=140, n_ctrl_Lo=6, n_ctrl_Ld=10):
        self.params = params if params is not None else MembraneParameters()
        self.n_points = int(n_points)
        self.n_ctrl_Lo = int(n_ctrl_Lo)
        self.n_ctrl_Ld = int(n_ctrl_Ld)

        self.S_tot_0 = self.params.S_tot
        self.s_split_0 = self.params.s_split
        self.S_Ld_0 = self.S_tot_0 - self.s_split_0

        self.S_tot = self.S_tot_0
        self.s_split = self.s_split_0
        self.S_Ld = self.S_Ld_0

        # Sub-grid node counts and normalized grids on [0, 1]
        self.n_pts_Lo = max(35, int(round(self.n_points * 0.35)))
        self.n_pts_Ld = max(65, self.n_points - self.n_pts_Lo + 1)

        self.u_Lo = np.linspace(0.0, 1.0, self.n_pts_Lo)
        self.u_Ld = np.linspace(0.0, 1.0, self.n_pts_Ld)

        self.s_Lo = self.u_Lo * self.s_split_0
        self.s_Ld = self.s_split_0 + self.u_Ld * self.S_Ld_0
        self.s_grid = np.concatenate([self.s_Lo, self.s_Ld[1:]])

        # Normalized spline control nodes on [0, 1] with boundary clustering near u=1 (Lo) and u=0,1 (Ld)
        u_lin_Lo = np.linspace(0.0, 1.0, self.n_ctrl_Lo + 2)[1:]
        self.u_ctrl_Lo = 0.5 * u_lin_Lo + 0.5 * (0.5 * (1.0 - np.cos(np.pi * u_lin_Lo)))
        self.u_ctrl_Lo[-1] = 1.0

        u_lin_Ld = np.linspace(0.0, 1.0, self.n_ctrl_Ld + 2)[1:-1]
        self.u_ctrl_Ld = 0.5 * u_lin_Ld + 0.5 * (0.5 * (1.0 - np.cos(np.pi * u_lin_Ld)))

        self.s_ctrl_Lo = self.u_ctrl_Lo * self.s_split_0
        self.s_ctrl_Ld = self.s_split_0 + self.u_ctrl_Ld * self.S_Ld_0
        self.s_ctrl = np.concatenate([self.s_ctrl_Lo, self.s_ctrl_Ld])
        self.n_angle_params = len(self.s_ctrl)
        self.n_params = self.n_angle_params + 2

        self._cache_x = None
        self._cache_theta = None
        self._cache_energy = None

    def _unpack_and_build_splines(self, x):
        """
        Builds piecewise C^2 splines on Lo [0, S_Lo] and Ld [S_Lo, S_tot]
        joined C^0 at psi_b = x[n_ctrl_Lo], with free meridian arc lengths S_Lo, S_Ld.
        """
        psi_Lo_ctrl = x[:self.n_ctrl_Lo + 1]
        psi_b = psi_Lo_ctrl[-1]
        psi_Ld_ctrl = x[self.n_ctrl_Lo + 1:self.n_angle_params]
        if len(x) >= self.n_angle_params + 2:
            S_Lo = float(x[self.n_angle_params])
            S_Ld = float(x[self.n_angle_params + 1])
        else:
            S_Lo = self.s_split_0
            S_Ld = self.S_Ld_0
        S_tot = S_Lo + S_Ld

        s_Lo = self.u_Lo * S_Lo
        s_Ld = S_Lo + self.u_Ld * S_Ld

        s_nodes_Lo = np.concatenate([[0.0], self.u_ctrl_Lo * S_Lo])
        psi_nodes_Lo = np.concatenate([[0.0], psi_Lo_ctrl])
        cs_Lo = CubicSpline(s_nodes_Lo, psi_nodes_Lo, bc_type=((2, 0.0), 'not-a-knot'))

        s_nodes_Ld = S_Lo + np.concatenate([[0.0], self.u_ctrl_Ld * S_Ld, [S_Ld]])
        psi_nodes_Ld = np.concatenate([[psi_b], psi_Ld_ctrl, [np.pi]])
        cs_Ld = CubicSpline(s_nodes_Ld, psi_nodes_Ld, bc_type=('not-a-knot', (2, 0.0)))

        psi_Lo = cs_Lo(s_Lo)
        c1_Lo = cs_Lo(s_Lo, 1)

        psi_Ld = cs_Ld(s_Ld)
        c1_Ld = cs_Ld(s_Ld, 1)

        return s_Lo, psi_Lo, c1_Lo, s_Ld, psi_Ld, c1_Ld

    def _param_to_psi(self, x):
        _, psi_Lo, _, _, psi_Ld, _ = self._unpack_and_build_splines(x)
        return np.concatenate([psi_Lo, psi_Ld[1:]])

    def _eval_cached(self, x, theta):
        if (self._cache_x is not None and np.array_equal(x, self._cache_x) and
                self._cache_theta == theta):
            return self._cache_energy

        s_Lo, psi_Lo, c1_Lo, s_Ld, psi_Ld, c1_Ld = self._unpack_and_build_splines(x)
        breakdown = compute_energy_breakdown_piecewise(
            s_Lo, psi_Lo, c1_Lo, s_Ld, psi_Ld, c1_Ld, self.params, theta
        )
        self._cache_x = np.copy(x)
        self._cache_theta = theta
        self._cache_energy = breakdown
        return breakdown

    def get_initial_seed(self, branch="inward"):
        """Generates physical initial seed profiles for designated morphology branch."""
        u_Lo = self.u_ctrl_Lo
        u_Ld = self.u_ctrl_Ld
        psi_b_sphere = (self.s_split_0 / self.S_tot_0) * np.pi
        scale = min(1.0, max(0.2, self.s_split_0 / 1.5))

        if branch == "sphere":
            psi_Lo_ctrl = u_Lo * psi_b_sphere
            psi_b = psi_Lo_ctrl[-1]
            psi_Ld_ctrl = psi_b + u_Ld * (np.pi - psi_b)
        elif branch == "inward":
            # Inverted cap seed: psi < 0 throughout Lo cap (including psi_b < 0 at boundary)
            psi_Lo_ctrl = -0.85 * u_Lo * psi_b_sphere - 0.35 * scale * np.sin(np.pi * u_Lo)
            psi_b = psi_Lo_ctrl[-1]
            psi_Ld_ctrl = psi_b + u_Ld * (np.pi - psi_b)
        elif branch == "outward":
            # Outward dome seed: steepened positive angle in Lo cap
            psi_Lo_ctrl = u_Lo * psi_b_sphere + 0.5 * scale * np.sin(np.pi * u_Lo)
            psi_b = psi_Lo_ctrl[-1]
            psi_Ld_ctrl = psi_b + u_Ld * (np.pi - psi_b)
        else:
            psi_Lo_ctrl = u_Lo * psi_b_sphere
            psi_b = psi_Lo_ctrl[-1]
            psi_Ld_ctrl = psi_b + u_Ld * (np.pi - psi_b)

        return np.concatenate([psi_Lo_ctrl, psi_Ld_ctrl, [self.s_split_0, self.S_Ld_0]])

    def optimize_shape(self, theta=0.0, branch="inward", maxiter=500, ftol=1e-9, method="SLSQP"):
        """
        Performs exact constrained minimization of physical energy.
        If branch == 'auto', evaluates inward, sphere, and outward seeds and returns the global minimum.
        """
        if branch == "auto":
            branches_to_try = ["inward", "sphere", "outward"] if theta > 0 else ["sphere", "outward"]
            best_res = None
            inward_res = None
            for b in branches_to_try:
                res = self._optimize_single_branch(theta=theta, branch=b, maxiter=maxiter, ftol=ftol, method=method)
                if b == "inward":
                    inward_res = res
                if best_res is None:
                    best_res = res
                elif res['geom_converged'] and not best_res['geom_converged']:
                    best_res = res
                elif res['geom_converged'] and best_res['geom_converged'] and res['total_energy'] < best_res['total_energy']:
                    best_res = res
                elif not best_res['geom_converged'] and not res['geom_converged'] and res['max_geom_error'] < best_res['max_geom_error']:
                    best_res = res
            best_res = dict(best_res)
            best_res["inward_branch_result"] = inward_res if inward_res is not None else best_res
            return best_res
        else:
            res = self._optimize_single_branch(theta=theta, branch=branch, maxiter=maxiter, ftol=ftol, method=method)
            res["inward_branch_result"] = res
            return res

    def _optimize_single_branch(self, theta=0.0, branch="inward", maxiter=500, ftol=1e-9, method="SLSQP"):
        if self.params.reduced_volume >= 0.9999:
            x_sph = self.get_initial_seed(branch="sphere")
            final_metrics = dict(self._eval_cached(x_sph, theta))
            final_metrics['res_V'] = 0.0
            final_metrics['res_A_Lo'] = 0.0
            final_metrics['res_A_tot'] = 0.0
            final_metrics['res_r_pole'] = 0.0
            final_metrics['converged'] = True
            final_metrics['opt_success'] = True
            final_metrics['geom_converged'] = True
            final_metrics['opt_result'] = None
            final_metrics['solver_message'] = 'Exact spherical state (v=1.00)'
            final_metrics['solver_status'] = 0
            final_metrics['solver_method'] = method
            final_metrics['max_geom_error'] = 0.0
            final_metrics['branch'] = 'sphere'
            final_metrics['theta'] = theta
            final_metrics['iterations'] = 0
            return final_metrics

        x0 = self.get_initial_seed(branch=branch)
        self._cache_x = None
        self._cache_theta = None
        self._cache_energy = None

        def objective(x):
            b = self._eval_cached(x, theta)
            return b['total_energy']

        def constr_volume(x):
            b = self._eval_cached(x, theta)
            return b['res_V']

        def constr_area_Lo(x):
            b = self._eval_cached(x, theta)
            return b['res_A_Lo']

        def constr_area_tot(x):
            b = self._eval_cached(x, theta)
            return b['res_A_tot']

        def constr_r_pole(x):
            b = self._eval_cached(x, theta)
            return b['res_r_pole']

        constraints = [
            {'type': 'eq', 'fun': constr_volume},
            {'type': 'eq', 'fun': constr_area_Lo},
            {'type': 'eq', 'fun': constr_area_tot},
            {'type': 'eq', 'fun': constr_r_pole}
        ]

        bounds = (
            [(-0.5 * np.pi, 1.5 * np.pi)] * self.n_angle_params +
            [(0.98 * self.s_split_0, 1.18 * self.s_split_0), (0.95 * self.S_Ld_0, 1.10 * self.S_Ld_0)]
        )

        options = {
            'maxiter': int(maxiter),
            'ftol': float(ftol),
            'disp': False
        }

        res_opt = minimize(
            objective,
            x0,
            method=method,
            bounds=bounds,
            constraints=constraints,
            options=options
        )

        final_metrics = dict(self._eval_cached(res_opt.x, theta))

        max_geom_error = max(
            abs(final_metrics['res_V']),
            abs(final_metrics['res_A_Lo']),
            abs(final_metrics['res_A_tot']),
            abs(final_metrics['res_r_pole'])
        )

        opt_success = bool(res_opt.success or res_opt.status == 0)
        geom_converged = bool(max_geom_error < 1.0e-3 and final_metrics['min_r_interior'] > 0.0)
        is_converged = bool(geom_converged and (opt_success or max_geom_error < 1.0e-5))

        final_metrics['converged'] = is_converged
        final_metrics['opt_success'] = opt_success
        final_metrics['geom_converged'] = geom_converged
        final_metrics['opt_result'] = res_opt
        final_metrics['solver_message'] = str(res_opt.message)
        final_metrics['solver_status'] = int(res_opt.status) if hasattr(res_opt, 'status') else 0
        final_metrics['solver_method'] = method
        final_metrics['max_geom_error'] = max_geom_error
        final_metrics['branch'] = branch
        final_metrics['theta'] = theta
        final_metrics['iterations'] = int(res_opt.nit) if hasattr(res_opt, 'nit') else 0

        return final_metrics
