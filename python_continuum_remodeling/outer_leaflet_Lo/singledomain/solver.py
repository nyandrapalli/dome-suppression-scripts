"""
Variational Shape Solver for HermesLoouter Single-Domain (PRL Standard).
Solves exact constrained optimization:
  Minimize E_total(psi)
  subject to:
    V(psi) = V_target
    A_Lo(psi) = A_Lo_target
    A_tot(psi) = A_tot_target
    psi(0) = 0, psi(S_tot) = pi
"""

import sys
import os
import numpy as np
from scipy.optimize import minimize
from scipy.interpolate import CubicSpline

sys.path.insert(0, os.path.dirname(__file__))
from mechanics import (
    MembraneParameters,
    reconstruct_profile,
    compute_curvatures,
    compute_energy_breakdown
)

class HermesSingleSolver:
    def __init__(self, params=None, n_points=140, n_ctrl_Lo=6, n_ctrl_Ld=10):
        self.params = params if params is not None else MembraneParameters()
        self.n_points = int(n_points)
        self.n_ctrl_Lo = int(n_ctrl_Lo)
        self.n_ctrl_Ld = int(n_ctrl_Ld)
        
        self.S_tot = self.params.S_tot
        self.s_split = self.params.s_split
        self.s_grid = np.linspace(0.0, self.S_tot, self.n_points)
        
        # Internal control nodes for B-spline representation
        self.s_ctrl_Lo = np.linspace(0.0, self.s_split, self.n_ctrl_Lo + 2)[1:]
        self.s_ctrl_Ld = np.linspace(self.s_split, self.S_tot, self.n_ctrl_Ld + 2)[1:-1]
        self.s_ctrl = np.unique(np.concatenate([self.s_ctrl_Lo, self.s_ctrl_Ld]))
        self.n_params = len(self.s_ctrl)
        
        # Cache
        self._cache_x = None
        self._cache_theta = None
        self._cache_energy = None

    def _param_to_psi(self, x):
        """Converts control node angles x to full arc-length profile psi(s)."""
        s_all = np.concatenate([[0.0], self.s_ctrl, [self.S_tot]])
        psi_all = np.concatenate([[0.0], x, [np.pi]])
        cs = CubicSpline(s_all, psi_all, bc_type='natural')
        return cs(self.s_grid)

    def _eval_cached(self, x, theta):
        if (self._cache_x is not None and np.array_equal(x, self._cache_x) and 
            self._cache_theta == theta):
            return self._cache_energy
        
        psi_grid = self._param_to_psi(x)
        breakdown = compute_energy_breakdown(
            self.s_grid, psi_grid, self.s_split, self.params, theta
        )
        self._cache_x = np.copy(x)
        self._cache_theta = theta
        self._cache_energy = breakdown
        return breakdown

    def get_initial_seed(self, branch="outward"):
        """Generates initial seed profiles for designated morphology branch."""
        psi_sphere = self.s_grid / self.params.R_ves
        
        if branch == "sphere":
            psi_seed = psi_sphere
        elif branch == "outward":
            # Outward dome / bud seed (promoted by positive c0 on Lo outer leaflet)
            scale = min(1.0, max(0.2, self.s_split / 1.5))
            psi_seed = np.copy(psi_sphere)
            mask_Lo = (self.s_grid <= self.s_split)
            mask_Ld = (self.s_grid > self.s_split)
            psi_seed[mask_Lo] = 0.5 * scale * np.sin(np.pi * self.s_grid[mask_Lo] / self.s_split) + (self.s_grid[mask_Lo] / self.s_split) * (self.s_split / self.S_tot * np.pi)
            psi_at_b = psi_seed[mask_Lo][-1]
            psi_seed[mask_Ld] = psi_at_b + ((self.s_grid[mask_Ld] - self.s_split) / (self.S_tot - self.s_split)) * (np.pi - psi_at_b)
        elif branch == "inward":
            # Inward dimple seed
            scale = min(1.0, max(0.2, self.s_split / 1.5))
            psi_seed = np.copy(psi_sphere)
            mask_Lo = (self.s_grid <= self.s_split)
            mask_Ld = (self.s_grid > self.s_split)
            psi_seed[mask_Lo] = -0.5 * scale * np.sin(np.pi * self.s_grid[mask_Lo] / self.s_split) + (self.s_grid[mask_Lo] / self.s_split) * (self.s_split / self.S_tot * np.pi)
            psi_at_b = psi_seed[mask_Lo][-1]
            psi_seed[mask_Ld] = psi_at_b + ((self.s_grid[mask_Ld] - self.s_split) / (self.S_tot - self.s_split)) * (np.pi - psi_at_b)
        else:
            psi_seed = psi_sphere
            
        cs_seed = CubicSpline(self.s_grid, psi_seed)
        return cs_seed(self.s_ctrl)

    def optimize_shape(self, theta=0.0, branch="auto", maxiter=350, ftol=1e-9, method="SLSQP"):
        """
        Performs exact constrained minimization.
        If branch == 'auto', runs both outward and sphere/inward seeds and selects the lowest energy solution.
        """
        if branch == "auto":
            branches_to_try = ["outward", "sphere", "inward"] if theta > 0 else ["sphere"]
            best_res = None
            for b in branches_to_try:
                res = self._optimize_single_branch(theta, branch=b, maxiter=maxiter, ftol=ftol, method=method)
                if best_res is None or (res['converged'] and (not best_res['converged'] or res['total_energy'] < best_res['total_energy'])):
                    best_res = res
                elif not best_res['converged'] and res['opt_success'] and res['total_energy'] < best_res['total_energy']:
                    best_res = res
            return best_res
        else:
            return self._optimize_single_branch(theta, branch=branch, maxiter=maxiter, ftol=ftol, method=method)

    def _optimize_single_branch(self, theta=0.0, branch="outward", maxiter=350, ftol=1e-9, method="SLSQP"):
        x0 = self.get_initial_seed(branch=branch)
        
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

        constraints = [
            {'type': 'eq', 'fun': constr_volume},
            {'type': 'eq', 'fun': constr_area_Lo},
            {'type': 'eq', 'fun': constr_area_tot}
        ]
        
        bounds = [(- 0.5 * np.pi, 1.5 * np.pi)] * self.n_params
        
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
        
        psi_final = self._param_to_psi(res_opt.x)
        final_metrics = compute_energy_breakdown(
            self.s_grid, psi_final, self.s_split, self.params, theta
        )
        
        max_geom_error = max(
            abs(final_metrics['res_V']),
            abs(final_metrics['res_A_Lo']),
            abs(final_metrics['res_A_tot'])
        )
        
        opt_success = bool(res_opt.success or res_opt.status == 0)
        geom_converged = bool(max_geom_error < 1.0e-3)
        is_converged = bool(opt_success and geom_converged)
        
        final_metrics['converged'] = is_converged
        final_metrics['opt_success'] = opt_success
        final_metrics['geom_converged'] = geom_converged
        final_metrics['psi'] = psi_final
        final_metrics['opt_result'] = res_opt
        final_metrics['solver_message'] = str(res_opt.message)
        final_metrics['solver_status'] = int(res_opt.status) if hasattr(res_opt, 'status') else 0
        final_metrics['solver_method'] = method
        final_metrics['max_geom_error'] = max_geom_error
        final_metrics['branch'] = branch
        final_metrics['theta'] = theta
        final_metrics['iterations'] = int(res_opt.nit) if hasattr(res_opt, 'nit') else 0
        
        return final_metrics
