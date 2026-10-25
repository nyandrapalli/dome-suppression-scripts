"""
Variational Shape Solver for Two-Phase Vesicles with Streptavidin Remodeling.
Minimizes the Jülicher-Lipowsky + ADE free energy functional across competing morphological branches.
"""

import numpy as np
from scipy.optimize import minimize
from scipy.interpolate import CubicSpline
from sim.mechanics import (
    MembraneParameters,
    reconstruct_profile,
    compute_curvatures,
    compute_energy_and_constraints
)

class VesicleShapeSolver:
    def __init__(self, params=None, n_points=120):
        self.params = params if params is not None else MembraneParameters()
        self.n_points = n_points
        
        # Estimate total arc length for sphere of radius R_ves
        self.S_tot = np.pi * self.params.R_ves
        
        # S_split for Lo domain with area A_Lo_target ~ pi * R_Lo^2
        theta_split = 2.0 * np.arcsin(self.params.R_Lo / (2.0 * self.params.R_ves))
        self.s_split = self.params.R_ves * theta_split
        
        # Discretization grid along arc-length s with high resolution at boundary and Lo cap
        s_Lo = np.linspace(0, self.s_split, 45)
        s_Ld = np.linspace(self.s_split, self.S_tot, self.n_points - 45 + 1)[1:]
        self.s_grid = np.concatenate([s_Lo, s_Ld])
        self.idx_split = len(s_Lo) - 1
        
    def generate_initial_guess(self, branch="sphere"):
        """
        Generates initial psi(s) profiles for:
        - "sphere": standard spherical profile psi(s) = (s / S_tot) * pi
        - "inward": Lo domain indented inward (invagination with negative psi)
        - "outward": Lo domain bulging outward (budded with steeper positive psi)
        """
        s = self.s_grid
        s_b = self.s_split
        mask_Lo = s <= s_b
        mask_Ld = s > s_b
        
        if branch == "sphere":
            psi = (s / self.S_tot) * np.pi
        elif branch == "inward":
            # Inward dimple: psi starts negative (inward curvature H < 0)
            psi = np.zeros_like(s)
            # In Lo domain: inverted sine wave dipping negative, then matching boundary slope
            psi[mask_Lo] = - 1.2 * np.sin(np.pi * s[mask_Lo] / s_b) + (s[mask_Lo] / s_b) * (s_b / self.S_tot * np.pi)
            psi_at_b = psi[mask_Lo][-1]
            psi[mask_Ld] = psi_at_b + ((s[mask_Ld] - s_b) / (self.S_tot - s_b)) * (np.pi - psi_at_b)
        elif branch == "outward":
            # Outward bud: steeper positive slope in Lo
            psi = np.zeros_like(s)
            psi[mask_Lo] = 1.4 * np.sin(np.pi * s[mask_Lo] / s_b) + (s[mask_Lo] / s_b) * (s_b / self.S_tot * np.pi)
            psi_at_b = psi[mask_Lo][-1]
            psi[mask_Ld] = psi_at_b + ((s[mask_Ld] - s_b) / (self.S_tot - s_b)) * (np.pi - psi_at_b)
        else:
            raise ValueError(f"Unknown branch: {branch}")
            
        return psi

    def _param_to_psi(self, x, branch="sphere"):
        """
        Converts control parameter vector into a smooth psi(s) profile via cubic spline.
        """
        n_ctrl_Lo = 7
        n_ctrl_Ld = 9
        
        s_ctrl_Lo = np.linspace(0, self.s_split, n_ctrl_Lo)
        s_ctrl_Ld = np.linspace(self.s_split, self.S_tot, n_ctrl_Ld + 1)[1:]
        s_ctrl = np.concatenate([s_ctrl_Lo, s_ctrl_Ld])
        
        psi_ctrl = np.zeros(len(s_ctrl))
        psi_ctrl[0] = 0.0
        psi_ctrl[-1] = np.pi
        psi_ctrl[1:-1] = x
        
        cs = CubicSpline(s_ctrl, psi_ctrl, bc_type='natural')
        psi_grid = cs(self.s_grid)
        return psi_grid

    def optimize_shape(self, theta=0.0, branch="sphere", maxiter=600):
        """
        Optimizes the vesicle shape for a given coverage theta and branch.
        """
        psi_init = self.generate_initial_guess(branch)
        
        n_ctrl_Lo = 7
        n_ctrl_Ld = 9
        s_ctrl_Lo = np.linspace(0, self.s_split, n_ctrl_Lo)
        s_ctrl_Ld = np.linspace(self.s_split, self.S_tot, n_ctrl_Ld + 1)[1:]
        s_ctrl = np.concatenate([s_ctrl_Lo, s_ctrl_Ld])
        
        x0 = np.interp(s_ctrl[1:-1], self.s_grid, psi_init)
        
        # Constraint penalty multipliers
        lambda_A_Lo = 1.0e5
        lambda_A_tot = 1.0e5
        lambda_V = 1.0e5
        lambda_pole = 1.0e5
        
        def objective(x):
            psi_grid = self._param_to_psi(x, branch)
            res = compute_energy_and_constraints(
                self.s_grid, psi_grid, self.s_split, self.params, theta
            )
            
            penalty = (
                lambda_A_Lo * (res['res_A_Lo']**2) +
                lambda_A_tot * (res['res_A_tot']**2) +
                lambda_V * (res['res_V']**2) +
                lambda_pole * (res['res_r_pole']**2)
            )
            
            # Branch topological guidance
            # Dimple depth = z[0] - z[idx_split] (inward when z[0] > z[idx_split], i.e., apex pushed inside)
            r_prof, z_prof, _ = reconstruct_profile(self.s_grid, psi_grid, self.s_split)
            # apex is at s=0 (z[0]=0). Domain boundary is at z[idx_split].
            # For sphere: z_b > z_apex = 0.
            # For inward dimple: the domain curves inward such that relative to boundary, apex is indented inside.
            # Thus, mean curvature in Lo should be negative on average for inward dimple.
            c1, c2, H, _ = compute_curvatures(self.s_grid, psi_grid, r_prof)
            H_Lo_mean = np.mean(H[:self.idx_split])
            
            if branch == "inward":
                # Guide inward branch to explore H_Lo < 0 (invagination)
                if H_Lo_mean > 0.05:
                    penalty += 5.0e3 * ((H_Lo_mean - 0.05)**2)
            elif branch == "outward":
                # Guide outward branch to explore H_Lo > 0.3 (budding)
                if H_Lo_mean < 0.25:
                    penalty += 5.0e3 * ((H_Lo_mean - 0.25)**2)
            elif branch == "sphere":
                # Sphere: H_Lo ~ 1/R_ves ~ 0.2
                if abs(H_Lo_mean - (1.0 / self.params.R_ves)) > 0.15:
                    penalty += 2.0e3 * ((H_Lo_mean - (1.0 / self.params.R_ves))**2)
                    
            return res['total_energy'] + penalty

        res_opt = minimize(
            objective,
            x0,
            method='L-BFGS-B',
            options={'maxiter': maxiter, 'ftol': 1e-10, 'gtol': 1e-8}
        )
        
        psi_final = self._param_to_psi(res_opt.x, branch)
        eval_final = compute_energy_and_constraints(
            self.s_grid, psi_final, self.s_split, self.params, theta
        )
        
        eval_final['branch'] = branch
        eval_final['theta'] = theta
        eval_final['opt_success'] = res_opt.success
        eval_final['opt_nit'] = res_opt.nit
        eval_final['s_grid'] = self.s_grid
        eval_final['psi_grid'] = psi_final
        return eval_final
