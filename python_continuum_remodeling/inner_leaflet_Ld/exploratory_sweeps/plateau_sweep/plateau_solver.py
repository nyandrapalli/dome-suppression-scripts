"""
Plateau Solver for Two-Phase GUVs across Large Domain Fractions (D_Lo = 8.0 um to 14.0 um).
Models deep stomatocyte invagination, neck pinching, and internal sphere-in-sphere plateau limits.
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

class PlateauVesicleSolver:
    def __init__(self, D_Lo=8.0, R_ves=5.0, n_points=150):
        self.D_Lo = D_Lo
        self.R_ves = R_ves
        self.n_points = n_points
        
        # Calculate domain area fraction: f_Lo = (D_Lo / (2*R_ves))^2 * 0.25 (flat disc approximation)
        # Or spherical cap area: A_Lo = pi * (D_Lo/2)^2
        self.params = MembraneParameters(
            D_Lo=D_Lo,
            R_ves=R_ves,
            kappa_Lo=80.0,
            kappa_Ld_0=20.0,
            kappa_prot_sat=70.0
        )
        
        # S_tot ~ pi * R_ves
        self.S_tot = np.pi * self.R_ves
        
        # Target area fraction
        f_Lo = min(0.50, self.params.A_Lo_target / self.params.A_tot_target)
        theta_split = np.arccos(1.0 - 2.0 * f_Lo)
        self.s_split = self.R_ves * theta_split
        
        # High resolution in Lo domain: 50-80 points
        frac_Lo = min(0.60, max(0.40, self.s_split / self.S_tot))
        n_Lo = int(self.n_points * frac_Lo)
        n_Ld = self.n_points - n_Lo
        
        s_Lo = np.linspace(0, self.s_split, n_Lo)
        s_Ld = np.linspace(self.s_split, self.S_tot, n_Ld + 1)[1:]
        self.s_grid = np.concatenate([s_Lo, s_Ld])
        self.idx_split = len(s_Lo) - 1
        self.n_Lo = n_Lo
        
    def generate_initial_guess(self, branch="sphere"):
        s = self.s_grid
        s_b = self.s_split
        mask_Lo = s <= s_b
        mask_Ld = s > s_b
        
        psi_sphere = (s / self.S_tot) * np.pi
        
        if branch == "sphere":
            return psi_sphere
        elif branch == "inward":
            psi = np.copy(psi_sphere)
            # Invert Lo cap to negative tangent slope for deep invagination
            psi[mask_Lo] = - 1.2 * np.sin(np.pi * s[mask_Lo] / s_b) + (s[mask_Lo] / s_b) * psi_sphere[mask_Lo][-1]
            return psi
        elif branch == "outward":
            psi = np.copy(psi_sphere)
            psi[mask_Lo] = 1.3 * np.sin(np.pi * s[mask_Lo] / s_b) + (s[mask_Lo] / s_b) * psi_sphere[mask_Lo][-1]
            return psi
        else:
            raise ValueError(f"Unknown branch: {branch}")

    def _param_to_psi(self, x):
        n_ctrl_Lo = 10
        n_ctrl_Ld = 10
        
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
        psi_init = self.generate_initial_guess(branch)
        
        n_ctrl_Lo = 10
        n_ctrl_Ld = 10
        s_ctrl_Lo = np.linspace(0, self.s_split, n_ctrl_Lo)
        s_ctrl_Ld = np.linspace(self.s_split, self.S_tot, n_ctrl_Ld + 1)[1:]
        s_ctrl = np.concatenate([s_ctrl_Lo, s_ctrl_Ld])
        
        x0 = np.interp(s_ctrl[1:-1], self.s_grid, psi_init)
        bounds = [(-np.pi/2, 3*np.pi/2) for _ in x0]
        
        lambda_A_Lo = 2.5e4
        lambda_A_tot = 2.5e4
        lambda_V = 2.5e4
        lambda_pole = 2.5e4
        
        def objective(x):
            psi_grid = self._param_to_psi(x)
            res = compute_energy_and_constraints(
                self.s_grid, psi_grid, self.s_split, self.params, theta
            )
            
            penalty = (
                lambda_A_Lo * (res['res_A_Lo']**2) +
                lambda_A_tot * (res['res_A_tot']**2) +
                lambda_V * (res['res_V']**2) +
                lambda_pole * (res['res_r_pole']**2)
            )
            
            min_r = np.min(res['r_profile'][1:-1])
            if min_r < 0:
                penalty += 1.0e6 * (min_r**2)
                
            r_prof, _, _ = reconstruct_profile(self.s_grid, psi_grid, self.s_split)
            _, _, H, _ = compute_curvatures(self.s_grid, psi_grid, r_prof)
            H_Lo_mean = np.mean(H[:self.idx_split])
            
            if branch == "inward":
                if H_Lo_mean > 0.05:
                    penalty += 3.0e3 * ((H_Lo_mean - 0.05)**2)
            elif branch == "outward":
                if H_Lo_mean < 0.25:
                    penalty += 3.0e3 * ((H_Lo_mean - 0.25)**2)
            elif branch == "sphere":
                if abs(H_Lo_mean - (1.0 / self.params.R_ves)) > 0.15:
                    penalty += 2.0e3 * ((H_Lo_mean - (1.0 / self.params.R_ves))**2)
                    
            return res['total_energy'] + penalty

        res_opt = minimize(
            objective,
            x0,
            bounds=bounds,
            method='L-BFGS-B',
            options={'maxiter': maxiter, 'ftol': 1e-9, 'gtol': 1e-7}
        )
        
        psi_final = self._param_to_psi(res_opt.x)
        eval_final = compute_energy_and_constraints(
            self.s_grid, psi_final, self.s_split, self.params, theta
        )
        
        eval_final['D_Lo'] = self.D_Lo
        eval_final['branch'] = branch
        eval_final['theta'] = theta
        eval_final['n_Lo'] = self.n_Lo
        eval_final['opt_success'] = res_opt.success
        eval_final['s_grid'] = self.s_grid
        eval_final['psi_grid'] = psi_final
        
        z_apex = eval_final['z_profile'][0]
        z_rim = eval_final['z_profile'][self.idx_split]
        eval_final['measured_dimple_depth'] = z_rim - z_apex
        
        return eval_final
