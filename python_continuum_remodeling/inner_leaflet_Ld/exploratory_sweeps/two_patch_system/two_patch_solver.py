"""
Two-Patch Antipodal Vesicle Solver (North & South Lo domains, Equatorial Ld belt).
Minimizes energy for symmetric biconcave discocyte remodeling across domain sizes (D_Lo = 2um to 7um).
"""

import numpy as np
from scipy.optimize import minimize
from scipy.interpolate import CubicSpline
from sim.two_patch_system.two_patch_mechanics import (
    TwoPatchMembraneParameters,
    reconstruct_profile_two_patch,
    compute_curvatures_two_patch,
    compute_two_patch_energy
)

class TwoPatchVesicleSolver:
    def __init__(self, D_Lo=3.0, R_ves=5.0, n_points=140):
        self.D_Lo = D_Lo
        self.R_ves = R_ves
        self.n_points = n_points
        
        self.params = TwoPatchMembraneParameters(
            D_Lo=D_Lo,
            R_ves=R_ves,
            kappa_Lo=80.0,
            kappa_Ld_0=20.0,
            kappa_prot_sat=70.0
        )
        
        self.S_tot = np.pi * self.R_ves
        
        # Arc length boundaries:
        # North patch: s in [0, s_b1]
        # South patch: s in [s_b2, S_tot]
        # Area fraction per patch: f_single = (D_Lo / (2*R_ves))^2 * 0.25 (or exact cap)
        f_single = min(0.24, self.params.A_Lo_single / self.params.A_tot_target)
        theta_split1 = np.arccos(1.0 - 2.0 * f_single)
        self.s_b1 = self.R_ves * theta_split1
        self.s_b2 = self.S_tot - self.s_b1
        
        # Discretize: allocate points to North Lo, Equatorial Ld, South Lo
        frac_Lo1 = max(0.18, min(0.35, self.s_b1 / self.S_tot))
        n_Lo1 = int(self.n_points * frac_Lo1)
        n_Lo2 = n_Lo1
        n_Ld = self.n_points - 2 * n_Lo1
        
        s_Lo1 = np.linspace(0, self.s_b1, n_Lo1)
        s_Ld = np.linspace(self.s_b1, self.s_b2, n_Ld + 2)[1:-1]
        s_Lo2 = np.linspace(self.s_b2, self.S_tot, n_Lo2)
        self.s_grid = np.concatenate([s_Lo1, s_Ld, s_Lo2])
        self.idx_b1 = len(s_Lo1) - 1
        self.idx_b2 = len(s_Lo1) + len(s_Ld) - 1
        self.n_Lo1 = n_Lo1
        self.n_Lo2 = n_Lo2
        
    def generate_initial_guess(self, branch="sphere"):
        s = self.s_grid
        s_b1 = self.s_b1
        s_b2 = self.s_b2
        S_tot = self.S_tot
        
        psi_sphere = (s / S_tot) * np.pi
        
        if branch == "sphere":
            return psi_sphere
        elif branch == "inward":
            # Dual inward dimple (biconcave discocyte)
            psi = np.copy(psi_sphere)
            mask_Lo1 = s <= s_b1
            mask_Lo2 = s >= s_b2
            
            # North cap inverts (negative slope)
            psi[mask_Lo1] = - 1.0 * np.sin(np.pi * s[mask_Lo1] / s_b1) + (s[mask_Lo1] / s_b1) * psi_sphere[mask_Lo1][-1]
            
            # South cap inverts (steeper slope flattening inward at south pole)
            s_rel_south = (s[mask_Lo2] - s_b2) / (S_tot - s_b2)
            psi[mask_Lo2] = psi_sphere[mask_Lo2][0] + s_rel_south * (np.pi - psi_sphere[mask_Lo2][0]) + 1.0 * np.sin(np.pi * s_rel_south)
            return psi
        elif branch == "outward":
            psi = np.copy(psi_sphere)
            mask_Lo1 = s <= s_b1
            mask_Lo2 = s >= s_b2
            psi[mask_Lo1] = 1.2 * np.sin(np.pi * s[mask_Lo1] / s_b1) + (s[mask_Lo1] / s_b1) * psi_sphere[mask_Lo1][-1]
            s_rel_south = (s[mask_Lo2] - s_b2) / (S_tot - s_b2)
            psi[mask_Lo2] = psi_sphere[mask_Lo2][0] + s_rel_south * (np.pi - psi_sphere[mask_Lo2][0]) - 1.0 * np.sin(np.pi * s_rel_south)
            return psi
        else:
            raise ValueError(f"Unknown branch: {branch}")

    def _param_to_psi(self, x):
        n_ctrl_Lo1 = 7
        n_ctrl_Ld = 8
        n_ctrl_Lo2 = 7
        
        s_ctrl_Lo1 = np.linspace(0, self.s_b1, n_ctrl_Lo1)
        s_ctrl_Ld = np.linspace(self.s_b1, self.s_b2, n_ctrl_Ld + 2)[1:-1]
        s_ctrl_Lo2 = np.linspace(self.s_b2, self.S_tot, n_ctrl_Lo2)
        s_ctrl = np.concatenate([s_ctrl_Lo1, s_ctrl_Ld, s_ctrl_Lo2])
        
        psi_ctrl = np.zeros(len(s_ctrl))
        psi_ctrl[0] = 0.0
        psi_ctrl[-1] = np.pi
        psi_ctrl[1:-1] = x
        
        cs = CubicSpline(s_ctrl, psi_ctrl, bc_type='natural')
        psi_grid = cs(self.s_grid)
        return psi_grid

    def optimize_shape(self, theta=0.0, branch="sphere", maxiter=600):
        psi_init = self.generate_initial_guess(branch)
        
        n_ctrl_Lo1 = 7
        n_ctrl_Ld = 8
        n_ctrl_Lo2 = 7
        s_ctrl_Lo1 = np.linspace(0, self.s_b1, n_ctrl_Lo1)
        s_ctrl_Ld = np.linspace(self.s_b1, self.s_b2, n_ctrl_Ld + 2)[1:-1]
        s_ctrl_Lo2 = np.linspace(self.s_b2, self.S_tot, n_ctrl_Lo2)
        s_ctrl = np.concatenate([s_ctrl_Lo1, s_ctrl_Ld, s_ctrl_Lo2])
        
        x0 = np.interp(s_ctrl[1:-1], self.s_grid, psi_init)
        bounds = [(-np.pi/2, 3*np.pi/2) for _ in x0]
        
        lambda_A_Lo = 2.0e4
        lambda_A_tot = 2.0e4
        lambda_V = 2.0e4
        lambda_pole = 2.0e4
        
        def objective(x):
            psi_grid = self._param_to_psi(x)
            res = compute_two_patch_energy(
                self.s_grid, psi_grid, self.s_b1, self.s_b2, self.params, theta
            )
            
            penalty = (
                lambda_A_Lo * (res['res_A_Lo1']**2) +
                lambda_A_Lo * (res['res_A_Lo2']**2) +
                lambda_A_tot * (res['res_A_tot']**2) +
                lambda_V * (res['res_V']**2) +
                lambda_pole * (res['res_r_pole']**2)
            )
            
            min_r = np.min(res['r_profile'][1:-1])
            if min_r < 0:
                penalty += 1.0e6 * (min_r**2)
                
            r_prof, _, _, _ = reconstruct_profile_two_patch(self.s_grid, psi_grid, self.s_b1, self.s_b2)
            _, _, H, _ = compute_curvatures_two_patch(self.s_grid, psi_grid, r_prof)
            H_Lo1_mean = np.mean(H[:self.idx_b1])
            H_Lo2_mean = np.mean(H[self.idx_b2:])
            
            if branch == "inward":
                if H_Lo1_mean > 0.05:
                    penalty += 3.0e3 * ((H_Lo1_mean - 0.05)**2)
                if H_Lo2_mean > 0.05:
                    penalty += 3.0e3 * ((H_Lo2_mean - 0.05)**2)
            elif branch == "sphere":
                if abs(H_Lo1_mean - (1.0 / self.params.R_ves)) > 0.15:
                    penalty += 2.0e3 * ((H_Lo1_mean - (1.0 / self.params.R_ves))**2)
                    
            return res['total_energy'] + penalty

        res_opt = minimize(
            objective,
            x0,
            bounds=bounds,
            method='L-BFGS-B',
            options={'maxiter': maxiter, 'ftol': 1e-9, 'gtol': 1e-7}
        )
        
        psi_final = self._param_to_psi(res_opt.x)
        eval_final = compute_two_patch_energy(
            self.s_grid, psi_final, self.s_b1, self.s_b2, self.params, theta
        )
        
        eval_final['D_Lo'] = self.D_Lo
        eval_final['branch'] = branch
        eval_final['theta'] = theta
        eval_final['n_Lo1'] = self.n_Lo1
        eval_final['n_Lo2'] = self.n_Lo2
        eval_final['opt_success'] = res_opt.success
        eval_final['s_grid'] = self.s_grid
        eval_final['psi_grid'] = psi_final
        
        # Depths of North and South dimples (depth into vesicle)
        z_north_apex = eval_final['z_profile'][0]
        z_north_rim = eval_final['z_profile'][self.idx_b1]
        eval_final['dimple_depth_north'] = z_north_rim - z_north_apex
        
        z_south_apex = eval_final['z_profile'][-1]
        z_south_rim = eval_final['z_profile'][self.idx_b2]
        eval_final['dimple_depth_south'] = z_south_apex - z_south_rim
        
        eval_final['dimple_depth_mean'] = 0.5 * (abs(eval_final['dimple_depth_north']) + abs(eval_final['dimple_depth_south']))
        
        return eval_final
