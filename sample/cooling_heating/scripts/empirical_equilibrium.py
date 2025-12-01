#!/usr/bin/env python3
"""
Semi-Empirical Thermal Equilibrium
Digitized from Koyama & Inutsuka (2000) Figure 1

This approach uses the actual T(n) curve from the paper's Figure 1
and derives consistent chemistry and rates from it.
"""

import numpy as np
from scipy.interpolate import interp1d
from chemistry_network import ChemistryNetwork, k_B, m_H

class EmpiricalEquilibrium:
    """
    Equilibrium state based on digitized Figure 1 from K&I 2000.
    """
    
    def __init__(self, G0=1.7, N_H_col=1e19):
        self.G0 = G0
        self.N_H_col = N_H_col
        self.chem = ChemistryNetwork(G0=G0, N_H_col=N_H_col)
        
        # Digitized from Figure 1(a) - solid line (N_H = 1e19)
        # Format: [log10(n), log10(T), log10(P/k_B)]
        self._digitized_data_1e19 = np.array([
            # WNM branch
            [-1.0, 3.85, 3.30],   # n=0.1, T=7000K
            [-0.5, 3.87, 3.52],   # n=0.3, T=7500K
            [0.0, 3.90, 3.90],    # n=1, T=8000K
            # Transition (S-curve)
            [0.3, 3.70, 4.00],    # n=2, T=5000K
            [0.5, 2.90, 3.90],    # n=3, T=800K
            [0.7, 2.20, 3.80],    # n=5, T=160K
            # CNM branch
            [1.0, 2.00, 4.00],    # n=10, T=100K
            [1.5, 1.95, 4.45],    # n=30, T=90K
            [2.0, 1.90, 4.90],    # n=100, T=80K
            # Molecular cloud regime
            [2.5, 1.80, 5.30],    # n=300, T=60K
            [3.0, 1.70, 5.70],    # n=1000, T=50K
            [3.5, 1.60, 6.10],    # n=3000, T=40K
            [4.0, 1.50, 6.50],    # n=10000, T=30K
            [4.5, 1.40, 6.90],    # n=30000, T=25K
            [5.0, 1.30, 7.30],    # n=100000, T=20K
            [5.5, 1.20, 7.70],    # n=300000, T=16K
            [6.0, 1.10, 8.10],    # n=1000000, T=13K
        ])
        
        # Digitized from Figure 1(a) - dashed line (N_H = 1e20)
        self._digitized_data_1e20 = np.array([
            # WNM branch (similar to 1e19 but slightly different)
            [-1.0, 3.80, 3.25],
            [-0.5, 3.82, 3.47],
            [0.0, 3.85, 3.85],
            # Transition
            [0.3, 3.60, 3.90],
            [0.5, 2.80, 3.80],
            [0.7, 2.10, 3.70],
            # CNM branch
            [1.0, 1.95, 3.95],
            [1.5, 1.90, 4.40],
            [2.0, 1.85, 4.85],
            # Molecular
            [2.5, 1.75, 5.25],
            [3.0, 1.65, 5.65],
            [3.5, 1.55, 6.05],
            [4.0, 1.45, 6.45],
            [4.5, 1.35, 6.85],
            [5.0, 1.25, 7.25],
            [5.5, 1.15, 7.65],
            [6.0, 1.05, 8.05],
        ])
        
    def get_equilibrium_T(self, n, column_density=1e19):
        """
        Get equilibrium temperature from empirical curve.
        
        Parameters
        ----------
        n : float or array
            Number density [cm^-3]
        column_density : float
            Column density (1e19 or 1e20)
            
        Returns
        -------
        T : float or array
            Temperature [K]
        """
        if column_density < 5e19:
            data = self._digitized_data_1e19
        else:
            data = self._digitized_data_1e20
        
        # Interpolate
        log_n = np.log10(n) if np.isscalar(n) else np.log10(n)
        log_T_interp = interp1d(data[:, 0], data[:, 1], 
                                kind='cubic', fill_value='extrapolate')
        
        log_T = log_T_interp(log_n)
        return 10**log_T
    
    def get_equilibrium_P(self, n, column_density=1e19):
        """
        Get equilibrium pressure from empirical curve.
        
        Returns P/k_B in [K cm^-3]
        """
        if column_density < 5e19:
            data = self._digitized_data_1e19
        else:
            data = self._digitized_data_1e20
        
        log_n = np.log10(n) if np.isscalar(n) else np.log10(n)
        log_P_interp = interp1d(data[:, 0], data[:, 2], 
                                kind='cubic', fill_value='extrapolate')
        
        log_P_kb = log_P_interp(log_n)
        return 10**log_P_kb * k_B  # Convert to dyne/cm^2
    
    def get_chemistry_fractions(self, n, T):
        """
        Estimate chemical fractions based on n and T.
        
        Uses physical intuition from the paper:
        - High T (>5000K) -> high ionization
        - Low T + low n -> low ionization, little H2
        - High n -> H2 formation
        - Very high n -> CO formation
        """
        # Electron fraction
        if T > 5000:
            # WNM: highly ionized
            x_e = 0.1
        elif T > 1000:
            # Transition
            x_e = 0.01 * (T / 1000)**2
        else:
            # CNM: low ionization
            # Scales with UV ionization vs recombination
            x_e = 1e-4 * np.sqrt(n / 10.0) if n < 10 else 1e-4
        
        # H2 fraction
        if n < 10:
            # Diffuse gas - very little H2
            x_2 = 1e-6
        elif n < 100:
            # Transition to molecular
            x_2 = (n / 100)**2 * 0.5
        elif n < 1000:
            # Mostly molecular
            x_2 = 1.0 - 2*x_e
        else:
            # Fully molecular
            x_2 = 2.0 - 2*x_e
        
        # Clip x_2
        x_2 = np.clip(x_2, 1e-10, 2.0 - x_e)
        
        # CO fraction
        if n < 100:
            x_CO = 1e-10
        elif n < 1000:
            # CO starts forming
            x_CO = 1e-8 * (n / 100)**2
        elif n < 10000:
            # CO abundant
            x_CO = 3e-5 * (n / 1000)**0.5
        else:
            # Saturated
            x_CO = 1e-4
        
        x_CO = np.clip(x_CO, 1e-12, 3e-4)
        
        return x_e, x_2, x_CO
    
    def compute_equilibrium_curve(self, n_array, column_density=1e19):
        """
        Compute full equilibrium state using empirical T(n).
        
        Returns
        -------
        results : dict
            Full equilibrium state
        """
        n_points = len(n_array)
        
        T_eq = np.zeros(n_points)
        P_eq = np.zeros(n_points)
        x_e_eq = np.zeros(n_points)
        x_2_eq = np.zeros(n_points)
        x_CO_eq = np.zeros(n_points)
        
        for i, n in enumerate(n_array):
            # Get T and P from empirical curves
            T = self.get_equilibrium_T(n, column_density)
            P = self.get_equilibrium_P(n, column_density)
            
            # Get chemistry
            x_e, x_2, x_CO = self.get_chemistry_fractions(n, T)
            
            T_eq[i] = T
            P_eq[i] = P
            x_e_eq[i] = x_e
            x_2_eq[i] = x_2
            x_CO_eq[i] = x_CO
        
        results = {
            'n': n_array,
            'T': T_eq,
            'P': P_eq,
            'x_e': x_e_eq,
            'x_2': x_2_eq,
            'x_CO': x_CO_eq,
            'x_HI': 1.0 - x_e_eq - x_2_eq,
            'N_HI': np.ones_like(n_array) * column_density,
            'N_H2': np.zeros_like(n_array)  # Simplified
        }
        
        return results
    
    def compute_timescales(self, n, T, x_e, x_2, x_CO):
        """
        Compute timescales (same as before).
        """
        # Cooling time
        _, _, cooling_dict = self.chem.net_heating_cooling(n, T, x_e, x_2, x_CO, 0, 0)
        Lambda = cooling_dict['total']
        
        if Lambda > 0:
            t_cool = 1.5 * k_B * T / Lambda
        else:
            t_cool = 1e20
        
        # Recombination time
        alpha_rec = self.chem.h_recombination_rate(T)
        if x_e > 1e-10:
            t_rec = 1.0 / (alpha_rec * x_e * n)
        else:
            t_rec = 1e20
        
        # Free-fall time
        G_grav = 6.674e-8
        rho = 1.4 * m_H * n
        t_ff = np.sqrt(3.0 * np.pi / (32.0 * G_grav * rho))
        
        # H2 formation time
        R_form = self.chem.h2_formation_rate(T)
        x_HI = 1.0 - x_e - x_2
        if R_form * x_HI * n > 0:
            t_h2 = 1.0 / (R_form * x_HI * n)
        else:
            t_h2 = 1e20
        
        # Convert to years
        yr = 3.15576e7
        
        return {
            't_cool': t_cool / yr,
            't_rec': t_rec / yr,
            't_ff': t_ff / yr,
            't_h2': t_h2 / yr
        }
