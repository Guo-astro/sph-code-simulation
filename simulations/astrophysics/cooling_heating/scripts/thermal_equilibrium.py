#!/usr/bin/env python3
"""
Thermal and Chemical Equilibrium Solver
Based on Koyama & Inutsuka (2000, ApJ 533, 2000)

Solves for equilibrium temperature and chemical abundances
at a given density and radiation field strength.
"""

import numpy as np
from scipy.optimize import fsolve, brentq
from chemistry_network import ChemistryNetwork, k_B, m_H

class ThermalEquilibriumSolver:
    """
    Solve for thermal and chemical equilibrium state.
    """
    
    def __init__(self, G0=1.7, N_H_col=1e19):
        """
        Initialize solver.
        
        Parameters
        ----------
        G0 : float
            UV radiation field (Habing units)
        N_H_col : float
            Total H column density [cm^-2]
        """
        self.chem = ChemistryNetwork(G0=G0, N_H_col=N_H_col)
        self.G0 = G0
        self.N_H_col = N_H_col
        
    def pressure_equilibrium(self, n, T, x_e, x_2):
        """
        Calculate gas pressure.
        
        Parameters
        ----------
        n : float
            Number density [cm^-3]
        T : float
            Temperature [K]
        x_e : float
            Electron fraction
        x_2 : float
            H2 fraction
            
        Returns
        -------
        P : float
            Pressure [dyne cm^-2]
        """
        # Equation of state from paper
        # P = (1.1 + x_e - x_2/2) * n * k_B * T
        # He contribution: 0.1, electron: x_e, H2 molecular: -x_2/2
        
        P = (1.1 + x_e - 0.5 * x_2) * n * k_B * T
        
        return P
    
    def solve_chemistry_equilibrium(self, n, T, N_HI=0, N_H2=0):
        """
        Solve chemical equilibrium at fixed n, T.
        
        Parameters
        ----------
        n : float
            Number density [cm^-3]
        T : float
            Temperature [K]
        N_HI : float
            HI column density [cm^-2] (for UV attenuation)
        N_H2 : float
            H2 column density [cm^-2] (for self-shielding)
            
        Returns
        -------
        x_e, x_2, x_CO : floats
            Chemical fractions at equilibrium
        """
        # Initial guess
        if T > 5000:
            x_e_init = 0.1
            x_2_init = 1e-6
        else:
            x_e_init = 1e-4
            x_2_init = 0.5
        x_CO_init = 1e-6
        
        def residuals(X):
            x_e, x_2, x_CO = X
            
            # Bounds checking
            x_e = np.clip(x_e, 1e-10, 1.0)
            x_2 = np.clip(x_2, 1e-10, 2.0)
            x_CO = np.clip(x_CO, 1e-20, 0.001)
            
            x_H = 1.0 - x_e - x_2
            if x_H < 0:
                x_H = 1e-10
            
            # H ionization equilibrium
            zeta_ion = self.chem.h_ionization_rate(n, T, x_e, N_HI)
            alpha_rec = self.chem.h_recombination_rate(T)
            
            # dx_e/dt = 0: ionization = recombination
            res_e = zeta_ion * x_H - alpha_rec * x_e**2 * n
            
            # H2 formation/destruction equilibrium
            R_form = self.chem.h2_formation_rate(T)
            R_ad = self.chem.h2_associative_detachment_rate(T)
            k_D = self.chem.h2_collisional_dissociation_rate(T)
            R_pump = self.chem.h2_photodissociation_rate(N_H2)
            R_CR = self.chem.h2_cosmic_ray_dissociation_rate()
            
            # Formation rate (grain + associative detachment)
            formation = 0.5 * (R_form * x_H * n + R_ad * x_e * x_H * n**2)
            
            # Destruction rate
            n_H2 = 0.5 * x_2 * n
            destruction = n_H2 * (k_D * x_H * n + R_pump + R_CR)
            
            # dx_2/dt = 0
            res_2 = formation - destruction
            
            # CO chemistry (simplified)
            x_Cp = 3e-4  # Assume all C is C+ initially
            x_OI = 4.6e-4  # Neutral O
            
            rate_form = self.chem.co_formation_rate(n, x_Cp, x_OI)
            rate_dest = self.chem.co_photodissociation_rate()
            
            # dxCO/dt = 0
            res_CO = rate_form - x_CO * n * rate_dest
            
            return [res_e, res_2, res_CO]
        
        try:
            solution = fsolve(residuals, [x_e_init, x_2_init, x_CO_init], 
                              full_output=False, xtol=1e-6)
            x_e, x_2, x_CO = solution
            
            # Sanity checks
            x_e = np.clip(x_e, 1e-10, 1.0)
            x_2 = np.clip(x_2, 1e-10, 2.0)
            x_CO = np.clip(x_CO, 1e-20, 3e-4)
            
            return x_e, x_2, x_CO
            
        except Exception as e:
            print(f"Warning: Chemistry solver failed at n={n:.2e}, T={T:.1f}: {e}")
            # Return reasonable defaults
            if T > 5000:
                return 0.1, 1e-6, 1e-10
            else:
                return 1e-4, 0.5, 1e-6
    
    def solve_thermal_equilibrium(self, n, N_HI=0, N_H2=0, 
                                   T_min=10, T_max=2e4):
        """
        Solve for thermal equilibrium at fixed density.
        
        Find ALL temperature roots where net heating = 0.
        Thermal bistability means there can be multiple solutions!
        
        Parameters
        ----------
        n : float
            Number density [cm^-3]
        N_HI : float
            HI column density [cm^-2]
        N_H2 : float
            H2 column density [cm^-2]
        T_min, T_max : float
            Temperature search range [K]
            
        Returns
        -------
        T_eq : float or array
            Equilibrium temperature(s) [K]
        x_e, x_2, x_CO : floats or arrays
            Chemical fractions at equilibrium
        """
        
        def net_heating(T):
            """Net heating at given T (with chemistry equilibrium)."""
            # Solve chemistry at this T
            x_e, x_2, x_CO = self.solve_chemistry_equilibrium(n, T, N_HI, N_H2)
            
            # Calculate net heating
            net, _, _ = self.chem.net_heating_cooling(n, T, x_e, x_2, x_CO, 
                                                      N_HI, N_H2)
            
            return net
        
        # Scan temperature range to find all roots
        # Use fine grid to capture thermal bistability
        T_scan = np.logspace(np.log10(T_min), np.log10(T_max), 200)
        net_scan = np.array([net_heating(T) for T in T_scan])
        
        # Find all sign changes (roots)
        sign_changes = np.where(np.diff(np.sign(net_scan)))[0]
        
        roots = []
        for idx in sign_changes:
            try:
                # Refine root in this interval
                T_root = brentq(net_heating, T_scan[idx], T_scan[idx+1], xtol=1.0)
                roots.append(T_root)
            except:
                pass
        
        if len(roots) == 0:
            # No root found - use minimum
            idx_min = np.argmin(np.abs(net_scan))
            T_eq = T_scan[idx_min]
            x_e, x_2, x_CO = self.solve_chemistry_equilibrium(n, T_eq, N_HI, N_H2)
            return T_eq, x_e, x_2, x_CO
        
        elif len(roots) == 1:
            # Single solution
            T_eq = roots[0]
            x_e, x_2, x_CO = self.solve_chemistry_equilibrium(n, T_eq, N_HI, N_H2)
            return T_eq, x_e, x_2, x_CO
        
        else:
            # Multiple solutions - thermal bistability!
            # Return the stable branches (coldest and hottest)
            # Middle root is thermally unstable
            
            # For K&I 2000, we need to trace BOTH stable branches
            # At low n: return WNM (hot branch)
            # At high n: return CNM (cold branch)
            # Transition happens around n ~ 1 cm^-3
            
            if n < 1.0:
                # Low density - prefer hot branch (WNM)
                T_eq = max(roots)
            else:
                # High density - prefer cold branch (CNM)
                T_eq = min(roots)
            
            x_e, x_2, x_CO = self.solve_chemistry_equilibrium(n, T_eq, N_HI, N_H2)
            return T_eq, x_e, x_2, x_CO
    
    def compute_equilibrium_curve(self, n_array, N_HI_array=None, N_H2_array=None):
        """
        Compute equilibrium T, P, chemistry along density sequence.
        
        Parameters
        ----------
        n_array : array
            Array of densities [cm^-3]
        N_HI_array : array, optional
            HI column densities [cm^-2]
        N_H2_array : array, optional
            H2 column densities [cm^-2]
            
        Returns
        -------
        results : dict
            Dictionary containing equilibrium quantities
        """
        n_points = len(n_array)
        
        if N_HI_array is None:
            N_HI_array = np.zeros(n_points)
        if N_H2_array is None:
            N_H2_array = np.zeros(n_points)
        
        T_eq = np.zeros(n_points)
        P_eq = np.zeros(n_points)
        x_e_eq = np.zeros(n_points)
        x_2_eq = np.zeros(n_points)
        x_CO_eq = np.zeros(n_points)
        
        for i, n in enumerate(n_array):
            print(f"Solving equilibrium at n = {n:.2e} cm^-3...")
            
            N_HI = N_HI_array[i]
            N_H2 = N_H2_array[i]
            
            # Solve thermal equilibrium
            T, x_e, x_2, x_CO = self.solve_thermal_equilibrium(n, N_HI, N_H2)
            
            # Calculate pressure
            P = self.pressure_equilibrium(n, T, x_e, x_2)
            
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
            'N_HI': N_HI_array,
            'N_H2': N_H2_array
        }
        
        return results
    
    def compute_timescales(self, n, T, x_e, x_2, x_CO):
        """
        Compute various timescales.
        
        Returns
        -------
        timescales : dict
            Dictionary of timescales [years]
        """
        # Cooling time: t_cool = E / Lambda
        # E = 1.5 * n * k_B * T (thermal energy per unit volume)
        # Lambda_vol = n * Lambda (cooling rate per volume)
        
        _, _, cooling_dict = self.chem.net_heating_cooling(n, T, x_e, x_2, x_CO)
        
        Lambda = cooling_dict['total']
        
        if Lambda > 0:
            t_cool = 1.5 * k_B * T / Lambda
        else:
            t_cool = 1e20  # Very long
        
        # Recombination time
        alpha_rec = self.chem.h_recombination_rate(T)
        if x_e > 1e-10:
            t_rec = 1.0 / (alpha_rec * x_e * n)
        else:
            t_rec = 1e20
        
        # Free-fall time
        G_grav = 6.674e-8  # cm^3 g^-1 s^-2
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
        yr = 3.15576e7  # seconds per year
        
        timescales = {
            't_cool': t_cool / yr,
            't_rec': t_rec / yr,
            't_ff': t_ff / yr,
            't_h2': t_h2 / yr
        }
        
        return timescales
