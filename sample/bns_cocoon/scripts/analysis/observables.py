#!/usr/bin/env python3
"""
Observable Quantities from Cocoon Shock Breakout

Computes predicted observables from simulation results:
- Radiated energy E_rad
- Observed temperature T_obs
- Peak luminosity L_peak
- Flash duration t_flash
- Isotropic equivalent energy E_γ,iso

For comparison with GRB 170817A and predictions for future events.
"""

import argparse
import json
from pathlib import Path
from dataclasses import dataclass, asdict
from typing import Optional

import numpy as np


# Physical constants (CGS)
C_CGS = 2.998e10        # Speed of light [cm/s]
SIGMA_SB = 5.67e-5      # Stefan-Boltzmann constant [erg/cm²/s/K⁴]
A_RAD = 4 * SIGMA_SB / C_CGS  # Radiation constant [erg/cm³/K⁴]
K_B = 1.38e-16          # Boltzmann constant [erg/K]
H_PLANCK = 6.63e-27     # Planck constant [erg·s]
KEV_ERG = 1.6e-9        # 1 keV in erg


@dataclass
class BreakoutObservables:
    """Observable quantities from shock breakout."""
    # Breakout parameters
    r_breakout_cm: float      # Breakout radius [cm]
    v_breakout_c: float       # Breakout velocity [c]
    gamma_breakout: float     # Lorentz factor at breakout
    
    # Energetics
    e_internal_erg: float     # Internal energy at breakout [erg]
    f_rad: float              # Radiation efficiency
    e_rad_erg: float          # Radiated energy [erg]
    e_iso_erg: float          # Isotropic equivalent energy [erg]
    
    # Temperatures
    t_comoving_K: float       # Comoving temperature [K]
    t_comoving_keV: float     # Comoving temperature [keV]
    t_observed_K: float       # Observer-frame temperature [K]
    t_observed_keV: float     # Observer-frame temperature [keV]
    
    # Timescales
    t_light_crossing_s: float # Light crossing time [s]
    t_curvature_s: float      # Curvature timescale [s]
    t_flash_s: float          # Flash duration [s]
    
    # Luminosity
    l_peak_erg_s: float       # Peak luminosity [erg/s]
    l_peak_iso_erg_s: float   # Isotropic equivalent luminosity [erg/s]
    
    # Geometry
    solid_angle_sr: float     # Emission solid angle [sr]
    opening_angle_deg: float  # Opening angle [degrees]


class ObservablesCalculator:
    """
    Calculate observable quantities from breakout state.
    
    Based on Nakar & Sari (2012), Nakar & Piran (2017) formalism.
    """
    
    def __init__(self, 
                 f_rad: float = 0.1,
                 opening_angle_deg: float = 30.0,
                 redshift: float = 0.0,
                 length_unit_cm: float = 1e10,
                 density_unit_cgs: float = 1e-10):
        """
        Initialize calculator.
        
        Args:
            f_rad: Radiation efficiency (fraction of internal energy radiated)
            opening_angle_deg: Emission cone half-angle [degrees]
            redshift: Source redshift
            length_unit_cm: Code length unit [cm]
            density_unit_cgs: Code density unit [g/cm³]
        """
        self.f_rad = f_rad
        self.opening_angle_deg = opening_angle_deg
        self.opening_angle_rad = np.radians(opening_angle_deg)
        self.redshift = redshift
        self.length_unit = length_unit_cm
        self.density_unit = density_unit_cgs
        
        # Solid angle of emission cone
        self.solid_angle = 2 * np.pi * (1 - np.cos(self.opening_angle_rad))
    
    def compute(self, r_breakout: float, v_breakout: float,
                internal_energy: float, shell_volume: float,
                gamma_breakout: Optional[float] = None) -> BreakoutObservables:
        """
        Compute observables from breakout state.
        
        All inputs in code units unless specified.
        
        Args:
            r_breakout: Breakout radius [code units]
            v_breakout: Breakout velocity [c]
            internal_energy: Internal energy [code units]
            shell_volume: Shocked shell volume [code units³]
            gamma_breakout: Lorentz factor (computed from v if not provided)
            
        Returns:
            BreakoutObservables with all computed quantities
        """
        # Convert to CGS
        r_bo_cm = r_breakout * self.length_unit
        v_bo = min(v_breakout, 0.9999)  # Limit to < c
        
        if gamma_breakout is None:
            gamma_bo = 1.0 / np.sqrt(1 - v_bo**2)
        else:
            gamma_bo = gamma_breakout
        
        # Internal energy
        # Code energy = ρ V u, where u is specific internal energy
        # E_int = ∫ ρ u dV ≈ <ρ u> V
        e_int_code = internal_energy
        e_int_erg = e_int_code * self.density_unit * self.length_unit**3 * C_CGS**2
        
        # Radiated energy
        e_rad_erg = self.f_rad * e_int_erg
        
        # Isotropic equivalent energy
        e_iso_erg = (4 * np.pi / self.solid_angle) * e_rad_erg
        
        # Shell volume in CGS
        V_shell_cm3 = shell_volume * self.length_unit**3
        
        # Comoving radiation energy density
        e_rad_density = e_int_erg / V_shell_cm3  # Approximate
        
        # Comoving temperature from Stefan-Boltzmann
        # e = a T^4
        T_comoving = (e_rad_density / A_RAD)**0.25
        T_comoving_keV = T_comoving * K_B / KEV_ERG
        
        # Observed temperature (Doppler boosted for head-on)
        # For isotropic expansion: T_obs ≈ Γ T / (1+z)
        doppler_factor = gamma_bo  # Approximate for head-on
        T_observed = doppler_factor * T_comoving / (1 + self.redshift)
        T_observed_keV = T_observed * K_B / KEV_ERG
        
        # Timescales
        # Light crossing time
        t_lc = r_bo_cm / C_CGS
        
        # Curvature timescale (angular spreading)
        # t_curv ≈ R_bo / (2 c Γ²)
        t_curv = r_bo_cm / (2 * C_CGS * gamma_bo**2)
        
        # Shell thickness crossing time
        # ΔR ≈ R_bo / Γ for relativistic shocks
        delta_R = r_bo_cm / gamma_bo
        t_shell = delta_R / C_CGS
        
        # Flash duration is longer of curvature and shell crossing
        t_flash = max(t_curv, t_shell)
        
        # Peak luminosity
        L_peak = e_rad_erg / t_flash
        
        # Isotropic equivalent luminosity
        L_peak_iso = (4 * np.pi / self.solid_angle) * L_peak
        
        return BreakoutObservables(
            r_breakout_cm=r_bo_cm,
            v_breakout_c=v_bo,
            gamma_breakout=gamma_bo,
            e_internal_erg=e_int_erg,
            f_rad=self.f_rad,
            e_rad_erg=e_rad_erg,
            e_iso_erg=e_iso_erg,
            t_comoving_K=T_comoving,
            t_comoving_keV=T_comoving_keV,
            t_observed_K=T_observed,
            t_observed_keV=T_observed_keV,
            t_light_crossing_s=t_lc,
            t_curvature_s=t_curv,
            t_flash_s=t_flash,
            l_peak_erg_s=L_peak,
            l_peak_iso_erg_s=L_peak_iso,
            solid_angle_sr=self.solid_angle,
            opening_angle_deg=self.opening_angle_deg,
        )
    
    def print_summary(self, obs: BreakoutObservables):
        """Print summary of observables."""
        print("\n" + "="*60)
        print("SHOCK BREAKOUT OBSERVABLES")
        print("="*60)
        
        print(f"\nBreakout Geometry:")
        print(f"  R_bo     = {obs.r_breakout_cm:.2e} cm")
        print(f"  v_bo     = {obs.v_breakout_c:.3f} c")
        print(f"  Γ_bo     = {obs.gamma_breakout:.2f}")
        print(f"  θ_open   = {obs.opening_angle_deg:.1f}°")
        
        print(f"\nEnergetics:")
        print(f"  E_int    = {obs.e_internal_erg:.2e} erg")
        print(f"  E_rad    = {obs.e_rad_erg:.2e} erg")
        print(f"  E_iso    = {obs.e_iso_erg:.2e} erg")
        
        print(f"\nTemperature:")
        print(f"  T_co     = {obs.t_comoving_K:.2e} K ({obs.t_comoving_keV:.1f} keV)")
        print(f"  T_obs    = {obs.t_observed_K:.2e} K ({obs.t_observed_keV:.1f} keV)")
        
        print(f"\nTimescales:")
        print(f"  t_lc     = {obs.t_light_crossing_s:.2e} s")
        print(f"  t_curv   = {obs.t_curvature_s:.2e} s")
        print(f"  t_flash  = {obs.t_flash_s:.2e} s")
        
        print(f"\nLuminosity:")
        print(f"  L_peak   = {obs.l_peak_erg_s:.2e} erg/s")
        print(f"  L_iso    = {obs.l_peak_iso_erg_s:.2e} erg/s")
        
        print("="*60)


def compare_with_grb170817a(obs: BreakoutObservables):
    """Compare computed observables with GRB 170817A."""
    print("\n" + "="*60)
    print("COMPARISON WITH GRB 170817A")
    print("="*60)
    
    # GRB 170817A observed values
    grb170817a = {
        'E_iso': 4.6e46,       # erg (Goldstein+2017)
        'E_iso_err': 1.5e46,
        'T_obs_keV': 185,      # keV (peak energy)
        'T_obs_err': 62,
        't_duration': 2.0,     # s (T90)
        'L_iso': 2.3e47,       # erg/s (peak)
    }
    
    print("\n         | Simulation | GRB 170817A | Ratio")
    print("-" * 55)
    
    ratio_E = obs.e_iso_erg / grb170817a['E_iso']
    print(f"E_iso    | {obs.e_iso_erg:.2e} | {grb170817a['E_iso']:.2e} | {ratio_E:.2f}")
    
    ratio_T = obs.t_observed_keV / grb170817a['T_obs_keV']
    print(f"T_obs    | {obs.t_observed_keV:.1f} keV | {grb170817a['T_obs_keV']} keV | {ratio_T:.2f}")
    
    ratio_t = obs.t_flash_s / grb170817a['t_duration']
    print(f"t_flash  | {obs.t_flash_s:.2f} s | {grb170817a['t_duration']} s | {ratio_t:.2f}")
    
    ratio_L = obs.l_peak_iso_erg_s / grb170817a['L_iso']
    print(f"L_iso    | {obs.l_peak_iso_erg_s:.2e} | {grb170817a['L_iso']:.2e} | {ratio_L:.2f}")
    
    print("-" * 55)
    
    # Assessment
    print("\nAssessment:")
    if 0.1 < ratio_E < 10 and 0.1 < ratio_T < 10:
        print("✓ Results are broadly consistent with GRB 170817A")
    elif ratio_E < 0.01 or ratio_E > 100:
        print("✗ Energetics differ significantly from GRB 170817A")
    else:
        print("~ Results are within order of magnitude of GRB 170817A")


def main():
    parser = argparse.ArgumentParser(description='Compute breakout observables')
    parser.add_argument('breakout_file', nargs='?', default=None,
                       help='Breakout state JSON file')
    parser.add_argument('--r-breakout', type=float, default=1e11,
                       help='Breakout radius [cm]')
    parser.add_argument('--v-breakout', type=float, default=0.5,
                       help='Breakout velocity [c]')
    parser.add_argument('--e-internal', type=float, default=1e50,
                       help='Internal energy [erg]')
    parser.add_argument('--f-rad', type=float, default=0.1,
                       help='Radiation efficiency')
    parser.add_argument('--opening-angle', type=float, default=30.0,
                       help='Opening angle [degrees]')
    parser.add_argument('--redshift', type=float, default=0.0,
                       help='Source redshift')
    parser.add_argument('--compare-grb', action='store_true',
                       help='Compare with GRB 170817A')
    parser.add_argument('--output', '-o', default=None,
                       help='Output JSON file')
    
    args = parser.parse_args()
    
    calc = ObservablesCalculator(
        f_rad=args.f_rad,
        opening_angle_deg=args.opening_angle,
        redshift=args.redshift,
        length_unit_cm=1.0,  # Using CGS directly for command-line
        density_unit_cgs=1.0,
    )
    
    # Compute
    gamma_bo = 1.0 / np.sqrt(1 - args.v_breakout**2)
    V_shell = (4/3) * np.pi * args.r_breakout**3 / gamma_bo
    
    obs = calc.compute(
        r_breakout=args.r_breakout,
        v_breakout=args.v_breakout,
        internal_energy=args.e_internal / C_CGS**2,  # Convert to mass units
        shell_volume=V_shell,
        gamma_breakout=gamma_bo,
    )
    
    calc.print_summary(obs)
    
    if args.compare_grb:
        compare_with_grb170817a(obs)
    
    if args.output:
        with open(args.output, 'w') as f:
            json.dump(asdict(obs), f, indent=2)
        print(f"\nSaved to {args.output}")


if __name__ == '__main__':
    main()
