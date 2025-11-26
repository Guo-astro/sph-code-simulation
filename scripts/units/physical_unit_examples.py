#!/usr/bin/env python3
"""
Physical Unit Conversion Examples for SR-GSPH

This script demonstrates how to convert between code units (c=1) and
physical units (CGS) for real astrophysical problems.

Usage:
    python physical_unit_examples.py
    
Examples shown:
    1. Neutron star merger: NS-NS inspiral and merger
    2. GRB afterglow: External shock deceleration
    3. Relativistic jet: AGN/microquasar jet dynamics
"""

import sys
from pathlib import Path

# Add scripts directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from units.unit_system import (
    RelativisticUnits, 
    SPEED_OF_LIGHT_CGS,
    GRAVITATIONAL_CONSTANT_CGS,
    MSUN_TO_G,
    KM_TO_CM,
    PC_TO_CM,
    PROTON_MASS_G
)
import numpy as np


def example_neutron_star_merger():
    """
    Example: Binary Neutron Star Merger
    
    Physical setup:
    - Two 1.4 M_sun neutron stars
    - Radius ~ 12 km each
    - Initial separation ~ 50 km
    - Nuclear density ~ 2.8 × 10^14 g/cm³
    - Merger time ~ 10 ms after start
    """
    print("\n" + "="*70)
    print("EXAMPLE 1: Binary Neutron Star Merger")
    print("="*70)
    
    # Create unit system for NS merger
    # Code unit L = 10 km, ρ = 10^14 g/cm³
    units = RelativisticUnits.create_neutron_star(length_km=10.0, density_scale=1e14)
    
    print(f"\nUnit System:")
    print(f"  Type: {units.get_type_name()}")
    print(f"  Speed of light in code: c = {units.c_code}")
    print(f"  Length label: {units.length_label}")
    print(f"  Time label: {units.time_label}")
    
    print(f"\nConversion Factors:")
    print(f"  1 code_length = {units.length_to_cgs:.3e} cm = {units.length_to_cgs/KM_TO_CM:.1f} km")
    print(f"  1 code_time = {units.time_to_cgs:.3e} s = {units.time_to_cgs*1e3:.4f} ms")
    print(f"  1 code_density = {units.density_to_cgs:.3e} g/cm³")
    print(f"  1 code_velocity = c = {units.velocity_to_cgs:.3e} cm/s")
    
    # Convert physical quantities to code units
    print(f"\nPhysical to Code Unit Conversion:")
    
    # NS radius
    R_ns_km = 12.0
    R_ns_code = units.from_physical_length(R_ns_km * KM_TO_CM)
    print(f"  NS radius: {R_ns_km} km → {R_ns_code:.2f} code_length")
    
    # Initial separation
    d_init_km = 50.0
    d_init_code = units.from_physical_length(d_init_km * KM_TO_CM)
    print(f"  Initial separation: {d_init_km} km → {d_init_code:.1f} code_length")
    
    # Nuclear density
    rho_nuc = 2.8e14  # g/cm³
    rho_nuc_code = units.from_physical_density(rho_nuc)
    print(f"  Nuclear density: {rho_nuc:.2e} g/cm³ → {rho_nuc_code:.2f} code_density")
    
    # Orbital velocity (approximate: v ~ sqrt(GM/r))
    M_ns = 1.4 * MSUN_TO_G
    v_orb = np.sqrt(GRAVITATIONAL_CONSTANT_CGS * M_ns / (d_init_km * KM_TO_CM))
    v_orb_c = v_orb / SPEED_OF_LIGHT_CGS
    print(f"  Orbital velocity: {v_orb:.3e} cm/s = {v_orb_c:.3f} c")
    
    # Merger timescale 
    print(f"\nSimulation Time Estimates:")
    t_merger_ms = 10.0
    t_merger_code = t_merger_ms * 1e-3 / units.time_to_cgs
    print(f"  Merger time: {t_merger_ms} ms → {t_merger_code:.0f} code_time")
    
    # Output interval
    dt_out_ms = 0.1
    dt_out_code = dt_out_ms * 1e-3 / units.time_to_cgs
    print(f"  Output interval: {dt_out_ms} ms → {dt_out_code:.1f} code_time")


def example_grb_afterglow():
    """
    Example: GRB Afterglow / External Shock
    
    Physical setup:
    - Isotropic energy E_iso ~ 10^52 erg
    - Initial Lorentz factor Γ₀ ~ 300
    - Circumburst density n₀ ~ 1 cm⁻³
    - Deceleration radius ~ 10^17 cm
    - Afterglow duration ~ days to weeks
    """
    print("\n" + "="*70)
    print("EXAMPLE 2: GRB Afterglow / External Shock")
    print("="*70)
    
    # Physical parameters
    E_iso = 1e52  # erg
    Gamma_0 = 300
    n_0 = 1.0  # cm^-3 (ambient density)
    rho_0 = n_0 * PROTON_MASS_G  # g/cm³
    
    # Deceleration radius: R_dec ~ (3E / 4π n m_p c²)^(1/3) / Gamma^(2/3)
    R_dec = (3 * E_iso / (4 * np.pi * n_0 * PROTON_MASS_G * SPEED_OF_LIGHT_CGS**2))**(1/3) / Gamma_0**(2/3)
    
    # Create unit system with L ~ R_dec
    L_scale = 1e17  # cm (order of deceleration radius)
    units = RelativisticUnits.create_relativistic(
        code_length_cm=L_scale,
        code_density_g_cm3=rho_0,
        length_label="10^17 cm",
        density_label="n_ISM"
    )
    
    print(f"\nPhysical Parameters:")
    print(f"  E_iso = {E_iso:.1e} erg")
    print(f"  Γ₀ = {Gamma_0}")
    print(f"  n₀ = {n_0} cm⁻³ = {rho_0:.3e} g/cm³")
    print(f"  R_dec ≈ {R_dec:.2e} cm")
    
    print(f"\nUnit System:")
    print(f"  1 code_length = {units.length_to_cgs:.3e} cm")
    print(f"  1 code_time = {units.time_to_cgs:.3e} s = {units.time_to_cgs/86400:.1f} days")
    print(f"  1 code_density = {units.density_to_cgs:.3e} g/cm³")
    
    print(f"\nCode Unit Conversions:")
    print(f"  R_dec in code: {R_dec/L_scale:.2f} code_length")
    
    # Deceleration time: t_dec ~ R_dec / (c * Gamma_0²)
    t_dec = R_dec / (SPEED_OF_LIGHT_CGS * Gamma_0**2)
    t_dec_code = t_dec / units.time_to_cgs
    print(f"  t_dec ≈ {t_dec:.2e} s = {t_dec/86400:.2f} days → {t_dec_code:.2f} code_time")
    
    # Lorentz factor (stored as velocity ~ sqrt(1 - 1/Γ²) )
    v_code = np.sqrt(1 - 1/Gamma_0**2)
    print(f"  Initial velocity: v = {v_code:.6f} c (Γ = {1/np.sqrt(1-v_code**2):.0f})")
    
    # Energy density
    e_internal = E_iso / (4/3 * np.pi * R_dec**3)  # erg/cm³
    P_internal = e_internal / 3  # For relativistic gas
    P_code = P_internal / units.pressure_to_cgs
    print(f"  Internal pressure scale: {P_internal:.2e} dyn/cm² → {P_code:.2f} code_pressure")


def example_relativistic_jet():
    """
    Example: Relativistic Jet (AGN / Microquasar)
    
    Physical setup:
    - Jet Lorentz factor Γ ~ 5-10
    - Jet radius ~ 0.01 pc
    - Length scale ~ 1 pc
    - Jet density contrast ~ 0.1 (jet/ambient)
    - Ambient ISM: n ~ 1 cm⁻³
    """
    print("\n" + "="*70)
    print("EXAMPLE 3: Relativistic Jet")
    print("="*70)
    
    # Physical parameters
    Gamma_jet = 7
    r_jet_pc = 0.01  # jet radius
    L_jet_pc = 1.0  # propagation length
    n_ism = 1.0  # cm^-3
    eta = 0.1  # density contrast (jet/ISM)
    
    rho_ism = n_ism * PROTON_MASS_G
    rho_jet = eta * rho_ism
    
    # Create unit system
    units = RelativisticUnits.create_relativistic_jet(
        length_pc=1.0,
        density_scale=rho_ism
    )
    
    print(f"\nPhysical Parameters:")
    print(f"  Γ_jet = {Gamma_jet}")
    print(f"  r_jet = {r_jet_pc} pc = {r_jet_pc * PC_TO_CM:.2e} cm")
    print(f"  L_jet = {L_jet_pc} pc")
    print(f"  n_ISM = {n_ism} cm⁻³")
    print(f"  η (jet/ISM density) = {eta}")
    
    print(f"\nUnit System:")
    print(f"  1 code_length = {units.length_to_cgs:.3e} cm = 1 pc")
    print(f"  1 code_time = {units.time_to_cgs:.3e} s = {units.time_to_cgs/(365.25*86400):.2f} years")
    print(f"  1 code_density = {units.density_to_cgs:.3e} g/cm³")
    
    print(f"\nCode Unit Conversions:")
    r_jet_code = r_jet_pc  # Since L_scale = 1 pc
    L_jet_code = L_jet_pc
    print(f"  Jet radius: {r_jet_code:.3f} code_length")
    print(f"  Jet length: {L_jet_code:.1f} code_length")
    
    rho_jet_code = units.from_physical_density(rho_jet)
    rho_ism_code = units.from_physical_density(rho_ism)
    print(f"  Jet density: {rho_jet_code:.2f} code_density")
    print(f"  ISM density: {rho_ism_code:.2f} code_density")
    
    # Jet velocity
    v_jet = np.sqrt(1 - 1/Gamma_jet**2)
    print(f"  Jet velocity: {v_jet:.4f} c")
    
    # Jet propagation time
    t_cross = L_jet_pc * PC_TO_CM / (v_jet * SPEED_OF_LIGHT_CGS)
    t_cross_years = t_cross / (365.25 * 86400)
    t_cross_code = t_cross / units.time_to_cgs
    print(f"  Crossing time: {t_cross_years:.2f} years → {t_cross_code:.2f} code_time")
    
    # Jet power estimate (kinetic)
    # L_jet ~ π r² ρ v³ Γ²
    L_jet_power = np.pi * (r_jet_pc * PC_TO_CM)**2 * rho_jet * (v_jet * SPEED_OF_LIGHT_CGS)**3 * Gamma_jet**2
    print(f"  Jet kinetic power: {L_jet_power:.2e} erg/s = {L_jet_power/1e44:.2f} × 10^44 erg/s")


def example_conversion_table():
    """Print a conversion table for quick reference"""
    print("\n" + "="*70)
    print("QUICK REFERENCE: Unit Conversion Table")
    print("="*70)
    
    print(f"\nPhysical Constants (CGS):")
    print(f"  c = {SPEED_OF_LIGHT_CGS:.6e} cm/s")
    print(f"  G = {GRAVITATIONAL_CONSTANT_CGS:.6e} cm³ g⁻¹ s⁻²")
    print(f"  M_sun = {MSUN_TO_G:.6e} g")
    print(f"  m_p = {PROTON_MASS_G:.6e} g")
    print(f"  pc = {PC_TO_CM:.6e} cm")
    
    print(f"\nFor c=1 relativistic units:")
    print(f"  Time = Length / c")
    print(f"  Velocity is dimensionless (in units of c)")
    print(f"  Energy = Mass × c² → [Energy] = [Mass]")
    print(f"  Pressure = Energy density → [Pressure] = [Mass] / [Length]³")
    
    print(f"\nTypical Scales:")
    print(f"  {'Scenario':<25} {'L (cm)':<15} {'t (s)':<15} {'ρ (g/cm³)':<15}")
    print(f"  {'-'*25} {'-'*15} {'-'*15} {'-'*15}")
    print(f"  {'NS merger':<25} {'1e6 (10 km)':<15} {'3e-5 (0.03 ms)':<15} {'1e14':<15}")
    print(f"  {'GRB afterglow':<25} {'1e17':<15} {'3e6 (40 days)':<15} {'2e-24':<15}")
    print(f"  {'AGN jet':<25} {'3e18 (1 pc)':<15} {'1e8 (3 yr)':<15} {'2e-24':<15}")
    print(f"  {'SNR (young)':<25} {'3e18 (1 pc)':<15} {'1e8 (3 yr)':<15} {'2e-24':<15}")


def main():
    print("="*70)
    print("Physical Unit Conversion Examples for SR-GSPH")
    print("="*70)
    print("\nAll simulations use natural units with c=1.")
    print("This script shows how to convert between code and physical units.")
    
    example_neutron_star_merger()
    example_grb_afterglow()
    example_relativistic_jet()
    example_conversion_table()
    
    print("\n" + "="*70)
    print("End of Examples")
    print("="*70)


if __name__ == '__main__':
    main()
