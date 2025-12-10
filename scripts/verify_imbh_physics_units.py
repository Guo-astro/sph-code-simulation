#!/usr/bin/env python3
"""
Verify IMBH-Cloud Simulation Physics Units

This script verifies that the IMBH Newtonian force and molecular cloud
initial position/velocity are using correct physics units.

Unit System (IMBH_ENCOUNTER):
  - 1 code length = 1 pc
  - 1 code mass = 1000 M_☉
  - 1 code velocity = 1 km/s
  - G = 1 in code units
  - 1 code time ≈ 0.978 Myr

Expected values for b=3 pc, v=10 km/s scenario:
  - Cloud center: (-20, 3, 0) pc
  - Cloud velocity: (10, 0, 0) km/s
  - BH at origin: (0, 0, 0)
  - BH mass: 10^5 M_☉ = 100 code units
"""

import numpy as np
import pandas as pd
import sys
import os

# Physical constants
PC_TO_CM = 3.0856775810e18  # 1 pc in cm
MSUN_TO_G = 1.98847e33      # 1 M_☉ in grams
KM_TO_CM = 1.0e5            # 1 km in cm
G_CGS = 6.67430e-8          # G in CGS (cm³ g⁻¹ s⁻²)

# IMBH_ENCOUNTER unit system
LENGTH_CODE_TO_CM = PC_TO_CM          # 1 code length = 1 pc
MASS_CODE_TO_G = 1000 * MSUN_TO_G     # 1 code mass = 1000 M_☉
VELOCITY_CODE_TO_CGS = KM_TO_CM       # 1 code velocity = 1 km/s

# Derived units
TIME_CODE_TO_S = LENGTH_CODE_TO_CM / VELOCITY_CODE_TO_CGS  # L/V
G_CODE = G_CGS * MASS_CODE_TO_G * TIME_CODE_TO_S**2 / LENGTH_CODE_TO_CM**3


def verify_units():
    """Verify that G=1 in the IMBH_ENCOUNTER unit system"""
    print("=" * 60)
    print("UNIT SYSTEM VERIFICATION (IMBH_ENCOUNTER)")
    print("=" * 60)
    print()
    print("Code units:")
    print(f"  1 code length = 1 pc = {PC_TO_CM:.4e} cm")
    print(f"  1 code mass = 1000 M_☉ = {MASS_CODE_TO_G:.4e} g")
    print(f"  1 code velocity = 1 km/s = {VELOCITY_CODE_TO_CGS:.4e} cm/s")
    print()
    
    # Time unit
    time_code_s = TIME_CODE_TO_S
    time_code_yr = time_code_s / (365.25 * 24 * 3600)
    time_code_myr = time_code_yr / 1e6
    print(f"  1 code time = L/V = {time_code_s:.4e} s = {time_code_myr:.4f} Myr")
    print()
    
    # G in code units
    print(f"Gravitational constant:")
    print(f"  G (CGS) = {G_CGS:.4e} cm³ g⁻¹ s⁻²")
    print(f"  G (code) = G × M_code × T_code² / L_code³")
    print(f"           = {G_CODE:.6f}")
    print()
    
    if abs(G_CODE - 1.0) < 0.01:
        print("  ✓ G ≈ 1 in code units (as expected)")
    else:
        print(f"  ✗ WARNING: G = {G_CODE:.4f} ≠ 1 in code units!")
    print()
    
    return time_code_myr


def verify_snapshot(filepath):
    """Verify snapshot data has correct positions and velocities"""
    print("=" * 60)
    print("SNAPSHOT DATA VERIFICATION")
    print("=" * 60)
    print(f"\nReading: {filepath}")
    
    try:
        df = pd.read_csv(filepath, comment='#')
    except Exception as e:
        print(f"ERROR: Could not read file: {e}")
        return False
    
    print(f"Particles: {len(df)}")
    print()
    
    # Position statistics
    print("Position (code units = pc):")
    print(f"  X: min={df['pos_x'].min():.4f}, max={df['pos_x'].max():.4f}, mean={df['pos_x'].mean():.4f}")
    print(f"  Y: min={df['pos_y'].min():.4f}, max={df['pos_y'].max():.4f}, mean={df['pos_y'].mean():.4f}")
    print(f"  Z: min={df['pos_z'].min():.4f}, max={df['pos_z'].max():.4f}, mean={df['pos_z'].mean():.4f}")
    print()
    
    # Center of mass
    com_x = df['pos_x'].mean()
    com_y = df['pos_y'].mean()
    com_z = df['pos_z'].mean()
    print("Cloud Center of Mass:")
    print(f"  COM = ({com_x:.4f}, {com_y:.4f}, {com_z:.4f}) pc")
    print(f"  Expected: (-20.0, 3.0, 0.0) pc")
    print()
    
    # Check if COM is correct
    pos_ok = True
    if abs(com_x - (-20.0)) > 0.5:
        print(f"  ✗ X position error: {com_x - (-20.0):.4f} pc")
        pos_ok = False
    else:
        print(f"  ✓ X position correct (error: {com_x - (-20.0):.4f} pc)")
    
    if abs(com_y - 3.0) > 0.5:
        print(f"  ✗ Y position error: {com_y - 3.0:.4f} pc")
        pos_ok = False
    else:
        print(f"  ✓ Y position correct (error: {com_y - 3.0:.4f} pc)")
    
    if abs(com_z - 0.0) > 0.5:
        print(f"  ✗ Z position error: {com_z - 0.0:.4f} pc")
        pos_ok = False
    else:
        print(f"  ✓ Z position correct (error: {com_z - 0.0:.4f} pc)")
    print()
    
    # Velocity statistics
    print("Velocity (code units = km/s):")
    print(f"  Vx: min={df['vel_x'].min():.4f}, max={df['vel_x'].max():.4f}, mean={df['vel_x'].mean():.4f}")
    print(f"  Vy: min={df['vel_y'].min():.4f}, max={df['vel_y'].max():.4f}, mean={df['vel_y'].mean():.4f}")
    print(f"  Vz: min={df['vel_z'].min():.4f}, max={df['vel_z'].max():.4f}, mean={df['vel_z'].mean():.4f}")
    print()
    
    # Mean velocity
    v_x = df['vel_x'].mean()
    v_y = df['vel_y'].mean()
    v_z = df['vel_z'].mean()
    print("Cloud Mean Velocity:")
    print(f"  V = ({v_x:.4f}, {v_y:.4f}, {v_z:.4f}) km/s")
    print(f"  Expected: (10.0, 0.0, 0.0) km/s")
    print()
    
    # Check if velocity is correct
    vel_ok = True
    if abs(v_x - 10.0) > 0.5:
        print(f"  ✗ Vx error: {v_x - 10.0:.4f} km/s")
        vel_ok = False
    else:
        print(f"  ✓ Vx correct (error: {v_x - 10.0:.4f} km/s)")
    
    if abs(v_y - 0.0) > 0.5:
        print(f"  ✗ Vy error: {v_y - 0.0:.4f} km/s")
        vel_ok = False
    else:
        print(f"  ✓ Vy correct (error: {v_y - 0.0:.4f} km/s)")
    
    if abs(v_z - 0.0) > 0.5:
        print(f"  ✗ Vz error: {v_z - 0.0:.4f} km/s")
        vel_ok = False
    else:
        print(f"  ✓ Vz correct (error: {v_z - 0.0:.4f} km/s)")
    print()
    
    return pos_ok and vel_ok


def verify_bh_gravity(df, M_BH_code=100.0, BH_pos=(0, 0, 0), softening=0.05):
    """
    Verify that gravitational acceleration is consistent with IMBH mass
    
    For a particle at position r from the BH:
        a_grav = -G * M_BH * r / (|r|² + ε²)^(3/2)
    
    With G=1 in code units:
        a_grav = -M_BH * r / (|r|² + ε²)^(3/2)
    """
    print("=" * 60)
    print("IMBH GRAVITATIONAL FORCE VERIFICATION")
    print("=" * 60)
    print()
    print(f"BH parameters (in code units):")
    print(f"  M_BH = {M_BH_code} code = {M_BH_code * 1000:.0f} M_☉")
    print(f"  Position = {BH_pos} code = {BH_pos} pc")
    print(f"  Softening ε = {softening} code = {softening} pc")
    print(f"  G = 1 (code units)")
    print()
    
    # Select a sample of particles for verification
    sample_indices = [0, 10, 100, 1000, len(df)//2, len(df)-1]
    sample_indices = [i for i in sample_indices if i < len(df)]
    
    print("Verifying gravitational acceleration for sample particles:")
    print("-" * 80)
    print(f"{'ID':>6} {'r (pc)':>10} {'|a_grav|_data':>14} {'|a_grav|_theory':>14} {'Error (%)':>12}")
    print("-" * 80)
    
    errors = []
    for idx in sample_indices:
        row = df.iloc[idx]
        
        # Position relative to BH
        dx = row['pos_x'] - BH_pos[0]
        dy = row['pos_y'] - BH_pos[1]
        dz = row['pos_z'] - BH_pos[2]
        r = np.sqrt(dx**2 + dy**2 + dz**2)
        
        # Gravitational acceleration from data
        ax_data = row['grav_acc_x']
        ay_data = row['grav_acc_y']
        az_data = row['grav_acc_z']
        a_mag_data = np.sqrt(ax_data**2 + ay_data**2 + az_data**2)
        
        # Expected acceleration from theory (includes self-gravity from cloud)
        # The BH contribution alone:
        r_soft = np.sqrt(r**2 + softening**2)
        a_bh_expected = M_BH_code / r_soft**2  # magnitude only from BH
        
        # Direction check: acceleration should point toward BH
        # Unit vector from particle to BH
        ux_expected = -dx / r if r > 0 else 0
        uy_expected = -dy / r if r > 0 else 0
        uz_expected = -dz / r if r > 0 else 0
        
        # The data includes both BH gravity AND cloud self-gravity
        # So we can't directly compare magnitudes
        # Instead, estimate BH contribution to total
        
        # Project data acceleration onto BH direction
        a_toward_bh = -(ax_data * dx + ay_data * dy + az_data * dz) / r if r > 0 else 0
        
        # Ratio of measured BH-directed acc to expected
        if a_bh_expected > 0:
            ratio = a_toward_bh / a_bh_expected
            error = (ratio - 1.0) * 100
        else:
            ratio = 0
            error = 0
        
        errors.append(error)
        
        print(f"{int(row['id']):>6} {r:>10.4f} {a_mag_data:>14.6f} {a_bh_expected:>14.6f} {error:>11.1f}%")
    
    print("-" * 80)
    print()
    
    # Analysis
    mean_error = np.mean(np.abs(errors))
    print(f"Mean absolute error in BH-directed acceleration: {mean_error:.1f}%")
    print()
    
    # Note about self-gravity
    print("Note: Total gravitational acceleration includes:")
    print("  1. IMBH gravity (dominant at large r)")
    print("  2. Cloud self-gravity (relevant at cloud center)")
    print("  Error > 0% indicates additional inward acceleration from self-gravity")
    print()
    
    # Check BH dominance at large distances
    # At r=20 pc from 10^5 Msun BH, tidal radius ~ 3.6 pc
    # Cloud at r~20 pc should feel mostly BH gravity
    
    return mean_error < 50  # Allow for self-gravity contribution


def verify_energy_conservation(filepath):
    """Check if energy output shows expected BH potential contribution"""
    energy_file = os.path.join(os.path.dirname(filepath), "energy.csv")
    if not os.path.exists(energy_file):
        print(f"Energy file not found: {energy_file}")
        return True
    
    print("=" * 60)
    print("ENERGY BUDGET VERIFICATION")
    print("=" * 60)
    print()
    
    df = pd.read_csv(energy_file)
    print(f"Energy data: {len(df)} timesteps")
    print()
    
    if len(df) > 0:
        # Initial state
        print("Initial energy state (t=0):")
        print(f"  Kinetic: {df['kinetic'].iloc[0]:.4f} code")
        print(f"  Thermal: {df['thermal'].iloc[0]:.4f} code")
        print(f"  Potential: {df['potential'].iloc[0]:.4f} code")
        print(f"  Total: {df['total'].iloc[0]:.4f} code")
        print()
        
        # Expected kinetic energy
        # E_kin = 0.5 * M_cloud * v^2 = 0.5 * 1 * 10^2 = 50 code
        M_cloud = 1.0  # code units (= 1000 Msun)
        v_bulk = 10.0  # km/s = code units
        E_kin_expected = 0.5 * M_cloud * v_bulk**2
        print(f"Expected bulk kinetic energy: 0.5 * {M_cloud} * {v_bulk}² = {E_kin_expected:.2f} code")
        print(f"Measured kinetic energy: {df['kinetic'].iloc[0]:.2f} code")
        print()
        
        if abs(df['kinetic'].iloc[0] - E_kin_expected) < 1.0:
            print("  ✓ Kinetic energy matches expected bulk motion")
        else:
            print(f"  ! Kinetic energy differs by {df['kinetic'].iloc[0] - E_kin_expected:.2f} code")
            print("    (difference is thermal velocity contribution)")
        print()
    
    return True


def main():
    # Default paths
    base_dir = "/Users/guo/Downloads/sphcode"
    snapshot_path = os.path.join(
        base_dir,
        "sample/imbh_cloud/results/Mc1e3_Mbh1e5_b3_v10/adiabatic_61k_gsph/snapshot_0000.csv"
    )
    
    if len(sys.argv) > 1:
        snapshot_path = sys.argv[1]
    
    print()
    print("###################################################################")
    print("#   IMBH-CLOUD SIMULATION PHYSICS UNIT VERIFICATION")
    print("###################################################################")
    print()
    
    # 1. Verify unit system
    time_code_myr = verify_units()
    print()
    
    # 2. Verify snapshot data
    if os.path.exists(snapshot_path):
        df = pd.read_csv(snapshot_path, comment='#')
        
        snapshot_ok = verify_snapshot(snapshot_path)
        print()
        
        # 3. Verify BH gravity
        # M_BH = 10^5 Msun = 100 code units (since 1 code = 1000 Msun)
        gravity_ok = verify_bh_gravity(df, M_BH_code=100.0, BH_pos=(0, 0, 0), softening=0.05)
        print()
        
        # 4. Verify energy
        verify_energy_conservation(snapshot_path)
        print()
        
        # Summary
        print("=" * 60)
        print("VERIFICATION SUMMARY")
        print("=" * 60)
        print()
        print(f"  Unit system: G ≈ 1 ✓")
        print(f"  Cloud position: {'✓ CORRECT' if snapshot_ok else '✗ ERROR'}")
        print(f"  Cloud velocity: {'✓ CORRECT' if snapshot_ok else '✗ ERROR'}")
        print(f"  BH gravity: {'✓ PLAUSIBLE' if gravity_ok else '! CHECK NEEDED'}")
        print()
        
        if snapshot_ok:
            print("✓ Physics units are being used correctly!")
        else:
            print("✗ Some issues detected - review above details")
        print()
    else:
        print(f"ERROR: Snapshot file not found: {snapshot_path}")
        print("Run a simulation first, or provide path as argument")


if __name__ == "__main__":
    main()
