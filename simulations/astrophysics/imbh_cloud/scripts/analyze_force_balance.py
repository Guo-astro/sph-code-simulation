#!/usr/bin/env python3
"""
Force Balance Analysis for Bonnor-Ebert Sphere

This script analyzes the force balance in a Bonnor-Ebert sphere simulation
to identify the root cause of central density deficit.

For hydrostatic equilibrium:
    Pressure gradient force = Gravitational force (inward)
    dP/dr = -ρ * G * M(r) / r²
    
The SPH acceleration should be:
    a_total = a_pressure + a_gravity ≈ 0

If a_total_radial > 0: Net outward force → expansion
If a_total_radial < 0: Net inward force → contraction

Author: Force balance diagnostic tool
"""

import math
import csv
import sys
from collections import defaultdict

def load_csv(filename):
    """Load CSV file, skipping comment lines"""
    data = defaultdict(list)
    with open(filename, 'r') as f:
        # Skip comment lines
        for line in f:
            if not line.startswith('#'):
                # This is the header
                headers = line.strip().split(',')
                break
        
        reader = csv.DictReader(f, fieldnames=headers)
        for row in reader:
            for col in headers:
                try:
                    data[col].append(float(row[col]))
                except (ValueError, KeyError):
                    data[col].append(0.0)
    return dict(data)

def analyze_snapshot(filename, verbose=True):
    """Comprehensive force balance analysis"""
    
    data = load_csv(filename)
    n = len(data['pos_x'])
    
    # Extract data
    x = data['pos_x']
    y = data['pos_y']
    z = data['pos_z']
    vx = data['vel_x']
    vy = data['vel_y']
    vz = data['vel_z']
    ax_total = data['acc_x']
    ay_total = data['acc_y']
    az_total = data['acc_z']
    ax_grav = data['grav_acc_x']
    ay_grav = data['grav_acc_y']
    az_grav = data['grav_acc_z']
    dens = data['dens']
    pres = data['pres']
    mass = data['mass']
    sml = data['sml']
    is_ghost = data.get('is_ghost', [0]*n)
    
    # Compute derived quantities for each particle
    results = {
        'r': [],           # radius
        'v_r': [],         # radial velocity
        'a_total': [],     # total acceleration magnitude
        'a_grav': [],      # gravitational acceleration magnitude
        'a_pres': [],      # pressure acceleration magnitude
        'a_total_r': [],   # radial component of total acc
        'a_grav_r': [],    # radial component of grav acc
        'a_pres_r': [],    # radial component of pressure acc
        'dens': [],
        'pres': [],
        'is_ghost': [],
        'sml': [],
    }
    
    for i in range(n):
        # Radius
        r = math.sqrt(x[i]**2 + y[i]**2 + z[i]**2)
        results['r'].append(r)
        
        # Skip very central particles (r ≈ 0)
        if r < 1e-6:
            r_unit = (0, 0, 0)
        else:
            r_unit = (x[i]/r, y[i]/r, z[i]/r)
        
        # Radial velocity (positive = outward)
        v_r = vx[i]*r_unit[0] + vy[i]*r_unit[1] + vz[i]*r_unit[2]
        results['v_r'].append(v_r)
        
        # Pressure acceleration = total - gravity
        ax_pres = ax_total[i] - ax_grav[i]
        ay_pres = ay_total[i] - ay_grav[i]
        az_pres = az_total[i] - az_grav[i]
        
        # Magnitudes
        a_total = math.sqrt(ax_total[i]**2 + ay_total[i]**2 + az_total[i]**2)
        a_grav = math.sqrt(ax_grav[i]**2 + ay_grav[i]**2 + az_grav[i]**2)
        a_pres = math.sqrt(ax_pres**2 + ay_pres**2 + az_pres**2)
        
        results['a_total'].append(a_total)
        results['a_grav'].append(a_grav)
        results['a_pres'].append(a_pres)
        
        # Radial components (positive = outward)
        a_total_r = ax_total[i]*r_unit[0] + ay_total[i]*r_unit[1] + az_total[i]*r_unit[2]
        a_grav_r = ax_grav[i]*r_unit[0] + ay_grav[i]*r_unit[1] + az_grav[i]*r_unit[2]
        a_pres_r = ax_pres*r_unit[0] + ay_pres*r_unit[1] + az_pres*r_unit[2]
        
        results['a_total_r'].append(a_total_r)
        results['a_grav_r'].append(a_grav_r)
        results['a_pres_r'].append(a_pres_r)
        
        results['dens'].append(dens[i])
        results['pres'].append(pres[i])
        results['is_ghost'].append(is_ghost[i])
        results['sml'].append(sml[i])
    
    return results

def print_analysis(results, title=""):
    """Print comprehensive analysis"""
    
    n = len(results['r'])
    
    # Separate cloud and ghost particles
    cloud_idx = [i for i in range(n) if results['is_ghost'][i] == 0]
    ghost_idx = [i for i in range(n) if results['is_ghost'][i] == 1]
    
    # Define radial bins
    bins = [
        (0.0, 0.05, "Very center (r<0.05)"),
        (0.05, 0.1, "Inner core (0.05-0.1)"),
        (0.1, 0.2, "Core (0.1-0.2)"),
        (0.2, 0.4, "Middle (0.2-0.4)"),
        (0.4, 0.6, "Outer (0.4-0.6)"),
        (0.6, 0.75, "Edge (0.6-0.75)"),
        (0.75, 1.0, "Boundary (0.75-1.0)"),
    ]
    
    print("=" * 80)
    print(f"FORCE BALANCE ANALYSIS {title}")
    print("=" * 80)
    print(f"Total particles: {n}")
    print(f"Cloud particles: {len(cloud_idx)}")
    print(f"Ghost particles: {len(ghost_idx)}")
    print()
    
    print("-" * 80)
    print("RADIAL ANALYSIS (Cloud particles only)")
    print("-" * 80)
    print(f"{'Region':<22} {'N':>6} {'ρ_mean':>10} {'a_grav_r':>12} {'a_pres_r':>12} {'a_net_r':>12} {'v_r':>10}")
    print(f"{'':22} {'':>6} {'M☉/pc³':>10} {'pc/Myr²':>12} {'pc/Myr²':>12} {'pc/Myr²':>12} {'pc/Myr':>10}")
    print("-" * 80)
    
    for r_min, r_max, label in bins:
        idx = [i for i in cloud_idx if r_min <= results['r'][i] < r_max]
        if len(idx) == 0:
            continue
            
        dens_mean = sum(results['dens'][i] for i in idx) / len(idx)
        a_grav_r_mean = sum(results['a_grav_r'][i] for i in idx) / len(idx)
        a_pres_r_mean = sum(results['a_pres_r'][i] for i in idx) / len(idx)
        a_net_r_mean = sum(results['a_total_r'][i] for i in idx) / len(idx)
        v_r_mean = sum(results['v_r'][i] for i in idx) / len(idx)
        
        # Interpretation
        direction = "→OUT" if a_net_r_mean > 0 else "←IN "
        
        print(f"{label:<22} {len(idx):>6} {dens_mean:>10.3f} {a_grav_r_mean:>12.6f} {a_pres_r_mean:>12.6f} {a_net_r_mean:>12.6f} {v_r_mean:>10.6f} {direction}")
    
    print()
    print("-" * 80)
    print("FORCE IMBALANCE INTERPRETATION")
    print("-" * 80)
    
    # Analyze center (r < 0.1)
    center_idx = [i for i in cloud_idx if results['r'][i] < 0.1]
    if center_idx:
        a_grav_r_center = [results['a_grav_r'][i] for i in center_idx]
        a_pres_r_center = [results['a_pres_r'][i] for i in center_idx]
        a_net_r_center = [results['a_total_r'][i] for i in center_idx]
        
        grav_mean = sum(a_grav_r_center) / len(a_grav_r_center)
        pres_mean = sum(a_pres_r_center) / len(a_pres_r_center)
        net_mean = sum(a_net_r_center) / len(a_net_r_center)
        
        print(f"CENTER (r < 0.1 pc, N={len(center_idx)}):")
        print(f"  Gravity (radial):  {grav_mean:+.6f} pc/Myr² (should be negative/inward)")
        print(f"  Pressure (radial): {pres_mean:+.6f} pc/Myr² (should be positive/outward)")
        print(f"  NET (radial):      {net_mean:+.6f} pc/Myr² (should be ~0 for equilibrium)")
        print()
        
        if net_mean > 1e-6:
            print("  ⚠️  NET OUTWARD FORCE → Particles will expand outward")
            if pres_mean > abs(grav_mean):
                print("  🔍 ROOT CAUSE: Pressure force EXCEEDS gravity")
                print("      Possible reasons:")
                print("      1. Initial density too low → pressure gradient too steep")
                print("      2. Ghost particles not providing enough confinement")
                print("      3. SPH E0 error at density gradients")
            else:
                print("  🔍 Gravity is dominant but net is still outward - check numerics")
        elif net_mean < -1e-6:
            print("  ⚠️  NET INWARD FORCE → Particles will collapse inward")
        else:
            print("  ✓  Approximately in equilibrium")
    
    # Analyze edge (r ~ 0.6-0.75)
    edge_idx = [i for i in cloud_idx if 0.6 <= results['r'][i] < 0.75]
    if edge_idx:
        a_grav_r_edge = [results['a_grav_r'][i] for i in edge_idx]
        a_pres_r_edge = [results['a_pres_r'][i] for i in edge_idx]
        a_net_r_edge = [results['a_total_r'][i] for i in edge_idx]
        
        grav_mean = sum(a_grav_r_edge) / len(a_grav_r_edge)
        pres_mean = sum(a_pres_r_edge) / len(a_pres_r_edge)
        net_mean = sum(a_net_r_edge) / len(a_net_r_edge)
        
        print()
        print(f"EDGE (r = 0.6-0.75 pc, N={len(edge_idx)}):")
        print(f"  Gravity (radial):  {grav_mean:+.6f} pc/Myr² (inward)")
        print(f"  Pressure (radial): {pres_mean:+.6f} pc/Myr² (outward)")
        print(f"  NET (radial):      {net_mean:+.6f} pc/Myr² (should be ~0)")
    
    # Ghost analysis
    if ghost_idx:
        print()
        print("-" * 80)
        print("GHOST PARTICLE ANALYSIS")
        print("-" * 80)
        
        ghost_r = [results['r'][i] for i in ghost_idx]
        ghost_dens = [results['dens'][i] for i in ghost_idx]
        ghost_pres = [results['pres'][i] for i in ghost_idx]
        
        print(f"Ghost particles: {len(ghost_idx)}")
        print(f"Ghost radius range: {min(ghost_r):.3f} - {max(ghost_r):.3f} pc")
        print(f"Ghost density: {sum(ghost_dens)/len(ghost_dens):.3f} M☉/pc³")
        print(f"Ghost pressure: {sum(ghost_pres)/len(ghost_pres):.3f}")
        
        # Check gap between cloud edge and ghost inner edge
        cloud_r_max = max(results['r'][i] for i in cloud_idx)
        ghost_r_min = min(ghost_r)
        gap = ghost_r_min - cloud_r_max
        
        print(f"\nCloud outer radius: {cloud_r_max:.4f} pc")
        print(f"Ghost inner radius: {ghost_r_min:.4f} pc")
        print(f"Gap: {gap:.4f} pc")
        
        if gap > 0.05:
            print("  ⚠️  Large gap between cloud and ghost envelope!")
            print("      This can cause pressure drop at boundary")
    
    return results

def compare_with_analytical(results):
    """Compare with analytical Bonnor-Ebert solution"""
    
    print()
    print("=" * 80)
    print("COMPARISON WITH ANALYTICAL SOLUTION")
    print("=" * 80)
    
    # BE sphere parameters (from config)
    n_center = 1800  # cm^-3
    T = 7  # K
    mu = 1.27
    xi_s = 6.0
    
    # Expected central density in M☉/pc³
    # From previous calculations: ~56.5 M☉/pc³
    rho_center_analytical = 56.5
    
    # Get actual central density
    cloud_idx = [i for i in range(len(results['r'])) if results['is_ghost'][i] == 0]
    center_idx = [i for i in cloud_idx if results['r'][i] < 0.05]
    
    if center_idx:
        rho_center_sim = sum(results['dens'][i] for i in center_idx) / len(center_idx)
        deficit = (rho_center_analytical - rho_center_sim) / rho_center_analytical * 100
        
        print(f"Analytical central density: {rho_center_analytical:.2f} M☉/pc³")
        print(f"Simulated central density:  {rho_center_sim:.2f} M☉/pc³")
        print(f"Deficit: {deficit:.1f}%")
        
        # Estimate pressure
        # P = ρ * c_s² where c_s = 0.213 km/s = 0.218 pc/Myr
        c_s = 0.218  # pc/Myr
        P_analytical = rho_center_analytical * c_s**2
        P_sim = sum(results['pres'][i] for i in center_idx) / len(center_idx)
        
        print(f"\nExpected central pressure (ρ*c_s²): {P_analytical:.4f}")
        print(f"Simulated central pressure: {P_sim:.4f}")


def main():
    if len(sys.argv) < 2:
        print("Usage: python analyze_force_balance.py <snapshot.csv> [snapshot2.csv ...]")
        print("\nAnalyzing default snapshots...")
        
        import os
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        results_dir = os.path.join(base_dir, "results", "phase2_hydrostatic")
        
        # Analyze first snapshot
        snap1 = os.path.join(results_dir, "snapshot_0001.csv")
        if os.path.exists(snap1):
            print(f"\nAnalyzing: {snap1}")
            results1 = analyze_snapshot(snap1)
            print_analysis(results1, title="(t=0)")
            compare_with_analytical(results1)
            
            # Compare with later snapshot
            snap5 = os.path.join(results_dir, "snapshot_0005.csv")
            if os.path.exists(snap5):
                print("\n" + "=" * 80)
                print("\nAnalyzing: {snap5}")
                results5 = analyze_snapshot(snap5)
                print_analysis(results5, title="(t~50 Myr)")
        else:
            print(f"File not found: {snap1}")
    else:
        for filename in sys.argv[1:]:
            print(f"\nAnalyzing: {filename}")
            results = analyze_snapshot(filename)
            print_analysis(results)
            compare_with_analytical(results)


if __name__ == "__main__":
    main()
