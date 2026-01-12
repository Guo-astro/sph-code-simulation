#!/usr/bin/env python3
"""
Detailed Force Analysis for Bonnor-Ebert Sphere Edge Region

This script analyzes why pressure exceeds gravity at the cloud boundary.

Key questions:
1. What is the density gradient at the edge?
2. How does SPH compute pressure gradient?
3. Is the ghost envelope providing correct boundary conditions?
"""

import math
import csv
from collections import defaultdict

def load_csv(filename):
    """Load CSV file, skipping comment lines"""
    data = defaultdict(list)
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith('#'):
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

def main():
    filename = "simulations/astrophysics/imbh_cloud/results/phase2_hydrostatic/snapshot_0001.csv"
    data = load_csv(filename)
    n = len(data['pos_x'])
    
    # Calculate radius and forces for each particle
    particles = []
    for i in range(n):
        x, y, z = data['pos_x'][i], data['pos_y'][i], data['pos_z'][i]
        r = math.sqrt(x**2 + y**2 + z**2)
        
        if r < 1e-6:
            r_unit = (0, 0, 0)
        else:
            r_unit = (x/r, y/r, z/r)
        
        # Accelerations
        ax_t = data['acc_x'][i]
        ay_t = data['acc_y'][i]
        az_t = data['acc_z'][i]
        ax_g = data['grav_acc_x'][i]
        ay_g = data['grav_acc_y'][i]
        az_g = data['grav_acc_z'][i]
        
        # Pressure accel = total - gravity
        ax_p = ax_t - ax_g
        ay_p = ay_t - ay_g
        az_p = az_t - az_g
        
        # Radial components
        a_grav_r = ax_g*r_unit[0] + ay_g*r_unit[1] + az_g*r_unit[2]
        a_pres_r = ax_p*r_unit[0] + ay_p*r_unit[1] + az_p*r_unit[2]
        a_total_r = ax_t*r_unit[0] + ay_t*r_unit[1] + az_t*r_unit[2]
        
        particles.append({
            'r': r,
            'dens': data['dens'][i],
            'pres': data['pres'][i],
            'a_grav_r': a_grav_r,
            'a_pres_r': a_pres_r,
            'a_total_r': a_total_r,
            'sml': data['sml'][i],
            'is_ghost': data['is_ghost'][i],
            'neighbor': data['neighbor'][i],
        })
    
    # Sort by radius
    particles.sort(key=lambda p: p['r'])
    
    # Analyze edge region in detail (0.6 < r < 0.85)
    print("=" * 100)
    print("DETAILED EDGE ANALYSIS (0.6 < r < 0.85 pc)")
    print("=" * 100)
    print()
    
    print(f"{'r':>8} {'dens':>10} {'pres':>10} {'a_grav_r':>12} {'a_pres_r':>12} {'a_net_r':>12} {'sml':>8} {'Nngb':>6} {'ghost':>6}")
    print("-" * 100)
    
    edge_particles = [p for p in particles if 0.6 < p['r'] < 0.85]
    
    # Sample every 50th particle
    for i in range(0, len(edge_particles), 50):
        p = edge_particles[i]
        ghost_str = "GHOST" if p['is_ghost'] == 1 else ""
        print(f"{p['r']:8.4f} {p['dens']:10.4f} {p['pres']:10.4f} {p['a_grav_r']:+12.6f} {p['a_pres_r']:+12.6f} {p['a_total_r']:+12.6f} {p['sml']:8.4f} {int(p['neighbor']):6d} {ghost_str}")
    
    # Detailed analysis of the cloud-ghost boundary
    print()
    print("=" * 100)
    print("CLOUD-GHOST BOUNDARY ANALYSIS")
    print("=" * 100)
    
    # Find cloud edge (last cloud particle)
    cloud_particles = [p for p in particles if p['is_ghost'] == 0]
    ghost_particles = [p for p in particles if p['is_ghost'] == 1]
    
    cloud_edge = max(p['r'] for p in cloud_particles)
    ghost_inner = min(p['r'] for p in ghost_particles)
    
    print(f"\nCloud outer radius: {cloud_edge:.4f} pc")
    print(f"Ghost inner radius: {ghost_inner:.4f} pc")
    print(f"Gap: {ghost_inner - cloud_edge:.4f} pc")
    
    # Particles near the gap
    print("\n--- Cloud particles near edge (r > 0.7) ---")
    near_edge_cloud = [p for p in cloud_particles if p['r'] > 0.7]
    near_edge_cloud.sort(key=lambda p: p['r'])
    
    print(f"{'r':>8} {'dens':>10} {'pres':>10} {'a_grav_r':>12} {'a_pres_r':>12} {'a_net_r':>12} {'sml':>8} {'ratio':>8}")
    print("-" * 100)
    
    for p in near_edge_cloud[-20:]:  # Last 20 cloud particles
        ratio = abs(p['a_pres_r'] / p['a_grav_r']) if abs(p['a_grav_r']) > 1e-10 else 0
        print(f"{p['r']:8.4f} {p['dens']:10.4f} {p['pres']:10.4f} {p['a_grav_r']:+12.6f} {p['a_pres_r']:+12.6f} {p['a_total_r']:+12.6f} {p['sml']:8.4f} {ratio:8.3f}")
    
    print("\n--- Ghost particles near edge (r < 0.85) ---")
    near_edge_ghost = [p for p in ghost_particles if p['r'] < 0.85]
    near_edge_ghost.sort(key=lambda p: p['r'])
    
    print(f"{'r':>8} {'dens':>10} {'pres':>10} {'a_grav_r':>12} {'a_pres_r':>12} {'a_net_r':>12} {'sml':>8}")
    print("-" * 100)
    
    for p in near_edge_ghost[:10]:  # First 10 ghost particles
        print(f"{p['r']:8.4f} {p['dens']:10.4f} {p['pres']:10.4f} {p['a_grav_r']:+12.6f} {p['a_pres_r']:+12.6f} {p['a_total_r']:+12.6f} {p['sml']:8.4f}")
    
    # Analyze pressure ratio (should be ~1.0 for equilibrium)
    print()
    print("=" * 100)
    print("PRESSURE/GRAVITY RATIO ANALYSIS")
    print("=" * 100)
    print("(Ratio = |a_pres_r| / |a_grav_r|, should be ~1.0 for equilibrium)")
    print()
    
    bins = [
        (0.0, 0.1, "Center"),
        (0.1, 0.2, "Inner"),
        (0.2, 0.4, "Middle"),
        (0.4, 0.6, "Outer"),
        (0.6, 0.7, "Edge-inner"),
        (0.7, 0.75, "Edge-outer"),
    ]
    
    for r_min, r_max, label in bins:
        bin_particles = [p for p in cloud_particles if r_min <= p['r'] < r_max]
        if not bin_particles:
            continue
        
        ratios = []
        for p in bin_particles:
            if abs(p['a_grav_r']) > 1e-10:
                ratios.append(abs(p['a_pres_r']) / abs(p['a_grav_r']))
        
        if ratios:
            mean_ratio = sum(ratios) / len(ratios)
            excess = (mean_ratio - 1.0) * 100
            status = "OK" if abs(excess) < 5 else ("PRESSURE EXCESS" if excess > 0 else "GRAVITY EXCESS")
            print(f"{label:15} (N={len(bin_particles):5}): ratio = {mean_ratio:.4f}  ({excess:+.1f}%)  {status}")
    
    # Calculate analytical expectations
    print()
    print("=" * 100)
    print("ANALYTICAL EXPECTATIONS")
    print("=" * 100)
    
    # For isothermal gas in hydrostatic equilibrium:
    # dP/dr = -ρ * g
    # P = ρ * c_s²  (isothermal)
    # So: c_s² * dρ/dr = -ρ * g
    # => a_pres = -c_s² * (1/ρ) * dρ/dr = -c_s² * d(ln ρ)/dr
    # => a_grav = -g = -G * M(r) / r²
    
    # At equilibrium: a_pres + a_grav = 0
    # => c_s² * d(ln ρ)/dr = G * M(r) / r²
    
    c_s = 0.218  # pc/Myr
    c_s_sq = c_s**2
    
    print(f"Sound speed: c_s = {c_s:.4f} pc/Myr")
    print(f"c_s² = {c_s_sq:.6f} (pc/Myr)²")
    print()
    
    # Check pressure formula: P = ρ * c_s² for isothermal
    center_particles = [p for p in cloud_particles if p['r'] < 0.1]
    dens_center = sum(p['dens'] for p in center_particles) / len(center_particles)
    pres_center = sum(p['pres'] for p in center_particles) / len(center_particles)
    expected_pres = dens_center * c_s_sq
    
    print(f"Center density: {dens_center:.4f} M☉/pc³")
    print(f"Center pressure (simulated): {pres_center:.4f}")
    print(f"Center pressure (expected ρ*c_s²): {expected_pres:.4f}")
    print(f"Ratio (sim/expected): {pres_center/expected_pres:.4f}")
    
    # Check if EOS is correct
    eos_ratio = pres_center / dens_center / c_s_sq
    print(f"\nEOS check: P/(ρ*c_s²) = {eos_ratio:.4f} (should be 1.0 for isothermal)")
    
    if abs(eos_ratio - 1.0) > 0.01:
        print("⚠️  EOS MISMATCH - Pressure is not following P = ρ * c_s²!")
        print("   This could be due to:")
        print("   1. Wrong gamma (using adiabatic instead of isothermal)")
        print("   2. Wrong internal energy initialization")
        print("   3. Wrong sound speed calculation")

if __name__ == "__main__":
    main()
