#!/usr/bin/env python3
"""
Hydrostatic Equilibrium Force Balance Diagnostic

This script performs a detailed analysis of the force imbalance in a 
Bonnor-Ebert sphere simulation to identify the root cause.

Physics Background:
==================
For isothermal hydrostatic equilibrium:
    ∇P = -ρ ∇Φ  (pressure gradient balances gravity)
    
In SPH, the pressure acceleration is:
    a_pres = -Σ_j m_j [(P_i/ρ_i² Ω_i) ∇W_i + (P_j/ρ_j² Ω_j) ∇W_j]
    
Where Ω is the grad-h correction factor.

The gravitational acceleration is:
    a_grav = -G Σ_j m_j r_ij / |r_ij|³  (or from tree code)

For equilibrium: a_pres + a_grav ≈ 0

Possible error sources:
1. SPH E0 error (∝ ∇²ρ, worst at density gradients)
2. Grad-h correction errors at boundaries
3. Kernel truncation at cloud edge (asymmetric neighbors)
4. Ghost particle boundary pressure mismatch
5. Initial condition deviation from analytical BE profile

This diagnostic will quantify each error source.
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


def analyze_be_equilibrium(filename):
    """Comprehensive equilibrium analysis"""
    
    print("=" * 80)
    print("BONNOR-EBERT HYDROSTATIC EQUILIBRIUM DIAGNOSTIC")
    print("=" * 80)
    print(f"Analyzing: {filename}")
    print()
    
    data = load_csv(filename)
    n = len(data['pos_x'])
    
    # Physical constants from config
    G = 0.00430091  # pc³ / (M☉ * Myr²)
    gamma = 1.0001
    c_s_target = 0.218  # pc/Myr (target sound speed)
    
    # Extract data
    x = data['pos_x']
    y = data['pos_y']
    z = data['pos_z']
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
    ene = data['ene']
    gradh = data['gradh']
    neighbor = data['neighbor']
    is_ghost = data.get('is_ghost', [0]*n)
    
    # Compute derived quantities
    particles = []
    for i in range(n):
        r = math.sqrt(x[i]**2 + y[i]**2 + z[i]**2)
        
        if r < 1e-6:
            r_unit = (0, 0, 0)
        else:
            r_unit = (x[i]/r, y[i]/r, z[i]/r)
        
        # Pressure acceleration = total - gravity
        ax_pres = ax_total[i] - ax_grav[i]
        ay_pres = ay_total[i] - ay_grav[i]
        az_pres = az_total[i] - az_grav[i]
        
        # Radial components (positive = outward)
        a_grav_r = ax_grav[i]*r_unit[0] + ay_grav[i]*r_unit[1] + az_grav[i]*r_unit[2]
        a_pres_r = ax_pres*r_unit[0] + ay_pres*r_unit[1] + az_pres*r_unit[2]
        a_total_r = ax_total[i]*r_unit[0] + ay_total[i]*r_unit[1] + az_total[i]*r_unit[2]
        
        # Sound speed from internal energy
        c_s_actual = math.sqrt((gamma - 1) * ene[i])
        
        particles.append({
            'r': r,
            'dens': dens[i],
            'pres': pres[i],
            'mass': mass[i],
            'sml': sml[i],
            'ene': ene[i],
            'gradh': gradh[i],
            'neighbor': neighbor[i],
            'is_ghost': is_ghost[i],
            'a_grav_r': a_grav_r,
            'a_pres_r': a_pres_r,
            'a_total_r': a_total_r,
            'c_s': c_s_actual,
        })
    
    # Separate cloud and ghost
    cloud = [p for p in particles if p['is_ghost'] == 0]
    ghost = [p for p in particles if p['is_ghost'] == 1]
    
    print(f"Total particles: {n}")
    print(f"Cloud particles: {len(cloud)}")
    print(f"Ghost particles: {len(ghost)}")
    print()
    
    # ========================================
    # DIAGNOSTIC 1: EOS Consistency Check
    # ========================================
    print("-" * 80)
    print("DIAGNOSTIC 1: EQUATION OF STATE CONSISTENCY")
    print("-" * 80)
    print("For isothermal gas: P = ρ * c_s²")
    print("In SPH: P = (γ-1) * ρ * u")
    print()
    
    # Check EOS for cloud particles
    center_cloud = [p for p in cloud if p['r'] < 0.2]
    if center_cloud:
        dens_avg = sum(p['dens'] for p in center_cloud) / len(center_cloud)
        pres_avg = sum(p['pres'] for p in center_cloud) / len(center_cloud)
        ene_avg = sum(p['ene'] for p in center_cloud) / len(center_cloud)
        c_s_avg = sum(p['c_s'] for p in center_cloud) / len(center_cloud)
        
        # Expected pressure
        P_expected = dens_avg * c_s_target**2
        P_from_eos = (gamma - 1) * dens_avg * ene_avg
        
        print(f"Center (r<0.2) statistics:")
        print(f"  Mean density: {dens_avg:.4f} M☉/pc³")
        print(f"  Mean pressure: {pres_avg:.4f}")
        print(f"  Mean internal energy: {ene_avg:.4f}")
        print(f"  Mean sound speed: {c_s_avg:.4f} pc/Myr (target: {c_s_target:.4f})")
        print()
        print(f"  Expected P (ρ*c_s_target²): {P_expected:.4f}")
        print(f"  Actual P from EOS: {P_from_eos:.4f}")
        print(f"  Ratio: {P_from_eos/P_expected:.4f}")
        
        if abs(P_from_eos/P_expected - 1) > 0.05:
            print("  ⚠️  PRESSURE MISMATCH: EOS not matching isothermal expectation")
        else:
            print("  ✓  EOS is consistent with isothermal")
    
    # ========================================
    # DIAGNOSTIC 2: Grad-h Correction Factor
    # ========================================
    print()
    print("-" * 80)
    print("DIAGNOSTIC 2: GRAD-H CORRECTION FACTOR")
    print("-" * 80)
    print("Ω = 1 / (1 + h/(D*ρ) * dρ/dh)")
    print("For uniform density: Ω = 1")
    print("At density gradients: Ω ≠ 1 (correction needed)")
    print()
    
    bins = [
        (0.0, 0.1, "Center"),
        (0.1, 0.3, "Inner"),
        (0.3, 0.5, "Middle"),
        (0.5, 0.7, "Outer"),
        (0.7, 0.75, "Edge"),
    ]
    
    print(f"{'Region':<12} {'N':>6} {'gradh_mean':>12} {'gradh_std':>12} {'gradh_min':>12} {'gradh_max':>12}")
    print("-" * 70)
    
    for r_min, r_max, label in bins:
        bin_particles = [p for p in cloud if r_min <= p['r'] < r_max]
        if not bin_particles:
            continue
        
        gradh_values = [p['gradh'] for p in bin_particles]
        mean_gradh = sum(gradh_values) / len(gradh_values)
        std_gradh = math.sqrt(sum((g - mean_gradh)**2 for g in gradh_values) / len(gradh_values))
        min_gradh = min(gradh_values)
        max_gradh = max(gradh_values)
        
        print(f"{label:<12} {len(bin_particles):>6} {mean_gradh:>12.4f} {std_gradh:>12.4f} {min_gradh:>12.4f} {max_gradh:>12.4f}")
    
    # ========================================
    # DIAGNOSTIC 3: Neighbor Count Analysis
    # ========================================
    print()
    print("-" * 80)
    print("DIAGNOSTIC 3: NEIGHBOR COUNT ANALYSIS")
    print("-" * 80)
    print("Target neighbors: ~50")
    print("Low neighbor count at boundaries can cause pressure gradient errors")
    print()
    
    print(f"{'Region':<12} {'N':>6} {'Nngb_mean':>12} {'Nngb_min':>12} {'Nngb_max':>12}")
    print("-" * 60)
    
    for r_min, r_max, label in bins:
        bin_particles = [p for p in cloud if r_min <= p['r'] < r_max]
        if not bin_particles:
            continue
        
        neighbor_values = [p['neighbor'] for p in bin_particles]
        mean_ngb = sum(neighbor_values) / len(neighbor_values)
        min_ngb = min(neighbor_values)
        max_ngb = max(neighbor_values)
        
        print(f"{label:<12} {len(bin_particles):>6} {mean_ngb:>12.1f} {int(min_ngb):>12d} {int(max_ngb):>12d}")
    
    # ========================================
    # DIAGNOSTIC 4: Force Balance Analysis
    # ========================================
    print()
    print("-" * 80)
    print("DIAGNOSTIC 4: FORCE BALANCE (Gravity vs Pressure)")
    print("-" * 80)
    print("For equilibrium: |a_pres_r| / |a_grav_r| = 1.0")
    print("Ratio > 1: Pressure excess → expansion")
    print("Ratio < 1: Gravity excess → collapse")
    print()
    
    print(f"{'Region':<12} {'N':>6} {'|a_grav_r|':>12} {'|a_pres_r|':>12} {'Ratio':>10} {'Excess':>10}")
    print("-" * 70)
    
    for r_min, r_max, label in bins:
        bin_particles = [p for p in cloud if r_min <= p['r'] < r_max]
        if not bin_particles:
            continue
        
        a_grav_r_avg = sum(abs(p['a_grav_r']) for p in bin_particles) / len(bin_particles)
        a_pres_r_avg = sum(abs(p['a_pres_r']) for p in bin_particles) / len(bin_particles)
        
        if a_grav_r_avg > 1e-10:
            ratio = a_pres_r_avg / a_grav_r_avg
            excess = (ratio - 1) * 100
        else:
            ratio = 0
            excess = 0
        
        excess_str = f"+{excess:.1f}%" if excess >= 0 else f"{excess:.1f}%"
        print(f"{label:<12} {len(bin_particles):>6} {a_grav_r_avg:>12.6f} {a_pres_r_avg:>12.6f} {ratio:>10.4f} {excess_str:>10}")
    
    # ========================================
    # DIAGNOSTIC 5: Net Radial Acceleration
    # ========================================
    print()
    print("-" * 80)
    print("DIAGNOSTIC 5: NET RADIAL ACCELERATION")
    print("-" * 80)
    print("a_net_r = a_pres_r + a_grav_r (should be ~0)")
    print("Positive = outward force, Negative = inward force")
    print()
    
    print(f"{'Region':<12} {'N':>6} {'a_net_r':>12} {'a_net_r/|a_g|':>14} {'Direction':>12}")
    print("-" * 70)
    
    for r_min, r_max, label in bins:
        bin_particles = [p for p in cloud if r_min <= p['r'] < r_max]
        if not bin_particles:
            continue
        
        a_net_r_avg = sum(p['a_total_r'] for p in bin_particles) / len(bin_particles)
        a_grav_r_avg = sum(abs(p['a_grav_r']) for p in bin_particles) / len(bin_particles)
        
        if a_grav_r_avg > 1e-10:
            relative_net = a_net_r_avg / a_grav_r_avg * 100
        else:
            relative_net = 0
        
        direction = "→ OUT" if a_net_r_avg > 0 else "← IN"
        print(f"{label:<12} {len(bin_particles):>6} {a_net_r_avg:>+12.6f} {relative_net:>+13.1f}% {direction:>12}")
    
    # ========================================
    # DIAGNOSTIC 6: Analytical BE Profile Check
    # ========================================
    print()
    print("-" * 80)
    print("DIAGNOSTIC 6: ANALYTICAL BE PROFILE COMPARISON")
    print("-" * 80)
    
    # BE sphere parameters
    xi_s = 6.0
    R_cloud = 0.75  # pc (cloud radius at ξ_s)
    rho_center_analytical = 56.5  # M☉/pc³
    
    # Lane-Emden solution points (approximate)
    # ξ: 0, 1, 2, 3, 4, 5, 6
    # θ: 1, 0.9157, 0.7033, 0.4684, 0.2721, 0.1206, 0.0000
    # ρ/ρ_c = exp(-θ) for isothermal
    xi_points = [0, 1, 2, 3, 4, 5, 6]
    theta_values = [0, 0.0880, 0.3582, 0.7579, 1.2984, 2.0116, 3.0000]  # BE solution
    
    # Map to radius
    r_scale = R_cloud / xi_s  # pc per unit ξ
    
    print(f"BE sphere parameters: ξ_s = {xi_s}, R = {R_cloud} pc, ρ_c = {rho_center_analytical} M☉/pc³")
    print()
    print(f"{'r (pc)':<10} {'ξ':>8} {'ρ_analytical':>14} {'ρ_simulated':>14} {'Error':>10}")
    print("-" * 60)
    
    for xi, theta in zip(xi_points, theta_values):
        r = xi * r_scale
        rho_analytical = rho_center_analytical * math.exp(-theta)
        
        # Find particles in this radial bin
        r_min = r - 0.03
        r_max = r + 0.03
        bin_particles = [p for p in cloud if r_min <= p['r'] < r_max]
        
        if bin_particles:
            rho_sim = sum(p['dens'] for p in bin_particles) / len(bin_particles)
            error = (rho_sim - rho_analytical) / rho_analytical * 100
            print(f"{r:<10.3f} {xi:>8.1f} {rho_analytical:>14.3f} {rho_sim:>14.3f} {error:>+9.1f}%")
        else:
            print(f"{r:<10.3f} {xi:>8.1f} {rho_analytical:>14.3f} {'N/A':>14} {'N/A':>10}")
    
    # ========================================
    # DIAGNOSTIC 7: Ghost Boundary Analysis
    # ========================================
    print()
    print("-" * 80)
    print("DIAGNOSTIC 7: GHOST BOUNDARY CONDITIONS")
    print("-" * 80)
    
    if ghost:
        cloud_r_max = max(p['r'] for p in cloud)
        ghost_r_min = min(p['r'] for p in ghost)
        ghost_r_max = max(p['r'] for p in ghost)
        gap = ghost_r_min - cloud_r_max
        
        # Ghost properties
        ghost_dens_avg = sum(p['dens'] for p in ghost) / len(ghost)
        ghost_pres_avg = sum(p['pres'] for p in ghost) / len(ghost)
        
        # Edge cloud properties  
        edge_cloud = [p for p in cloud if p['r'] > cloud_r_max - 0.05]
        edge_dens_avg = sum(p['dens'] for p in edge_cloud) / len(edge_cloud)
        edge_pres_avg = sum(p['pres'] for p in edge_cloud) / len(edge_cloud)
        
        # Expected external pressure for BE sphere at ξ_s = 6.0
        P_external_analytical = rho_center_analytical * math.exp(-3.0) * c_s_target**2  # θ(6) ≈ 3
        
        print(f"Cloud outer radius: {cloud_r_max:.4f} pc")
        print(f"Ghost inner radius: {ghost_r_min:.4f} pc")
        print(f"Ghost outer radius: {ghost_r_max:.4f} pc")
        print(f"Gap between cloud and ghost: {gap:.4f} pc")
        print()
        print(f"Edge cloud density: {edge_dens_avg:.4f} M☉/pc³")
        print(f"Ghost density: {ghost_dens_avg:.4f} M☉/pc³")
        print(f"Density ratio (ghost/edge): {ghost_dens_avg/edge_dens_avg:.4f}")
        print()
        print(f"Edge cloud pressure: {edge_pres_avg:.4f}")
        print(f"Ghost pressure: {ghost_pres_avg:.4f}")
        print(f"Expected external P (BE at ξ_s): {P_external_analytical:.4f}")
        print()
        
        if gap > 0.03:
            print("⚠️  LARGE GAP: Cloud edge particles don't feel ghost confinement!")
            print("    This causes massive pressure excess at edge.")
        
        if ghost_pres_avg < edge_pres_avg * 0.5:
            print("⚠️  GHOST PRESSURE TOO LOW: Not providing enough boundary confinement")
    
    # ========================================
    # ROOT CAUSE SUMMARY
    # ========================================
    print()
    print("=" * 80)
    print("ROOT CAUSE SUMMARY")
    print("=" * 80)
    print()
    
    # Collect evidence
    center_excess = None
    edge_excess = None
    
    for r_min, r_max, label in bins:
        bin_particles = [p for p in cloud if r_min <= p['r'] < r_max]
        if not bin_particles:
            continue
        
        a_grav_r_avg = sum(abs(p['a_grav_r']) for p in bin_particles) / len(bin_particles)
        a_pres_r_avg = sum(abs(p['a_pres_r']) for p in bin_particles) / len(bin_particles)
        
        if a_grav_r_avg > 1e-10:
            excess = (a_pres_r_avg / a_grav_r_avg - 1) * 100
        else:
            excess = 0
        
        if label == "Center":
            center_excess = excess
        elif label == "Edge":
            edge_excess = excess
    
    print("IDENTIFIED ISSUES:")
    print()
    
    if center_excess and center_excess > 3:
        print(f"1. ❌ UNIFORM PRESSURE EXCESS ({center_excess:.1f}% at center)")
        print("   Root cause: SPH pressure gradient has systematic positive bias")
        print("   This is likely due to:")
        print("   - SPH E0 (zeroth-order) error in density gradient regions")
        print("   - The pressure term Σ m_j (P_i/ρ_i² + P_j/ρ_j²) ∇W")
        print("     overestimates force when there's a density gradient")
        print()
    
    if edge_excess and edge_excess > 50:
        print(f"2. ❌ MASSIVE EDGE PRESSURE EXCESS ({edge_excess:.1f}%)")
        print("   Root cause: Gap between cloud and ghost envelope")
        print("   Edge particles only see cloud neighbors on one side")
        print("   This creates asymmetric pressure gradient → outward force")
        print()
    
    if ghost:
        gap_val = ghost_r_min - cloud_r_max
        if gap_val > 0.03:
            print(f"3. ❌ GHOST ENVELOPE GAP ({gap_val:.3f} pc)")
            print("   Ghost particles don't overlap with cloud edge")
            print("   No pressure confinement at boundary")
            print()
    
    print()
    print("RECOMMENDED FIXES:")
    print()
    print("Fix 1: REDUCE GHOST GAP")
    print("  - Move ghost envelope closer to cloud (inner radius at 0.73-0.75 pc)")
    print("  - Or use larger smoothing lengths at boundary")
    print()
    print("Fix 2: MATCH GHOST PRESSURE TO BE BOUNDARY")
    print("  - Ghost pressure should match P_external = ρ_edge * c_s²")
    print("  - Currently ghost pressure may be too low")
    print()
    print("Fix 3: USE GSPH INSTEAD OF SSPH")
    print("  - GSPH has better pressure gradient estimation")
    print("  - Uses Riemann solver for pairwise interactions")
    print()
    print("Fix 4: INCREASE PARTICLE RESOLUTION")
    print("  - More particles = smaller SPH E0 error")
    print("  - Error scales as O(h²) where h is smoothing length")
    print()


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = "simulations/astrophysics/imbh_cloud/results/phase2_hydrostatic/snapshot_0001.csv"
    
    analyze_be_equilibrium(filename)
