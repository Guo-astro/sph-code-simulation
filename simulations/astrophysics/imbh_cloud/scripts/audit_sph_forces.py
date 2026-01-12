#!/usr/bin/env python3
"""
First-Principles SPH Force Audit

This script analyzes the force imbalance in a Bonnor-Ebert sphere simulation
to identify the root cause of central density deficit.

NO NUMPY REQUIRED - pure Python implementation.
"""

import math
import csv
from collections import defaultdict


def load_csv(filename):
    """Load CSV file, skipping comment lines"""
    data = defaultdict(list)
    with open(filename, 'r') as f:
        headers = None
        for line in f:
            if not line.startswith('#'):
                headers = line.strip().split(',')
                break
        
        if not headers:
            raise ValueError("No header found in file")
        
        reader = csv.DictReader(f, fieldnames=headers)
        for row in reader:
            for col in headers:
                try:
                    data[col].append(float(row[col]))
                except (ValueError, KeyError):
                    data[col].append(0.0)
    return dict(data)


def solve_lane_emden(xi_s, n_points=5000):
    """
    Solve Lane-Emden equation for isothermal case (n=∞).
    
    d²ψ/dξ² + (2/ξ)dψ/dξ = exp(-ψ)
    where ρ/ρ_c = exp(-ψ)
    
    Returns: xi, rho_ratio, dpsi arrays
    """
    dxi = xi_s / n_points
    psi = [0.0] * (n_points + 1)
    dpsi = [0.0] * (n_points + 1)
    xi = [i * dxi for i in range(n_points + 1)]
    
    # Initial conditions: ψ(0) = 0, ψ'(0) = 0
    psi[0] = 0.0
    dpsi[0] = 0.0
    
    # RK4 integration
    for i in range(n_points):
        xi_i = max(xi[i], 1e-10)
        
        def f2(p, dp, x):
            if x > 1e-10:
                return math.exp(-p) - (2.0/x)*dp
            else:
                return math.exp(-p)  # At origin, term vanishes
        
        # RK4 step for dψ/dξ and ψ
        k1_p = dpsi[i]
        k1_dp = f2(psi[i], dpsi[i], xi_i)
        
        k2_p = dpsi[i] + 0.5*dxi*k1_dp
        k2_dp = f2(psi[i] + 0.5*dxi*k1_p, dpsi[i] + 0.5*dxi*k1_dp, xi_i + 0.5*dxi)
        
        k3_p = dpsi[i] + 0.5*dxi*k2_dp
        k3_dp = f2(psi[i] + 0.5*dxi*k2_p, dpsi[i] + 0.5*dxi*k2_dp, xi_i + 0.5*dxi)
        
        k4_p = dpsi[i] + dxi*k3_dp
        k4_dp = f2(psi[i] + dxi*k3_p, dpsi[i] + dxi*k3_dp, xi_i + dxi)
        
        psi[i+1] = psi[i] + dxi/6.0 * (k1_p + 2*k2_p + 2*k3_p + k4_p)
        dpsi[i+1] = dpsi[i] + dxi/6.0 * (k1_dp + 2*k2_dp + 2*k3_dp + k4_dp)
    
    # ρ/ρ_c = exp(-ψ)
    rho_ratio = [math.exp(-p) for p in psi]
    
    return xi, rho_ratio, dpsi


def interpolate_be(xi_val, xi, rho_ratio, dpsi):
    """Interpolate BE solution at given xi value"""
    if xi_val <= 0:
        return 1.0, 0.0
    if xi_val >= xi[-1]:
        return rho_ratio[-1], dpsi[-1]
    
    # Find index
    idx = int(xi_val / xi[-1] * (len(xi) - 1))
    idx = min(idx, len(xi) - 2)
    
    # Linear interpolation
    t = (xi_val - xi[idx]) / (xi[idx+1] - xi[idx]) if xi[idx+1] != xi[idx] else 0
    rho_r = rho_ratio[idx] * (1-t) + rho_ratio[idx+1] * t
    dpsi_r = dpsi[idx] * (1-t) + dpsi[idx+1] * t
    
    return rho_r, dpsi_r


def audit_sph_forces(filename):
    """
    Comprehensive audit of SPH forces.
    """
    
    print("=" * 80)
    print("FIRST-PRINCIPLES SPH FORCE AUDIT")
    print("=" * 80)
    print(f"\nAnalyzing: {filename}\n")
    
    data = load_csv(filename)
    n = len(data['pos_x'])
    
    # Physical parameters (from config)
    G = 0.00430091  # pc³ / (M☉ * Myr²)
    gamma = 1.0001  # Nearly isothermal
    c_s_target = 0.218  # pc/Myr
    xi_s = 6.0
    R_cloud = 0.75  # pc (approximate)
    
    # Solve BE equation once
    print("Solving Lane-Emden equation for isothermal sphere...")
    xi_be, rho_ratio_be, dpsi_be = solve_lane_emden(xi_s)
    print(f"  ξ_s = {xi_s}, ρ(ξ_s)/ρ_c = {rho_ratio_be[-1]:.6f}")
    print()
    
    # Extract particle data
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
    is_ghost = data.get('is_ghost', [0]*n)
    
    # Compute particle properties
    particles = []
    for i in range(n):
        r = math.sqrt(x[i]**2 + y[i]**2 + z[i]**2)
        
        # Radial unit vector
        if r > 1e-10:
            r_unit = (x[i]/r, y[i]/r, z[i]/r)
        else:
            r_unit = (0, 0, 0)
        
        # Pressure acceleration = total - gravity
        ax_pres = ax_total[i] - ax_grav[i]
        ay_pres = ay_total[i] - ay_grav[i]
        az_pres = az_total[i] - az_grav[i]
        
        # Radial components (positive = outward)
        a_grav_r = ax_grav[i]*r_unit[0] + ay_grav[i]*r_unit[1] + az_grav[i]*r_unit[2]
        a_pres_r = ax_pres*r_unit[0] + ay_pres*r_unit[1] + az_pres*r_unit[2]
        a_total_r = ax_total[i]*r_unit[0] + ay_total[i]*r_unit[1] + az_total[i]*r_unit[2]
        
        # Sound speed from internal energy
        c_s = math.sqrt((gamma - 1) * ene[i]) if ene[i] > 0 else 0
        
        particles.append({
            'r': r,
            'dens': dens[i],
            'pres': pres[i],
            'mass': mass[i],
            'sml': sml[i],
            'ene': ene[i],
            'gradh': gradh[i],
            'is_ghost': is_ghost[i],
            'a_grav_r': a_grav_r,
            'a_pres_r': a_pres_r,
            'a_total_r': a_total_r,
            'c_s': c_s,
        })
    
    # Separate cloud and ghost
    cloud = [p for p in particles if p['is_ghost'] == 0]
    ghost = [p for p in particles if p['is_ghost'] == 1]
    
    print(f"Total particles: {n}")
    print(f"Cloud particles: {len(cloud)}")
    print(f"Ghost particles: {len(ghost)}")
    
    # Find central density
    center_particles = [p for p in cloud if p['r'] < 0.1]
    if center_particles:
        rho_center = sum(p['dens'] for p in center_particles) / len(center_particles)
        c_s_actual = sum(p['c_s'] for p in center_particles) / len(center_particles)
    else:
        rho_center = max(p['dens'] for p in cloud)
        c_s_actual = c_s_target
    
    print(f"Central density: {rho_center:.4f} M☉/pc³")
    print(f"Sound speed: {c_s_actual:.4f} pc/Myr (target: {c_s_target})")
    print()
    
    # ========================================
    # ANALYSIS 1: Compare with Analytical BE
    # ========================================
    print("-" * 80)
    print("ANALYSIS 1: COMPARISON WITH ANALYTICAL BONNOR-EBERT SOLUTION")
    print("-" * 80)
    print()
    
    print(f"{'Radius':<10} {'|a_grav|_code':>15} {'|a_grav|_BE':>15} {'|a_pres|_code':>15} {'|a_pres|_BE':>15} {'Ratio_p':>10}")
    print("-" * 95)
    
    radial_bins = [(0.0, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.4), (0.4, 0.5), (0.5, 0.6), (0.6, 0.7)]
    
    analysis_data = []
    
    for r_min, r_max in radial_bins:
        bin_particles = [p for p in cloud if r_min <= p['r'] < r_max]
        if not bin_particles:
            continue
        
        r_avg = sum(p['r'] for p in bin_particles) / len(bin_particles)
        
        # Code forces (magnitudes)
        a_grav_r_avg = sum(abs(p['a_grav_r']) for p in bin_particles) / len(bin_particles)
        a_pres_r_avg = sum(abs(p['a_pres_r']) for p in bin_particles) / len(bin_particles)
        
        # Analytical BE forces at this radius
        xi_r = r_avg * xi_s / R_cloud
        if xi_r < xi_s and r_avg > 0.01:
            rho_r, dpsi_r = interpolate_be(xi_r, xi_be, rho_ratio_be, dpsi_be)
            
            # a_pres_analytical = c_s² * d(ln ρ)/dr = c_s² * (-dψ/dξ) * (ξ_s/R)
            a_pres_be = abs(c_s_target**2 * (-dpsi_r) * (xi_s / R_cloud))
            
            # M(r) from BE solution: ξ² dψ/dξ
            M_enclosed = 4.0 * math.pi * rho_center * (R_cloud/xi_s)**3 * (xi_r**2 * abs(dpsi_r))
            a_grav_be = G * M_enclosed / (r_avg**2)
        else:
            a_pres_be = 0
            a_grav_be = 0
        
        ratio_p = a_pres_r_avg / a_pres_be if a_pres_be > 1e-15 else 0
        
        print(f"{r_min:.1f}-{r_max:.1f}      {a_grav_r_avg:>15.6e} {a_grav_be:>15.6e} {a_pres_r_avg:>15.6e} {a_pres_be:>15.6e} {ratio_p:>10.3f}")
        
        analysis_data.append({
            'r_avg': r_avg,
            'a_grav_code': a_grav_r_avg,
            'a_grav_be': a_grav_be,
            'a_pres_code': a_pres_r_avg,
            'a_pres_be': a_pres_be,
            'ratio_pres': ratio_p,
            'n_particles': len(bin_particles),
        })
    
    print()
    
    # ========================================
    # ANALYSIS 2: Force Balance by Radius
    # ========================================
    print("-" * 80)
    print("ANALYSIS 2: FORCE BALANCE (PRESSURE/GRAVITY RATIO)")
    print("-" * 80)
    print("For equilibrium: |a_pres| / |a_grav| = 1.0")
    print("Values > 1 indicate net outward force (expansion)")
    print("Values < 1 indicate net inward force (collapse)")
    print()
    
    print(f"{'Radius':<10} {'|a_grav|':>12} {'|a_pres|':>12} {'P/G Ratio':>12} {'Net_a_r':>12} {'Direction':<10}")
    print("-" * 80)
    
    force_balance_data = []
    
    for r_min, r_max in radial_bins:
        bin_particles = [p for p in cloud if r_min <= p['r'] < r_max]
        if not bin_particles:
            continue
        
        a_grav_r_avg = sum(abs(p['a_grav_r']) for p in bin_particles) / len(bin_particles)
        a_pres_r_avg = sum(abs(p['a_pres_r']) for p in bin_particles) / len(bin_particles)
        a_total_r_avg = sum(p['a_total_r'] for p in bin_particles) / len(bin_particles)
        
        pg_ratio = a_pres_r_avg / a_grav_r_avg if a_grav_r_avg > 1e-15 else 0
        direction = "OUTWARD" if a_total_r_avg > 0 else "INWARD"
        
        print(f"{r_min:.1f}-{r_max:.1f}      {a_grav_r_avg:>12.6e} {a_pres_r_avg:>12.6e} {pg_ratio:>12.4f} {a_total_r_avg:>+12.6e} {direction:<10}")
        
        force_balance_data.append({
            'r_bin': f"{r_min:.1f}-{r_max:.1f}",
            'pg_ratio': pg_ratio,
            'excess_pct': (pg_ratio - 1) * 100,
        })
    
    print()
    
    # ========================================
    # ANALYSIS 3: Density Profile Comparison
    # ========================================
    print("-" * 80)
    print("ANALYSIS 3: DENSITY PROFILE VS ANALYTICAL BE")
    print("-" * 80)
    print()
    
    print(f"{'Radius':<10} {'ρ_SPH':>12} {'ρ_BE':>12} {'Ratio':>10} {'Error%':>10}")
    print("-" * 60)
    
    for r_min, r_max in radial_bins:
        bin_particles = [p for p in cloud if r_min <= p['r'] < r_max]
        if not bin_particles:
            continue
        
        r_avg = sum(p['r'] for p in bin_particles) / len(bin_particles)
        rho_sph = sum(p['dens'] for p in bin_particles) / len(bin_particles)
        
        # Analytical BE density
        xi_r = r_avg * xi_s / R_cloud
        if xi_r < xi_s:
            rho_r, _ = interpolate_be(xi_r, xi_be, rho_ratio_be, dpsi_be)
            rho_be = rho_center * rho_r
        else:
            rho_be = 0
        
        ratio = rho_sph / rho_be if rho_be > 0 else 0
        error_pct = (ratio - 1) * 100
        
        print(f"{r_min:.1f}-{r_max:.1f}      {rho_sph:>12.4f} {rho_be:>12.4f} {ratio:>10.4f} {error_pct:>+10.1f}%")
    
    print()
    
    # ========================================
    # ANALYSIS 4: Neighbor Count Analysis
    # ========================================
    print("-" * 80)
    print("ANALYSIS 4: NEIGHBOR COUNT BY RADIUS")
    print("-" * 80)
    print("Low neighbor count at edge indicates kernel truncation")
    print()
    
    neighbor_data = data.get('neighbor', [0]*n)
    
    print(f"{'Radius':<10} {'Mean N':>10} {'Min N':>10} {'Max N':>10} {'Std':>10}")
    print("-" * 60)
    
    for r_min, r_max in radial_bins:
        bin_indices = [i for i in range(n) if is_ghost[i] == 0 and r_min <= particles[i]['r'] < r_max]
        if not bin_indices:
            continue
        
        neighbors = [neighbor_data[i] for i in bin_indices]
        mean_n = sum(neighbors) / len(neighbors)
        min_n = min(neighbors)
        max_n = max(neighbors)
        
        # Standard deviation
        variance = sum((nb - mean_n)**2 for nb in neighbors) / len(neighbors)
        std_n = math.sqrt(variance)
        
        print(f"{r_min:.1f}-{r_max:.1f}      {mean_n:>10.1f} {min_n:>10.0f} {max_n:>10.0f} {std_n:>10.1f}")
    
    print()
    
    # ========================================
    # ANALYSIS 5: Ghost Particle Analysis
    # ========================================
    print("-" * 80)
    print("ANALYSIS 5: GHOST BOUNDARY ANALYSIS")
    print("-" * 80)
    print()
    
    if ghost:
        cloud_r_max = max(p['r'] for p in cloud)
        ghost_r_min = min(p['r'] for p in ghost)
        ghost_r_max = max(p['r'] for p in ghost)
        gap = ghost_r_min - cloud_r_max
        
        # Edge cloud smoothing length
        edge_cloud = [p for p in cloud if p['r'] > cloud_r_max - 0.05]
        edge_sml = sum(p['sml'] for p in edge_cloud) / len(edge_cloud) if edge_cloud else 0
        
        print(f"Cloud outer radius:  {cloud_r_max:.4f} pc")
        print(f"Ghost inner radius:  {ghost_r_min:.4f} pc")
        print(f"Gap:                 {gap:.4f} pc")
        print(f"Edge smoothing len:  {edge_sml:.4f} pc")
        print()
        
        if gap > edge_sml:
            print("⚠️  CRITICAL: Gap > smoothing length!")
            print("   Edge particles don't feel ghost confinement")
            print("   This causes MASSIVE outward pressure force at edge")
        elif gap > 0:
            print(f"⚠️  WARNING: Positive gap ({gap:.4f} pc)")
            print("   Ghost envelope should overlap with cloud edge")
        else:
            print("✓  Ghost envelope overlaps with cloud (good)")
        
        # Ghost density vs edge cloud density
        edge_dens = sum(p['dens'] for p in edge_cloud) / len(edge_cloud) if edge_cloud else 0
        ghost_dens = sum(p['dens'] for p in ghost) / len(ghost)
        
        print()
        print(f"Edge cloud density:  {edge_dens:.4f} M☉/pc³")
        print(f"Ghost density:       {ghost_dens:.4f} M☉/pc³")
        print(f"Ratio (ghost/edge):  {ghost_dens/edge_dens:.4f}" if edge_dens > 0 else "")
    else:
        print("No ghost particles found!")
    
    print()
    
    # ========================================
    # SUMMARY: ROOT CAUSE IDENTIFICATION
    # ========================================
    print("=" * 80)
    print("ROOT CAUSE IDENTIFICATION SUMMARY")
    print("=" * 80)
    print()
    
    # Check for uniform pressure excess
    interior_data = [d for d in force_balance_data if float(d['r_bin'].split('-')[1]) <= 0.5]
    edge_data_fb = [d for d in force_balance_data if float(d['r_bin'].split('-')[0]) >= 0.5]
    
    if interior_data:
        avg_interior_excess = sum(d['excess_pct'] for d in interior_data) / len(interior_data)
        
        print(f"INTERIOR REGION (r < 0.5 pc):")
        print(f"  Average pressure excess: {avg_interior_excess:+.1f}%")
        
        if avg_interior_excess > 3:
            print(f"  ⚠️  DIAGNOSIS: Systematic SPH pressure overestimate")
            print(f"     This is the SPH 'E0 error' - inherent in symmetric SPH formulation")
            print(f"     The pressure gradient term has O(h²) truncation error")
    
    print()
    
    if edge_data_fb:
        avg_edge_excess = sum(d['excess_pct'] for d in edge_data_fb) / len(edge_data_fb)
        
        print(f"EDGE REGION (r > 0.5 pc):")
        print(f"  Average pressure excess: {avg_edge_excess:+.1f}%")
        
        if avg_edge_excess > 50:
            print(f"  ⚠️  DIAGNOSIS: Massive edge pressure excess")
            print(f"     Cause: Ghost boundary gap - edge particles have asymmetric neighbors")
    
    print()
    
    # Sound speed check
    c_s_error = (c_s_actual / c_s_target - 1) * 100
    print(f"EQUATION OF STATE:")
    print(f"  Sound speed error: {c_s_error:+.1f}%")
    
    if abs(c_s_error) > 2:
        u_needed = c_s_target**2 / (gamma - 1)
        u_actual = center_particles[0]['ene'] if center_particles else 0
        print(f"  ⚠️  Internal energy mismatch!")
        print(f"     Current u = {u_actual:.2f}")
        print(f"     Needed u = {u_needed:.2f} for c_s = {c_s_target}")
    
    print()
    print("=" * 80)
    print("RECOMMENDED FIXES")
    print("=" * 80)
    print()
    print("1. FIX GHOST BOUNDARY GAP (most impactful for edge):")
    print("   - Move ghost inner radius to ≤ cloud outer radius")
    print("   - Ideally: ghost_r_min = cloud_r_max - 2*h_edge")
    print()
    print("2. FIX INTERNAL ENERGY (for correct EOS):")
    print(f"   - Set u = c_s²/(γ-1) = {c_s_target**2 / (gamma - 1):.2f}")
    print()
    print("3. FOR SPH E0 ERROR:")
    print("   - Use GSPH instead of SSPH (Riemann solver has better accuracy)")
    print("   - Or increase resolution (error ∝ h²)")
    print()
    print("4. VERIFY INITIAL BE PROFILE:")
    print("   - Check sample generator against analytical BE solution")
    print()


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = "results/phase2_hydrostatic/snapshot_0001.csv"
    
    audit_sph_forces(filename)
