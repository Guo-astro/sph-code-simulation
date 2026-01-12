#!/usr/bin/env python3
"""
Plot radial density profile and force balance for Bonnor-Ebert sphere.

Generates two plots:
1. Radial density profile: SPH vs analytical BE
2. Force balance profile: |a_pres|, |a_grav|, and P/G ratio
"""

import math
import csv
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


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
    """Solve Lane-Emden equation for isothermal sphere"""
    dxi = xi_s / n_points
    psi = [0.0] * (n_points + 1)
    dpsi = [0.0] * (n_points + 1)
    xi = [i * dxi for i in range(n_points + 1)]
    
    psi[0] = 0.0
    dpsi[0] = 0.0
    
    for i in range(n_points):
        xi_i = max(xi[i], 1e-10)
        
        def f2(p, dp, x):
            if x > 1e-10:
                return math.exp(-p) - (2.0/x)*dp
            else:
                return math.exp(-p)
        
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
    
    rho_ratio = [math.exp(-p) for p in psi]
    return xi, rho_ratio, dpsi


def plot_profiles(filename, output_file=None):
    """Generate density and force balance plots"""
    
    print(f"Loading data from: {filename}")
    data = load_csv(filename)
    n = len(data['pos_x'])
    
    # Physical parameters
    G = 0.00430091
    gamma = 1.0001
    c_s_target = 0.218
    xi_s = 6.0
    R_cloud = 0.75
    
    # Solve BE equation
    xi_be, rho_ratio_be, dpsi_be = solve_lane_emden(xi_s)
    
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
    is_ghost = data.get('is_ghost', [0]*n)
    
    # Compute particle properties
    particles = []
    for i in range(n):
        if is_ghost[i] == 1:
            continue
            
        r = math.sqrt(x[i]**2 + y[i]**2 + z[i]**2)
        if r < 1e-10:
            continue
        
        r_unit = (x[i]/r, y[i]/r, z[i]/r)
        
        ax_pres = ax_total[i] - ax_grav[i]
        ay_pres = ay_total[i] - ay_grav[i]
        az_pres = az_total[i] - az_grav[i]
        
        a_grav_r = ax_grav[i]*r_unit[0] + ay_grav[i]*r_unit[1] + az_grav[i]*r_unit[2]
        a_pres_r = ax_pres*r_unit[0] + ay_pres*r_unit[1] + az_pres*r_unit[2]
        
        particles.append({
            'r': r,
            'dens': dens[i],
            'a_grav_r': a_grav_r,
            'a_pres_r': a_pres_r,
        })
    
    # Find central density
    center_p = [p for p in particles if p['r'] < 0.1]
    rho_center = sum(p['dens'] for p in center_p) / len(center_p) if center_p else 50.0
    
    # Bin data
    n_bins = 30
    r_max = 0.75
    bin_edges = [i * r_max / n_bins for i in range(n_bins + 1)]
    
    bin_r = []
    bin_dens_sph = []
    bin_dens_be = []
    bin_a_grav = []
    bin_a_pres = []
    bin_pg_ratio = []
    bin_count = []
    
    for i in range(n_bins):
        r_min_bin = bin_edges[i]
        r_max_bin = bin_edges[i + 1]
        r_mid = (r_min_bin + r_max_bin) / 2
        
        bin_particles = [p for p in particles if r_min_bin <= p['r'] < r_max_bin]
        if len(bin_particles) < 3:
            continue
        
        # SPH averages
        dens_avg = sum(p['dens'] for p in bin_particles) / len(bin_particles)
        a_grav_avg = sum(abs(p['a_grav_r']) for p in bin_particles) / len(bin_particles)
        a_pres_avg = sum(abs(p['a_pres_r']) for p in bin_particles) / len(bin_particles)
        
        # Analytical BE
        xi_r = r_mid * xi_s / R_cloud
        if xi_r < xi_s:
            idx = int(xi_r / xi_s * (len(xi_be) - 1))
            idx = min(idx, len(xi_be) - 2)
            t = (xi_r - xi_be[idx]) / (xi_be[idx+1] - xi_be[idx])
            rho_r = rho_ratio_be[idx] * (1-t) + rho_ratio_be[idx+1] * t
            dens_be = rho_center * rho_r
        else:
            dens_be = 0
        
        pg_ratio = a_pres_avg / a_grav_avg if a_grav_avg > 1e-15 else 1.0
        
        bin_r.append(r_mid)
        bin_dens_sph.append(dens_avg)
        bin_dens_be.append(dens_be)
        bin_a_grav.append(a_grav_avg)
        bin_a_pres.append(a_pres_avg)
        bin_pg_ratio.append(pg_ratio)
        bin_count.append(len(bin_particles))
    
    # Create figure
    fig = plt.figure(figsize=(14, 10))
    gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1], hspace=0.3, wspace=0.3)
    
    # Plot 1: Density Profile
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.scatter([p['r'] for p in particles[::5]], [p['dens'] for p in particles[::5]], 
                s=1, alpha=0.3, c='blue', label='SPH particles')
    ax1.plot(bin_r, bin_dens_sph, 'b-', linewidth=2, label='SPH (binned)')
    ax1.plot(bin_r, bin_dens_be, 'r--', linewidth=2, label='Analytical BE')
    ax1.set_xlabel('Radius [pc]', fontsize=12)
    ax1.set_ylabel('Density [M☉/pc³]', fontsize=12)
    ax1.set_title('Radial Density Profile', fontsize=14)
    ax1.legend(loc='upper right')
    ax1.set_xlim(0, 0.8)
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Density Ratio
    ax2 = fig.add_subplot(gs[0, 1])
    ratio = [s/b if b > 0 else 1 for s, b in zip(bin_dens_sph, bin_dens_be)]
    ax2.plot(bin_r, ratio, 'g-', linewidth=2, marker='o', markersize=4)
    ax2.axhline(y=1.0, color='k', linestyle='--', alpha=0.5, label='Perfect match')
    ax2.fill_between(bin_r, [0.95]*len(bin_r), [1.05]*len(bin_r), alpha=0.2, color='green', label='±5% band')
    ax2.set_xlabel('Radius [pc]', fontsize=12)
    ax2.set_ylabel('ρ_SPH / ρ_BE', fontsize=12)
    ax2.set_title('Density: SPH vs Analytical', fontsize=14)
    ax2.set_xlim(0, 0.8)
    ax2.set_ylim(0.9, 1.15)
    ax2.legend(loc='upper right')
    ax2.grid(True, alpha=0.3)
    
    # Add text annotation for average error
    avg_ratio = sum(ratio) / len(ratio)
    ax2.text(0.05, 0.95, f'Average: {avg_ratio:.3f}\n(+{(avg_ratio-1)*100:.1f}%)', 
             transform=ax2.transAxes, fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Plot 3: Force Magnitudes
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.plot(bin_r, bin_a_grav, 'r-', linewidth=2, label='|a_gravity| (inward)')
    ax3.plot(bin_r, bin_a_pres, 'b-', linewidth=2, label='|a_pressure| (outward)')
    ax3.set_xlabel('Radius [pc]', fontsize=12)
    ax3.set_ylabel('Acceleration [pc/Myr²]', fontsize=12)
    ax3.set_title('Force Balance: Pressure vs Gravity', fontsize=14)
    ax3.legend(loc='upper right')
    ax3.set_xlim(0, 0.8)
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: P/G Ratio
    ax4 = fig.add_subplot(gs[1, 1])
    colors = ['red' if r > 1.03 else 'orange' if r > 1.01 else 'green' for r in bin_pg_ratio]
    ax4.bar(bin_r, bin_pg_ratio, width=0.02, color=colors, alpha=0.7, edgecolor='black')
    ax4.axhline(y=1.0, color='k', linestyle='-', linewidth=2, label='Equilibrium (P/G = 1)')
    ax4.axhline(y=1.05, color='r', linestyle='--', alpha=0.7, label='+5% excess')
    ax4.set_xlabel('Radius [pc]', fontsize=12)
    ax4.set_ylabel('|a_pres| / |a_grav|', fontsize=12)
    ax4.set_title('Pressure/Gravity Ratio (1.0 = equilibrium)', fontsize=14)
    ax4.set_xlim(0, 0.8)
    ax4.set_ylim(0.95, 1.15)
    ax4.legend(loc='upper right')
    ax4.grid(True, alpha=0.3)
    
    # Add text annotation
    avg_pg = sum(bin_pg_ratio) / len(bin_pg_ratio)
    ax4.text(0.05, 0.95, f'Average P/G: {avg_pg:.3f}\nExcess: +{(avg_pg-1)*100:.1f}%', 
             transform=ax4.transAxes, fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Main title
    fig.suptitle('Bonnor-Ebert Sphere: Hydrostatic Equilibrium Analysis\n(GSPH with grad-h correction)', 
                 fontsize=16, fontweight='bold')
    
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    
    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Plot saved to: {output_file}")
    
    plt.show()
    
    return fig


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = "results/phase2_hydrostatic/snapshot_0001.csv"
    
    output_file = "results/force_balance_analysis.png"
    if len(sys.argv) > 2:
        output_file = sys.argv[2]
    
    plot_profiles(filename, output_file)
