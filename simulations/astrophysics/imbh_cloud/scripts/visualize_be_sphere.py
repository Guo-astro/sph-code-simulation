#!/usr/bin/env python3
"""
Bonnor-Ebert Sphere Visualization

Generates:
1. Radial density profile comparison (SPH vs analytical BE)
2. Force balance profile (pressure vs gravity)
3. GIF animation of time evolution
"""

import math
import csv
import os
import glob
from collections import defaultdict

# Check for matplotlib
try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("ERROR: matplotlib not found. Install with: pip3 install matplotlib")
    exit(1)


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


def get_time_from_csv(filename):
    """Extract time from CSV header"""
    with open(filename, 'r') as f:
        for line in f:
            if '# Time (code):' in line:
                try:
                    return float(line.split(':')[1].strip())
                except:
                    pass
    return 0.0


def solve_lane_emden(xi_s, n_points=5000):
    """Solve Lane-Emden equation for isothermal case"""
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


def compute_profiles(filename, xi_s=6.0, R_cloud=0.75):
    """Compute radial density and force profiles from CSV"""
    data = load_csv(filename)
    n = len(data['pos_x'])
    
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
    
    # Compute for each particle
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
    
    # Bin particles
    r_bins = [0.05 + i*0.05 for i in range(15)]
    
    density_profile = {'r': [], 'rho': [], 'rho_std': []}
    force_profile = {'r': [], 'a_grav': [], 'a_pres': [], 'ratio': []}
    
    for r_center in r_bins:
        r_min, r_max = r_center - 0.025, r_center + 0.025
        bin_particles = [p for p in particles if r_min <= p['r'] < r_max]
        
        if len(bin_particles) < 5:
            continue
        
        rho_vals = [p['dens'] for p in bin_particles]
        rho_mean = sum(rho_vals) / len(rho_vals)
        rho_std = math.sqrt(sum((r - rho_mean)**2 for r in rho_vals) / len(rho_vals))
        
        a_grav_mean = sum(abs(p['a_grav_r']) for p in bin_particles) / len(bin_particles)
        a_pres_mean = sum(abs(p['a_pres_r']) for p in bin_particles) / len(bin_particles)
        ratio = a_pres_mean / a_grav_mean if a_grav_mean > 1e-15 else 1.0
        
        density_profile['r'].append(r_center)
        density_profile['rho'].append(rho_mean)
        density_profile['rho_std'].append(rho_std)
        
        force_profile['r'].append(r_center)
        force_profile['a_grav'].append(a_grav_mean)
        force_profile['a_pres'].append(a_pres_mean)
        force_profile['ratio'].append(ratio)
    
    # Get central density for BE scaling
    center_particles = [p for p in particles if p['r'] < 0.1]
    rho_center = sum(p['dens'] for p in center_particles) / len(center_particles) if center_particles else 50.0
    
    return density_profile, force_profile, rho_center


def plot_density_profile(filename, output_file='density_profile.png'):
    """Plot radial density profile comparison"""
    print(f"Generating density profile plot...")
    
    density_profile, _, rho_center = compute_profiles(filename)
    
    # Analytical BE solution
    xi_s = 6.0
    R_cloud = 0.75
    xi_be, rho_ratio_be, _ = solve_lane_emden(xi_s)
    r_be = [xi * R_cloud / xi_s for xi in xi_be]
    rho_be = [rho_center * rho_r for rho_r in rho_ratio_be]
    
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Plot SPH data
    ax.errorbar(density_profile['r'], density_profile['rho'], 
                yerr=density_profile['rho_std'],
                fmt='o', markersize=8, capsize=3, 
                color='blue', label='SPH Simulation', alpha=0.8)
    
    # Plot analytical BE
    ax.plot(r_be, rho_be, 'r-', linewidth=2, label='Analytical BE')
    
    ax.set_xlabel('Radius (pc)', fontsize=14)
    ax.set_ylabel('Density (M☉/pc³)', fontsize=14)
    ax.set_title('Radial Density Profile: SPH vs Analytical Bonnor-Ebert', fontsize=16)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 0.8)
    ax.set_ylim(0, max(density_profile['rho']) * 1.2)
    
    # Add text box with parameters
    time = get_time_from_csv(filename)
    textstr = f't = {time:.1f} Myr\nρ_center = {rho_center:.1f} M☉/pc³\nξ_s = 6.0'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.95, 0.95, textstr, transform=ax.transAxes, fontsize=11,
            verticalalignment='top', horizontalalignment='right', bbox=props)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    plt.close()
    print(f"  Saved: {output_file}")


def plot_force_profile(filename, output_file='force_profile.png'):
    """Plot force balance profile"""
    print(f"Generating force profile plot...")
    
    _, force_profile, _ = compute_profiles(filename)
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True)
    
    # Top: Force magnitudes
    ax1.plot(force_profile['r'], force_profile['a_grav'], 'b-o', 
             markersize=6, linewidth=2, label='|a_gravity|')
    ax1.plot(force_profile['r'], force_profile['a_pres'], 'r-s', 
             markersize=6, linewidth=2, label='|a_pressure|')
    
    ax1.set_ylabel('Acceleration Magnitude', fontsize=14)
    ax1.set_title('Force Balance Analysis', fontsize=16)
    ax1.legend(fontsize=12)
    ax1.grid(True, alpha=0.3)
    ax1.set_yscale('log')
    
    # Bottom: Pressure/Gravity ratio
    ax2.plot(force_profile['r'], force_profile['ratio'], 'g-^', 
             markersize=8, linewidth=2, label='|a_pres|/|a_grav|')
    ax2.axhline(y=1.0, color='k', linestyle='--', linewidth=1.5, label='Equilibrium (ratio=1)')
    
    ax2.set_xlabel('Radius (pc)', fontsize=14)
    ax2.set_ylabel('Pressure/Gravity Ratio', fontsize=14)
    ax2.legend(fontsize=12)
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(0.9, 1.2)
    
    # Add annotation for excess
    avg_ratio = sum(force_profile['ratio']) / len(force_profile['ratio'])
    excess_pct = (avg_ratio - 1) * 100
    ax2.text(0.5, 0.85, f'Average excess: {excess_pct:+.1f}%', 
             transform=ax2.transAxes, fontsize=12,
             bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    plt.close()
    print(f"  Saved: {output_file}")


def create_animation(result_dir, output_file='be_evolution.gif'):
    """Create GIF animation of density evolution"""
    print(f"Creating animation...")
    
    # Find all snapshots
    snapshots = sorted(glob.glob(os.path.join(result_dir, 'snapshot_*.csv')))
    if not snapshots:
        print(f"  No snapshots found in {result_dir}")
        return
    
    print(f"  Found {len(snapshots)} snapshots")
    
    # Precompute analytical BE
    xi_s = 6.0
    R_cloud = 0.75
    xi_be, rho_ratio_be, _ = solve_lane_emden(xi_s)
    r_be = [xi * R_cloud / xi_s for xi in xi_be]
    
    # Setup figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Initialize plots
    line_sph, = ax1.plot([], [], 'bo', markersize=4, alpha=0.6, label='SPH')
    line_be, = ax1.plot([], [], 'r-', linewidth=2, label='Analytical BE')
    line_ratio, = ax2.plot([], [], 'g-^', markersize=6, linewidth=2)
    line_eq = ax2.axhline(y=1.0, color='k', linestyle='--', linewidth=1.5)
    
    ax1.set_xlim(0, 0.8)
    ax1.set_ylim(0, 60)
    ax1.set_xlabel('Radius (pc)', fontsize=12)
    ax1.set_ylabel('Density (M☉/pc³)', fontsize=12)
    ax1.set_title('Density Profile', fontsize=14)
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    ax2.set_xlim(0, 0.8)
    ax2.set_ylim(0.9, 1.2)
    ax2.set_xlabel('Radius (pc)', fontsize=12)
    ax2.set_ylabel('Pressure/Gravity Ratio', fontsize=12)
    ax2.set_title('Force Balance', fontsize=14)
    ax2.grid(True, alpha=0.3)
    
    time_text = ax1.text(0.02, 0.95, '', transform=ax1.transAxes, fontsize=12,
                         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    excess_text = ax2.text(0.02, 0.95, '', transform=ax2.transAxes, fontsize=12,
                           bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.8))
    
    def init():
        line_sph.set_data([], [])
        line_be.set_data([], [])
        line_ratio.set_data([], [])
        time_text.set_text('')
        excess_text.set_text('')
        return line_sph, line_be, line_ratio, time_text, excess_text
    
    def animate(frame_idx):
        if frame_idx >= len(snapshots):
            return line_sph, line_be, line_ratio, time_text, excess_text
        
        filename = snapshots[frame_idx]
        time = get_time_from_csv(filename)
        
        try:
            density_profile, force_profile, rho_center = compute_profiles(filename)
        except Exception as e:
            print(f"  Error processing {filename}: {e}")
            return line_sph, line_be, line_ratio, time_text, excess_text
        
        # Update density plot
        line_sph.set_data(density_profile['r'], density_profile['rho'])
        rho_be = [rho_center * rho_ratio_be[i] for i in range(len(r_be))]
        line_be.set_data(r_be, rho_be)
        
        # Update force ratio plot
        line_ratio.set_data(force_profile['r'], force_profile['ratio'])
        
        # Update text
        time_text.set_text(f't = {time:.1f} Myr\nρ_c = {rho_center:.1f}')
        
        if force_profile['ratio']:
            avg_ratio = sum(force_profile['ratio']) / len(force_profile['ratio'])
            excess_pct = (avg_ratio - 1) * 100
            excess_text.set_text(f'Excess: {excess_pct:+.1f}%')
        
        print(f"  Frame {frame_idx+1}/{len(snapshots)}: t={time:.1f} Myr")
        
        return line_sph, line_be, line_ratio, time_text, excess_text
    
    anim = FuncAnimation(fig, animate, init_func=init, 
                         frames=len(snapshots), interval=200, blit=True)
    
    # Save as GIF
    writer = PillowWriter(fps=5)
    anim.save(output_file, writer=writer)
    plt.close()
    print(f"  Saved: {output_file}")


def main():
    import sys
    
    if len(sys.argv) < 2:
        # Default paths
        result_dir = "results/phase2_hydrostatic"
        snapshot = "results/phase2_hydrostatic/snapshot_0001.csv"
    else:
        arg = sys.argv[1]
        if os.path.isdir(arg):
            result_dir = arg
            snapshots = sorted(glob.glob(os.path.join(arg, 'snapshot_*.csv')))
            snapshot = snapshots[-1] if snapshots else None
        else:
            snapshot = arg
            result_dir = os.path.dirname(arg)
    
    if snapshot and os.path.exists(snapshot):
        plot_density_profile(snapshot, 'density_profile.png')
        plot_force_profile(snapshot, 'force_profile.png')
    
    if result_dir and os.path.isdir(result_dir):
        create_animation(result_dir, 'be_evolution.gif')
    
    print("\nDone!")


if __name__ == "__main__":
    main()
