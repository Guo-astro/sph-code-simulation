#!/usr/bin/env python3
"""
Radial Profile Animation for Bonnor-Ebert Sphere Relaxation
============================================================

Creates animated GIF showing:
- Top panel: Radial density profile vs analytical BE solution
- Bottom panel: Force balance (SPH force vs gravity)

Author: SPH Simulation Team
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import argparse
from pathlib import Path
import glob

# Dark theme colors
DARK_BG = '#0d0d12'
DARK_PANEL = '#15151f'
TEXT_COLOR = '#f0f0f0'
GRID_COLOR = '#3a3a4a'

def solve_lane_emden_isothermal(xi_max=10.0, n_points=1000):
    """
    Solve the isothermal Lane-Emden equation:
    d²ψ/dξ² + (2/ξ) dψ/dξ = exp(-ψ)
    
    with boundary conditions: ψ(0) = 0, ψ'(0) = 0
    
    Returns:
        xi: dimensionless radius array
        psi: dimensionless potential array
        dpsi: derivative dψ/dξ array
    """
    from scipy.integrate import odeint
    
    def derivatives(y, xi):
        psi, dpsi = y
        if xi < 1e-10:
            # At center: d²ψ/dξ² = exp(-ψ)/3 (L'Hopital's rule)
            d2psi = np.exp(-psi) / 3.0
        else:
            d2psi = np.exp(-psi) - 2.0 * dpsi / xi
        return [dpsi, d2psi]
    
    # Initial conditions at center
    y0 = [0.0, 0.0]
    
    # Solve from small xi to avoid singularity
    xi = np.linspace(1e-6, xi_max, n_points)
    solution = odeint(derivatives, y0, xi)
    
    # Prepend center point
    xi = np.concatenate([[0.0], xi])
    psi = np.concatenate([[0.0], solution[:, 0]])
    dpsi = np.concatenate([[0.0], solution[:, 1]])
    
    return xi, psi, dpsi


def load_snapshot(filepath):
    """Load snapshot CSV and return relevant data."""
    import pandas as pd
    
    df = pd.read_csv(filepath, comment='#')
    
    # Extract positions
    x = df['pos_x'].values
    y = df['pos_y'].values
    z = df['pos_z'].values
    r = np.sqrt(x**2 + y**2 + z**2)
    
    # Extract density and acceleration
    dens = df['dens'].values
    
    # Check for ghost flag
    if 'is_ghost' in df.columns:
        is_ghost = df['is_ghost'].values.astype(bool)
    else:
        is_ghost = np.zeros(len(df), dtype=bool)
    
    # Extract accelerations
    ax = df['acc_x'].values if 'acc_x' in df.columns else np.zeros_like(x)
    ay = df['acc_y'].values if 'acc_y' in df.columns else np.zeros_like(x)
    az = df['acc_z'].values if 'acc_z' in df.columns else np.zeros_like(x)
    a_mag = np.sqrt(ax**2 + ay**2 + az**2)
    
    # Radial acceleration (positive = outward)
    with np.errstate(divide='ignore', invalid='ignore'):
        a_r = (ax*x + ay*y + az*z) / r
        a_r = np.nan_to_num(a_r, nan=0.0)
    
    # Gravity acceleration
    gx = df['grav_acc_x'].values if 'grav_acc_x' in df.columns else np.zeros_like(x)
    gy = df['grav_acc_y'].values if 'grav_acc_y' in df.columns else np.zeros_like(x)
    gz = df['grav_acc_z'].values if 'grav_acc_z' in df.columns else np.zeros_like(x)
    g_mag = np.sqrt(gx**2 + gy**2 + gz**2)
    
    # Radial gravity (should be negative = inward)
    with np.errstate(divide='ignore', invalid='ignore'):
        g_r = (gx*x + gy*y + gz*z) / r
        g_r = np.nan_to_num(g_r, nan=0.0)
    
    return {
        'r': r,
        'dens': dens,
        'is_ghost': is_ghost,
        'a_r': a_r,
        'a_mag': a_mag,
        'g_r': g_r,
        'g_mag': g_mag,
    }


def get_time_from_snapshot(filepath):
    """Extract time from snapshot header."""
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('# Time (code):'):
                return float(line.split(':')[1].strip())
    return 0.0


def create_animation(results_dir, output_path, xi_s=6.0, R_cloud=0.73):
    """Create radial profile animation."""
    
    # Find all snapshots
    snapshots = sorted(glob.glob(f"{results_dir}/snapshot_*.csv"))
    if not snapshots:
        print(f"No snapshots found in {results_dir}")
        return
    
    print(f"Found {len(snapshots)} snapshots")
    
    # Solve Lane-Emden for analytical profile
    xi_arr, psi_arr, dpsi_arr = solve_lane_emden_isothermal(xi_max=xi_s*1.5)
    
    # Convert to physical coordinates
    r_0 = R_cloud / xi_s  # Scale length
    r_analytical = xi_arr * r_0
    rho_ratio = np.exp(-psi_arr)  # ρ/ρ_c
    
    # Load first snapshot to get normalization
    data0 = load_snapshot(snapshots[0])
    cloud_mask = ~data0['is_ghost']
    
    # Get central density from first snapshot (particles near r=0)
    center_mask = (data0['r'] < 0.1) & cloud_mask
    if np.sum(center_mask) > 0:
        rho_c = np.median(data0['dens'][center_mask])
    else:
        rho_c = np.max(data0['dens'][cloud_mask])
    
    print(f"Central density (from sim): {rho_c:.2f}")
    print(f"R_cloud: {R_cloud}, r_0: {r_0:.4f}, xi_s: {xi_s}")
    
    # Setup figure with dark theme
    plt.style.use('dark_background')
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), facecolor=DARK_BG)
    fig.patch.set_facecolor(DARK_BG)
    
    for ax in [ax1, ax2]:
        ax.set_facecolor(DARK_PANEL)
        ax.tick_params(colors=TEXT_COLOR)
        ax.spines['bottom'].set_color(GRID_COLOR)
        ax.spines['top'].set_color(GRID_COLOR)
        ax.spines['left'].set_color(GRID_COLOR)
        ax.spines['right'].set_color(GRID_COLOR)
        ax.grid(True, alpha=0.3, color=GRID_COLOR)
    
    # Create bin edges for radial averaging
    r_bins = np.linspace(0, R_cloud * 1.5, 50)
    r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])
    
    def update(frame):
        ax1.clear()
        ax2.clear()
        
        for ax in [ax1, ax2]:
            ax.set_facecolor(DARK_PANEL)
            ax.tick_params(colors=TEXT_COLOR)
            ax.grid(True, alpha=0.3, color=GRID_COLOR)
        
        # Load data
        data = load_snapshot(snapshots[frame])
        time = get_time_from_snapshot(snapshots[frame])
        cloud_mask = ~data['is_ghost']
        
        r = data['r'][cloud_mask]
        dens = data['dens'][cloud_mask]
        a_r = data['a_r'][cloud_mask]
        g_r = data['g_r'][cloud_mask]
        
        # Compute radial averages
        dens_avg = np.zeros(len(r_centers))
        a_r_avg = np.zeros(len(r_centers))
        g_r_avg = np.zeros(len(r_centers))
        
        for i in range(len(r_centers)):
            mask = (r >= r_bins[i]) & (r < r_bins[i+1])
            if np.sum(mask) > 0:
                dens_avg[i] = np.median(dens[mask])
                a_r_avg[i] = np.median(a_r[mask])
                g_r_avg[i] = np.median(g_r[mask])
        
        # Top panel: Density profile
        ax1.scatter(r, dens, s=1, alpha=0.3, c='#00ffff', label='Particles')
        ax1.plot(r_centers, dens_avg, 'w-', lw=2, label='Radial average')
        ax1.plot(r_analytical, rho_c * rho_ratio, 'y--', lw=2, label='Analytical BE')
        ax1.axvline(R_cloud, color='#ff4444', ls=':', lw=1.5, label=f'R_cloud={R_cloud:.2f}')
        
        ax1.set_xlabel('Radius [pc]', color=TEXT_COLOR, fontsize=12)
        ax1.set_ylabel('Density [M☉/pc³]', color=TEXT_COLOR, fontsize=12)
        ax1.set_xlim(0, R_cloud * 1.5)
        ax1.set_ylim(0, rho_c * 1.5)
        ax1.legend(loc='upper right', fontsize=9, facecolor=DARK_PANEL, edgecolor=GRID_COLOR)
        ax1.set_title(f'Radial Density Profile (GSPH 1st Order)', color=TEXT_COLOR, fontsize=14)
        
        # Bottom panel: Force balance
        # Pressure gradient acceleration (a_r - g_r if gravity is included in a_r)
        # Or directly show a_r and g_r
        
        ax2.scatter(r, a_r, s=1, alpha=0.3, c='#00ff00', label='Hydro acc (radial)')
        ax2.scatter(r, g_r, s=1, alpha=0.3, c='#ff00ff', label='Gravity acc (radial)')
        ax2.plot(r_centers, a_r_avg, 'g-', lw=2, label='Hydro avg')
        ax2.plot(r_centers, g_r_avg, 'm-', lw=2, label='Gravity avg')
        ax2.axhline(0, color='white', ls='-', lw=0.5)
        ax2.axvline(R_cloud, color='#ff4444', ls=':', lw=1.5)
        
        ax2.set_xlabel('Radius [pc]', color=TEXT_COLOR, fontsize=12)
        ax2.set_ylabel('Radial Acceleration', color=TEXT_COLOR, fontsize=12)
        ax2.set_xlim(0, R_cloud * 1.5)
        
        # Auto-scale y-axis
        all_acc = np.concatenate([a_r, g_r])
        y_max = np.percentile(np.abs(all_acc[np.isfinite(all_acc)]), 99) * 1.2
        ax2.set_ylim(-y_max, y_max)
        
        ax2.legend(loc='upper right', fontsize=9, facecolor=DARK_PANEL, edgecolor=GRID_COLOR)
        ax2.set_title(f'Force Balance', color=TEXT_COLOR, fontsize=14)
        
        # Time annotation
        fig.suptitle(f'Bonnor-Ebert Relaxation  |  Step {frame}/{len(snapshots)-1}  |  t = {time:.3f} code',
                     color=TEXT_COLOR, fontsize=14, y=0.98)
        
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        
        return ax1, ax2
    
    # Create animation
    print(f"Creating animation with {len(snapshots)} frames...")
    anim = FuncAnimation(fig, update, frames=len(snapshots), interval=200, blit=False)
    
    # Save
    print(f"Saving to {output_path}...")
    writer = PillowWriter(fps=5)
    anim.save(output_path, writer=writer, dpi=100)
    
    plt.close()
    print(f"✓ Animation saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(description='Create radial profile animation')
    parser.add_argument('results_dir', help='Directory containing snapshot CSVs')
    parser.add_argument('-o', '--output', default='radial_profile.gif', help='Output GIF path')
    parser.add_argument('--xi_s', type=float, default=6.0, help='BE dimensionless radius')
    parser.add_argument('--R_cloud', type=float, default=0.73, help='Cloud radius in pc')
    
    args = parser.parse_args()
    
    create_animation(args.results_dir, args.output, args.xi_s, args.R_cloud)


if __name__ == '__main__':
    main()
