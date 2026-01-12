#!/usr/bin/env python3
"""
XY scatter plot of SPH particles colored by various quantities.

Shows particle distribution with envelope (ghost) particles distinguished
from cloud (real) particles.

Usage:
    python plot_xy_scatter.py <snapshot> -o output.png [options]
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize


def load_snapshot(filepath):
    """Load snapshot CSV with comment handling."""
    return pd.read_csv(filepath, comment='#')


def get_time_from_snapshot(filepath):
    """Extract simulation time from snapshot header."""
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('# time:'):
                return float(line.split(':')[1].strip())
    return 0.0


def main():
    parser = argparse.ArgumentParser(description='XY scatter plot of SPH particles')
    parser.add_argument('snapshot', help='Path to snapshot CSV')
    parser.add_argument('-o', '--output', default='xy_scatter.png', help='Output filename')
    parser.add_argument('--color_by', default='ene', 
                        choices=['dens', 'ene', 'temp', 'pres', 'sound'],
                        help='Quantity to color by')
    parser.add_argument('--density_to_nH', type=float, default=31.86,
                        help='Code density to number density [cm^-3]')
    parser.add_argument('--u_to_cgs', type=float, default=1.0e10,
                        help='Code energy to CGS [erg/g]')
    parser.add_argument('--R_cloud', type=float, default=0.75, help='Cloud radius [pc]')
    parser.add_argument('--log_color', action='store_true', help='Use log scale for color')
    parser.add_argument('--show_envelope', action='store_true', 
                        help='Show envelope particles')
    args = parser.parse_args()
    
    # Physical constants for temperature conversion
    k_B = 1.3807e-16  # erg/K
    m_p = 1.6726e-24  # g
    mu = 1.27
    gamma = 5.0/3.0  # For temperature conversion (not for simulation)
    
    # Load data
    df = load_snapshot(args.snapshot)
    time = get_time_from_snapshot(args.snapshot)
    
    # Separate real cloud particles from envelope/ghost particles
    real_mask = df['is_ghost'] == 0
    ghost_mask = df['is_ghost'] == 1
    
    df_cloud = df[real_mask]
    df_envelope = df[ghost_mask]
    
    n_cloud = len(df_cloud)
    n_envelope = len(df_envelope)
    
    print(f"Snapshot time: {time:.2f} code time")
    print(f"Cloud particles (is_ghost=0): {n_cloud}")
    print(f"Envelope particles (is_ghost=1): {n_envelope}")
    print(f"Total particles: {len(df)}")
    
    # Calculate color quantity
    def get_color_quantity(df_in, name):
        if name == 'dens':
            return df_in['dens'].values * args.density_to_nH, 'Number Density [cm⁻³]'
        elif name == 'ene':
            return df_in['ene'].values, 'Internal Energy [code]'
        elif name == 'temp':
            u_cgs = df_in['ene'].values * args.u_to_cgs
            T = (gamma - 1) * mu * m_p * u_cgs / k_B
            return T, 'Temperature [K]'
        elif name == 'pres':
            return df_in['pres'].values, 'Pressure [code]'
        elif name == 'sound':
            return df_in['sound'].values, 'Sound Speed [code]'
    
    cloud_color, color_label = get_color_quantity(df_cloud, args.color_by)
    
    # Statistics
    print(f"\n=== Cloud Particle Statistics ({args.color_by}) ===")
    print(f"  Min: {cloud_color.min():.4g}")
    print(f"  Max: {cloud_color.max():.4g}")
    print(f"  Mean: {cloud_color.mean():.4g}")
    print(f"  Median: {np.median(cloud_color):.4g}")
    
    if n_envelope > 0:
        envelope_color, _ = get_color_quantity(df_envelope, args.color_by)
        print(f"\n=== Envelope Particle Statistics ({args.color_by}) ===")
        print(f"  Min: {envelope_color.min():.4g}")
        print(f"  Max: {envelope_color.max():.4g}")
        print(f"  Mean: {envelope_color.mean():.4g}")
        print(f"  Median: {np.median(envelope_color):.4g}")
    
    # Create figure with multiple views
    fig = plt.figure(figsize=(16, 12))
    
    # XY projection (z=0 slice)
    ax1 = fig.add_subplot(2, 2, 1)
    # XZ projection
    ax2 = fig.add_subplot(2, 2, 2)
    # Color vs radius for cloud
    ax3 = fig.add_subplot(2, 2, 3)
    # Histogram of color quantity
    ax4 = fig.add_subplot(2, 2, 4)
    
    # Set up color normalization
    vmin, vmax = np.percentile(cloud_color, [1, 99])
    if args.log_color and vmin > 0:
        norm = LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = Normalize(vmin=vmin, vmax=vmax)
    
    # Plot XY projection
    if args.show_envelope and n_envelope > 0:
        ax1.scatter(df_envelope['pos_x'], df_envelope['pos_y'], 
                   c='lightgray', s=1, alpha=0.3, label='Envelope')
    
    sc1 = ax1.scatter(df_cloud['pos_x'], df_cloud['pos_y'], 
                     c=cloud_color, s=2, alpha=0.7, cmap='viridis', norm=norm)
    
    # Add cloud boundary circle
    theta = np.linspace(0, 2*np.pi, 100)
    ax1.plot(args.R_cloud * np.cos(theta), args.R_cloud * np.sin(theta), 
             'r--', lw=2, label=f'R_cloud = {args.R_cloud} pc')
    
    ax1.set_xlabel('X [pc]')
    ax1.set_ylabel('Y [pc]')
    ax1.set_aspect('equal')
    ax1.set_title(f'XY Projection (t = {time:.2f})')
    ax1.legend(loc='upper right')
    cbar1 = plt.colorbar(sc1, ax=ax1)
    cbar1.set_label(color_label)
    
    # Plot XZ projection
    if args.show_envelope and n_envelope > 0:
        ax2.scatter(df_envelope['pos_x'], df_envelope['pos_z'], 
                   c='lightgray', s=1, alpha=0.3, label='Envelope')
    
    sc2 = ax2.scatter(df_cloud['pos_x'], df_cloud['pos_z'], 
                     c=cloud_color, s=2, alpha=0.7, cmap='viridis', norm=norm)
    
    ax2.plot(args.R_cloud * np.cos(theta), args.R_cloud * np.sin(theta), 
             'r--', lw=2, label=f'R_cloud = {args.R_cloud} pc')
    
    ax2.set_xlabel('X [pc]')
    ax2.set_ylabel('Z [pc]')
    ax2.set_aspect('equal')
    ax2.set_title(f'XZ Projection (t = {time:.2f})')
    ax2.legend(loc='upper right')
    cbar2 = plt.colorbar(sc2, ax=ax2)
    cbar2.set_label(color_label)
    
    # Color vs radius scatter
    r_cloud = np.sqrt(df_cloud['pos_x']**2 + df_cloud['pos_y']**2 + df_cloud['pos_z']**2)
    ax3.scatter(r_cloud, cloud_color, s=1, alpha=0.3, c='blue', label='Cloud')
    
    if args.show_envelope and n_envelope > 0:
        r_envelope = np.sqrt(df_envelope['pos_x']**2 + df_envelope['pos_y']**2 + df_envelope['pos_z']**2)
        envelope_color_vals, _ = get_color_quantity(df_envelope, args.color_by)
        ax3.scatter(r_envelope, envelope_color_vals, s=1, alpha=0.3, c='red', label='Envelope')
    
    ax3.axvline(args.R_cloud, color='r', ls='--', lw=2, label=f'R_cloud = {args.R_cloud} pc')
    ax3.set_xlabel('Radius [pc]')
    ax3.set_ylabel(color_label)
    ax3.set_title(f'{color_label} vs Radius')
    ax3.legend()
    if args.log_color and vmin > 0:
        ax3.set_yscale('log')
    
    # Histogram
    ax4.hist(cloud_color, bins=50, alpha=0.7, label='Cloud', color='blue', density=True)
    if args.show_envelope and n_envelope > 0:
        ax4.hist(envelope_color_vals, bins=50, alpha=0.5, label='Envelope', color='red', density=True)
    ax4.set_xlabel(color_label)
    ax4.set_ylabel('Probability Density')
    ax4.set_title(f'Distribution of {color_label}')
    ax4.legend()
    if args.log_color and vmin > 0:
        ax4.set_xscale('log')
    
    plt.suptitle(f'SPH Particle Distribution - {args.snapshot.split("/")[-1]}\n'
                 f'Cloud: {n_cloud} particles, Envelope: {n_envelope} particles', 
                 fontsize=12, y=1.02)
    
    plt.tight_layout()
    plt.savefig(args.output, dpi=150, bbox_inches='tight')
    print(f"\n✓ Saved: {args.output}")
    plt.close()


if __name__ == '__main__':
    main()
