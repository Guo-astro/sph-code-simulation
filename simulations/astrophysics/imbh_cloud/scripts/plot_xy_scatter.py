#!/usr/bin/env python3
"""
3D and 2D scatter plots of SPH particles colored by various quantities.

Shows particle distribution with envelope (ghost) particles distinguished
from cloud (real) particles. Supports both CSV and HDF5 snapshots.

Usage:
    python plot_xy_scatter.py <snapshot> -o output.png [options]
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from mpl_toolkits.mplot3d import Axes3D
import os


def load_snapshot(filepath):
    """Load snapshot (CSV or HDF5) and return as dict."""
    if filepath.endswith('.h5') or filepath.endswith('.hdf5'):
        import h5py
        with h5py.File(filepath, 'r') as f:
            data = {}
            for key in f['particles'].keys():
                data[key] = np.array(f['particles'][key])
            # Rename for compatibility
            if 'pos_x' in data:
                data['pos_x'] = data['pos_x']
                data['pos_y'] = data['pos_y']
                data['pos_z'] = data['pos_z']
            if 'dens' in data and 'dens' not in data:
                pass  # already named correctly
            return data
    else:
        import pandas as pd
        df = pd.read_csv(filepath, comment='#')
        return {col: df[col].values for col in df.columns}


def get_time_from_snapshot(filepath):
    """Extract simulation time from snapshot header."""
    if filepath.endswith('.h5') or filepath.endswith('.hdf5'):
        import h5py
        with h5py.File(filepath, 'r') as f:
            if 'metadata' in f and 'time' in f['metadata']:
                return float(np.array(f['metadata/time'])[0])
            elif 'time' in f:
                return float(np.array(f['time'])[0])
        return 0.0
    else:
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('# time:'):
                    return float(line.split(':')[1].strip())
        return 0.0


def main():
    parser = argparse.ArgumentParser(description='3D and 2D scatter plots of SPH particles')
    parser.add_argument('snapshot', help='Path to snapshot (CSV or HDF5)')
    parser.add_argument('-o', '--output', default='xy_scatter.png', help='Output filename')
    parser.add_argument('--color_by', default='dens',
                        choices=['dens', 'ene', 'temp', 'pres', 'sound'],
                        help='Quantity to color by')
    parser.add_argument('--density_to_nH', type=float, default=31.86,
                        help='Code density to number density [cm^-3]')
    parser.add_argument('--u_to_cgs', type=float, default=1.0e10,
                        help='Code energy to CGS [erg/g]')
    parser.add_argument('--R_cloud', type=float, default=1.55, help='Cloud radius [pc]')
    parser.add_argument('--log_color', action='store_true', help='Use log scale for color')
    parser.add_argument('--show_envelope', action='store_true',
                        help='Show envelope particles')
    parser.add_argument('--max_particles', type=int, default=20000,
                        help='Max particles to plot (subsample if more)')
    args = parser.parse_args()

    # Physical constants for temperature conversion
    k_B = 1.3807e-16  # erg/K
    m_p = 1.6726e-24  # g
    mu = 1.27
    gamma = 5.0/3.0  # For temperature conversion (not for simulation)

    # Load data
    data = load_snapshot(args.snapshot)
    time = get_time_from_snapshot(args.snapshot)

    # Separate real cloud particles from envelope/ghost particles
    is_ghost = data.get('is_ghost', np.zeros(len(data['pos_x'])))
    real_mask = is_ghost == 0
    ghost_mask = is_ghost == 1

    # Extract arrays for cloud and envelope
    cloud = {k: v[real_mask] for k, v in data.items()}
    envelope = {k: v[ghost_mask] for k, v in data.items()}

    n_cloud = len(cloud['pos_x'])
    n_envelope = len(envelope['pos_x']) if 'pos_x' in envelope else 0
    
    print(f"Snapshot: {os.path.basename(args.snapshot)}")
    print(f"Time: {time:.2f} code units")
    print(f"Cloud particles: {n_cloud:,}")
    print(f"Envelope particles: {n_envelope:,}")
    print(f"Total: {n_cloud + n_envelope:,}")

    # Calculate color quantity
    def get_color_quantity(d, name):
        if name == 'dens':
            return d['dens'] * args.density_to_nH, 'Number Density [cm⁻³]'
        elif name == 'ene':
            return d['ene'], 'Internal Energy [code]'
        elif name == 'temp':
            u_cgs = d['ene'] * args.u_to_cgs
            T = (gamma - 1) * mu * m_p * u_cgs / k_B
            return T, 'Temperature [K]'
        elif name == 'pres':
            return d['pres'], 'Pressure [code]'
        elif name == 'sound':
            return d['sound'], 'Sound Speed [code]'

    cloud_color, color_label = get_color_quantity(cloud, args.color_by)

    # Statistics
    print(f"\n=== Cloud Statistics ({args.color_by}) ===")
    print(f"  Min: {cloud_color.min():.4g}")
    print(f"  Max: {cloud_color.max():.4g}")
    print(f"  Mean: {cloud_color.mean():.4g}")

    envelope_color = None
    if n_envelope > 0:
        envelope_color, _ = get_color_quantity(envelope, args.color_by)
        print(f"\n=== Envelope Statistics ({args.color_by}) ===")
        print(f"  Min: {envelope_color.min():.4g}")
        print(f"  Max: {envelope_color.max():.4g}")

    # Subsample for plotting if needed
    step = max(1, n_cloud // args.max_particles)
    if step > 1:
        print(f"\nSubsampling: showing 1/{step} particles")
    idx = slice(None, None, step)
    
    # Create figure with 3D + 2D projections
    fig = plt.figure(figsize=(16, 12))

    # 3D view
    ax3d = fig.add_subplot(2, 2, 1, projection='3d')
    # XY projection
    ax_xy = fig.add_subplot(2, 2, 2)
    # XZ projection
    ax_xz = fig.add_subplot(2, 2, 3)
    # Radial profile
    ax_rad = fig.add_subplot(2, 2, 4)

    # Set up color normalization
    vmin, vmax = np.percentile(cloud_color, [1, 99])
    if args.log_color and vmin > 0:
        norm = LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = Normalize(vmin=vmin, vmax=vmax)

    # Extract subsampled coordinates
    x, y, z = cloud['pos_x'][idx], cloud['pos_y'][idx], cloud['pos_z'][idx]
    c = cloud_color[idx]

    # 3D scatter plot
    if args.show_envelope and n_envelope > 0:
        env_step = max(1, n_envelope // (args.max_particles // 4))
        ax3d.scatter(envelope['pos_x'][::env_step], envelope['pos_y'][::env_step],
                    envelope['pos_z'][::env_step], c='lightgray', s=0.5, alpha=0.2)

    sc3d = ax3d.scatter(x, y, z, c=c, s=1, alpha=0.6, cmap='viridis', norm=norm)
    ax3d.set_xlabel('X [pc]')
    ax3d.set_ylabel('Y [pc]')
    ax3d.set_zlabel('Z [pc]')
    ax3d.set_title(f'3D View (t = {time:.2f})')
    lim = args.R_cloud * 1.5
    ax3d.set_xlim(-lim, lim)
    ax3d.set_ylim(-lim, lim)
    ax3d.set_zlim(-lim, lim)

    # Draw sphere wireframe for cloud boundary
    u = np.linspace(0, 2 * np.pi, 30)
    v = np.linspace(0, np.pi, 20)
    xs = args.R_cloud * np.outer(np.cos(u), np.sin(v))
    ys = args.R_cloud * np.outer(np.sin(u), np.sin(v))
    zs = args.R_cloud * np.outer(np.ones(np.size(u)), np.cos(v))
    ax3d.plot_wireframe(xs, ys, zs, color='red', alpha=0.1, linewidth=0.5)

    # XY projection
    theta = np.linspace(0, 2*np.pi, 100)
    if args.show_envelope and n_envelope > 0:
        ax_xy.scatter(envelope['pos_x'][::env_step], envelope['pos_y'][::env_step],
                     c='lightgray', s=0.5, alpha=0.3)

    sc_xy = ax_xy.scatter(x, y, c=c, s=1, alpha=0.6, cmap='viridis', norm=norm)
    ax_xy.plot(args.R_cloud * np.cos(theta), args.R_cloud * np.sin(theta),
               'r--', lw=1.5, label=f'R = {args.R_cloud} pc')
    ax_xy.set_xlabel('X [pc]')
    ax_xy.set_ylabel('Y [pc]')
    ax_xy.set_aspect('equal')
    ax_xy.set_title('XY Projection')
    ax_xy.set_xlim(-lim, lim)
    ax_xy.set_ylim(-lim, lim)
    ax_xy.legend(loc='upper right', fontsize=8)
    plt.colorbar(sc_xy, ax=ax_xy, label=color_label)

    # XZ projection
    if args.show_envelope and n_envelope > 0:
        ax_xz.scatter(envelope['pos_x'][::env_step], envelope['pos_z'][::env_step],
                     c='lightgray', s=0.5, alpha=0.3)

    sc_xz = ax_xz.scatter(x, z, c=c, s=1, alpha=0.6, cmap='viridis', norm=norm)
    ax_xz.plot(args.R_cloud * np.cos(theta), args.R_cloud * np.sin(theta),
               'r--', lw=1.5, label=f'R = {args.R_cloud} pc')
    ax_xz.set_xlabel('X [pc]')
    ax_xz.set_ylabel('Z [pc]')
    ax_xz.set_aspect('equal')
    ax_xz.set_title('XZ Projection')
    ax_xz.set_xlim(-lim, lim)
    ax_xz.set_ylim(-lim, lim)
    ax_xz.legend(loc='upper right', fontsize=8)
    plt.colorbar(sc_xz, ax=ax_xz, label=color_label)

    # Radial profile
    r_cloud = np.sqrt(cloud['pos_x']**2 + cloud['pos_y']**2 + cloud['pos_z']**2)
    ax_rad.scatter(r_cloud[idx], cloud_color[idx], s=0.5, alpha=0.3, c='blue', label='Cloud')

    if args.show_envelope and n_envelope > 0:
        r_env = np.sqrt(envelope['pos_x']**2 + envelope['pos_y']**2 + envelope['pos_z']**2)
        ax_rad.scatter(r_env[::env_step], envelope_color[::env_step],
                      s=0.5, alpha=0.3, c='red', label='Envelope')

    ax_rad.axvline(args.R_cloud, color='r', ls='--', lw=1.5, label=f'R_cloud')
    ax_rad.set_xlabel('Radius [pc]')
    ax_rad.set_ylabel(color_label)
    ax_rad.set_title(f'{args.color_by.capitalize()} vs Radius')
    ax_rad.legend(fontsize=8)
    ax_rad.set_xlim(0, lim)
    if args.log_color and vmin > 0:
        ax_rad.set_yscale('log')

    plt.suptitle(f'{os.path.basename(args.snapshot)} | Cloud: {n_cloud:,} | Envelope: {n_envelope:,}',
                 fontsize=11)

    plt.tight_layout()
    plt.savefig(args.output, dpi=150, bbox_inches='tight')
    print(f"\nSaved: {args.output}")
    plt.close()


if __name__ == '__main__':
    main()
