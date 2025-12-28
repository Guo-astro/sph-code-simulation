#!/usr/bin/env python3
"""
Plot MHD Shock Tube results from SPH simulation.
Dai-Woodward shock tube test (Iwasaki & Inutsuka 2011 Fig. 5)

Usage:
    python plot_mhd_results.py                    # Plot final state
    python plot_mhd_results.py --animate          # Create animation GIF
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from pathlib import Path
import glob
import argparse


def read_sph_csv(filepath):
    """Read SPH CSV output file, skipping comment lines."""
    with open(filepath, 'r') as f:
        skip_lines = 0
        for line in f:
            if line.startswith('#'):
                skip_lines += 1
            else:
                break
    return pd.read_csv(filepath, skiprows=skip_lines)


def get_initial_condition():
    """Return the analytic initial condition for Dai-Woodward test."""
    sqrt_4pi = np.sqrt(4.0 * np.pi)

    # Left state (x < 0)
    left = {
        'rho': 1.08, 'P': 0.95, 'vx': 1.2, 'vy': 0.01, 'vz': 0.5,
        'Bx': 2.0/sqrt_4pi, 'By': 3.6/sqrt_4pi, 'Bz': 2.0/sqrt_4pi
    }

    # Right state (x > 0)
    right = {
        'rho': 1.0, 'P': 1.0, 'vx': 0.0, 'vy': 0.0, 'vz': 0.0,
        'Bx': 2.0/sqrt_4pi, 'By': 4.0/sqrt_4pi, 'Bz': 2.0/sqrt_4pi
    }

    x = np.linspace(-0.75, 0.51, 1000)
    ic = {}
    for key in left.keys():
        ic[key] = np.where(x < 0, left[key], right[key])
    ic['x'] = x

    return ic


def create_final_plot(result_dir):
    """Create final comparison plot with initial condition overlay."""
    snapshot_files = sorted(glob.glob(str(result_dir / 'snapshot_*.csv')))
    if not snapshot_files:
        print("No snapshot files found!")
        return

    data = read_sph_csv(snapshot_files[-1])
    data = data[data['is_ghost'] == 0]
    ic = get_initial_condition()

    idx = int(Path(snapshot_files[-1]).stem.split('_')[1])
    time = idx * 0.02

    fig, axes = plt.subplots(3, 3, figsize=(15, 12))
    fig.patch.set_facecolor('white')

    plots = [
        ('pos_x', 'dens', 'Density', 'rho', r'$\rho$'),
        ('pos_x', 'pres', 'Pressure', 'P', r'$P$'),
        ('pos_x', 'vel_x', 'Velocity (x)', 'vx', r'$v_x$'),
        ('pos_x', 'vy_mhd', 'Velocity (y)', 'vy', r'$v_y$'),
        ('pos_x', 'vz_mhd', 'Velocity (z)', 'vz', r'$v_z$'),
        ('pos_x', 'B_x', 'Magnetic Field (x)', 'Bx', r'$B_x$'),
        ('pos_x', 'B_y', 'Magnetic Field (y)', 'By', r'$B_y$'),
        ('pos_x', 'B_z', 'Magnetic Field (z)', 'Bz', r'$B_z$'),
        ('pos_x', 'ene', 'Internal Energy', None, r'$u$'),
    ]

    for ax, (xcol, ycol, title, ic_key, ylabel) in zip(axes.flat, plots):
        ax.scatter(data[xcol], data[ycol], s=2, c='blue', alpha=0.7, label='GSPMHD (t=0.2)')
        if ic_key and ic_key in ic:
            ax.plot(ic['x'], ic[ic_key], 'r--', lw=1.5, alpha=0.6, label='Initial (t=0)')

        ax.set_xlabel('x', fontsize=10)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_title(title, fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.legend(loc='best', fontsize=8)

        y = data[ycol].values
        q1, q99 = np.percentile(y, [2, 98])
        margin = (q99 - q1) * 0.15
        if margin > 0:
            ax.set_ylim(q1 - margin, q99 + margin)
        ax.set_xlim(-0.75, 0.55)

    fig.suptitle(f'MHD Shock Tube (Dai-Woodward) at t = {time:.2f}\n'
                 'Iwasaki & Inutsuka (2011) - GSPMHD Method',
                 fontsize=13, fontweight='bold')

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    output_path = result_dir / 'mhd_final_with_ic.png'
    plt.savefig(str(output_path), dpi=150, bbox_inches='tight', facecolor='white')
    print(f"Saved: {output_path}")
    plt.close()


def create_animation(result_dir):
    """Create animation of MHD shock tube evolution."""
    snapshot_files = sorted(glob.glob(str(result_dir / 'snapshot_*.csv')))
    if not snapshot_files:
        print("No snapshot files found!")
        return

    print(f"Found {len(snapshot_files)} snapshots")
    ic = get_initial_condition()

    snapshots = []
    times = []
    for f in snapshot_files:
        data = read_sph_csv(f)
        data = data[data['is_ghost'] == 0]
        snapshots.append(data)
        idx = int(Path(f).stem.split('_')[1])
        times.append(idx * 0.02)

    fig, axes = plt.subplots(3, 3, figsize=(14, 11))
    fig.patch.set_facecolor('white')

    plots = [
        ('pos_x', 'dens', 'Density', 'rho', r'$\rho$'),
        ('pos_x', 'pres', 'Pressure', 'P', r'$P$'),
        ('pos_x', 'vel_x', 'Velocity (x)', 'vx', r'$v_x$'),
        ('pos_x', 'vy_mhd', 'Velocity (y)', 'vy', r'$v_y$'),
        ('pos_x', 'vz_mhd', 'Velocity (z)', 'vz', r'$v_z$'),
        ('pos_x', 'B_x', 'Magnetic Field (x)', 'Bx', r'$B_x$'),
        ('pos_x', 'B_y', 'Magnetic Field (y)', 'By', r'$B_y$'),
        ('pos_x', 'B_z', 'Magnetic Field (z)', 'Bz', r'$B_z$'),
        ('pos_x', 'ene', 'Internal Energy', None, r'$u$'),
    ]

    xlim = (-0.75, 0.55)
    ylims = {}
    for xcol, ycol, title, ic_key, ylabel in plots:
        all_y = np.concatenate([s[ycol].values for s in snapshots])
        q1, q99 = np.percentile(all_y, [1, 99])
        margin = (q99 - q1) * 0.1
        ylims[ycol] = (q1 - margin, q99 + margin) if margin > 0 else (q1 - 0.1, q99 + 0.1)

    scatters = []
    for idx, (ax, (xcol, ycol, title, ic_key, ylabel)) in enumerate(zip(axes.flat, plots)):
        sc = ax.scatter([], [], s=1, c='blue', alpha=0.7, label='GSPMHD')
        scatters.append(sc)
        if ic_key and ic_key in ic:
            ax.plot(ic['x'], ic[ic_key], 'r-', lw=1.5, alpha=0.7, label='Initial')
        ax.set_xlim(xlim)
        ax.set_ylim(ylims[ycol])
        ax.set_xlabel('x')
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        if ic_key:
            ax.legend(loc='best', fontsize=8)

    title_text = fig.suptitle('', fontsize=14, fontweight='bold')

    def animate(frame):
        data = snapshots[frame]
        for sc, (xcol, ycol, _, _, _) in zip(scatters, plots):
            sc.set_offsets(np.c_[data[xcol], data[ycol]])
        title_text.set_text(f'MHD Shock Tube (Dai-Woodward) - t = {times[frame]:.3f}')
        return scatters + [title_text]

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    anim = animation.FuncAnimation(fig, animate, frames=len(snapshots), interval=500, blit=False)

    output_path = result_dir / 'mhd_shock_tube_animation.gif'
    print(f"Saving animation to {output_path}...")
    anim.save(str(output_path), writer='pillow', fps=2, dpi=100)
    print("Animation saved!")
    plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot MHD shock tube results')
    parser.add_argument('--animate', action='store_true', help='Create animation GIF')
    parser.add_argument('--dir', type=str, default='simulations/mhd/results/mhd_shock_tube_1',
                        help='Results directory')
    args = parser.parse_args()

    result_dir = Path(args.dir)
    if args.animate:
        create_animation(result_dir)
    create_final_plot(result_dir)
