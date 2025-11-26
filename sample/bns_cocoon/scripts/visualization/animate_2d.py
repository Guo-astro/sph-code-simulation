#!/usr/bin/env python3
"""
BNS Cocoon 2D Animation

Creates animated GIF from simulation snapshots.
"""

import argparse
import glob
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LogNorm
from matplotlib.animation import FuncAnimation, PillowWriter


def load_csv(filename):
    """Load CSV snapshot with header comments."""
    import pandas as pd
    # Read skipping comment lines
    data = pd.read_csv(filename, comment='#')
    return data


def make_animation(input_dir, output_file, pattern='snapshot_*.csv'):
    """Create animation from snapshots."""
    
    # Find all snapshot files
    files = sorted(glob.glob(os.path.join(input_dir, pattern)))
    if not files:
        print(f"No files found matching {pattern}")
        sys.exit(1)
    
    print(f"Found {len(files)} snapshots")
    
    # Load first snapshot to get bounds
    data0 = load_csv(files[0])
    x0 = data0['pos_x']
    y0 = data0['pos_y']
    x_lim = [x0.min() * 1.1, x0.max() * 1.1]
    y_lim = [y0.min() * 1.1, y0.max() * 1.1]
    
    # Get density range across all snapshots (sample a few)
    sample_indices = np.linspace(0, len(files)-1, min(10, len(files)), dtype=int)
    all_dens = []
    for idx in sample_indices:
        d = load_csv(files[idx])
        all_dens.extend(d['dens'])
    dens_min = max(np.percentile(all_dens, 1), 1e-12)
    dens_max = np.percentile(all_dens, 99)
    
    # Setup figure
    fig, ax = plt.subplots(figsize=(8, 8))
    
    def update(frame):
        ax.clear()
        data = load_csv(files[frame])
        x = data['pos_x']
        y = data['pos_y']
        dens = data['dens']
        
        sc = ax.scatter(x, y, c=dens, s=2, cmap='viridis',
                        norm=LogNorm(vmin=dens_min, vmax=dens_max))
        
        ax.set_xlim(x_lim)
        ax.set_ylim(y_lim)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title(f'BNS Cocoon 2D - Frame {frame:04d}')
        ax.set_aspect('equal')
        
        return sc,
    
    print("Generating animation...")
    anim = FuncAnimation(fig, update, frames=len(files), 
                         interval=100, blit=False)
    
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    anim.save(output_file, writer=PillowWriter(fps=10))
    print(f"Animation saved: {output_file}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Animate 2D BNS Cocoon simulation')
    parser.add_argument('--input', '-i', required=True,
                        help='Input directory with CSV files')
    parser.add_argument('--output', '-o', required=True,
                        help='Output GIF file')
    parser.add_argument('--pattern', default='snapshot_*.csv',
                        help='File pattern to match')
    args = parser.parse_args()
    
    make_animation(args.input, args.output, args.pattern)


if __name__ == '__main__':
    main()
