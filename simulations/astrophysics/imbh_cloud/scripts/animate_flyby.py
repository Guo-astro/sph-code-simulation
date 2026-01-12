#!/usr/bin/env python3
"""
Animate IMBH-Cloud Flyby Simulation
Shows cloud particles and BH position over time
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from pathlib import Path
import sys

# Configuration
RESULTS_DIR = Path("/Users/guo-opt-p148/sph-code-simulation/simulations/astrophysics/imbh_cloud/results/phase3_flyby")
OUTPUT_GIF = RESULTS_DIR / "imbh_flyby.gif"

# BH parameters (from config)
BH_X0, BH_Y0 = 15.0, 4.07  # Initial position
BH_VX = -80.0  # Velocity in x

def load_snapshot(filepath):
    """Load snapshot and extract time from header."""
    time = 0.0
    with open(filepath) as f:
        for line in f:
            if 'Time (code):' in line:
                time = float(line.split(':')[1].strip())
                break
    df = pd.read_csv(filepath, comment='#')
    return df, time

def main():
    # Find all snapshots
    snapshots = sorted(RESULTS_DIR.glob("snapshot_*.csv"))
    print(f"Found {len(snapshots)} snapshots")
    
    if len(snapshots) == 0:
        print("No snapshots found!")
        return
    
    # Load first snapshot to get limits
    df0, _ = load_snapshot(snapshots[0])
    
    # Set up figure
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle("IMBH-Cloud Tidal Flyby", fontsize=14, fontweight='bold')
    
    # Determine plot limits
    x_range = 20  # pc
    y_range = 12  # pc
    
    def animate(frame_idx):
        for ax in axes:
            ax.clear()
        
        snap_path = snapshots[frame_idx]
        df, t = load_snapshot(snap_path)
        
        # Calculate BH position
        bh_x = BH_X0 + BH_VX * t
        bh_y = BH_Y0
        bh_r = np.sqrt(bh_x**2 + bh_y**2)
        
        # Cloud properties
        cm_x = df['pos_x'].mean()
        cm_y = df['pos_y'].mean()
        v_bulk = np.sqrt(df['vel_x'].mean()**2 + df['vel_y'].mean()**2)
        sigma_v = np.std(np.sqrt(df['vel_x']**2 + df['vel_y']**2))
        
        # Left panel: X-Y view (top-down)
        ax1 = axes[0]
        
        # Color by velocity magnitude
        v_mag = np.sqrt(df['vel_x']**2 + df['vel_y']**2)
        scatter = ax1.scatter(df['pos_x'], df['pos_y'], c=v_mag, 
                             cmap='plasma', s=0.5, alpha=0.6, vmin=0, vmax=3)
        
        # Plot BH
        ax1.scatter([bh_x], [bh_y], c='black', s=200, marker='*', 
                   edgecolors='yellow', linewidths=1, zorder=10, label='IMBH')
        
        # Draw BH trajectory
        t_traj = np.linspace(0, 0.6, 100)
        bh_traj_x = BH_X0 + BH_VX * t_traj
        ax1.plot(bh_traj_x, [BH_Y0]*len(t_traj), 'k--', alpha=0.3, lw=1)
        
        # Mark pericenter
        ax1.axvline(x=0, color='gray', linestyle=':', alpha=0.5)
        ax1.scatter([0], [BH_Y0], c='red', s=50, marker='x', zorder=5, label='Pericenter')
        
        # Cloud center
        ax1.scatter([cm_x], [cm_y], c='cyan', s=100, marker='+', 
                   linewidths=2, zorder=8, label='Cloud CoM')
        
        ax1.set_xlim(-x_range, x_range)
        ax1.set_ylim(-y_range/2, y_range)
        ax1.set_xlabel('x [pc]', fontsize=11)
        ax1.set_ylabel('y [pc]', fontsize=11)
        ax1.set_aspect('equal')
        ax1.legend(loc='upper right', fontsize=9)
        ax1.set_title(f'X-Y Plane (t = {t:.3f})', fontsize=12)
        ax1.grid(True, alpha=0.3)
        
        # Right panel: Statistics
        ax2 = axes[1]
        
        # Velocity histogram
        ax2.hist(v_mag, bins=50, color='steelblue', alpha=0.7, density=True)
        ax2.axvline(x=v_bulk, color='red', linestyle='-', lw=2, label=f'Bulk: {v_bulk:.2f} km/s')
        ax2.axvline(x=sigma_v, color='orange', linestyle='--', lw=2, label=f'σ_v: {sigma_v:.2f} km/s')
        
        ax2.set_xlabel('Velocity magnitude [km/s]', fontsize=11)
        ax2.set_ylabel('Probability density', fontsize=11)
        ax2.set_xlim(0, 4)
        ax2.legend(loc='upper right', fontsize=10)
        ax2.set_title('Velocity Distribution', fontsize=12)
        ax2.grid(True, alpha=0.3)
        
        # Add info text
        info_text = (
            f"Frame: {frame_idx+1}/{len(snapshots)}\n"
            f"Time: {t:.4f} code = {t*0.978:.3f} Myr\n"
            f"BH position: ({bh_x:.1f}, {bh_y:.1f}) pc\n"
            f"BH distance: {bh_r:.2f} pc\n"
            f"Cloud CoM: ({cm_x:.3f}, {cm_y:.3f}) pc\n"
            f"Bulk velocity: {v_bulk:.3f} km/s\n"
            f"Dispersion σ: {sigma_v:.3f} km/s"
        )
        ax2.text(0.98, 0.02, info_text, transform=ax2.transAxes, 
                fontsize=9, verticalalignment='bottom', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
                family='monospace')
        
        plt.tight_layout()
        return axes
    
    print("Creating animation...")
    anim = FuncAnimation(fig, animate, frames=len(snapshots), interval=100)
    
    print(f"Saving GIF to {OUTPUT_GIF}...")
    writer = PillowWriter(fps=10)
    anim.save(OUTPUT_GIF, writer=writer, dpi=100)
    
    print(f"✓ GIF saved: {OUTPUT_GIF}")
    plt.close()

if __name__ == '__main__':
    main()
