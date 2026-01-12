#!/usr/bin/env python3
"""
Generate XY slice GIF for IMBH-cloud simulation.

Usage:
    python make_xy_gif.py <results_dir>
    
Example:
    python make_xy_gif.py results/phase2.5_cooling_gamma53
"""
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import Normalize

# Constants for unit conversion
GAMMA = 5.0 / 3.0
MU = 1.27
K_B = 1.38e-16  # erg/K
M_H = 1.67e-24  # g
U_TO_CGS = 1.0e10  # code energy to erg/g
N_TO_CGS = 31.86   # code density to n_H cm^-3
T_TO_MYR = 0.978   # code time to Myr

def u_to_T(u):
    """Convert internal energy to temperature."""
    u_cgs = u * U_TO_CGS
    return (GAMMA - 1) * MU * M_H * u_cgs / K_B

def read_snapshot(filepath):
    """Read snapshot and return particle data."""
    t = 0.0
    pos_x, pos_y, pos_z, ene, dens = [], [], [], [], []
    
    with open(filepath) as f:
        for line in f:
            if "Time (code):" in line:
                t = float(line.split(":")[-1].strip())
            if line.startswith('#') or 'id,' in line:
                continue
            parts = line.strip().split(',')
            if len(parts) > 13:
                try:
                    x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                    e, d = float(parts[13]), float(parts[11])
                    if not (np.isnan(x) or np.isnan(e)):
                        pos_x.append(x)
                        pos_y.append(y)
                        pos_z.append(z)
                        ene.append(e)
                        dens.append(d)
                except:
                    pass
    
    return {
        't': t,
        'x': np.array(pos_x),
        'y': np.array(pos_y),
        'z': np.array(pos_z),
        'T': u_to_T(np.array(ene)),
        'n_H': np.array(dens) * N_TO_CGS
    }

def main(results_dir):
    snapshots = sorted(glob.glob(f"{results_dir}/snapshot_*.csv"))
    
    if not snapshots:
        print(f"❌ No snapshots found in {results_dir}")
        sys.exit(1)
    
    # Filter out NaN snapshots
    valid_snapshots = []
    for snap in snapshots:
        data = read_snapshot(snap)
        if len(data['x']) > 0 and not np.isnan(data['T']).all():
            valid_snapshots.append(snap)
    
    print(f"Found {len(valid_snapshots)} valid snapshots")
    
    if len(valid_snapshots) < 2:
        print("❌ Need at least 2 valid snapshots for GIF")
        sys.exit(1)
    
    # Create figure with 2 panels: Temperature and Density
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Get initial data
    data0 = read_snapshot(valid_snapshots[0])
    
    # XY slice (|z| < 0.1)
    z_slice = np.abs(data0['z']) < 0.1
    
    # Temperature plot
    sc_T = axes[0].scatter(data0['x'][z_slice], data0['y'][z_slice], 
                           c=data0['T'][z_slice], s=3, cmap='coolwarm',
                           vmin=5, vmax=12)
    axes[0].set_xlabel('X [pc]')
    axes[0].set_ylabel('Y [pc]')
    axes[0].set_title('Temperature (|z| < 0.1 pc)')
    axes[0].set_xlim(-1, 1)
    axes[0].set_ylim(-1, 1)
    axes[0].set_aspect('equal')
    cbar_T = plt.colorbar(sc_T, ax=axes[0], label='T [K]')
    
    # Density plot
    sc_n = axes[1].scatter(data0['x'][z_slice], data0['y'][z_slice], 
                           c=np.log10(data0['n_H'][z_slice]), s=3, cmap='viridis',
                           vmin=1, vmax=4)
    axes[1].set_xlabel('X [pc]')
    axes[1].set_ylabel('Y [pc]')
    axes[1].set_title('Density (|z| < 0.1 pc)')
    axes[1].set_xlim(-1, 1)
    axes[1].set_ylim(-1, 1)
    axes[1].set_aspect('equal')
    cbar_n = plt.colorbar(sc_n, ax=axes[1], label='log₁₀(n_H [cm⁻³])')
    
    time_text = fig.suptitle(f't = 0.00 Myr  |  T_mean = {data0["T"].mean():.1f} K', 
                             fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    def update(frame):
        data = read_snapshot(valid_snapshots[frame])
        z_slice = np.abs(data['z']) < 0.1
        
        # Clear and redraw
        axes[0].clear()
        axes[1].clear()
        
        # Temperature
        sc_T = axes[0].scatter(data['x'][z_slice], data['y'][z_slice], 
                               c=data['T'][z_slice], s=3, cmap='coolwarm',
                               vmin=5, vmax=12)
        axes[0].set_xlabel('X [pc]')
        axes[0].set_ylabel('Y [pc]')
        axes[0].set_title('Temperature (|z| < 0.1 pc)')
        axes[0].set_xlim(-1, 1)
        axes[0].set_ylim(-1, 1)
        axes[0].set_aspect('equal')
        
        # Density
        sc_n = axes[1].scatter(data['x'][z_slice], data['y'][z_slice], 
                               c=np.log10(data['n_H'][z_slice]), s=3, cmap='viridis',
                               vmin=1, vmax=4)
        axes[1].set_xlabel('X [pc]')
        axes[1].set_ylabel('Y [pc]')
        axes[1].set_title('Density (|z| < 0.1 pc)')
        axes[1].set_xlim(-1, 1)
        axes[1].set_ylim(-1, 1)
        axes[1].set_aspect('equal')
        
        t_myr = data['t'] * T_TO_MYR
        time_text.set_text(f't = {t_myr:.2f} Myr  |  T_mean = {data["T"].mean():.1f} K')
        
        return time_text,
    
    anim = FuncAnimation(fig, update, frames=len(valid_snapshots), 
                         interval=500, blit=False)
    
    output_path = f"{results_dir}/xy_slice.gif"
    anim.save(output_path, writer='pillow', fps=2)
    print(f"✓ Saved: {output_path}")
    plt.close()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python make_xy_gif.py <results_dir>")
        sys.exit(1)
    main(sys.argv[1])
