#!/usr/bin/env python3
"""
Generate radial profile GIF for IMBH-cloud simulation.

Usage:
    python make_radial_gif.py <results_dir>
    
Example:
    python make_radial_gif.py results/phase2.5_cooling_gamma53
"""
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

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

def compute_radial_profile(data, n_bins=30):
    """Compute radial profiles of T and n_H."""
    r = np.sqrt(data['x']**2 + data['y']**2 + data['z']**2)
    r_bins = np.linspace(0, 0.8, n_bins + 1)
    r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])
    
    T_profile = np.zeros(n_bins)
    n_profile = np.zeros(n_bins)
    
    for i in range(n_bins):
        mask = (r >= r_bins[i]) & (r < r_bins[i+1])
        if mask.sum() > 0:
            T_profile[i] = data['T'][mask].mean()
            n_profile[i] = data['n_H'][mask].mean()
        else:
            T_profile[i] = np.nan
            n_profile[i] = np.nan
    
    return r_centers, T_profile, n_profile

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
    
    # Create figure
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # Initialize plots
    data0 = read_snapshot(valid_snapshots[0])
    r_centers, T_profile, n_profile = compute_radial_profile(data0)
    
    line_T, = axes[0].plot(r_centers, T_profile, 'b-o', markersize=4)
    axes[0].axhline(7.0, color='green', linestyle='--', alpha=0.7, label='T_eq ≈ 7K')
    axes[0].axhline(5.0, color='orange', linestyle=':', alpha=0.7, label='T_floor = 5K')
    axes[0].set_xlabel('Radius [pc]')
    axes[0].set_ylabel('Temperature [K]')
    axes[0].set_title('Temperature Radial Profile')
    axes[0].set_ylim(4, 12)
    axes[0].legend(loc='upper right')
    axes[0].grid(True, alpha=0.3)
    
    line_n, = axes[1].semilogy(r_centers, n_profile, 'r-o', markersize=4)
    axes[1].set_xlabel('Radius [pc]')
    axes[1].set_ylabel('n_H [cm⁻³]')
    axes[1].set_title('Density Radial Profile')
    axes[1].set_ylim(10, 1e4)
    axes[1].grid(True, alpha=0.3)
    
    scatter = axes[2].scatter(data0['n_H'], data0['T'], s=1, alpha=0.3, c='blue')
    axes[2].axhline(7.0, color='green', linestyle='--', alpha=0.7)
    axes[2].axhline(5.0, color='orange', linestyle=':', alpha=0.7)
    axes[2].set_xlabel('n_H [cm⁻³]')
    axes[2].set_ylabel('Temperature [K]')
    axes[2].set_title('T-n Phase Diagram')
    axes[2].set_xscale('log')
    axes[2].set_xlim(10, 1e4)
    axes[2].set_ylim(4, 12)
    axes[2].grid(True, alpha=0.3)
    
    time_text = fig.suptitle(f't = 0.00 Myr', fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    def update(frame):
        data = read_snapshot(valid_snapshots[frame])
        r_centers, T_profile, n_profile = compute_radial_profile(data)
        
        line_T.set_ydata(T_profile)
        line_n.set_ydata(n_profile)
        
        # Update scatter (need to recreate)
        axes[2].clear()
        axes[2].scatter(data['n_H'], data['T'], s=1, alpha=0.3, c='blue')
        axes[2].axhline(7.0, color='green', linestyle='--', alpha=0.7)
        axes[2].axhline(5.0, color='orange', linestyle=':', alpha=0.7)
        axes[2].set_xlabel('n_H [cm⁻³]')
        axes[2].set_ylabel('Temperature [K]')
        axes[2].set_title('T-n Phase Diagram')
        axes[2].set_xscale('log')
        axes[2].set_xlim(10, 1e4)
        axes[2].set_ylim(4, 12)
        axes[2].grid(True, alpha=0.3)
        
        t_myr = data['t'] * T_TO_MYR
        time_text.set_text(f't = {t_myr:.2f} Myr  |  T_mean = {data["T"].mean():.1f} K')
        
        return line_T, line_n, time_text
    
    anim = FuncAnimation(fig, update, frames=len(valid_snapshots), 
                         interval=500, blit=False)
    
    output_path = f"{results_dir}/radial_profile.gif"
    anim.save(output_path, writer='pillow', fps=2)
    print(f"✓ Saved: {output_path}")
    plt.close()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python make_radial_gif.py <results_dir>")
        sys.exit(1)
    main(sys.argv[1])
