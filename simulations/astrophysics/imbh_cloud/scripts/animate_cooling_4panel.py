#!/usr/bin/env python3
"""
5-Panel Animation for Phase 2.5 Cooling Test

Creates animated GIF showing:
1. n-T phase diagram with KI2000 equilibrium curve
2. Radial density profile
3. Radial temperature profile
4. Radial internal energy profile
5. Radial cooling rate profile

Usage:
    python animate_cooling_4panel.py [results_dir] [output.gif]
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from pathlib import Path
import argparse
import sys

# Physical constants
K_B = 1.380649e-16  # Boltzmann constant [erg/K]
M_H = 1.6726e-24    # Hydrogen mass [g]
MU = 1.27           # Mean molecular weight for molecular gas
GAMMA = 1.0001      # Near-isothermal

# Unit conversions for IMBH cloud
DENSITY_TO_N_H = 31.9      # cm^-3 per code density unit
U_TO_CGS = 1e10            # erg/g per code energy unit
T_TO_CGS = 3.086e13        # seconds per code time unit

# KI2000 equilibrium curve data for N_H = 10^20 cm^-2
N_GRID_KI2000 = np.array([
    0.01, 0.0158, 0.0251, 0.0398, 0.0631, 0.1, 0.158, 0.251, 0.398, 0.631,
    1.0, 1.58, 2.51, 3.98, 6.31, 10.0, 15.8, 25.1, 39.8, 63.1,
    100.0, 158.0, 251.0, 398.0, 631.0, 1000.0, 1585.0, 2512.0, 3981.0, 6310.0,
    10000.0, 15849.0, 25119.0, 39811.0, 63096.0, 100000.0, 158489.0, 251189.0,
    398107.0, 630957.0, 1000000.0, 1584893.0, 2511886.0, 3981072.0, 6309573.0,
    10000000.0, 15848932.0, 25118864.0, 39810717.0, 63095734.0, 100000000.0,
    158489319.0, 251188643.0, 398107171.0, 631578947.0, 1000000000.0,
    1584893192.0, 2511886432.0, 3981071706.0, 6309573445.0, 10000000000.0,
    15848931925.0, 25118864315.0, 39810717055.0, 63095734448.0
])

T_EQ_GRID_KI2000 = np.array([
    5756.37, 5721.56, 5672.93, 5604.63, 5508.00, 5371.37, 5180.03, 4916.63,
    4561.24, 4094.75, 3503.19, 2784.86, 2017.36, 1349.45, 893.60, 601.08,
    413.12, 287.83, 201.29, 140.48, 97.57, 67.40, 46.27, 31.53, 21.28,
    14.20, 9.37, 6.12, 6.12, 6.33, 6.57, 6.84, 7.13, 7.47, 7.84,
    8.24, 8.70, 9.21, 9.77, 10.41, 11.12, 11.94, 12.87, 13.95, 15.20,
    16.69, 18.46, 20.59, 23.16, 26.28, 30.04, 34.58, 40.04, 46.46, 54.02,
    62.91, 73.35, 85.63, 100.08, 117.20, 137.39, 161.12, 189.21, 222.23, 260.81
])


def internal_energy_to_temperature(u_code):
    """Convert code internal energy to temperature [K]."""
    u_cgs = u_code * U_TO_CGS
    T = u_cgs * MU * M_H * (GAMMA - 1.0) / K_B
    return T


def get_ki2000_equilibrium_temp(n_H):
    """Get KI2000 equilibrium temperature for given density."""
    if n_H <= 0:
        return np.nan
    log_n = np.log10(n_H)
    log_n_grid = np.log10(N_GRID_KI2000)
    log_T_grid = np.log10(T_EQ_GRID_KI2000)
    return 10**np.interp(log_n, log_n_grid, log_T_grid)


def compute_cooling_rate(n_H, T):
    """
    Compute approximate cooling rate based on KI2000 relaxation.
    Returns net cooling rate in erg/g/s (negative = heating).
    """
    T_eq = np.array([get_ki2000_equilibrium_temp(n) for n in n_H])
    # Relaxation timescale (approximate)
    t_cool = 1e6 * 3.154e7  # 1 Myr in seconds
    # Net cooling = (T - T_eq) / t_cool * specific heat
    c_v = K_B / (MU * M_H * (GAMMA - 1.0))  # specific heat
    dT = T - T_eq
    cooling_rate = dT * c_v / t_cool  # erg/g/s
    return cooling_rate, T_eq


def load_snapshot(filepath):
    """Load a snapshot CSV file with comment lines."""
    # Find header line (first non-comment line)
    with open(filepath, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                header = line.strip().split(',')
                break
    
    # Load data skipping comment lines
    data_lines = []
    with open(filepath, 'r') as f:
        for line in f:
            if not line.startswith('#') and line.strip():
                # Skip header line too
                if line.strip().startswith('id,'):
                    continue
                data_lines.append(line.strip())
    
    # Parse numeric data
    arr = np.array([list(map(float, line.split(','))) for line in data_lines])
    return {col: arr[:, i] for i, col in enumerate(header)}


def load_energy_file(filepath):
    """Load energy.dat file."""
    return np.loadtxt(filepath, skiprows=1)


def scan_data_ranges(snapshots):
    """Pre-scan all snapshots to find proper data ranges."""
    print("Scanning snapshots for data ranges...")
    
    T_all, u_all, cooling_all, n_all, r_all = [], [], [], [], []
    
    for i, snapshot in enumerate(snapshots):
        if i % 10 == 0:
            print(f"  Scanning {i+1}/{len(snapshots)}...")
        
        data = load_snapshot(snapshot)
        
        # Filter for cloud particles (non-ghost)
        if 'is_ghost' in data:
            mask = data['is_ghost'] == 0
        elif 'particle_type' in data:
            mask = data['particle_type'] == 0
        else:
            mask = data['dens'] > 10
        
        # Compute radius
        x = data['pos_x'][mask]
        y = data['pos_y'][mask]
        z = data['pos_z'][mask]
        r = np.sqrt(x**2 + y**2 + z**2)
        
        # Get density and temperature
        dens = data['dens'][mask]
        n_H = dens * DENSITY_TO_N_H
        
        if 'ene' in data:
            u = data['ene'][mask]
        elif 'eng' in data:
            u = data['eng'][mask]
        else:
            u = np.full_like(dens, 454.97)
        
        T = internal_energy_to_temperature(u)
        
        # Compute cooling rate
        cooling_rate, _ = compute_cooling_rate(n_H, T)
        
        T_all.extend(T)
        u_all.extend(u)
        cooling_all.extend(cooling_rate)
        n_all.extend(n_H)
        r_all.extend(r)
    
    T_all = np.array(T_all)
    u_all = np.array(u_all)
    cooling_all = np.array(cooling_all)
    n_all = np.array(n_all)
    r_all = np.array(r_all)
    
    # Compute ranges with 10% padding
    def get_range(arr, padding=0.1):
        vmin, vmax = np.nanmin(arr), np.nanmax(arr)
        span = vmax - vmin
        if span == 0:
            span = abs(vmin) * 0.1 if vmin != 0 else 1.0
        return vmin - span * padding, vmax + span * padding
    
    ranges = {
        'T': get_range(T_all),
        'u': get_range(u_all),
        'cooling': get_range(cooling_all),
        'n': (max(1, np.nanmin(n_all) * 0.5), np.nanmax(n_all) * 2),
        'r': (0, np.nanmax(r_all) * 1.1),
    }
    
    print(f"  Temperature range: {ranges['T'][0]:.2f} - {ranges['T'][1]:.2f} K")
    print(f"  Energy range: {ranges['u'][0]:.2f} - {ranges['u'][1]:.2f} code units")
    print(f"  Cooling range: {ranges['cooling'][0]:.2e} - {ranges['cooling'][1]:.2e} erg/g/s")
    print(f"  Density range: {ranges['n'][0]:.1f} - {ranges['n'][1]:.1f} cm^-3")
    print(f"  Radius range: {ranges['r'][0]:.2f} - {ranges['r'][1]:.2f} pc")
    
    return ranges


def main():
    parser = argparse.ArgumentParser(description='4-panel cooling animation')
    parser.add_argument('results_dir', nargs='?',
                       default='simulations/astrophysics/imbh_cloud/results/phase2.5_cooling',
                       help='Results directory')
    parser.add_argument('output', nargs='?', default='be_cooling_4panel.gif',
                       help='Output GIF filename')
    args = parser.parse_args()
    
    results_dir = Path(args.results_dir)
    output_path = results_dir / args.output
    
    # Find all snapshots
    snapshots = sorted(results_dir.glob('snapshot_*.csv'))
    if len(snapshots) < 2:
        print(f"Need at least 2 snapshots, found {len(snapshots)}")
        sys.exit(1)
    
    print(f"Found {len(snapshots)} snapshots")
    
    # Pre-scan all snapshots to find data ranges
    data_ranges = scan_data_ranges(snapshots)
    
    # Load energy data for timing
    energy_file = results_dir / 'energy.dat'
    if energy_file.exists():
        energy_data = load_energy_file(energy_file)
        times = energy_data[:, 0] if len(energy_data.shape) > 1 else [0]
    else:
        times = np.arange(len(snapshots)) * 2.0  # Assume 2 Myr intervals
    
    # Set up figure with 2x3 layout
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    ax_phase = axes[0, 0]
    ax_dens = axes[0, 1]
    ax_temp = axes[0, 2]
    ax_energy = axes[1, 0]
    ax_cooling = axes[1, 1]
    ax_empty = axes[1, 2]
    ax_empty.axis('off')  # Hide the 6th panel
    
    fig.suptitle('Bonnor-Ebert Sphere: KI2000 Cooling Test', fontsize=14, fontweight='bold')
    
    # Pre-plot KI2000 equilibrium curve
    n_range = np.logspace(-1, 5, 500)
    T_eq_curve = np.array([get_ki2000_equilibrium_temp(n) for n in n_range])
    
    # Initialize plots with detected ranges
    ax_phase.loglog(n_range, T_eq_curve, 'k-', lw=2, label='KI2000 T_eq', zorder=1)
    scatter_phase = ax_phase.scatter([], [], c=[], cmap='viridis', s=1, alpha=0.5, zorder=2)
    ax_phase.set_xlim(data_ranges['n'][0], data_ranges['n'][1])
    ax_phase.set_ylim(max(1, data_ranges['T'][0] * 0.5), data_ranges['T'][1] * 2)
    ax_phase.set_xlabel(r'$n_H$ [cm$^{-3}$]', fontsize=12)
    ax_phase.set_ylabel('T [K]', fontsize=12)
    ax_phase.legend(loc='upper right')
    ax_phase.set_title('Phase Diagram')
    ax_phase.axhline(7, color='red', ls='--', alpha=0.5, label='T=7K')
    ax_phase.axvline(1800, color='blue', ls='--', alpha=0.5, label='n=1800')
    
    # Density profile panel
    line_dens, = ax_dens.plot([], [], 'c.', ms=2, alpha=0.3)
    ax_dens.set_xlim(0, data_ranges['r'][1])
    ax_dens.set_ylim(data_ranges['n'][0], data_ranges['n'][1])
    ax_dens.set_yscale('log')
    ax_dens.set_xlabel('r [pc]', fontsize=12)
    ax_dens.set_ylabel(r'$n_H$ [cm$^{-3}$]', fontsize=12)
    ax_dens.set_title('Radial Density Profile')
    ax_dens.axhline(1800, color='red', ls='--', alpha=0.5, label=r'$n_{center}$=1800')
    ax_dens.axhline(162, color='blue', ls='--', alpha=0.5, label=r'$n_{edge}$=162')
    ax_dens.legend(loc='upper right')
    
    line_temp, = ax_temp.plot([], [], 'b.', ms=2, alpha=0.3)
    line_temp_eq, = ax_temp.plot([], [], 'r-', lw=2, label=r'$T_{eq}$ (KI2000)')
    ax_temp.set_xlim(0, data_ranges['r'][1])
    ax_temp.set_ylim(data_ranges['T'][0], data_ranges['T'][1])
    ax_temp.set_xlabel('r [pc]', fontsize=12)
    ax_temp.set_ylabel('T [K]', fontsize=12)
    ax_temp.legend(loc='upper right')
    ax_temp.set_title('Radial Temperature Profile')
    ax_temp.axhline(7, color='gray', ls='--', alpha=0.5, label='T=7K')
    ax_temp.legend(loc='upper right')
    
    line_energy, = ax_energy.plot([], [], 'g.', ms=2, alpha=0.3)
    ax_energy.set_xlim(0, data_ranges['r'][1])
    ax_energy.set_ylim(data_ranges['u'][0], data_ranges['u'][1])
    ax_energy.set_xlabel('r [pc]', fontsize=12)
    ax_energy.set_ylabel('u [code units]', fontsize=12)
    ax_energy.set_title('Radial Internal Energy Profile')
    ax_energy.axhline(454.97, color='red', ls='--', alpha=0.5, label='u(T=7K)=454.97')
    ax_energy.legend(loc='upper right')
    
    line_cooling, = ax_cooling.plot([], [], 'm.', ms=2, alpha=0.3)
    ax_cooling.set_xlim(0, data_ranges['r'][1])
    ax_cooling.set_ylim(data_ranges['cooling'][0], data_ranges['cooling'][1])
    ax_cooling.set_xlabel('r [pc]', fontsize=12)
    ax_cooling.set_ylabel(r'$\dot{u}$ [erg/g/s]', fontsize=12)
    ax_cooling.set_title('Net Cooling Rate (+ = cooling, - = heating)')
    ax_cooling.axhline(0, color='gray', ls='--', alpha=0.5, label='Equilibrium')
    ax_cooling.legend(loc='upper right')
    
    time_text = fig.text(0.5, 0.02, '', ha='center', fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.93, bottom=0.08)
    
    def init():
        scatter_phase.set_offsets(np.c_[[], []])
        line_dens.set_data([], [])
        line_temp.set_data([], [])
        line_temp_eq.set_data([], [])
        line_energy.set_data([], [])
        line_cooling.set_data([], [])
        time_text.set_text('')
        return scatter_phase, line_dens, line_temp, line_temp_eq, line_energy, line_cooling, time_text
    
    def update(frame):
        snapshot = snapshots[frame]
        data = load_snapshot(snapshot)
        
        # Get time
        t_myr = times[frame] if frame < len(times) else frame * 2.0
        
        # Filter for cloud particles (non-ghost)
        if 'particle_type' in data:
            mask = data['particle_type'] == 0
        else:
            # Assume particles with density > 10 are cloud
            mask = data['dens'] > 10
        
        # Compute radius
        x = data['pos_x'][mask]
        y = data['pos_y'][mask]
        z = data['pos_z'][mask]
        r = np.sqrt(x**2 + y**2 + z**2)
        
        # Get density and temperature
        dens = data['dens'][mask]
        n_H = dens * DENSITY_TO_N_H
        
        if 'eng' in data:
            u = data['eng'][mask]
        else:
            u = np.full_like(dens, 454.97)
        
        T = internal_energy_to_temperature(u)
        
        # Compute cooling rate
        cooling_rate, T_eq = compute_cooling_rate(n_H, T)
        
        # Update phase diagram
        scatter_phase.set_offsets(np.c_[n_H, T])
        scatter_phase.set_array(r)
        
        # Sort by radius for profiles
        sort_idx = np.argsort(r)
        r_sorted = r[sort_idx]
        T_sorted = T[sort_idx]
        u_sorted = u[sort_idx]
        cooling_sorted = cooling_rate[sort_idx]
        n_sorted = n_H[sort_idx]
        T_eq_sorted = T_eq[sort_idx]
        
        # Update density profile
        line_dens.set_data(r_sorted, n_sorted)
        
        # Update temperature profile
        line_temp.set_data(r_sorted, T_sorted)
        line_temp_eq.set_data(r_sorted, T_eq_sorted)
        
        # Update energy profile
        line_energy.set_data(r_sorted, u_sorted)
        
        # Update cooling rate profile
        line_cooling.set_data(r_sorted, cooling_sorted)
        
        # Update time text
        time_text.set_text(f't = {t_myr:.1f} Myr')
        
        # Keep fixed y-limits for consistent visualization
        # (ranges set during initialization)
        
        print(f"  Frame {frame+1}/{len(snapshots)}: t={t_myr:.1f} Myr, "
              f"T_mean={T_sorted.mean():.2f} K, u_mean={u_sorted.mean():.2f}, n_center={n_sorted[0]:.0f}")
        
        return scatter_phase, line_dens, line_temp, line_temp_eq, line_energy, line_cooling, time_text
    
    print("Creating animation...")
    anim = FuncAnimation(fig, update, frames=len(snapshots),
                        init_func=init, blit=False, interval=500)
    
    print(f"Saving to {output_path}...")
    anim.save(output_path, writer='pillow', fps=2, dpi=100)
    print(f"✓ Animation saved: {output_path}")
    
    plt.close()


if __name__ == '__main__':
    main()
