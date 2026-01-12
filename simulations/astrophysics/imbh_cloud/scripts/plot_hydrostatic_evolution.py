#!/usr/bin/env python3
"""
Plot density profile evolution during hydrostatic test.
Monitors cloud stability by comparing SPH density to analytical BE profile.
"""

import os
import sys
import glob
import numpy as np

# Bonnor-Ebert parameters from config
MU = 1.27
T_CLOUD = 7.0  # K
N_CENTER = 1800.0  # cm^-3
XI_S = 6.0
G_CODE = 0.00430091

# Physical constants
K_B_CGS = 1.380649e-16  # erg/K
M_PROTON_CGS = 1.6726219e-24  # g
MSUN_CGS = 1.989e33  # g
PC_CGS = 3.086e18  # cm
KMS_CGS = 1.0e5  # cm/s

def solve_lane_emden(xi_max, n_points=1000):
    """Solve isothermal Lane-Emden equation."""
    xi = np.zeros(n_points)
    psi = np.zeros(n_points)
    dpsi = np.zeros(n_points)
    
    dxi = xi_max / (n_points - 1)
    xi[0] = 1e-6
    psi[0] = xi[0]**2 / 6.0
    dpsi[0] = xi[0] / 3.0
    
    for i in range(1, n_points):
        x, y, z = xi[i-1], psi[i-1], dpsi[i-1]
        
        def f2(x_, y_, z_):
            if x_ < 1e-10:
                return 0.0
            return np.exp(-y_) - 2.0 * z_ / x_
        
        k1_y = dxi * z
        k1_z = dxi * f2(x, y, z)
        k2_y = dxi * (z + 0.5*k1_z)
        k2_z = dxi * f2(x + 0.5*dxi, y + 0.5*k1_y, z + 0.5*k1_z)
        k3_y = dxi * (z + 0.5*k2_z)
        k3_z = dxi * f2(x + 0.5*dxi, y + 0.5*k2_y, z + 0.5*k2_z)
        k4_y = dxi * (z + k3_z)
        k4_z = dxi * f2(x + dxi, y + k3_y, z + k3_z)
        
        xi[i] = x + dxi
        psi[i] = y + (k1_y + 2*k2_y + 2*k3_y + k4_y) / 6.0
        dpsi[i] = z + (k1_z + 2*k2_z + 2*k3_z + k4_z) / 6.0
    
    return xi, psi, dpsi

def compute_be_profile():
    """Compute analytical Bonnor-Ebert profile."""
    # Sound speed
    c_s_sq_cgs = K_B_CGS * T_CLOUD / (MU * M_PROTON_CGS)
    c_s = np.sqrt(c_s_sq_cgs) / KMS_CGS  # km/s
    
    # Density conversion
    density_code_to_cgs = MSUN_CGS / PC_CGS**3
    density_to_n = density_code_to_cgs / (MU * M_PROTON_CGS)
    rho_c = N_CENTER / density_to_n  # code units
    
    # Scale length
    r_0 = c_s / np.sqrt(4 * np.pi * G_CODE * rho_c)
    R_cloud = XI_S * r_0
    
    # Solve Lane-Emden
    xi_arr, psi_arr, _ = solve_lane_emden(XI_S * 1.1)
    
    # Physical radius and density
    r_arr = xi_arr * r_0
    rho_arr = rho_c * np.exp(-psi_arr)
    
    return r_arr, rho_arr, R_cloud, rho_c

def load_snapshot(filepath):
    """Load snapshot CSV file."""
    # Read header to get column names (skip comment lines)
    with open(filepath, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                header = line.strip().split(',')
                break
    
    # Load data, skipping comment lines
    data = np.genfromtxt(filepath, delimiter=',', names=header, skip_header=0, comments='#')
    
    # Filter real particles (not ghost)
    if 'is_ghost' in data.dtype.names:
        mask = data['is_ghost'] == 0
        data = data[mask]
    
    return data

def plot_evolution(results_dir, output_dir=None):
    """Plot density profile evolution."""
    if output_dir is None:
        output_dir = results_dir
    
    # Find all snapshots
    snapshots = sorted(glob.glob(os.path.join(results_dir, 'snapshot_*.csv')))
    if not snapshots:
        print(f"No snapshots found in {results_dir}")
        return
    
    print(f"Found {len(snapshots)} snapshots")
    
    # Get analytical profile
    r_be, rho_be, R_cloud, rho_c = compute_be_profile()
    
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        
        # Create figure with subplots
        n_plots = min(len(snapshots), 6)  # Max 6 panels
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        axes = axes.flatten()
        
        # Select snapshots to plot (evenly spaced)
        indices = np.linspace(0, len(snapshots)-1, n_plots, dtype=int)
        
        for idx, (ax, snap_idx) in enumerate(zip(axes, indices)):
            filepath = snapshots[snap_idx]
            data = load_snapshot(filepath)
            
            # Compute radius for each particle (numpy structured array)
            if 'pos_x' in data.dtype.names:
                x = data['pos_x']
                y = data['pos_y']
                z = data['pos_z'] if 'pos_z' in data.dtype.names else np.zeros_like(x)
            else:
                x = data['x']
                y = data['y']
                z = data['z'] if 'z' in data.dtype.names else np.zeros_like(x)
            
            r = np.sqrt(x**2 + y**2 + z**2)
            rho = data['dens'] if 'dens' in data.dtype.names else data['density']
            
            # Bin data for profile
            r_bins = np.linspace(0, R_cloud * 1.2, 50)
            r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])
            rho_mean = np.zeros(len(r_centers))
            rho_std = np.zeros(len(r_centers))
            
            for i in range(len(r_centers)):
                mask = (r >= r_bins[i]) & (r < r_bins[i+1])
                if np.sum(mask) > 0:
                    rho_mean[i] = np.mean(rho[mask])
                    rho_std[i] = np.std(rho[mask])
            
            # Extract time from filename
            snap_num = int(os.path.basename(filepath).split('_')[1].split('.')[0])
            time = snap_num * 0.05  # outputTime from config
            
            # Plot
            ax.scatter(r, rho, s=1, alpha=0.3, c='blue', label='SPH')
            ax.plot(r_be, rho_be, 'r-', lw=2, label='Analytical BE')
            ax.axvline(R_cloud, color='green', ls='--', alpha=0.5, label=f'R={R_cloud:.2f} pc')
            
            ax.set_xlabel('r [pc]')
            ax.set_ylabel('ρ [M☉/pc³]')
            ax.set_title(f't = {time:.2f} Myr')
            ax.set_xlim(0, R_cloud * 1.3)
            ax.set_ylim(0, rho_c * 1.2)
            if idx == 0:
                ax.legend(loc='upper right', fontsize=8)
        
        plt.suptitle('Hydrostatic Test: Density Profile Evolution', fontsize=14)
        plt.tight_layout()
        
        outfile = os.path.join(output_dir, 'hydrostatic_evolution.png')
        plt.savefig(outfile, dpi=150)
        print(f"✓ Saved: {outfile}")
        plt.close()
        
        # Also create a single comparison plot (initial vs final)
        fig, ax = plt.subplots(figsize=(10, 7))
        
        # Initial
        data0 = load_snapshot(snapshots[0])
        if 'pos_x' in data0.dtype.names:
            r0 = np.sqrt(data0['pos_x']**2 + data0['pos_y']**2 + data0['pos_z']**2)
        else:
            r0 = np.sqrt(data0['x']**2 + data0['y']**2 + data0['z']**2)
        rho0 = data0['dens'] if 'dens' in data0.dtype.names else data0['density']
        ax.scatter(r0, rho0, s=2, alpha=0.3, c='blue', label='t=0 (Initial)')
        
        # Final
        data_f = load_snapshot(snapshots[-1])
        if 'pos_x' in data_f.dtype.names:
            rf = np.sqrt(data_f['pos_x']**2 + data_f['pos_y']**2 + data_f['pos_z']**2)
        else:
            rf = np.sqrt(data_f['x']**2 + data_f['y']**2 + data_f['z']**2)
        rhof = data_f['dens'] if 'dens' in data_f.dtype.names else data_f['density']
        ax.scatter(rf, rhof, s=2, alpha=0.3, c='orange', label=f't={len(snapshots)*0.05:.1f} Myr (Final)')
        
        # Analytical
        ax.plot(r_be, rho_be, 'r-', lw=2, label='Analytical BE')
        ax.axvline(R_cloud, color='green', ls='--', alpha=0.7, label=f'R_cloud={R_cloud:.2f} pc')
        
        ax.set_xlabel('r [pc]', fontsize=12)
        ax.set_ylabel('ρ [M☉/pc³]', fontsize=12)
        ax.set_title('Hydrostatic Test: Initial vs Final Density Profile', fontsize=14)
        ax.set_xlim(0, R_cloud * 1.3)
        ax.legend(loc='upper right')
        ax.grid(True, alpha=0.3)
        
        outfile2 = os.path.join(output_dir, 'hydrostatic_comparison.png')
        plt.savefig(outfile2, dpi=150)
        print(f"✓ Saved: {outfile2}")
        plt.close()
        
    except ImportError:
        print("Matplotlib not available, printing ASCII summary")
        for filepath in snapshots[::max(1, len(snapshots)//5)]:
            data = load_snapshot(filepath)
            snap_num = int(os.path.basename(filepath).split('_')[1].split('.')[0])
            rho = data['dens'] if 'dens' in data.dtype.names else data['density']
            print(f"Snapshot {snap_num}: rho_mean={np.mean(rho):.2f}, rho_max={np.max(rho):.2f}")

if __name__ == '__main__':
    if len(sys.argv) > 1:
        results_dir = sys.argv[1]
    else:
        results_dir = 'simulations/astrophysics/imbh_cloud/results/phase2_hydrostatic'
    
    plot_evolution(results_dir)
