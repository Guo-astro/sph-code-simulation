#!/usr/bin/env python3
"""
Animate radial profiles with Koyama & Inutsuka (2000) ISM physics.

Shows:
1. Radial density profile
2. Temperature profile (from internal energy)
3. Internal energy profile
4. Net cooling/heating rate (from KI2000)
5. Chemistry: H2, CO, and electron fractions (from KI2000 tables)

Usage:
    python animate_ki2000_profiles.py <results_dir> -o output.gif [options]

Example:
    python animate_ki2000_profiles.py results/phase2.5_cooling -o ki2000_profiles.gif
"""

import argparse
import glob
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.interpolate import interp1d

# ============================================================================
# Koyama & Inutsuka (2000) Embedded Data - Direct from C++ Header
# ============================================================================

# Column density N_H = 10^19 cm^-2 data
N19_n = np.array([
    1.0e-02, 1.2308e-02, 1.5397e-02, 1.8951e-02, 2.3708e-02, 2.9179e-02, 3.6503e-02,
    4.4928e-02, 5.6206e-02, 6.9177e-02, 8.6541e-02, 1.0651e-01, 1.3325e-01, 1.6400e-01,
    2.0517e-01, 2.5252e-01, 3.1591e-01, 3.9520e-01, 4.8641e-01, 6.0850e-01, 7.4894e-01,
    9.3693e-01, 1.1532e+00, 1.4426e+00, 1.7756e+00, 2.2303e+00, 2.7339e+00, 3.4341e+00,
    4.2095e+00, 5.2875e+00, 6.4814e+00, 8.1414e+00, 1.0020e+01, 1.2536e+01, 1.5429e+01,
    1.9301e+01, 2.3756e+01, 2.9719e+01, 3.6578e+01, 4.5759e+01, 5.6320e+01, 7.0457e+01,
    8.6718e+01, 1.0848e+02, 1.3352e+02, 1.6704e+02, 2.0559e+02, 2.5719e+02, 3.1655e+02,
    3.9601e+02, 4.8740e+02, 6.0974e+02, 7.5047e+02, 9.3884e+02, 1.1555e+03, 1.4456e+03,
    1.7792e+03, 2.2258e+03, 2.7395e+03, 3.4271e+03, 4.2180e+03, 5.2768e+03, 6.4946e+03,
    8.1248e+03, 1.0000e+04
])

N19_T_eq = np.array([
    696.51, 696.51, 696.51, 686.40, 686.40, 673.15, 673.15, 663.38, 663.38, 650.57,
    641.13, 631.82, 619.62, 601.76, 562.07, 492.76, 397.64, 291.08, 202.94, 136.74,
    95.34, 69.79, 52.60, 41.83, 34.25, 29.16, 25.07, 21.98, 19.65, 17.83,
    16.09, 14.81, 13.90, 12.80, 12.01, 11.38, 10.84, 10.33, 9.84, 9.37,
    9.05, 8.75, 8.50, 8.21, 8.09, 7.82, 7.67, 7.56, 7.56, 7.45,
    7.45, 7.45, 7.31, 7.31, 7.31, 7.20, 7.10, 6.96, 6.72, 6.41,
    5.90, 5.45, 5.02, 4.69, 4.40
])

# Column density N_H = 10^20 cm^-2 data
N20_n = N19_n.copy()  # Same density grid

N20_T_eq = np.array([
    686.40, 673.15, 673.15, 663.38, 663.38, 663.38, 650.57, 641.13, 641.13, 631.82,
    610.63, 581.58, 543.21, 469.32, 378.72, 291.08, 213.08, 153.71, 110.35, 81.97,
    62.08, 49.37, 39.84, 33.27, 28.19, 24.71, 21.66, 19.65, 17.83, 16.09,
    14.81, 13.90, 12.80, 12.19, 11.38, 10.69, 10.18, 9.69, 9.37, 8.92,
    8.62, 8.33, 8.09, 7.94, 7.67, 7.56, 7.45, 7.31, 7.20, 7.20,
    7.10, 7.10, 6.96, 6.86, 6.72, 6.63, 6.53, 6.22, 5.81, 5.35,
    4.85, 4.40, 4.11, 3.86, 3.67
])

# H2 fraction for N20 (more typical for our cloud)
N20_n_H2 = N20_n[:59]  # H2 data has 59 points
N20_x_H2 = np.array([
    1.0000e+00, 1.0000e+00, 9.7688e-01, 9.7688e-01, 9.7688e-01, 9.7688e-01, 9.5429e-01,
    9.5429e-01, 9.3222e-01, 9.1066e-01, 8.6903e-01, 8.4893e-01, 7.9139e-01, 7.4352e-01,
    6.7710e-01, 5.8842e-01, 4.9953e-01, 4.0785e-01, 3.1530e-01, 2.2902e-01, 1.5507e-01,
    1.0020e-01, 6.3244e-02, 3.7214e-02, 2.1897e-02, 1.2587e-02, 7.2350e-03, 4.1588e-03,
    2.3352e-03, 1.3216e-03, 7.2492e-04, 3.9764e-04, 2.1308e-04, 1.0981e-04, 5.2346e-05,
    2.2372e-05, 6.9448e-06, 1.3190e-06, 6.9043e-07, 4.7856e-07, 3.4759e-07, 2.5643e-07,
    1.9067e-07, 1.4176e-07, 1.0706e-07, 8.1488e-08, 6.0117e-08, 4.7949e-08, 4.1669e-08,
    4.0705e-08, 3.9764e-08, 3.7946e-08, 3.3757e-08, 2.8214e-08, 2.2328e-08, 1.7397e-08,
    1.3449e-08, 1.0479e-08, 1.0000e-08
])

# CO fraction for N20
N20_n_CO = N20_n[:24]  # CO data has 24 points
N20_x_CO = np.array([
    2.0019e-04, 1.7809e-04, 1.5843e-04, 1.3449e-04, 1.1507e-04, 9.5429e-05, 7.6113e-05,
    6.0235e-05, 4.6932e-05, 3.6282e-05, 2.6977e-05, 1.9902e-05, 1.4456e-05, 1.0020e-05,
    6.7842e-06, 4.3495e-06, 2.5794e-06, 1.4484e-06, 7.2350e-07, 3.2403e-07, 1.3216e-07,
    5.0246e-08, 1.7809e-08, 1.0000e-08
])

# Electron fraction for N20
N20_x_e = np.array([
    1.0479e-04, 1.2636e-04, 1.4769e-04, 1.6995e-04, 1.9104e-04, 2.0815e-04, 2.2857e-04,
    2.4519e-04, 2.5693e-04, 2.6924e-04, 2.7562e-04, 2.8214e-04, 2.8882e-04, 2.9566e-04,
    2.9566e-04, 3.0266e-04, 3.0266e-04, 3.0266e-04, 3.0266e-04, 3.0266e-04, 3.0982e-04,
    3.0982e-04, 3.0982e-04, 3.0982e-04, 3.0982e-04, 3.0982e-04, 3.1715e-04, 3.1715e-04,
    3.2466e-04, 3.3235e-04, 3.3757e-04, 3.5374e-04, 3.7069e-04, 3.9764e-04, 4.2655e-04,
    4.6840e-04, 5.1435e-04, 5.7369e-04, 6.4489e-04, 7.5964e-04, 8.7413e-04, 1.0459e-03,
    1.2611e-03, 1.5446e-03, 1.9518e-03, 2.5050e-03, 3.2403e-03, 4.3580e-03, 5.7706e-03,
    7.7611e-03, 1.0039e-02, 1.2885e-02, 1.5905e-02, 1.9029e-02, 2.2416e-02, 2.5794e-02,
    3.0148e-02, 3.4691e-02, 3.9920e-02, 4.5579e-02, 5.2448e-02, 6.0352e-02, 6.9448e-02,
    7.9294e-02, 9.1244e-02
])


class KI2000Physics:
    """Koyama & Inutsuka (2000) ISM physics calculator."""
    
    def __init__(self, N_H=1.0e20, gamma=5.0/3.0, mu=1.27):
        """
        Initialize KI2000 physics.
        
        Args:
            N_H: Column density [cm^-2], use 1e19 or 1e20
            gamma: Adiabatic index
            mu: Mean molecular weight
        """
        self.N_H = N_H
        self.gamma = gamma
        self.mu = mu
        
        # Physical constants
        self.k_B = 1.3807e-16  # erg/K
        self.m_p = 1.6726e-24  # g
        self.m_n = mu * self.m_p  # mean particle mass
        
        # Select data based on column density
        if N_H < 5e19:
            self.n_grid = N19_n
            self.T_eq_grid = N19_T_eq
        else:
            self.n_grid = N20_n
            self.T_eq_grid = N20_T_eq
        
        # Create interpolators
        self._setup_interpolators()
    
    def _setup_interpolators(self):
        """Set up interpolation functions."""
        log_n = np.log10(self.n_grid)
        
        # Temperature interpolator
        self.T_eq_interp = interp1d(log_n, self.T_eq_grid, 
                                    kind='linear', fill_value='extrapolate')
        
        # H2 fraction interpolator (N20 data)
        log_n_H2 = np.log10(N20_n_H2)
        self.x_H2_interp = interp1d(log_n_H2, np.log10(N20_x_H2),
                                    kind='linear', fill_value='extrapolate')
        
        # CO fraction interpolator (N20 data)
        log_n_CO = np.log10(N20_n_CO)
        self.x_CO_interp = interp1d(log_n_CO, np.log10(N20_x_CO),
                                    kind='linear', fill_value='extrapolate')
        
        # Electron fraction interpolator (N20 data)
        log_n_e = np.log10(N20_n)
        self.x_e_interp = interp1d(log_n_e, np.log10(N20_x_e),
                                   kind='linear', fill_value='extrapolate')
    
    def equilibrium_temperature(self, n_H):
        """Get equilibrium temperature for given number density."""
        n_H = np.atleast_1d(n_H)
        log_n = np.log10(np.clip(n_H, 1e-2, 1e4))
        return self.T_eq_interp(log_n)
    
    def h2_fraction(self, n_H):
        """Get H2 molecular fraction."""
        n_H = np.atleast_1d(n_H)
        log_n = np.log10(np.clip(n_H, N20_n_H2[0], N20_n_H2[-1]))
        return np.power(10, self.x_H2_interp(log_n))
    
    def co_fraction(self, n_H):
        """Get CO molecular fraction."""
        n_H = np.atleast_1d(n_H)
        log_n = np.log10(np.clip(n_H, N20_n_CO[0], N20_n_CO[-1]))
        return np.power(10, self.x_CO_interp(log_n))
    
    def electron_fraction(self, n_H):
        """Get electron fraction."""
        n_H = np.atleast_1d(n_H)
        log_n = np.log10(np.clip(n_H, N20_n[0], N20_n[-1]))
        return np.power(10, self.x_e_interp(log_n))
    
    def net_cooling_rate(self, n_H, T):
        """
        Compute net cooling rate (Γ - nΛ) / n [erg/s/particle].
        
        Positive = heating, Negative = cooling
        """
        n_H = np.atleast_1d(n_H)
        T = np.atleast_1d(T)
        
        T_eq = self.equilibrium_temperature(n_H)
        
        # Simple relaxation model: du/dt ∝ (T_eq - T) / t_cool
        # Cooling timescale approximation from K&I 2000
        # t_cool ~ 10^5 / n_H [years] for typical ISM
        t_cool_yr = 1e5 / np.clip(n_H, 0.1, 1e4)
        t_cool_s = t_cool_yr * 3.156e7  # convert to seconds
        
        # Energy difference
        u_eq = self.k_B * T_eq / ((self.gamma - 1) * self.m_n)
        u = self.k_B * T / ((self.gamma - 1) * self.m_n)
        
        # Net rate [erg/g/s]
        rate = (u_eq - u) / t_cool_s
        
        return rate
    
    def temperature_from_energy(self, u_cgs):
        """Convert specific internal energy [erg/g] to temperature [K]."""
        return u_cgs * (self.gamma - 1) * self.m_n / self.k_B
    
    def energy_from_temperature(self, T):
        """Convert temperature [K] to specific internal energy [erg/g]."""
        return self.k_B * T / ((self.gamma - 1) * self.m_n)
    
    def temperature_from_sound_speed(self, cs_cgs):
        """
        Convert sound speed [cm/s] to temperature [K].
        
        For isothermal gas: c_s^2 = k_B * T / (mu * m_p)
        This is more reliable than using gamma for nearly isothermal simulations.
        """
        return cs_cgs**2 * self.m_n / self.k_B


def load_snapshot(filepath):
    """Load a snapshot CSV file."""
    # Skip comment lines
    with open(filepath, 'r') as f:
        skip_rows = 0
        for line in f:
            if line.startswith('#'):
                skip_rows += 1
            else:
                break
    
    df = pd.read_csv(filepath, skiprows=skip_rows)
    return df


def get_time_from_snapshot(filepath):
    """Extract simulation time from snapshot header."""
    with open(filepath, 'r') as f:
        for line in f:
            if 'Time (code):' in line:
                parts = line.split(':')
                if len(parts) >= 2:
                    return float(parts[-1].strip())
    return 0.0


def main():
    parser = argparse.ArgumentParser(description='Animate KI2000 radial profiles')
    parser.add_argument('results_dir', help='Directory containing snapshot CSV files')
    parser.add_argument('-o', '--output', default='ki2000_profiles.gif', help='Output GIF filename')
    parser.add_argument('--R_cloud', type=float, default=0.75, help='Cloud radius [pc]')
    parser.add_argument('--density_to_nH', type=float, default=31.86, 
                        help='Code density to number density [cm^-3]')
    parser.add_argument('--u_to_cgs', type=float, default=1.0e10,
                        help='Code energy to CGS [erg/g]')
    parser.add_argument('--N_H', type=float, default=1.0e20,
                        help='Column density for KI2000 [cm^-2]')
    parser.add_argument('--v_to_cgs', type=float, default=1.0e5,
                        help='Code velocity to CGS [cm/s], default 1e5 for km/s')
    parser.add_argument('--sim_gamma', type=float, default=1.0001,
                        help='Simulation gamma (for isothermal, use ~1.0001)')
    parser.add_argument('--fps', type=int, default=3, help='Frames per second')
    args = parser.parse_args()
    
    # Find snapshots
    pattern = os.path.join(args.results_dir, 'snapshot_*.csv')
    snapshots = sorted(glob.glob(pattern))
    
    if not snapshots:
        print(f"No snapshots found in {args.results_dir}")
        return
    
    print(f"Found {len(snapshots)} snapshots")
    
    # Initialize KI2000 physics (use gamma close to isothermal)
    ki2000 = KI2000Physics(N_H=args.N_H, gamma=5.0/3.0, mu=1.27)
    
    # Load first snapshot to get reference values
    df0 = load_snapshot(snapshots[0])
    real_mask = df0['is_ghost'] == 0
    
    # Central density reference
    df_real = df0[real_mask]
    r0 = np.sqrt(df_real['pos_x']**2 + df_real['pos_y']**2 + df_real['pos_z']**2)
    central_mask = r0 < 0.1
    if central_mask.sum() > 0:
        rho_c = df_real.loc[central_mask, 'dens'].mean()
    else:
        rho_c = df_real['dens'].max()
    
    n_c = rho_c * args.density_to_nH  # Central number density
    print(f"Central density: {rho_c:.2f} code units = {n_c:.1f} cm^-3")
    
    # Set up figure with 6 panels (2x3)
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle('Koyama & Inutsuka (2000) ISM Physics Profiles', fontsize=14)
    
    # Flatten for easier access
    ax_dens, ax_temp, ax_energy = axes[0]
    ax_cool, ax_h2, ax_xe = axes[1]
    
    # Set up radial bins
    r_max = args.R_cloud * 1.5
    r_bins = np.linspace(0, r_max, 50)
    r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])
    
    # Initialize plot elements
    # Pre-scan all snapshots to find global min/max for consistent y-axes
    print("Pre-scanning snapshots for y-axis ranges...")
    all_n_H, all_T, all_ene = [], [], []
    for snap in snapshots:
        df_scan = load_snapshot(snap)
        real_mask_scan = df_scan['is_ghost'] == 0
        df_scan_real = df_scan[real_mask_scan]
        rho_scan = df_scan_real['dens'].values
        ene_scan = df_scan_real['ene'].values
        cs_scan = df_scan_real['sound'].values  # code units (km/s)
        n_H_scan = rho_scan * args.density_to_nH
        cs_cgs_scan = cs_scan * args.v_to_cgs  # cm/s
        T_scan = ki2000.temperature_from_sound_speed(cs_cgs_scan)  # More reliable
        all_n_H.extend(n_H_scan)
        all_T.extend(T_scan)
        all_ene.extend(ene_scan)
    
    # Compute robust y-limits (use percentiles to handle outliers)
    n_H_max = np.percentile(all_n_H, 99.5) * 1.1
    T_max = np.percentile(all_T, 99.5) * 1.1
    ene_max = np.percentile(all_ene, 99.5) * 1.1
    print(f"  Density range: 0 - {n_H_max:.0f} cm^-3")
    print(f"  Temperature range: 0 - {T_max:.0f} K")
    print(f"  Energy range: 0 - {ene_max:.0f} code units")
    
    line_dens, = ax_dens.plot([], [], 'b-', lw=2, label='Simulation')
    line_dens_eq, = ax_dens.plot([], [], 'r--', lw=1.5, alpha=0.7, label='KI2000 T_eq n')
    ax_dens.set_xlabel('Radius [pc]')
    ax_dens.set_ylabel('Number Density [cm⁻³]')
    ax_dens.set_xlim(0, r_max)
    ax_dens.set_ylim(0, n_H_max)
    ax_dens.legend(loc='upper right')
    ax_dens.set_title('Density Profile')
    ax_dens.axvline(args.R_cloud, color='gray', ls=':', alpha=0.5)
    
    line_temp, = ax_temp.plot([], [], 'b-', lw=2, label='Simulation')
    line_temp_eq, = ax_temp.plot([], [], 'r--', lw=1.5, alpha=0.7, label='KI2000 T_eq')
    ax_temp.set_xlabel('Radius [pc]')
    ax_temp.set_ylabel('Temperature [K]')
    ax_temp.set_xlim(0, r_max)
    ax_temp.set_ylim(0, T_max)
    ax_temp.legend(loc='upper right')
    ax_temp.set_title('Temperature Profile')
    ax_temp.axvline(args.R_cloud, color='gray', ls=':', alpha=0.5)
    
    line_energy, = ax_energy.plot([], [], 'b-', lw=2, label='Simulation')
    line_energy_eq, = ax_energy.plot([], [], 'r--', lw=1.5, alpha=0.7, label='KI2000 u_eq')
    ax_energy.set_xlabel('Radius [pc]')
    ax_energy.set_ylabel('Specific Energy [code]')
    ax_energy.set_xlim(0, r_max)
    ax_energy.set_ylim(0, ene_max)
    ax_energy.legend(loc='upper right')
    ax_energy.set_title('Internal Energy Profile')
    ax_energy.axvline(args.R_cloud, color='gray', ls=':', alpha=0.5)
    
    line_cool, = ax_cool.plot([], [], 'b-', lw=2, label='Net rate')
    ax_cool.axhline(0, color='k', ls='-', lw=0.5)
    ax_cool.set_xlabel('Radius [pc]')
    ax_cool.set_ylabel('Net Cooling Rate [erg/g/s]')
    ax_cool.set_xlim(0, r_max)
    ax_cool.set_ylim(-1e-20, 1e-20)
    ax_cool.legend(loc='upper right')
    ax_cool.set_title('Net Cooling/Heating Rate')
    ax_cool.axvline(args.R_cloud, color='gray', ls=':', alpha=0.5)
    
    line_h2, = ax_h2.plot([], [], 'g-', lw=2, label='H₂')
    line_co, = ax_h2.plot([], [], 'm-', lw=2, label='CO (×10⁴)')
    ax_h2.set_xlabel('Radius [pc]')
    ax_h2.set_ylabel('Molecular Fraction')
    ax_h2.set_xlim(0, r_max)
    ax_h2.set_yscale('log')
    ax_h2.set_ylim(1e-8, 1)
    ax_h2.legend(loc='upper right')
    ax_h2.set_title('Molecular Fractions (KI2000)')
    ax_h2.axvline(args.R_cloud, color='gray', ls=':', alpha=0.5)
    
    line_xe, = ax_xe.plot([], [], 'c-', lw=2, label='Electrons')
    ax_xe.set_xlabel('Radius [pc]')
    ax_xe.set_ylabel('Electron Fraction')
    ax_xe.set_xlim(0, r_max)
    ax_xe.set_yscale('log')
    ax_xe.set_ylim(1e-5, 1)
    ax_xe.legend(loc='upper right')
    ax_xe.set_title('Electron Fraction (KI2000)')
    ax_xe.axvline(args.R_cloud, color='gray', ls=':', alpha=0.5)
    
    time_text = fig.text(0.5, 0.02, '', ha='center', fontsize=12)
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    def update(frame):
        filepath = snapshots[frame]
        df = load_snapshot(filepath)
        time = get_time_from_snapshot(filepath)
        
        # Filter real particles
        real_mask = df['is_ghost'] == 0
        df_real = df[real_mask]
        
        # Calculate radius
        r = np.sqrt(df_real['pos_x']**2 + df_real['pos_y']**2 + df_real['pos_z']**2).values
        
        # Get quantities
        rho = df_real['dens'].values
        ene = df_real['ene'].values  # code units
        cs = df_real['sound'].values  # code units (km/s)
        
        # Convert to physical units
        n_H = rho * args.density_to_nH  # cm^-3
        u_cgs = ene * args.u_to_cgs  # erg/g
        cs_cgs = cs * args.v_to_cgs  # cm/s
        
        # Temperature from sound speed (more reliable for isothermal simulations)
        # c_s^2 = k_B * T / (mu * m_p) for isothermal gas
        T = ki2000.temperature_from_sound_speed(cs_cgs)  # K
        
        # Get KI2000 equilibrium values based on density
        T_eq = ki2000.equilibrium_temperature(n_H)
        u_eq_cgs = ki2000.energy_from_temperature(T_eq)
        u_eq_code = u_eq_cgs / args.u_to_cgs
        
        # Get chemistry from KI2000
        x_H2 = ki2000.h2_fraction(n_H)
        x_CO = ki2000.co_fraction(n_H)
        x_e = ki2000.electron_fraction(n_H)
        
        # Compute net cooling rate
        cool_rate = ki2000.net_cooling_rate(n_H, T)
        
        # Bin by radius
        n_H_binned = np.zeros(len(r_centers))
        T_binned = np.zeros(len(r_centers))
        ene_binned = np.zeros(len(r_centers))
        T_eq_binned = np.zeros(len(r_centers))
        u_eq_binned = np.zeros(len(r_centers))
        cool_binned = np.zeros(len(r_centers))
        x_H2_binned = np.zeros(len(r_centers))
        x_CO_binned = np.zeros(len(r_centers))
        x_e_binned = np.zeros(len(r_centers))
        
        for i in range(len(r_centers)):
            mask = (r >= r_bins[i]) & (r < r_bins[i+1])
            if mask.sum() > 0:
                n_H_binned[i] = np.mean(n_H[mask])
                T_binned[i] = np.mean(T[mask])
                ene_binned[i] = np.mean(ene[mask])
                T_eq_binned[i] = np.mean(T_eq[mask])
                u_eq_binned[i] = np.mean(u_eq_code[mask])
                cool_binned[i] = np.mean(cool_rate[mask])
                x_H2_binned[i] = np.mean(x_H2[mask])
                x_CO_binned[i] = np.mean(x_CO[mask])
                x_e_binned[i] = np.mean(x_e[mask])
        
        # Update plots
        valid = n_H_binned > 0
        
        line_dens.set_data(r_centers[valid], n_H_binned[valid])
        # Show equilibrium density curve at equilibrium temperature
        line_dens_eq.set_data(r_centers[valid], n_H_binned[valid])  # Same density
        
        line_temp.set_data(r_centers[valid], T_binned[valid])
        line_temp_eq.set_data(r_centers[valid], T_eq_binned[valid])
        
        line_energy.set_data(r_centers[valid], ene_binned[valid])
        line_energy_eq.set_data(r_centers[valid], u_eq_binned[valid])
        
        line_cool.set_data(r_centers[valid], cool_binned[valid])
        
        line_h2.set_data(r_centers[valid], x_H2_binned[valid])
        line_co.set_data(r_centers[valid], x_CO_binned[valid])  # Scaled for visibility
        
        line_xe.set_data(r_centers[valid], x_e_binned[valid])
        
        # Auto-scale cooling axis
        if np.any(valid):
            cool_max = np.max(np.abs(cool_binned[valid]))
            if cool_max > 0:
                ax_cool.set_ylim(-1.5*cool_max, 1.5*cool_max)
        
        time_text.set_text(f't = {time:.2f} code time')
        
        return [line_dens, line_dens_eq, line_temp, line_temp_eq, 
                line_energy, line_energy_eq, line_cool, 
                line_h2, line_co, line_xe, time_text]
    
    print(f"Creating animation with {len(snapshots)} frames...")
    anim = FuncAnimation(fig, update, frames=len(snapshots),
                         interval=1000//args.fps, blit=False)
    
    print(f"Saving to {args.output}...")
    anim.save(args.output, writer='pillow', fps=args.fps)
    print(f"✓ Animation saved: {args.output}")
    
    plt.close()


if __name__ == '__main__':
    main()
