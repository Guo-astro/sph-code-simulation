#!/usr/bin/env python3
"""
Density profile plot with analytical Bonnor-Ebert overlay and envelope particles.

Solves the isothermal Lane-Emden equation:
    (1/xi^2) d/dxi(xi^2 dpsi/dxi) = exp(-psi)

where:
    xi = r / r_0        (dimensionless radius)
    r_0 = c_s / sqrt(4*pi*G*rho_c)
    psi = -ln(rho/rho_c)
    rho = rho_c * exp(-psi)
"""
import sys
import os
import glob
import math
import json

# Physical constants in CGS
K_B_CGS = 1.380649e-16       # erg/K
M_PROTON_CGS = 1.6726219e-24 # g
MSUN_CGS = 1.989e33          # g
PC_CGS = 3.086e18            # cm
KMS_CGS = 1.0e5              # cm/s
G_CGS = 6.674e-8             # cm^3/g/s^2
G_CODE = 0.00430091          # pc^3 / (M_sun * Myr^2)


def solve_lane_emden_simple(xi_max, n_points=500):
    """
    Solve isothermal Lane-Emden equation using simple RK4.
    
    Returns:
        xi: Dimensionless radius array
        psi: Dimensionless potential array
        dpsi: dpsi/dxi array
    """
    # Initial conditions from series expansion near origin
    xi0 = 1e-6
    psi0 = xi0**2 / 6.0
    dpsi0 = xi0 / 3.0
    
    dxi = (xi_max - xi0) / n_points
    
    xi_arr = [xi0]
    psi_arr = [psi0]
    dpsi_arr = [dpsi0]
    
    xi = xi0
    psi = psi0
    dpsi = dpsi0
    
    for _ in range(n_points):
        # RK4 integration
        def f1(xi, psi, dpsi):
            return dpsi
        
        def f2(xi, psi, dpsi):
            if xi < 1e-10:
                return 0.0
            return math.exp(-psi) - 2.0 * dpsi / xi
        
        k1_psi = dxi * f1(xi, psi, dpsi)
        k1_dpsi = dxi * f2(xi, psi, dpsi)
        
        k2_psi = dxi * f1(xi + 0.5*dxi, psi + 0.5*k1_psi, dpsi + 0.5*k1_dpsi)
        k2_dpsi = dxi * f2(xi + 0.5*dxi, psi + 0.5*k1_psi, dpsi + 0.5*k1_dpsi)
        
        k3_psi = dxi * f1(xi + 0.5*dxi, psi + 0.5*k2_psi, dpsi + 0.5*k2_dpsi)
        k3_dpsi = dxi * f2(xi + 0.5*dxi, psi + 0.5*k2_psi, dpsi + 0.5*k2_dpsi)
        
        k4_psi = dxi * f1(xi + dxi, psi + k3_psi, dpsi + k3_dpsi)
        k4_dpsi = dxi * f2(xi + dxi, psi + k3_psi, dpsi + k3_dpsi)
        
        xi = xi + dxi
        psi = psi + (k1_psi + 2*k2_psi + 2*k3_psi + k4_psi) / 6.0
        dpsi = dpsi + (k1_dpsi + 2*k2_dpsi + 2*k3_dpsi + k4_dpsi) / 6.0
        
        xi_arr.append(xi)
        psi_arr.append(psi)
        dpsi_arr.append(dpsi)
    
    return xi_arr, psi_arr, dpsi_arr


def compute_be_profile(params, n_points=500):
    """
    Compute physical BE profile from parameters.
    
    Returns dict with:
        r: Radius array [pc]
        rho: Density array [M_sun/pc^3]
        n: Number density array [cm^-3]
    """
    mu = params.get('mu', 1.27)
    T = params.get('T_cloud', 7.0)
    n_center = params.get('n_center', 1800.0)
    xi_s = params.get('xi_s', 6.0)
    
    # Sound speed: c_s = sqrt(k_B * T / (mu * m_p))
    c_s_cgs = math.sqrt(K_B_CGS * T / (mu * M_PROTON_CGS))  # cm/s
    c_s = c_s_cgs / KMS_CGS  # km/s
    
    # Central density: n_center [cm^-3] -> rho_c [M_sun/pc^3]
    rho_c_cgs = n_center * mu * M_PROTON_CGS  # g/cm^3
    rho_c = rho_c_cgs * PC_CGS**3 / MSUN_CGS  # M_sun/pc^3
    
    # Scale length: r_0 = c_s / sqrt(4*pi*G*rho_c)
    r_0 = c_s / math.sqrt(4.0 * math.pi * G_CODE * rho_c)  # pc
    
    # Solve Lane-Emden
    xi_arr, psi_arr, dpsi_arr = solve_lane_emden_simple(xi_s * 1.1, n_points)
    
    # Physical quantities
    r = [x * r_0 for x in xi_arr]
    rho = [rho_c * math.exp(-p) for p in psi_arr]
    
    # Number density
    density_to_n = (MSUN_CGS / PC_CGS**3) / (mu * M_PROTON_CGS)
    n = [rh * density_to_n for rh in rho]
    
    # Enclosed mass: M(r) = 4*pi*rho_c*r_0^3 * xi^2 * dpsi
    M_enc = [4.0 * math.pi * rho_c * r_0**3 * xi**2 * dp 
             for xi, dp in zip(xi_arr, dpsi_arr)]
    
    # Truncate at xi_s
    mask_indices = [i for i, xi in enumerate(xi_arr) if xi <= xi_s]
    
    return {
        'r': [r[i] for i in mask_indices],
        'rho': [rho[i] for i in mask_indices],
        'n': [n[i] for i in mask_indices],
        'M_enc': [M_enc[i] for i in mask_indices],
        'r_0': r_0,
        'R_cloud': xi_s * r_0,
        'rho_c': rho_c,
        'n_center': n_center,
        'c_s': c_s,
        'xi_s': xi_s,
        'mu': mu,
        'T': T,
    }


def load_csv_data(filepath):
    """Load CSV data, separating real and ghost particles."""
    real_particles = {'x': [], 'y': [], 'z': [], 'rho': []}
    ghost_particles = {'x': [], 'y': [], 'z': [], 'rho': []}
    
    with open(filepath, 'r') as f:
        # Skip comment lines
        header = None
        for line in f:
            if not line.startswith('#'):
                header = line.strip().split(',')
                break
        
        if not header:
            return real_particles, ghost_particles
        
        # Find column indices
        try:
            x_idx = header.index('pos_x')
            y_idx = header.index('pos_y')
            z_idx = header.index('pos_z')
            rho_idx = header.index('dens')
            ghost_idx = header.index('is_ghost') if 'is_ghost' in header else None
        except ValueError:
            return real_particles, ghost_particles
        
        # Read data
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split(',')
            if len(parts) <= max(x_idx, y_idx, z_idx, rho_idx):
                continue
            
            x = float(parts[x_idx])
            y = float(parts[y_idx])
            z = float(parts[z_idx])
            rho = float(parts[rho_idx])
            
            is_ghost = int(parts[ghost_idx]) if ghost_idx is not None else 0
            
            if is_ghost == 1:
                ghost_particles['x'].append(x)
                ghost_particles['y'].append(y)
                ghost_particles['z'].append(z)
                ghost_particles['rho'].append(rho)
            else:
                real_particles['x'].append(x)
                real_particles['y'].append(y)
                real_particles['z'].append(z)
                real_particles['rho'].append(rho)
    
    return real_particles, ghost_particles


def bin_radial_profile(x, y, z, rho, n_bins=30):
    """Compute binned radial density profile."""
    radii = [math.sqrt(xi**2 + yi**2 + zi**2) for xi, yi, zi in zip(x, y, z)]
    
    if not radii:
        return [], [], []
    
    r_min, r_max = 0, max(radii) * 1.05
    bin_width = (r_max - r_min) / n_bins
    
    bins = [[] for _ in range(n_bins)]
    for r, rh in zip(radii, rho):
        idx = min(int((r - r_min) / bin_width), n_bins - 1)
        bins[idx].append(rh)
    
    bin_centers = []
    bin_means = []
    bin_stds = []
    
    for i, bin_data in enumerate(bins):
        if len(bin_data) > 0:
            r_center = r_min + (i + 0.5) * bin_width
            mean = sum(bin_data) / len(bin_data)
            
            if len(bin_data) > 1:
                variance = sum((x - mean)**2 for x in bin_data) / (len(bin_data) - 1)
                std = math.sqrt(variance)
            else:
                std = 0
            
            bin_centers.append(r_center)
            bin_means.append(mean)
            bin_stds.append(std)
    
    return bin_centers, bin_means, bin_stds


def main():
    # Find the results directory
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    results_dir = os.path.join(base_dir, "results", "phase1_relaxation")
    config_file = os.path.join(base_dir, "config", "presets", "phase1_relaxation.json")
    
    # Load config parameters
    params = {
        'mu': 1.27,
        'T_cloud': 7.0,
        'n_center': 1800.0,
        'xi_s': 6.0,
        'R_cloud': 0.75,
    }
    
    if os.path.exists(config_file):
        try:
            with open(config_file, 'r') as f:
                config = json.load(f)
            params.update({
                'mu': config.get('mu', params['mu']),
                'T_cloud': config.get('T_cloud', params['T_cloud']),
                'n_center': config.get('n_center', params['n_center']),
                'xi_s': config.get('xi_s', params['xi_s']),
                'R_cloud': config.get('R_cloud', params['R_cloud']),
            })
            print(f"Loaded parameters from config: mu={params['mu']}, T={params['T_cloud']}K, "
                  f"n_center={params['n_center']} cm^-3, xi_s={params['xi_s']}")
        except Exception as e:
            print(f"Warning: Could not load config: {e}")
    
    # Find CSV files
    csv_files = sorted(glob.glob(os.path.join(results_dir, "snapshot_*.csv")))
    if not csv_files:
        csv_files = sorted(glob.glob(os.path.join(results_dir, "*.csv")))
    
    if not csv_files:
        print(f"No CSV files found in {results_dir}")
        sys.exit(1)
    
    print(f"Found {len(csv_files)} CSV files")
    latest_file = csv_files[-1]
    print(f"Using latest: {latest_file}")
    
    # Load data
    real_data, ghost_data = load_csv_data(latest_file)
    print(f"Loaded {len(real_data['x'])} real particles, {len(ghost_data['x'])} ghost particles")
    
    # Compute analytical BE profile
    print("\nComputing analytical Bonnor-Ebert profile...")
    be = compute_be_profile(params, n_points=500)
    print(f"  Scale length r_0 = {be['r_0']:.4f} pc")
    print(f"  Cloud radius R = {be['R_cloud']:.4f} pc")
    print(f"  Central density = {be['rho_c']:.2f} M_sun/pc^3")
    print(f"  Sound speed = {be['c_s']:.4f} km/s")
    
    # Compute radial profiles
    real_r, real_rho, real_std = bin_radial_profile(
        real_data['x'], real_data['y'], real_data['z'], real_data['rho'], n_bins=30)
    ghost_r, ghost_rho, ghost_std = bin_radial_profile(
        ghost_data['x'], ghost_data['y'], ghost_data['z'], ghost_data['rho'], n_bins=20)
    
    # Calculate radii for scatter plots
    real_radii = [math.sqrt(x**2 + y**2 + z**2) 
                  for x, y, z in zip(real_data['x'], real_data['y'], real_data['z'])]
    ghost_radii = [math.sqrt(x**2 + y**2 + z**2) 
                   for x, y, z in zip(ghost_data['x'], ghost_data['y'], ghost_data['z'])]
    
    # Try matplotlib
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        
        fig = plt.figure(figsize=(16, 12))
        
        # ===== Plot 1: Density scatter with analytical overlay =====
        ax1 = fig.add_subplot(2, 2, 1)
        
        # Scatter plot of all particles (downsampled for visibility)
        step = max(1, len(real_radii) // 2000)
        ax1.scatter(real_radii[::step], [real_data['rho'][i] for i in range(0, len(real_radii), step)],
                   s=3, alpha=0.3, c='blue', label='SPH particles')
        
        # Ghost particles
        if ghost_radii:
            step_g = max(1, len(ghost_radii) // 1000)
            ax1.scatter(ghost_radii[::step_g], [ghost_data['rho'][i] for i in range(0, len(ghost_radii), step_g)],
                       s=3, alpha=0.3, c='red', label='Envelope (ghost)')
        
        # Analytical BE profile
        ax1.plot(be['r'], be['rho'], 'k-', lw=2.5, label=f'Analytical BE (ξₛ={be["xi_s"]})', zorder=10)
        
        ax1.axvline(be['R_cloud'], color='green', ls='--', alpha=0.7, 
                   label=f'R_cloud = {be["R_cloud"]:.2f} pc')
        
        ax1.set_xlabel('Radius [pc]', fontsize=12)
        ax1.set_ylabel(r'Density [M$_\odot$/pc$^3$]', fontsize=12)
        ax1.set_title('Density vs Radius (All Particles)', fontsize=14)
        ax1.legend(loc='upper right', fontsize=10)
        ax1.set_xlim(0, max(ghost_radii) * 1.1 if ghost_radii else be['R_cloud'] * 1.5)
        ax1.grid(True, alpha=0.3)
        ax1.set_yscale('log')
        
        # ===== Plot 2: Binned profile with error bars =====
        ax2 = fig.add_subplot(2, 2, 2)
        
        # Analytical profile
        ax2.plot(be['r'], be['rho'], 'k-', lw=2.5, label='Analytical BE', zorder=10)
        
        # Real particle profile with error bars
        if real_r:
            ax2.errorbar(real_r, real_rho, yerr=real_std, fmt='o', color='blue',
                        markersize=6, capsize=3, label='SPH (binned mean ± σ)')
        
        # Ghost particle profile
        if ghost_r:
            ax2.errorbar(ghost_r, ghost_rho, yerr=ghost_std, fmt='s', color='red',
                        markersize=5, capsize=3, alpha=0.7, label='Envelope (binned)')
        
        ax2.axvline(be['R_cloud'], color='green', ls='--', alpha=0.7)
        ax2.axhline(ghost_rho[0] if ghost_rho else params.get('n_edge', 162) * params['mu'] * M_PROTON_CGS * PC_CGS**3 / MSUN_CGS,
                   color='red', ls=':', alpha=0.5, label='P_ext target')
        
        ax2.set_xlabel('Radius [pc]', fontsize=12)
        ax2.set_ylabel(r'Density [M$_\odot$/pc$^3$]', fontsize=12)
        ax2.set_title('Radial Density Profile (Binned)', fontsize=14)
        ax2.legend(loc='upper right', fontsize=10)
        ax2.set_xlim(0, max(ghost_radii) * 1.1 if ghost_radii else be['R_cloud'] * 1.5)
        ax2.grid(True, alpha=0.3)
        ax2.set_yscale('log')
        
        # ===== Plot 3: Linear scale density (cloud region only) =====
        ax3 = fig.add_subplot(2, 2, 3)
        
        ax3.plot(be['r'], be['rho'], 'k-', lw=2.5, label='Analytical BE', zorder=10)
        
        if real_r:
            ax3.scatter(real_r, real_rho, s=50, c='blue', marker='o', 
                       label='SPH (binned)', zorder=5)
            ax3.fill_between(real_r, 
                            [m - s for m, s in zip(real_rho, real_std)],
                            [m + s for m, s in zip(real_rho, real_std)],
                            alpha=0.3, color='blue')
        
        ax3.axvline(be['R_cloud'], color='green', ls='--', alpha=0.7, 
                   label=f'R_cloud = {be["R_cloud"]:.2f} pc')
        ax3.axvline(be['r_0'], color='orange', ls=':', alpha=0.7,
                   label=f'r_0 = {be["r_0"]:.3f} pc')
        
        ax3.set_xlabel('Radius [pc]', fontsize=12)
        ax3.set_ylabel(r'Density [M$_\odot$/pc$^3$]', fontsize=12)
        ax3.set_title('Cloud Region (Linear Scale)', fontsize=14)
        ax3.legend(loc='upper right', fontsize=10)
        ax3.set_xlim(0, be['R_cloud'] * 1.2)
        ax3.set_ylim(0, be['rho_c'] * 1.2)
        ax3.grid(True, alpha=0.3)
        
        # ===== Plot 4: XY projection showing cloud and envelope =====
        ax4 = fig.add_subplot(2, 2, 4)
        
        # Filter to a thin slice in z
        z_slice = 0.1  # pc
        
        real_xy_x = [real_data['x'][i] for i in range(len(real_data['x'])) 
                     if abs(real_data['z'][i]) < z_slice]
        real_xy_y = [real_data['y'][i] for i in range(len(real_data['y'])) 
                     if abs(real_data['z'][i]) < z_slice]
        real_xy_rho = [real_data['rho'][i] for i in range(len(real_data['rho'])) 
                       if abs(real_data['z'][i]) < z_slice]
        
        ghost_xy_x = [ghost_data['x'][i] for i in range(len(ghost_data['x'])) 
                      if abs(ghost_data['z'][i]) < z_slice]
        ghost_xy_y = [ghost_data['y'][i] for i in range(len(ghost_data['y'])) 
                      if abs(ghost_data['z'][i]) < z_slice]
        
        # Plot ghost particles first (background)
        if ghost_xy_x:
            ax4.scatter(ghost_xy_x, ghost_xy_y, s=5, c='red', alpha=0.3, 
                       label=f'Envelope ({len(ghost_xy_x)} particles)')
        
        # Plot real particles with density coloring
        if real_xy_x:
            import matplotlib.colors as mcolors
            norm = mcolors.LogNorm(vmin=min(real_xy_rho), vmax=max(real_xy_rho))
            sc = ax4.scatter(real_xy_x, real_xy_y, s=10, c=real_xy_rho, 
                            cmap='viridis', norm=norm, alpha=0.7,
                            label=f'Cloud ({len(real_xy_x)} particles)')
            plt.colorbar(sc, ax=ax4, label=r'Density [M$_\odot$/pc$^3$]')
        
        # Draw R_cloud circle
        theta = [i * 2 * math.pi / 100 for i in range(101)]
        circle_x = [be['R_cloud'] * math.cos(t) for t in theta]
        circle_y = [be['R_cloud'] * math.sin(t) for t in theta]
        ax4.plot(circle_x, circle_y, 'g--', lw=2, label=f'R_cloud = {be["R_cloud"]:.2f} pc')
        
        ax4.set_xlabel('x [pc]', fontsize=12)
        ax4.set_ylabel('y [pc]', fontsize=12)
        ax4.set_title(f'XY Projection (|z| < {z_slice} pc slice)', fontsize=14)
        ax4.set_aspect('equal')
        ax4.legend(loc='upper right', fontsize=9)
        ax4.grid(True, alpha=0.3)
        
        # Set axis limits
        max_extent = max(ghost_radii) * 1.1 if ghost_radii else be['R_cloud'] * 1.5
        ax4.set_xlim(-max_extent, max_extent)
        ax4.set_ylim(-max_extent, max_extent)
        
        plt.tight_layout()
        
        # Add title
        fig.suptitle(f'Phase 1 Relaxation: Bonnor-Ebert Sphere (ξₛ={params["xi_s"]}, T={params["T_cloud"]}K)',
                    fontsize=16, y=1.02)
        
        output_path = os.path.join(base_dir, "results", "density_with_analytic.png")
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"\n✓ Plot saved to: {output_path}")
        
        # Also create a summary plot
        fig2, ax = plt.subplots(figsize=(10, 7))
        
        # Scatter all data
        ax.scatter(real_radii, real_data['rho'], s=2, alpha=0.2, c='blue', label='SPH Cloud')
        if ghost_radii:
            ax.scatter(ghost_radii, ghost_data['rho'], s=2, alpha=0.2, c='red', label='Envelope')
        
        # Analytical
        ax.plot(be['r'], be['rho'], 'k-', lw=3, label=f'Analytical BE (ξₛ={be["xi_s"]})', zorder=10)
        
        # Binned profile
        ax.plot(real_r, real_rho, 'b-', lw=2, marker='o', markersize=5, label='SPH binned mean')
        
        ax.axvline(be['R_cloud'], color='green', ls='--', lw=2, alpha=0.8)
        ax.text(be['R_cloud']*1.02, ax.get_ylim()[1]*0.8, f'R={be["R_cloud"]:.2f} pc', 
               fontsize=10, color='green')
        
        ax.set_xlabel('Radius [pc]', fontsize=14)
        ax.set_ylabel(r'Density [M$_\odot$/pc$^3$]', fontsize=14)
        ax.set_title(f'BE Sphere Density Profile: SPH vs Analytical\n'
                    f'(T={params["T_cloud"]}K, n_c={params["n_center"]} cm⁻³, ξₛ={params["xi_s"]})',
                    fontsize=14)
        ax.legend(loc='upper right', fontsize=11)
        ax.set_yscale('log')
        ax.grid(True, alpha=0.3)
        
        output_path2 = os.path.join(base_dir, "results", "density_profile_analytic.png")
        plt.savefig(output_path2, dpi=150, bbox_inches='tight')
        print(f"✓ Summary plot saved to: {output_path2}")
        
    except ImportError as e:
        print(f"\nMatplotlib not available: {e}")
        print("ASCII profile shown above.")
    except Exception as e:
        print(f"\nError creating plot: {e}")
        import traceback
        traceback.print_exc()


if __name__ == '__main__':
    main()
