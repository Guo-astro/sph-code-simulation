#!/usr/bin/env python3
"""
Bonnor-Ebert Sphere Visualization

Generates:
1. Radial density profile comparison (SPH vs analytical BE)
2. Force balance profile (pressure vs gravity)
3. GIF animation of time evolution

Supports both CSV and HDF5 snapshot formats.
"""

import math
import os
import glob
import json
from collections import defaultdict

# Check for required packages
try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    print("ERROR: numpy not found. Install with: pip3 install numpy")
    exit(1)

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("ERROR: matplotlib not found. Install with: pip3 install matplotlib")
    exit(1)

try:
    import h5py
    HAS_H5PY = True
except ImportError:
    HAS_H5PY = False
    print("WARNING: h5py not found. HDF5 support disabled.")

try:
    from scipy.integrate import odeint
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


def load_snapshot(filename):
    """Load snapshot file (CSV or HDF5)"""
    if filename.endswith('.h5') or filename.endswith('.hdf5'):
        return load_hdf5(filename)
    else:
        return load_csv(filename)


def load_hdf5(filename):
    """Load HDF5 snapshot file"""
    if not HAS_H5PY:
        raise ImportError("h5py required for HDF5 files")

    data = {}
    with h5py.File(filename, 'r') as f:
        # Map HDF5 field names to expected names
        field_map = {
            'pos_x': 'particles/pos_x',
            'pos_y': 'particles/pos_y',
            'pos_z': 'particles/pos_z',
            'vel_x': 'particles/vel_x',
            'vel_y': 'particles/vel_y',
            'vel_z': 'particles/vel_z',
            'acc_x': 'particles/acc_x',
            'acc_y': 'particles/acc_y',
            'acc_z': 'particles/acc_z',
            'dens': 'particles/dens',
            'mass': 'particles/mass',
            'sml': 'particles/sml',
            'pres': 'particles/pres',
            'ene': 'particles/ene',
        }

        for key, path in field_map.items():
            if path in f:
                data[key] = np.array(f[path])

        # Get time from metadata if available
        if 'metadata' in f and 'time' in f['metadata']:
            data['_time'] = float(f['metadata/time'][()])
        else:
            data['_time'] = 0.0

    return data


def load_csv(filename):
    """Load CSV file, skipping comment lines"""
    import csv
    data = defaultdict(list)
    time = 0.0

    with open(filename, 'r') as f:
        headers = None
        for line in f:
            if line.startswith('# Time (code):'):
                try:
                    time = float(line.split(':')[1].strip())
                except:
                    pass
            elif not line.startswith('#'):
                headers = line.strip().split(',')
                break

        if not headers:
            raise ValueError("No header found in file")

        reader = csv.DictReader(f, fieldnames=headers)
        for row in reader:
            for col in headers:
                try:
                    data[col].append(float(row[col]))
                except (ValueError, KeyError):
                    data[col].append(0.0)

    result = {k: np.array(v) for k, v in data.items()}
    result['_time'] = time
    return result


def get_time_from_snapshot(filename):
    """Extract time from snapshot"""
    try:
        data = load_snapshot(filename)
        return data.get('_time', 0.0)
    except:
        return 0.0


def solve_lane_emden_scipy(xi_s, n_points=1000):
    """Solve Lane-Emden equation using scipy"""
    def lane_emden(y, xi):
        psi, dpsi = y
        if xi < 1e-10:
            return [dpsi, 0]
        d2psi = np.exp(-psi) - 2*dpsi/xi
        return [dpsi, d2psi]

    xi_arr = np.linspace(1e-6, xi_s, n_points)
    sol = odeint(lane_emden, [0, 0], xi_arr)
    psi_arr = sol[:, 0]
    dpsi_arr = sol[:, 1]
    rho_ratio = np.exp(-psi_arr)

    return xi_arr, rho_ratio, dpsi_arr


def solve_lane_emden_rk4(xi_s, n_points=5000):
    """Solve Lane-Emden equation using RK4 (fallback)"""
    dxi = xi_s / n_points
    psi = np.zeros(n_points + 1)
    dpsi = np.zeros(n_points + 1)
    xi = np.linspace(0, xi_s, n_points + 1)

    def f2(p, dp, x):
        if x > 1e-10:
            return math.exp(-p) - (2.0/x)*dp
        else:
            return math.exp(-p)

    for i in range(n_points):
        xi_i = max(xi[i], 1e-10)

        k1_p = dpsi[i]
        k1_dp = f2(psi[i], dpsi[i], xi_i)

        k2_p = dpsi[i] + 0.5*dxi*k1_dp
        k2_dp = f2(psi[i] + 0.5*dxi*k1_p, dpsi[i] + 0.5*dxi*k1_dp, xi_i + 0.5*dxi)

        k3_p = dpsi[i] + 0.5*dxi*k2_dp
        k3_dp = f2(psi[i] + 0.5*dxi*k2_p, dpsi[i] + 0.5*dxi*k2_dp, xi_i + 0.5*dxi)

        k4_p = dpsi[i] + dxi*k3_dp
        k4_dp = f2(psi[i] + dxi*k3_p, dpsi[i] + dxi*k3_dp, xi_i + dxi)

        psi[i+1] = psi[i] + dxi/6.0 * (k1_p + 2*k2_p + 2*k3_p + k4_p)
        dpsi[i+1] = dpsi[i] + dxi/6.0 * (k1_dp + 2*k2_dp + 2*k3_dp + k4_dp)

    rho_ratio = np.exp(-psi)
    return xi, rho_ratio, dpsi


def solve_lane_emden(xi_s, n_points=1000):
    """Solve Lane-Emden equation for isothermal case"""
    if HAS_SCIPY:
        return solve_lane_emden_scipy(xi_s, n_points)
    else:
        return solve_lane_emden_rk4(xi_s, n_points)


def load_config(result_dir):
    """Try to load config from result directory or parent"""
    config_paths = [
        os.path.join(result_dir, 'config.json'),
        os.path.join(os.path.dirname(result_dir), 'config.json'),
    ]

    # Default values
    config = {
        'xi_s': 6.0,
        'T_cloud': 10.0,
        'n_center': 100.0,
        'mu': 1.27,
    }

    for path in config_paths:
        if os.path.exists(path):
            try:
                with open(path, 'r') as f:
                    loaded = json.load(f)
                    config.update(loaded)
                print(f"  Loaded config from {path}")
                break
            except:
                pass

    return config


def compute_analytic_profile(config):
    """Compute analytical BE density profile"""
    # Physical constants (CGS)
    k_B = 1.380649e-16
    m_H = 1.6735575e-24
    G_cgs = 6.67430e-8
    pc_to_cm = 3.0857e18
    Msun_to_g = 1.989e33

    T_cloud = config.get('T_cloud', 10.0)
    n_center = config.get('n_center', 100.0)
    xi_s = config.get('xi_s', 6.0)
    mu = config.get('mu', 1.27)

    # Sound speed and scale length
    c_s = np.sqrt(k_B * T_cloud / (mu * m_H))
    rho_c_cgs = n_center * mu * m_H
    r_0_cgs = c_s / np.sqrt(4 * np.pi * G_cgs * rho_c_cgs)
    r_0_pc = r_0_cgs / pc_to_cm

    # Solve Lane-Emden
    xi_arr, rho_ratio, _ = solve_lane_emden(xi_s)

    # Convert to code units
    rho_c_code = rho_c_cgs * (pc_to_cm**3) / Msun_to_g
    r_analytic = xi_arr * r_0_pc
    rho_analytic = rho_c_code * rho_ratio

    return r_analytic, rho_analytic, r_0_pc * xi_s


def compute_profiles(filename, config=None):
    """Compute radial density profile from snapshot"""
    data = load_snapshot(filename)

    x = data['pos_x']
    y = data['pos_y']
    z = data.get('pos_z', np.zeros_like(x))
    dens = data['dens']
    time = data.get('_time', 0.0)

    # Compute radius
    radius = np.sqrt(x**2 + y**2 + z**2)

    return radius, dens, time


def plot_density_profile(filename, output_file='density_profile.png', config=None):
    """Plot radial density profile comparison"""
    print(f"Generating density profile plot...")

    if config is None:
        config = {}

    radius, dens, time = compute_profiles(filename, config)
    r_analytic, rho_analytic, R_cloud = compute_analytic_profile(config)

    fig, ax = plt.subplots(figsize=(10, 7))

    # Plot SPH scatter
    ax.scatter(radius, dens, s=0.5, alpha=0.3, c='blue', label='SPH')

    # Plot analytical BE
    ax.plot(r_analytic, rho_analytic, 'r-', linewidth=2, label='Analytical BE')

    ax.set_xlabel('Radius (pc)', fontsize=14)
    ax.set_ylabel('Density (code units)', fontsize=14)
    ax.set_title('Radial Density Profile: SPH vs Analytical Bonnor-Ebert', fontsize=16)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, R_cloud * 1.3)
    ax.set_ylim(0, rho_analytic.max() * 1.1)

    # Add text box
    textstr = f't = {time:.2f}\nR_cloud = {R_cloud:.2f} pc'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.95, 0.95, textstr, transform=ax.transAxes, fontsize=11,
            verticalalignment='top', horizontalalignment='right', bbox=props)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    plt.close()
    print(f"  Saved: {output_file}")


def create_animation(result_dir, output_file='be_evolution.gif', config=None):
    """Create GIF animation of density evolution"""
    print(f"Creating animation...")

    # Find all snapshots (CSV or HDF5)
    snapshots = sorted(glob.glob(os.path.join(result_dir, 'snapshot_*.h5')))
    if not snapshots:
        snapshots = sorted(glob.glob(os.path.join(result_dir, 'snapshot_*.hdf5')))
    if not snapshots:
        snapshots = sorted(glob.glob(os.path.join(result_dir, 'snapshot_*.csv')))

    if not snapshots:
        print(f"  No snapshots found in {result_dir}")
        return

    print(f"  Found {len(snapshots)} snapshots")

    if config is None:
        config = load_config(result_dir)

    # Compute analytical profile
    r_analytic, rho_analytic, R_cloud = compute_analytic_profile(config)

    # Generate frames
    frames = []
    for i, snap in enumerate(snapshots):
        radius, dens, time = compute_profiles(snap, config)

        fig, ax = plt.subplots(figsize=(10, 7))
        ax.scatter(radius, dens, s=0.5, alpha=0.3, c='blue', label='SPH')
        ax.plot(r_analytic, rho_analytic, 'r-', linewidth=2, label='Analytical BE')
        ax.set_xlabel('Radius (pc)', fontsize=12)
        ax.set_ylabel('Density (code units)', fontsize=12)
        ax.set_title(f'Relaxation - Frame {i}', fontsize=14)
        ax.set_xlim(0, R_cloud * 1.3)
        ax.set_ylim(0, max(rho_analytic.max(), np.percentile(dens, 99)) * 1.1)
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper right')

        frame_path = f'/tmp/be_frame_{i:04d}.png'
        plt.savefig(frame_path, dpi=100)
        plt.close()
        frames.append(frame_path)

        if (i + 1) % 10 == 0:
            print(f"  Frame {i+1}/{len(snapshots)}")

    # Create GIF using PIL
    try:
        from PIL import Image
        images = [Image.open(f) for f in frames]
        images[0].save(output_file, save_all=True, append_images=images[1:],
                       duration=150, loop=0)
        print(f"  Saved: {output_file}")
    except ImportError:
        print("  ERROR: PIL/Pillow not found. Install with: pip3 install pillow")

    # Cleanup temp files
    for f in frames:
        try:
            os.remove(f)
        except:
            pass


def main():
    import sys

    if len(sys.argv) < 2:
        print("Usage: python visualize_be_sphere.py <result_dir_or_snapshot>")
        print("Example: python visualize_be_sphere.py results/phase1_compact")
        return

    arg = sys.argv[1]

    if os.path.isdir(arg):
        result_dir = arg
        config = load_config(result_dir)

        # Find latest snapshot for static plot
        snapshots = sorted(glob.glob(os.path.join(arg, 'snapshot_*.h5')))
        if not snapshots:
            snapshots = sorted(glob.glob(os.path.join(arg, 'snapshot_*.csv')))

        if snapshots:
            snapshot = snapshots[-1]
            plot_density_profile(snapshot,
                                 os.path.join(result_dir, 'density_profile.png'),
                                 config)

        create_animation(result_dir,
                        os.path.join(result_dir, 'be_evolution.gif'),
                        config)
    else:
        snapshot = arg
        result_dir = os.path.dirname(arg)
        config = load_config(result_dir)
        plot_density_profile(snapshot, 'density_profile.png', config)

    print("\nDone!")


if __name__ == "__main__":
    main()
