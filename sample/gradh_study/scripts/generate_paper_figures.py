#!/usr/bin/env python3
"""
Generate Publication-Ready Figures for Grad-H Study Paper

Creates high-quality figures suitable for academic publication demonstrating
the importance of the grad-h correction term in SPH hydrostatic simulations.

Usage:
    python3 generate_paper_figures.py --results-dir <dir> --output <outdir>
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import glob
import os
import argparse

# Publication-quality settings
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'lines.linewidth': 1.5,
})

parser = argparse.ArgumentParser(description='Generate paper figures')
parser.add_argument('--results-dir', required=True, help='Base results directory')
parser.add_argument('--output', required=True, help='Output directory')
args = parser.parse_args()

results_dir = args.results_dir
output_dir = args.output

print('=' * 70)
print('Generating Publication-Ready Figures for Grad-H Study')
print('=' * 70)

os.makedirs(output_dir, exist_ok=True)

# Colors for publication (colorblind-safe)
COLOR_WITH_GRADH = '#0072B2'    # Blue
COLOR_NO_GRADH = '#D55E00'      # Vermillion


def load_snapshot(filename):
    """Load snapshot data"""
    data = {}
    with open(filename, 'r') as f:
        metadata = {}
        for line in f:
            if line.startswith('#'):
                if ':' in line:
                    key, value = line[1:].strip().split(':', 1)
                    metadata[key.strip()] = value.strip()
            else:
                break
        
        f.seek(0)
        lines = [ln for ln in f.readlines() if not ln.startswith('#')]
        header = lines[0].strip().split(',')
        
        for col_name in header:
            data[col_name] = []
        
        for line in lines[1:]:
            values = line.strip().split(',')
            for i, col_name in enumerate(header):
                try:
                    data[col_name].append(float(values[i]))
                except (ValueError, IndexError):
                    pass
    
    for key in data:
        data[key] = np.array(data[key])
    
    return data, metadata


def get_snapshot_list(directory):
    """Get sorted list of snapshot numbers"""
    snapshots = sorted(glob.glob(os.path.join(directory, "snapshot_*.csv")))
    numbers = []
    for snap in snapshots:
        num_str = os.path.basename(snap).replace('snapshot_', '').replace('.csv', '')
        try:
            numbers.append(int(num_str))
        except ValueError:
            pass
    return numbers


def compute_central_density(data, dim=2):
    """Compute central density"""
    x = data.get('pos_x', np.array([]))
    y = data.get('pos_y', np.array([]))
    z = data.get('pos_z', np.array([]))
    dens = data.get('dens', np.array([]))
    
    if len(x) == 0:
        return None
    
    if dim == 3 and len(z) > 0:
        r = np.sqrt(x**2 + y**2 + z**2)
    else:
        r = np.sqrt(x**2 + y**2)
    
    r_max = np.percentile(r, 90)
    central_mask = r < 0.1 * r_max
    
    if np.sum(central_mask) > 0:
        return np.mean(dens[central_mask])
    return None


# ================================================================================
# Figure 1: Two-panel comparison (particles + density evolution)
# ================================================================================
print("\nFigure 1: Main comparison figure...")

with_dir = os.path.join(results_dir, 'with_gradh')
no_dir = os.path.join(results_dir, 'no_gradh')

if os.path.exists(with_dir) and os.path.exists(no_dir):
    fig = plt.figure(figsize=(7.0, 4.5))  # Single column width
    gs = gridspec.GridSpec(2, 3, figure=fig, height_ratios=[1, 1], 
                          width_ratios=[1, 1, 1.2], wspace=0.3, hspace=0.35)
    
    # Get snapshot list
    snaps_with = get_snapshot_list(with_dir)
    snaps_no = get_snapshot_list(no_dir)
    common = sorted(set(snaps_with) & set(snaps_no))
    
    if common:
        # Select initial and final snapshots
        snap_init = common[0]
        snap_final = common[-1]
        
        # Top row: Initial state
        for col, (case_dir, label) in enumerate([(with_dir, 'With Grad-H'), (no_dir, 'Without Grad-H')]):
            ax = fig.add_subplot(gs[0, col])
            snap_file = os.path.join(case_dir, f"snapshot_{snap_init:04d}.csv")
            if os.path.exists(snap_file):
                data, meta = load_snapshot(snap_file)
                x, y, dens = data.get('pos_x', []), data.get('pos_y', []), data.get('dens', [])
                ax.scatter(x, y, c=dens, s=0.5, cmap='viridis', rasterized=True)
                ax.set_xlim(-0.55, 0.55)
                ax.set_ylim(-0.55, 0.55)
                ax.set_aspect('equal')
                ax.set_title(f'{label}', fontsize=10)
                if col == 0:
                    ax.set_ylabel('y', fontsize=10)
                    ax.text(-0.48, 0.42, f't = {float(meta.get("time", 0)):.1f}', fontsize=8,
                           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
                ax.set_xticklabels([])
        
        # Bottom row: Final state
        for col, (case_dir, label) in enumerate([(with_dir, 'With Grad-H'), (no_dir, 'Without Grad-H')]):
            ax = fig.add_subplot(gs[1, col])
            snap_file = os.path.join(case_dir, f"snapshot_{snap_final:04d}.csv")
            if os.path.exists(snap_file):
                data, meta = load_snapshot(snap_file)
                x, y, dens = data.get('pos_x', []), data.get('pos_y', []), data.get('dens', [])
                ax.scatter(x, y, c=dens, s=0.5, cmap='viridis', rasterized=True)
                ax.set_xlim(-0.55, 0.55)
                ax.set_ylim(-0.55, 0.55)
                ax.set_aspect('equal')
                ax.set_xlabel('x', fontsize=10)
                if col == 0:
                    ax.set_ylabel('y', fontsize=10)
                    ax.text(-0.48, 0.42, f't = {float(meta.get("time", 0)):.1f}', fontsize=8,
                           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        # Right column: Density evolution
        ax_dens = fig.add_subplot(gs[:, 2])
        
        for case_dir, label, color, ls in [(with_dir, r'$\Omega \neq 1$', COLOR_WITH_GRADH, '-'),
                                           (no_dir, r'$\Omega = 1$', COLOR_NO_GRADH, '--')]:
            times, central_dens = [], []
            for snap_num in common:
                snap_file = os.path.join(case_dir, f"snapshot_{snap_num:04d}.csv")
                if os.path.exists(snap_file):
                    data, meta = load_snapshot(snap_file)
                    t = float(meta.get('time', 0))
                    rho_c = compute_central_density(data)
                    if rho_c is not None:
                        times.append(t)
                        central_dens.append(rho_c)
            
            if times:
                times = np.array(times)
                central_dens = np.array(central_dens)
                if central_dens[0] > 0:
                    central_dens /= central_dens[0]
                ax_dens.plot(times, central_dens, ls, color=color, label=label, linewidth=1.5)
        
        ax_dens.axhline(y=1.0, color='gray', linestyle=':', linewidth=0.8, alpha=0.7)
        ax_dens.set_xlabel('Time', fontsize=10)
        ax_dens.set_ylabel(r'$\rho_c / \rho_{c,0}$', fontsize=10)
        ax_dens.legend(frameon=True, fancybox=False, edgecolor='gray')
        ax_dens.grid(True, alpha=0.3, linewidth=0.5)
        ax_dens.set_ylim(0, None)
    
    # Add panel labels
    for i, label in enumerate(['(a)', '(b)', '(c)', '(d)', '(e)']):
        pass  # Could add panel labels if needed
    
    plt.savefig(os.path.join(output_dir, 'fig1_gradh_comparison.pdf'), format='pdf')
    plt.savefig(os.path.join(output_dir, 'fig1_gradh_comparison.png'), format='png', dpi=300)
    plt.close()
    print("  Saved: fig1_gradh_comparison.pdf/png")
else:
    print("  Skipped: Periodic box results not found")

# ================================================================================
# Figure 2: Self-gravitating comparison (if available)
# ================================================================================
print("\nFigure 2: Self-gravitating comparison...")

sg_with_dir = os.path.join(results_dir, 'selfgrav_with_gradh')
sg_no_dir = os.path.join(results_dir, 'selfgrav_no_gradh')

if os.path.exists(sg_with_dir) and os.path.exists(sg_no_dir):
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0))
    
    snaps_with = get_snapshot_list(sg_with_dir)
    snaps_no = get_snapshot_list(sg_no_dir)
    common = sorted(set(snaps_with) & set(snaps_no))
    
    if common:
        # Panel (a): Central density evolution
        ax = axes[0]
        for case_dir, label, color, ls in [(sg_with_dir, r'With $\Omega$', COLOR_WITH_GRADH, '-'),
                                           (sg_no_dir, r'Without $\Omega$', COLOR_NO_GRADH, '--')]:
            times, central_dens = [], []
            for snap_num in common:
                snap_file = os.path.join(case_dir, f"snapshot_{snap_num:04d}.csv")
                if os.path.exists(snap_file):
                    data, meta = load_snapshot(snap_file)
                    t = float(meta.get('time', 0))
                    rho_c = compute_central_density(data, dim=3)
                    if rho_c is not None:
                        times.append(t)
                        central_dens.append(rho_c)
            
            if times:
                times = np.array(times)
                central_dens = np.array(central_dens)
                if central_dens[0] > 0:
                    central_dens /= central_dens[0]
                ax.plot(times, central_dens, ls, color=color, label=label)
        
        ax.axhline(y=1.0, color='gray', linestyle=':', linewidth=0.8)
        ax.set_xlabel('Time', fontsize=10)
        ax.set_ylabel(r'$\rho_c / \rho_{c,0}$', fontsize=10)
        ax.set_title('(a) Central density', fontsize=10)
        ax.legend(frameon=True, fancybox=False, edgecolor='gray')
        ax.grid(True, alpha=0.3, linewidth=0.5)
        ax.set_yscale('log')
        
        # Panel (b): Final density profile
        ax = axes[1]
        snap_final = common[-1]
        
        for case_dir, label, color, ls in [(sg_with_dir, r'With $\Omega$', COLOR_WITH_GRADH, '-'),
                                           (sg_no_dir, r'Without $\Omega$', COLOR_NO_GRADH, '--')]:
            snap_file = os.path.join(case_dir, f"snapshot_{snap_final:04d}.csv")
            if os.path.exists(snap_file):
                data, meta = load_snapshot(snap_file)
                x = data.get('pos_x', np.array([]))
                y = data.get('pos_y', np.array([]))
                z = data.get('pos_z', np.array([]))
                dens = data.get('dens', np.array([]))
                
                if len(x) > 0:
                    if len(z) > 0:
                        r = np.sqrt(x**2 + y**2 + z**2)
                    else:
                        r = np.sqrt(x**2 + y**2)
                    
                    r_max = np.percentile(r, 95)
                    r_bins = np.linspace(0, r_max, 40)
                    r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])
                    dens_profile = []
                    
                    for i in range(len(r_bins) - 1):
                        mask = (r >= r_bins[i]) & (r < r_bins[i+1])
                        if np.sum(mask) > 0:
                            dens_profile.append(np.mean(dens[mask]))
                        else:
                            dens_profile.append(np.nan)
                    
                    ax.plot(r_centers, dens_profile, ls, color=color, label=label)
        
        ax.set_xlabel('Radius', fontsize=10)
        ax.set_ylabel('Density', fontsize=10)
        ax.set_title(f'(b) Final profile (t = {float(meta.get("time", 0)):.1f})', fontsize=10)
        ax.legend(frameon=True, fancybox=False, edgecolor='gray')
        ax.grid(True, alpha=0.3, linewidth=0.5)
        ax.set_yscale('log')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'fig2_selfgrav_gradh.pdf'), format='pdf')
    plt.savefig(os.path.join(output_dir, 'fig2_selfgrav_gradh.png'), format='png', dpi=300)
    plt.close()
    print("  Saved: fig2_selfgrav_gradh.pdf/png")
else:
    print("  Skipped: Self-gravitating results not found")

print()
print('=' * 70)
print('Paper figures generated!')
print(f'Output directory: {output_dir}')
print('=' * 70)
