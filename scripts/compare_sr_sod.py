#!/usr/bin/env python3
"""
Compare two SR Sod simulations side-by-side
Shows before/after density blip fix
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def load_snapshot(filename):
    """Load CSV snapshot"""
    data = []
    time = 0.0

    with open(filename, 'r') as f:
        for line in f:
            if '# Time (code):' in line:
                time = float(line.split(':')[1].strip())
            if line.startswith('#') or line.strip().startswith('id'):
                continue
            parts = line.strip().split(',')
            if len(parts) >= 15:
                try:
                    data.append([float(x) for x in parts])
                except ValueError:
                    continue

    data = np.array(data)
    result = {
        'time': time,
        'pos_x': data[:, 1],
        'vel_x': data[:, 4],
        'nu': data[:, 10],
        'dens': data[:, 11],
        'pres': data[:, 12],
    }

    idx = np.argsort(result['pos_x'])
    for key in ['pos_x', 'vel_x', 'nu', 'dens', 'pres']:
        result[key] = result[key][idx]

    # Compute lab-frame density
    c = 1.0
    v = result['vel_x']
    gamma = 1.0 / np.sqrt(1.0 - v**2 / c**2)
    result['N'] = gamma * result['dens']
    result['gamma'] = gamma

    return result

def plot_comparison(data1, data2, label1, label2, output_file):
    """Create side-by-side comparison plot"""
    fig, axes = plt.subplots(3, 2, figsize=(14, 12))

    t1, t2 = data1['time'], data2['time']

    # Lab-frame density N
    ax = axes[0, 0]
    ax.plot(data1['pos_x'], data1['N'], 'r-', linewidth=2, label=label1)
    ax.set_xlabel('Position x')
    ax.set_ylabel('Lab-frame Density N')
    ax.set_title(f'{label1} (t={t1:.3f})')
    ax.legend()
    ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    ax.plot(data2['pos_x'], data2['N'], 'b-', linewidth=2, label=label2)
    ax.set_xlabel('Position x')
    ax.set_ylabel('Lab-frame Density N')
    ax.set_title(f'{label2} (t={t2:.3f})')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Pressure
    ax = axes[1, 0]
    ax.plot(data1['pos_x'], data1['pres'], 'r-', linewidth=2)
    ax.set_xlabel('Position x')
    ax.set_ylabel('Pressure P')
    ax.set_title(f'Pressure: {label1}')
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    ax.plot(data2['pos_x'], data2['pres'], 'b-', linewidth=2)
    ax.set_xlabel('Position x')
    ax.set_ylabel('Pressure P')
    ax.set_title(f'Pressure: {label2}')
    ax.grid(True, alpha=0.3)

    # Velocity
    ax = axes[2, 0]
    ax.plot(data1['pos_x'], data1['vel_x'], 'r-', linewidth=2)
    ax.axhline(0, color='k', linestyle='--', alpha=0.3)
    ax.set_xlabel('Position x')
    ax.set_ylabel('Velocity v')
    ax.set_title(f'Velocity: {label1}')
    ax.grid(True, alpha=0.3)

    ax = axes[2, 1]
    ax.plot(data2['pos_x'], data2['vel_x'], 'b-', linewidth=2)
    ax.axhline(0, color='k', linestyle='--', alpha=0.3)
    ax.set_xlabel('Position x')
    ax.set_ylabel('Velocity v')
    ax.set_title(f'Velocity: {label2}')
    ax.grid(True, alpha=0.3)

    plt.suptitle('SR Sod Comparison: Before vs After', fontsize=14, y=0.995)
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Saved comparison plot: {output_file}")

def main():
    if len(sys.argv) < 3:
        print("Usage: python compare_sr_sod.py <snapshot1.csv> <snapshot2.csv> [output.png]")
        print("\nExample:")
        print("  python compare_sr_sod.py \\")
        print("    sample/sr_sod/results/sharp/snapshot_0131.csv \\")
        print("    sample/sr_sod/results/sharp_fixed/snapshot_0175.csv \\")
        print("    comparison.png")
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]

    if len(sys.argv) > 3:
        output_file = sys.argv[3]
    else:
        output_file = "sr_sod_comparison.png"

    print(f"Loading {file1}...")
    data1 = load_snapshot(file1)

    print(f"Loading {file2}...")
    data2 = load_snapshot(file2)

    # Extract labels from paths
    label1 = "Before (C_cd=0.2 default)"
    label2 = "After (C_cd=1.0 fixed)"

    if "sharp/" in file1 and "sharp_fixed" not in file1:
        label1 = "Before: 100p, C_cd=0.2, 1st order"
    if "sharp_fixed/" in file2:
        label2 = "After: 400p, C_cd=1.0, 2nd order"

    plot_comparison(data1, data2, label1, label2, output_file)

if __name__ == "__main__":
    main()
