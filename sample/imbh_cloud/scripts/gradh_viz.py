#!/usr/bin/env python3
"""
Generate grad-h comparison plots for IMBH 3D tests.
Usage: python3 gradh_viz.py <results_dir> <output_dir>
"""

import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 gradh_viz.py <results_dir> <output_dir>")
        sys.exit(1)
    
    results_base = sys.argv[1]
    output_dir = sys.argv[2]
    os.makedirs(output_dir, exist_ok=True)

    methods = {
        'gsph_gradh': {'label': 'GSPH + grad-h', 'color': '#0066CC', 'linestyle': '-'},
        'gsph_nogradh': {'label': 'GSPH - grad-h', 'color': '#CC0000', 'linestyle': '--'},
        'gsph_csmooth2': {'label': 'GSPH + C_smooth=2', 'color': '#009933', 'linestyle': '-.'},
    }

    # Read energy data
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    for method, props in methods.items():
        energy_file = os.path.join(results_base, method, 'energy.dat')
        if os.path.exists(energy_file):
            try:
                data = np.loadtxt(energy_file, skiprows=1)
                if data.ndim >= 2 and len(data) > 0:
                    t = data[:, 0]
                    E_kin = data[:, 1]
                    E_th = data[:, 2]
                    E_pot = data[:, 3]
                    E_total = data[:, 4] if data.shape[1] > 4 else E_kin + E_th + E_pot
                    
                    # Total energy
                    axes[0, 0].plot(t, E_total, color=props['color'], linestyle=props['linestyle'],
                                   label=props['label'], linewidth=2)
                    # Kinetic energy
                    axes[0, 1].plot(t, E_kin, color=props['color'], linestyle=props['linestyle'],
                                   label=props['label'], linewidth=2)
                    # Energy drift
                    E0 = E_total[0]
                    drift = (E_total - E0) / np.abs(E0) * 100
                    axes[1, 0].plot(t, drift, color=props['color'], linestyle=props['linestyle'],
                                   label=props['label'], linewidth=2)
                    # Potential energy
                    axes[1, 1].plot(t, E_pot, color=props['color'], linestyle=props['linestyle'],
                                   label=props['label'], linewidth=2)
                    print(f"  ✓ Loaded {method}: {len(t)} time steps")
            except Exception as e:
                print(f"  ⚠️  Failed to load {method}: {e}")
        else:
            print(f"  ⚠️  No energy file for {method}")

    axes[0, 0].set_xlabel('Time')
    axes[0, 0].set_ylabel('Total Energy')
    axes[0, 0].set_title('Total Energy Evolution')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)

    axes[0, 1].set_xlabel('Time')
    axes[0, 1].set_ylabel('Kinetic Energy')
    axes[0, 1].set_title('Kinetic Energy (Stability Indicator)')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)

    axes[1, 0].set_xlabel('Time')
    axes[1, 0].set_ylabel('Energy Drift (%)')
    axes[1, 0].set_title('Energy Conservation Error')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)

    axes[1, 1].set_xlabel('Time')
    axes[1, 1].set_ylabel('Potential Energy')
    axes[1, 1].set_title('Gravitational Potential Energy')
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)

    plt.suptitle('Grad-h Comparison: 3 GSPH Variants\n(IMBH 10k particles, 3D kernel gravity)', fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(os.path.join(output_dir, 'gradh_energy_comparison.png'), dpi=150)
    plt.close()
    print(f"\n  📊 Saved: {output_dir}/gradh_energy_comparison.png")

if __name__ == "__main__":
    main()
