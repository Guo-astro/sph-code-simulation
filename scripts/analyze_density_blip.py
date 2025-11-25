#!/usr/bin/env python3
"""
Analyze density blip in SR Sod shock tube results
Identifies the contact discontinuity and measures the density overshoot
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def load_snapshot(filename):
    """Load CSV snapshot and extract particle data"""
    print(f"Loading {filename}...")

    # Read CSV, skipping header lines that start with # and column names
    data = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            # Skip column header line (starts with "id")
            if line.strip().startswith('id'):
                continue
            parts = line.strip().split(',')
            if len(parts) >= 15:  # Ensure we have all columns
                try:
                    data.append([float(x) for x in parts])
                except ValueError:
                    continue  # Skip non-numeric lines

    data = np.array(data)

    # Extract columns based on header description:
    # id, pos_x, pos_y, pos_z, vel_x, vel_y, vel_z, acc_x, acc_y, acc_z,
    # mass, dens, pres, ene, sml, sound, alpha, balsara, gradh, phi, neighbor

    result = {
        'id': data[:, 0],
        'pos_x': data[:, 1],
        'vel_x': data[:, 4],
        'nu': data[:, 10],      # mass = baryon number
        'dens': data[:, 11],    # rest-frame density n
        'pres': data[:, 12],
        'ene': data[:, 13],
        'sml': data[:, 14],
        'sound': data[:, 15]
    }

    # Sort by position
    idx = np.argsort(result['pos_x'])
    for key in result:
        result[key] = result[key][idx]

    print(f"Loaded {len(result['id'])} particles")
    return result

def compute_gamma_and_N(data):
    """Compute Lorentz factor and lab-frame density N"""
    c = 1.0  # Speed of light in code units
    v = data['vel_x']

    # Lorentz factor
    gamma = 1.0 / np.sqrt(1.0 - v**2 / c**2)

    # Lab-frame density N = γ * n
    N = gamma * data['dens']

    return gamma, N

def find_contact_discontinuity(data, N):
    """Find contact discontinuity location (pressure constant, density jump)"""
    x = data['pos_x']
    P = data['pres']

    # Contact discontinuity is where pressure is roughly constant
    # but density changes significantly
    # Look for maximum density gradient
    dN_dx = np.abs(np.diff(N) / np.diff(x))

    # Ignore edges
    idx_max = np.argmax(dN_dx[5:-5]) + 5

    x_contact = (x[idx_max] + x[idx_max+1]) / 2

    print(f"\nContact discontinuity located at x ≈ {x_contact:.4f}")

    return x_contact, idx_max

def analyze_density_blip(data, N, x_contact, idx_contact):
    """Analyze the density blip around contact discontinuity"""
    x = data['pos_x']

    # Get densities near contact
    window = 15  # particles on each side
    idx_left = max(0, idx_contact - window)
    idx_right = min(len(x), idx_contact + window + 1)

    x_local = x[idx_left:idx_right]
    N_local = N[idx_left:idx_right]
    n_local = data['dens'][idx_left:idx_right]
    P_local = data['pres'][idx_left:idx_right]
    v_local = data['vel_x'][idx_left:idx_right]

    # Find left and right plateau values (away from discontinuity)
    N_left_plateau = np.mean(N[idx_contact-window:idx_contact-5])
    N_right_plateau = np.mean(N[idx_contact+5:idx_contact+window])

    # Find maximum overshoot/undershoot
    N_max_local = np.max(N_local)
    N_min_local = np.min(N_local)

    # Overshoot on left side (above left plateau)
    overshoot_left = (N_max_local - N_left_plateau) / N_left_plateau * 100

    # Undershoot on right side (below right plateau)
    undershoot_right = (N_right_plateau - N_min_local) / N_right_plateau * 100

    print(f"\n=== DENSITY BLIP ANALYSIS ===")
    print(f"Left plateau N:  {N_left_plateau:.6f}")
    print(f"Right plateau N: {N_right_plateau:.6f}")
    print(f"Max N in region: {N_max_local:.6f}")
    print(f"Min N in region: {N_min_local:.6f}")
    print(f"Overshoot (left):  {overshoot_left:.2f}%")
    print(f"Undershoot (right): {undershoot_right:.2f}%")

    # Check pressure constancy
    P_left = np.mean(P_local[P_local < np.median(P_local)])
    P_right = np.mean(P_local[P_local > np.median(P_local)])
    P_variation = abs(P_left - P_right) / np.mean([P_left, P_right]) * 100
    print(f"\nPressure variation across contact: {P_variation:.2f}%")

    return {
        'x': x_local,
        'N': N_local,
        'n': n_local,
        'P': P_local,
        'v': v_local,
        'overshoot_pct': overshoot_left,
        'undershoot_pct': undershoot_right
    }

def plot_results(data, N, x_contact, blip_data, output_file):
    """Create diagnostic plots"""
    fig, axes = plt.subplots(3, 1, figsize=(12, 10))

    x = data['pos_x']

    # Plot 1: Full domain - lab-frame density N
    ax = axes[0]
    ax.plot(x, N, 'ko-', markersize=4, label='Lab-frame N = γn')
    ax.axvline(x_contact, color='r', linestyle='--', alpha=0.5, label='Contact disc.')
    ax.set_xlabel('Position x')
    ax.set_ylabel('Lab-frame Density N')
    ax.set_title('SR Sod Test: Lab-Frame Density N = γn')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 2: Zoom on contact discontinuity
    ax = axes[1]
    ax.plot(blip_data['x'], blip_data['N'], 'ko-', markersize=6, linewidth=2)
    ax.axvline(x_contact, color='r', linestyle='--', alpha=0.5, label='Contact')
    ax.set_xlabel('Position x')
    ax.set_ylabel('Lab-frame Density N')
    ax.set_title(f'Density Blip at Contact Discontinuity '
                 f'(Overshoot: {blip_data["overshoot_pct"]:.1f}%, '
                 f'Undershoot: {blip_data["undershoot_pct"]:.1f}%)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 3: Pressure across contact
    ax = axes[2]
    ax.plot(blip_data['x'], blip_data['P'], 'bo-', markersize=6, linewidth=2)
    ax.axvline(x_contact, color='r', linestyle='--', alpha=0.5, label='Contact')
    ax.set_xlabel('Position x')
    ax.set_ylabel('Pressure P')
    ax.set_title('Pressure Across Contact Discontinuity (should be constant)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\nPlot saved to: {output_file}")

    return fig

def main():
    if len(sys.argv) < 2:
        print("Usage: python analyze_density_blip.py <snapshot_file.csv>")
        print("\nExample:")
        print("  python analyze_density_blip.py sample/sr_sod/results/sharp/snapshot_0131.csv")
        sys.exit(1)

    snapshot_file = sys.argv[1]

    if not os.path.exists(snapshot_file):
        print(f"Error: File not found: {snapshot_file}")
        sys.exit(1)

    # Load data
    data = load_snapshot(snapshot_file)

    # Compute derived quantities
    gamma, N = compute_gamma_and_N(data)

    # Find contact discontinuity
    x_contact, idx_contact = find_contact_discontinuity(data, N)

    # Analyze density blip
    blip_data = analyze_density_blip(data, N, x_contact, idx_contact)

    # Create diagnostic plots
    output_dir = os.path.dirname(snapshot_file)
    output_file = os.path.join(output_dir, "density_blip_analysis.png")
    plot_results(data, N, x_contact, blip_data, output_file)

    print("\n" + "="*60)
    if blip_data['overshoot_pct'] > 5.0 or blip_data['undershoot_pct'] > 5.0:
        print("⚠️  LARGE DENSITY BLIP DETECTED!")
        print("\nLikely causes:")
        print("  1. C_cd parameter too small (should be 1.0, not 0.2)")
        print("  2. V²_ij interpolation issues")
        print("  3. Insufficient smoothing (C_smooth too small)")
    else:
        print("✓ Density blip is within acceptable range (<5%)")
    print("="*60)

if __name__ == "__main__":
    main()
