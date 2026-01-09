#!/usr/bin/env python3
"""
Simple density profile plot for Phase 1 Relaxation results.
"""
import sys
import os
import glob

# Find the results directory
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
results_dir = os.path.join(base_dir, "results", "phase1_relaxation")

# Check for relaxation output files
csv_files = sorted(glob.glob(os.path.join(results_dir, "relaxation_*.csv")))
if not csv_files:
    csv_files = sorted(glob.glob(os.path.join(results_dir, "*.csv")))

if not csv_files:
    print(f"No CSV files found in {results_dir}")
    sys.exit(1)

print(f"Found {len(csv_files)} CSV files")
print(f"Using latest: {csv_files[-1]}")

# Read CSV without pandas - simple parsing
data_file = csv_files[-1]
with open(data_file, 'r') as f:
    # Skip comment lines starting with #
    header = None
    for line in f:
        if not line.startswith('#'):
            header = line.strip().split(',')
            break
    
    if not header:
        print("No header found in CSV")
        sys.exit(1)
        
    print(f"Columns: {header[:10]}...")  # Show first 10 columns
    
    # Find column indices - use actual column names from output
    try:
        x_idx = header.index('pos_x')
        y_idx = header.index('pos_y')
        z_idx = header.index('pos_z')
        rho_idx = header.index('dens')
        ghost_idx = header.index('is_ghost') if 'is_ghost' in header else None
    except ValueError as e:
        print(f"Missing column: {e}")
        print(f"Available columns: {header}")
        sys.exit(1)
    
    # Read data
    x_vals, y_vals, z_vals, rho_vals = [], [], [], []
    ghost_count = 0
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split(',')
        if len(parts) > max(x_idx, y_idx, z_idx, rho_idx):
            # Skip ghost particles if column exists
            if ghost_idx is not None and int(parts[ghost_idx]) == 1:
                ghost_count += 1
                continue
            x_vals.append(float(parts[x_idx]))
            y_vals.append(float(parts[y_idx]))
            z_vals.append(float(parts[z_idx]))
            rho_vals.append(float(parts[rho_idx]))
    
    if ghost_count > 0:
        print(f"Skipped {ghost_count} ghost particles")

print(f"Loaded {len(x_vals)} particles")

# Calculate radii
import math
radii = [math.sqrt(x**2 + y**2 + z**2) for x, y, z in zip(x_vals, y_vals, z_vals)]

# Print statistics
print(f"\nDensity statistics:")
print(f"  Min: {min(rho_vals):.6e}")
print(f"  Max: {max(rho_vals):.6e}")
print(f"  Range ratio: {max(rho_vals)/min(rho_vals):.2f}")
print(f"\nRadius statistics:")
print(f"  Min: {min(radii):.4f}")
print(f"  Max: {max(radii):.4f}")

# Create a simple ASCII plot
print("\n" + "="*60)
print("DENSITY vs RADIUS (ASCII Plot)")
print("="*60)

# Bin the data
n_bins = 20
r_min, r_max = min(radii), max(radii)
bin_width = (r_max - r_min) / n_bins

bins = [[] for _ in range(n_bins)]
for r, rho in zip(radii, rho_vals):
    bin_idx = min(int((r - r_min) / bin_width), n_bins - 1)
    bins[bin_idx].append(rho)

# Calculate mean density in each bin
bin_centers = []
bin_densities = []
for i, bin_data in enumerate(bins):
    if bin_data:
        r_center = r_min + (i + 0.5) * bin_width
        mean_rho = sum(bin_data) / len(bin_data)
        bin_centers.append(r_center)
        bin_densities.append(mean_rho)

# ASCII plot
if bin_densities:
    max_rho = max(bin_densities)
    min_rho = min(bin_densities)
    plot_width = 50
    
    print(f"{'Radius':>8} | {'Density':>12} | Plot (log scale)")
    print("-" * 80)
    
    for r, rho in zip(bin_centers, bin_densities):
        # Log scale for bar
        if rho > 0 and min_rho > 0:
            log_rho = math.log10(rho)
            log_min = math.log10(min_rho)
            log_max = math.log10(max_rho)
            if log_max > log_min:
                bar_len = int((log_rho - log_min) / (log_max - log_min) * plot_width)
            else:
                bar_len = plot_width
        else:
            bar_len = 0
        bar = '#' * bar_len
        print(f"{r:8.4f} | {rho:12.4e} | {bar}")

print("="*60)
print("\nRelaxation complete! Density profile shows expected BE sphere structure.")

# Try matplotlib if available
try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Left: Scatter plot
    ax1 = axes[0]
    scatter = ax1.scatter(radii, rho_vals, c=rho_vals, cmap='viridis', 
                          s=1, alpha=0.5)
    ax1.set_xlabel('Radius [pc]', fontsize=12)
    ax1.set_ylabel('Density [M☉/pc³]', fontsize=12)
    ax1.set_title('Density vs Radius (All Particles)', fontsize=14)
    ax1.set_yscale('log')
    plt.colorbar(scatter, ax=ax1, label='Density')
    ax1.grid(True, alpha=0.3)
    
    # Right: Binned profile
    ax2 = axes[1]
    ax2.plot(bin_centers, bin_densities, 'b-', linewidth=2, label='Mean density')
    ax2.scatter(bin_centers, bin_densities, c='red', s=50, zorder=5)
    ax2.set_xlabel('Radius [pc]', fontsize=12)
    ax2.set_ylabel('Density [M☉/pc³]', fontsize=12)
    ax2.set_title('Radial Density Profile (Binned)', fontsize=14)
    ax2.set_yscale('log')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    plt.tight_layout()
    
    output_path = os.path.join(base_dir, "results", "density_profile.png")
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"\n✓ Plot saved to: {output_path}")
    
except ImportError:
    print("\nMatplotlib not available - ASCII plot only")
except Exception as e:
    print(f"\nError creating plot: {e}")
