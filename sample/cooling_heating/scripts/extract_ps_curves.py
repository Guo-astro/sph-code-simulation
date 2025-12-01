#!/usr/bin/env python3
"""
Extract actual curve data from PostScript files.
Uses moveto (M) and relative lineto (V) commands.
"""

import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def extract_curves_from_ps(filename):
    """Extract all curves from PostScript file."""
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    curves = []
    current_curve = []
    current_x, current_y = 0, 0
    in_curve = False
    
    for line in lines:
        line = line.strip()
        
        # Start of new curve with absolute position
        if re.match(r'^\d+ \d+ M$', line):
            if current_curve:
                curves.append(np.array(current_curve))
                current_curve = []
            x, y, _ = line.split()
            current_x, current_y = int(x), int(y)
            current_curve.append([current_x, current_y])
            in_curve = True
        
        # Relative movement (V for vertical-ish, H for horizontal)
        elif in_curve and re.match(r'^-?\d+ -?\d+ V$', line):
            dx, dy, _ = line.split()
            current_x += int(dx)
            current_y += int(dy)
            current_curve.append([current_x, current_y])
        
        # End curve at line type change
        elif line.startswith('LT') and current_curve:
            curves.append(np.array(current_curve))
            current_curve = []
            in_curve = False
    
    if current_curve:
        curves.append(np.array(current_curve))
    
    return curves

def find_axis_range(filename):
    """Extract axis range from PostScript comments/labels."""
    
    with open(filename, 'r') as f:
        content = f.read()
    
    # Look for axis labels in the PostScript
    # Gnuplot typically embeds these as text strings
    x_labels = re.findall(r'\((-?\d+(?:\.\d+)?(?:e[+-]?\d+)?)\) Rshow', content)
    y_labels = re.findall(r'\((-?\d+(?:\.\d+)?(?:e[+-]?\d+)?)\) Cshow', content)
    
    print(f"Found X-axis labels: {x_labels}")
    print(f"Found Y-axis labels: {y_labels}")
    
    return x_labels, y_labels

def convert_ps_to_data(curves, x_range, y_range, x_canvas_range, y_canvas_range):
    """
    Convert PostScript canvas coordinates to actual data coordinates.
    
    Parameters:
    -----------
    curves : list of arrays
        PostScript canvas coordinates
    x_range : tuple
        (log10(x_min), log10(x_max)) in data coordinates
    y_range : tuple
        (log10(y_min), log10(y_max)) in data coordinates
    x_canvas_range : tuple
        (x_min_canvas, x_max_canvas)
    y_canvas_range : tuple
        (y_min_canvas, y_max_canvas)
    """
    
    converted_curves = []
    
    for curve in curves:
        if len(curve) < 5:  # Skip short paths (likely labels/ticks)
            continue
        
        # Convert canvas to data coordinates (log scale)
        x_canvas = curve[:, 0]
        y_canvas = curve[:, 1]
        
        # Linear interpolation in log space
        log_x = x_range[0] + (x_canvas - x_canvas_range[0]) / (x_canvas_range[1] - x_canvas_range[0]) * (x_range[1] - x_range[0])
        log_y = y_range[0] + (y_canvas - y_canvas_range[0]) / (y_canvas_range[1] - y_canvas_range[0]) * (y_range[1] - y_range[0])
        
        # Convert from log to linear
        x_data = 10**log_x
        y_data = 10**log_y
        
        converted_curves.append(np.column_stack([x_data, y_data]))
    
    return converted_curves

def main():
    ps_file = '/Users/guo/Downloads/sphcode/docs/papers/cooling-heating/f1a.ps'
    
    print("=" * 70)
    print("Extracting Curve Data from f1a.ps")
    print("=" * 70)
    
    curves = extract_curves_from_ps(ps_file)
    
    print(f"\nFound {len(curves)} curves")
    for i, curve in enumerate(curves):
        print(f"Curve {i}: {len(curve)} points")
        if len(curve) > 5:
            print(f"  X range (canvas): {curve[:, 0].min():.0f} - {curve[:, 0].max():.0f}")
            print(f"  Y range (canvas): {curve[:, 1].min():.0f} - {curve[:, 1].max():.0f}")
    
    # Find axis labels
    print("\nExtracting axis information...")
    x_labels, y_labels = find_axis_range(ps_file)
    
    # Plot raw canvas coordinates
    fig, ax = plt.subplots(figsize=(12, 8))
    
    for i, curve in enumerate(curves):
        if len(curve) > 5:
            ax.plot(curve[:, 0], curve[:, 1], 'o-', linewidth=2, markersize=3,
                   label=f'Curve {i} ({len(curve)} pts)')
    
    ax.set_xlabel('Canvas X')
    ax.set_ylabel('Canvas Y')
    ax.set_title('Raw PostScript Coordinates from f1a.ps')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('../results/ps_raw_curves.png', dpi=200, bbox_inches='tight')
    print("\n✓ Raw curves plot saved: ../results/ps_raw_curves.png")
    
    # Try to convert to data coordinates
    # Figure 1(a) should be: x: log(n) from -1 to 6, y: log(T) from 1 to 5
    print("\nAttempting coordinate conversion...")
    
    # Estimate canvas range from the curves
    all_points = np.vstack([c for c in curves if len(c) > 5])
    x_canvas_range = (all_points[:, 0].min(), all_points[:, 0].max())
    y_canvas_range = (all_points[:, 1].min(), all_points[:, 1].max())
    
    print(f"Canvas X range: {x_canvas_range}")
    print(f"Canvas Y range: {y_canvas_range}")
    
    # Known data ranges for Figure 1(a)
    x_data_range = (-1, 6)  # log10(n) from 0.1 to 1e6
    y_data_range = (1, 5)   # log10(T) from 10 to 100000
    
    data_curves = convert_ps_to_data(curves, x_data_range, y_data_range,
                                     x_canvas_range, y_canvas_range)
    
    # Plot converted data
    fig2, ax2 = plt.subplots(figsize=(12, 8))
    
    for i, curve in enumerate(data_curves):
        ax2.loglog(curve[:, 0], curve[:, 1], 'o-', linewidth=2.5, markersize=4,
                  label=f'Curve {i} ({len(curve)} pts)')
    
    ax2.set_xlabel('n [cm$^{-3}$]', fontsize=13)
    ax2.set_ylabel('T [K] or P/k$_B$ [K cm$^{-3}$]', fontsize=13)
    ax2.set_title('Converted Data from f1a.ps', fontsize=14, fontweight='bold')
    ax2.set_xlim(0.1, 1e6)
    ax2.set_ylim(10, 1e5)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('../results/ps_converted_curves.png', dpi=200, bbox_inches='tight')
    print("✓ Converted curves plot saved: ../results/ps_converted_curves.png")
    
    # Save the main temperature curve as digitized data
    if len(data_curves) > 0:
        # The longest curve is likely the temperature curve
        main_curve_idx = max(range(len(data_curves)), key=lambda i: len(data_curves[i]))
        main_curve = data_curves[main_curve_idx]
        
        print(f"\nSaving main curve (index {main_curve_idx}) as digitized data...")
        print(f"  {len(main_curve)} data points")
        
        # Save to file
        output_file = '../results/f1a_digitized.txt'
        np.savetxt(output_file, main_curve, 
                  header='n[cm^-3]  T[K] or P/k_B[K/cm^3] - Digitized from f1a.ps',
                  fmt='%.6e')
        print(f"✓ Digitized data saved: {output_file}")
        
        # Print first/last few points
        print("\nFirst 10 points (n, value):")
        for i in range(min(10, len(main_curve))):
            print(f"  {main_curve[i, 0]:.3e}  {main_curve[i, 1]:.3e}")
        
        print(f"\nLast 10 points:")
        for i in range(max(0, len(main_curve)-10), len(main_curve)):
            print(f"  {main_curve[i, 0]:.3e}  {main_curve[i, 1]:.3e}")
    
    print("\n" + "=" * 70)
    print("✓ PostScript data extraction complete!")
    print("=" * 70)

if __name__ == '__main__':
    main()
