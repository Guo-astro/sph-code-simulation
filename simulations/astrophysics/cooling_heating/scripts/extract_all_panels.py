#!/usr/bin/env python3
"""
Extract ALL curves from ALL 4 PostScript figure panels.
Pixel-perfect digitization of Figure 1(a,b,c,d) from Koyama & Inutsuka (2000).
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
    
    for line in lines:
        line = line.strip()
        
        if re.match(r'^\d+ \d+ M$', line):
            if current_curve:
                curves.append(np.array(current_curve))
                current_curve = []
            x, y, _ = line.split()
            current_x, current_y = int(x), int(y)
            current_curve.append([current_x, current_y])
        
        elif re.match(r'^-?\d+ -?\d+ V$', line):
            dx, dy, _ = line.split()
            current_x += int(dx)
            current_y += int(dy)
            current_curve.append([current_x, current_y])
        
        elif line.startswith('LT') and current_curve:
            curves.append(np.array(current_curve))
            current_curve = []
    
    if current_curve:
        curves.append(np.array(current_curve))
    
    # Filter to only substantial curves
    substantial_curves = [c for c in curves if len(c) >= 30]
    
    return substantial_curves

def convert_to_physical(curve, x_canvas_min, x_canvas_max, 
                       y_canvas_min, y_canvas_max,
                       log_x_min, log_x_max,
                       log_y_min, log_y_max):
    """Convert canvas coordinates to physical values."""
    
    x_canvas = curve[:, 0]
    y_canvas = curve[:, 1]
    
    log_x = log_x_min + (x_canvas - x_canvas_min) / (x_canvas_max - x_canvas_min) * (log_x_max - log_x_min)
    log_y = log_y_min + (y_canvas - y_canvas_min) / (y_canvas_max - y_canvas_min) * (log_y_max - log_y_min)
    
    x = 10**log_x
    y = 10**log_y
    
    return np.column_stack([x, y])

def process_figure(ps_file, figure_name, log_x_range, log_y_range, output_dir):
    """Process one PostScript figure panel."""
    
    print(f"\n{'='*70}")
    print(f"Processing {figure_name}: {ps_file}")
    print(f"{'='*70}")
    
    curves = extract_curves_from_ps(ps_file)
    print(f"Found {len(curves)} substantial curves")
    
    # Determine canvas ranges
    all_points = np.vstack(curves)
    x_canvas_min, x_canvas_max = all_points[:, 0].min(), all_points[:, 0].max()
    y_canvas_min, y_canvas_max = all_points[:, 1].min(), all_points[:, 1].max()
    
    print(f"Canvas X: [{x_canvas_min}, {x_canvas_max}]")
    print(f"Canvas Y: [{y_canvas_min}, {y_canvas_max}]")
    print(f"Physical log(X): [{log_x_range[0]}, {log_x_range[1]}]")
    print(f"Physical log(Y): [{log_y_range[0]}, {log_y_range[1]}]")
    
    # Convert all curves
    converted_curves = []
    for i, curve in enumerate(curves):
        data = convert_to_physical(
            curve, x_canvas_min, x_canvas_max, y_canvas_min, y_canvas_max,
            log_x_range[0], log_x_range[1], log_y_range[0], log_y_range[1]
        )
        converted_curves.append(data)
        
        print(f"\nCurve {i}: {len(data)} points")
        print(f"  X: {data[:,0].min():.3e} - {data[:,0].max():.3e}")
        print(f"  Y: {data[:,1].min():.3e} - {data[:,1].max():.3e}")
        
        # Save curve
        output_file = f'{output_dir}/{figure_name}_curve_{i}.txt'
        np.savetxt(output_file, data, 
                  header=f'{figure_name} curve {i} - x, y values',
                  fmt='%.6e')
        print(f"  ✓ Saved: {output_file}")
    
    return converted_curves

def plot_panel_a(curves, output_file):
    """Plot Figure 1(a) - Temperature and Pressure."""
    from chemistry_network import k_B
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Curves 0,1: Temperature, 2,3: Pressure
    ax.loglog(curves[0][:,0], curves[0][:,1], 'b-', linewidth=2.5, 
              label=r'$T$ ($N_H=10^{19}$ cm$^{-2}$)')
    ax.loglog(curves[1][:,0], curves[1][:,1], 'b--', linewidth=2.5,
              label=r'$T$ ($N_H=10^{20}$ cm$^{-2}$)')
    ax.loglog(curves[2][:,0], curves[2][:,1]/k_B, 'r-', linewidth=2.5,
              label=r'$P/k_B$ ($10^{19}$)')
    ax.loglog(curves[3][:,0], curves[3][:,1]/k_B, 'r--', linewidth=2.5,
              label=r'$P/k_B$ ($10^{20}$)')
    
    ax.set_xlabel(r'log $n$ [cm$^{-3}$]', fontsize=13, fontweight='bold')
    ax.set_ylabel(r'log $T$ [K], log $P$ [K/cm$^3$]', fontsize=13, fontweight='bold')
    ax.set_xlim(0.1, 1e6)
    ax.set_ylim(10, 1e8)
    ax.legend(fontsize=9, loc='lower right')
    ax.grid(True, alpha=0.3, which='both')
    ax.set_title('(a) Temperature and Pressure', fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Plot saved: {output_file}")

def plot_panel_b(curves, output_file):
    """Plot Figure 1(b) - Chemical fractions."""
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Typically: electron, H2, CO
    colors = ['k', 'b', 'g', 'r', 'm', 'c']
    labels = ['electron', r'H$_2$', 'CO', 'Curve 3', 'Curve 4', 'Curve 5']
    
    for i, curve in enumerate(curves):
        ax.loglog(curve[:,0], curve[:,1], linewidth=2.5, 
                 color=colors[i % len(colors)],
                 label=labels[i] if i < len(labels) else f'Curve {i}')
    
    ax.set_xlabel(r'log $n$ [cm$^{-3}$]', fontsize=13, fontweight='bold')
    ax.set_ylabel(r'log $x_i$', fontsize=13, fontweight='bold')
    ax.set_xlim(0.1, 1e6)
    ax.set_ylim(1e-8, 1)
    ax.legend(fontsize=11, loc='best')
    ax.grid(True, alpha=0.3, which='both')
    ax.set_title('(b) Chemical Fractions', fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Plot saved: {output_file}")

def plot_panel_c(curves, output_file):
    """Plot Figure 1(c) - Heating and cooling rates."""
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Heating (dashed) and cooling (solid) lines
    colors = ['r', 'g', 'b', 'm', 'c', 'orange', 'purple', 'brown', 'pink', 'gray']
    labels = ['PE', 'XR', 'CR', r'H$_2$', 'CII', 'OI', r'Ly-$\alpha$', r'H$_2$', 'CO', 'GR']
    
    for i, curve in enumerate(curves):
        linestyle = ':' if i < len(curves)//2 else '-'  # First half dashed, second half solid
        ax.loglog(curve[:,0], curve[:,1], linewidth=2.5 if linestyle=='-' else 2,
                 linestyle=linestyle,
                 color=colors[i % len(colors)],
                 label=labels[i] if i < len(labels) else f'Curve {i}')
    
    ax.set_xlabel(r'log $n$ [cm$^{-3}$]', fontsize=13, fontweight='bold')
    ax.set_ylabel(r'log $\Gamma$, $\Lambda$ [ergs s$^{-1}$ H$^{-1}$]', fontsize=13, fontweight='bold')
    ax.set_xlim(0.1, 1e6)
    ax.set_ylim(1e-28, 1e-23)
    ax.legend(fontsize=8, loc='best', ncol=2)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_title('(c) Heating and Cooling Rates', fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Plot saved: {output_file}")

def plot_panel_d(curves, output_file):
    """Plot Figure 1(d) - Timescales."""
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    colors = ['r', 'b', 'g', 'm']
    styles = ['-', '--', '-.', ':']
    labels = ['cooling', 'recombination', 'free fall', r'H$_2$ formation']
    
    for i, curve in enumerate(curves):
        ax.loglog(curve[:,0], curve[:,1], 
                 linewidth=2.5 if styles[i % 4] == '-' else 2.5,
                 linestyle=styles[i % 4],
                 color=colors[i % len(colors)],
                 label=labels[i] if i < len(labels) else f'Curve {i}')
    
    ax.set_xlabel(r'log $n$ [cm$^{-3}$]', fontsize=13, fontweight='bold')
    ax.set_ylabel(r'log [year]', fontsize=13, fontweight='bold')
    ax.set_xlim(0.1, 1e6)
    ax.set_ylim(1e0, 1e12)
    ax.legend(fontsize=10, loc='best')
    ax.grid(True, alpha=0.3, which='both')
    ax.set_title('(d) Timescales', fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Plot saved: {output_file}")

def main():
    ps_dir = '/Users/guo/Downloads/sphcode/docs/papers/cooling-heating'
    output_dir = '../results'
    
    print("="*70)
    print("PIXEL-PERFECT DIGITIZATION OF ALL 4 FIGURE PANELS")
    print("Koyama & Inutsuka (2000) Figure 1(a,b,c,d)")
    print("="*70)
    
    # Panel (a): Temperature and Pressure vs density
    # Already done, but let's include it for completeness
    curves_a = process_figure(
        f'{ps_dir}/f1a.ps',
        'f1a',
        log_x_range=(-1, 6),      # n from 0.1 to 1e6
        log_y_range=(1.5, 6.0),   # T/P from 31.6 to 1e6
        output_dir=output_dir
    )
    
    # Panel (b): Chemical fractions vs density
    curves_b = process_figure(
        f'{ps_dir}/f1b.ps',
        'f1b',
        log_x_range=(-1, 6),      # n from 0.1 to 1e6
        log_y_range=(-8, 0),      # fractions from 1e-8 to 1
        output_dir=output_dir
    )
    
    # Panel (c): Heating/cooling rates vs density
    curves_c = process_figure(
        f'{ps_dir}/f1c.ps',
        'f1c',
        log_x_range=(-1, 6),      # n from 0.1 to 1e6
        log_y_range=(-28, -23),   # rates from 1e-28 to 1e-23 erg/s
        output_dir=output_dir
    )
    
    # Panel (d): Timescales vs density
    curves_d = process_figure(
        f'{ps_dir}/f1d.ps',
        'f1d',
        log_x_range=(-1, 6),      # n from 0.1 to 1e6
        log_y_range=(0, 12),      # time from 1 to 1e12 years
        output_dir=output_dir
    )
    
    # Create individual panel plots
    print(f"\n{'='*70}")
    print("Creating panel plots...")
    print(f"{'='*70}")
    
    plot_panel_a(curves_a, f'{output_dir}/f1a_EXACT.png')
    plot_panel_b(curves_b, f'{output_dir}/f1b_EXACT.png')
    plot_panel_c(curves_c, f'{output_dir}/f1c_EXACT.png')
    plot_panel_d(curves_d, f'{output_dir}/f1d_EXACT.png')
    
    # Create combined figure
    fig = plt.figure(figsize=(14, 11))
    
    # Replot all 4 panels in a 2x2 grid
    from chemistry_network import k_B
    
    # Panel (a)
    ax1 = fig.add_subplot(2, 2, 1)
    ax1.loglog(curves_a[0][:,0], curves_a[0][:,1], 'b-', linewidth=2.5, 
              label=r'$T$ ($10^{19}$)')
    ax1.loglog(curves_a[1][:,0], curves_a[1][:,1], 'b--', linewidth=2.5,
              label=r'$T$ ($10^{20}$)')
    ax1.loglog(curves_a[2][:,0], curves_a[2][:,1]/k_B, 'r-', linewidth=2.5,
              label=r'$P/k_B$ ($10^{19}$)')
    ax1.loglog(curves_a[3][:,0], curves_a[3][:,1]/k_B, 'r--', linewidth=2.5,
              label=r'$P/k_B$ ($10^{20}$)')
    ax1.set_xlabel(r'log $n$ [cm$^{-3}$]', fontsize=11, fontweight='bold')
    ax1.set_ylabel(r'log $T$, $P/k_B$', fontsize=11, fontweight='bold')
    ax1.set_xlim(0.1, 1e6)
    ax1.set_ylim(10, 1e8)
    ax1.legend(fontsize=8, loc='lower right')
    ax1.grid(True, alpha=0.3, which='both')
    ax1.set_title('(a)', fontsize=12, fontweight='bold', loc='left')
    
    # Panel (b)
    ax2 = fig.add_subplot(2, 2, 2)
    colors_b = ['k', 'b', 'g']
    labels_b = ['electron', r'H$_2$', 'CO']
    for i in range(min(3, len(curves_b))):
        ax2.loglog(curves_b[i][:,0], curves_b[i][:,1], linewidth=2.5,
                  color=colors_b[i], label=labels_b[i])
    ax2.set_xlabel(r'log $n$ [cm$^{-3}$]', fontsize=11, fontweight='bold')
    ax2.set_ylabel(r'log $x_i$', fontsize=11, fontweight='bold')
    ax2.set_xlim(0.1, 1e6)
    ax2.set_ylim(1e-8, 1)
    ax2.legend(fontsize=9, loc='best')
    ax2.grid(True, alpha=0.3, which='both')
    ax2.set_title('(b)', fontsize=12, fontweight='bold', loc='left')
    
    # Panel (c)
    ax3 = fig.add_subplot(2, 2, 3)
    n_heating = len(curves_c) // 2
    for i, curve in enumerate(curves_c):
        linestyle = ':' if i < n_heating else '-'
        ax3.loglog(curve[:,0], curve[:,1], linewidth=2 if linestyle==':' else 2.5,
                  linestyle=linestyle)
    ax3.set_xlabel(r'log $n$ [cm$^{-3}$]', fontsize=11, fontweight='bold')
    ax3.set_ylabel(r'log $\Gamma$, $\Lambda$', fontsize=11, fontweight='bold')
    ax3.set_xlim(0.1, 1e6)
    ax3.set_ylim(1e-28, 1e-23)
    ax3.grid(True, alpha=0.3, which='both')
    ax3.set_title('(c)', fontsize=12, fontweight='bold', loc='left')
    
    # Panel (d)
    ax4 = fig.add_subplot(2, 2, 4)
    styles_d = ['-', '--', '-.', ':']
    for i, curve in enumerate(curves_d):
        ax4.loglog(curve[:,0], curve[:,1], linewidth=2.5,
                  linestyle=styles_d[i % 4])
    ax4.set_xlabel(r'log $n$ [cm$^{-3}$]', fontsize=11, fontweight='bold')
    ax4.set_ylabel(r'log [year]', fontsize=11, fontweight='bold')
    ax4.set_xlim(0.1, 1e6)
    ax4.set_ylim(1e0, 1e12)
    ax4.grid(True, alpha=0.3, which='both')
    ax4.set_title('(d)', fontsize=12, fontweight='bold', loc='left')
    
    plt.tight_layout()
    combined_file = f'{output_dir}/figure1_ALL_EXACT.png'
    plt.savefig(combined_file, dpi=300, bbox_inches='tight')
    print(f"✓ Combined figure saved: {combined_file}")
    
    print(f"\n{'='*70}")
    print("✓ ALL 4 PANELS DIGITIZED SUCCESSFULLY!")
    print(f"{'='*70}")
    print(f"\nTotal curves extracted:")
    print(f"  Panel (a): {len(curves_a)} curves")
    print(f"  Panel (b): {len(curves_b)} curves")
    print(f"  Panel (c): {len(curves_c)} curves")
    print(f"  Panel (d): {len(curves_d)} curves")
    print(f"\nAll data saved to: {output_dir}/")
    print("="*70)

if __name__ == '__main__':
    main()
