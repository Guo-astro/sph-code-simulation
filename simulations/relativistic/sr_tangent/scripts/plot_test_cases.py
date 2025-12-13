#!/usr/bin/env python3
"""
Generate visualization plots for SR tangent velocity test cases.

This script creates various plots comparing the tabulated test case
solutions from Pons et al. (2000) and Rezzolla et al. (2003).

Usage:
    python plot_test_cases.py                    # Generate all plots
    python plot_test_cases.py --plot pressure   # Generate specific plot
    python plot_test_cases.py --list            # List available plots
"""

import sys
import argparse
from pathlib import Path

# Add script directory to path
SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))

import matplotlib.pyplot as plt
import numpy as np
from sr_tangent_test_cases import (
    PONS2000_TESTS, REZZOLLA2003_TESTS, 
    get_all_tests, WavePattern
)


def setup_plot_style():
    """Set up consistent plot style."""
    plt.style.use('default')
    plt.rcParams.update({
        'font.size': 11,
        'axes.labelsize': 12,
        'axes.titlesize': 13,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
        'figure.dpi': 100,
    })


def get_output_dir():
    """Get output directory for plots."""
    output_dir = SCRIPT_DIR.parent / "results" / "plots"
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def plot_pressure_comparison(output_dir: Path) -> None:
    """Plot P* (contact pressure) for all test cases."""
    tests = get_all_tests()
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    P_star = [t.expected.P_star for t in tests]
    x = range(len(tests))
    
    # Color by source
    colors = ['blue'] * len(PONS2000_TESTS) + ['green'] * len(REZZOLLA2003_TESTS)
    
    ax.semilogy(x, P_star, 'o-', markersize=8, color='gray', alpha=0.3)
    for i, (xi, Pi, c) in enumerate(zip(x, P_star, colors)):
        ax.semilogy(xi, Pi, 'o', markersize=8, color=c)
    
    ax.set_xlabel('Test Case Index')
    ax.set_ylabel(r'$P^*$ (Contact Pressure)')
    ax.set_title('Contact Pressure Across All Test Cases')
    ax.set_xticks(x)
    labels = [t.name.replace('pons2000_', 'P_').replace('rezzolla2003_', 'R_') 
              for t in tests]
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
    
    # Add legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', 
               markersize=10, label='Pons et al. (2000)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='green', 
               markersize=10, label='Rezzolla et al. (2003)')
    ]
    ax.legend(handles=legend_elements, loc='upper right')
    
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    
    output_path = output_dir / "pressure_comparison.png"
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"✓ Saved: {output_path}")


def plot_velocity_comparison(output_dir: Path) -> None:
    """Plot v^x* (contact velocity) for all test cases."""
    tests = get_all_tests()
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    vx_star = [t.expected.vx_star for t in tests]
    x = range(len(tests))
    
    colors = ['blue'] * len(PONS2000_TESTS) + ['green'] * len(REZZOLLA2003_TESTS)
    
    ax.plot(x, vx_star, 'o-', markersize=8, color='gray', alpha=0.3)
    for i, (xi, vi, c) in enumerate(zip(x, vx_star, colors)):
        ax.plot(xi, vi, 'o', markersize=8, color=c)
    
    ax.set_xlabel('Test Case Index')
    ax.set_ylabel(r'$v^x_*$ (Contact Velocity)')
    ax.set_title('Contact Velocity Across All Test Cases')
    ax.set_xticks(x)
    labels = [t.name.replace('pons2000_', 'P_').replace('rezzolla2003_', 'R_') 
              for t in tests]
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
    ax.axhline(y=0, color='k', linestyle='--', alpha=0.3)
    
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    
    output_path = output_dir / "velocity_comparison.png"
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"✓ Saved: {output_path}")


def plot_density_comparison(output_dir: Path) -> None:
    """Plot ρ'_L and ρ'_R (post-wave densities) for all test cases."""
    tests = get_all_tests()
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    rho_L = [t.expected.rho_L_prime for t in tests]
    rho_R = [t.expected.rho_R_prime for t in tests]
    x = range(len(tests))
    
    ax.semilogy(x, rho_L, '^-', markersize=7, color='blue', label=r"$\rho'_L$ (left of contact)")
    ax.semilogy(x, rho_R, 'v-', markersize=7, color='red', label=r"$\rho'_R$ (right of contact)")
    
    ax.set_xlabel('Test Case Index')
    ax.set_ylabel('Post-Wave Density')
    ax.set_title('Post-Wave Densities Across All Test Cases')
    ax.legend()
    ax.set_xticks(x)
    labels = [t.name.replace('pons2000_', 'P_').replace('rezzolla2003_', 'R_') 
              for t in tests]
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
    
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    
    output_path = output_dir / "density_comparison.png"
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"✓ Saved: {output_path}")


def plot_wave_pattern(output_dir: Path) -> None:
    """Plot wave pattern distribution."""
    tests = get_all_tests()
    
    # Count patterns
    patterns = {}
    for t in tests:
        p = t.wave_pattern.name
        patterns[p] = patterns.get(p, 0) + 1
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    colors = {'RS': 'blue', 'SR': 'green', 'SS': 'red', 'RR': 'orange'}
    pattern_names = list(patterns.keys())
    counts = list(patterns.values())
    bar_colors = [colors.get(k, 'gray') for k in pattern_names]
    
    bars = ax.bar(pattern_names, counts, color=bar_colors, edgecolor='black')
    
    ax.set_xlabel('Wave Pattern')
    ax.set_ylabel('Count')
    ax.set_title('Distribution of Wave Patterns in Test Suite')
    
    # Add count labels on bars
    for bar, count in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.2, 
                str(count), ha='center', fontweight='bold', fontsize=12)
    
    # Add pattern descriptions
    pattern_desc = {
        'RS': 'Rarefaction-Shock',
        'SR': 'Shock-Rarefaction', 
        'SS': 'Shock-Shock',
        'RR': 'Rarefaction-Rarefaction'
    }
    legend_text = '\n'.join([f'{k}: {pattern_desc.get(k, k)}' for k in pattern_names])
    ax.text(0.98, 0.98, legend_text, transform=ax.transAxes, 
            verticalalignment='top', horizontalalignment='right',
            fontsize=10, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    ax.set_ylim(0, max(counts) + 2)
    plt.tight_layout()
    
    output_path = output_dir / "wave_pattern_distribution.png"
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"✓ Saved: {output_path}")


def plot_tangent_effect(output_dir: Path) -> None:
    """Plot effect of tangential velocity on P* (Pons2000 tests)."""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Extract data
    vt_L = [t.left.vt for t in PONS2000_TESTS]
    vt_R = [t.right.vt for t in PONS2000_TESTS]
    P_star = [t.expected.P_star for t in PONS2000_TESTS]
    
    labels = [f'({vL:.2f},{vR:.2f})' for vL, vR in zip(vt_L, vt_R)]
    x = range(len(PONS2000_TESTS))
    
    ax.semilogy(x, P_star, 'go-', markersize=12, linewidth=2)
    
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_xlabel(r'$(v^t_L, v^t_R)$')
    ax.set_ylabel(r'$P^*$ (Contact Pressure)')
    ax.set_title('Effect of Tangential Velocity on Contact Pressure\n(Pons et al. 2000)')
    
    # Add annotations for key observations
    ax.annotate('High $v^t_L$ reduces $P^*$', 
                xy=(3, P_star[3]), xytext=(4.5, P_star[3]*5),
                arrowprops=dict(arrowstyle='->', color='red'),
                fontsize=10, color='red')
    
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    
    output_path = output_dir / "tangent_velocity_effect.png"
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"✓ Saved: {output_path}")


def plot_rezzolla_vt_sweep(output_dir: Path) -> None:
    """Plot Rezzolla tests showing transition between wave patterns."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # First set: v^x_L=0.5, v^x_R=0, varying v^t_R
    tests_1 = [t for t in REZZOLLA2003_TESTS if 'vx05_' in t.name]
    vt_R = [t.right.vt for t in tests_1]
    P_star_1 = [t.expected.P_star for t in tests_1]
    vx_star_1 = [t.expected.vx_star for t in tests_1]
    patterns_1 = [t.wave_pattern.name for t in tests_1]
    
    ax1 = axes[0]
    colors_1 = ['green' if p == 'SR' else 'red' for p in patterns_1]
    ax1.scatter(vt_R, P_star_1, c=colors_1, s=100, zorder=3)
    ax1.plot(vt_R, P_star_1, 'k--', alpha=0.3, zorder=2)
    ax1.set_xlabel(r'$v^t_R$')
    ax1.set_ylabel(r'$P^*$')
    ax1.set_title(r'$v^x_L=0.5, v^x_R=0$: Varying Right Tangent Velocity')
    ax1.grid(True, alpha=0.3)
    
    # Add pattern transition annotation
    ax1.axvline(x=0.85, color='gray', linestyle=':', alpha=0.5)
    ax1.text(0.6, max(P_star_1)*0.9, 'SR', fontsize=12, color='green', fontweight='bold')
    ax1.text(0.92, max(P_star_1)*0.9, '2S', fontsize=12, color='red', fontweight='bold')
    
    # Second set: v^x_L=0, v^x_R=0.5, varying v^t_L
    tests_2 = [t for t in REZZOLLA2003_TESTS if 'vx00_vx05_' in t.name]
    vt_L = [t.left.vt for t in tests_2]
    P_star_2 = [t.expected.P_star for t in tests_2]
    patterns_2 = [t.wave_pattern.name for t in tests_2]
    
    ax2 = axes[1]
    colors_2 = ['green' if p == 'SR' else 'orange' for p in patterns_2]
    ax2.scatter(vt_L, P_star_2, c=colors_2, s=100, zorder=3)
    ax2.plot(vt_L, P_star_2, 'k--', alpha=0.3, zorder=2)
    ax2.set_xlabel(r'$v^t_L$')
    ax2.set_ylabel(r'$P^*$')
    ax2.set_title(r'$v^x_L=0, v^x_R=0.5$: Varying Left Tangent Velocity')
    ax2.grid(True, alpha=0.3)
    
    # Add pattern transition annotation
    ax2.axvline(x=0.6, color='gray', linestyle=':', alpha=0.5)
    ax2.text(0.2, max(P_star_2)*0.9, 'SR', fontsize=12, color='green', fontweight='bold')
    ax2.text(0.8, max(P_star_2)*0.9, '2R', fontsize=12, color='orange', fontweight='bold')
    
    plt.tight_layout()
    
    output_path = output_dir / "rezzolla_vt_sweep.png"
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"✓ Saved: {output_path}")


AVAILABLE_PLOTS = {
    'pressure': plot_pressure_comparison,
    'velocity': plot_velocity_comparison,
    'density': plot_density_comparison,
    'wave_pattern': plot_wave_pattern,
    'tangent_effect': plot_tangent_effect,
    'rezzolla_sweep': plot_rezzolla_vt_sweep,
}


def main():
    parser = argparse.ArgumentParser(
        description="Generate visualization plots for SR tangent velocity tests"
    )
    parser.add_argument(
        '--plot', '-p',
        choices=list(AVAILABLE_PLOTS.keys()) + ['all'],
        default='all',
        help='Plot to generate (default: all)'
    )
    parser.add_argument(
        '--list', '-l',
        action='store_true',
        help='List available plots'
    )
    args = parser.parse_args()
    
    if args.list:
        print("Available plots:")
        for name, func in AVAILABLE_PLOTS.items():
            print(f"  {name}: {func.__doc__.strip().split(chr(10))[0]}")
        return
    
    setup_plot_style()
    output_dir = get_output_dir()
    
    print(f"Output directory: {output_dir}\n")
    
    if args.plot == 'all':
        for name, func in AVAILABLE_PLOTS.items():
            func(output_dir)
    else:
        AVAILABLE_PLOTS[args.plot](output_dir)
    
    print(f"\n✓ All plots saved to: {output_dir}")


if __name__ == "__main__":
    main()
