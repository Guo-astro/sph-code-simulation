#!/usr/bin/env python3
"""
Kernel-Convolved Gravity Comparison
Analyzes energy conservation across 4 test cases:
1. Grad-H ON + Hernquist-Katz (baseline)
2. Grad-H OFF + Hernquist-Katz (expected to have energy drift)
3. Grad-H ON + Wendland C4
4. Grad-H OFF + Wendland C4 (KEY TEST)
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path

# Configuration
RESULTS_DIR = Path("sample/gradh_study_3d/results/kernel_test")
OUTPUT_DIR = Path("sample/gradh_study_3d/figures_v3")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

cases = {
    'gradh_hk': {
        'label': 'Ω ON + H-K',
        'color': 'C0',
        'linestyle': '-',
        'description': 'Grad-h correction ON, Hernquist-Katz softening'
    },
    'no_gradh_hk': {
        'label': 'Ω OFF + H-K',
        'color': 'C1',
        'linestyle': '--',
        'description': 'Grad-h correction OFF, Hernquist-Katz softening'
    },
    'gradh_wendland': {
        'label': 'Ω ON + Wendland',
        'color': 'C2',
        'linestyle': '-',
        'description': 'Grad-h correction ON, Wendland C4 softening'
    },
    'no_gradh_wendland': {
        'label': 'Ω OFF + Wendland',
        'color': 'C3',
        'linestyle': '--',
        'description': 'Grad-h correction OFF, Wendland C4 softening (KEY TEST)'
    }
}

def load_energy_file(case_dir):
    """Load energy data from energy.dat file."""
    energy_file = RESULTS_DIR / case_dir / "energy.dat"
    if not energy_file.exists():
        print(f"Warning: {energy_file} not found")
        return None
    
    data = np.loadtxt(energy_file, comments='#')
    return {
        'time': data[:, 0],
        'kinetic': data[:, 1],
        'thermal': data[:, 2],
        'potential': data[:, 3],
        'total': data[:, 4]
    }

def main():
    # Load all data
    data = {}
    for case_name in cases:
        d = load_energy_file(case_name)
        if d is not None:
            data[case_name] = d
            print(f"Loaded {case_name}: {len(d['time'])} data points")
    
    if not data:
        print("No data loaded!")
        return
    
    # Create figure with 4 subplots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # =====================
    # Plot 1: Total Energy vs Time
    # =====================
    ax1 = axes[0, 0]
    for case_name, case_info in cases.items():
        if case_name in data:
            d = data[case_name]
            ax1.plot(d['time'], d['total'], 
                    label=case_info['label'],
                    color=case_info['color'],
                    linestyle=case_info['linestyle'],
                    linewidth=2)
    
    ax1.set_xlabel('Time [code units]', fontsize=12)
    ax1.set_ylabel('Total Energy [code units]', fontsize=12)
    ax1.set_title('Total Energy Evolution', fontsize=14)
    ax1.legend(loc='best', fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # =====================
    # Plot 2: Energy Conservation Error (relative to initial)
    # =====================
    ax2 = axes[0, 1]
    for case_name, case_info in cases.items():
        if case_name in data:
            d = data[case_name]
            E0 = d['total'][0]
            error = (d['total'] - E0) / abs(E0) * 100  # percent
            ax2.plot(d['time'], error, 
                    label=case_info['label'],
                    color=case_info['color'],
                    linestyle=case_info['linestyle'],
                    linewidth=2)
    
    ax2.set_xlabel('Time [code units]', fontsize=12)
    ax2.set_ylabel('Energy Error [%]', fontsize=12)
    ax2.set_title('Energy Conservation Error', fontsize=14)
    ax2.legend(loc='best', fontsize=10)
    ax2.grid(True, alpha=0.3)
    ax2.axhline(y=0, color='k', linestyle=':', alpha=0.5)
    
    # =====================
    # Plot 3: Energy Components
    # =====================
    ax3 = axes[1, 0]
    case_name = 'gradh_hk'  # Show baseline case
    if case_name in data:
        d = data[case_name]
        ax3.plot(d['time'], d['kinetic'], label='Kinetic', color='C0', linewidth=2)
        ax3.plot(d['time'], d['thermal'], label='Thermal', color='C1', linewidth=2)
        ax3.plot(d['time'], d['potential'], label='Gravitational', color='C2', linewidth=2)
        ax3.plot(d['time'], d['total'], label='Total', color='C3', linewidth=2, linestyle='--')
    
    ax3.set_xlabel('Time [code units]', fontsize=12)
    ax3.set_ylabel('Energy [code units]', fontsize=12)
    ax3.set_title(f'Energy Components ({cases[case_name]["label"]})', fontsize=14)
    ax3.legend(loc='best', fontsize=10)
    ax3.grid(True, alpha=0.3)
    
    # =====================
    # Plot 4: Bar chart of final energy error
    # =====================
    ax4 = axes[1, 1]
    
    case_names = []
    errors = []
    colors = []
    
    for case_name, case_info in cases.items():
        if case_name in data:
            d = data[case_name]
            E0 = d['total'][0]
            Ef = d['total'][-1]
            error = abs((Ef - E0) / E0) * 100
            case_names.append(case_info['label'])
            errors.append(error)
            colors.append(case_info['color'])
    
    bars = ax4.bar(case_names, errors, color=colors, alpha=0.7, edgecolor='black')
    ax4.set_ylabel('Energy Error [%]', fontsize=12)
    ax4.set_title('Final Energy Conservation Error', fontsize=14)
    ax4.grid(True, alpha=0.3, axis='y')
    
    # Add value labels on bars
    for bar, error in zip(bars, errors):
        height = bar.get_height()
        ax4.annotate(f'{error:.3f}%',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=10)
    
    plt.tight_layout()
    
    # Save figure
    output_path = OUTPUT_DIR / "kernel_gravity_comparison.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"\n✓ Saved figure to: {output_path}")
    
    # Print summary
    print("\n" + "="*60)
    print("KERNEL-CONVOLVED GRAVITY TEST RESULTS")
    print("="*60)
    print("\nTest Hypothesis: Using true kernel-convolved gravity")
    print("(Wendland C4) should improve energy conservation even")
    print("without grad-h correction (Ω).")
    print("\nResults Summary:")
    print("-"*60)
    
    for case_name, case_info in cases.items():
        if case_name in data:
            d = data[case_name]
            E0 = d['total'][0]
            Ef = d['total'][-1]
            error = (Ef - E0) / abs(E0) * 100
            print(f"{case_info['label']:20s}: ΔE/|E₀| = {error:+.4f}%")
    
    print("-"*60)
    print("\nConclusion:")
    
    # Analyze results
    if 'gradh_hk' in data and 'no_gradh_hk' in data:
        e_gradh_hk = abs((data['gradh_hk']['total'][-1] - data['gradh_hk']['total'][0]) / data['gradh_hk']['total'][0])
        e_no_gradh_hk = abs((data['no_gradh_hk']['total'][-1] - data['no_gradh_hk']['total'][0]) / data['no_gradh_hk']['total'][0])
        
    if 'gradh_wendland' in data and 'no_gradh_wendland' in data:
        e_gradh_w = abs((data['gradh_wendland']['total'][-1] - data['gradh_wendland']['total'][0]) / data['gradh_wendland']['total'][0])
        e_no_gradh_w = abs((data['no_gradh_wendland']['total'][-1] - data['no_gradh_wendland']['total'][0]) / data['no_gradh_wendland']['total'][0])
        
        if e_no_gradh_w < e_no_gradh_hk:
            print("✓ Wendland C4 softening shows BETTER energy conservation")
            print("  without grad-h correction than Hernquist-Katz.")
        else:
            print("✗ Wendland C4 softening does NOT improve energy conservation")
            print("  without grad-h correction compared to Hernquist-Katz.")
            
        if e_no_gradh_w < e_gradh_hk * 2:
            print("✓ Wendland C4 without Ω achieves comparable accuracy to")
            print("  Hernquist-Katz with Ω.")
        else:
            print("✗ Wendland C4 without Ω still shows significant energy drift.")
    
    plt.show()

if __name__ == "__main__":
    main()
