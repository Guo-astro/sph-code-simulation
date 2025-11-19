#!/usr/bin/env python3
"""
Visualize comparison between C++ and Python Riemann solver results

This script parses the output files from both solvers and creates
detailed comparison plots showing convergence behavior and final results.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import re
import os


def parse_cpp_results(filename):
    """Parse C++ test results file"""
    tests = []
    
    with open(filename, 'r') as f:
        content = f.read()
    
    # Split into test cases
    test_blocks = re.split(r'={80,}', content)
    
    for block in test_blocks:
        if not block.strip():
            continue
        
        # Extract test name
        name_match = re.search(r'TEST:\s+(\S+)', block)
        if not name_match:
            continue
        
        name = name_match.group(1)
        
        # Extract parameters
        gamma_match = re.search(r'gamma=([0-9.e+-]+)', block)
        c_match = re.search(r'c=([0-9.e+-]+)', block)
        
        # Extract states
        left_match = re.search(r'LEFT:\s+v=([0-9.e+-]+),\s+n=([0-9.e+-]+),\s+P=([0-9.e+-]+)', block)
        right_match = re.search(r'RIGHT:\s+v=([0-9.e+-]+),\s+n=([0-9.e+-]+),\s+P=([0-9.e+-]+)', block)
        
        # Extract results
        pstar_match = re.search(r'pstar=([0-9.e+-]+)', block)
        vstar_match = re.search(r'vstar=([0-9.e+-]+)', block)
        iter_match = re.search(r'iterations=([0-9]+)', block)
        conv_match = re.search(r'converged=(true|false)', block)
        resid_match = re.search(r'final_residual=([0-9.e+-]+)', block)
        
        # Extract iteration history
        history_section = re.search(r'ITERATION_HISTORY:\s*\n(.*?)(?:={80,}|$)', block, re.DOTALL)
        p_history = []
        f_history = []
        
        if history_section:
            for line in history_section.group(1).strip().split('\n'):
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        p_history.append(float(parts[1]))
                        f_history.append(float(parts[2]))
                    except ValueError:
                        pass
        
        if all([gamma_match, c_match, left_match, right_match, pstar_match, vstar_match]):
            tests.append({
                'name': name,
                'gamma': float(gamma_match.group(1)),
                'c': float(c_match.group(1)),
                'left': [float(x) for x in left_match.groups()],
                'right': [float(x) for x in right_match.groups()],
                'pstar': float(pstar_match.group(1)),
                'vstar': float(vstar_match.group(1)),
                'iterations': int(iter_match.group(1)) if iter_match else 0,
                'converged': conv_match.group(1) == 'true' if conv_match else False,
                'residual': float(resid_match.group(1)) if resid_match else 0.0,
                'p_history': np.array(p_history),
                'f_history': np.array(f_history)
            })
    
    return tests


def parse_python_results(filename):
    """Parse Python test results file"""
    tests = []
    
    with open(filename, 'r') as f:
        content = f.read()
    
    # Split into test cases
    test_blocks = re.split(r'={80,}', content)
    
    for block in test_blocks:
        if not block.strip():
            continue
        
        # Extract test name
        name_match = re.search(r'TEST:\s+(\S+)', block)
        if not name_match:
            continue
        
        name = name_match.group(1)
        
        # Extract parameters
        gamma_match = re.search(r'gamma=([0-9.e+-]+)', block)
        c_match = re.search(r'c=([0-9.e+-]+)', block)
        
        # Extract states
        left_match = re.search(r'LEFT:\s+v=([0-9.e+-]+),\s+n=([0-9.e+-]+),\s+P=([0-9.e+-]+)', block)
        right_match = re.search(r'RIGHT:\s+v=([0-9.e+-]+),\s+n=([0-9.e+-]+),\s+P=([0-9.e+-]+)', block)
        
        # Extract results
        pstar_match = re.search(r'pstar=([0-9.e+-]+)', block)
        vstar_match = re.search(r'vstar=([0-9.e+-]+)', block)
        
        if all([gamma_match, c_match, left_match, right_match, pstar_match, vstar_match]):
            tests.append({
                'name': name,
                'gamma': float(gamma_match.group(1)),
                'c': float(c_match.group(1)),
                'left': [float(x) for x in left_match.groups()],
                'right': [float(x) for x in right_match.groups()],
                'pstar': float(pstar_match.group(1)),
                'vstar': float(vstar_match.group(1))
            })
    
    return tests


def compare_and_visualize(cpp_tests, python_tests, output_dir='comparison_plots'):
    """Create comparison visualizations"""
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Create summary table
    print("\n" + "="*100)
    print("COMPARISON SUMMARY")
    print("="*100)
    print(f"{'Test Name':<30} {'C++ P*':<18} {'Py P*':<18} {'ΔP* (%)':<12} {'C++ v*':<18} {'Py v*':<18} {'Δv*':<12}")
    print("-"*100)
    
    comparison_data = []
    
    for cpp_test in cpp_tests:
        # Find matching Python test
        py_test = next((t for t in python_tests if t['name'] == cpp_test['name']), None)
        
        if py_test is None:
            print(f"{cpp_test['name']:<30} NO PYTHON MATCH")
            continue
        
        # Calculate differences
        pstar_diff = abs(cpp_test['pstar'] - py_test['pstar'])
        pstar_rel = 100 * pstar_diff / max(abs(cpp_test['pstar']), 1e-10)
        
        vstar_diff = abs(cpp_test['vstar'] - py_test['vstar'])
        
        print(f"{cpp_test['name']:<30} "
              f"{cpp_test['pstar']:<18.10e} "
              f"{py_test['pstar']:<18.10e} "
              f"{pstar_rel:<12.6e} "
              f"{cpp_test['vstar']:<18.10e} "
              f"{py_test['vstar']:<18.10e} "
              f"{vstar_diff:<12.6e}")
        
        comparison_data.append({
            'cpp': cpp_test,
            'python': py_test,
            'pstar_diff': pstar_diff,
            'pstar_rel': pstar_rel,
            'vstar_diff': vstar_diff
        })
    
    print("="*100)
    
    # Create detailed plots for each test
    for i, comp in enumerate(comparison_data):
        cpp_test = comp['cpp']
        py_test = comp['python']
        
        fig = plt.figure(figsize=(16, 10))
        gs = GridSpec(3, 2, figure=fig, hspace=0.3, wspace=0.3)
        
        test_name = cpp_test['name']
        
        # Main title
        fig.suptitle(f"Test: {test_name}\n"
                    f"C++ P*={cpp_test['pstar']:.6e}, Python P*={py_test['pstar']:.6e}, "
                    f"ΔP*={comp['pstar_rel']:.3e}%\n"
                    f"C++ v*={cpp_test['vstar']:.6e}, Python v*={py_test['vstar']:.6e}, "
                    f"Δv*={comp['vstar_diff']:.3e}",
                    fontsize=12, fontweight='bold')
        
        # Plot 1: Convergence history (pressure)
        ax1 = fig.add_subplot(gs[0, 0])
        if len(cpp_test['p_history']) > 0:
            ax1.plot(cpp_test['p_history'], 'b.-', label='C++ P*', linewidth=2, markersize=8)
            ax1.axhline(py_test['pstar'], color='r', linestyle='--', label='Python P*', linewidth=2)
            ax1.set_xlabel('Iteration', fontsize=10)
            ax1.set_ylabel('Pressure P*', fontsize=10)
            ax1.set_title('Convergence: Pressure', fontsize=11, fontweight='bold')
            ax1.legend()
            ax1.grid(True, alpha=0.3)
        else:
            ax1.text(0.5, 0.5, 'No convergence history', 
                    ha='center', va='center', transform=ax1.transAxes)
            ax1.set_title('Convergence: Pressure', fontsize=11, fontweight='bold')
        
        # Plot 2: Convergence history (residual)
        ax2 = fig.add_subplot(gs[0, 1])
        if len(cpp_test['f_history']) > 0:
            ax2.semilogy(np.abs(cpp_test['f_history']), 'b.-', 
                        label='C++ |f(P*)|', linewidth=2, markersize=8)
            ax2.axhline(1e-8, color='g', linestyle='--', label='Tolerance', linewidth=2)
            ax2.set_xlabel('Iteration', fontsize=10)
            ax2.set_ylabel('|Velocity Residual|', fontsize=10)
            ax2.set_title('Convergence: Residual', fontsize=11, fontweight='bold')
            ax2.legend()
            ax2.grid(True, alpha=0.3)
        else:
            ax2.text(0.5, 0.5, 'No convergence history', 
                    ha='center', va='center', transform=ax2.transAxes)
            ax2.set_title('Convergence: Residual', fontsize=11, fontweight='bold')
        
        # Plot 3: Initial conditions comparison
        ax3 = fig.add_subplot(gs[1, :])
        
        states_table = [
            ['State', 'Source', 'v', 'n', 'P'],
            ['Left', 'C++', f"{cpp_test['left'][0]:.6e}", 
             f"{cpp_test['left'][1]:.6e}", f"{cpp_test['left'][2]:.6e}"],
            ['Left', 'Python', f"{py_test['left'][0]:.6e}", 
             f"{py_test['left'][1]:.6e}", f"{py_test['left'][2]:.6e}"],
            ['Right', 'C++', f"{cpp_test['right'][0]:.6e}", 
             f"{cpp_test['right'][1]:.6e}", f"{cpp_test['right'][2]:.6e}"],
            ['Right', 'Python', f"{py_test['right'][0]:.6e}", 
             f"{py_test['right'][1]:.6e}", f"{py_test['right'][2]:.6e}"],
        ]
        
        ax3.axis('tight')
        ax3.axis('off')
        table = ax3.table(cellText=states_table, cellLoc='center', loc='center',
                         colWidths=[0.15, 0.15, 0.23, 0.23, 0.23])
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1, 2)
        
        # Header row formatting
        for i in range(5):
            table[(0, i)].set_facecolor('#40466e')
            table[(0, i)].set_text_props(weight='bold', color='white')
        
        ax3.set_title('Initial Conditions', fontsize=11, fontweight='bold', pad=20)
        
        # Plot 4: Results comparison
        ax4 = fig.add_subplot(gs[2, 0])
        
        results_table = [
            ['Variable', 'C++', 'Python', 'Difference'],
            ['P* (pstar)', f"{cpp_test['pstar']:.10e}", f"{py_test['pstar']:.10e}", 
             f"{comp['pstar_diff']:.3e}"],
            ['v* (vstar)', f"{cpp_test['vstar']:.10e}", f"{py_test['vstar']:.10e}", 
             f"{comp['vstar_diff']:.3e}"],
            ['Iterations', f"{cpp_test['iterations']}", "N/A", "N/A"],
            ['Converged', f"{cpp_test['converged']}", "N/A", "N/A"],
        ]
        
        ax4.axis('tight')
        ax4.axis('off')
        table2 = ax4.table(cellText=results_table, cellLoc='center', loc='center',
                          colWidths=[0.25, 0.25, 0.25, 0.25])
        table2.auto_set_font_size(False)
        table2.set_fontsize(9)
        table2.scale(1, 2)
        
        # Header row formatting
        for i in range(4):
            table2[(0, i)].set_facecolor('#40466e')
            table2[(0, i)].set_text_props(weight='bold', color='white')
        
        ax4.set_title('Results Comparison', fontsize=11, fontweight='bold', pad=20)
        
        # Plot 5: Error metrics
        ax5 = fig.add_subplot(gs[2, 1])
        
        metrics = ['P* Rel Error (%)', 'v* Abs Error']
        values = [comp['pstar_rel'], comp['vstar_diff']]
        colors = ['#3498db' if v < 1e-6 else '#e74c3c' for v in [comp['pstar_rel'], comp['vstar_diff']]]
        
        bars = ax5.bar(metrics, values, color=colors, alpha=0.7, edgecolor='black', linewidth=1.5)
        ax5.set_ylabel('Error Magnitude', fontsize=10)
        ax5.set_title('Error Metrics', fontsize=11, fontweight='bold')
        ax5.set_yscale('log')
        ax5.grid(True, alpha=0.3, axis='y')
        
        # Add value labels on bars
        for bar, val in zip(bars, values):
            height = bar.get_height()
            ax5.text(bar.get_x() + bar.get_width()/2., height,
                    f'{val:.2e}', ha='center', va='bottom', fontsize=9, fontweight='bold')
        
        # Save figure
        output_file = os.path.join(output_dir, f'{test_name}_comparison.png')
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Saved: {output_file}")
        plt.close()
    
    # Create summary comparison plot
    fig, axes = plt.subplots(2, 1, figsize=(14, 10))
    
    test_names = [c['cpp']['name'] for c in comparison_data]
    pstar_errors = [c['pstar_rel'] for c in comparison_data]
    vstar_errors = [c['vstar_diff'] for c in comparison_data]
    
    # P* relative errors
    axes[0].bar(range(len(test_names)), pstar_errors, color='#3498db', alpha=0.7, 
               edgecolor='black', linewidth=1.5)
    axes[0].set_xticks(range(len(test_names)))
    axes[0].set_xticklabels(test_names, rotation=45, ha='right')
    axes[0].set_ylabel('Relative Error (%)', fontsize=11)
    axes[0].set_title('Pressure P* Relative Error: C++ vs Python', fontsize=12, fontweight='bold')
    axes[0].set_yscale('log')
    axes[0].grid(True, alpha=0.3, axis='y')
    
    # v* absolute errors
    axes[1].bar(range(len(test_names)), vstar_errors, color='#e74c3c', alpha=0.7,
               edgecolor='black', linewidth=1.5)
    axes[1].set_xticks(range(len(test_names)))
    axes[1].set_xticklabels(test_names, rotation=45, ha='right')
    axes[1].set_ylabel('Absolute Error', fontsize=11)
    axes[1].set_title('Velocity v* Absolute Error: C++ vs Python', fontsize=12, fontweight='bold')
    axes[1].set_yscale('log')
    axes[1].grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    summary_file = os.path.join(output_dir, 'summary_comparison.png')
    plt.savefig(summary_file, dpi=150, bbox_inches='tight')
    print(f"\nSaved summary: {summary_file}")
    plt.close()
    
    print(f"\nAll plots saved to: {output_dir}/")


def main():
    """Main comparison and visualization routine"""
    print("Parsing C++ results...")
    cpp_tests = parse_cpp_results('test_results_cpp.txt')
    print(f"Found {len(cpp_tests)} C++ test results")
    
    print("\nParsing Python results...")
    python_tests = parse_python_results('test_results_python.txt')
    print(f"Found {len(python_tests)} Python test results")
    
    if not cpp_tests or not python_tests:
        print("\nError: Could not find test results.")
        print("Please run the following first:")
        print("1. cd test && g++ -std=c++17 -O3 test_iterative_riemann_solver.cpp -o test_solver && ./test_solver")
        print("2. python3 compare_riemann_solvers.py")
        return
    
    print("\nCreating comparison visualizations...")
    compare_and_visualize(cpp_tests, python_tests)
    
    print("\n✓ Comparison complete!")


if __name__ == '__main__':
    main()
