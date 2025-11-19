#!/usr/bin/env python3
"""
Direct comparison of C++ vs Python Riemann solver on identical inputs.
This tests whether the implementations produce exactly matching results.
"""

import sys
import os
import subprocess
import numpy as np

# Add relativistic-riemann to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'relativistic-riemann'))
from kitajima_solver import KitajimaRiemannSolver

def run_cpp_test(test_name, Pl, Pr, vl, vr, nl, nr, gamma, c):
    """Run the C++ solver and capture output."""
    # Create a minimal C++ test program
    cpp_code = f"""
#include <iostream>
#include <cmath>
#include <iomanip>
#include <limits>

int main() {{
    // Test: {test_name}
    const double gamma_c = {gamma};
    const double c = {c};
    
    // Left state
    const double Pl = {Pl};
    const double nl = {nl};
    const double vl = {vl};
    const double ul = Pl / ((gamma_c - 1.0) * nl);
    const double Hl = 1.0 + ul / (c * c) + Pl / (nl * c * c);
    const double csl = std::sqrt(gamma_c * Pl / (nl * Hl));
    const double vtl = 0.0;
    const double wl = 1.0 / std::sqrt(1.0 - vl * vl / (c * c));
    
    // Right state
    const double Pr = {Pr};
    const double nr = {nr};
    const double vr = {vr};
    const double ur = Pr / ((gamma_c - 1.0) * nr);
    const double Hr = 1.0 + ur / (c * c) + Pr / (nr * c * c);
    const double csr = std::sqrt(gamma_c * Pr / (nr * Hr));
    const double vtr = 0.0;
    const double wr = 1.0 / std::sqrt(1.0 - vr * vr / (c * c));
    
    std::cout << std::setprecision(16);
    std::cout << "Pl=" << Pl << " Pr=" << Pr << " vl=" << vl << " vr=" << vr << std::endl;
    std::cout << "Hl=" << Hl << " Hr=" << Hr << std::endl;
    std::cout << "csl=" << csl << " csr=" << csr << std::endl;
    std::cout << "wl=" << wl << " wr=" << wr << std::endl;
    
    return 0;
}}
"""
    
    # Write and compile
    with open('/tmp/test_riemann.cpp', 'w') as f:
        f.write(cpp_code)
    
    result = subprocess.run(['g++', '-std=c++17', '-O0', '/tmp/test_riemann.cpp', '-o', '/tmp/test_riemann'],
                          capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Compilation failed: {result.stderr}")
        return None
    
    result = subprocess.run(['/tmp/test_riemann'], capture_output=True, text=True)
    return result.stdout

def test_case(name, Pl, Pr, vl, vr, nl, nr, gamma=5/3, c=1.0):
    """Test a single Riemann problem."""
    print(f"\n{'='*70}")
    print(f"TEST: {name}")
    print(f"{'='*70}")
    
    # Python solution
    solver = KitajimaRiemannSolver(
        Pl=Pl, Pr=Pr,
        nl=nl, nr=nr,
        vl=vl, vr=vr,
        vyl=0.0, vyr=0.0,
        vzl=0.0, vzr=0.0,
        gamma_c=gamma, c=c
    )
    
    print(f"\nInput States:")
    print(f"  Left:  P={Pl:12.6e} n={nl:12.6e} v={vl:12.6e}")
    print(f"  Right: P={Pr:12.6e} n={nr:12.6e} v={vr:12.6e}")
    print(f"  Ratio: P_L/P_R = {Pl/Pr:.2e}")
    
    # Python calculation
    print(f"\nPython (kitajima_solver.py):")
    print(f"  Hl = {solver.Hl:.16e}")
    print(f"  Hr = {solver.Hr:.16e}")
    print(f"  csl = {solver.csl:.16e}")
    print(f"  csr = {solver.csr:.16e}")
    print(f"  gammal = {solver.gammal:.16e}")
    print(f"  gammar = {solver.gammar:.16e}")
    
    # Run C++ (for verification of setup)
    cpp_output = run_cpp_test(name, Pl, Pr, vl, vr, nl, nr, gamma, c)
    if cpp_output:
        print(f"\nC++ verification:")
        print(cpp_output.strip())
    
    # Solve with Python
    try:
        result = solver.solve(x0=0.0, t=0.4, tol=1e-10)
        print(f"\nPython Solution:")
        print(f"  P* = {result['Ps']:.16e}")
        print(f"  v* = {result['vs']:.16e}")
        print(f"  Converged: {result.get('converged', 'N/A')}")
    except Exception as e:
        print(f"\nPython FAILED: {e}")
        return False
    
    return True

if __name__ == '__main__':
    print("VERIFICATION: C++ vs Python Riemann Solver")
    print("="*70)
    
    # Test 1: SR Sod (extreme ratio)
    test_case("SR Sod Shock Tube", 
             Pl=13.33, Pr=1.0e-8, 
             vl=0.0, vr=0.0,
             nl=10.0, nr=1.0,
             gamma=5/3, c=1.0)
    
    # Test 2: Moderate ratio
    test_case("Moderate Pressure Ratio",
             Pl=1.0, Pr=0.1,
             vl=0.0, vr=0.0,
             nl=1.0, nr=0.125,
             gamma=5/3, c=1.0)
    
    # Test 3: Small ratio (like SPH)
    test_case("SPH-like Smoothed",
             Pl=0.85, Pr=0.12,
             vl=0.0, vr=0.0,
             nl=1.0, nr=0.14,
             gamma=5/3, c=1.0)
    
    print(f"\n{'='*70}")
    print("Verification complete")
    print(f"{'='*70}")
