#!/usr/bin/env python3
"""
Compare C++ iterative Riemann solver with Python Kitajima solver

This script runs the same test cases through both the C++ and Python
implementations and compares the results to identify any discrepancies.
"""

import numpy as np
import sys
import os

# Add relativistic-riemann to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'relativistic-riemann'))

from kitajima_solver import KitajimaRiemannSolver


class TestCase:
    """Test case structure matching C++ version"""
    def __init__(self, name, gamma, c, left, right):
        self.name = name
        self.gamma = gamma
        self.c = c
        self.left = left  # [v^n, n, P, cs, |v^t|]
        self.right = right


def compute_sound_speed(P, n, gamma, c):
    """Compute sound speed from primitive variables"""
    u = P / ((gamma - 1.0) * n)
    H = 1.0 + u/(c*c) + P/(n*c*c)
    cs = np.sqrt(gamma * P / (n * H))
    return cs


def run_python_solver(tc):
    """Run Kitajima solver on test case"""
    print(f"\n{'='*80}")
    print(f"PYTHON TEST: {tc.name}")
    print(f"{'='*80}")
    
    # Extract state vectors
    vl, nl, Pl, csl, vtl = tc.left
    vr, nr, Pr, csr, vtr = tc.right
    
    print(f"gamma={tc.gamma}, c={tc.c}")
    print(f"LEFT:  v={vl}, n={nl}, P={Pl}, cs={csl}, vt={vtl}")
    print(f"RIGHT: v={vr}, n={nr}, P={Pr}, cs={csr}, vt={vtr}")
    
    # Create solver
    solver = KitajimaRiemannSolver(tc.gamma, tc.c)
    
    # Set initial states (note: Python solver takes (P, n, v) order)
    # For 1D problem, tangential velocities are in y,z components
    # We'll set vy = vt, vz = 0 for simplicity
    solver.set_initial_states(
        Pl, nl, vl,  # Left state
        Pr, nr, vr,  # Right state
        vyl=vtl, vzl=0.0,  # Left tangential
        vyr=vtr, vzr=0.0   # Right tangential
    )
    
    # Find pressure bounds (matching C++ logic)
    pmin = (Pl + Pr) / 2.0
    pmax = pmin
    
    for _ in range(100):
        pmin = 0.5 * max(pmin, 0.0)
        pmax = 2.0 * pmax
        
        dvel1 = solver.get_dvel(pmin)
        dvel2 = solver.get_dvel(pmax)
        
        if dvel1 * dvel2 <= 0.0:
            break
    
    print(f"\nPressure bounds: [{pmin}, {pmax}]")
    
    # Solve for star pressure
    pstar = solver.get_pressure(pmin, pmax, tol=0.0)
    vstar = 0.5 * (solver.vls + solver.vrs)
    
    print(f"\nRESULT:")
    print(f"pstar={pstar}")
    print(f"vstar={vstar}")
    print(f"vls={solver.vls}, vrs={solver.vrs}")
    print(f"nls={solver.nls}, nrs={solver.nrs}")
    
    return {
        'pstar': pstar,
        'vstar': vstar,
        'vls': solver.vls,
        'vrs': solver.vrs,
        'nls': solver.nls,
        'nrs': solver.nrs,
        'Hls': solver.Hls,
        'Hrs': solver.Hrs
    }


def main():
    """Run comparison tests"""
    print("Relativistic Riemann Solver Comparison: C++ vs Python")
    print("="*80)
    
    # Define test cases (matching C++ version)
    test_cases = []
    
    # Test 1: SR Sod shock tube
    tc = TestCase(
        "SR_Sod_Shock_Tube",
        5.0/3.0,  # gamma
        1.0,      # c
        [0.0, 10.0, 13.33, 0.0, 0.0],  # left: v, n, P, cs, vt
        [0.0, 1.0, 1.0e-8, 0.0, 0.0]   # right
    )
    # Compute sound speeds
    tc.left[3] = compute_sound_speed(tc.left[2], tc.left[1], tc.gamma, tc.c)
    tc.right[3] = compute_sound_speed(tc.right[2], tc.right[1], tc.gamma, tc.c)
    test_cases.append(tc)
    
    # Test 2: Two rarefactions
    tc = TestCase(
        "Two_Rarefactions",
        5.0/3.0,
        1.0,
        [-0.5, 1.0, 1.0, 0.0, 0.0],
        [0.5, 1.0, 1.0, 0.0, 0.0]
    )
    tc.left[3] = compute_sound_speed(tc.left[2], tc.left[1], tc.gamma, tc.c)
    tc.right[3] = compute_sound_speed(tc.right[2], tc.right[1], tc.gamma, tc.c)
    test_cases.append(tc)
    
    # Test 3: Two shocks
    tc = TestCase(
        "Two_Shocks",
        5.0/3.0,
        1.0,
        [0.5, 1.0, 1.0, 0.0, 0.0],
        [-0.5, 1.0, 1.0, 0.0, 0.0]
    )
    tc.left[3] = compute_sound_speed(tc.left[2], tc.left[1], tc.gamma, tc.c)
    tc.right[3] = compute_sound_speed(tc.right[2], tc.right[1], tc.gamma, tc.c)
    test_cases.append(tc)
    
    # Test 4: Left shock, right rarefaction
    tc = TestCase(
        "Left_Shock_Right_Rarefaction",
        5.0/3.0,
        1.0,
        [0.3, 1.0, 10.0, 0.0, 0.0],
        [0.0, 1.0, 1.0, 0.0, 0.0]
    )
    tc.left[3] = compute_sound_speed(tc.left[2], tc.left[1], tc.gamma, tc.c)
    tc.right[3] = compute_sound_speed(tc.right[2], tc.right[1], tc.gamma, tc.c)
    test_cases.append(tc)
    
    # Test 5: With tangential velocity
    tc = TestCase(
        "With_Tangential_Velocity",
        5.0/3.0,
        1.0,
        [0.0, 1.0, 1.0, 0.0, 0.3],  # vt = 0.3
        [0.0, 1.0, 0.1, 0.0, 0.0]
    )
    tc.left[3] = compute_sound_speed(tc.left[2], tc.left[1], tc.gamma, tc.c)
    tc.right[3] = compute_sound_speed(tc.right[2], tc.right[1], tc.gamma, tc.c)
    test_cases.append(tc)
    
    # Run Python tests and save results
    results = []
    with open('test_results_python.txt', 'w') as f:
        for tc in test_cases:
            result = run_python_solver(tc)
            results.append({'name': tc.name, 'result': result})
            
            # Write to file
            f.write(f"\n{'='*80}\n")
            f.write(f"TEST: {tc.name}\n")
            f.write(f"gamma={tc.gamma}, c={tc.c}\n")
            f.write(f"LEFT:  v={tc.left[0]}, n={tc.left[1]}, P={tc.left[2]}, "
                   f"cs={tc.left[3]}, vt={tc.left[4]}\n")
            f.write(f"RIGHT: v={tc.right[0]}, n={tc.right[1]}, P={tc.right[2]}, "
                   f"cs={tc.right[3]}, vt={tc.right[4]}\n")
            f.write(f"\nRESULT:\n")
            f.write(f"pstar={result['pstar']}\n")
            f.write(f"vstar={result['vstar']}\n")
            f.write(f"vls={result['vls']}\n")
            f.write(f"vrs={result['vrs']}\n")
            f.write(f"nls={result['nls']}\n")
            f.write(f"nrs={result['nrs']}\n")
            f.write(f"{'='*80}\n")
    
    print("\n\nResults written to test_results_python.txt")
    print("\nTo compare with C++ results:")
    print("1. Compile and run: cd test && g++ -std=c++17 -O3 test_iterative_riemann_solver.cpp -o test_solver && ./test_solver")
    print("2. Compare: diff test_results_cpp.txt test_results_python.txt")
    print("3. Or run: python3 visualize_comparison.py")


if __name__ == '__main__':
    main()
