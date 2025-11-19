#!/usr/bin/env python3
"""
Extract actual Riemann solver inputs from SPH simulation
and compare C++ vs Python results
"""

import sys
import numpy as np
sys.path.insert(0, '/Users/guo/Downloads/sphcode/relativistic-riemann')
from kitajima_solver import KitajimaRiemannSolver

def test_riemann_solver_on_sph_data():
    """
    Test with actual SPH initial conditions from sr_sod.cpp
    
    From the code:
    - Left:  P=1.0, n=1.0, v=0.0
    - Right: P=0.1, n=0.125, v=0.0
    - gamma = 5/3, c = 1.0
    
    But these are SMOOTHED near x=0 over ~20h width
    """
    
    print("=" * 80)
    print("RIEMANN SOLVER COMPARISON: C++ vs Python on SR Sod States")
    print("=" * 80)
    
    # Sharp discontinuity (theoretical problem)
    print("\n1. SHARP DISCONTINUITY (Theoretical SR Sod):")
    print("-" * 80)
    test_case_sharp = {
        'name': 'Sharp_SR_Sod',
        'gamma': 5.0/3.0,
        'c': 1.0,
        'left': {'v': 0.0, 'n': 1.0, 'P': 1.0, 'vt': 0.0},
        'right': {'v': 0.0, 'n': 0.125, 'P': 0.1, 'vt': 0.0}
    }
    run_comparison(test_case_sharp)
    
    # Smoothed states (what SPH actually uses)
    print("\n2. SMOOTHED STATES (SPH Initial Conditions):")
    print("-" * 80)
    print("Near x=0, states are smoothly interpolated over ~20h width")
    
    # Example: particle pair with smoothed states
    # At x ~ -0.05: mostly left state
    # At x ~ +0.05: mostly right state
    # Effective ratio much smaller than theoretical 10x
    
    # Typical smoothed pair (interpolated)
    s = 0.3  # interpolation factor (example)
    P_smooth_left = 1.0 * (1-s) + 0.1 * s      # ~0.73
    n_smooth_left = 1.0 * (1-s) + 0.125 * s    # ~0.74
    
    s = 0.7
    P_smooth_right = 1.0 * (1-s) + 0.1 * s     # ~0.37
    n_smooth_right = 1.0 * (1-s) + 0.125 * s   # ~0.39
    
    test_case_smooth = {
        'name': 'Smoothed_Near_Contact',
        'gamma': 5.0/3.0,
        'c': 1.0,
        'left': {'v': 0.0, 'n': n_smooth_left, 'P': P_smooth_left, 'vt': 0.0},
        'right': {'v': 0.0, 'n': n_smooth_right, 'P': P_smooth_right, 'vt': 0.0}
    }
    run_comparison(test_case_smooth)
    
    # Moderate ratio example
    print("\n3. MODERATE RATIO (Pressure ratio ~2x):")
    print("-" * 80)
    test_case_moderate = {
        'name': 'Moderate_Ratio',
        'gamma': 5.0/3.0,
        'c': 1.0,
        'left': {'v': 0.0, 'n': 0.6, 'P': 0.8, 'vt': 0.0},
        'right': {'v': 0.0, 'n': 0.4, 'P': 0.4, 'vt': 0.0}
    }
    run_comparison(test_case_moderate)

def run_comparison(test):
    """Run C++ and Python on same inputs and compare"""
    
    # Python solver
    solver = KitajimaRiemannSolver(gamma_c=test['gamma'], c=test['c'])
    solver.set_initial_states(
        Pl=test['left']['P'],
        nl=test['left']['n'],
        vl=test['left']['v'],
        vyl=test['left']['vt'],
        vzl=0.0,
        Pr=test['right']['P'],
        nr=test['right']['n'],
        vr=test['right']['v'],
        vyr=test['right']['vt'],
        vzr=0.0
    )
    
    # Bracketing
    pmin = (solver.Pl + solver.Pr) / 2.0
    pmax = pmin
    
    for i in range(100):
        pmin = 0.5 * max(pmin, 0.0)
        pmax = 2.0 * pmax
        
        dvel1 = solver.get_dvel(pmin)
        dvel2 = solver.get_dvel(pmax)
        
        if dvel1 * dvel2 <= 0.0:
            break
    
    # Solve
    Ps_python = solver.get_pressure(pmin, pmax, tol=1.0e-10)
    _ = solver.get_dvel(Ps_python)  # Update velocities
    vs_python = 0.5 * (solver.vls + solver.vrs)
    
    # For C++, we'll write the input and read expected output
    # (C++ test program should be run separately)
    
    # Calculate pressure ratio
    p_ratio = max(test['left']['P'], test['right']['P']) / min(test['left']['P'], test['right']['P'])
    n_ratio = max(test['left']['n'], test['right']['n']) / min(test['left']['n'], test['right']['n'])
    
    print(f"\n{test['name']}:")
    print(f"  Input states:")
    print(f"    Left:  P={test['left']['P']:.6f}, n={test['left']['n']:.6f}, v={test['left']['v']:.6f}")
    print(f"    Right: P={test['right']['P']:.6f}, n={test['right']['n']:.6f}, v={test['right']['v']:.6f}")
    print(f"  Ratios: P_max/P_min={p_ratio:.2f}, n_max/n_min={n_ratio:.2f}")
    print(f"\n  Python Results:")
    print(f"    Bracket: [{pmin:.6e}, {pmax:.6e}]")
    print(f"    P* = {Ps_python:.16e}")
    print(f"    v* = {vs_python:.16e}")
    print(f"    vls = {solver.vls:.16e}")
    print(f"    vrs = {solver.vrs:.16e}")
    
    # Sound speeds for reference
    csl = np.sqrt(solver.gamma_c * solver.Pl / (solver.nl * solver.Hl))
    csr = np.sqrt(solver.gamma_c * solver.Pr / (solver.nr * solver.Hr))
    print(f"    Sound speeds: cs_left={csl:.6f}, cs_right={csr:.6f}")
    
    return {
        'pstar': Ps_python,
        'vstar': vs_python,
        'input': test
    }

if __name__ == '__main__':
    test_riemann_solver_on_sph_data()
