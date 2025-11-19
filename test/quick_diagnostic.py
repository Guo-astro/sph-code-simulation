#!/usr/bin/env python3
"""
Quick diagnostic script to identify the bug in C++ iterative solver

This script runs a single test case with maximum debug output
to help identify where the C++ and Python solvers diverge.
"""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'relativistic-riemann'))
from kitajima_solver import KitajimaRiemannSolver


def debug_single_iteration():
    """Debug SR Sod test case step-by-step"""
    
    print("="*80)
    print("DETAILED DEBUG: SR Sod Shock Tube")
    print("="*80)
    
    # Test parameters
    gamma = 5.0/3.0
    c = 1.0
    
    # Initial states
    vl, nl, Pl = 0.0, 10.0, 13.33
    vr, nr, Pr = 0.0, 1.0, 1.0e-8
    
    print(f"\nINITIAL STATES:")
    print(f"Left:  v={vl}, n={nl}, P={Pl}")
    print(f"Right: v={vr}, n={nr}, P={Pr}")
    print(f"Pressure ratio: P_L/P_R = {Pl/Pr:.3e}")
    
    # Compute derived quantities
    ul = Pl / ((gamma - 1.0) * nl)
    ur = Pr / ((gamma - 1.0) * nr)
    Hl = 1.0 + ul/(c*c) + Pl/(nl*c*c)
    Hr = 1.0 + ur/(c*c) + Pr/(nr*c*c)
    
    print(f"\nDERIVED QUANTITIES:")
    print(f"u_L = {ul:.6e}, u_R = {ur:.6e}")
    print(f"H_L = {Hl:.10f}, H_R = {Hr:.15f}")
    
    # Compute sound speeds
    csl = np.sqrt(gamma * Pl / (nl * Hl))
    csr = np.sqrt(gamma * Pr / (nr * Hr))
    
    print(f"cs_L = {csl:.10f}, cs_R = {csr:.15f}")
    
    # Lorentz factors (for 1D, no tangential velocity)
    wl = 1.0  # γ = 1 for v=0
    wr = 1.0
    
    print(f"γ_L = {wl:.10f}, γ_R = {wr:.10f}")
    
    # Initial guess (two-rarefaction approximation)
    exponent = (gamma - 1.0) / (2.0 * gamma)
    num = csl + csr - 0.5 * (gamma - 1.0) * (vr - vl)
    denom = csl / (Pl**exponent) + csr / (Pr**exponent)
    pguess = (num / denom)**(2.0 * gamma / (gamma - 1.0))
    
    print(f"\nINITIAL GUESS (Two-Rarefaction):")
    print(f"exponent = {exponent:.10f}")
    print(f"numerator = {num:.10f}")
    print(f"denominator = {denom:.15e}")
    print(f"P_guess = {pguess:.10e}")
    
    # Check if guess is reasonable
    print(f"\nSANITY CHECKS:")
    print(f"P_guess > P_R? {pguess > Pr} ({pguess/Pr:.3e} times)")
    print(f"P_guess < P_L? {pguess < Pl} ({pguess/Pl:.3f} times)")
    
    # Run Python solver for reference
    print(f"\n{'='*80}")
    print("PYTHON SOLVER (Reference)")
    print(f"{'='*80}")
    
    solver = KitajimaRiemannSolver(gamma, c)
    solver.set_initial_states(Pl, nl, vl, Pr, nr, vr)
    
    # Find pressure bounds
    pmin = (Pl + Pr) / 2.0
    pmax = pmin
    
    for _ in range(100):
        pmin = 0.5 * max(pmin, 0.0)
        pmax = 2.0 * pmax
        
        dvel1 = solver.get_dvel(pmin)
        dvel2 = solver.get_dvel(pmax)
        
        if dvel1 * dvel2 <= 0.0:
            break
    
    print(f"Pressure bounds: [{pmin:.6e}, {pmax:.6e}]")
    
    # Solve
    pstar = solver.get_pressure(pmin, pmax, tol=0.0)
    vstar = 0.5 * (solver.vls + solver.vrs)
    
    print(f"\nPYTHON RESULT:")
    print(f"P* = {pstar:.10e}")
    print(f"v* = {vstar:.10e}")
    print(f"v_L* = {solver.vls:.10e}")
    print(f"v_R* = {solver.vrs:.10e}")
    print(f"n_L* = {solver.nls:.10e}")
    print(f"n_R* = {solver.nrs:.10e}")
    
    # Analyze wave types
    print(f"\nWAVE TYPES:")
    print(f"Left wave:  {'SHOCK' if pstar > Pl else 'RAREFACTION'}")
    print(f"Right wave: {'SHOCK' if pstar > Pr else 'RAREFACTION'}")
    
    # Compare with initial guess
    print(f"\nCOMPARISON:")
    print(f"P_guess / P* = {pguess/pstar:.6e}")
    print(f"The initial guess is {'too high' if pguess > pstar else 'too low'}")
    print(f"  by a factor of {abs(pguess - pstar)/pstar * 100:.2f}%")
    
    # Recommend fix
    print(f"\n{'='*80}")
    print("RECOMMENDATION:")
    print(f"{'='*80}")
    
    if pguess/pstar > 10:
        print("❌ Initial guess is WAY too high!")
        print("   Two-rarefaction approximation is inappropriate for this case.")
        print("   Consider using:")
        print(f"   1. Arithmetic mean: ({Pl} + {Pr})/2 = {(Pl+Pr)/2:.6e}")
        print(f"   2. Geometric mean: sqrt({Pl} * {Pr}) = {np.sqrt(Pl*Pr):.6e}")
        print(f"   3. PVRS method (Toro's book)")
    elif pguess/pstar < 0.1:
        print("❌ Initial guess is too low!")
        print("   Consider starting closer to the actual solution.")
    else:
        print("✓ Initial guess is reasonable")
        print("  Bug must be in iteration logic or derivative calculation")
    
    # Test derivative at pguess
    print(f"\n{'='*80}")
    print("DERIVATIVE TEST AT INITIAL GUESS:")
    print(f"{'='*80}")
    
    f_guess = solver.get_dvel(pguess)
    print(f"f(P_guess) = {f_guess:.6e}")
    print(f"This should be close to zero if guess is good")
    
    # Test derivative with finite difference
    dp = 1.0e-6 * pguess
    f_plus = solver.get_dvel(pguess + dp)
    dfdp = (f_plus - f_guess) / dp
    
    print(f"df/dP ≈ {dfdp:.6e}")
    print(f"Newton step: dP = -f/f' = {-f_guess/dfdp:.6e}")
    print(f"New pressure: P_new = {pguess - f_guess/dfdp:.6e}")


if __name__ == '__main__':
    debug_single_iteration()
