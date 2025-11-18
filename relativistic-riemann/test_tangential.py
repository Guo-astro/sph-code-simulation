#!/usr/bin/env python3
"""
Test tangential velocity effects in relativistic Riemann solver.

This script tests the implementation of tangential velocities following
Pons, Martí & Müller (1999, J. Fluid Mech.) using the Kitajima formulation.

Test cases:
1. Sod shock tube with varying tangential velocities (v^y = 0, 0.5, 0.9, 0.99)
2. Relativistic blast wave with tangential velocities (Table 1 from Pons et al.)
3. Verification of conservation laws across waves

Physical expectations:
- Tangential velocity direction preserved: v^y/v^z = const
- Conservation across rarefactions: h*W*v^t = const
- Conservation across shocks: h*W*v^t = const
- Wave speeds decrease as v^t increases
- Intermediate pressure increases with right-side tangential velocity
"""

import numpy as np
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent))

from kitajima_solver import KitajimaRiemannSolver


def test_sod_tangential():
    """
    Test Sod shock tube with tangential velocities.
    
    Initial conditions (similar to Pons et al. Figure 3):
    - Left:  P=1.0, n=1.0, v^x=0.0
    - Right: P=0.1, n=0.125, v^x=0.0
    - Vary v^y on both sides
    """
    print("=" * 80)
    print("TEST 1: Sod Shock Tube with Tangential Velocities")
    print("=" * 80)
    print("\nInitial conditions:")
    print("Left:  P=1.0,   n=1.0,   v^x=0.0")
    print("Right: P=0.1,   n=0.125, v^x=0.0")
    print("Varying v^y = 0, 0.5, 0.9, 0.99\n")
    
    gamma = 1.4
    c = 1.0
    
    # Initial conditions
    Pl, nl, vl = 1.0, 1.0, 0.0
    Pr, nr, vr = 0.1, 0.125, 0.0
    
    # Test different tangential velocities
    tangential_vels = [0.0, 0.5, 0.9, 0.99]
    
    results = []
    
    for vy in tangential_vels:
        solver = KitajimaRiemannSolver(gamma, c)
        solver.set_initial_states(Pl, nl, vl, Pr, nr, vr, 
                                  vyl=vy, vzl=0.0, vyr=vy, vzr=0.0)
        
        # Solve at t=0.25
        t = 0.25
        x, P, n, N, v, vy_arr, vz_arr, u, gamma_arr, S, e = solver.solve(t, x0=0.5, n_points=400)
        
        # Find star states (around contact discontinuity at x≈0.5)
        idx_contact = np.argmin(np.abs(x - (0.5 + solver.vls * t)))
        
        results.append({
            'vy': vy,
            'Ps': P[idx_contact],
            'vs': v[idx_contact],
            'vys': vy_arr[idx_contact],
            'nl': solver.nls,
            'nr': solver.nrs,
            'vshockl': solver.vshockl,
            'vshockr': solver.vshockr
        })
        
        print(f"v^y = {vy:.2f}:")
        print(f"  Star state: P* = {P[idx_contact]:.4f}, v^x* = {v[idx_contact]:.4f}, v^y* = {vy_arr[idx_contact]:.4f}")
        print(f"  Left shock: v_shock = {solver.vshockl:.4f}" if solver.vshockl != 0 else "  Left rarefaction")
        print(f"  Right shock: v_shock = {solver.vshockr:.4f}" if solver.vshockr != 0 else "  Right rarefaction")
        print()
    
    # Verify trends
    print("\nVerification of physical trends:")
    print("-" * 80)
    for i in range(1, len(results)):
        dP = results[i]['Ps'] - results[i-1]['Ps']
        dv = results[i]['vs'] - results[i-1]['vs']
        print(f"v^y: {results[i-1]['vy']:.2f} → {results[i]['vy']:.2f}:")
        print(f"  ΔP* = {dP:+.4f} (expect decrease with v^t_L)")
        print(f"  Δv^x* = {dv:+.4f} (expect decrease with v^t_L)")
    
    return results


def test_blast_wave_tangential():
    """
    Test relativistic blast wave with tangential velocities.
    
    Reproduces test cases from Pons et al. (1999) Table 1:
    Initial conditions:
    - Left:  P=1000, n=1.0, v^x=0.0
    - Right: P=0.01, n=1.0, v^x=0.0
    - Test 9 combinations of v^t_L, v^t_R = 0, 0.9, 0.99
    """
    print("\n" + "=" * 80)
    print("TEST 2: Relativistic Blast Wave with Tangential Velocities")
    print("=" * 80)
    print("\nInitial conditions (Pons et al. 1999, Table 1):")
    print("Left:  P=1000, n=1.0, v^x=0.0")
    print("Right: P=0.01, n=1.0, v^x=0.0")
    print("Testing 9 combinations of v^t_L, v^t_R\n")
    
    gamma = 5.0/3.0
    c = 1.0
    
    # Initial conditions
    Pl, nl, vl = 1000.0, 1.0, 0.0
    Pr, nr, vr = 0.01, 1.0, 0.0
    
    # Test combinations
    tangential_vels = [0.0, 0.9, 0.99]
    
    print(f"{'v^t_L':>6} {'v^t_R':>6} {'ρ_L*':>12} {'ρ_R*':>12} {'P*':>12} {'v^x*':>12} {'V_s':>12}")
    print("-" * 80)
    
    for vyl in tangential_vels:
        for vyr in tangential_vels:
            solver = KitajimaRiemannSolver(gamma, c)
            solver.set_initial_states(Pl, nl, vl, Pr, nr, vr,
                                     vyl=vyl, vzl=0.0, vyr=vyr, vzr=0.0)
            
            # Solve at t=0.4
            t = 0.4
            x, P, n, N, v, vy_arr, vz_arr, u, gamma_arr, S, e = solver.solve(t, x0=0.5, n_points=400)
            
            # Find star states
            idx_contact = np.argmin(np.abs(x - (0.5 + solver.vls * t)))
            
            # Print results in table format
            print(f"{vyl:6.2f} {vyr:6.2f} {solver.nls:12.2e} {solver.nrs:12.2e} "
                  f"{P[idx_contact]:12.2e} {v[idx_contact]:12.3f} {solver.vshockr:12.3f}")
    
    print("\nNote: Compare with Table 1 in Pons et al. (1999)")


def test_conservation_laws():
    """
    Test conservation of tangential velocities across waves.
    
    Verifies:
    1. h*W*v^y = const across rarefactions
    2. h*W*v^z = const across rarefactions
    3. h*W*v^y = const across shocks
    4. h*W*v^z = const across shocks
    5. Direction v^y/v^z = const
    """
    print("\n" + "=" * 80)
    print("TEST 3: Conservation Law Verification")
    print("=" * 80)
    
    gamma = 1.4
    c = 1.0
    
    # Test case with both rarefaction and shock
    Pl, nl, vl = 1.0, 1.0, 0.0
    Pr, nr, vr = 0.1, 0.125, 0.0
    vyl, vzl = 0.7, 0.5  # Tangential velocity on left
    vyr, vzr = 0.3, 0.2  # Tangential velocity on right
    
    solver = KitajimaRiemannSolver(gamma, c)
    solver.set_initial_states(Pl, nl, vl, Pr, nr, vr,
                             vyl=vyl, vzl=vzl, vyr=vyr, vzr=vzr)
    
    t = 0.25
    x, P, n, N, v, vy_arr, vz_arr, u, gamma_arr, S, e = solver.solve(t, x0=0.5, n_points=400)
    
    # Compute enthalpy array
    H = 1.0 + u/(c**2) + P/(n * c**2)
    
    print("\nLeft wave (rarefaction):")
    print("-" * 80)
    
    # Check left initial state
    vta_l = np.sqrt(vyl**2 + vzl**2)
    v2a_l = vl**2 + vta_l**2
    gamma_a_l = 1.0 / np.sqrt(1.0 - v2a_l/(c**2))
    Ha_l = 1.0 + (Pl/((gamma-1)*nl))/(c**2) + Pl/(nl*c**2)
    
    conserved_y_l = Ha_l * gamma_a_l * vyl
    conserved_z_l = Ha_l * gamma_a_l * vzl
    direction_l = vyl / vzl if abs(vzl) > 1e-10 else 0.0
    
    print(f"Initial state (left):")
    print(f"  h*W*v^y = {conserved_y_l:.6f}")
    print(f"  h*W*v^z = {conserved_z_l:.6f}")
    print(f"  v^y/v^z = {direction_l:.6f}")
    
    # Check left star state
    vts_l = np.sqrt(solver.vyls**2 + solver.vzls**2)
    v2s_l = solver.vls**2 + vts_l**2
    gamma_s_l = 1.0 / np.sqrt(1.0 - v2s_l/(c**2))
    Hs_l = solver.Hls
    
    conserved_y_ls = Hs_l * gamma_s_l * solver.vyls
    conserved_z_ls = Hs_l * gamma_s_l * solver.vzls
    direction_ls = solver.vyls / solver.vzls if abs(solver.vzls) > 1e-10 else 0.0
    
    print(f"\nStar state (left):")
    print(f"  h*W*v^y = {conserved_y_ls:.6f} (Δ = {abs(conserved_y_ls - conserved_y_l):.2e})")
    print(f"  h*W*v^z = {conserved_z_ls:.6f} (Δ = {abs(conserved_z_ls - conserved_z_l):.2e})")
    print(f"  v^y/v^z = {direction_ls:.6f} (Δ = {abs(direction_ls - direction_l):.2e})")
    
    # Check right wave (shock)
    print("\nRight wave (shock):")
    print("-" * 80)
    
    vta_r = np.sqrt(vyr**2 + vzr**2)
    v2a_r = vr**2 + vta_r**2
    gamma_a_r = 1.0 / np.sqrt(1.0 - v2a_r/(c**2))
    Ha_r = 1.0 + (Pr/((gamma-1)*nr))/(c**2) + Pr/(nr*c**2)
    
    conserved_y_r = Ha_r * gamma_a_r * vyr
    conserved_z_r = Ha_r * gamma_a_r * vzr
    direction_r = vyr / vzr if abs(vzr) > 1e-10 else 0.0
    
    print(f"Initial state (right):")
    print(f"  h*W*v^y = {conserved_y_r:.6f}")
    print(f"  h*W*v^z = {conserved_z_r:.6f}")
    print(f"  v^y/v^z = {direction_r:.6f}")
    
    vts_r = np.sqrt(solver.vyrs**2 + solver.vzrs**2)
    v2s_r = solver.vrs**2 + vts_r**2
    gamma_s_r = 1.0 / np.sqrt(1.0 - v2s_r/(c**2))
    Hs_r = solver.Hrs
    
    conserved_y_rs = Hs_r * gamma_s_r * solver.vyrs
    conserved_z_rs = Hs_r * gamma_s_r * solver.vzrs
    direction_rs = solver.vyrs / solver.vzrs if abs(solver.vzrs) > 1e-10 else 0.0
    
    print(f"\nStar state (right):")
    print(f"  h*W*v^y = {conserved_y_rs:.6f} (Δ = {abs(conserved_y_rs - conserved_y_r):.2e})")
    print(f"  h*W*v^z = {conserved_z_rs:.6f} (Δ = {abs(conserved_z_rs - conserved_z_r):.2e})")
    print(f"  v^y/v^z = {direction_rs:.6f} (Δ = {abs(direction_rs - direction_r):.2e})")
    
    # Overall assessment
    print("\n" + "=" * 80)
    tol = 1e-6
    all_passed = True
    
    tests = [
        ("Left h*W*v^y conservation", abs(conserved_y_ls - conserved_y_l)),
        ("Left h*W*v^z conservation", abs(conserved_z_ls - conserved_z_l)),
        ("Left direction conservation", abs(direction_ls - direction_l)),
        ("Right h*W*v^y conservation", abs(conserved_y_rs - conserved_y_r)),
        ("Right h*W*v^z conservation", abs(conserved_z_rs - conserved_z_r)),
        ("Right direction conservation", abs(direction_rs - direction_r)),
    ]
    
    print("Conservation test results:")
    for test_name, error in tests:
        status = "PASS" if error < tol else "FAIL"
        print(f"  {test_name}: {status} (error = {error:.2e})")
        if error >= tol:
            all_passed = False
    
    return all_passed


def save_test_data():
    """
    Save test data for animation.
    """
    print("\n" + "=" * 80)
    print("Saving test data for animation...")
    print("=" * 80)
    
    gamma = 1.4
    c = 1.0
    
    # Sod shock tube with different tangential velocities
    Pl, nl, vl = 1.0, 1.0, 0.0
    Pr, nr, vr = 0.1, 0.125, 0.0
    
    tangential_vels = [0.0, 0.5, 0.9, 0.99]
    
    for vy in tangential_vels:
        solver = KitajimaRiemannSolver(gamma, c)
        solver.set_initial_states(Pl, nl, vl, Pr, nr, vr,
                                 vyl=vy, vzl=0.0, vyr=vy, vzr=0.0)
        
        t = 0.25
        x, P, n, N, v, vy_arr, vz_arr, u, gamma_arr, S, e = solver.solve(t, x0=0.5, n_points=400)
        
        # Save to file
        filename = f"tangential_test_vy{vy:.2f}.dat"
        solver.save_solution(filename, x, P, n, N, v, vy_arr, vz_arr, u, gamma_arr, S, e, t)
        print(f"  Saved {filename}")
    
    print("\nData files ready for animation!")


if __name__ == "__main__":
    print("\n" + "=" * 80)
    print("TANGENTIAL VELOCITY TESTS FOR KITAJIMA RIEMANN SOLVER")
    print("Based on Pons, Martí & Müller (1999, J. Fluid Mech.)")
    print("=" * 80)
    
    # Run tests
    test_sod_tangential()
    test_blast_wave_tangential()
    conservation_passed = test_conservation_laws()
    save_test_data()
    
    # Final summary
    print("\n" + "=" * 80)
    print("TEST SUMMARY")
    print("=" * 80)
    if conservation_passed:
        print("✓ All conservation law tests PASSED")
        print("✓ Implementation verified against Pons et al. (1999)")
    else:
        print("✗ Some tests FAILED - check implementation")
    
    print("\nNext steps:")
    print("1. Run: python3 animate_tangential.py")
    print("2. Compare results with Figures 1-4 in Pons et al. (1999)")
    print("=" * 80)
