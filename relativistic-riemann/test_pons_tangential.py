#!/usr/bin/env python3
"""
Test script for the updated Kitajima solver with tangential velocities.

Verifies against Table 1 from Pons, Martí & Müller (1999):
"The exact solution of the Riemann problem with non-zero tangential velocities
in relativistic hydrodynamics"
"""

import numpy as np
from kitajima_solver import KitajimaRiemannSolver


def test_pons_blast_wave():
    """
    Test relativistic blast wave problem from Pons et al. (1999), Table 1.
    
    Initial data:
    - Left:  P_L = 1000, rho_L = 1, v^x_L = 0
    - Right: P_R = 0.01, rho_R = 1, v^x_R = 0
    - gamma = 5/3
    - Time t = 0.4
    
    Test 9 combinations of tangential velocities: v^t_{L,R} = 0, 0.9, 0.99
    """
    
    print("=" * 80)
    print("Testing Pons et al. (1999) Relativistic Blast Wave")
    print("=" * 80)
    
    # Initial conditions (in natural units c=1)
    gamma_c = 5.0/3.0
    c = 1.0
    
    Pl = 1000.0
    nl = 1.0
    vxl = 0.0
    
    Pr = 0.01
    nr = 1.0
    vxr = 0.0
    
    t = 0.4
    
    # Expected results from Table 1
    expected = {
        (0.00, 0.00): {'rho_ls': 9.16e-2, 'rho_rs': 1.04e1, 'p_s': 1.86e1, 'v_s': 0.960},
        (0.00, 0.90): {'rho_ls': 1.51e-1, 'rho_rs': 1.46e1, 'p_s': 4.28e1, 'v_s': 0.913},
        (0.00, 0.99): {'rho_ls': 2.89e-1, 'rho_rs': 4.36e1, 'p_s': 1.27e2, 'v_s': 0.767},
        (0.90, 0.00): {'rho_ls': 5.83e-3, 'rho_rs': 3.44e0, 'p_s': 1.89e-1, 'v_s': 0.328},
        (0.90, 0.90): {'rho_ls': 1.49e-2, 'rho_rs': 4.46e0, 'p_s': 9.04e-1, 'v_s': 0.319},
        (0.90, 0.99): {'rho_ls': 5.72e-2, 'rho_rs': 7.83e0, 'p_s': 8.48e0, 'v_s': 0.292},
        (0.99, 0.00): {'rho_ls': 1.99e-3, 'rho_rs': 1.91e0, 'p_s': 3.16e-2, 'v_s': 0.099},
        (0.99, 0.90): {'rho_ls': 3.80e-3, 'rho_rs': 2.90e0, 'p_s': 9.27e-2, 'v_s': 0.098},
        (0.99, 0.99): {'rho_ls': 1.29e-2, 'rho_rs': 4.29e0, 'p_s': 7.06e-1, 'v_s': 0.095},
    }
    
    print("\nTesting 9 combinations of tangential velocities:")
    print(f"{'v^t_L':>6} {'v^t_R':>6} | {'rho_L*':>10} {'rho_R*':>10} {'P*':>10} {'v^x*':>8} | Status")
    print("-" * 80)
    
    all_passed = True
    
    for vtl in [0.0, 0.9, 0.99]:
        for vtr in [0.0, 0.9, 0.99]:
            # Create solver
            solver = KitajimaRiemannSolver(gamma_c, c)
            
            # Set initial states (tangential velocity in y-direction)
            solver.set_initial_states(Pl, nl, vxl, Pr, nr, vxr, 
                                     vyl=vtl, vzl=0.0, vyr=vtr, vzr=0.0)
            
            # Solve
            try:
                x, P, n, N, v, vy, vz, u, gamma, S, e = solver.solve(t)
                
                # Find intermediate states (at contact discontinuity x ~ 0.5)
                idx_contact = np.argmin(np.abs(x - 0.5))
                
                # Get values from solution
                rho_ls = solver.nls  # Left star density
                rho_rs = solver.nrs  # Right star density
                p_s = P[idx_contact]
                v_s = v[idx_contact]
                
                # Compare with expected
                exp = expected.get((vtl, vtr))
                if exp is not None:
                    err_rho_ls = abs(rho_ls - exp['rho_ls']) / exp['rho_ls']
                    err_rho_rs = abs(rho_rs - exp['rho_rs']) / exp['rho_rs']
                    err_p_s = abs(p_s - exp['p_s']) / exp['p_s']
                    err_v_s = abs(v_s - exp['v_s']) / max(exp['v_s'], 0.01)
                    
                    # Tolerance (10% for this initial test)
                    tol = 0.10
                    passed = (err_rho_ls < tol and err_rho_rs < tol and 
                             err_p_s < tol and err_v_s < tol)
                    
                    status = "✓ PASS" if passed else "✗ FAIL"
                    if not passed:
                        all_passed = False
                else:
                    status = "? N/A"
                    passed = True
                
                print(f"{vtl:6.2f} {vtr:6.2f} | {rho_ls:10.3e} {rho_rs:10.3e} "
                      f"{p_s:10.3e} {v_s:8.3f} | {status}")
                
                if exp and not passed:
                    print(f"  Expected: {exp['rho_ls']:10.3e} {exp['rho_rs']:10.3e} "
                          f"{exp['p_s']:10.3e} {exp['v_s']:8.3f}")
                    print(f"  Errors:   {err_rho_ls:10.1%} {err_rho_rs:10.1%} "
                          f"{err_p_s:10.1%} {err_v_s:10.1%}")
                
            except Exception as ex:
                print(f"{vtl:6.2f} {vtr:6.2f} | {'ERROR':>10} | {str(ex)[:40]}")
                all_passed = False
    
    print("-" * 80)
    if all_passed:
        print("✓ All tests PASSED")
    else:
        print("✗ Some tests FAILED")
    print()
    
    return all_passed


def test_simple_case():
    """Test simple case without tangential velocities (Sod problem)."""
    
    print("=" * 80)
    print("Testing Simple Sod Problem (no tangential velocities)")
    print("=" * 80)
    
    gamma_c = 1.4
    c = 1.0
    
    # Sod shock tube
    solver = KitajimaRiemannSolver(gamma_c, c)
    solver.set_initial_states(Pl=1.0, nl=1.0, vl=0.0, 
                             Pr=0.1, nr=0.125, vr=0.0)
    
    try:
        x, P, n, N, v, vy, vz, u, gamma, S, e = solver.solve(t=0.25)
        
        print(f"Solution computed successfully")
        print(f"  Pressure range: [{P.min():.3f}, {P.max():.3f}]")
        print(f"  Density range:  [{n.min():.3f}, {n.max():.3f}]")
        print(f"  Velocity range: [{v.min():.3f}, {v.max():.3f}]")
        print(f"  Star pressure:  {P[len(P)//2]:.3f}")
        print(f"  Star velocity:  {v[len(v)//2]:.3f}")
        
        # Check tangential velocities are zero
        assert np.allclose(vy, 0.0), "Tangential vy should be zero"
        assert np.allclose(vz, 0.0), "Tangential vz should be zero"
        
        print("✓ Test PASSED")
        return True
        
    except Exception as ex:
        print(f"✗ Test FAILED: {ex}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    # Run tests
    test1_passed = test_simple_case()
    print()
    test2_passed = test_pons_blast_wave()
    
    # Summary
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Simple Sod test:         {'✓ PASSED' if test1_passed else '✗ FAILED'}")
    print(f"Pons blast wave test:    {'✓ PASSED' if test2_passed else '✗ FAILED'}")
    print()
