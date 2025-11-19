#!/usr/bin/env python3
"""
Run Python Riemann solver on same test cases and compare with C++
"""
import sys
sys.path.insert(0, '/Users/guo/Downloads/sphcode/relativistic-riemann')
from kitajima_solver import KitajimaRiemannSolver

test_cases = [
    {
        'name': 'SR_Sod_Shock_Tube',
        'gamma': 5.0/3.0,
        'c': 1.0,
        'left': {'v': 0.0, 'n': 10.0, 'P': 13.33, 'vt': 0.0},
        'right': {'v': 0.0, 'n': 1.0, 'P': 1.0e-8, 'vt': 0.0},
        'cpp': {'pstar': 1.447682689862212, 'vstar': 0.7139906460205292}
    },
    {
        'name': 'Two_Rarefactions',
        'gamma': 5.0/3.0,
        'c': 1.0,
        'left': {'v': -0.5, 'n': 1.0, 'P': 1.0, 'vt': 0.0},
        'right': {'v': 0.5, 'n': 1.0, 'P': 1.0, 'vt': 0.0},
        'cpp': {'pstar': 0.2497071955739552, 'vstar': -5.551115123141402e-17}
    },
    {
        'name': 'Two_Shocks',
        'gamma': 5.0/3.0,
        'c': 1.0,
        'left': {'v': 0.5, 'n': 1.0, 'P': 1.0, 'vt': 0.0},
        'right': {'v': -0.5, 'n': 1.0, 'P': 1.0, 'vt': 0.0},
        'cpp': {'pstar': 3.591598452917558, 'vstar': 0.0}
    },
    {
        'name': 'Left_Shock_Right_Rarefaction',
        'gamma': 5.0/3.0,
        'c': 1.0,
        'left': {'v': 0.3, 'n': 1.0, 'P': 10.0, 'vt': 0.0},
        'right': {'v': 0.0, 'n': 1.0, 'P': 1.0, 'vt': 0.0},
        'cpp': {'pstar': 4.705161581292606, 'vstar': 0.5853329022543065}
    },
    {
        'name': 'With_Tangential_Velocity',
        'gamma': 5.0/3.0,
        'c': 1.0,
        'left': {'v': 0.0, 'n': 1.0, 'P': 1.0, 'vt': 0.3},
        'right': {'v': 0.0, 'n': 1.0, 'P': 0.1, 'vt': 0.0},
        'cpp': {'pstar': 0.4156788143096682, 'vstar': 0.3396854041208034}
    }
]

print("=" * 80)
print("EXACT MATCHING TEST: C++ vs Python Brent's Method")
print("=" * 80)

all_match = True
max_p_err = 0.0
max_v_err = 0.0

for test in test_cases:
    # Create solver
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
    
    # Bracketing (matching Python code exactly)
    pmin = (solver.Pl + solver.Pr) / 2.0
    pmax = pmin
    
    for i in range(100):
        pmin = 0.5 * max(pmin, 0.0)
        pmax = 2.0 * pmax
        
        dvel1 = solver.get_dvel(pmin)
        dvel2 = solver.get_dvel(pmax)
        
        print(f"  Bracket iter {i}: pmin={pmin:.3e}, pmax={pmax:.3e}, f(pmin)={dvel1:.3e}, f(pmax)={dvel2:.3e}, product={dvel1*dvel2:.3e}")
        
        if dvel1 * dvel2 <= 0.0:
            print(f"  ✓ Bracketing complete: iterations={i+1}")
            break
    
    # Brent's method
    print(f"  Calling get_pressure with pmin={pmin:.16e}, pmax={pmax:.16e}")
    
    # Check initial function values
    fa_init = solver.get_dvel(pmin)
    fb_init = solver.get_dvel(pmax)
    print(f"  Initial: f(pmin)={fa_init:.16e}, f(pmax)={fb_init:.16e}")
    
    Ps = solver.get_pressure(pmin, pmax, tol=1.0e-10)
    print(f"  get_pressure returned: Ps={Ps:.16e}")
    
    # CRITICAL: get_pressure() leaves vls/vrs at the last evaluation point (not Ps!)
    # Must re-evaluate at Ps to get correct velocities
    _ = solver.get_dvel(Ps)
    vs = 0.5 * (solver.vls + solver.vrs)
    print(f"  After Brent: Ps={Ps:.16e}, vls={solver.vls:.16e}, vrs={solver.vrs:.16e}, vs={vs:.16e}")
    
    # Compare with C++
    Ps_cpp = test['cpp']['pstar']
    vs_cpp = test['cpp']['vstar']
    
    dp = abs(Ps - Ps_cpp)
    dv = abs(vs - vs_cpp)
    
    p_rel = dp / abs(Ps) if Ps != 0 else 0.0
    v_rel = dv / abs(vs) if vs != 0 else dv
    
    max_p_err = max(max_p_err, p_rel)
    max_v_err = max(max_v_err, v_rel)
    
    print(f"\n{test['name']}:")
    print(f"  Python: P* = {Ps:.16e}, v* = {vs:.16e}")
    print(f"  C++:    P* = {Ps_cpp:.16e}, v* = {vs_cpp:.16e}")
    print(f"  ΔP = {dp:.3e} (rel: {p_rel:.3e})")
    print(f"  Δv = {dv:.3e} (rel: {v_rel:.3e})")
    
    # Machine precision check
    match = (p_rel < 1e-14) and (v_rel < 1e-14 or (vs == 0 and vs_cpp == 0))
    
    if match:
        print(f"  ✅ EXACT MATCH")
    else:
        print(f"  ❌ MISMATCH")
        all_match = False

print("\n" + "=" * 80)
print(f"Maximum relative errors:")
print(f"  Pressure:  {max_p_err:.3e}")
print(f"  Velocity:  {max_v_err:.3e}")
print("=" * 80)

if all_match:
    print("✅ SUCCESS: C++ matches Python EXACTLY (within machine precision)")
    print("=" * 80)
    sys.exit(0)
else:
    print("❌ FAILURE: Differences exceed machine precision")
    print("=" * 80)
    sys.exit(1)
