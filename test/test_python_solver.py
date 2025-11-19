#!/usr/bin/env python3
"""Quick test of Python solver"""
import sys
sys.path.insert(0, '/Users/guo/Downloads/sphcode/relativistic-riemann')
from kitajima_solver import KitajimaRiemannSolver

# Test SR Sod
solver = KitajimaRiemannSolver(gamma_c=5.0/3.0, c=1.0)
solver.set_initial_states(
    Pl=13.33, nl=10.0, vl=0.0, vyl=0.0, vzl=0.0,
    Pr=1.0e-8, nr=1.0, vr=0.0, vyr=0.0, vzr=0.0
)

# Solve at t=0.4
x, P, n, N, v, vy, vz, u, gamma, S, e = solver.solve(t=0.4, tol=1.0e-10)

# Find contact discontinuity (where v is constant = vstar)
import numpy as np
v_unique = np.unique(np.round(v, 10))
print(f"Unique velocities: {v_unique}")

# Star region should have constant v and P
mid = len(x) // 2
print(f"Mid-region: P[{mid}] = {P[mid]:.16e}, v[{mid}] = {v[mid]:.16e}")

# Now test get_pressure directly
pmin = (solver.Pl + solver.Pr) / 2.0
pmax = pmin

for _ in range(100):
    pmin = 0.5 * max(pmin, 0.0)
    pmax = 2.0 * pmax
    
    dvel1 = solver.get_dvel(pmin)
    dvel2 = solver.get_dvel(pmax)
    
    if dvel1 * dvel2 <= 0.0:
        break

print(f"\nBracketing: pmin={pmin:.6e}, pmax={pmax:.6e}")
print(f"f(pmin)={dvel1:.6e}, f(pmax)={dvel2:.6e}")

Ps_direct = solver.get_pressure(pmin, pmax, tol=1.0e-10)
vs_direct = 0.5 * (solver.vls + solver.vrs)

print(f"\nDirect get_pressure(): Ps = {Ps_direct:.16e}")
print(f"Direct velocities: vs = {vs_direct:.16e}")
