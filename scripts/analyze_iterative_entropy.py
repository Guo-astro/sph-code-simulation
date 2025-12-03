#!/usr/bin/env python3
"""Analyze entropy for iterative solver hydrostatic test."""

import pandas as pd
import numpy as np

# Load snapshots
s0 = pd.read_csv('sample/imbh_cloud/results/hydrostatic/10k/GSPH/snapshot_0000.csv')
s1 = pd.read_csv('sample/imbh_cloud/results/hydrostatic/10k/GSPH/snapshot_0001.csv')

# Calculate radius
s0['r'] = np.sqrt(s0['x']**2 + s0['y']**2 + s0['z']**2)
s1['r'] = np.sqrt(s1['x']**2 + s1['y']**2 + s1['z']**2)

# Calculate entropy proxy s = P/rho^gamma (gamma = 5/3)
gamma = 5/3
s0['entropy'] = s0['pres'] / (s0['dens']**(gamma))
s1['entropy'] = s1['pres'] / (s1['dens']**(gamma))

# Get central particles (r < 0.1)
c0 = s0[s0['r'] < 0.1]
c1 = s1[s1['r'] < 0.1]

print('=== ITERATIVE SOLVER - Early Entropy Check ===')
print(f't=0: Central particles = {len(c0)}')
print(f'  Max density: {c0["dens"].max():.6f}')
print(f'  Mean entropy: {c0["entropy"].mean():.6f}')
print(f'  Min entropy: {c0["entropy"].min():.6f}')
print()
print(f't=1: Central particles = {len(c1)}')  
print(f'  Max density: {c1["dens"].max():.6f}')
print(f'  Mean entropy: {c1["entropy"].mean():.6f}')
print(f'  Min entropy: {c1["entropy"].min():.6f}')
print()
entropy_change = (c1["entropy"].mean() - c0["entropy"].mean()) / c0["entropy"].mean() * 100
density_change = (c1["dens"].max() - c0["dens"].max()) / c0["dens"].max() * 100
print(f'Entropy change: {entropy_change:.3f}%')
print(f'Density change: {density_change:.3f}%')

# Compare with HLL at t=81 (for reference)
print()
print('=== For Reference: HLL at t=81 ===')
print('  Entropy change at center: -11% (from 0.408 to 0.363)')
print('  Density change at center: +8% (from 1.457 to 1.572)')

# Estimate: if entropy loss is linear over time
if entropy_change != 0:
    estimated_t81_entropy_change = entropy_change * 81
    print()
    print(f'=== Estimate for t=81 (linear extrapolation) ===')
    print(f'  Estimated entropy change: {estimated_t81_entropy_change:.3f}%')
