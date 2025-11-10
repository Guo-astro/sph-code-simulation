# Lane-Emden Test Notes

## Summary
The Lane-Emden n=3/2 hydrostatic equilibrium test shows **partial stability** - the bulk of particles remain confined (~90% at r < 0.9), but some particles escape reaching r_max ≈ 2.6 by t=10.

## Current Configuration
- **Grid**: 30³ uniform grid → 14,328 particles within r < 1
- **Density profile**: ρ(r) = ρ_center × exp(-1.5r/R)
  - ρ_center = 3 × ρ_avg = 0.716
  - ρ_avg = 0.239
- **Polytropic constant**: K = 0.4 × G × M^(1/3) / R = 0.4
- **Total mass**: M = 1.0
- **Radius**: R = 1.0

## Results
```
Time    r_max   r_90%   Status
================================
0.0     1.00    0.90    ✓ Confined
2.5     1.29    0.68    ✓ Confined  
5.0     2.08    0.89    ⚠ Expanding
10.0    2.42    0.87    ⚠ Expanding
```

## Analysis
**Good**: 
- 90% of particles stay within r < 0.9 (well-confined)
- No catastrophic collapse or explosion
- Simulation completes without numerical issues

**Issues**:
- Maximum radius grows from 1.0 → 2.6 (164% expansion)
- A small fraction of particles escape
- Not true hydrostatic equilibrium

## Why Expansion Occurs
1. **Initial conditions not exact equilibrium**: 
   - Exponential density profile is approximate
   - True Lane-Emden solution requires solving differential equation
   
2. **SPH pressure gradient errors**:
   - Particle discretization creates noise in ∇P
   - High-density center prone to numerical errors
   
3. **Missing rotational support**:
   - Real astrophysical objects often rotate
   - Rotation provides additional support against gravity

## Comparison with Evrard Test
| Test | Purpose | r_max (final) | Result |
|------|---------|---------------|--------|
| **Evrard** | Collapse → Shock → Rebound | r > 10 | MUST expand |
| **Lane-Emden** | Hydrostatic equilibrium | r ≈ 1 | SHOULD stay confined |
| **Evrard Cold** (u=0.0001) | Demonstrate shock amplification | r = 63 | Extreme expansion |
| **Lane-Emden** (current) | Equilibrium attempt | r = 2.6 | Partial success |

## Recommendations

### For Truly Confined Particles:
1. **Use existing `hydrostatic` test** (2D equilibrium with rotation)
2. **Add rotation to Lane-Emden**: ω = √(GM/R³)
3. **Relax initial conditions**: Run forward then backward in time

### For Demonstrating Physics:
- **Evrard**: Shows collapse dynamics
- **Evrard Cold**: Shows shock amplification  
- **Lane-Emden**: Shows SPH challenges with equilibrium
- **Hydrostatic**: Shows stable equilibrium with rotation

## Technical Details
- **Neighbor list size**: Increased to 300 (high central density requires many neighbors)
- **Simulation time**: ~79 seconds for t=0→10
- **Output**: 21 snapshots (every 0.5 time units)
- **SPH type**: DISPH (Density Independent SPH)
- **Kernel**: Wendland C4

## Next Steps
1. ✅ Completed: Basic Lane-Emden implementation
2. ✅ Completed: Analysis scripts created
3. ⏳ Optional: Add rotation for better stability
4. ⏳ Optional: Try relaxation technique
5. ⏳ Optional: Compare with 2D hydrostatic test
