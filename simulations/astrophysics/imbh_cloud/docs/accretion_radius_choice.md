# Choosing the Accretion Radius for IMBH-Cloud Simulations

## Unit System

| Quantity | Code Unit | Physical Unit |
|----------|-----------|---------------|
| Length   | 1 pc      | 3.086 × 10¹⁸ cm |
| Mass     | 1 M☉     | 1.989 × 10³³ g |
| Velocity | 1 km/s    | 10⁵ cm/s |
| Time     | 1 code    | 0.978 Myr |
| G        | 0.00430091 | 6.674 × 10⁻⁸ cgs |

## Physical Setup

- **Black Hole Mass**: M_BH = 10⁵ M☉ = 1.989 × 10³⁸ g
- **Cloud Mass**: M_cloud = 127.5 M☉ = 2.54 × 10³⁵ g
- **Cloud Radius**: R_cloud = 1.55 pc = 4.78 × 10¹⁸ cm
- **Cloud Temperature**: T = 35 K
- **Central Density**: n_center = 500 cm⁻³
- **Pericenter Distance**: r_p = 0.4 pc = 1.23 × 10¹⁸ cm
- **Pericenter Velocity**: v_p = 46.4 km/s = 4.64 × 10⁶ cm/s

## The Timestep Problem

### CFL Constraints in SPH

The simulation timestep is limited by multiple CFL conditions:

1. **Sound CFL** (wave propagation):
   ```
   dt_sound = C_sound × h / c_s
   ```

   For isothermal gas at T = 35 K:
   ```
   c_s = √(k_B T / μ m_H) = 0.48 km/s = 4.8 × 10⁴ cm/s
   ```

2. **Force CFL** (acceleration limit):
   ```
   dt_force = C_force × √(h / |a|)
   ```

### Why Force CFL Dominates Near the BH

The gravitational acceleration from the BH at distance r is:
```
a_BH = G × M_BH / r²
```

| Distance r | r (cm) | a_BH (cm/s²) | a_BH (code) | dt_force |
|------------|--------|--------------|-------------|----------|
| 1.0 pc | 3.1 × 10¹⁸ | 4.4 × 10⁻³ | 4.3 × 10³ | 0.15 Myr |
| 0.1 pc | 3.1 × 10¹⁷ | 4.4 × 10⁻¹ | 4.3 × 10⁵ | 15 kyr |
| 0.05 pc | 1.5 × 10¹⁷ | 1.8 | 1.7 × 10⁶ | 7.3 kyr |
| 0.02 pc | 6.2 × 10¹⁶ | 11 | 1.1 × 10⁷ | 2.9 kyr |
| 0.01 pc | 3.1 × 10¹⁶ | 44 | 4.3 × 10⁷ | 1.5 kyr |

**Assumptions**: h = 0.01 pc = 3.1 × 10¹⁶ cm, C_force = 0.3

**Key insight**: Timestep scales as dt ∝ r, so particles at 0.01 pc require 10× smaller timesteps than at 0.1 pc.

## Relevant Physical Scales

### 1. Schwarzschild Radius
```
r_s = 2 G M_BH / c²
    = 2 × (6.67 × 10⁻⁸) × (10⁵ × 2 × 10³³) / (3 × 10¹⁰)²
    = 2.95 × 10¹⁰ cm = 9.6 × 10⁻⁹ pc
```
Far below our resolution - we use a sink particle instead.

### 2. Tidal Radius
The distance where BH tidal forces overcome cloud self-gravity:
```
r_tidal = R_cloud × (M_BH / M_cloud)^(1/3)
        = 1.55 pc × (10⁵ / 127.5)^(1/3)
        = 1.55 pc × 9.3
        = 14.4 pc = 4.4 × 10¹⁹ cm
```
Cloud disruption begins well outside any reasonable sink radius.

### 3. Bondi Radius
The radius where BH gravity dominates over thermal pressure:
```
r_B = G M_BH / c_s²
    = (6.67 × 10⁻⁸) × (10⁵ × 2 × 10³³) / (4.8 × 10⁴)²
    = 5.8 × 10¹⁸ cm = 1.9 pc
```
Relevant for steady-state accretion, less so for tidal disruption.

### 4. Smoothing Length
Typical SPH resolution scale:
```
h ≈ 0.01 - 0.02 pc = (3 - 6) × 10¹⁶ cm
```
Sink radius should be > h for proper resolution.

## Accretion Radius Selection

### Physical Constraints

| Constraint | Requirement | Value |
|------------|-------------|-------|
| Inside orbit | r_acc < r_p | < 0.4 pc |
| Resolved by SPH | r_acc > few × h | > 0.02 pc |
| Above softening | r_acc > ε | > 0.005 pc |
| Capture dynamics | r_acc << r_tidal | << 14.4 pc |

### Computational Constraints

Simulation time to reach t_end = 6 code units = 5.87 Myr:

| r_acc | dt_min | Physical dt | Loops needed | Wall time |
|-------|--------|-------------|--------------|-----------|
| 0.01 pc | 3 × 10⁻⁵ | 29 yr | ~200,000 | Days |
| 0.05 pc | 1.5 × 10⁻⁴ | 150 yr | ~40,000 | Hours |
| 0.1 pc | 3 × 10⁻⁴ | 290 yr | ~20,000 | Hours |
| 0.2 pc | 6 × 10⁻⁴ | 590 yr | ~10,000 | Hour |

### Our Choice: r_acc = 0.1 pc = 3.1 × 10¹⁷ cm

**Rationale:**

1. **Resolves tidal disruption**: r_tidal = 14.4 pc >> r_acc = 0.1 pc
   - All tidal physics captured before accretion

2. **Captures pericenter dynamics**: r_acc = r_p / 4
   - Stream formation visible
   - Partial accretion measured

3. **SPH-resolved**: r_acc = 10h (with h ~ 0.01 pc)
   - Sink boundary well-resolved
   - No numerical artifacts

4. **Computationally feasible**: dt ~ 290 yr per step
   - ~20,000 loops to t = 5.87 Myr
   - Tractable on single workstation

## Physical Interpretation

At r = 0.1 pc from a 10⁵ M☉ BH:
- **Orbital velocity**: v_orb = √(GM/r) = 65 km/s
- **Dynamical time**: t_dyn = √(r³/GM) = 4700 yr
- **Free-fall time**: t_ff = π/2 × t_dyn = 7400 yr

Gas reaching r_acc = 0.1 pc is essentially committed to accretion on timescales << simulation time.

## Summary

| Property | Value | Physical |
|----------|-------|----------|
| Sink radius | 0.1 pc | 3.1 × 10¹⁷ cm |
| Min timestep | ~3 × 10⁻⁴ code | ~290 yr |
| Tidal radius | 14.4 pc | 4.4 × 10¹⁹ cm |
| Pericenter | 0.4 pc | 1.2 × 10¹⁸ cm |
| r_acc / r_p | 0.25 | - |
| r_acc / r_tidal | 0.007 | - |

The choice r_acc = 0.1 pc captures all tidal disruption physics while maintaining computational efficiency with timesteps of ~290 years.

## References

- Monaghan, J.J. (1992). Smoothed particle hydrodynamics. ARAA, 30, 543.
- Rees, M.J. (1988). Tidal disruption of stars by black holes. Nature, 333, 523.
- Guillochon, J. & Ramirez-Ruiz, E. (2013). Hydrodynamical simulations to determine the feeding rate of black holes by the tidal disruption of stars. ApJ, 767, 25.
