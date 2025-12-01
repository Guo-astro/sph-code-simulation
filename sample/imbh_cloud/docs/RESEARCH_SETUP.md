# IMBH-Molecular Cloud Interaction: Research Setup

## Scientific Objective

Simulate intermediate-mass black hole (IMBH) tidal disruption of molecular clouds with the following parameters:

### Physical Parameters

- **IMBH Mass**: M_BH = 10⁵ M_☉ (100,000 solar masses)
- **Impact Parameters**: b = 3-6 parsecs
- **Cloud Structure**: 3D Lane-Emden sphere with polytropic index n = 5/3 (γ = 5/3)
- **Thermal Physics**: Koyama & Inutsuka (2000) thermal equilibrium curve
  - Assumption: Thermal timescale << dynamical timescale
  - Cloud remains in thermal equilibrium throughout simulation

## Resolution Analysis

### 1. Jeans Mass Criterion

The Jeans mass for a self-gravitating cloud is:

```
M_J = (π^(5/2) / 6) * (c_s^3 / (G^(3/2) * ρ^(1/2)))
```

For thermal equilibrium in CNM (Cold Neutral Medium):
- n_H ~ 100-1000 cm⁻³ → T_eq ~ 10-50 K (from K&I thermal curve)
- Sound speed: c_s ~ 0.2-0.3 km/s

**Jeans mass**: M_J ~ 1-10 M_☉ (typical for molecular clouds)

**Requirement**: To resolve fragmentation and thermal pressure support:
- **N_J ≥ 50-100 particles per Jeans mass**
- For cloud mass M_cloud ~ 10³-10⁴ M_☉:
  - **Minimum N ≥ 5×10⁴ - 10⁶ particles**

### 2. Tidal Disruption Timescale

The tidal (pancake) disruption occurs when tidal force exceeds self-gravity:

```
t_tidal ~ √(R_cloud³ / (G M_BH))
```

For M_BH = 10⁵ M_☉ and R_cloud ~ 5-10 pc:
- **t_tidal ~ 10⁴ - 10⁵ years**

Compare to cloud crossing time:
```
t_cross ~ R_cloud / v_rel
```

For v_rel ~ 10-50 km/s:
- **t_cross ~ 10⁵ - 10⁶ years**

**Timescale ordering**: t_thermal << t_tidal ~ t_cross

### 3. Spatial Resolution Requirements

**Minimum spatial resolution needed**:

a) **Tidal radius** (where BH gravity ~ cloud self-gravity):
   ```
   r_t ~ R_cloud * (M_BH / M_cloud)^(1/3)
   ```
   For M_cloud ~ 10⁴ M_☉, M_BH = 10⁵ M_☉:
   - r_t ~ 2-5 R_cloud
   - **Need to resolve ~0.01-0.1 R_cloud** to capture tidal effects

b) **Hill radius** (tidal truncation):
   ```
   r_H ~ b * (M_cloud / (3 M_BH))^(1/3)
   ```
   For b = 3-6 pc:
   - r_H ~ 0.5-1 pc
   - **Need N_neighbor ~ 50-100 within r_H**

c) **Smoothing length constraint**:
   - h ~ R_cloud / √N for uniform distribution
   - For N = 10⁵: h ~ 0.01 R_cloud ✓
   - For N = 10⁶: h ~ 0.003 R_cloud ✓

### 4. Recommended Particle Numbers

| Cloud Mass | N_particles | h/R_cloud | N_J | Purpose |
|------------|-------------|-----------|-----|---------|
| 10³ M_☉   | 5×10⁴       | ~0.014    | 50  | Quick test |
| 10⁴ M_☉   | 2×10⁵       | ~0.007    | 100 | Standard resolution |
| 10⁴ M_☉   | 10⁶         | ~0.003    | 500 | High resolution |

**Recommended**: Start with **N = 2×10⁵** for production runs, **N = 5×10⁴** for parameter studies.

### 5. Thermal Equilibrium Constraint

From Koyama & Inutsuka (2000) cooling timescale:

```
t_cool ~ u / |du/dt| ~ u / |(u_eq - u) / τ_relax|
```

For τ_relax ~ 0.1 Myr (short compared to t_tidal):
- **Thermal equilibrium maintained**: T → T_eq(n) within ~1% of t_tidal

**Implementation**: Use `KoyamaInutsukaCooling` with relaxation time:
```cpp
real tau_relax = 0.01 * t_tidal;  // 1% of tidal timescale
real du_dt = cooling.cooling_rate(n_H, T_current, tau_relax);
```

## Physics Modules Required

### 1. External Force: Point-Mass Black Hole

**Header**: `include/external_forces/point_mass_bh.hpp`

```cpp
namespace sph {
namespace external_forces {

class PointMassBlackHole {
    vec_t m_position;      // BH position [pc]
    real m_mass;           // BH mass [M_☉]
    real m_softening;      // Gravitational softening [pc]
    
public:
    vec_t acceleration(const vec_t& r_particle) const;
    real potential(const vec_t& r_particle) const;
};

}
}
```

**Force**:
```
F_BH = -G M_BH (r - r_BH) / (|r - r_BH|² + ε²)^(3/2)
```

where ε = softening length ~ 0.01 R_cloud (prevents singularity).

### 2. Thermal Physics: K&I Cooling

**Existing module**: `include/thermal/koyama_inutsuka_cooling.hpp`

- Already implemented ✓
- Provides T_eq(n_H), P_eq(n_H)
- Cooling rate: du/dt = (u_eq - u) / τ_relax

### 3. Initial Conditions: Lane-Emden + Thermal Equilibrium

**Based on**: `src/sample/lane_emden.cpp`

Modifications needed:
1. Set initial T = T_eq(n) from K&I curve (not polytro

pic P = K ρ^γ)
2. Add IMBH at offset position
3. Optionally add relative velocity

### 4. SPH Method: GDISPH (Recommended)

**Why GDISPH**:
- Best for density gradients (cloud disruption creates steep gradients)
- Hybrid Godunov + pressure-energy formulation
- Built-in Riemann solver handles shocks from tidal compression
- Already tested in `sample/sedov/`, `sample/khi/`

**Configuration**:
```json
{
  "numerical": {
    "sph_type": "gdisph",
    "kernel": "wendland",
    "neighbor_number": 50,
    "use_gravity": true
  },
  "artificial_viscosity": {
    "alpha": 1.0,
    "use_balsara_switch": true,
    "use_time_dependent_av": true
  }
}
```

## Code Units and Conversion

### Recommended Unit System

**Length unit**: 1 code unit = 1 parsec (pc)
**Mass unit**: 1 code unit = 1 M_☉
**Time unit**: Derived from G

```
G = 4.302 × 10⁻³ pc (km/s)² M_☉⁻¹
t_unit = √(L³ / (G M)) = √(pc³ / (G M_☉))
       ≈ 0.978 Myr
```

**Velocity unit**: v_unit = L/t_unit ≈ 1.02 km/s

### Physical Parameters in Code Units

| Quantity | Physical Value | Code Units |
|----------|----------------|------------|
| M_BH | 10⁵ M_☉ | 10⁵ |
| M_cloud | 10⁴ M_☉ | 10⁴ |
| R_cloud | 5 pc | 5.0 |
| b (impact) | 3-6 pc | 3.0-6.0 |
| v_rel | 20 km/s | ~20 |
| T_CNM | 50 K | — |
| n_H | 100 cm⁻³ | — |
| t_end | 1 Myr | ~1.0 |

### Density Conversion

From code density ρ [M_☉/pc³] to number density n_H [cm⁻³]:

```
n_H = ρ * (M_☉/pc³) / (μ m_H)
    = ρ * 40.5  [for μ = 1, neutral H]
```

where:
- μ = 1 (neutral hydrogen)
- m_H = 1.67 × 10⁻²⁴ g
- 1 pc³ = 2.938 × 10⁵⁷ cm³

## Implementation Checklist

### Phase 1: External Force Module
- [ ] `include/external_forces/point_mass_bh.hpp`
- [ ] `src/external_forces/point_mass_bh.cpp`
- [ ] Add to CMakeLists.txt
- [ ] Unit tests

### Phase 2: Initial Conditions
- [ ] `src/sample/imbh_cloud.cpp`
  - [ ] Lane-Emden sphere with n=5/3
  - [ ] K&I thermal equilibrium T_eq(n)
  - [ ] IMBH placement at distance b
  - [ ] Optional relative velocity
  
### Phase 3: Configuration Presets
- [ ] `sample/imbh_cloud/config/presets/imbh_cloud_b3pc_gdisph.json`
- [ ] `sample/imbh_cloud/config/presets/imbh_cloud_b6pc_gdisph.json`
- [ ] Parameter scan configs (b = 3, 4, 5, 6 pc)

### Phase 4: Visualization & Analysis
- [ ] `sample/imbh_cloud/scripts/visualize_disruption.py`
  - [ ] Density evolution
  - [ ] Tidal deformation (axis ratios)
  - [ ] Mass loss rate
  - [ ] Temperature vs density phase diagram
  
### Phase 5: Makefile
- [ ] `sample/imbh_cloud/Makefile.imbh_cloud`
  - [ ] Single run targets
  - [ ] Parameter scan targets
  - [ ] Visualization targets

## Expected Physics

### Tidal Disruption Sequence

1. **Initial approach** (t = 0 - 0.3 t_tidal):
   - Cloud feels increasing tidal gradient
   - Elongation along radial direction
   - Compression perpendicular to orbital plane (pancake)

2. **Maximum compression** (t ~ 0.5 t_tidal):
   - Peak density reached at pericenter passage
   - Shocks from rapid compression (if v_rel large)
   - Potential triggered star formation (future work)

3. **Disruption** (t > t_tidal):
   - Leading/trailing tidal tails form
   - Bound vs unbound material separation
   - Self-gravity determines remnant mass

4. **Thermal response**:
   - Compression → n increases → T_eq decreases (CNM branch)
   - Cooling maintains thermal equilibrium
   - No shock heating runaway (unlike adiabatic case)

### Diagnostics to Track

1. **Morphology**:
   - Axis ratios (a:b:c) from moment of inertia tensor
   - Elongation parameter: e = 1 - c/a

2. **Energetics**:
   - Kinetic energy
   - Thermal energy (should track T_eq)
   - Gravitational potential (cloud + BH)
   - Total energy (check conservation)

3. **Mass Budget**:
   - Bound mass (E < 0)
   - Unbound mass (E > 0)
   - Accretion onto BH (within R_acc ~ 0.01 pc)

4. **Thermal State**:
   - n-T diagram (should follow K&I curve)
   - Departures from equilibrium (measure τ_cool)

## References

1. Koyama & Inutsuka (2000), ApJ 532, 980: Thermal equilibrium curves
2. Guillochon & Loeb (2015), ApJ 811, 20: Tidal disruption events
3. Burkert & Naab (2013), MNRAS 434, 36: Cloud-BH interactions
4. Stone et al. (1998), ApJS 114, 345: Numerical methods for tidal disruption

## Next Steps

1. Review existing gravity force implementation in `src/gravity_force.cpp`
2. Implement external point-mass force module
3. Adapt Lane-Emden setup to include K&I thermal equilibrium
4. Create test case with N = 5×10⁴ particles
5. Validate conservation and thermal equilibrium before production runs
