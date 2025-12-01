# Complete IMBH-Cloud Workflow with Hydrostatic Testing

## Quick Reference: 3-Step Process

```bash
# From repository root: /Users/guo/Downloads/sphcode/

# STEP 1: Relaxation (create equilibrium initial conditions)
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_relax_2k
# Runtime: ~5-10 minutes
# Output: sample/imbh_cloud/results/lane_emden_2k_relax/

# STEP 2: Hydrostatic Test (verify self-gravity implementation)
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_hydrostatic_2k
# Runtime: ~5-10 minutes
# Output: sample/imbh_cloud/results/lane_emden_2k_hydrostatic/

# STEP 3: IMBH Disruption (proceed with science runs)
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_b3pc_run
# Ready for production simulations!
```

## Why Each Step is Necessary

### Step 1: Relaxation
**Purpose:** Remove initial transients and settle particles into hydrostatic equilibrium.

**What happens:**
- Particles start from analytical Lane-Emden positions
- Small random velocities cause initial motion
- Analytical gravity force (no self-gravity yet) brings system to rest
- After ~20k-400k steps, velocities → 0

**Output:** Relaxed snapshot with particles at rest in correct density profile

### Step 2: Hydrostatic Test (NEW!)
**Purpose:** Verify self-gravity works correctly before expensive IMBH runs.

**What happens:**
- Load relaxed particles from Step 1
- Turn ON self-gravity (Barnes-Hut tree)
- Evolve for 10 crossing times with DISPH hydrodynamics
- Monitor density profile, velocities, energy

**Critical validation:**
- ✓ Density stays constant → gravity balances pressure
- ✓ Velocities stay << 1% c_s → no spurious forces
- ✓ Energy conserved → numerical stability
- ✓ Duration >> tidal timescale → relevant for IMBH physics

**Why this matters:**
If self-gravity is incorrect (wrong sign, wrong normalization, tree bugs), the sphere will:
- Collapse (gravity too strong)
- Expand (gravity too weak)
- Oscillate (spurious forces)
- Drift energy (poor conservation)

**Catching bugs here saves:**
- Days of wasted IMBH simulation time
- Incorrect science results
- Debugging expensive production runs

### Step 3: IMBH Disruption
**Purpose:** Science! Study tidal disruption with verified gravity.

**Uses:**
- Relaxed initial conditions (Step 1)
- Verified self-gravity (Step 2)
- Add IMBH external force
- Watch tidal disruption unfold

## Testing vs Production

### Quick Testing (2k particles)

**Full workflow (~15-25 minutes):**
```bash
# Step 1: Relaxation (5-10 min)
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_relax_2k

# Step 2: Hydrostatic test (5-10 min)
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_hydrostatic_2k

# Step 3: Quick IMBH test (5-10 min)
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_b3pc_run
```

**Use cases:**
- Method development
- Code testing
- Parameter exploration
- Bug fixing

### Production (200k particles)

**Full workflow (~20-32 hours):**
```bash
# Step 1: High-resolution relaxation (4-8 hours)
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_relax_200k

# Step 2: Production hydrostatic test (8-16 hours)
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_hydrostatic_200k

# Step 3: Production IMBH runs (8-16 hours each)
# Now safe to run expensive production simulations!
```

**Use cases:**
- Publication-quality results
- High-resolution studies
- Jeans mass resolution (N_J > 50-100)

## Visualization Workflow

### Check Each Step

```bash
# After Step 1: Verify relaxation quality
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_relax_2k_viz
open sample/imbh_cloud/results/lane_emden_2k_relax/relaxation_animation.gif

# After Step 2: Verify hydrostatic stability
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_hydrostatic_2k_viz
open sample/imbh_cloud/results/lane_emden_2k_hydrostatic/hydrostatic_animation.gif

# After Step 3: Analyze IMBH disruption
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_visualize
```

## Quality Checklist

### ✅ Step 1 Complete When:
- [ ] Velocities converged to << 0.001 code units
- [ ] Density profile matches Lane-Emden analytic
- [ ] Energy file shows smooth evolution
- [ ] Animation shows particles settling to rest

### ✅ Step 2 Complete When:
- [ ] Test ran for full 10.0 code units (10 t_cross)
- [ ] Density RMS error < 5%
- [ ] Max velocity < 1% sound speed
- [ ] Median velocity < 0.1% sound speed
- [ ] Energy drift < 1%
- [ ] Animation shows stable density profile
- [ ] Summary statistics show PASS for all criteria

### ✅ Step 3 Ready When:
- [ ] Steps 1 & 2 both passed
- [ ] Hydrostatic test verified self-gravity
- [ ] Parameters chosen (impact parameter, velocity)
- [ ] Sufficient disk space for output

## Troubleshooting Decision Tree

### Problem: Step 1 fails (relaxation)

**Symptoms:**
- Velocities not converging
- Density profile wrong
- Particles escaping

**Actions:**
1. Increase relaxation steps: Edit `lane_emden_2k_relax.json`
2. Check gamma = 5/3 (not 1.666...)
3. Check polytropic constant K
4. Verify Lane-Emden parameters (ξ₁, α, ρ_c)

### Problem: Step 2 fails (hydrostatic test)

**Symptoms:**
- Density error > 10%
- Velocities > 5% c_s
- Energy drift > 5%
- Sphere expanding or collapsing

**Actions:**
1. **First:** Check Step 1 relaxation quality
2. Verify gravity constant G = 1.0 in code units
3. Check tree opening angle θ
4. Reduce timestep (lower CFL)
5. See detailed guide: `docs/HYDROSTATIC_TEST.md`

**DO NOT proceed to Step 3 if Step 2 fails!**

### Problem: Step 3 fails (IMBH disruption)

**Symptoms:**
- Unphysical behavior
- Particle explosions
- Energy non-conservation

**Actions:**
1. **First:** Verify Steps 1 & 2 passed
2. Check IMBH force parameters
3. Verify impact parameter and velocity
4. Check boundary conditions

## File Locations

### Presets (Configuration Files)
```
sample/imbh_cloud/config/presets/
├── lane_emden_2k_relax.json         # Step 1 (testing)
├── lane_emden_200k_relax.json       # Step 1 (production)
├── lane_emden_2k_hydrostatic.json   # Step 2 (testing)
├── lane_emden_200k_hydrostatic.json # Step 2 (production)
└── imbh_cloud_b3pc_gdisph.json      # Step 3 (IMBH runs)
```

### Results
```
sample/imbh_cloud/results/
├── lane_emden_2k_relax/        # Step 1 output
├── lane_emden_2k_hydrostatic/  # Step 2 output
└── disruption_runs/            # Step 3 output
```

### Scripts
```
sample/imbh_cloud/scripts/
├── visualize_relaxation.py     # Step 1 visualization
├── visualize_hydrostatic.py    # Step 2 visualization
└── visualize_disruption.py     # Step 3 visualization
```

## Recommended Development Cycle

### Initial Development (Use 2k)
```bash
# Quick iteration during development
while [ code_being_developed ]; do
    make clean
    make -j8
    make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_relax_2k
    make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_hydrostatic_2k
    # Check results, fix bugs, repeat
done
```

### Pre-Production Validation (Use 2k)
```bash
# Once code stable, run full 2k validation
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_relax_2k
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_hydrostatic_2k

# Thoroughly check results
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_relax_2k_viz
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_hydrostatic_2k_viz

# Verify all pass criteria met
```

### Production Runs (Use 200k)
```bash
# Only after 2k validation passes!
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_relax_200k
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_hydrostatic_200k

# Then proceed with expensive science runs
```

## Time Estimates

### 2k Particle Timeline
| Step | Time | Cumulative |
|------|------|------------|
| Relaxation | 5-10 min | 10 min |
| Hydrostatic test | 5-10 min | 20 min |
| Visualization | 2-5 min | 25 min |
| **Total** | **~25 min** | Ready for IMBH! |

### 200k Particle Timeline
| Step | Time | Cumulative |
|------|------|------------|
| Relaxation | 4-8 hours | 8 hours |
| Hydrostatic test | 8-16 hours | 24 hours |
| Visualization | 5-10 min | 24 hours |
| **Total** | **~24 hours** | Ready for production! |

**Recommendation:** Run 200k overnight or over weekend.

## Help & Documentation

```bash
# Quick help
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_help

# Full documentation
cat sample/imbh_cloud/README.md
cat sample/imbh_cloud/docs/HYDROSTATIC_TEST.md

# Examples
ls sample/imbh_cloud/config/presets/
```

## Summary

**The hydrostatic test is NOT optional!**

It's a critical validation step that:
- Catches gravity implementation bugs
- Verifies numerical stability
- Ensures accuracy for IMBH physics
- Saves time by failing fast

**Workflow:**
1. Relax → particles at rest
2. Test → gravity verified
3. Science → confident results

---

**Last updated:** 2025-12-01  
**Author:** Guo
