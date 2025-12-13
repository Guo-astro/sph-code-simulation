# Integration Checklist for IMBH-Cloud Simulation

## Status: ✅ Setup Complete, ⏳ Integration Pending

All code files have been created following the repository architecture. The following integration steps are needed before the code can be compiled and run.

---

## Phase 1: Build System Integration

### 1. Add External Forces to Build

**File**: `src/CMakeLists.txt`

**Action**: Add the following line to the list of subdirectories:

```cmake
add_subdirectory(external_forces)
```

**Location**: After other subdirectories like `add_subdirectory(gdisph)`, `add_subdirectory(thermal)`, etc.

**Verification**:
```bash
grep "add_subdirectory(external_forces)" src/CMakeLists.txt
```

---

## Phase 2: Sample Registration

### 2. Register IMBH-Cloud Sample

**File**: `src/solver.cpp`

**Action**: Add to the sample type dispatcher in `Solver::initialize()`:

Find the section that looks like:
```cpp
} else if (sample_type == "lane_emden") {
    make_lane_emden();
} else if (sample_type == "khi") {
    make_khi();
```

Add after it:
```cpp
} else if (sample_type == "imbh_cloud") {
    make_imbh_cloud();
```

**Verification**:
```bash
grep "make_imbh_cloud" src/solver.cpp
```

---

### 3. Declare IMBH-Cloud Function

**File**: `include/solver.hpp`

**Action**: Add declaration in the `Solver` class:

```cpp
void make_imbh_cloud();
```

**Location**: In the private section with other sample functions like `make_lane_emden()`, `make_khi()`, etc.

**Verification**:
```bash
grep "make_imbh_cloud" include/solver.hpp
```

---

## Phase 3: Build Configuration

### 4. Set Dimension to 3D

**File**: `include/defines.hpp`

**Action**: Ensure DIM is set to 3:

```cpp
#define DIM 3
```

**Note**: IMBH-cloud simulation requires 3D. You may want to keep 2D/1D configs for other tests.

**Verification**:
```bash
grep "define DIM" include/defines.hpp
```

---

## Phase 4: Build & Test

### 5. Clean Build

```bash
cd /Users/guo/Downloads/sphcode
make clean
make -j8
```

**Expected**: No compilation errors. The `build/sph` executable should be created.

**Troubleshooting**:
- If "external_forces not found": Check step 1 (CMakeLists.txt)
- If "make_imbh_cloud undefined": Check steps 2-3 (solver.cpp and solver.hpp)
- If dimension errors: Check step 4 (defines.hpp)

---

### 6. Test Initial Conditions Generation

```bash
make -f simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud imbh_check_dim
```

**Expected**: Pass dimension check (DIM=3)

---

### 7. Run Test Simulation (Small N)

**Action**: Create a test preset with fewer particles for quick validation:

**File**: `simulations/astrophysics/imbh_cloud/config/presets/imbh_cloud_test.json`

```json
{
  "name": "imbh_cloud_test",
  "description": "Quick test with N=1000",
  "sample": {
    "type": "imbh_cloud",
    "N": 1000,
    "M_cloud": 10000.0,
    "R_cloud": 5.0,
    "M_BH": 100000.0,
    "b": 3.0
  },
  "simulation": {
    "end_time": 0.1,
    "output_time": 0.01,
    "output_directory": "simulations/astrophysics/imbh_cloud/results/test"
  },
  "numerical": {
    "neighbor_number": 50,
    "kernel": "wendland",
    "sph_type": "gdisph",
    "use_gravity": true
  }
}
```

**Run**:
```bash
cp simulations/astrophysics/imbh_cloud/config/presets/imbh_cloud_test.json simulations/astrophysics/imbh_cloud/imbh_cloud.json
./build/sph imbh_cloud
```

**Expected**: Should complete in <1 minute and create snapshots in `simulations/astrophysics/imbh_cloud/results/test/`

---

## Phase 5: Validation (Optional but Recommended)

### 8. Energy Conservation Test

**Check**: `simulations/astrophysics/imbh_cloud/results/test/energy.dat`

**Validation**:
```python
import numpy as np
data = np.loadtxt('simulations/astrophysics/imbh_cloud/results/test/energy.dat')
time = data[:, 0]
E_total = data[:, 4]  # Total energy column
drift = (E_total[-1] - E_total[0]) / E_total[0]
print(f"Energy drift: {drift * 100:.2f}%")
```

**Expected**: < 1% drift

---

### 9. Thermal Equilibrium Check

**Load snapshot**:
```python
import numpy as np
data = np.loadtxt('simulations/astrophysics/imbh_cloud/results/test/snapshot_0001.csv', delimiter=',')
dens = data[:, 11]  # Density column
ene = data[:, 13]   # Energy column
n_H = dens * 40.5   # Convert to number density

# Plot n-T diagram and verify it follows K&I curve
```

---

## Phase 6: Production Run (After Validation)

### 10. Run Full Resolution

```bash
make -f simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud imbh_b3pc_run
```

**Expected**: ~12 hours for N=200,000 particles

**Monitor**:
```bash
tail -f simulations/astrophysics/imbh_cloud/results/b3pc_gdisph/energy.dat
```

---

### 11. Generate Visualizations

```bash
make -f simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud imbh_visualize
```

**Expected**:
- Plots in `simulations/astrophysics/imbh_cloud/results/plots/`
- Animations in `simulations/astrophysics/imbh_cloud/results/animations/`

---

## Phase 7: Parameter Exploration (Research Phase)

### 12. Impact Parameter Scan

```bash
make -f simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud imbh_scan_run
```

**Expected**: 4 runs (b = 3, 4, 5, 6 pc), ~2 days total

---

### 13. Resolution Convergence

**Manual**: Run with N = 50k, 200k, 1M and compare:
- Tidal disruption threshold
- Mass loss fraction
- Energy conservation quality

---

## Troubleshooting Common Issues

### Compilation Errors

| Error | Fix |
|-------|-----|
| `external_forces: No such file` | Add to src/CMakeLists.txt (step 1) |
| `make_imbh_cloud undefined` | Declare in solver.hpp and solver.cpp (steps 2-3) |
| `DIM != 3` assertion | Edit defines.hpp (step 4) |
| `KoyamaInutsukaCooling not found` | Thermal module already exists, check includes |

### Runtime Errors

| Error | Fix |
|-------|-----|
| `Cannot open Lane-Emden file` | Run from repository root, not subdirectory |
| `sample type imbh_cloud not recognized` | Register sample (step 2) |
| Segfault | Check N is reasonable (< 10M for 8GB RAM) |

### Physics Issues

| Issue | Check |
|-------|-------|
| Energy drift > 1% | Reduce timestep (CFL_sound, CFL_force) |
| No tidal elongation | Verify BH position (should be at x=b, y=0, z=0) |
| Cloud explodes | Check artificial viscosity (alpha ~ 1.0) |

---

## Completion Checklist

- [ ] Step 1: Add external_forces to CMakeLists.txt
- [ ] Step 2: Register sample in solver.cpp
- [ ] Step 3: Declare in solver.hpp
- [ ] Step 4: Set DIM=3 in defines.hpp
- [ ] Step 5: Clean build (make clean && make)
- [ ] Step 6: Pass dimension check
- [ ] Step 7: Test run with N=1000
- [ ] Step 8: Verify energy conservation
- [ ] Step 9: Check thermal equilibrium
- [ ] Step 10: Production run (optional)
- [ ] Step 11: Generate visualizations (optional)

---

## Quick Command Reference

```bash
# Build
make clean && make -j8

# Test
make -f simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud imbh_check_dim

# Run
make -f simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud imbh_b3pc_run

# Visualize
make -f simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud imbh_visualize

# Clean
make -f simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud imbh_clean
```

---

## Files Modified During Integration

| File | Action | Status |
|------|--------|--------|
| `src/CMakeLists.txt` | Add subdirectory | ⏳ Pending |
| `src/solver.cpp` | Register sample | ⏳ Pending |
| `include/solver.hpp` | Declare function | ⏳ Pending |
| `include/defines.hpp` | Set DIM=3 | ⏳ Pending |

**All other files are ready (created by setup).**

---

## Next Steps After Integration

1. **Validate**: Run test simulation (N=1000) and check conservation
2. **Benchmark**: Run N=50k for performance baseline
3. **Explore**: Parameter scan (b = 3-6 pc)
4. **Analyze**: Generate plots and interpret physics
5. **Publish**: Write up results with references to K&I (2000)

---

**Questions?** See:
- `simulations/astrophysics/imbh_cloud/README.md` for usage
- `simulations/astrophysics/imbh_cloud/docs/RESEARCH_SETUP.md` for physics
- `simulations/astrophysics/imbh_cloud/docs/IMPLEMENTATION_SUMMARY.md` for architecture
