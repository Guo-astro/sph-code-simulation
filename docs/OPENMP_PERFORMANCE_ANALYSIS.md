# OpenMP Parallelization Analysis for SPH Code

## Current Status: ✅ **GOOD** (but can be optimized)

### What's Already Parallelized

The code **IS** using OpenMP parallelization across all critical sections:

#### 1. **Relaxation Loop** (`src/solver.cpp`)
- ✅ Position integration (line 1097): `#pragma omp parallel for`
- ✅ Velocity zeroing (same loop)
- ✅ Periodic boundary application

#### 2. **Force Calculations**
- ✅ Lane-Emden relaxation forces (`src/relaxation/lane_emden_relaxation.cpp:121`)
- ✅ GDISPH pre-interaction (`src/gdisph/gd_pre_interaction.cpp:57`)
- ✅ GDISPH fluid forces (`src/gdisph/gd_fluid_force.cpp:77`)
- ✅ Gravity forces (`src/gravity_force.cpp:66`)

#### 3. **Tree Building** (`src/bhtree.cpp`)
- ✅ Node initialization (line 53)
- ✅ Bounding box calculation (line 75)
- ✅ Morton key computation (line 105)

#### 4. **Other Operations**
- ✅ Timestep calculation (`src/timestep.cpp:24`)
- ✅ Energy calculations (solver.cpp lines 1764, 1789)
- ✅ Output manager energy calculation (`src/output_manager.cpp`)

### Current Configuration

**From simulation output:**
```
Open MP is valid.
the number of threads = 8
```

**System specs:**
- Physical cores: 8
- Logical cores: 8
- Current OMP_NUM_THREADS: **not set** (defaults to system maximum)

## Performance Analysis

### ✅ What's Working Well

1. **All major loops are parallelized**
   - Force calculations
   - Position updates
   - Tree operations
   - Energy calculations

2. **Thread count is optimal**
   - Using 8 threads on 8-core system

3. **No obvious serial bottlenecks** in computational sections

### ⚠️ Potential Bottlenecks

#### 1. **Sequential Tree Rebuilding**
Every relaxation step (line 1063-1073):
```cpp
for(int step = start_step; step < target_step; ++step) {
    tree->resize(num_p);        // Sequential
    tree->make(particles, num_p);  // Partially parallel
    m_pre->calculation(m_sim);   // Parallel ✓
    m_fforce->calculation(m_sim); // Parallel ✓
    // ... position update (parallel) ✓
}
```

**The main loop itself is SEQUENTIAL** - each iteration must complete before the next begins. This is correct for physics (time-stepping), but means parallelization is within each step, not across steps.

#### 2. **I/O Operations**
- Snapshot writing (every 1000 steps): Sequential, single-threaded
- Energy file writing: Sequential

#### 3. **Tree Construction**
While Morton key computation is parallel, tree node linking may have serial sections.

#### 4. **Memory Bandwidth**
With 10,648 particles and 8 threads, each thread processes ~1,331 particles.
- Good for computation (enough work per thread)
- Possible cache contention if particles accessed randomly

## Performance Metrics

### Expected Speedup

**Theoretical (Amdahl's Law):**
- If 95% of work is parallelizable across 8 cores:
  - Speedup = 1 / (0.05 + 0.95/8) = **5.9x**
- If 90% parallelizable:
  - Speedup = 1 / (0.10 + 0.90/8) = **4.7x**

**Actual speedup depends on:**
- Memory bandwidth (shared among cores)
- Cache efficiency
- Load balancing
- Tree construction overhead

### CPU Utilization

**Expected**: 600-750% CPU usage (on macOS Activity Monitor)
- 8 threads × ~75-95% efficiency = 6-7.6 cores actively used

**Factors reducing efficiency:**
- I/O waiting
- Memory bandwidth saturation
- Load imbalance (some particles have more neighbors)
- Tree rebuilding overhead

## Optimization Recommendations

### 1. **Set OMP Environment Variables** (EASY, IMMEDIATE IMPACT)

Add to your shell profile (`~/.zshrc`):
```bash
export OMP_NUM_THREADS=8           # Use all cores
export OMP_SCHEDULE="dynamic,100"  # Dynamic load balancing
export OMP_PROC_BIND=close         # Pin threads to cores
export OMP_PLACES=cores            # Thread placement strategy
```

Or set before running:
```bash
OMP_NUM_THREADS=8 OMP_SCHEDULE="dynamic" make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_relax
```

### 2. **Reduce Tree Rebuilding Frequency** (MODERATE IMPACT)

Current: Tree rebuilt EVERY step (line 1068-1072)

**Optimization**: Rebuild only when needed
```cpp
if(step % 10 == 0) {  // Rebuild every 10 steps instead of every step
    tree->resize(num_p);
    tree->make(particles, num_p);
}
```

**Tradeoff**: Slightly less accurate neighbor lists, but 10x fewer tree builds

### 3. **Batch Snapshot Writes** (SMALL IMPACT)

Current: Snapshot written inside parallel section

**Optimization**: Accumulate snapshots, write in batch at end

### 4. **Use NUMA-Aware Allocation** (ADVANCED, macOS limited)

macOS doesn't have traditional NUMA, but first-touch allocation can help:
```cpp
#pragma omp parallel for
for(int i = 0; i < num_p; ++i) {
    particles[i].pos[0] = initial_pos[i][0];  // Initialize in parallel
}
```

### 5. **Reduce Output Frequency During Debugging**

Change `relaxationOutputFreq` from 1000 to 5000 or 10000 to reduce I/O overhead.

## Verification Commands

### Check CPU Usage During Run

**Terminal 1:**
```bash
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_relax_resume \
     SNAPSHOT=sample/imbh_cloud/results/lane_emden_50k_relax/snapshot_0021.csv \
     STEPS=10000 FREQ=5000
```

**Terminal 2 (while running):**
```bash
# Check CPU usage
top -pid $(pgrep sph) -stats pid,cpu,threads,mem

# Or use Activity Monitor.app and search for "sph"
```

**Expected**: 
- CPU: 600-750% (on 8-core system)
- Threads: 9-10 (1 main + 8 worker threads)

### Profile with Instruments (macOS)

```bash
# Build with debug symbols
cd build
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -DSPH_DIM=3 ..
make -j8

# Profile
instruments -t "Time Profiler" ./sph ../sample/imbh_cloud/config/presets/lane_emden_50k_relax.json
```

## Benchmark: Single vs Multi-threaded

Create test script:
```bash
#!/bin/bash
# Quick benchmark
CONFIG="sample/imbh_cloud/config/presets/lane_emden_50k_relax.json"

echo "Single-threaded (OMP_NUM_THREADS=1):"
time OMP_NUM_THREADS=1 ./build/sph $CONFIG

echo -e "\nMulti-threaded (OMP_NUM_THREADS=8):"
time OMP_NUM_THREADS=8 ./build/sph $CONFIG
```

**Expected result**: 4-6x speedup with 8 threads

## Current Performance Estimate

For 10,000 relaxation steps with 10,648 particles:

- **Single-threaded**: ~20-30 minutes
- **8 threads (current)**: ~4-6 minutes ✓
- **Optimized 8 threads**: ~3-4 minutes (with recommendations)

## Conclusion

### ✅ **Current Status: GOOD**

The code is **properly parallelized** and using **all 8 CPU cores**. The main bottlenecks are:

1. **Physics-imposed serial time-stepping** (unavoidable)
2. **Tree rebuilding every step** (can be optimized)
3. **Memory bandwidth** (hardware limitation)

### 🚀 **Quick Win Recommendations**

1. **Set OMP environment variables** (2 min, 10-15% improvement)
2. **Reduce tree rebuild frequency** (10 min coding, 30-50% improvement)
3. **Increase output frequency** (1 min, 5-10% improvement for I/O-heavy runs)

### 📊 **Verify Performance**

Run with monitoring:
```bash
OMP_NUM_THREADS=8 OMP_SCHEDULE=dynamic \
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_relax_resume \
     SNAPSHOT=sample/imbh_cloud/results/lane_emden_50k_relax/snapshot_0021.csv \
     STEPS=1000 FREQ=1000
```

Check in Activity Monitor:
- Should see ~600-750% CPU usage
- Memory should be stable (~200-300 MB)

---

**Last updated**: 2025-12-01  
**Analysis by**: Guo  
**System**: 8-core macOS
