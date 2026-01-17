# Tree Optimization Implementation Plan

## Overview

**Branch:** `feature/optimize-tree-remove-exhaustive-search`
**Goal:** Achieve 4-10x speedup through three major optimizations:

1. **Morton Code Particle Reordering** - 2-3x speedup via cache-friendly memory access
2. **Iterative Tree Traversal** - 10-20% speedup via stack elimination
3. **SIMD Batch Distance Calculations** - 2-4x speedup via vectorized operations

**Approach:** TDD - Use exhaustive search as ground truth to validate tree implementation before removal

---

## Current Codebase Analysis

### Tree Structure (`include/bhtree.hpp`)
- `BHNode`: ~160-200 bytes per node with `children[NCHILD]` pointers
- NCHILD = 8 (3D), 4 (2D), 2 (1D)
- Recursive traversal in `neighbor_search()` and `calc_force()`

### Particle Layout (`include/particle.hpp`)
- AoS (Array-of-Structures): ~200-400 bytes per particle
- Critical fields for neighbor search: `pos[DIM]`, `sml`, `id`
- Many unused fields during traversal cause cache pollution

### Distance Calculations
- Scalar calculations in tight loops
- ~50-500 neighbors per particle × N particles = millions of distance calcs/timestep
- Periodic boundary handling with `minimum_image()`

---

## Architecture Design

### 1. Morton Code Particle Reordering

**Goal**: Improve cache locality by reordering particles along space-filling curve

**Design**:
```
┌─────────────────────────────────────────────────────────────────┐
│                    Morton Code Architecture                      │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  ┌──────────────┐    ┌───────────────┐    ┌─────────────────┐  │
│  │ Particle     │───▶│ Morton Code   │───▶│ Sort Particles  │  │
│  │ Positions    │    │ Calculator    │    │ by Morton Code  │  │
│  └──────────────┘    └───────────────┘    └─────────────────┘  │
│                                                  │               │
│                                                  ▼               │
│  ┌──────────────┐    ┌───────────────┐    ┌─────────────────┐  │
│  │ Reorder      │◀───│ Build Index   │◀───│ Sorted Morton   │  │
│  │ Particle     │    │ Mapping       │    │ Keys            │  │
│  │ Arrays       │    │ old→new       │    │                 │  │
│  └──────────────┘    └───────────────┘    └─────────────────┘  │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

**New Files**:
- `include/morton.hpp` - Morton code computation
- `src/morton.cpp` - Morton sorting and reordering

**Key Functions**:
```cpp
namespace sph {
namespace morton {
    // Compute Morton code from normalized position
    uint64_t encode(const real pos[DIM], const real domain_min[DIM],
                    const real domain_size[DIM]);

    // Sort particles by Morton code (returns permutation)
    std::vector<int> sort_by_morton(const std::vector<SPHParticle>& particles,
                                    const real domain_min[DIM],
                                    const real domain_size[DIM]);

    // Apply permutation to reorder particles
    void reorder_particles(std::vector<SPHParticle>& particles,
                          const std::vector<int>& permutation);
}
}
```

**Integration Points**:
- Call after tree construction in `Simulation::build_tree()`
- Update neighbor lists to use new particle indices
- Reorder once per tree rebuild (not every timestep)

---

### 2. Iterative Tree Traversal

**Goal**: Eliminate recursion overhead with explicit stack

**Design**:
```
┌─────────────────────────────────────────────────────────────────┐
│               Iterative Traversal Architecture                   │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  Current (Recursive):           Proposed (Iterative):           │
│  ┌─────────────────┐            ┌─────────────────┐             │
│  │ neighbor_search │            │ neighbor_search │             │
│  │   if leaf:      │            │   stack.push(r) │             │
│  │     check       │    ───▶    │   while !empty: │             │
│  │   else:         │            │     node=pop()  │             │
│  │     recurse(c)  │            │     if leaf:    │             │
│  └─────────────────┘            │       check     │             │
│                                 │     else:       │             │
│                                 │       push kids │             │
│                                 └─────────────────┘             │
│                                                                  │
│  Stack Implementation:                                           │
│  ┌─────────────────────────────────────────────────────────┐    │
│  │ thread_local std::array<BHNode*, 64> stack;             │    │
│  │ // 64 levels = 2^64 nodes max (more than enough)        │    │
│  └─────────────────────────────────────────────────────────┘    │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

**Modified Functions** (in `src/bhtree.cpp`):
- `neighbor_search()` → `neighbor_search_iterative()`
- `calc_force()` → `calc_force_iterative()`

**Key Implementation**:
```cpp
void BHTree::neighbor_search_iterative(
    const SPHParticle& p_i,
    std::vector<int>& neighbor_list,
    int& n_neighbor,
    int max_neighbors,
    const bool is_ij,
    const Periodic* periodic) const
{
    thread_local std::array<BHNode*, 64> stack;
    int stack_top = 0;
    stack[stack_top++] = &m_root;

    while (stack_top > 0) {
        BHNode* node = stack[--stack_top];

        // AABB overlap test
        if (!overlaps_search_region(node, p_i.pos, kernel_size)) continue;

        if (node->is_leaf) {
            // Check particles in leaf
            for (auto* p = node->first; p != nullptr; p = p->next) {
                // Distance check and add to neighbor list
            }
        } else {
            // Push children (reverse order for depth-first)
            for (int c = NCHILD - 1; c >= 0; --c) {
                if (node->childs[c]) {
                    stack[stack_top++] = node->childs[c];
                }
            }
        }
    }
}
```

---

### 3. SIMD Batch Distance Calculations

**Goal**: Vectorize distance calculations using AVX2/AVX-512

**Design**:
```
┌─────────────────────────────────────────────────────────────────┐
│                  SIMD Distance Architecture                      │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  ┌──────────────────────────────────────────────────────────┐   │
│  │ Batch Buffer (process 4/8 candidates at once)            │   │
│  ├──────────────────────────────────────────────────────────┤   │
│  │ candidate_x[8]: [x0, x1, x2, x3, x4, x5, x6, x7]        │   │
│  │ candidate_y[8]: [y0, y1, y2, y3, y4, y5, y6, y7]        │   │
│  │ candidate_z[8]: [z0, z1, z2, z3, z4, z5, z6, z7]        │   │
│  │ candidate_id[8]: [id0, id1, ...]                         │   │
│  └──────────────────────────────────────────────────────────┘   │
│                              │                                   │
│                              ▼                                   │
│  ┌──────────────────────────────────────────────────────────┐   │
│  │ SIMD Distance Calculation                                │   │
│  ├──────────────────────────────────────────────────────────┤   │
│  │ __m256d dx = _mm256_sub_pd(target_x, candidate_x);      │   │
│  │ __m256d dy = _mm256_sub_pd(target_y, candidate_y);      │   │
│  │ __m256d dz = _mm256_sub_pd(target_z, candidate_z);      │   │
│  │ __m256d r2 = dx*dx + dy*dy + dz*dz;                     │   │
│  │ __m256d mask = _mm256_cmp_pd(r2, h2, _CMP_LT_OQ);       │   │
│  └──────────────────────────────────────────────────────────┘   │
│                              │                                   │
│                              ▼                                   │
│  ┌──────────────────────────────────────────────────────────┐   │
│  │ Compact Results (only neighbors within kernel)           │   │
│  ├──────────────────────────────────────────────────────────┤   │
│  │ int mask_bits = _mm256_movemask_pd(mask);               │   │
│  │ for each set bit: neighbor_list[count++] = candidate_id │   │
│  └──────────────────────────────────────────────────────────┘   │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

**New Files**:
- `include/simd_distance.hpp` - SIMD distance calculation utilities

**Key Implementation**:
```cpp
#include <immintrin.h>

namespace sph {
namespace simd {

#ifdef SPH_HAS_AVX2
// Process 4 double-precision distances at once
inline void batch_distance_check_avx2(
    const double target_x, const double target_y, const double target_z,
    const double* __restrict__ cand_x,  // aligned arrays of 4
    const double* __restrict__ cand_y,
    const double* __restrict__ cand_z,
    const double h2,  // kernel radius squared
    int* __restrict__ results,  // output: which candidates are neighbors
    int& count)
{
    __m256d tx = _mm256_set1_pd(target_x);
    __m256d ty = _mm256_set1_pd(target_y);
    __m256d tz = _mm256_set1_pd(target_z);
    __m256d h2_vec = _mm256_set1_pd(h2);

    __m256d cx = _mm256_load_pd(cand_x);
    __m256d cy = _mm256_load_pd(cand_y);
    __m256d cz = _mm256_load_pd(cand_z);

    __m256d dx = _mm256_sub_pd(tx, cx);
    __m256d dy = _mm256_sub_pd(ty, cy);
    __m256d dz = _mm256_sub_pd(tz, cz);

    __m256d r2 = _mm256_fmadd_pd(dx, dx,
                  _mm256_fmadd_pd(dy, dy,
                   _mm256_mul_pd(dz, dz)));

    __m256d mask = _mm256_cmp_pd(r2, h2_vec, _CMP_LT_OQ);
    int mask_bits = _mm256_movemask_pd(mask);

    // Compact results
    for (int i = 0; i < 4; ++i) {
        if (mask_bits & (1 << i)) {
            results[count++] = i;
        }
    }
}
#endif

// Fallback scalar implementation
inline void batch_distance_check_scalar(...) { /* existing code */ }

}  // namespace simd
}  // namespace sph
```

**CMake Integration**:
```cmake
# Detect SIMD capabilities
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-mavx2" COMPILER_SUPPORTS_AVX2)
check_cxx_compiler_flag("-mavx512f" COMPILER_SUPPORTS_AVX512)

if(COMPILER_SUPPORTS_AVX2)
    target_compile_definitions(sph PRIVATE SPH_HAS_AVX2)
    target_compile_options(sph PRIVATE -mavx2 -mfma)
endif()
```

---

## Implementation Plan

### Phase 1: Morton Code Reordering (Priority: High)

| Step | Task | Files | Test |
|------|------|-------|------|
| 1.1 | Create morton.hpp with encode function | `include/morton.hpp` | Unit test Morton codes |
| 1.2 | Implement sort_by_morton() | `src/morton.cpp` | Test sorting preserves particles |
| 1.3 | Add reorder_particles() | `src/morton.cpp` | Test reorder/inverse |
| 1.4 | Integrate into tree build | `src/bhtree.cpp` | Existing neighbor tests pass |
| 1.5 | Benchmark cache performance | N/A | Measure L1/L2 cache misses |

### Phase 2: Iterative Traversal (Priority: Medium)

| Step | Task | Files | Test |
|------|------|-------|------|
| 2.1 | Add neighbor_search_iterative() | `src/bhtree.cpp` | Compare with recursive |
| 2.2 | Add calc_force_iterative() | `src/bhtree.cpp` | Compare gravity results |
| 2.3 | Add compile-time switch | `include/defines.hpp` | Both paths work |
| 2.4 | Benchmark recursion overhead | N/A | Measure stack depth |
| 2.5 | Make iterative default | `include/defines.hpp` | All tests pass |

### Phase 3: SIMD Distance Calculations (Priority: Medium)

| Step | Task | Files | Test |
|------|------|-------|------|
| 3.1 | Create simd_distance.hpp | `include/simd_distance.hpp` | Unit test distance calc |
| 3.2 | Add CMake SIMD detection | `CMakeLists.txt` | Build on different CPUs |
| 3.3 | Integrate into leaf checks | `src/bhtree.cpp` | Existing tests pass |
| 3.4 | Add aligned batch buffers | `src/bhtree.cpp` | Valgrind alignment check |
| 3.5 | Benchmark SIMD speedup | N/A | Compare scalar vs SIMD |

### Phase 4: SoA Particle Layout (Future/Optional)

| Step | Task | Files | Test |
|------|------|-------|------|
| 4.1 | Design ParticleSoA struct | `include/particle_soa.hpp` | N/A |
| 4.2 | Add conversion functions | `src/particle_soa.cpp` | Round-trip test |
| 4.3 | Update hot loops to use SoA | Multiple | All physics tests pass |

---

## Progress Checklist

### Phase 0: Remove Exhaustive Search (COMPLETED)
- [x] Create `test_neighbor_search.cpp` with validation tests
- [x] Run tests to confirm both methods agree (9/9 tests pass)
- [x] Comment out EXHAUSTIVE_SEARCH in defines.hpp
- [x] Verify all tests still pass
- [x] Run physics simulation for validation

### Phase 1: Morton Code Reordering
- [ ] Create `include/morton.hpp` with Morton code encoding
- [ ] Create `src/morton.cpp` with sorting implementation
- [ ] Add `tests/test_morton.cpp` unit tests
- [ ] Integrate Morton reordering into `BHTree::make()`
- [ ] Update CMakeLists.txt with new files
- [ ] Verify all existing tests still pass
- [ ] Benchmark: measure cache miss reduction
- [ ] Update this documentation with results

### Phase 2: Iterative Tree Traversal
- [ ] Add `neighbor_search_iterative()` to bhtree.cpp
- [ ] Add `calc_force_iterative()` to bhtree.cpp
- [ ] Add `#define USE_ITERATIVE_TRAVERSAL` switch
- [ ] Create comparison tests (iterative vs recursive)
- [ ] Benchmark function call overhead reduction
- [ ] Make iterative the default
- [ ] Update documentation

### Phase 3: SIMD Batch Distance
- [ ] Create `include/simd_distance.hpp`
- [ ] Add CMake AVX2/AVX-512 detection
- [ ] Implement AVX2 batch distance check
- [ ] Add scalar fallback for non-AVX systems
- [ ] Integrate into `neighbor_search_iterative()` leaf processing
- [ ] Add alignment to batch buffers (32-byte for AVX2)
- [ ] Benchmark SIMD vs scalar speedup
- [ ] Update documentation

### Validation (Each Phase)
- [ ] All existing unit tests pass
- [ ] Neighbor search validation tests pass (tree vs exhaustive)
- [ ] Physics simulation produces expected results (Sod shock tube)
- [ ] No memory leaks (AddressSanitizer)
- [ ] No undefined behavior (UBSanitizer)

---

## Expected Performance Gains

| Optimization | Mechanism | Expected Speedup |
|-------------|-----------|------------------|
| Morton reordering | Cache locality | 2-3x |
| Iterative traversal | Stack elimination | 10-20% |
| SIMD distance | Vectorization | 2-4x |
| **Combined** | All optimizations | **4-10x** |

### Target Performance

| Particle Count | Current (est.) | Target | Speedup |
|----------------|----------------|--------|---------|
| 100,000 | ~5s/step | <1s/step | 5x |
| 1,000,000 | ~60s/step | <10s/step | 6x |
| 10,000,000 | ~700s/step | <100s/step | 7x |

---

## Risk Mitigation

1. **Maintain backwards compatibility** - All optimizations behind compile-time switches
2. **TDD validation** - Each optimization validated against exhaustive search
3. **Incremental rollout** - One optimization at a time, fully tested
4. **Portable fallbacks** - SIMD has scalar fallback, Morton has unsorted fallback
5. **Benchmark at each step** - Verify actual speedup before moving on

---

## Files to Modify/Create

| File | Action | Description |
|------|--------|-------------|
| `include/morton.hpp` | Create | Morton code computation |
| `src/morton.cpp` | Create | Morton sorting/reordering |
| `include/simd_distance.hpp` | Create | SIMD distance utilities |
| `src/bhtree.cpp` | Modify | Add iterative traversal, SIMD integration |
| `include/bhtree.hpp` | Modify | Add iterative function declarations |
| `include/defines.hpp` | Modify | Add optimization switches |
| `CMakeLists.txt` | Modify | Add new files, SIMD detection |
| `tests/test_morton.cpp` | Create | Morton code unit tests |
| `docs/TREE_OPTIMIZATION_IMPLEMENTATION.md` | Update | This document |

---

## Test Specifications (Phase 0 - Completed)

### Test 1: NeighborCountComparison
```
GIVEN: N particles in random distribution
WHEN: Finding neighbors with both exhaustive and tree search
THEN: Neighbor counts match exactly for each particle
```

### Test 2: NeighborIDSetEquality
```
GIVEN: N particles in random distribution
WHEN: Finding neighbors with both methods
THEN: Set of neighbor IDs is identical (order may differ)
```

### Test 3: DistanceSorting
```
GIVEN: N particles
WHEN: Neighbors found and sorted by distance
THEN: Both methods produce same distance ordering
```

### Test 4: PeriodicBoundaryCorrectness
```
GIVEN: Particles near periodic boundaries
WHEN: Finding neighbors
THEN: Neighbors across boundary are found correctly
```

### Test 5: SymmetricSearchMode
```
GIVEN: Particles with varying smoothing lengths
WHEN: is_ij=true (use max of h_i, h_kernel)
THEN: Both methods find same neighbors
```

### Test 6: EdgeCases
```
GIVEN: Particles at exact same position, or exactly at boundary
WHEN: Finding neighbors
THEN: No crashes, consistent results
```

---

## Rollback Plan

If issues are found after an optimization is enabled:
1. Disable the optimization switch in `defines.hpp`
2. All code paths remain functional
3. Investigate failing test cases
4. For exhaustive search: Uncomment `#define EXHAUSTIVE_SEARCH` in `defines.hpp`

---

## Implementation Progress Log

### 2024-XX-XX: Phase 0 Complete
- Removed O(N²) exhaustive search
- Tree-based O(N log N) neighbor search enabled by default
- All 9 validation tests passing
- ~2x speedup observed at N=10K

### Next: Phase 1 - Morton Code Reordering
- Implementation pending
