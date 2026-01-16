# Tree Optimization Implementation Plan

## Overview

**Branch:** `feature/optimize-tree-remove-exhaustive-search`
**Goal:** Remove O(N²) exhaustive search and enable O(N log N) tree-based neighbor search for 1M+ particle simulations
**Approach:** TDD - Use exhaustive search as ground truth to validate tree implementation before removal

---

## Definition of Done

### Phase 1: Validation Infrastructure ✅ COMPLETE WHEN:
- [ ] Test file `test_neighbor_search.cpp` created
- [ ] Tests compare tree search vs exhaustive search results
- [ ] Tests cover:
  - [ ] Neighbor count matches within tolerance
  - [ ] Neighbor IDs match (order-independent)
  - [ ] Distance ordering is correct
  - [ ] Edge cases: particles at boundaries, overlapping particles
  - [ ] Periodic boundary conditions
  - [ ] Variable smoothing lengths (h_i ≠ h_j)
- [ ] Tests pass with current `EXHAUSTIVE_SEARCH` enabled
- [ ] Tests can run with tree search and validate against exhaustive

### Phase 2: Tree Search Validation ✅ COMPLETE WHEN:
- [ ] All validation tests pass comparing tree vs exhaustive
- [ ] Performance test shows tree search is faster for N > 1000
- [ ] No memory leaks (run with sanitizers)
- [ ] Works in 1D, 2D, and 3D configurations

### Phase 3: Remove EXHAUSTIVE_SEARCH ✅ COMPLETE WHEN:
- [ ] `#define EXHAUSTIVE_SEARCH` commented out in `defines.hpp`
- [ ] All existing unit tests pass
- [ ] Validation tests confirm tree search produces correct results
- [ ] At least one physics simulation runs and produces expected output
- [ ] Code compiles without warnings

### Phase 4: Cleanup ✅ COMPLETE WHEN:
- [ ] Exhaustive search code remains but is disabled by default
- [ ] Documentation updated
- [ ] Performance benchmarks recorded

---

## Implementation Progress

### Phase 1: Validation Infrastructure

| Task | Status | Notes |
|------|--------|-------|
| Create `test_neighbor_search.cpp` | ✅ Complete | 9 test cases |
| Test: basic neighbor count comparison | ✅ Complete | Pass |
| Test: neighbor ID set equality | ✅ Complete | Pass |
| Test: periodic boundary conditions | ✅ Complete | Pass |
| Test: variable smoothing length | ✅ Complete | Pass (tree conservative) |
| Test: is_ij=true mode (symmetric search) | ✅ Complete | Pass |
| CMake integration | ✅ Complete | test_neighbor_search target |

### Phase 2: Tree Search Validation

| Task | Status | Notes |
|------|--------|-------|
| Run validation tests with tree search | ✅ Complete | 9/9 tests pass |
| Fix any discrepancies found | ✅ Complete | None found |
| Memory leak check with sanitizers | ⬜ Not Started | |
| Test all dimensions (1D, 2D, 3D) | 🔄 Tested 2D | Need 1D, 3D |

### Phase 3: Remove EXHAUSTIVE_SEARCH

| Task | Status | Notes |
|------|--------|-------|
| Comment out EXHAUSTIVE_SEARCH define | ✅ Complete | defines.hpp:39-40 |
| Fix missing tree variable in srmhd | ✅ Complete | srmhd_fluid_force.cpp |
| Verify all tests pass | ✅ Complete | |
| Build main sph executable | ✅ Complete | |
| Run physics simulation (e.g., Sod shock tube) | ⬜ Not Started | |
| Verify simulation output correctness | ⬜ Not Started | |

### Phase 4: Cleanup

| Task | Status | Notes |
|------|--------|-------|
| Add compile-time switch documentation | ✅ Complete | In defines.hpp |
| Performance benchmark (N=10K, 100K, 1M) | 🔄 Partial | 2x speedup at N=10K |
| Update README if needed | ⬜ Not Started | |

---

## Test Specifications

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

## Risk Mitigation

1. **Keep exhaustive search code** - Don't delete, just disable with `#ifdef`
2. **Test incrementally** - Validate each dimension separately
3. **Use deterministic tests** - Fixed random seeds for reproducibility
4. **CI integration** - Ensure tests run on all builds

---

## Files Modified

| File | Change Type | Description |
|------|-------------|-------------|
| `include/defines.hpp` | Modified | Comment out EXHAUSTIVE_SEARCH |
| `CMakeLists.txt` | Modified | Add test_neighbor_search target |
| `tests/test_neighbor_search.cpp` | New | Validation test suite (9 tests) |
| `src/srmhd/srmhd_fluid_force.cpp` | Fixed | Add missing tree variable |
| `docs/TREE_OPTIMIZATION_IMPLEMENTATION.md` | New | This tracking document |

---

## Rollback Plan

If issues are found after removing EXHAUSTIVE_SEARCH:
1. Uncomment `#define EXHAUSTIVE_SEARCH` in `defines.hpp`
2. All code paths remain functional
3. Investigate failing test cases

---

## Performance Targets

| Particle Count | Exhaustive (est.) | Tree (target) | Speedup |
|----------------|-------------------|---------------|---------|
| 10,000 | ~1s | <0.1s | 10x |
| 100,000 | ~100s | <1s | 100x |
| 1,000,000 | ~10000s | <10s | 1000x |

---

## Next Steps

1. ✅ Create this tracking document
2. 🔄 Create `test_neighbor_search.cpp` with validation tests
3. ⬜ Run tests to confirm both methods agree
4. ⬜ Comment out EXHAUSTIVE_SEARCH
5. ⬜ Verify all tests still pass
6. ⬜ Run a physics simulation for final validation
