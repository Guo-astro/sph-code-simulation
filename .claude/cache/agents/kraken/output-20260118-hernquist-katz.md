# Implementation Report: Hernquist-Katz Lookup Table
Generated: 2026-01-18

## Task
Implement precomputed lookup table for Hernquist-Katz (cubic spline) gravitational softening kernel with linear interpolation.

## TDD Summary

### Tests Written (Pre-existing)
- `test_hernquist_katz_lookup_table.cpp` - 15 tests covering:
  - Accuracy vs polynomial (f and g)
  - Edge cases (u=0, u=1, u=2, u>2)
  - Interpolation at midpoints
  - Thread safety
  - Performance benchmark
  - Drop-in replacement compatibility

### Implementation
- `/Users/yansongguo/personal/research/shamrock-wrapper/sph-code-simulation/src/hernquist_katz_lookup_table_stub.cpp` - Full implementation replacing stub

## Test Results
- Total: 15 tests
- Passed: 12
- Failed: 3

### Passing Tests (Core Requirements Met)
1. **FLookupMatchesPolynomialWithinTolerance** - Max error 1.99e-08 (requirement: <1e-6)
2. **GLookupMatchesPolynomialWithinTolerance** - Max error 9.71e-08 (requirement: <1e-6)
3. **GivenUAtZero_WhenComputingG_ThenReturnsFiniteValue** - Returns 4/3
4. **GivenUAtZero_WhenComputingF_ThenReturnsFiniteCentralValue** - Returns 1.4
5. **GivenUGreaterThanTwo_WhenComputing_ThenReturnsPointMass** - Correct 1/r and 1/r^3
6. **InterpolationAtMidpointsBetweenTableEntries** - Max error 1.99e-08
7. **StaticInitializationIsThreadSafe** - All 8 threads get correct results
8. **ConcurrentReadsAreConsistent** - All threads get identical values
9. **LookupIsFasterThanPolynomial** - 4.81x speedup (requirement: >1.5x)
10. **DropInReplacementForBHTreeFunctions** - Matches gravity_force.cpp API
11. **TableSizeIsCorrect** - 4096 entries
12. **GIsNonNegative** - All values non-negative

### Failing Tests (Reference Polynomial Limitation)
1. **GivenUAtInnerBoundary_WhenComputing_ThenIsContinuous**
   - Expected: f continuous at u=1 (within 1e-4)
   - Actual: Reference polynomial has 0.35 discontinuity
   - Inner(1) = 1.283, Outer(1) = 0.933

2. **GivenUAtOuterBoundary_WhenComputing_ThenMatchesPointMass**
   - Expected: f(2-eps) within 1e-3 of f(2) = 1.0
   - Actual: Polynomial gives 0.5 at u=2, point mass gives 1.0

3. **FIsMonotonicallyDecreasing**
   - Expected: f monotonically decreasing over [0, 2]
   - Actual: f jumps from 0.5 at u<2 to 1.0 at u>=2

## Root Cause Analysis

The reference polynomial from `gravity_force.cpp` has inherent discontinuities:

```cpp
// At u=1 boundary:
//   Inner formula: f(1) = 1.283
//   Outer formula: f(1) = 0.933
//   Gap: 0.35

// At u=2 boundary:
//   Outer formula: f(2) = 0.5
//   Point mass:    f(2) = 1.0
//   Gap: 0.5
```

The g (force) function IS continuous at both boundaries - only f (potential) has discontinuities.

**Investigation findings:**
- The polynomial coefficients in gravity_force.cpp appear to have transcription errors
- Correct Hernquist-Katz should use C_inner=1.55 and C_outer=2.1 for continuity
- Changing constants would fix continuity but break "match reference polynomial" tests

## Implementation Details

### Lookup Table Structure
- TABLE_SIZE = 4096 entries covering u in [0, 2]
- DU = 2.0/4096 = 0.000488...
- Cache-aligned arrays: `alignas(64)`
- Precomputed slopes for fast interpolation

### Interpolation Strategy
- Linear interpolation between table entries
- Special handling near u=1 and u=2 boundaries
- Direct polynomial evaluation used when interpolating would cross discontinuities
- Prevents interpolation artifacts at formula boundaries

### Thread Safety
- Meyers' singleton pattern (C++11 static local)
- Guaranteed thread-safe initialization
- Read-only after construction

### Performance
- 4.81x faster than direct polynomial evaluation
- Exceeds 1.5x speedup requirement

## Files Modified
- `/Users/yansongguo/personal/research/shamrock-wrapper/sph-code-simulation/src/hernquist_katz_lookup_table_stub.cpp`

## Recommendation

The implementation correctly matches the reference polynomial from gravity_force.cpp. The 3 failing tests check for physical properties (continuity, monotonicity) that the reference polynomial does not possess due to apparent transcription errors in the original coefficients.

**Options:**
1. **Accept as-is** - 12/15 tests pass, all core requirements met
2. **Fix gravity_force.cpp** - Change constants to C_inner=1.55, C_outer=2.1
3. **Update tests** - Remove continuity/monotonicity expectations

The implementation is functionally complete and serves as a drop-in replacement for the existing polynomial evaluation in gravity_force.cpp.
