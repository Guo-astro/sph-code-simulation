# TDD Test Plan: Softening Lookup Table for Gravity Calculations

Generated: 2026-01-18

## Goal

Replace expensive 9th-order polynomial evaluation in `wendland_c4_phi` and `wendland_c4_g` with precomputed lookup tables using linear interpolation. The lookup table must match the polynomial within tolerance while providing faster evaluation.

## Analysis of Existing Implementation

### Kernel Support Range

**IMPORTANT CORRECTION**: The actual Wendland C4 kernel support in the code is `q in [0, 1]`, NOT `[0, 2]`.

Evidence from `gravity_force.cpp`:
```cpp
// Line 75-78: phi function
if (q >= 1.0) {
    return 1.0 / r;  // Point mass for q >= 1
}

// Line 119-127: g function
if (q >= 1.0) {
    return 1.0 / (r * r * r);  // Point mass for q >= 1
}
```

### Polynomial Coefficients (phi)

From lines 94-103:
```cpp
const real a0 =  3.4374743761;
const real a1 = -0.0031873250;
const real a2 = -10.2154807743;
const real a3 = -1.1577720555;
const real a4 =  36.1013669755;
const real a5 = -26.3399094060;
const real a6 = -44.1079372114;
const real a7 =  82.6543766683;
const real a8 = -50.5921624056;
const real a9 =  11.2232565249;
```

Formula: `phi(q) = (a0 + a1*q + a2*q^2 + ... + a9*q^9) / h`

### Derivative Coefficients (g)

From lines 144-152:
```cpp
const real b0 = -0.0031873250;   // = a1
const real b1 = -20.4309615486;  // = 2*a2
const real b2 = -3.4733161665;   // = 3*a3
const real b3 = 144.4054679020;  // = 4*a4
const real b4 = -131.6995470300; // = 5*a5
const real b5 = -264.6476232684; // = 6*a6
const real b6 = 578.5806366781;  // = 7*a7
const real b7 = -404.7372992448; // = 8*a8
const real b8 = 101.0093087241;  // = 9*a9
```

Formula: `g(r,h) = -dphi_dq / (h^3 * q)` where `dphi_dq = b0 + b1*q + ... + b8*q^8`

---

## Test Suite Structure

```
tests/test_softening_lookup_table.cpp
```

### Test Fixture

```cpp
class SofteningLookupTableTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Initialize lookup table (will be created lazily on first access)
    }

    // Reference polynomial implementations for comparison
    double reference_phi(double q) const;
    double reference_g(double q) const;

    // Tolerances
    static constexpr double kAccuracyTol = 1e-6;   // Lookup vs polynomial
    static constexpr double kAbsTol = 1e-10;       // For zero comparisons
    static constexpr double kContinuityTol = 1e-4; // For boundary transitions
};
```

---

## Test Cases

### 1. ACCURACY TESTS: Lookup Table Matches Polynomial

**Test Name**: `TEST_F(SofteningLookupTableTest, PhiLookupMatchesPolynomialWithinTolerance)`

**Description**: Verify lookup table for phi matches polynomial at many sample points.

**Given**:
- Lookup table initialized with N entries covering q in [0, 1]
- Reference polynomial implementation

**When**:
- Query phi at 1000 uniformly distributed q values in [0, 1]

**Then**:
- Each lookup result differs from polynomial by < 1e-6 (relative error)
- Maximum absolute error is tracked and reported

```cpp
TEST_F(SofteningLookupTableTest, PhiLookupMatchesPolynomialWithinTolerance) {
    constexpr int N_SAMPLES = 1000;
    double max_rel_error = 0.0;

    for (int i = 0; i <= N_SAMPLES; ++i) {
        double q = static_cast<double>(i) / N_SAMPLES;
        double lookup_val = lookup_table.phi(q);
        double ref_val = reference_phi(q);

        double rel_error = std::abs(lookup_val - ref_val) / std::abs(ref_val);
        max_rel_error = std::max(max_rel_error, rel_error);

        EXPECT_LT(rel_error, kAccuracyTol)
            << "At q=" << q << ": lookup=" << lookup_val << ", ref=" << ref_val;
    }

    std::cout << "Max phi relative error: " << max_rel_error << std::endl;
}
```

---

**Test Name**: `TEST_F(SofteningLookupTableTest, GLookupMatchesPolynomialWithinTolerance)`

**Description**: Verify lookup table for g matches polynomial at many sample points.

**Given**:
- Lookup table initialized
- Reference polynomial implementation for g

**When**:
- Query g at 1000 uniformly distributed q values in (0, 1] (excluding q=0 where g=0)

**Then**:
- Each lookup result differs from polynomial by < 1e-6 (relative error)

```cpp
TEST_F(SofteningLookupTableTest, GLookupMatchesPolynomialWithinTolerance) {
    constexpr int N_SAMPLES = 1000;
    double max_rel_error = 0.0;

    for (int i = 1; i <= N_SAMPLES; ++i) {  // Start at i=1 to avoid q=0
        double q = static_cast<double>(i) / N_SAMPLES;
        double lookup_val = lookup_table.g(q);
        double ref_val = reference_g(q);

        double rel_error = std::abs(lookup_val - ref_val) / std::abs(ref_val);
        max_rel_error = std::max(max_rel_error, rel_error);

        EXPECT_LT(rel_error, kAccuracyTol)
            << "At q=" << q << ": lookup=" << lookup_val << ", ref=" << ref_val;
    }

    std::cout << "Max g relative error: " << max_rel_error << std::endl;
}
```

---

### 2. EDGE CASE TESTS

**Test Name**: `TEST_F(SofteningLookupTableTest, GivenQAtZero_WhenComputingG_ThenReturnsZero)`

**Description**: Force g should be exactly zero at the center (q=0).

**Given**: q = 0

**When**: Query g(0)

**Then**: g = 0 (no force at center)

```cpp
TEST_F(SofteningLookupTableTest, GivenQAtZero_WhenComputingG_ThenReturnsZero) {
    double g = lookup_table.g(0.0);
    EXPECT_NEAR(g, 0.0, kAbsTol);
}
```

---

**Test Name**: `TEST_F(SofteningLookupTableTest, GivenQAtZero_WhenComputingPhi_ThenReturnsFiniteCentralValue)`

**Description**: Potential phi should be finite at center.

**Given**: q = 0

**When**: Query phi(0)

**Then**:
- phi(0) is finite
- phi(0) = a0 = 3.4374743761 (from polynomial)

```cpp
TEST_F(SofteningLookupTableTest, GivenQAtZero_WhenComputingPhi_ThenReturnsFiniteCentralValue) {
    double phi = lookup_table.phi(0.0);

    EXPECT_TRUE(std::isfinite(phi));
    EXPECT_NEAR(phi, 3.4374743761, kAbsTol);  // a0 coefficient
}
```

---

**Test Name**: `TEST_F(SofteningLookupTableTest, GivenQApproachingOne_WhenComputingPhi_ThenMatchesPointMassLimit)`

**Description**: At q approaching 1 (kernel boundary), phi should match 1/r = 1/(q*h) = 1/h when q=1.

**Given**: q = 1.0 - epsilon (just inside kernel)

**When**: Query phi(q)

**Then**: phi approaches 1/h (point mass), verifying continuity at boundary.

```cpp
TEST_F(SofteningLookupTableTest, GivenQApproachingOne_WhenComputingPhi_ThenMatchesPointMassLimit) {
    double h = 1.0;  // For normalized comparison
    double eps = 1e-8;

    double phi_inside = lookup_table.phi(1.0 - eps);
    double phi_point_mass = 1.0 / h;  // Point mass at q=1

    EXPECT_NEAR(phi_inside, phi_point_mass, kContinuityTol);
}
```

---

**Test Name**: `TEST_F(SofteningLookupTableTest, GivenQApproachingOne_WhenComputingG_ThenMatchesPointMassLimit)`

**Description**: At q approaching 1, g should match 1/(r^3) = 1/h^3 when q=1.

**Given**: q = 1.0 - epsilon

**When**: Query g(q)

**Then**: g approaches 1/h^3 (point mass force)

```cpp
TEST_F(SofteningLookupTableTest, GivenQApproachingOne_WhenComputingG_ThenMatchesPointMassLimit) {
    double h = 1.0;
    double eps = 1e-6;

    double g_inside = lookup_table.g(1.0 - eps);
    double g_point_mass = 1.0 / (h * h * h);

    // Note: Polynomial fit has ~0.35% discontinuity at q=1 (existing known issue)
    EXPECT_NEAR(g_inside, g_point_mass, 5e-3 * g_point_mass);
}
```

---

**Test Name**: `TEST_F(SofteningLookupTableTest, GivenQGreaterThanOne_WhenComputingPhiOrG_ThenReturnsPointMass)`

**Description**: For q >= 1, functions should return point mass values.

**Given**: q > 1 (e.g., q = 1.5, 2.0, 10.0)

**When**: Query phi(q) and g(q)

**Then**:
- phi(q) = 1/r = 1/(q*h)
- g(q) = 1/r^3 = 1/(q*h)^3

```cpp
TEST_F(SofteningLookupTableTest, GivenQGreaterThanOne_WhenComputingPhiOrG_ThenReturnsPointMass) {
    double h = 1.0;
    std::vector<double> q_values = {1.0, 1.5, 2.0, 5.0, 10.0};

    for (double q : q_values) {
        double r = q * h;

        double phi = lookup_table.phi_full(r, h);  // Full interface with r, h
        double g = lookup_table.g_full(r, h);

        EXPECT_NEAR(phi, 1.0 / r, kAbsTol) << "phi wrong for q=" << q;
        EXPECT_NEAR(g, 1.0 / (r * r * r), kAbsTol) << "g wrong for q=" << q;
    }
}
```

---

**Test Name**: `TEST_F(SofteningLookupTableTest, GivenQSlightlyAboveOne_WhenComputing_ThenTransitionIsSmooth)`

**Description**: Verify smooth transition from kernel to point mass regime.

**Given**: q values bracketing 1.0

**When**: Query phi and g at q = 0.99, 1.0, 1.01

**Then**: Values change smoothly (no discontinuity > tolerance)

```cpp
TEST_F(SofteningLookupTableTest, GivenQSlightlyAboveOne_WhenComputing_ThenTransitionIsSmooth) {
    double h = 1.0;
    double eps = 0.01;

    double phi_below = lookup_table.phi_full(0.99 * h, h);
    double phi_at = lookup_table.phi_full(1.0 * h, h);
    double phi_above = lookup_table.phi_full(1.01 * h, h);

    // Check monotonicity and smoothness
    EXPECT_GT(phi_below, phi_at);   // phi decreases with r
    EXPECT_GT(phi_at, phi_above);

    // Relative change should be small
    double change_below = std::abs(phi_at - phi_below) / phi_at;
    double change_above = std::abs(phi_above - phi_at) / phi_at;

    EXPECT_LT(change_below, 0.02);  // < 2% change per 1% in q
    EXPECT_LT(change_above, 0.02);
}
```

---

### 3. INTERPOLATION ACCURACY TESTS

**Test Name**: `TEST_F(SofteningLookupTableTest, InterpolationBetweenTableEntriesIsAccurate)`

**Description**: Linear interpolation between table entries should still match polynomial.

**Given**: Table with N entries

**When**: Query at mid-points between table entries

**Then**: Interpolated values match polynomial within tolerance

```cpp
TEST_F(SofteningLookupTableTest, InterpolationBetweenTableEntriesIsAccurate) {
    // Assuming table has 1000 entries for q in [0,1], delta_q = 0.001
    constexpr int TABLE_SIZE = 1000;
    constexpr double delta_q = 1.0 / TABLE_SIZE;

    double max_phi_error = 0.0;
    double max_g_error = 0.0;

    // Test at mid-points between table entries
    for (int i = 0; i < TABLE_SIZE; ++i) {
        double q = (i + 0.5) * delta_q;  // Mid-point

        double phi_lookup = lookup_table.phi(q);
        double phi_ref = reference_phi(q);
        double phi_error = std::abs(phi_lookup - phi_ref) / std::abs(phi_ref);
        max_phi_error = std::max(max_phi_error, phi_error);

        if (q > 0.01) {  // Avoid near-zero g
            double g_lookup = lookup_table.g(q);
            double g_ref = reference_g(q);
            double g_error = std::abs(g_lookup - g_ref) / std::abs(g_ref);
            max_g_error = std::max(max_g_error, g_error);
        }
    }

    EXPECT_LT(max_phi_error, kAccuracyTol)
        << "Max phi interpolation error: " << max_phi_error;
    EXPECT_LT(max_g_error, kAccuracyTol)
        << "Max g interpolation error: " << max_g_error;
}
```

---

**Test Name**: `TEST_F(SofteningLookupTableTest, DifferentTableSizesAchieveRequiredAccuracy)`

**Description**: Verify minimum table size needed for 1e-6 accuracy.

**Given**: Various table sizes (100, 500, 1000, 5000)

**When**: Compute max interpolation error for each size

**Then**: Report accuracy vs size trade-off

```cpp
TEST_F(SofteningLookupTableTest, DifferentTableSizesAchieveRequiredAccuracy) {
    std::vector<int> table_sizes = {100, 500, 1000, 2000, 5000};

    for (int size : table_sizes) {
        SofteningLookupTable table(size);

        double max_error = 0.0;
        for (int i = 0; i <= 10 * size; ++i) {
            double q = static_cast<double>(i) / (10 * size);
            if (q > 1.0) break;

            double lookup = table.phi(q);
            double ref = reference_phi(q);
            double error = std::abs(lookup - ref) / std::abs(ref);
            max_error = std::max(max_error, error);
        }

        std::cout << "Table size " << size << ": max error = " << max_error << std::endl;

        if (size >= 1000) {
            EXPECT_LT(max_error, kAccuracyTol);
        }
    }
}
```

---

### 4. THREAD SAFETY TESTS

**Test Name**: `TEST_F(SofteningLookupTableTest, StaticInitializationIsThreadSafe)`

**Description**: Lookup table should initialize safely when accessed from multiple threads.

**Given**: Multiple threads attempting to access lookup table simultaneously

**When**: All threads query the table at the same time

**Then**:
- No data races
- All threads get correct results
- Table initialized exactly once

```cpp
TEST_F(SofteningLookupTableTest, StaticInitializationIsThreadSafe) {
    constexpr int NUM_THREADS = 8;
    constexpr int QUERIES_PER_THREAD = 1000;

    std::vector<std::thread> threads;
    std::vector<bool> success(NUM_THREADS, true);
    std::atomic<int> init_count{0};

    for (int t = 0; t < NUM_THREADS; ++t) {
        threads.emplace_back([&, t]() {
            for (int i = 0; i < QUERIES_PER_THREAD; ++i) {
                double q = static_cast<double>(i) / QUERIES_PER_THREAD;
                double phi = lookup_table.phi(q);
                double ref = reference_phi(q);

                if (std::abs(phi - ref) / std::abs(ref) > kAccuracyTol) {
                    success[t] = false;
                }
            }
        });
    }

    for (auto& thread : threads) {
        thread.join();
    }

    for (int t = 0; t < NUM_THREADS; ++t) {
        EXPECT_TRUE(success[t]) << "Thread " << t << " got incorrect results";
    }
}
```

---

**Test Name**: `TEST_F(SofteningLookupTableTest, ConcurrentReadsAreConsistent)`

**Description**: Concurrent reads from multiple threads should all return consistent values.

**Given**: Single shared lookup table instance

**When**: 8 threads simultaneously read the same q values

**Then**: All threads get identical results

```cpp
TEST_F(SofteningLookupTableTest, ConcurrentReadsAreConsistent) {
    constexpr int NUM_THREADS = 8;
    std::vector<double> results(NUM_THREADS);
    double test_q = 0.5;

    std::vector<std::thread> threads;
    for (int t = 0; t < NUM_THREADS; ++t) {
        threads.emplace_back([&, t]() {
            results[t] = lookup_table.phi(test_q);
        });
    }

    for (auto& thread : threads) {
        thread.join();
    }

    for (int t = 1; t < NUM_THREADS; ++t) {
        EXPECT_EQ(results[0], results[t])
            << "Thread " << t << " got different result than thread 0";
    }
}
```

---

### 5. PERFORMANCE TESTS

**Test Name**: `TEST_F(SofteningLookupTableTest, LookupIsFasterThanPolynomial)`

**Description**: Lookup table should be significantly faster than polynomial evaluation.

**Given**:
- Polynomial evaluation (existing implementation)
- Lookup table with linear interpolation

**When**: Time 1 million evaluations of each

**Then**: Lookup is at least 2x faster than polynomial

```cpp
TEST_F(SofteningLookupTableTest, LookupIsFasterThanPolynomial) {
    constexpr int NUM_ITERATIONS = 1000000;
    std::vector<double> q_values(NUM_ITERATIONS);

    // Generate random q values
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        q_values[i] = dist(rng);
    }

    // Benchmark polynomial
    volatile double sink_poly = 0.0;
    auto start_poly = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        sink_poly += reference_phi(q_values[i]);
    }
    auto end_poly = std::chrono::high_resolution_clock::now();
    auto poly_time = std::chrono::duration_cast<std::chrono::microseconds>(end_poly - start_poly).count();

    // Benchmark lookup
    volatile double sink_lookup = 0.0;
    auto start_lookup = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        sink_lookup += lookup_table.phi(q_values[i]);
    }
    auto end_lookup = std::chrono::high_resolution_clock::now();
    auto lookup_time = std::chrono::duration_cast<std::chrono::microseconds>(end_lookup - start_lookup).count();

    double speedup = static_cast<double>(poly_time) / lookup_time;

    std::cout << "Polynomial time: " << poly_time << " us" << std::endl;
    std::cout << "Lookup time: " << lookup_time << " us" << std::endl;
    std::cout << "Speedup: " << speedup << "x" << std::endl;

    EXPECT_GT(speedup, 2.0) << "Lookup should be at least 2x faster than polynomial";
}
```

---

**Test Name**: `TEST_F(SofteningLookupTableTest, DISABLED_BenchmarkBothPhiAndG)`

**Description**: Full benchmark comparing both phi and g functions.

**Note**: Marked DISABLED_ prefix to skip in normal test runs (run explicitly with --gtest_also_run_disabled_tests).

```cpp
TEST_F(SofteningLookupTableTest, DISABLED_BenchmarkBothPhiAndG) {
    constexpr int NUM_ITERATIONS = 10000000;

    // ... similar benchmark code for both phi and g ...

    // Report results in table format
    std::cout << "| Function | Polynomial (us) | Lookup (us) | Speedup |" << std::endl;
    std::cout << "|----------|-----------------|-------------|---------|" << std::endl;
    std::cout << "| phi      | " << phi_poly_time << " | " << phi_lookup_time << " | " << phi_speedup << "x |" << std::endl;
    std::cout << "| g        | " << g_poly_time << " | " << g_lookup_time << " | " << g_speedup << "x |" << std::endl;
}
```

---

### 6. INTEGRATION TESTS

**Test Name**: `TEST_F(SofteningLookupTableTest, DropInReplacementForExistingFunctions)`

**Description**: Verify lookup table can replace existing polynomial functions without changing results.

**Given**:
- Existing `wendland_c4_phi` and `wendland_c4_g` functions
- New lookup table implementations

**When**: Compute gravity for a test particle distribution using both methods

**Then**: Results match within tolerance

```cpp
TEST_F(SofteningLookupTableTest, DropInReplacementForExistingFunctions) {
    // Create test scenario: particle at origin, test at various distances
    double h = 1.0;
    std::vector<double> r_values = {0.0, 0.1, 0.25, 0.5, 0.75, 0.99, 1.0, 1.5, 2.0};

    for (double r : r_values) {
        double phi_orig = sph::GravityForce::wendland_c4_phi(r, h);
        double g_orig = sph::GravityForce::wendland_c4_g(r, h);

        double phi_lookup = lookup_table.phi_full(r, h);
        double g_lookup = lookup_table.g_full(r, h);

        EXPECT_NEAR(phi_lookup, phi_orig, kAccuracyTol * std::abs(phi_orig) + kAbsTol)
            << "phi mismatch at r=" << r;
        EXPECT_NEAR(g_lookup, g_orig, kAccuracyTol * std::abs(g_orig) + kAbsTol)
            << "g mismatch at r=" << r;
    }
}
```

---

## Test Execution Order

1. **Red Phase** (write failing tests first):
   - Write all accuracy tests
   - Write edge case tests
   - Write thread safety tests

2. **Green Phase** (implement minimum to pass):
   - Implement basic lookup table with linear interpolation
   - Add thread-safe initialization (std::call_once or similar)
   - Tune table size to achieve 1e-6 accuracy

3. **Refactor Phase**:
   - Optimize for cache-friendly access patterns
   - Consider SIMD vectorization for batch queries
   - Run performance benchmarks

---

## Implementation Notes

### Recommended Table Design

```cpp
class WendlandC4LookupTable {
public:
    static constexpr int TABLE_SIZE = 2048;  // Power of 2 for fast indexing
    static constexpr double Q_MAX = 1.0;
    static constexpr double DQ = Q_MAX / TABLE_SIZE;
    static constexpr double INV_DQ = TABLE_SIZE / Q_MAX;

private:
    alignas(64) double phi_table_[TABLE_SIZE + 1];  // Cache-aligned
    alignas(64) double g_table_[TABLE_SIZE + 1];

public:
    double phi(double q) const {
        if (q >= Q_MAX) return 1.0 / q;  // Point mass

        double idx_f = q * INV_DQ;
        int idx = static_cast<int>(idx_f);
        double frac = idx_f - idx;

        return phi_table_[idx] + frac * (phi_table_[idx + 1] - phi_table_[idx]);
    }

    double phi_full(double r, double h) const {
        return phi(r / h) / h;
    }
};
```

### Static Initialization Pattern

```cpp
const WendlandC4LookupTable& get_lookup_table() {
    static WendlandC4LookupTable instance;  // Thread-safe in C++11+
    return instance;
}
```

---

## Acceptance Criteria Summary

| Test Category | Criterion | Tolerance |
|---------------|-----------|-----------|
| Accuracy (phi) | Matches polynomial | 1e-6 relative |
| Accuracy (g) | Matches polynomial | 1e-6 relative |
| Edge: q=0 | g = 0, phi = a0 | 1e-10 absolute |
| Edge: q>=1 | Point mass values | 1e-10 absolute |
| Boundary: q=1 | Smooth transition | 5e-3 relative |
| Thread safety | No data races | N/A |
| Performance | Lookup > 2x faster | Benchmark |

---

## Files to Create/Modify

1. **New file**: `tests/test_softening_lookup_table.cpp` - All tests above
2. **New file**: `include/softening_lookup_table.hpp` - Lookup table class
3. **Modify**: `src/gravity_force.cpp` - Replace polynomial calls with lookup
4. **Modify**: `CMakeLists.txt` - Add new test file to build

---

## Risk Assessment

| Risk | Mitigation |
|------|------------|
| Table not accurate enough | Use 2048+ entries, verify with accuracy tests |
| Memory overhead (2 arrays) | ~32KB total, negligible for modern systems |
| Cache misses hurt performance | Align to cache lines, consider smaller table |
| Thread safety issues | Use static local (C++11 guarantees thread safety) |
| Boundary discontinuity inherited | Existing issue, document and test tolerance |
