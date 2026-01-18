# TDD Test Plan: Gravity Optimizations for SPH Code

Generated: 2026-01-18

## Overview

This document defines comprehensive TDD test plans for five gravity optimizations, ordered by implementation priority. Each feature includes test cases for accuracy, edge cases, and performance.

---

## Feature 1: Cubic Spline (Hernquist-Katz) Lookup Table

**Priority:** HIGH - Direct reuse of existing Wendland C4 lookup table pattern

### Background

The Hernquist & Katz (1989) softening kernel is already implemented in `bhtree.cpp` (lines 351-377) and `gravity_force.cpp` (lines 26-52). It uses:
- Softening length: `e = h/2`
- Normalized distance: `u = r/e`
- Support: `u in [0, 2]`
- Three piecewise regions: `u < 1`, `1 <= u < 2`, `u >= 2`

The existing `WendlandC4LookupTable` pattern in `softening_lookup_table.hpp` can be directly reused.

### Mathematical Definition

**Potential f(r, h):**
```cpp
// u < 1
f = (-0.5 * u^2 * (1/3 - 3/20 * u^2 + u^3/20) + 1.4) / e

// 1 <= u < 2
f = -1/(15r) + (-u^2 * (4/3 - u + 0.3*u^2 - u^3/30) + 1.6) / e

// u >= 2 (point mass)
f = 1/r
```

**Force g(r, h):**
```cpp
// u < 1
g = (4/3 - 1.2*u^2 + 0.5*u^3) / e^3

// 1 <= u < 2
g = (-1/15 + 8/3*u^3 - 3*u^4 + 1.2*u^5 - u^6/6) / r^3

// u >= 2 (point mass)
g = 1/r^3
```

### Test File

```
tests/test_hernquist_katz_lookup_table.cpp
```

### Test Fixture

```cpp
class HernquistKatzLookupTableTest : public ::testing::Test {
protected:
    // Reference polynomial implementations (from bhtree.cpp)
    double reference_f(double u) const;
    double reference_g(double u) const;

    // Tolerances
    static constexpr double kAccuracyTol = 1e-6;   // Lookup vs polynomial
    static constexpr double kAbsTol = 1e-10;       // For zero comparisons
    static constexpr double kContinuityTol = 1e-4; // For boundary transitions

    const HernquistKatzLookupTable& lookup_table = HernquistKatzLookupTable::get_instance();
};
```

### Test Cases

#### 1.1 Accuracy Tests

```cpp
TEST_F(HernquistKatzLookupTableTest, FLookupMatchesPolynomialWithinTolerance) {
    // GIVEN: Lookup table covering u in [0, 2]
    // WHEN: Query f at 2000 uniformly distributed u values
    // THEN: Each lookup result differs from polynomial by < 1e-6

    constexpr int N_SAMPLES = 2000;
    double max_rel_error = 0.0;

    for (int i = 0; i <= N_SAMPLES; ++i) {
        double u = 2.0 * static_cast<double>(i) / N_SAMPLES;
        double lookup_val = lookup_table.f(u);
        double ref_val = reference_f(u);

        double rel_error = std::abs(lookup_val - ref_val) / (std::abs(ref_val) + kAbsTol);
        max_rel_error = std::max(max_rel_error, rel_error);

        EXPECT_LT(rel_error, kAccuracyTol)
            << "At u=" << u << ": lookup=" << lookup_val << ", ref=" << ref_val;
    }
}

TEST_F(HernquistKatzLookupTableTest, GLookupMatchesPolynomialWithinTolerance) {
    // GIVEN: Lookup table covering u in [0, 2]
    // WHEN: Query g at 2000 uniformly distributed u values (excluding u=0)
    // THEN: Each lookup result differs from polynomial by < 1e-6

    constexpr int N_SAMPLES = 2000;
    double max_rel_error = 0.0;

    for (int i = 1; i <= N_SAMPLES; ++i) {
        double u = 2.0 * static_cast<double>(i) / N_SAMPLES;
        double lookup_val = lookup_table.g(u);
        double ref_val = reference_g(u);

        double rel_error = std::abs(lookup_val - ref_val) / (std::abs(ref_val) + kAbsTol);
        max_rel_error = std::max(max_rel_error, rel_error);

        EXPECT_LT(rel_error, kAccuracyTol)
            << "At u=" << u << ": lookup=" << lookup_val << ", ref=" << ref_val;
    }
}
```

#### 1.2 Edge Case Tests

```cpp
TEST_F(HernquistKatzLookupTableTest, GivenUAtZero_WhenComputingG_ThenReturnsFiniteValue) {
    // GIVEN: u = 0 (at source position)
    // WHEN: Query g(0)
    // THEN: g = 4/3 / e^3 (finite, not singular)

    double g = lookup_table.g(0.0);
    double expected = 4.0 / 3.0;  // With e = 1 for normalized lookup

    EXPECT_TRUE(std::isfinite(g));
    EXPECT_NEAR(g, expected, kAbsTol);
}

TEST_F(HernquistKatzLookupTableTest, GivenUAtZero_WhenComputingF_ThenReturnsFiniteCentralValue) {
    // GIVEN: u = 0
    // WHEN: Query f(0)
    // THEN: f(0) = 1.4 / e (finite central potential)

    double f = lookup_table.f(0.0);
    EXPECT_TRUE(std::isfinite(f));
    EXPECT_NEAR(f, 1.4, kAbsTol);  // With e = 1
}

TEST_F(HernquistKatzLookupTableTest, GivenUAtInnerBoundary_WhenComputing_ThenIsContinuous) {
    // GIVEN: u = 1.0 (boundary between inner and outer regions)
    // WHEN: Query f and g at u = 1 - eps and u = 1 + eps
    // THEN: Values are continuous across boundary

    double eps = 1e-8;

    double f_below = lookup_table.f(1.0 - eps);
    double f_above = lookup_table.f(1.0 + eps);
    EXPECT_NEAR(f_below, f_above, kContinuityTol);

    double g_below = lookup_table.g(1.0 - eps);
    double g_above = lookup_table.g(1.0 + eps);
    EXPECT_NEAR(g_below, g_above, kContinuityTol);
}

TEST_F(HernquistKatzLookupTableTest, GivenUAtOuterBoundary_WhenComputing_ThenMatchesPointMass) {
    // GIVEN: u = 2.0 (boundary to point mass regime)
    // WHEN: Query f and g at u = 2 - eps and u = 2 + eps
    // THEN: Transition to point mass is smooth

    double eps = 1e-6;
    double e = 1.0;  // Normalized

    double f_inside = lookup_table.f(2.0 - eps);
    double f_point_mass = 1.0 / (2.0 * e);  // 1/r where r = u*e = 2
    EXPECT_NEAR(f_inside, f_point_mass, 1e-3);

    double g_inside = lookup_table.g(2.0 - eps);
    double g_point_mass = 1.0 / (8.0 * e * e * e);  // 1/r^3 where r = 2
    EXPECT_NEAR(g_inside, g_point_mass, 1e-3);
}

TEST_F(HernquistKatzLookupTableTest, GivenUGreaterThanTwo_WhenComputing_ThenReturnsPointMass) {
    // GIVEN: u > 2 (outside kernel support)
    // WHEN: Query f(u) and g(u)
    // THEN: Returns point mass values

    std::vector<double> u_values = {2.0, 3.0, 5.0, 10.0};

    for (double u : u_values) {
        double f = lookup_table.f_full(u, 1.0);  // e = 0.5 for h = 1.0
        double g = lookup_table.g_full(u, 1.0);

        double r = u * 0.5;  // r = u * e where e = h/2
        EXPECT_NEAR(f, 1.0 / r, kAbsTol) << "f wrong for u=" << u;
        EXPECT_NEAR(g, 1.0 / (r * r * r), kAbsTol) << "g wrong for u=" << u;
    }
}
```

#### 1.3 Interpolation Tests

```cpp
TEST_F(HernquistKatzLookupTableTest, InterpolationAtMidpointsBetweenTableEntries) {
    // GIVEN: Table with N entries for u in [0, 2]
    // WHEN: Query at midpoints between table entries
    // THEN: Interpolated values match polynomial within tolerance

    constexpr int TABLE_SIZE = 4096;  // Expected table size
    constexpr double delta_u = 2.0 / TABLE_SIZE;

    double max_f_error = 0.0;
    double max_g_error = 0.0;

    for (int i = 0; i < TABLE_SIZE; ++i) {
        double u = (i + 0.5) * delta_u;  // Midpoint

        double f_lookup = lookup_table.f(u);
        double f_ref = reference_f(u);
        max_f_error = std::max(max_f_error, std::abs(f_lookup - f_ref) / (std::abs(f_ref) + kAbsTol));

        if (u > 0.01) {
            double g_lookup = lookup_table.g(u);
            double g_ref = reference_g(u);
            max_g_error = std::max(max_g_error, std::abs(g_lookup - g_ref) / (std::abs(g_ref) + kAbsTol));
        }
    }

    EXPECT_LT(max_f_error, kAccuracyTol);
    EXPECT_LT(max_g_error, kAccuracyTol);
}
```

#### 1.4 Performance Tests

```cpp
TEST_F(HernquistKatzLookupTableTest, LookupIsFasterThanPolynomial) {
    // GIVEN: Polynomial evaluation vs lookup table
    // WHEN: Time 1 million evaluations of each
    // THEN: Lookup is at least 1.5x faster

    constexpr int NUM_ITERATIONS = 1000000;
    std::vector<double> u_values(NUM_ITERATIONS);

    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(0.0, 2.0);
    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        u_values[i] = dist(rng);
    }

    // Benchmark polynomial
    volatile double sink = 0.0;
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        sink += reference_g(u_values[i]);
    }
    auto poly_time = std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::high_resolution_clock::now() - start).count();

    // Benchmark lookup
    sink = 0.0;
    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        sink += lookup_table.g(u_values[i]);
    }
    auto lookup_time = std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::high_resolution_clock::now() - start).count();

    double speedup = static_cast<double>(poly_time) / lookup_time;
    EXPECT_GT(speedup, 1.5);
}
```

#### 1.5 Drop-in Replacement Tests

```cpp
TEST_F(HernquistKatzLookupTableTest, DropInReplacementForBHTreeFunctions) {
    // GIVEN: Existing f() and g() functions from bhtree.cpp
    // WHEN: Compute gravity at various r, h combinations
    // THEN: Lookup table results match existing implementation

    std::vector<std::pair<double, double>> test_cases = {
        {0.0, 1.0}, {0.1, 1.0}, {0.5, 1.0}, {0.99, 1.0}, {1.01, 1.0},
        {1.5, 1.0}, {1.99, 1.0}, {2.01, 1.0}, {5.0, 1.0},
        {0.5, 0.5}, {0.5, 2.0}  // Different h values
    };

    for (const auto& [r, h] : test_cases) {
        double f_orig = sph::f(r, h);  // Existing function
        double g_orig = sph::g(r, h);

        double f_lookup = lookup_table.f_full(r, h);
        double g_lookup = lookup_table.g_full(r, h);

        EXPECT_NEAR(f_lookup, f_orig, kAccuracyTol * std::abs(f_orig) + kAbsTol)
            << "f mismatch at r=" << r << ", h=" << h;
        EXPECT_NEAR(g_lookup, g_orig, kAccuracyTol * std::abs(g_orig) + kAbsTol)
            << "g mismatch at r=" << r << ", h=" << h;
    }
}
```

### Expected API

```cpp
class HernquistKatzLookupTable {
public:
    static constexpr int TABLE_SIZE = 4096;      // Larger than Wendland due to 2x support
    static constexpr double U_MAX = 2.0;         // Kernel support [0, 2]
    static constexpr double DU = U_MAX / TABLE_SIZE;
    static constexpr double INV_DU = static_cast<double>(TABLE_SIZE) / U_MAX;

    static const HernquistKatzLookupTable& get_instance();

    // Lookup normalized functions (u = r/e where e = h/2)
    double f(double u) const;
    double g(double u) const;

    // Full interface matching existing API
    double f_full(double r, double h) const;
    double g_full(double r, double h) const;

    int size() const { return TABLE_SIZE; }

private:
    HernquistKatzLookupTable();

    alignas(64) double f_table_[TABLE_SIZE + 1];
    alignas(64) double f_slope_[TABLE_SIZE];
    alignas(64) double g_table_[TABLE_SIZE + 1];
    alignas(64) double g_slope_[TABLE_SIZE];
};
```

### Integration with Existing Code

**Files to modify:**
- `include/softening_lookup_table.hpp` - Add `HernquistKatzLookupTable` class
- `src/softening_lookup_table.cpp` - Add implementation (or new file)
- `src/bhtree.cpp` - Replace inline `f()` and `g()` with lookup
- `src/gravity_force.cpp` - Replace inline `f()` and `g()` with lookup

---

## Feature 2: 3D Analytic Gravity Test (Uniform Sphere)

**Priority:** HIGH - Essential for validating gravity implementation

### Background

A uniform density sphere has an exact analytic gravity solution:
- **Inside** (r < R): `g = (4/3) * pi * G * rho * r` (linear in r)
- **Outside** (r > R): `g = G * M / r^2` (inverse square law)

This provides a ground-truth test for kernel-softened gravity implementations.

### Test File

```
tests/test_analytic_gravity.cpp
```

### Test Fixture

```cpp
class AnalyticGravityTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Physical parameters
        G = 1.0;        // Gravitational constant
        rho = 1.0;      // Uniform density
        R = 1.0;        // Sphere radius
        M = (4.0/3.0) * M_PI * rho * R * R * R;  // Total mass
    }

    double G, rho, R, M;

    // Analytic gravity inside sphere: g = (4/3) * pi * G * rho * r
    double analytic_g_inside(double r) const {
        return (4.0/3.0) * M_PI * G * rho * r;
    }

    // Analytic gravity outside sphere: g = G * M / r^2
    double analytic_g_outside(double r) const {
        return G * M / (r * r);
    }

    // Create uniform sphere particle distribution
    std::vector<SPHParticle> create_uniform_sphere(int N, double R, double rho);

    // Compute softened gravity at position r using particle distribution
    double compute_gravity(const vec_t& r, const std::vector<SPHParticle>& particles,
                          double softening, GravitySofteningType type);
};
```

### Test Cases

#### 2.1 Analytic Reference Tests

```cpp
TEST_F(AnalyticGravityTest, GivenUniformSphere_WhenAtCenter_ThenGravityIsZero) {
    // GIVEN: Uniform density sphere
    // WHEN: Compute gravity at r = 0
    // THEN: g = 0 (by symmetry)

    double g = analytic_g_inside(0.0);
    EXPECT_NEAR(g, 0.0, kAbsTol);
}

TEST_F(AnalyticGravityTest, GivenUniformSphere_WhenInsideSphere_ThenGravityIsLinear) {
    // GIVEN: Uniform density sphere of radius R
    // WHEN: Compute gravity at r < R
    // THEN: g proportional to r

    for (double r = 0.1; r < R; r += 0.1) {
        double g = analytic_g_inside(r);
        double expected = (4.0/3.0) * M_PI * G * rho * r;
        EXPECT_NEAR(g, expected, kAbsTol);
    }
}

TEST_F(AnalyticGravityTest, GivenUniformSphere_WhenOutsideSphere_ThenGravityIsInverseSquare) {
    // GIVEN: Uniform density sphere
    // WHEN: Compute gravity at r > R
    // THEN: g = GM/r^2

    for (double r = R + 0.1; r < 5.0 * R; r += 0.5) {
        double g = analytic_g_outside(r);
        double expected = G * M / (r * r);
        EXPECT_NEAR(g, expected, kAbsTol);
    }
}

TEST_F(AnalyticGravityTest, GivenUniformSphere_WhenAtSurface_ThenGravityIsContinuous) {
    // GIVEN: Uniform density sphere
    // WHEN: Compute gravity just inside and outside surface
    // THEN: g is continuous at r = R

    double eps = 1e-8;
    double g_inside = analytic_g_inside(R - eps);
    double g_outside = analytic_g_outside(R + eps);

    EXPECT_NEAR(g_inside, g_outside, kContinuityTol);
}
```

#### 2.2 Softening Convergence Tests

```cpp
TEST_F(AnalyticGravityTest, GivenSoftenedGravity_WhenSofteningDecreases_ThenConvergesToAnalytic) {
    // GIVEN: Uniform sphere particle distribution
    // WHEN: Compute gravity with decreasing softening length
    // THEN: Results converge to analytic solution

    const int N = 10000;  // Number of particles
    auto particles = create_uniform_sphere(N, R, rho);

    vec_t test_pos = {0.5 * R, 0.0, 0.0};  // Inside sphere
    double analytic = analytic_g_inside(0.5 * R);

    std::vector<double> softenings = {0.5, 0.2, 0.1, 0.05, 0.02};
    double prev_error = 1.0;

    for (double eps : softenings) {
        double computed = compute_gravity(test_pos, particles, eps,
                                          GravitySofteningType::WENDLAND_C4);
        double error = std::abs(computed - analytic) / analytic;

        // Error should decrease as softening decreases
        EXPECT_LT(error, prev_error);
        prev_error = error;

        std::cout << "eps=" << eps << ": error=" << error << std::endl;
    }
}

TEST_F(AnalyticGravityTest, GivenHernquistKatzSoftening_WhenSofteningDecreases_ThenConverges) {
    // GIVEN: Uniform sphere with Hernquist-Katz softening
    // WHEN: Softening epsilon -> 0
    // THEN: Converges to analytic gravity

    const int N = 10000;
    auto particles = create_uniform_sphere(N, R, rho);

    // Test multiple positions
    std::vector<double> radii = {0.0, 0.25 * R, 0.5 * R, 0.75 * R, R, 1.5 * R, 2.0 * R};

    for (double r : radii) {
        vec_t test_pos = {r, 0.0, 0.0};
        double analytic = (r < R) ? analytic_g_inside(r) : analytic_g_outside(r);

        std::vector<double> softenings = {0.2, 0.1, 0.05};
        for (double eps : softenings) {
            double computed = compute_gravity(test_pos, particles, eps,
                                              GravitySofteningType::HERNQUIST_KATZ);
            double error = std::abs(computed - analytic) / (analytic + kAbsTol);

            std::cout << "r=" << r << ", eps=" << eps << ": error=" << error << std::endl;
        }
    }
}
```

#### 2.3 Resolution Convergence Tests

```cpp
TEST_F(AnalyticGravityTest, GivenFixedSoftening_WhenResolutionIncreases_ThenErrorDecreases) {
    // GIVEN: Fixed softening epsilon
    // WHEN: Increase number of particles (N)
    // THEN: Error converges to softening-limited error floor

    double eps = 0.1;  // Fixed softening
    std::vector<int> particle_counts = {100, 500, 1000, 5000, 10000};

    vec_t test_pos = {0.5 * R, 0.0, 0.0};
    double analytic = analytic_g_inside(0.5 * R);

    for (int N : particle_counts) {
        auto particles = create_uniform_sphere(N, R, rho);
        double computed = compute_gravity(test_pos, particles, eps,
                                          GravitySofteningType::WENDLAND_C4);
        double error = std::abs(computed - analytic) / analytic;

        std::cout << "N=" << N << ": error=" << error << std::endl;
    }
}
```

#### 2.4 Comparison Tests (Hernquist-Katz vs Wendland C4)

```cpp
TEST_F(AnalyticGravityTest, GivenSameSoftening_WhenComparingKernels_ThenWendlandIsMoreAccurate) {
    // GIVEN: Same softening length for both kernels
    // WHEN: Compute gravity inside sphere
    // THEN: Wendland C4 should have smaller error (smoother kernel)

    const int N = 10000;
    auto particles = create_uniform_sphere(N, R, rho);
    double eps = 0.1;

    vec_t test_pos = {0.5 * R, 0.0, 0.0};
    double analytic = analytic_g_inside(0.5 * R);

    double g_hk = compute_gravity(test_pos, particles, eps, GravitySofteningType::HERNQUIST_KATZ);
    double g_wc4 = compute_gravity(test_pos, particles, eps, GravitySofteningType::WENDLAND_C4);

    double error_hk = std::abs(g_hk - analytic) / analytic;
    double error_wc4 = std::abs(g_wc4 - analytic) / analytic;

    std::cout << "HK error: " << error_hk << ", WC4 error: " << error_wc4 << std::endl;

    // Both should be reasonably accurate
    EXPECT_LT(error_hk, 0.1);
    EXPECT_LT(error_wc4, 0.1);
}
```

### Expected API

```cpp
namespace sph {

// Helper function for creating uniform sphere test distributions
std::vector<SPHParticle> create_uniform_sphere_distribution(
    int N, double R, double rho, unsigned seed = 42);

// Direct summation gravity computation for testing
struct DirectGravityResult {
    vec_t acceleration;
    real potential;
};

DirectGravityResult compute_direct_gravity(
    const vec_t& pos,
    const std::vector<SPHParticle>& particles,
    real softening,
    GravitySofteningType type);

} // namespace sph
```

---

## Feature 3: SoA (Structure of Arrays) Gravity Arrays

**Priority:** MEDIUM - Performance optimization for cache efficiency

### Background

The current implementation uses AoS (Array of Structures) layout via `SPHParticle` class:
```cpp
// Current AoS: particles[i].pos, particles[i].mass (scattered in memory)
for (int i = 0; i < N; ++i) {
    real m = particles[i].mass;      // Cache miss
    vec_t p = particles[i].pos;      // Another cache miss
}
```

SoA layout places same-type data contiguously for better cache utilization:
```cpp
// SoA: pos_x, pos_y, pos_z, mass are contiguous arrays
for (int i = 0; i < N; ++i) {
    real m = mass[i];          // Cache-friendly sequential access
    real x = pos_x[i];
}
```

### Test File

```
tests/test_soa_gravity.cpp
```

### Test Fixture

```cpp
class SoAGravityTest : public ::testing::Test {
protected:
    // Reference AoS implementation
    vec_t compute_gravity_aos(const SPHParticle& target,
                               const std::vector<SPHParticle>& sources,
                               real softening);

    // New SoA implementation
    vec_t compute_gravity_soa(int target_idx,
                               const GravityDataSoA& data,
                               real softening);

    // Convert between layouts
    GravityDataSoA convert_to_soa(const std::vector<SPHParticle>& particles);
};
```

### Test Cases

#### 3.1 Correctness Tests

```cpp
TEST_F(SoAGravityTest, GivenSameParticles_WhenComputingGravity_ThenSoAMatchesAoS) {
    // GIVEN: Particle distribution in both AoS and SoA layouts
    // WHEN: Compute gravity for each particle
    // THEN: Results are identical (bitwise or within epsilon)

    const int N = 1000;
    auto particles = create_random_particles(N);
    auto soa_data = convert_to_soa(particles);
    real softening = 0.1;

    for (int i = 0; i < N; ++i) {
        vec_t g_aos = compute_gravity_aos(particles[i], particles, softening);
        vec_t g_soa = compute_gravity_soa(i, soa_data, softening);

        for (int d = 0; d < DIM; ++d) {
            EXPECT_NEAR(g_aos[d], g_soa[d], kAbsTol)
                << "Mismatch for particle " << i << " dimension " << d;
        }
    }
}

TEST_F(SoAGravityTest, GivenUniformSphere_WhenUsingSoA_ThenMatchesAnalytic) {
    // GIVEN: Uniform sphere in SoA layout
    // WHEN: Compute gravity at test positions
    // THEN: Matches analytic solution

    const int N = 10000;
    double R = 1.0, rho = 1.0;
    auto particles = create_uniform_sphere(N, R, rho);
    auto soa_data = convert_to_soa(particles);

    // Test at multiple radii
    for (double r : {0.25, 0.5, 0.75, 1.0, 1.5, 2.0}) {
        vec_t test_pos = {r, 0.0, 0.0};
        double analytic = (r < R) ? analytic_g_inside(r) : analytic_g_outside(r);

        double computed = compute_gravity_soa_at_position(test_pos, soa_data, 0.05);
        double error = std::abs(computed - analytic) / (analytic + kAbsTol);

        EXPECT_LT(error, 0.05);
    }
}
```

#### 3.2 Layout Conversion Tests

```cpp
TEST_F(SoAGravityTest, ConversionFromAoSToSoAPreservesData) {
    // GIVEN: AoS particle array
    // WHEN: Convert to SoA
    // THEN: All data preserved exactly

    const int N = 1000;
    auto particles = create_random_particles(N);
    auto soa_data = convert_to_soa(particles);

    for (int i = 0; i < N; ++i) {
        EXPECT_EQ(soa_data.pos_x[i], particles[i].pos[0]);
        EXPECT_EQ(soa_data.pos_y[i], particles[i].pos[1]);
        EXPECT_EQ(soa_data.pos_z[i], particles[i].pos[2]);
        EXPECT_EQ(soa_data.mass[i], particles[i].mass);
        EXPECT_EQ(soa_data.sml[i], particles[i].sml);
    }
}

TEST_F(SoAGravityTest, RoundTripConversionPreservesData) {
    // GIVEN: Original AoS particles
    // WHEN: Convert AoS -> SoA -> AoS
    // THEN: Data unchanged

    auto particles = create_random_particles(1000);
    auto soa = convert_to_soa(particles);
    auto particles2 = convert_to_aos(soa);

    for (int i = 0; i < 1000; ++i) {
        EXPECT_EQ(particles[i].pos, particles2[i].pos);
        EXPECT_EQ(particles[i].mass, particles2[i].mass);
    }
}
```

#### 3.3 Memory Alignment Tests

```cpp
TEST_F(SoAGravityTest, SoAArraysAreCacheAligned) {
    // GIVEN: SoA gravity data structure
    // WHEN: Check memory alignment
    // THEN: Arrays are cache-line aligned (64 bytes)

    GravityDataSoA data(1000);

    EXPECT_EQ(reinterpret_cast<uintptr_t>(data.pos_x.data()) % 64, 0);
    EXPECT_EQ(reinterpret_cast<uintptr_t>(data.pos_y.data()) % 64, 0);
    EXPECT_EQ(reinterpret_cast<uintptr_t>(data.pos_z.data()) % 64, 0);
    EXPECT_EQ(reinterpret_cast<uintptr_t>(data.mass.data()) % 64, 0);
    EXPECT_EQ(reinterpret_cast<uintptr_t>(data.sml.data()) % 64, 0);
}
```

#### 3.4 Performance Tests

```cpp
TEST_F(SoAGravityTest, SoAIsFasterThanAoSForDirectSum) {
    // GIVEN: Same particle distribution in AoS and SoA
    // WHEN: Compute gravity via direct sum
    // THEN: SoA is faster due to better cache utilization

    const int N = 5000;  // O(N^2) so keep small
    auto particles = create_random_particles(N);
    auto soa_data = convert_to_soa(particles);
    real softening = 0.1;

    // Benchmark AoS
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i) {
        volatile vec_t g = compute_gravity_aos(particles[i], particles, softening);
    }
    auto aos_time = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - start).count();

    // Benchmark SoA
    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i) {
        volatile vec_t g = compute_gravity_soa(i, soa_data, softening);
    }
    auto soa_time = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - start).count();

    std::cout << "AoS: " << aos_time << " ms, SoA: " << soa_time << " ms" << std::endl;

    // SoA should be faster (at least 1.2x on most hardware)
    EXPECT_LT(soa_time, aos_time);
}
```

### Expected API

```cpp
namespace sph {

// SoA layout for gravity-critical data
struct GravityDataSoA {
    // Cache-aligned contiguous arrays
    alignas(64) std::vector<real> pos_x;
    alignas(64) std::vector<real> pos_y;
    alignas(64) std::vector<real> pos_z;
    alignas(64) std::vector<real> mass;
    alignas(64) std::vector<real> sml;

    // Optional output arrays
    alignas(64) std::vector<real> grav_acc_x;
    alignas(64) std::vector<real> grav_acc_y;
    alignas(64) std::vector<real> grav_acc_z;
    alignas(64) std::vector<real> phi;

    explicit GravityDataSoA(size_t n);

    // Conversion utilities
    static GravityDataSoA from_aos(const std::vector<SPHParticle>& particles);
    void copy_results_to_aos(std::vector<SPHParticle>& particles) const;
};

// SoA-based gravity computation
void compute_gravity_soa(GravityDataSoA& data, real G,
                         GravitySofteningType type, real softening);

} // namespace sph
```

---

## Feature 4: Batch Tree Queries

**Priority:** MEDIUM - Performance optimization for tree traversal

### Background

Current implementation traverses tree independently for each particle:
```cpp
for (int i = 0; i < N; ++i) {
    tree.tree_force(particles[i]);  // Independent traversal
}
```

Particles in the same tree cell can share node opening decisions:
```cpp
void tree_force_batch(const std::vector<int>& indices) {
    // All particles share same tree traversal path
    // Node opening decisions apply to entire batch
}
```

### Test File

```
tests/test_batch_tree_queries.cpp
```

### Test Fixture

```cpp
class BatchTreeQueryTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create test particle distribution
        particles = create_random_particles(10000);

        // Build tree
        tree = std::make_unique<BHTree>();
        tree->initialize(create_test_params());
        tree->resize(particles.size());
        tree->make(particles, particles.size());
    }

    std::vector<SPHParticle> particles;
    std::unique_ptr<BHTree> tree;

    // Get particles in same leaf cell
    std::vector<std::vector<int>> group_particles_by_leaf_cell();
};
```

### Test Cases

#### 4.1 Correctness Tests

```cpp
TEST_F(BatchTreeQueryTest, BatchQueryMatchesIndividualQueries) {
    // GIVEN: Particles grouped by leaf cell
    // WHEN: Compute gravity via batch query
    // THEN: Results match individual queries exactly

    auto groups = group_particles_by_leaf_cell();

    for (const auto& group : groups) {
        if (group.size() < 2) continue;

        // Compute individual results
        std::vector<vec_t> individual_results(group.size());
        for (size_t i = 0; i < group.size(); ++i) {
            tree->tree_force(particles[group[i]]);
            individual_results[i] = particles[group[i]].grav_acc;
        }

        // Reset and compute batch result
        for (int idx : group) {
            particles[idx].grav_acc = vec_t(0.0);
        }
        tree->tree_force_batch(group, particles);

        // Compare
        for (size_t i = 0; i < group.size(); ++i) {
            for (int d = 0; d < DIM; ++d) {
                EXPECT_NEAR(particles[group[i]].grav_acc[d],
                           individual_results[i][d], kAbsTol);
            }
        }
    }
}

TEST_F(BatchTreeQueryTest, BatchQueryWorksForSingleParticle) {
    // GIVEN: Batch with single particle
    // WHEN: Compute batch query
    // THEN: Result matches individual query

    std::vector<int> single = {0};

    tree->tree_force(particles[0]);
    vec_t individual = particles[0].grav_acc;

    particles[0].grav_acc = vec_t(0.0);
    tree->tree_force_batch(single, particles);

    for (int d = 0; d < DIM; ++d) {
        EXPECT_EQ(particles[0].grav_acc[d], individual[d]);
    }
}

TEST_F(BatchTreeQueryTest, BatchQueryHandlesEmptyBatch) {
    // GIVEN: Empty batch
    // WHEN: Call batch query
    // THEN: No crash, no side effects

    std::vector<int> empty;
    EXPECT_NO_THROW(tree->tree_force_batch(empty, particles));
}
```

#### 4.2 Node Visit Counting Tests

```cpp
TEST_F(BatchTreeQueryTest, BatchQueryVisitsFewerNodes) {
    // GIVEN: Particles in same leaf cell
    // WHEN: Compare node visits for batch vs individual
    // THEN: Batch visits fewer nodes

    auto groups = group_particles_by_leaf_cell();

    // Find a group with multiple particles
    std::vector<int> test_group;
    for (const auto& g : groups) {
        if (g.size() >= 4) {
            test_group = g;
            break;
        }
    }

    if (test_group.empty()) {
        GTEST_SKIP() << "No suitable particle groups found";
    }

    // Count individual traversal visits
    int individual_visits = 0;
    for (int idx : test_group) {
        individual_visits += tree->tree_force_counting(particles[idx]);
    }

    // Count batch traversal visits
    int batch_visits = tree->tree_force_batch_counting(test_group, particles);

    std::cout << "Individual visits: " << individual_visits
              << ", Batch visits: " << batch_visits << std::endl;

    EXPECT_LT(batch_visits, individual_visits);
}
```

#### 4.3 Performance Tests

```cpp
TEST_F(BatchTreeQueryTest, BatchQueryIsFaster) {
    // GIVEN: All particles grouped by leaf cell
    // WHEN: Compare batch vs individual computation time
    // THEN: Batch is faster

    auto groups = group_particles_by_leaf_cell();

    // Benchmark individual
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < particles.size(); ++i) {
        tree->tree_force(particles[i]);
    }
    auto individual_time = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - start).count();

    // Reset
    for (auto& p : particles) {
        p.grav_acc = vec_t(0.0);
    }

    // Benchmark batch
    start = std::chrono::high_resolution_clock::now();
    for (const auto& group : groups) {
        tree->tree_force_batch(group, particles);
    }
    auto batch_time = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - start).count();

    std::cout << "Individual: " << individual_time << " ms, Batch: " << batch_time << " ms" << std::endl;

    EXPECT_LT(batch_time, individual_time);
}
```

### Expected API

```cpp
class BHTree {
public:
    // Existing individual query
    void tree_force(SPHParticle& p_i);

    // New batch query (particles in same leaf cell)
    void tree_force_batch(const std::vector<int>& particle_indices,
                          std::vector<SPHParticle>& particles);

    // Debug versions with node visit counting
    int tree_force_counting(SPHParticle& p_i);
    int tree_force_batch_counting(const std::vector<int>& indices,
                                   std::vector<SPHParticle>& particles);

    // Group particles by leaf cell for batching
    std::vector<std::vector<int>> group_particles_by_leaf() const;
};
```

---

## Feature 5: Quadrupole Moments

**Priority:** LOW - Accuracy improvement allowing larger opening angle

### Background

Current tree uses monopole approximation only:
- Mass `M` and center of mass `r_cm`
- Error scales as `(l/r)^2` where `l` is cell size

Quadrupole tensor captures mass distribution shape:
- Allows larger opening angle `theta` for same accuracy
- Error scales as `(l/r)^4` with quadrupole correction

### Quadrupole Tensor Definition

For a mass distribution:
```
Q_ij = sum_k m_k * (3 * r_ki * r_kj - |r_k|^2 * delta_ij)
```

Where `r_k = pos_k - center_of_mass`.

The gravitational potential including quadrupole:
```
phi = -GM/r - (G/2r^5) * sum_{ij} Q_ij * r_i * r_j
```

### Test File

```
tests/test_quadrupole_moments.cpp
```

### Test Fixture

```cpp
class QuadrupoleMomentTest : public ::testing::Test {
protected:
    // Compute quadrupole tensor for particle set
    QuadrupoleTensor compute_quadrupole(const std::vector<SPHParticle>& particles,
                                        const vec_t& center_of_mass);

    // Analytic quadrupole for test cases
    QuadrupoleTensor analytic_rod_quadrupole(double M, double L);
    QuadrupoleTensor analytic_disk_quadrupole(double M, double R);
};
```

### Test Cases

#### 5.1 Quadrupole Tensor Computation

```cpp
TEST_F(QuadrupoleMomentTest, GivenPointMass_WhenComputingQuadrupole_ThenIsZero) {
    // GIVEN: Single point mass at origin
    // WHEN: Compute quadrupole tensor
    // THEN: All components are zero

    std::vector<SPHParticle> particles(1);
    particles[0].pos = vec_t(0.0);
    particles[0].mass = 1.0;

    auto Q = compute_quadrupole(particles, vec_t(0.0));

    EXPECT_NEAR(Q.xx, 0.0, kAbsTol);
    EXPECT_NEAR(Q.xy, 0.0, kAbsTol);
    EXPECT_NEAR(Q.xz, 0.0, kAbsTol);
    EXPECT_NEAR(Q.yy, 0.0, kAbsTol);
    EXPECT_NEAR(Q.yz, 0.0, kAbsTol);
    EXPECT_NEAR(Q.zz, 0.0, kAbsTol);
}

TEST_F(QuadrupoleMomentTest, GivenSymmetricSphere_WhenComputingQuadrupole_ThenIsZero) {
    // GIVEN: Symmetric spherical mass distribution
    // WHEN: Compute quadrupole tensor
    // THEN: All components are zero (by symmetry)

    auto particles = create_uniform_sphere(1000, 1.0, 1.0);
    vec_t com = compute_center_of_mass(particles);

    auto Q = compute_quadrupole(particles, com);

    // Should be small (not exactly zero due to discrete sampling)
    EXPECT_NEAR(Q.xx, 0.0, 0.1);
    EXPECT_NEAR(Q.xy, 0.0, 0.1);
    EXPECT_NEAR(Q.yy, 0.0, 0.1);
}

TEST_F(QuadrupoleMomentTest, GivenRodAlongZ_WhenComputingQuadrupole_ThenHasCorrectStructure) {
    // GIVEN: Mass distributed along z-axis (rod)
    // WHEN: Compute quadrupole tensor
    // THEN: Q_zz = -2*Q_xx = -2*Q_yy, off-diagonal = 0

    // Place particles along z-axis
    std::vector<SPHParticle> particles(100);
    double L = 2.0;  // Total length
    double m = 1.0 / 100;
    for (int i = 0; i < 100; ++i) {
        double z = L * (i - 49.5) / 100.0;
        particles[i].pos = vec_t(0.0, 0.0, z);
        particles[i].mass = m;
    }

    vec_t com = {0.0, 0.0, 0.0};
    auto Q = compute_quadrupole(particles, com);

    // For rod along z: Q_zz = M*L^2/12 * 2 = M*L^2/6
    // Q_xx = Q_yy = -Q_zz/2
    double M = 1.0;
    double expected_Qzz = M * L * L / 6.0;  // Approximate

    // Check traceless condition: Q_xx + Q_yy + Q_zz = 0
    EXPECT_NEAR(Q.xx + Q.yy + Q.zz, 0.0, kAbsTol);

    // Check symmetry: Q_xx = Q_yy
    EXPECT_NEAR(Q.xx, Q.yy, kAbsTol);

    // Check off-diagonal are zero
    EXPECT_NEAR(Q.xy, 0.0, kAbsTol);
    EXPECT_NEAR(Q.xz, 0.0, kAbsTol);
    EXPECT_NEAR(Q.yz, 0.0, kAbsTol);
}
```

#### 5.2 Tree Node Quadrupole Tests

```cpp
TEST_F(QuadrupoleMomentTest, TreeNodeQuadrupoleIsAccumulated) {
    // GIVEN: Tree built from particle distribution
    // WHEN: Traverse tree
    // THEN: Each node has accumulated quadrupole from children

    auto particles = create_random_particles(1000);
    BHTreeWithQuadrupole tree;
    tree.initialize(create_test_params());
    tree.make(particles);

    // Check root node quadrupole matches direct computation
    vec_t com = tree.get_root_center_of_mass();
    auto Q_direct = compute_quadrupole(particles, com);
    auto Q_root = tree.get_root_quadrupole();

    EXPECT_NEAR(Q_root.xx, Q_direct.xx, 1e-4);
    EXPECT_NEAR(Q_root.xy, Q_direct.xy, 1e-4);
    EXPECT_NEAR(Q_root.xz, Q_direct.xz, 1e-4);
    EXPECT_NEAR(Q_root.yy, Q_direct.yy, 1e-4);
    EXPECT_NEAR(Q_root.yz, Q_direct.yz, 1e-4);
    EXPECT_NEAR(Q_root.zz, Q_direct.zz, 1e-4);
}
```

#### 5.3 Gravity Accuracy with Quadrupole

```cpp
TEST_F(QuadrupoleMomentTest, QuadrupoleImprovesAccuracyForAsymmetricDistribution) {
    // GIVEN: Asymmetric mass distribution (e.g., elongated)
    // WHEN: Compute gravity with monopole-only vs monopole+quadrupole
    // THEN: Quadrupole correction is more accurate

    // Create elongated distribution along z-axis
    auto particles = create_elongated_distribution(1000, 3.0);  // 3:1 aspect ratio

    // Test position far from distribution
    vec_t test_pos = {5.0, 0.0, 0.0};

    // Direct summation (ground truth)
    vec_t g_direct = compute_direct_gravity(test_pos, particles);

    // Monopole approximation
    double M = total_mass(particles);
    vec_t com = compute_center_of_mass(particles);
    vec_t r = test_pos - com;
    double r_mag = std::abs(r);
    vec_t g_monopole = -r * (G * M / (r_mag * r_mag * r_mag));

    // Monopole + quadrupole
    auto Q = compute_quadrupole(particles, com);
    vec_t g_quad = g_monopole + compute_quadrupole_correction(test_pos, com, Q, M);

    double error_monopole = std::abs(g_monopole - g_direct) / std::abs(g_direct);
    double error_quad = std::abs(g_quad - g_direct) / std::abs(g_direct);

    std::cout << "Monopole error: " << error_monopole
              << ", Quadrupole error: " << error_quad << std::endl;

    EXPECT_LT(error_quad, error_monopole);
}

TEST_F(QuadrupoleMomentTest, QuadrupoleAllowsLargerOpeningAngle) {
    // GIVEN: Tree with quadrupole moments
    // WHEN: Use larger opening angle theta
    // THEN: Same accuracy as monopole-only with smaller theta

    auto particles = create_random_particles(10000);

    // Baseline: monopole with theta=0.5
    BHTree tree_mono;
    tree_mono.initialize(create_params_with_theta(0.5));
    tree_mono.make(particles);

    vec_t g_baseline = compute_tree_gravity(particles[0], tree_mono);

    // Quadrupole with theta=0.7 (40% larger)
    BHTreeWithQuadrupole tree_quad;
    tree_quad.initialize(create_params_with_theta(0.7));
    tree_quad.make(particles);

    vec_t g_direct = compute_direct_gravity(particles[0].pos, particles);
    vec_t g_quad = compute_tree_gravity_with_quad(particles[0], tree_quad);

    double error_mono = std::abs(g_baseline - g_direct) / std::abs(g_direct);
    double error_quad = std::abs(g_quad - g_direct) / std::abs(g_direct);

    // Quadrupole with larger theta should be at least as accurate
    EXPECT_LE(error_quad, error_mono * 1.1);  // Allow 10% margin
}
```

#### 5.4 Performance Tests

```cpp
TEST_F(QuadrupoleMomentTest, QuadrupoleBuildTimeOverhead) {
    // GIVEN: Particle distribution
    // WHEN: Build tree with and without quadrupole
    // THEN: Quadrupole adds acceptable overhead (<50%)

    auto particles = create_random_particles(100000);

    // Time monopole-only build
    BHTree tree_mono;
    auto start = std::chrono::high_resolution_clock::now();
    tree_mono.initialize(create_test_params());
    tree_mono.make(particles);
    auto mono_build = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - start).count();

    // Time quadrupole build
    BHTreeWithQuadrupole tree_quad;
    start = std::chrono::high_resolution_clock::now();
    tree_quad.initialize(create_test_params());
    tree_quad.make(particles);
    auto quad_build = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - start).count();

    double overhead = static_cast<double>(quad_build - mono_build) / mono_build;
    std::cout << "Build overhead: " << overhead * 100 << "%" << std::endl;

    EXPECT_LT(overhead, 0.5);  // Less than 50% overhead
}
```

### Expected API

```cpp
// Symmetric quadrupole tensor (6 independent components)
struct QuadrupoleTensor {
    real xx, xy, xz;
    real     yy, yz;
    real         zz;

    QuadrupoleTensor() : xx(0), xy(0), xz(0), yy(0), yz(0), zz(0) {}

    // Add contribution from particle at position r relative to COM
    void add_particle(const vec_t& r, real mass);

    // Combine with another tensor (for tree accumulation)
    QuadrupoleTensor& operator+=(const QuadrupoleTensor& other);

    // Shift tensor when moving reference point (parallel axis theorem)
    void shift_origin(const vec_t& displacement, real total_mass);
};

// Extended tree node with quadrupole
class BHNodeWithQuadrupole : public BHNode {
public:
    QuadrupoleTensor quadrupole;

    // Build quadrupole during tree construction
    void compute_quadrupole();
};

// Tree with quadrupole support
class BHTreeWithQuadrupole : public BHTree {
public:
    void tree_force_with_quadrupole(SPHParticle& p_i);
};

// Compute quadrupole gravity correction
vec_t quadrupole_gravity_correction(const vec_t& r, const QuadrupoleTensor& Q);
```

---

## Summary: Test File Structure

| Feature | Test File | Priority |
|---------|-----------|----------|
| Cubic Spline Lookup | `tests/test_hernquist_katz_lookup_table.cpp` | HIGH |
| Analytic Gravity | `tests/test_analytic_gravity.cpp` | HIGH |
| SoA Gravity | `tests/test_soa_gravity.cpp` | MEDIUM |
| Batch Tree Queries | `tests/test_batch_tree_queries.cpp` | MEDIUM |
| Quadrupole Moments | `tests/test_quadrupole_moments.cpp` | LOW |

## Files to Create/Modify

### New Files
1. `tests/test_hernquist_katz_lookup_table.cpp`
2. `tests/test_analytic_gravity.cpp`
3. `tests/test_soa_gravity.cpp`
4. `tests/test_batch_tree_queries.cpp`
5. `tests/test_quadrupole_moments.cpp`

### Modified Files
1. `include/softening_lookup_table.hpp` - Add HernquistKatzLookupTable
2. `src/softening_lookup_table.cpp` - Implement HernquistKatzLookupTable
3. `src/bhtree.cpp` - Add batch queries, quadrupole support
4. `include/bhtree.hpp` - Add QuadrupoleTensor, batch API

## Acceptance Criteria Summary

| Feature | Accuracy Criterion | Performance Criterion |
|---------|-------------------|----------------------|
| HK Lookup | 1e-6 vs polynomial | 1.5x speedup |
| Analytic Test | Converge as softening -> 0 | N/A |
| SoA | Exact match to AoS | >1.0x speedup |
| Batch | Exact match to individual | Fewer node visits |
| Quadrupole | Allows 40% larger theta | <50% build overhead |

---

## Implementation Order

1. **Feature 1 (HK Lookup)**: Direct port of existing WendlandC4LookupTable pattern
2. **Feature 2 (Analytic Test)**: Provides validation for all gravity implementations
3. **Feature 3 (SoA)**: Can be implemented independently, test against Feature 2
4. **Feature 4 (Batch)**: Requires tree understanding, test against Feature 2
5. **Feature 5 (Quadrupole)**: Most complex, depends on understanding Features 1-4
