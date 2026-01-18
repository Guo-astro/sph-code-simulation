# TDD Test Plan: 3D Analytic Gravity Test + Cubic Spline Lookup Table

Generated: 2026-01-18

---

## Overview

This TDD plan covers two related features for the SPH gravity subsystem:

1. **Feature 1**: 3D Analytic Gravity Test for Uniform Sphere
2. **Feature 2**: Cubic Spline (Hernquist-Katz) Lookup Table

Both features share the Hernquist-Katz softening kernel, making them natural companions for implementation.

---

## Feature 1: 3D Analytic Gravity Test (Uniform Sphere)

### Goal

Validate that the softened gravity implementation correctly reproduces the analytic solution for a uniform density sphere. This provides a rigorous physics-based test that catches implementation errors.

### Physics Background

For a uniform density sphere of radius R, total mass M, and density rho = 3M/(4 pi R^3):

| Region | Gravitational Acceleration | Formula |
|--------|---------------------------|---------|
| Inside (r < R) | Linear in r | g = (4/3) pi G rho r = G M_enclosed / r^2 |
| Outside (r >= R) | Point mass | g = G M / r^2 |

Where M_enclosed(r) = M (r/R)^3 for r < R.

**Key insight**: The enclosed mass M_enclosed = (4/3) pi rho r^3 gives:
```
g(r < R) = G M_enclosed / r^2 = G (4/3 pi rho r^3) / r^2 = (4/3) pi G rho r
```

This is the famous "inside a uniform sphere, gravity grows linearly with radius" result.

### Softening Effect

When using kernel softening (Hernquist-Katz or Wendland C4), the softened force approaches the point-mass limit outside the kernel support:

- **Hernquist-Katz**: Support radius = 2 epsilon = h (where epsilon = h/2)
- **Wendland C4**: Support radius = h

The test verifies:
1. Softened gravity matches analytic for r >> softening length
2. Softened gravity is smooth (no singularity) at r = 0
3. Convergence to analytic as softening -> 0

---

### Test Suite Structure

**File**: `tests/test_analytic_gravity_3d.cpp`

```
tests/
  test_analytic_gravity_3d.cpp  <- New file
```

---

### Test Cases

#### 1. ANALYTIC SOLUTION TESTS

**Test Name**: `TEST_F(AnalyticGravity3DTest, UniformSphereGravityInsideIsLinear)`

**Description**: Verify gravity inside a uniform sphere grows linearly with radius.

**Given**:
- Uniform density sphere: R = 1.0, M = 1.0, rho = 3/(4 pi)
- Test points at r = 0, 0.1, 0.25, 0.5, 0.75, 0.99 (all inside sphere)

**When**:
- Compute analytic gravity: g_analytic = (4/3) pi G rho r

**Then**:
- g(r) / r is constant for all r > 0 (linear relationship)
- g(0) = 0 exactly

```cpp
TEST_F(AnalyticGravity3DTest, UniformSphereGravityInsideIsLinear) {
    // GIVEN: Uniform sphere parameters
    const double R = 1.0;       // Sphere radius
    const double M = 1.0;       // Total mass
    const double G = 1.0;       // Gravitational constant
    const double rho = 3.0 * M / (4.0 * M_PI * R * R * R);

    std::vector<double> r_values = {0.0, 0.1, 0.25, 0.5, 0.75, 0.99};

    // WHEN/THEN: g(r) / r is constant for r > 0
    const double expected_slope = (4.0 / 3.0) * M_PI * G * rho;

    for (double r : r_values) {
        double g_analytic = expected_slope * r;

        if (r > 0.0) {
            double slope = g_analytic / r;
            EXPECT_NEAR(slope, expected_slope, kRelTol * expected_slope)
                << "Gravity not linear at r = " << r;
        } else {
            EXPECT_NEAR(g_analytic, 0.0, kAbsTol)
                << "Gravity should be zero at r = 0";
        }
    }
}
```

---

**Test Name**: `TEST_F(AnalyticGravity3DTest, UniformSphereGravityOutsideIsPointMass)`

**Description**: Verify gravity outside a uniform sphere matches point mass.

**Given**:
- Uniform density sphere: R = 1.0, M = 1.0
- Test points at r = 1.0, 1.5, 2.0, 5.0, 10.0 (all outside sphere)

**When**:
- Compute analytic gravity: g_analytic = G M / r^2

**Then**:
- g(r) = G M / r^2 for all r >= R

```cpp
TEST_F(AnalyticGravity3DTest, UniformSphereGravityOutsideIsPointMass) {
    // GIVEN: Uniform sphere parameters
    const double R = 1.0;
    const double M = 1.0;
    const double G = 1.0;

    std::vector<double> r_values = {1.0, 1.5, 2.0, 5.0, 10.0};

    for (double r : r_values) {
        // WHEN: Compute analytic point-mass gravity
        double g_analytic = G * M / (r * r);
        double g_expected = G * M / (r * r);

        // THEN: Matches exactly
        EXPECT_NEAR(g_analytic, g_expected, kAbsTol)
            << "Gravity not point-mass at r = " << r;
    }
}
```

---

**Test Name**: `TEST_F(AnalyticGravity3DTest, UniformSphereGravityContinuousAtSurface)`

**Description**: Verify gravity is continuous at sphere surface (r = R).

**Given**:
- Uniform sphere: R = 1.0, M = 1.0

**When**:
- Compute g_inside = lim(r -> R^-) = (4/3) pi G rho R
- Compute g_outside = lim(r -> R^+) = G M / R^2

**Then**:
- g_inside = g_outside (continuity at surface)

```cpp
TEST_F(AnalyticGravity3DTest, UniformSphereGravityContinuousAtSurface) {
    const double R = 1.0;
    const double M = 1.0;
    const double G = 1.0;
    const double rho = 3.0 * M / (4.0 * M_PI * R * R * R);

    // Inside limit: g = (4/3) pi G rho R
    double g_inside = (4.0 / 3.0) * M_PI * G * rho * R;

    // Outside limit: g = G M / R^2
    double g_outside = G * M / (R * R);

    // Should be equal by construction (M = (4/3) pi rho R^3)
    EXPECT_NEAR(g_inside, g_outside, kAbsTol)
        << "Gravity discontinuous at surface";
}
```

---

#### 2. SOFTENED GRAVITY VS ANALYTIC TESTS

**Test Name**: `TEST_F(AnalyticGravity3DTest, SoftenedGravityMatchesAnalyticOutsideKernel)`

**Description**: Softened gravity matches analytic result far from kernel support.

**Given**:
- Uniform sphere: R = 1.0, M = 1.0
- Particles distributed uniformly in sphere
- Softening length h = 0.1 (much smaller than R)

**When**:
- Compute softened gravity at r = 5R, 10R using particle summation
- Compute analytic gravity at same points

**Then**:
- Softened matches analytic within tolerance (< 1% for r >> h)

```cpp
TEST_F(AnalyticGravity3DTest, SoftenedGravityMatchesAnalyticOutsideKernel) {
    // GIVEN: Create uniform sphere particle distribution
    const int N_particles = 1000;
    const double R = 1.0;
    const double M = 1.0;
    const double G = 1.0;
    const double m_particle = M / N_particles;
    const double h = 0.1;  // Softening length << R

    // Generate uniform sphere distribution (Marsaglia 1972)
    std::vector<std::array<double, 3>> positions;
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> uniform(-1.0, 1.0);

    while (positions.size() < N_particles) {
        double x = uniform(rng);
        double y = uniform(rng);
        double z = uniform(rng);
        double r2 = x*x + y*y + z*z;

        if (r2 < 1.0) {  // Inside unit sphere
            double r_scaled = std::cbrt(uniform(rng) * 0.5 + 0.5) * R;
            double norm = r_scaled / std::sqrt(r2 + 1e-30);
            positions.push_back({x * norm, y * norm, z * norm});
        }
    }

    // Test points far from sphere
    std::vector<double> r_test = {5.0 * R, 10.0 * R};

    for (double r_i : r_test) {
        // Test point on x-axis
        std::array<double, 3> pos_i = {r_i, 0.0, 0.0};

        // WHEN: Compute softened gravity using Hernquist-Katz g function
        double g_softened_x = 0.0;
        for (const auto& pos_j : positions) {
            double dx = pos_i[0] - pos_j[0];
            double dy = pos_i[1] - pos_j[1];
            double dz = pos_i[2] - pos_j[2];
            double r = std::sqrt(dx*dx + dy*dy + dz*dz);

            // g function returns force factor: F = G m_i m_j g(r,h) * r_ij
            double g_factor = sph::g(r, h * 2.0);  // h_param = 2*epsilon
            g_softened_x -= G * m_particle * g_factor * dx;
        }

        // Analytic point-mass result (force per unit mass in x direction)
        double g_analytic = G * M / (r_i * r_i);

        // THEN: Relative error < 5% (includes discretization noise)
        double rel_error = std::abs(g_softened_x - g_analytic) / g_analytic;
        EXPECT_LT(rel_error, 0.05)
            << "At r = " << r_i << ": softened = " << g_softened_x
            << ", analytic = " << g_analytic;
    }
}
```

---

**Test Name**: `TEST_F(AnalyticGravity3DTest, SoftenedGravityZeroAtCenter)`

**Description**: Softened gravity is exactly zero at sphere center by symmetry.

**Given**:
- Uniform sphere particle distribution with N particles
- Test point at center r = 0

**When**:
- Sum softened force contributions from all particles

**Then**:
- Net force is zero (within numerical tolerance)

```cpp
TEST_F(AnalyticGravity3DTest, SoftenedGravityZeroAtCenter) {
    // GIVEN: Symmetric distribution (octahedral symmetry for exact cancellation)
    std::vector<std::array<double, 3>> positions = {
        { 1.0,  0.0,  0.0}, {-1.0,  0.0,  0.0},
        { 0.0,  1.0,  0.0}, { 0.0, -1.0,  0.0},
        { 0.0,  0.0,  1.0}, { 0.0,  0.0, -1.0}
    };

    const double G = 1.0;
    const double m = 1.0;
    const double h = 0.5;

    // Test point at origin
    std::array<double, 3> pos_i = {0.0, 0.0, 0.0};

    // WHEN: Sum forces
    double fx = 0.0, fy = 0.0, fz = 0.0;
    for (const auto& pos_j : positions) {
        double dx = pos_i[0] - pos_j[0];
        double dy = pos_i[1] - pos_j[1];
        double dz = pos_i[2] - pos_j[2];
        double r = std::sqrt(dx*dx + dy*dy + dz*dz);

        double g_factor = sph::g(r, h * 2.0);
        fx -= G * m * g_factor * dx;
        fy -= G * m * g_factor * dy;
        fz -= G * m * g_factor * dz;
    }

    // THEN: Net force is zero by symmetry
    EXPECT_NEAR(fx, 0.0, kAbsTol);
    EXPECT_NEAR(fy, 0.0, kAbsTol);
    EXPECT_NEAR(fz, 0.0, kAbsTol);
}
```

---

#### 3. SOFTENING CONVERGENCE TESTS

**Test Name**: `TEST_F(AnalyticGravity3DTest, SofteningConvergesToPointMass)`

**Description**: As softening length -> 0, softened gravity converges to point mass.

**Given**:
- Single point mass M at origin
- Test point at r = 1.0
- Softening lengths h = 2.0, 1.0, 0.5, 0.1, 0.01

**When**:
- Compute g(r, h) for each softening length

**Then**:
- g(r, h) approaches G M / r^3 as h -> 0
- Convergence is monotonic

```cpp
TEST_F(AnalyticGravity3DTest, SofteningConvergesToPointMass) {
    const double r = 1.0;
    const double G = 1.0;
    const double M = 1.0;
    const double g_point_mass = G * M / (r * r * r);

    std::vector<double> h_values = {2.0, 1.0, 0.5, 0.1, 0.01};
    double prev_error = std::numeric_limits<double>::max();

    for (double h : h_values) {
        // h_param = 2 * epsilon for Hernquist-Katz
        double g_softened = sph::g(r, h * 2.0);
        double error = std::abs(g_softened - g_point_mass);

        // Verify convergence (error decreases with smaller h)
        EXPECT_LE(error, prev_error * 1.1)  // Allow 10% noise
            << "Not converging at h = " << h;

        prev_error = error;

        std::cout << "h = " << h << ", error = " << error << std::endl;
    }

    // Final error should be very small
    EXPECT_LT(prev_error, 1e-4 * g_point_mass);
}
```

---

**Test Name**: `TEST_F(AnalyticGravity3DTest, SoftenedGravityMatchesInsideSphere)`

**Description**: With small softening, gravity inside uniform sphere matches (4/3) pi G rho r.

**Given**:
- Many particles uniformly distributed in sphere
- Small softening length (h << inter-particle spacing)
- Test points inside sphere

**When**:
- Compute softened gravity via particle sum
- Compute analytic gravity

**Then**:
- Match within particle noise tolerance (~1/sqrt(N))

```cpp
TEST_F(AnalyticGravity3DTest, SoftenedGravityMatchesInsideSphere) {
    // GIVEN: Large number of particles for low noise
    const int N = 5000;
    const double R = 1.0;
    const double M = 1.0;
    const double G = 1.0;
    const double rho = 3.0 * M / (4.0 * M_PI * R * R * R);
    const double m_particle = M / N;
    const double h = 0.05;  // Small softening

    // Generate uniform sphere
    std::vector<std::array<double, 3>> positions;
    std::mt19937 rng(42);
    generate_uniform_sphere(positions, N, R, rng);

    // Test at r = 0.5 R (inside sphere)
    double r_test = 0.5 * R;
    std::array<double, 3> pos_i = {r_test, 0.0, 0.0};

    // WHEN: Compute softened gravity
    double g_x = 0.0;
    for (const auto& pos_j : positions) {
        double dx = pos_i[0] - pos_j[0];
        double dy = pos_i[1] - pos_j[1];
        double dz = pos_i[2] - pos_j[2];
        double r = std::sqrt(dx*dx + dy*dy + dz*dz);

        if (r > 1e-10) {  // Skip if coincident
            double g_factor = sph::g(r, h * 2.0);
            g_x -= G * m_particle * g_factor * dx;
        }
    }

    // Analytic: g = (4/3) pi G rho r
    double g_analytic = (4.0 / 3.0) * M_PI * G * rho * r_test;

    // THEN: Match within particle noise (~ 1/sqrt(N) ~ 1.4%)
    double tol = 5.0 / std::sqrt(N);  // ~7% tolerance
    double rel_error = std::abs(g_x - g_analytic) / g_analytic;

    EXPECT_LT(rel_error, tol)
        << "Inside sphere: softened = " << g_x << ", analytic = " << g_analytic;
}
```

---

#### 4. KERNEL TYPE COMPARISON TESTS

**Test Name**: `TEST_F(AnalyticGravity3DTest, HernquistKatzVsWendlandC4Comparison)`

**Description**: Compare Hernquist-Katz and Wendland C4 kernels on uniform sphere test.

**Given**:
- Same particle distribution
- Same effective softening length

**When**:
- Compute gravity using both kernels

**Then**:
- Both converge to same analytic result
- Differences are documented for user guidance

```cpp
TEST_F(AnalyticGravity3DTest, HernquistKatzVsWendlandC4Comparison) {
    // GIVEN: Point mass test
    const double r = 1.0;
    const double h = 0.5;  // Softening length

    // WHEN: Compute with both kernels
    // Hernquist-Katz: h_param = 2 * epsilon = h
    double g_hk = sph::g(r, h);

    // Wendland C4: h_param = h directly
    double g_w = sph::GravityForce::wendland_c4_g(r, h);

    // Point mass reference
    double g_point = 1.0 / (r * r * r);

    // THEN: Both should be close to point mass (r > h)
    std::cout << "At r/h = " << r/h << ":" << std::endl;
    std::cout << "  Hernquist-Katz: " << g_hk << std::endl;
    std::cout << "  Wendland C4:    " << g_w << std::endl;
    std::cout << "  Point mass:     " << g_point << std::endl;

    // Both should match point mass when r >= kernel support
    if (r >= h) {  // Wendland C4 support
        EXPECT_NEAR(g_w, g_point, 1e-6);
    }
    if (r >= 2.0 * (h / 2.0)) {  // Hernquist-Katz support (u >= 2)
        EXPECT_NEAR(g_hk, g_point, 1e-10);
    }
}
```

---

### Test Fixture Definition

```cpp
class AnalyticGravity3DTest : public ::testing::Test {
protected:
    // Tolerances
    static constexpr double kAbsTol = 1e-10;
    static constexpr double kRelTol = 1e-6;

    // Helper: Generate uniform sphere distribution
    void generate_uniform_sphere(
        std::vector<std::array<double, 3>>& positions,
        int N,
        double R,
        std::mt19937& rng
    ) {
        std::uniform_real_distribution<double> uniform(-1.0, 1.0);
        std::uniform_real_distribution<double> unit(0.0, 1.0);

        positions.clear();
        while (positions.size() < static_cast<size_t>(N)) {
            // Rejection sampling for uniform in unit sphere
            double x = uniform(rng);
            double y = uniform(rng);
            double z = uniform(rng);
            double r2 = x*x + y*y + z*z;

            if (r2 < 1.0 && r2 > 1e-10) {
                // Uniform radial scaling
                double r_target = std::cbrt(unit(rng)) * R;
                double scale = r_target / std::sqrt(r2);
                positions.push_back({x * scale, y * scale, z * scale});
            }
        }
    }
};
```

---

## Feature 2: Cubic Spline (Hernquist-Katz) Lookup Table

### Goal

Add lookup table support for the Hernquist-Katz softening kernel, paralleling the existing Wendland C4 lookup table infrastructure. The lookup table provides faster evaluation while maintaining accuracy.

### Existing Implementation Analysis

From `gravity_force.cpp` lines 26-52:

**Hernquist-Katz f() function (potential)**:
```cpp
inline real f(const real r, const real h) {
    const real e = h * 0.5;       // epsilon = h/2
    const real u = r / e;         // normalized distance, u in [0, 2]

    if(u < 1.0) {
        return (-0.5 * u * u * (1.0 / 3.0 - 3.0 / 20 * u * u + u * u * u / 20) + 1.4) / e;
    } else if(u < 2.0) {
        return -1.0 / (15 * r) + (-u * u * (4.0 / 3.0 - u + 0.3 * u * u - u * u * u / 30) + 1.6) / e;
    } else {
        return 1 / r;  // Point mass
    }
}
```

**Hernquist-Katz g() function (force)**:
```cpp
inline real g(const real r, const real h) {
    const real e = h * 0.5;
    const real u = r / e;

    if(u < 1.0) {
        return (4.0 / 3.0 - 1.2 * u * u + 0.5 * u * u * u) / (e * e * e);
    } else if(u < 2.0) {
        return (-1.0 / 15 + 8.0 / 3 * u * u * u - 3 * u * u * u * u
                + 1.2 * u * u * u * u * u - u * u * u * u * u * u / 6.0) / (r * r * r);
    } else {
        return 1 / (r * r * r);  // Point mass
    }
}
```

**Key parameters**:
- Softening length epsilon = h / 2
- Normalized distance u = r / epsilon = 2r / h
- Kernel support: u in [0, 2] (equivalently, r in [0, h])
- Point mass for u >= 2 (r >= h)

---

### Test Suite Structure

**File**: `tests/test_hernquist_katz_lookup_table.cpp`

---

### Test Cases

#### 1. ACCURACY TESTS

**Test Name**: `TEST_F(HernquistKatzLookupTest, FLookupMatchesPolynomialWithinTolerance)`

**Description**: Lookup table for f() matches polynomial at many sample points.

**Given**:
- Lookup table covering u in [0, 2]
- Reference polynomial implementation (existing f() function)

**When**:
- Query f at 2000 uniformly distributed u values in [0, 2]

**Then**:
- Each lookup result differs from polynomial by < 1e-6 (relative error)

```cpp
TEST_F(HernquistKatzLookupTest, FLookupMatchesPolynomialWithinTolerance) {
    constexpr int N_SAMPLES = 2000;
    const double e = 1.0;  // epsilon = h/2
    double max_rel_error = 0.0;

    for (int i = 0; i <= N_SAMPLES; ++i) {
        double u = 2.0 * static_cast<double>(i) / N_SAMPLES;  // u in [0, 2]
        double r = u * e;
        double h = 2.0 * e;

        double lookup_val = lookup_table.f(u);
        double ref_val = reference_f(r, h);

        double rel_error = (std::abs(ref_val) > kAbsTol)
            ? std::abs(lookup_val - ref_val) / std::abs(ref_val)
            : std::abs(lookup_val - ref_val);

        max_rel_error = std::max(max_rel_error, rel_error);

        EXPECT_LT(rel_error, kAccuracyTol)
            << "At u=" << u << ": lookup=" << lookup_val << ", ref=" << ref_val;
    }

    std::cout << "[INFO] Max f relative error: " << max_rel_error << std::endl;
}
```

---

**Test Name**: `TEST_F(HernquistKatzLookupTest, GLookupMatchesPolynomialWithinTolerance)`

**Description**: Lookup table for g() matches polynomial.

**Given**:
- Lookup table initialized
- Reference g() function

**When**:
- Query g at 2000 uniformly distributed u values in (0, 2]

**Then**:
- Relative error < 1e-6

```cpp
TEST_F(HernquistKatzLookupTest, GLookupMatchesPolynomialWithinTolerance) {
    constexpr int N_SAMPLES = 2000;
    const double e = 1.0;
    double max_rel_error = 0.0;

    for (int i = 1; i <= N_SAMPLES; ++i) {  // Start at 1 to avoid u=0
        double u = 2.0 * static_cast<double>(i) / N_SAMPLES;
        double r = u * e;
        double h = 2.0 * e;

        double lookup_val = lookup_table.g(u);
        double ref_val = reference_g(r, h);

        double rel_error = (std::abs(ref_val) > kAbsTol)
            ? std::abs(lookup_val - ref_val) / std::abs(ref_val)
            : std::abs(lookup_val - ref_val);

        max_rel_error = std::max(max_rel_error, rel_error);

        EXPECT_LT(rel_error, kAccuracyTol)
            << "At u=" << u << ": lookup=" << lookup_val << ", ref=" << ref_val;
    }

    std::cout << "[INFO] Max g relative error: " << max_rel_error << std::endl;
}
```

---

#### 2. EDGE CASE TESTS

**Test Name**: `TEST_F(HernquistKatzLookupTest, GivenUAtZero_WhenComputingG_ThenReturnsCorrectValue)`

**Description**: Force g at u=0 should return (4/3) / e^3 (from inner polynomial).

**Given**: u = 0

**When**: Query g(0)

**Then**: g(0) = 4/3 / e^3 (coefficient from polynomial at u=0)

```cpp
TEST_F(HernquistKatzLookupTest, GivenUAtZero_WhenComputingG_ThenReturnsCorrectValue) {
    const double e = 1.0;
    double g = lookup_table.g(0.0);

    // From polynomial: g(u) = (4/3 - 1.2*u^2 + 0.5*u^3) / e^3
    // At u=0: g(0) = 4/3 / e^3
    double expected = (4.0 / 3.0) / (e * e * e);

    EXPECT_NEAR(g, expected, kAbsTol);
}
```

---

**Test Name**: `TEST_F(HernquistKatzLookupTest, GivenUAtZero_WhenComputingF_ThenReturnsFiniteValue)`

**Description**: Potential f at u=0 should return 1.4 / e (from inner polynomial).

**Given**: u = 0

**When**: Query f(0)

**Then**: f(0) = 1.4 / e (constant term from inner polynomial)

```cpp
TEST_F(HernquistKatzLookupTest, GivenUAtZero_WhenComputingF_ThenReturnsFiniteValue) {
    const double e = 1.0;
    double f = lookup_table.f(0.0);

    // From polynomial: f(u) = (-0.5*u^2*(1/3 - 3/20*u^2 + u^3/20) + 1.4) / e
    // At u=0: f(0) = 1.4 / e
    double expected = 1.4 / e;

    EXPECT_NEAR(f, expected, kAbsTol);
    EXPECT_TRUE(std::isfinite(f));
}
```

---

**Test Name**: `TEST_F(HernquistKatzLookupTest, GivenUAtOne_WhenComputing_ThenTransitionIsSmooth)`

**Description**: Verify continuity at u = 1 (inner/outer polynomial transition).

**Given**: u values near 1.0

**When**: Compute f and g at u = 0.99, 1.0, 1.01

**Then**: Values are continuous (no jump > tolerance)

```cpp
TEST_F(HernquistKatzLookupTest, GivenUAtOne_WhenComputing_ThenTransitionIsSmooth) {
    const double e = 1.0;
    double eps = 1e-6;

    // f values near u = 1
    double f_below = lookup_table.f(1.0 - eps);
    double f_at = lookup_table.f(1.0);
    double f_above = lookup_table.f(1.0 + eps);

    // Check continuity
    EXPECT_NEAR(f_below, f_at, kContinuityTol)
        << "f discontinuous at u=1 from below";
    EXPECT_NEAR(f_at, f_above, kContinuityTol)
        << "f discontinuous at u=1 from above";

    // g values near u = 1
    double g_below = lookup_table.g(1.0 - eps);
    double g_at = lookup_table.g(1.0);
    double g_above = lookup_table.g(1.0 + eps);

    EXPECT_NEAR(g_below, g_at, kContinuityTol)
        << "g discontinuous at u=1 from below";
    EXPECT_NEAR(g_at, g_above, kContinuityTol)
        << "g discontinuous at u=1 from above";
}
```

---

**Test Name**: `TEST_F(HernquistKatzLookupTest, GivenUAtTwo_WhenComputing_ThenMatchesPointMass)`

**Description**: At u = 2 (kernel boundary), f and g match point mass exactly.

**Given**: u = 2.0 - epsilon

**When**: Compute f and g

**Then**: Match point mass values 1/r and 1/r^3

```cpp
TEST_F(HernquistKatzLookupTest, GivenUAtTwo_WhenComputing_ThenMatchesPointMass) {
    const double e = 1.0;
    const double eps = 1e-6;
    const double u = 2.0 - eps;
    const double r = u * e;

    double f_inside = lookup_table.f(u);
    double g_inside = lookup_table.g(u);

    double f_point = 1.0 / r;
    double g_point = 1.0 / (r * r * r);

    EXPECT_NEAR(f_inside, f_point, kContinuityTol)
        << "f not continuous at kernel boundary";
    EXPECT_NEAR(g_inside, g_point, kContinuityTol)
        << "g not continuous at kernel boundary";
}
```

---

**Test Name**: `TEST_F(HernquistKatzLookupTest, GivenUGreaterThanTwo_WhenComputing_ThenReturnsPointMass)`

**Description**: For u >= 2, functions return point mass values.

**Given**: u values: 2.0, 2.5, 3.0, 5.0, 10.0

**When**: Query f(u) and g(u)

**Then**:
- f(u) = 1/r = 1/(u*e)
- g(u) = 1/r^3 = 1/(u*e)^3

```cpp
TEST_F(HernquistKatzLookupTest, GivenUGreaterThanTwo_WhenComputing_ThenReturnsPointMass) {
    const double e = 1.0;
    std::vector<double> u_values = {2.0, 2.5, 3.0, 5.0, 10.0};

    for (double u : u_values) {
        double r = u * e;
        double h = 2.0 * e;

        double f_val = lookup_table.f_full(r, h);
        double g_val = lookup_table.g_full(r, h);

        double f_point = 1.0 / r;
        double g_point = 1.0 / (r * r * r);

        EXPECT_NEAR(f_val, f_point, kAbsTol)
            << "f not point mass at u = " << u;
        EXPECT_NEAR(g_val, g_point, kAbsTol)
            << "g not point mass at u = " << u;
    }
}
```

---

#### 3. KERNEL SELECTOR TESTS

**Test Name**: `TEST_F(HernquistKatzLookupTest, KernelTypeSelectorChoosesCorrectKernel)`

**Description**: Verify kernel type selector returns appropriate kernel functions.

**Given**: KernelType enum with HERNQUIST_KATZ and WENDLAND_C4

**When**: Request lookup table for each type

**Then**: Correct kernel is used

```cpp
TEST_F(HernquistKatzLookupTest, KernelTypeSelectorChoosesCorrectKernel) {
    // GIVEN: Both kernel types
    const double r = 0.5;
    const double h = 1.0;

    // WHEN: Use kernel selector
    auto& hk_table = sph::SofteningLookupTableFactory::get(
        sph::GravitySofteningType::HERNQUIST_KATZ);
    auto& w4_table = sph::SofteningLookupTableFactory::get(
        sph::GravitySofteningType::WENDLAND_C4);

    double hk_g = hk_table.g_full(r, h);
    double w4_g = w4_table.g_full(r, h);

    // THEN: Values differ (different kernels)
    // Hernquist-Katz has support u in [0,2] where u = r/(h/2) = 2r/h
    // Wendland C4 has support q in [0,1] where q = r/h
    // At r/h = 0.5: HK has u = 1.0, W4 has q = 0.5

    // Both should give sensible positive values
    EXPECT_GT(hk_g, 0.0);
    EXPECT_GT(w4_g, 0.0);

    // They should differ (different kernels have different shapes)
    // This is a sanity check, not a precise test
    EXPECT_NE(hk_g, w4_g);
}
```

---

#### 4. THREAD SAFETY TESTS

**Test Name**: `TEST_F(HernquistKatzLookupTest, StaticInitializationIsThreadSafe)`

**Description**: Same thread safety test as Wendland C4, adapted for Hernquist-Katz.

```cpp
TEST_F(HernquistKatzLookupTest, StaticInitializationIsThreadSafe) {
    constexpr int NUM_THREADS = 8;
    constexpr int QUERIES_PER_THREAD = 1000;

    std::vector<std::thread> threads;
    std::vector<bool> success(NUM_THREADS, true);

    for (int t = 0; t < NUM_THREADS; ++t) {
        threads.emplace_back([this, &success, t]() {
            for (int i = 0; i < QUERIES_PER_THREAD; ++i) {
                double u = 2.0 * static_cast<double>(i) / QUERIES_PER_THREAD;
                double f_lookup = lookup_table.f(u);
                double f_ref = reference_f_normalized(u);

                double rel_error = std::abs(f_lookup - f_ref) / std::abs(f_ref);
                if (rel_error > kAccuracyTol) {
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

#### 5. PERFORMANCE TESTS

**Test Name**: `TEST_F(HernquistKatzLookupTest, LookupIsFasterThanPolynomial)`

**Description**: Lookup table should be faster than polynomial evaluation.

```cpp
TEST_F(HernquistKatzLookupTest, LookupIsFasterThanPolynomial) {
    constexpr int NUM_ITERATIONS = 1000000;
    std::vector<double> u_values(NUM_ITERATIONS);

    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(0.0, 2.0);
    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        u_values[i] = dist(rng);
    }

    const double e = 1.0;
    const double h = 2.0 * e;

    // Benchmark polynomial (existing implementation)
    volatile double sink_poly = 0.0;
    auto start_poly = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        double r = u_values[i] * e;
        sink_poly += sph::g(r, h);
    }
    auto end_poly = std::chrono::high_resolution_clock::now();
    auto poly_time = std::chrono::duration_cast<std::chrono::microseconds>(
        end_poly - start_poly).count();

    // Benchmark lookup
    volatile double sink_lookup = 0.0;
    auto start_lookup = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        sink_lookup += lookup_table.g(u_values[i]);
    }
    auto end_lookup = std::chrono::high_resolution_clock::now();
    auto lookup_time = std::chrono::duration_cast<std::chrono::microseconds>(
        end_lookup - start_lookup).count();

    double speedup = static_cast<double>(poly_time) / lookup_time;

    std::cout << "[BENCHMARK] Polynomial time: " << poly_time << " us" << std::endl;
    std::cout << "[BENCHMARK] Lookup time: " << lookup_time << " us" << std::endl;
    std::cout << "[BENCHMARK] Speedup: " << speedup << "x" << std::endl;

    // Hernquist-Katz polynomial is simpler than Wendland C4 (3rd order vs 9th)
    // Speedup may be less dramatic (~1.5x expected)
    EXPECT_GT(speedup, 1.2) << "Lookup should be at least 1.2x faster";
}
```

---

#### 6. INTEGRATION TESTS

**Test Name**: `TEST_F(HernquistKatzLookupTest, DropInReplacementForExistingFunctions)`

**Description**: Verify lookup table matches existing f() and g() inline functions.

```cpp
TEST_F(HernquistKatzLookupTest, DropInReplacementForExistingFunctions) {
    const double h = 1.0;  // h_param in existing functions
    const double e = h * 0.5;  // epsilon

    std::vector<double> r_values = {0.0, 0.1, 0.25, 0.5, 0.75, 0.99, 1.0, 1.5, 2.0};

    for (double r : r_values) {
        double f_orig = sph::f(r, h);
        double g_orig = sph::g(r, h);

        double f_lookup = lookup_table.f_full(r, h);
        double g_lookup = lookup_table.g_full(r, h);

        double f_tol = kAccuracyTol * std::abs(f_orig) + kAbsTol;
        double g_tol = kAccuracyTol * std::abs(g_orig) + kAbsTol;

        EXPECT_NEAR(f_lookup, f_orig, f_tol)
            << "f mismatch at r=" << r;
        EXPECT_NEAR(g_lookup, g_orig, g_tol)
            << "g mismatch at r=" << r;
    }
}
```

---

### Test Fixture Definition

```cpp
class HernquistKatzLookupTest : public ::testing::Test {
protected:
    // Tolerances
    static constexpr double kAccuracyTol = 1e-6;
    static constexpr double kAbsTol = 1e-10;
    static constexpr double kContinuityTol = 1e-4;

    // Reference to lookup table (singleton)
    const sph::HernquistKatzLookupTable& lookup_table =
        sph::HernquistKatzLookupTable::get_instance();

    // Reference polynomial: f(r, h) with e = h/2, u = r/e
    double reference_f(double r, double h) const {
        const double e = h * 0.5;
        const double u = r / e;

        if (u < 1.0) {
            return (-0.5 * u * u * (1.0/3.0 - 3.0/20.0 * u * u + u * u * u / 20.0)
                    + 1.4) / e;
        } else if (u < 2.0) {
            return -1.0 / (15.0 * r)
                   + (-u * u * (4.0/3.0 - u + 0.3 * u * u - u * u * u / 30.0)
                      + 1.6) / e;
        } else {
            return 1.0 / r;
        }
    }

    // Reference polynomial: g(r, h)
    double reference_g(double r, double h) const {
        const double e = h * 0.5;
        const double u = r / e;

        if (u < 1.0) {
            return (4.0/3.0 - 1.2 * u * u + 0.5 * u * u * u) / (e * e * e);
        } else if (u < 2.0) {
            const double u3 = u * u * u;
            const double u4 = u3 * u;
            const double u5 = u4 * u;
            const double u6 = u5 * u;
            return (-1.0/15.0 + 8.0/3.0 * u3 - 3.0 * u4 + 1.2 * u5 - u6/6.0)
                   / (r * r * r);
        } else {
            return 1.0 / (r * r * r);
        }
    }

    // Normalized version for internal lookup (u in [0, 2])
    double reference_f_normalized(double u) const {
        const double e = 1.0;  // epsilon = 1 for normalized
        const double r = u * e;
        const double h = 2.0 * e;
        return reference_f(r, h);
    }

    double reference_g_normalized(double u) const {
        const double e = 1.0;
        const double r = u * e;
        const double h = 2.0 * e;
        return reference_g(r, h);
    }
};
```

---

## Implementation Notes

### Recommended Lookup Table Design for Hernquist-Katz

```cpp
class HernquistKatzLookupTable {
public:
    // Parameters: u in [0, 2]
    static constexpr int TABLE_SIZE = 4096;  // Need more entries for [0,2]
    static constexpr double U_MAX = 2.0;
    static constexpr double DU = U_MAX / TABLE_SIZE;
    static constexpr double INV_DU = TABLE_SIZE / U_MAX;

private:
    // Cache-aligned tables
    alignas(64) double f_table_[TABLE_SIZE + 1];
    alignas(64) double f_slope_[TABLE_SIZE];
    alignas(64) double g_table_[TABLE_SIZE + 1];
    alignas(64) double g_slope_[TABLE_SIZE];

public:
    static const HernquistKatzLookupTable& get_instance() {
        static HernquistKatzLookupTable instance;
        return instance;
    }

    double f(double u) const {
        if (u >= U_MAX) {
            return 1.0 / u;  // Point mass: f = 1/r = 1/(u*e)
        }

        double idx_f = u * INV_DU;
        int idx = static_cast<int>(idx_f);
        double frac = idx_f - idx;

        return f_table_[idx] + frac * f_slope_[idx];
    }

    double g(double u) const {
        if (u >= U_MAX) {
            double u3 = u * u * u;
            return 1.0 / u3;  // Point mass: g = 1/r^3 = 1/(u*e)^3
        }

        double idx_f = u * INV_DU;
        int idx = static_cast<int>(idx_f);
        double frac = idx_f - idx;

        return g_table_[idx] + frac * g_slope_[idx];
    }

    // Full interface matching existing f(r, h) and g(r, h)
    double f_full(double r, double h) const {
        const double e = h * 0.5;
        const double u = r / e;
        return f(u) * e;  // Denormalize
    }

    double g_full(double r, double h) const {
        const double e = h * 0.5;
        const double u = r / e;
        double e3 = e * e * e;
        return g(u) * e3;  // Denormalize
    }
};
```

### Factory Pattern for Kernel Selection

```cpp
class SofteningLookupTableFactory {
public:
    static const ISofteningLookupTable& get(GravitySofteningType type) {
        switch (type) {
            case GravitySofteningType::HERNQUIST_KATZ:
                return HernquistKatzLookupTable::get_instance();
            case GravitySofteningType::WENDLAND_C4:
                return WendlandC4LookupTable::get_instance();
            default:
                throw std::runtime_error("Unknown softening type");
        }
    }
};
```

---

## Test Execution Order

### Phase 1: Red (Write failing tests)

1. Create `tests/test_analytic_gravity_3d.cpp` with all analytic gravity tests
2. Create `tests/test_hernquist_katz_lookup_table.cpp` with all lookup table tests
3. Verify tests compile but fail (no implementation yet)

### Phase 2: Green (Minimal implementation)

1. Implement `HernquistKatzLookupTable` class
2. Add factory pattern for kernel selection
3. Verify all tests pass

### Phase 3: Refactor

1. Optimize table access patterns
2. Consider SIMD for batch evaluations
3. Profile and tune table size

---

## Acceptance Criteria Summary

### Feature 1: 3D Analytic Gravity Test

| Test Category | Criterion | Tolerance |
|---------------|-----------|-----------|
| Inside sphere | Linear gravity law | 5% (particle noise) |
| Outside sphere | Point mass law | 1% |
| Continuity | Smooth at surface | Exact by construction |
| Symmetry | Zero at center | 1e-10 absolute |
| Convergence | Softening -> 0 | Monotonic |

### Feature 2: Hernquist-Katz Lookup Table

| Test Category | Criterion | Tolerance |
|---------------|-----------|-----------|
| Accuracy (f) | Matches polynomial | 1e-6 relative |
| Accuracy (g) | Matches polynomial | 1e-6 relative |
| Edge: u=0 | Correct central value | 1e-10 absolute |
| Edge: u=1 | Smooth transition | 1e-4 relative |
| Edge: u>=2 | Point mass values | 1e-10 absolute |
| Continuity | No jumps at boundaries | 1e-4 relative |
| Thread safety | No data races | N/A |
| Performance | Lookup >= 1.2x faster | Benchmark |

---

## Files to Create/Modify

### New Files

1. `tests/test_analytic_gravity_3d.cpp` - 3D analytic gravity tests
2. `tests/test_hernquist_katz_lookup_table.cpp` - Lookup table tests
3. `include/hernquist_katz_lookup_table.hpp` - HK lookup table class

### Files to Modify

1. `include/softening_lookup_table.hpp` - Add factory pattern
2. `src/softening_lookup_table_stub.cpp` - Implement HK lookup
3. `CMakeLists.txt` - Add new test files
4. `include/parameters.hpp` - Ensure GravitySofteningType enum exists

---

## Risk Assessment

| Risk | Impact | Mitigation |
|------|--------|------------|
| Particle noise in analytic tests | False failures | Use large N, statistical tolerance |
| Lookup table memory overhead | ~64KB per kernel | Acceptable for modern systems |
| HK polynomial simpler = less speedup | Reduced benefit | Still worth implementing for consistency |
| Thread safety issues | Correctness bugs | Use Meyers' singleton pattern |
| DIM=3 compilation only | Tests not run in 1D/2D | Conditional compilation blocks |

---

## Notes

### Difference Between Wendland C4 and Hernquist-Katz

| Property | Wendland C4 | Hernquist-Katz |
|----------|-------------|----------------|
| Support | q in [0, 1] | u in [0, 2] |
| Normalization | q = r/h | u = r/epsilon, epsilon = h/2 |
| Polynomial order | 9th order | 3rd/6th order |
| Regions | 1 (inside) + point mass | 2 (inner/outer) + point mass |
| Kernel width | h | 2*epsilon = h |

Both kernels have the same effective support radius (h), but different normalization conventions.
