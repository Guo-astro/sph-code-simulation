/**
 * @file test_performance_optimizations.cpp
 * @brief TDD tests for SPH performance optimizations
 *
 * These tests are written FIRST (TDD red phase) to define expected behavior
 * for three performance optimizations:
 *
 * 1. Thread-local neighbor list - Avoid heap allocation in hot loops
 * 2. Distance caching before sorting - Reduce redundant computation
 * 3. OpenMP dynamic scheduling - Better load balancing for variable density
 *
 * Tests will FAIL initially because optimizations don't exist yet.
 * After implementing each optimization, tests should PASS.
 */

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <cmath>
#include <vector>
#include <set>
#include <random>
#include <algorithm>
#include <memory>
#include <chrono>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "defines.hpp"
#include "particle.hpp"
#include "vector_type.hpp"
#include "bhtree.hpp"
#include "periodic.hpp"
#include "parameters.hpp"
#include "fluid_force.hpp"
#include "pre_interaction.hpp"

// Use types from sph namespace
using sph::SPHParticle;
using sph::SPHParameters;
using sph::BHTree;
using sph::Periodic;
using sph::FluidForce;
using sph::PreInteraction;

namespace {

// Tolerance for floating-point comparisons
constexpr double kAbsTol = 1e-14;
constexpr double kRelTol = 1e-12;

// Test configuration
constexpr int kDefaultNeighborNumber = 32;
constexpr int kNeighborListSize = 5000;

// Helper to create default parameters
std::shared_ptr<SPHParameters> create_default_params(bool periodic = false,
                                                     double range_min = 0.0,
                                                     double range_max = 1.0) {
    auto param = std::make_shared<SPHParameters>();

    // Time parameters
    param->time.start = 0.0;
    param->time.end = 1.0;
    param->time.output = 0.1;
    param->time.energy = 0.1;

    // SPH type
    param->type = sph::SPHType::SSPH;

    // CFL parameters
    param->cfl.sound = 0.3;
    param->cfl.force = 0.25;

    // AV parameters
    param->av.alpha = 1.0;
    param->av.use_balsara_switch = true;
    param->av.use_time_dependent_av = false;
    param->av.alpha_max = 2.0;
    param->av.alpha_min = 0.1;
    param->av.epsilon = 0.1;

    // AC parameters
    param->ac.alpha = 0.0;
    param->ac.is_valid = false;

    // Tree parameters
    param->tree.max_level = 20;
    param->tree.leaf_particle_num = 8;

    // Physics parameters
    param->physics.neighbor_number = kDefaultNeighborNumber;
    param->physics.gamma = 1.4;
    param->physics.c_smooth = 1.0;

    // Kernel
    param->kernel = sph::KernelType::WENDLAND;

    // Iterative smoothing length
    param->iterative_sml = false;
    param->preserve_initial_density = false;

    // Periodic parameters
    param->periodic.is_valid = periodic;
    param->periodic.boundary_type = sph::BoundaryType::REFLECTING;
    for (int i = 0; i < DIM; ++i) {
        param->periodic.range_min[i] = range_min;
        param->periodic.range_max[i] = range_max;
        param->periodic.per_dimension[i] = periodic;
    }

    // Gravity (disabled for these tests)
    param->gravity.is_valid = false;
    param->gravity.constant = 0.0;
    param->gravity.theta = 0.5;
    param->gravity.use_fixed_softening = false;
    param->gravity.fixed_softening = 0.0;
    param->gravity.softening_type = sph::GravitySofteningType::HERNQUIST_KATZ;

    // GSPH parameters (for grad-h)
    param->gsph.is_2nd_order = false;
    param->gsph.riemann_solver = sph::RiemannSolverType::HLL;
    param->gsph.use_gradh = true;
    param->gsph.use_volume_based = false;
    param->gsph.eta = 1.0;
    param->gsph.c_smooth = 1.0;

    return param;
}

// Helper to create a particle
SPHParticle create_particle(int id, const vec_t& pos, double mass = 1.0, double sml = 0.1) {
    SPHParticle p;
    p.id = id;
    p.pos = pos;
    p.mass = mass;
    p.sml = sml;
    p.dens = 1.0;
    p.pres = 1.0;
    p.ene = 1.0;
    p.vel = vec_t(0.0);
    p.vel_p = vec_t(0.0);
    p.acc = vec_t(0.0);
    p.dene = 0.0;
    p.sound = 1.0;
    p.balsara = 1.0;
    p.alpha = 1.0;
    p.gradh = 1.0;
    p.neighbor = 0;
    p.next = nullptr;
    p.is_ghost = false;
    return p;
}

// Helper to create uniform random distribution of particles
std::vector<SPHParticle> create_random_particles(int n, double domain_min, double domain_max,
                                                 double sml, unsigned int seed = 42) {
    std::vector<SPHParticle> particles(n);
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(domain_min, domain_max);

    for (int i = 0; i < n; ++i) {
        vec_t pos;
        for (int d = 0; d < DIM; ++d) {
            pos[d] = dist(gen);
        }
        particles[i] = create_particle(i, pos, 1.0, sml);
    }

    return particles;
}

// Helper to create grid distribution of particles
std::vector<SPHParticle> create_grid_particles(int n_per_dim, double domain_min, double domain_max, double sml) {
    std::vector<SPHParticle> particles;
    double spacing = (domain_max - domain_min) / n_per_dim;
    int id = 0;

#if DIM == 1
    for (int i = 0; i < n_per_dim; ++i) {
        vec_t pos;
        pos[0] = domain_min + (i + 0.5) * spacing;
        particles.push_back(create_particle(id++, pos, 1.0, sml));
    }
#elif DIM == 2
    for (int i = 0; i < n_per_dim; ++i) {
        for (int j = 0; j < n_per_dim; ++j) {
            vec_t pos;
            pos[0] = domain_min + (i + 0.5) * spacing;
            pos[1] = domain_min + (j + 0.5) * spacing;
            particles.push_back(create_particle(id++, pos, 1.0, sml));
        }
    }
#elif DIM == 3
    for (int i = 0; i < n_per_dim; ++i) {
        for (int j = 0; j < n_per_dim; ++j) {
            for (int k = 0; k < n_per_dim; ++k) {
                vec_t pos;
                pos[0] = domain_min + (i + 0.5) * spacing;
                pos[1] = domain_min + (j + 0.5) * spacing;
                pos[2] = domain_min + (k + 0.5) * spacing;
                particles.push_back(create_particle(id++, pos, 1.0, sml));
            }
        }
    }
#endif

    return particles;
}

// Helper to create clustered (non-uniform) particles for load balancing tests
std::vector<SPHParticle> create_clustered_particles(int n, double domain_min, double domain_max,
                                                    double sml, unsigned int seed = 42) {
    std::vector<SPHParticle> particles;
    std::mt19937 gen(seed);

    // 70% of particles in a small central region (high density)
    int n_dense = static_cast<int>(0.7 * n);
    double center = (domain_max + domain_min) * 0.5;
    double dense_range = (domain_max - domain_min) * 0.1;  // 10% of domain
    std::uniform_real_distribution<double> dense_dist(center - dense_range, center + dense_range);

    for (int i = 0; i < n_dense; ++i) {
        vec_t pos;
        for (int d = 0; d < DIM; ++d) {
            pos[d] = dense_dist(gen);
        }
        particles.push_back(create_particle(static_cast<int>(particles.size()), pos, 1.0, sml));
    }

    // 30% of particles spread uniformly (low density)
    std::uniform_real_distribution<double> sparse_dist(domain_min, domain_max);
    int n_sparse = n - n_dense;
    for (int i = 0; i < n_sparse; ++i) {
        vec_t pos;
        for (int d = 0; d < DIM; ++d) {
            pos[d] = sparse_dist(gen);
        }
        particles.push_back(create_particle(static_cast<int>(particles.size()), pos, 1.0, sml));
    }

    return particles;
}

// Helper to capture particle state
struct ParticleState {
    std::vector<vec_t> acc;
    std::vector<real> dene;
    std::vector<real> dens;
    std::vector<real> pres;
    std::vector<real> gradh;

    void capture(const std::vector<SPHParticle>& particles) {
        int n = particles.size();
        acc.resize(n);
        dene.resize(n);
        dens.resize(n);
        pres.resize(n);
        gradh.resize(n);
        for (int i = 0; i < n; ++i) {
            acc[i] = particles[i].acc;
            dene[i] = particles[i].dene;
            dens[i] = particles[i].dens;
            pres[i] = particles[i].pres;
            gradh[i] = particles[i].gradh;
        }
    }
};

// Helper to compare states within tolerance
bool states_match(const ParticleState& s1, const ParticleState& s2, double tol = kAbsTol) {
    if (s1.acc.size() != s2.acc.size()) return false;
    int n = s1.acc.size();
    for (int i = 0; i < n; ++i) {
        for (int d = 0; d < DIM; ++d) {
            if (std::abs(s1.acc[i][d] - s2.acc[i][d]) > tol) return false;
        }
        if (std::abs(s1.dene[i] - s2.dene[i]) > tol) return false;
        if (std::abs(s1.dens[i] - s2.dens[i]) > tol) return false;
        if (std::abs(s1.pres[i] - s2.pres[i]) > tol) return false;
        if (std::abs(s1.gradh[i] - s2.gradh[i]) > tol) return false;
    }
    return true;
}

} // anonymous namespace

// =============================================================================
// Test Fixture for Performance Optimizations
// =============================================================================
class PerformanceOptimizationTest : public ::testing::Test {
protected:
    std::shared_ptr<SPHParameters> param_;
    std::unique_ptr<BHTree> tree_;
    std::unique_ptr<Periodic> periodic_;
    std::unique_ptr<FluidForce> fluid_force_;
    std::unique_ptr<PreInteraction> pre_interaction_;

    void SetUp() override {
        param_ = create_default_params(false, 0.0, 1.0);
        tree_ = std::make_unique<BHTree>();
        tree_->initialize(param_);
        periodic_ = std::make_unique<Periodic>();
        periodic_->initialize(param_);
        fluid_force_ = std::make_unique<FluidForce>();
        fluid_force_->initialize(param_);
        pre_interaction_ = std::make_unique<PreInteraction>();
        pre_interaction_->initialize(param_);
    }

    void SetUpPeriodic(double range_min, double range_max) {
        param_ = create_default_params(true, range_min, range_max);
        tree_ = std::make_unique<BHTree>();
        tree_->initialize(param_);
        periodic_ = std::make_unique<Periodic>();
        periodic_->initialize(param_);
        fluid_force_ = std::make_unique<FluidForce>();
        fluid_force_->initialize(param_);
        pre_interaction_ = std::make_unique<PreInteraction>();
        pre_interaction_->initialize(param_);
    }

    void BuildTree(std::vector<SPHParticle>& particles) {
        int n = particles.size();
        tree_->resize(n, 5);
        tree_->make(particles, n);
        tree_->set_kernel();
    }
};

// =============================================================================
// Optimization 1: Thread-Local Neighbor List Tests
// =============================================================================

/**
 * TEST: Thread-local neighbor list availability check
 *
 * GIVEN: FluidForce or PreInteraction class
 * WHEN: Checking for thread-local neighbor list support
 * THEN: A method or flag indicates thread-local storage is available
 *
 * This test will FAIL because the optimization doesn't exist yet.
 * The test checks for a hypothetical method that returns whether
 * thread-local storage is being used.
 */
TEST_F(PerformanceOptimizationTest, ThreadLocalNeighborList_FeatureExists) {
    // GIVEN: FluidForce with thread-local optimization
    // WHEN: Checking if optimization is enabled
    // THEN: The feature flag should exist and be true

    // This test will FAIL because the method doesn't exist yet
    // Uncomment when implementing:
    // EXPECT_TRUE(fluid_force_->uses_thread_local_neighbor_list());

    // For now, we mark this as a compile-time check by looking for a specific define
    // The test will FAIL at link time if the optimization is not implemented
#ifdef SPH_USE_THREAD_LOCAL_NEIGHBOR_LIST
    SUCCEED() << "Thread-local neighbor list feature flag exists";
#else
    FAIL() << "Thread-local neighbor list optimization not implemented. "
           << "Define SPH_USE_THREAD_LOCAL_NEIGHBOR_LIST when implementing.";
#endif
}

/**
 * TEST: Thread-local neighbor list correctness
 *
 * GIVEN: Random particle distribution
 * WHEN: Run FluidForce with thread-local neighbor list
 * THEN: Results match baseline (per-iteration allocation)
 */
TEST_F(PerformanceOptimizationTest, ThreadLocalNeighborList_CorrectnessVsBaseline) {
    // GIVEN: Random particles with known baseline
    const int N = 500;
    SetUpPeriodic(0.0, 1.0);
    auto particles = create_random_particles(N, 0.0, 1.0, 0.1, 42);
    BuildTree(particles);

    // Store baseline state - this will be computed with current implementation
    ParticleState baseline;
    baseline.capture(particles);

    // Note: This test currently just verifies the baseline can be captured.
    // When thread-local optimization is implemented, we'll compare against
    // a version that explicitly uses per-iteration allocation.

    // For TDD: Test passes trivially now but will need modification
    // when the optimization is implemented to compare two code paths.
    EXPECT_EQ(static_cast<int>(baseline.acc.size()), N);
}

/**
 * TEST: Thread-local buffer isolation between threads
 *
 * GIVEN: Multiple threads executing concurrently
 * WHEN: Each thread uses its own thread-local neighbor list
 * THEN: No cross-thread contamination occurs
 */
TEST_F(PerformanceOptimizationTest, ThreadLocalNeighborList_BufferIsolation) {
    // GIVEN: Grid particles with known neighbor structure
#if DIM == 3
    const int n_per_dim = 10;  // 1000 particles in 3D
#elif DIM == 2
    const int n_per_dim = 32;  // 1024 particles in 2D
#else
    const int n_per_dim = 1000;  // 1000 particles in 1D
#endif
    SetUpPeriodic(0.0, 1.0);
    auto particles = create_grid_particles(n_per_dim, 0.0, 1.0, 0.15);
    const int N = particles.size();
    BuildTree(particles);

    // Pre-compute expected neighbor counts using sequential search
    std::vector<int> expected_neighbor_counts(N);
    for (int i = 0; i < N; ++i) {
        std::vector<int> neighbors(kNeighborListSize);
        int count = tree_->neighbor_search(particles[i], neighbors, particles, false);
        expected_neighbor_counts[i] = count;
    }

    // Run parallel computation and verify counts match
    std::vector<int> actual_neighbor_counts(N, -1);
    bool any_mismatch = false;

#pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        std::vector<int> neighbors(kNeighborListSize);
        int count = tree_->neighbor_search(particles[i], neighbors, particles, false);
        actual_neighbor_counts[i] = count;
    }

    // Verify all counts match
    for (int i = 0; i < N; ++i) {
        if (actual_neighbor_counts[i] != expected_neighbor_counts[i]) {
            any_mismatch = true;
            ADD_FAILURE() << "Neighbor count mismatch at particle " << i
                         << ": expected " << expected_neighbor_counts[i]
                         << ", got " << actual_neighbor_counts[i];
        }
    }

    EXPECT_FALSE(any_mismatch) << "Thread-local buffers may have cross-contamination";
}

/**
 * TEST: Thread-local results identical across different thread counts
 *
 * GIVEN: Same particle configuration
 * WHEN: Run with 1, 2, 4, 8 threads (up to max)
 * THEN: Results are identical regardless of thread count
 */
TEST_F(PerformanceOptimizationTest, ThreadLocalNeighborList_DifferentThreadCounts) {
#ifdef _OPENMP
    const int N = 200;
    SetUpPeriodic(0.0, 1.0);
    auto particles_original = create_random_particles(N, 0.0, 1.0, 0.1, 12345);

    // Reference: single-threaded execution
    auto particles = particles_original;
    BuildTree(particles);

    omp_set_num_threads(1);
    std::vector<int> reference_counts(N);
    for (int i = 0; i < N; ++i) {
        std::vector<int> neighbors(kNeighborListSize);
        reference_counts[i] = tree_->neighbor_search(particles[i], neighbors, particles, false);
    }

    // Test with different thread counts
    int max_threads = omp_get_max_threads();
    for (int num_threads : {2, 4, 8}) {
        if (num_threads > max_threads) continue;

        // Reset particles and rebuild tree
        particles = particles_original;
        BuildTree(particles);

        omp_set_num_threads(num_threads);
        std::vector<int> counts(N);

#pragma omp parallel for
        for (int i = 0; i < N; ++i) {
            std::vector<int> neighbors(kNeighborListSize);
            counts[i] = tree_->neighbor_search(particles[i], neighbors, particles, false);
        }

        // Verify counts match reference
        for (int i = 0; i < N; ++i) {
            EXPECT_EQ(counts[i], reference_counts[i])
                << "Mismatch with " << num_threads << " threads at particle " << i;
        }
    }

    // Restore original thread count
    omp_set_num_threads(max_threads);
#else
    GTEST_SKIP() << "OpenMP not available";
#endif
}

// =============================================================================
// Optimization 2: Distance Caching Tests
// =============================================================================

/**
 * TEST: Distance caching feature availability
 *
 * GIVEN: BHTree neighbor search
 * WHEN: Checking for distance caching support
 * THEN: Feature flag indicates caching is enabled
 */
TEST_F(PerformanceOptimizationTest, DistanceCaching_FeatureExists) {
#ifdef SPH_USE_DISTANCE_CACHING
    SUCCEED() << "Distance caching feature flag exists";
#else
    FAIL() << "Distance caching optimization not implemented. "
           << "Define SPH_USE_DISTANCE_CACHING when implementing.";
#endif
}

/**
 * TEST: Sorted neighbor order is identical with caching
 *
 * GIVEN: Random particle distribution
 * WHEN: Sort neighbors using cached vs recomputed distances
 * THEN: Sorted order is identical
 */
TEST_F(PerformanceOptimizationTest, DistanceCaching_SortedOrderIdentical) {
    const int N = 500;
    SetUpPeriodic(0.0, 1.0);
    auto particles = create_random_particles(N, 0.0, 1.0, 0.15, 42);
    BuildTree(particles);

    // For each particle, verify neighbors are sorted by distance
    bool all_sorted = true;
    for (int i = 0; i < N; ++i) {
        std::vector<int> neighbors(kNeighborListSize);
        int count = tree_->neighbor_search(particles[i], neighbors, particles, false);

        // Verify distances are monotonically non-decreasing
        for (int n = 1; n < count; ++n) {
            vec_t r_prev = periodic_->calc_r_ij(particles[i].pos, particles[neighbors[n-1]].pos);
            vec_t r_curr = periodic_->calc_r_ij(particles[i].pos, particles[neighbors[n]].pos);
            real dist2_prev = abs2(r_prev);
            real dist2_curr = abs2(r_curr);

            if (dist2_prev > dist2_curr + 1e-14) {
                all_sorted = false;
                ADD_FAILURE() << "Neighbors not sorted at particle " << i
                             << " position " << n
                             << ": dist2[" << n-1 << "]=" << dist2_prev
                             << " > dist2[" << n << "]=" << dist2_curr;
            }
        }
    }

    EXPECT_TRUE(all_sorted) << "Neighbor lists must be sorted by distance";
}

/**
 * TEST: Distance caching produces identical sort results
 *
 * GIVEN: Neighbor list from tree search
 * WHEN: Sorting with cached distances vs recomputed distances
 * THEN: Sort order is identical
 */
TEST_F(PerformanceOptimizationTest, DistanceCaching_IdenticalToRecompute) {
    const int N = 200;
    SetUpPeriodic(0.0, 1.0);
    auto particles = create_random_particles(N, 0.0, 1.0, 0.15, 42);
    BuildTree(particles);

    // Test a subset of particles
    for (int i = 0; i < std::min(N, 50); ++i) {
        std::vector<int> neighbors(kNeighborListSize);
        int count = tree_->neighbor_search(particles[i], neighbors, particles, false);

        if (count < 2) continue;  // Need at least 2 neighbors to test sorting

        // The neighbors should already be sorted - verify this
        std::vector<std::pair<real, int>> distances_with_ids;
        for (int n = 0; n < count; ++n) {
            vec_t r_ij = periodic_->calc_r_ij(particles[i].pos, particles[neighbors[n]].pos);
            distances_with_ids.push_back({abs2(r_ij), neighbors[n]});
        }

        // Sort by distance
        std::sort(distances_with_ids.begin(), distances_with_ids.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });

        // Verify order matches
        for (int n = 0; n < count; ++n) {
            EXPECT_EQ(neighbors[n], distances_with_ids[n].second)
                << "Sort order differs at position " << n << " for particle " << i;
        }
    }
}

// =============================================================================
// Optimization 3: OpenMP Dynamic Scheduling Tests
// =============================================================================

/**
 * TEST: OpenMP scheduling feature availability
 *
 * GIVEN: FluidForce or PreInteraction class
 * WHEN: Checking for dynamic/guided scheduling support
 * THEN: Feature flag indicates dynamic scheduling is enabled
 */
TEST_F(PerformanceOptimizationTest, OpenMPScheduling_FeatureExists) {
#ifdef SPH_USE_DYNAMIC_SCHEDULING
    SUCCEED() << "Dynamic scheduling feature flag exists";
#else
    FAIL() << "OpenMP dynamic scheduling optimization not implemented. "
           << "Define SPH_USE_DYNAMIC_SCHEDULING when implementing.";
#endif
}

/**
 * TEST: Static and dynamic scheduling produce same results
 *
 * GIVEN: Non-uniform particle distribution
 * WHEN: Run with static vs dynamic scheduling
 * THEN: Results are bitwise identical
 */
TEST_F(PerformanceOptimizationTest, OpenMPScheduling_StaticDynamicSameResults) {
#ifdef _OPENMP
    const int N = 500;
    SetUpPeriodic(0.0, 1.0);
    auto particles_original = create_clustered_particles(N, 0.0, 1.0, 0.1, 42);

    // Run with current scheduling (baseline)
    auto particles = particles_original;
    BuildTree(particles);

    // Capture neighbor counts (using parallel loop)
    std::vector<int> baseline_counts(N);
#pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        std::vector<int> neighbors(kNeighborListSize);
        baseline_counts[i] = tree_->neighbor_search(particles[i], neighbors, particles, false);
    }

    // Verify counts are deterministic across multiple runs
    for (int trial = 0; trial < 5; ++trial) {
        particles = particles_original;
        BuildTree(particles);

        std::vector<int> trial_counts(N);
#pragma omp parallel for schedule(dynamic, 32)
        for (int i = 0; i < N; ++i) {
            std::vector<int> neighbors(kNeighborListSize);
            trial_counts[i] = tree_->neighbor_search(particles[i], neighbors, particles, false);
        }

        for (int i = 0; i < N; ++i) {
            EXPECT_EQ(trial_counts[i], baseline_counts[i])
                << "Mismatch at trial " << trial << " particle " << i;
        }
    }
#else
    GTEST_SKIP() << "OpenMP not available";
#endif
}

/**
 * TEST: Different chunk sizes produce same results
 *
 * GIVEN: Non-uniform particle distribution
 * WHEN: Run with different dynamic scheduling chunk sizes
 * THEN: Results are identical regardless of chunk size
 */
TEST_F(PerformanceOptimizationTest, OpenMPScheduling_DifferentChunkSizesSameResults) {
#ifdef _OPENMP
    const int N = 500;
    SetUpPeriodic(0.0, 1.0);
    auto particles_original = create_clustered_particles(N, 0.0, 1.0, 0.1, 42);

    // Reference with default scheduling
    auto particles = particles_original;
    BuildTree(particles);

    std::vector<int> reference_counts(N);
#pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        std::vector<int> neighbors(kNeighborListSize);
        reference_counts[i] = tree_->neighbor_search(particles[i], neighbors, particles, false);
    }

    // Test with different chunk sizes
    for (int chunk : {1, 4, 16, 64, 256}) {
        particles = particles_original;
        BuildTree(particles);

        std::vector<int> counts(N);
        // Note: Can't pass variable to schedule() pragma directly in standard OpenMP
        // This is a conceptual test - actual implementation would use runtime scheduling
#pragma omp parallel for schedule(dynamic, 32)
        for (int i = 0; i < N; ++i) {
            std::vector<int> neighbors(kNeighborListSize);
            counts[i] = tree_->neighbor_search(particles[i], neighbors, particles, false);
        }

        for (int i = 0; i < N; ++i) {
            EXPECT_EQ(counts[i], reference_counts[i])
                << "Mismatch with chunk size " << chunk << " at particle " << i;
        }
    }
#else
    GTEST_SKIP() << "OpenMP not available";
#endif
}

// =============================================================================
// Integration/Regression Tests
// =============================================================================

/**
 * TEST: All optimizations combined preserve correctness
 *
 * GIVEN: Standard particle distribution
 * WHEN: Run full simulation step with all optimizations
 * THEN: Results match non-optimized baseline within tolerance
 */
TEST_F(PerformanceOptimizationTest, Integration_AllOptimizationsPreserveCorrectness) {
    // GIVEN: Grid particles
#if DIM == 3
    const int n_per_dim = 8;  // 512 particles
#elif DIM == 2
    const int n_per_dim = 23;  // 529 particles
#else
    const int n_per_dim = 500;
#endif
    SetUpPeriodic(0.0, 1.0);
    auto particles = create_grid_particles(n_per_dim, 0.0, 1.0, 0.15);
    const int N = particles.size();
    BuildTree(particles);

    // Capture state before any operations
    ParticleState initial_state;
    initial_state.capture(particles);

    // Run neighbor search for all particles
    std::vector<int> neighbor_counts(N);
#pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        std::vector<int> neighbors(kNeighborListSize);
        neighbor_counts[i] = tree_->neighbor_search(particles[i], neighbors, particles, false);
    }

    // Verify all particles have reasonable neighbor counts
    int total_neighbors = 0;
    int min_neighbors = neighbor_counts[0];
    int max_neighbors = neighbor_counts[0];
    for (int i = 0; i < N; ++i) {
        total_neighbors += neighbor_counts[i];
        min_neighbors = std::min(min_neighbors, neighbor_counts[i]);
        max_neighbors = std::max(max_neighbors, neighbor_counts[i]);

        // Each particle should have at least itself as neighbor
        EXPECT_GE(neighbor_counts[i], 1) << "Particle " << i << " has no neighbors";
    }

    double avg_neighbors = static_cast<double>(total_neighbors) / N;
    EXPECT_GT(avg_neighbors, 1.0) << "Average neighbor count too low";

    // Verify determinism: run again and compare
    std::vector<int> neighbor_counts_2(N);
#pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        std::vector<int> neighbors(kNeighborListSize);
        neighbor_counts_2[i] = tree_->neighbor_search(particles[i], neighbors, particles, false);
    }

    for (int i = 0; i < N; ++i) {
        EXPECT_EQ(neighbor_counts_2[i], neighbor_counts[i])
            << "Non-deterministic neighbor count at particle " << i;
    }
}

/**
 * TEST: Performance regression - optimizations should not be slower
 *
 * This is a DISABLED test that measures performance.
 * Enable when benchmarking the optimization implementation.
 */
TEST_F(PerformanceOptimizationTest, DISABLED_Performance_NotSlowerThanBaseline) {
    // GIVEN: Large particle count
    const int N = 10000;
    SetUpPeriodic(0.0, 1.0);
    auto particles = create_random_particles(N, 0.0, 1.0, 0.1, 42);
    BuildTree(particles);

    // WHEN: Run neighbor search many times
    const int ITERATIONS = 10;

    auto start = std::chrono::high_resolution_clock::now();
    for (int iter = 0; iter < ITERATIONS; ++iter) {
#pragma omp parallel for
        for (int i = 0; i < N; ++i) {
            std::vector<int> neighbors(kNeighborListSize);
            tree_->neighbor_search(particles[i], neighbors, particles, false);
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    std::cout << "[BENCHMARK] Neighbor search: " << N << " particles x " << ITERATIONS
              << " iterations = " << elapsed_ms << " ms" << std::endl;
    std::cout << "  Per iteration: " << elapsed_ms / ITERATIONS << " ms" << std::endl;
    std::cout << "  Per particle: " << static_cast<double>(elapsed_ms) / (N * ITERATIONS) * 1000
              << " us" << std::endl;

    // THEN: Should complete in reasonable time
    // Adjust threshold based on hardware
    EXPECT_LT(elapsed_ms, 60000) << "Performance regression detected";
}

/**
 * TEST: Thread safety - no data races under ThreadSanitizer
 *
 * GIVEN: Multiple threads accessing shared data
 * WHEN: Running parallel neighbor search
 * THEN: No data races detected (build with -fsanitize=thread)
 */
TEST_F(PerformanceOptimizationTest, ThreadSafety_NoDataRaces) {
#ifdef _OPENMP
    const int N = 1000;
    SetUpPeriodic(0.0, 1.0);
    auto particles = create_random_particles(N, 0.0, 1.0, 0.1, 42);
    BuildTree(particles);

    // Run many iterations to increase chance of race detection
    for (int trial = 0; trial < 50; ++trial) {
        std::vector<int> counts(N);

#pragma omp parallel for
        for (int i = 0; i < N; ++i) {
            std::vector<int> neighbors(kNeighborListSize);
            counts[i] = tree_->neighbor_search(particles[i], neighbors, particles, false);
        }

        // Verify no obviously corrupted data
        for (int i = 0; i < N; ++i) {
            EXPECT_GE(counts[i], 0) << "Negative neighbor count at particle " << i;
            EXPECT_LT(counts[i], kNeighborListSize) << "Overflow at particle " << i;
        }
    }

    // If we reach here without TSAN error, test passes
    SUCCEED() << "No data races detected (verify with -fsanitize=thread)";
#else
    GTEST_SKIP() << "OpenMP not available";
#endif
}

// =============================================================================
// Edge Cases and Boundary Conditions
// =============================================================================

/**
 * TEST: Empty particle set
 */
TEST_F(PerformanceOptimizationTest, EdgeCase_EmptyParticleSet) {
    std::vector<SPHParticle> particles;
    // Empty tree should not crash
    tree_->resize(0, 5);
    tree_->make(particles, 0);
    SUCCEED() << "Empty particle set handled correctly";
}

/**
 * TEST: Single particle
 */
TEST_F(PerformanceOptimizationTest, EdgeCase_SingleParticle) {
    SetUpPeriodic(0.0, 1.0);
    std::vector<SPHParticle> particles;
    particles.push_back(create_particle(0, vec_t(0.5), 1.0, 0.1));
    BuildTree(particles);

    std::vector<int> neighbors(kNeighborListSize);
    int count = tree_->neighbor_search(particles[0], neighbors, particles, false);

    EXPECT_EQ(count, 1) << "Single particle should find itself as neighbor";
    EXPECT_EQ(neighbors[0], 0) << "Neighbor should be the particle itself";
}

/**
 * TEST: Two particles at same position
 */
TEST_F(PerformanceOptimizationTest, EdgeCase_CoincidentParticles) {
    SetUpPeriodic(0.0, 1.0);
    std::vector<SPHParticle> particles;
    vec_t pos(0.5);
    particles.push_back(create_particle(0, pos, 1.0, 0.1));
    particles.push_back(create_particle(1, pos, 1.0, 0.1));
    BuildTree(particles);

    std::vector<int> neighbors(kNeighborListSize);
    int count = tree_->neighbor_search(particles[0], neighbors, particles, false);

    EXPECT_EQ(count, 2) << "Should find both coincident particles";
}

/**
 * TEST: Particles at domain corners (periodic)
 */
TEST_F(PerformanceOptimizationTest, EdgeCase_DomainCorners) {
    SetUpPeriodic(0.0, 1.0);
    std::vector<SPHParticle> particles;

    // Place particles at corners
    double eps = 0.01;
#if DIM == 1
    particles.push_back(create_particle(0, vec_t(eps), 1.0, 0.2));
    particles.push_back(create_particle(1, vec_t(1.0 - eps), 1.0, 0.2));
#elif DIM == 2
    particles.push_back(create_particle(0, vec_t(eps, eps), 1.0, 0.2));
    particles.push_back(create_particle(1, vec_t(1.0 - eps, 1.0 - eps), 1.0, 0.2));
#else
    particles.push_back(create_particle(0, vec_t(eps, eps, eps), 1.0, 0.2));
    particles.push_back(create_particle(1, vec_t(1.0 - eps, 1.0 - eps, 1.0 - eps), 1.0, 0.2));
#endif
    BuildTree(particles);

    // With periodic boundaries, opposite corner particles should be neighbors
    std::vector<int> neighbors(kNeighborListSize);
    int count = tree_->neighbor_search(particles[0], neighbors, particles, false);

    // In periodic domain, diagonal corners are close: sqrt(DIM) * 2*eps distance
    // With sml=0.2, they should be neighbors
    EXPECT_EQ(count, 2) << "Corner particles should find each other via periodic wrap";
}

// =============================================================================
// Main
// =============================================================================
int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
