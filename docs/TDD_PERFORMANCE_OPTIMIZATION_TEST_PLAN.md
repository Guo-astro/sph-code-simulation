# TDD Test Plan: SPH Performance Optimizations

**Date**: 2026-01-22
**Author**: Generated for Claude Code session
**Target Files**: `fluid_force.cpp`, `pre_interaction.cpp`, `bhtree.cpp`
**Test Framework**: Google Test (gtest)

---

## Overview

This document defines TDD test cases for three performance optimizations in the SPH code:

1. **Thread-local neighbor list** - Avoid heap allocation in hot loops
2. **Distance caching before sorting** - Reduce redundant computation in tree neighbor search
3. **OpenMP dynamic scheduling** - Better load balancing for variable particle density

---

## Optimization 1: Thread-Local Neighbor List

### Problem Statement

Current code allocates a `std::vector<int>` inside each iteration of `#pragma omp parallel for`:

```cpp
// fluid_force.cpp:34
#pragma omp parallel for
for(int i = 0; i < num; ++i) {
    std::vector<int> neighbor_list(m_neighbor_number * neighbor_list_size);  // HEAP ALLOCATION
    // ...
}
```

This causes:
- ~10,000 malloc/free calls per timestep (for 10K particles)
- Memory allocator contention between threads
- Cache pollution from scattered allocations

### Solution

Use thread-local storage to reuse pre-allocated buffers:

```cpp
// Option A: OpenMP threadprivate
static std::vector<int> tls_neighbor_list;
#pragma omp threadprivate(tls_neighbor_list)

// Option B: Explicit thread-local with resize-once pattern
thread_local std::vector<int> neighbor_list;
neighbor_list.resize(max_neighbors);  // No-op after first call
```

### Test Cases

#### 1.1 Correctness Tests

```cpp
// File: tests/test_thread_local_neighbor_list.cpp

class ThreadLocalNeighborListTest : public ::testing::Test {
protected:
    std::shared_ptr<SPHParameters> param_;
    std::unique_ptr<BHTree> tree_;
    std::unique_ptr<Periodic> periodic_;
    std::unique_ptr<FluidForce> fluid_force_;
    std::unique_ptr<PreInteraction> pre_interaction_;

    void SetUp() override;
    void SetUpParticles(std::vector<SPHParticle>& particles, int n, unsigned seed = 42);
    void CaptureState(const std::vector<SPHParticle>& particles,
                      std::vector<vec_t>& acc, std::vector<real>& dene);
};

TEST_F(ThreadLocalNeighborListTest, FluidForceOutputMatchesBaseline) {
    /**
     * GIVEN: Random particle distribution with known baseline results
     * WHEN: Run FluidForce::calculation() with thread-local neighbor list
     * THEN: Particle acc and dene match baseline within floating-point epsilon
     *
     * Rationale: Optimization must not change numerical results
     */

    const int N = 1000;
    auto particles = CreateRandomParticles(N);

    // Store baseline state (before optimization)
    std::vector<vec_t> baseline_acc(N);
    std::vector<real> baseline_dene(N);
    CaptureFluidForceOutput(particles, baseline_acc, baseline_dene);

    // Reset particles to initial state
    ResetParticles(particles);

    // Run with thread-local optimization
    fluid_force_->calculation(sim_);  // Uses thread-local

    // Compare
    for (int i = 0; i < N; ++i) {
        for (int d = 0; d < DIM; ++d) {
            EXPECT_NEAR(particles[i].acc[d], baseline_acc[i][d], 1e-14)
                << "acc mismatch at particle " << i << " dim " << d;
        }
        EXPECT_NEAR(particles[i].dene, baseline_dene[i], 1e-14)
            << "dene mismatch at particle " << i;
    }
}

TEST_F(ThreadLocalNeighborListTest, PreInteractionDensityMatchesBaseline) {
    /**
     * GIVEN: Random particle distribution
     * WHEN: Run PreInteraction::calculation() with thread-local neighbor list
     * THEN: Particle dens, pres, gradh match baseline
     */

    const int N = 1000;
    auto particles = CreateRandomParticles(N);

    std::vector<real> baseline_dens(N), baseline_pres(N), baseline_gradh(N);
    // ... capture baseline ...

    // Reset and run with optimization
    pre_interaction_->calculation(sim_);

    for (int i = 0; i < N; ++i) {
        EXPECT_NEAR(particles[i].dens, baseline_dens[i], 1e-14);
        EXPECT_NEAR(particles[i].pres, baseline_pres[i], 1e-14);
        EXPECT_NEAR(particles[i].gradh, baseline_gradh[i], 1e-14);
    }
}

TEST_F(ThreadLocalNeighborListTest, MultipleTimestepsAccumulateCorrectly) {
    /**
     * GIVEN: Particles evolved over multiple timesteps
     * WHEN: Thread-local buffers are reused across timesteps
     * THEN: No stale data from previous timesteps affects results
     *
     * Rationale: Buffer reuse must clear/overwrite properly
     */

    const int N = 500;
    const int TIMESTEPS = 10;
    auto particles = CreateRandomParticles(N);

    for (int t = 0; t < TIMESTEPS; ++t) {
        // Capture expected state
        auto expected_particles = particles;
        RunBaselineFluidForce(expected_particles);

        // Run optimized
        fluid_force_->calculation(sim_);

        // Compare
        for (int i = 0; i < N; ++i) {
            EXPECT_NEAR(particles[i].dene, expected_particles[i].dene, 1e-14)
                << "Mismatch at timestep " << t << " particle " << i;
        }

        // Evolve particles
        IntegrateOneTimestep(particles);
    }
}
```

#### 1.2 Thread Safety Tests

```cpp
TEST_F(ThreadLocalNeighborListTest, ThreadSafety_NoDataRaces) {
    /**
     * GIVEN: Multiple threads executing concurrently
     * WHEN: Each thread uses its own thread-local neighbor list
     * THEN: No data races occur (detected by ThreadSanitizer)
     *
     * Build with: -fsanitize=thread
     * Run with: TSAN_OPTIONS=detect_deadlocks=0
     */

    const int N = 2000;  // Enough to saturate threads
    auto particles = CreateRandomParticles(N);

    // Run 100 times to increase chance of race detection
    for (int trial = 0; trial < 100; ++trial) {
        fluid_force_->calculation(sim_);
        pre_interaction_->calculation(sim_);
    }

    // If we reach here without TSAN error, test passes
    SUCCEED();
}

TEST_F(ThreadLocalNeighborListTest, ThreadSafety_BufferIsolation) {
    /**
     * GIVEN: Thread-local buffers initialized per-thread
     * WHEN: Threads write different neighbor counts to their buffers
     * THEN: No cross-thread contamination occurs
     *
     * Test strategy: Use deterministic particle layout where
     * each thread processes particles with known neighbor counts.
     * Verify counts match expected values.
     */

    const int N = 1000;
    const int num_threads = omp_get_max_threads();
    std::vector<int> expected_neighbor_counts(N);

    // Create particles in grid with known neighbor structure
    auto particles = CreateGridParticles(10, 10, 10);  // 1000 particles

    // Pre-compute expected neighbor counts using exhaustive search
    for (int i = 0; i < N; ++i) {
        expected_neighbor_counts[i] = ExhaustiveNeighborCount(particles, i);
    }

    // Run parallel computation
    std::vector<int> actual_neighbor_counts(N, -1);

    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        // This uses thread-local neighbor list internally
        actual_neighbor_counts[i] = tree_->neighbor_search(...);
    }

    for (int i = 0; i < N; ++i) {
        EXPECT_EQ(actual_neighbor_counts[i], expected_neighbor_counts[i])
            << "Neighbor count mismatch at particle " << i;
    }
}

TEST_F(ThreadLocalNeighborListTest, ThreadSafety_DifferentThreadCounts) {
    /**
     * GIVEN: Same particle configuration
     * WHEN: Run with 1, 2, 4, 8 threads (up to max)
     * THEN: Results are identical regardless of thread count
     *
     * Rationale: Thread-local storage must produce deterministic results
     */

    const int N = 1000;
    auto particles = CreateRandomParticles(N);

    // Reference: single-threaded execution
    omp_set_num_threads(1);
    fluid_force_->calculation(sim_);

    std::vector<vec_t> reference_acc(N);
    std::vector<real> reference_dene(N);
    CaptureState(particles, reference_acc, reference_dene);

    // Test with different thread counts
    for (int num_threads : {2, 4, 8}) {
        if (num_threads > omp_get_max_threads()) continue;

        ResetParticles(particles);
        omp_set_num_threads(num_threads);
        fluid_force_->calculation(sim_);

        for (int i = 0; i < N; ++i) {
            for (int d = 0; d < DIM; ++d) {
                EXPECT_NEAR(particles[i].acc[d], reference_acc[i][d], 1e-14)
                    << "Mismatch with " << num_threads << " threads";
            }
        }
    }
}
```

#### 1.3 Performance Tests

```cpp
TEST_F(ThreadLocalNeighborListTest, DISABLED_Performance_NoMallocInHotLoop) {
    /**
     * GIVEN: Large particle count
     * WHEN: Run multiple timesteps
     * THEN: Thread-local version is faster than per-iteration allocation
     *
     * Measurement: Wall clock time for 10 timesteps
     * Expected speedup: 10-30% for malloc-heavy workloads
     */

    const int N = 10000;
    const int TIMESTEPS = 10;
    auto particles = CreateRandomParticles(N);

    // Benchmark: per-iteration allocation (baseline)
    auto start_baseline = std::chrono::high_resolution_clock::now();
    for (int t = 0; t < TIMESTEPS; ++t) {
        fluid_force_->calculation_baseline(sim_);  // Original code path
    }
    auto end_baseline = std::chrono::high_resolution_clock::now();
    auto baseline_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        end_baseline - start_baseline).count();

    // Benchmark: thread-local allocation
    ResetParticles(particles);
    auto start_optimized = std::chrono::high_resolution_clock::now();
    for (int t = 0; t < TIMESTEPS; ++t) {
        fluid_force_->calculation(sim_);  // Optimized code path
    }
    auto end_optimized = std::chrono::high_resolution_clock::now();
    auto optimized_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        end_optimized - start_optimized).count();

    std::cout << "[BENCHMARK] Thread-local neighbor list:" << std::endl;
    std::cout << "  Baseline:  " << baseline_ms << " ms" << std::endl;
    std::cout << "  Optimized: " << optimized_ms << " ms" << std::endl;
    std::cout << "  Speedup:   " << static_cast<double>(baseline_ms) / optimized_ms << "x" << std::endl;

    // Expect at least 10% speedup
    EXPECT_LT(optimized_ms, baseline_ms * 0.9)
        << "Thread-local should be at least 10% faster";
}

TEST_F(ThreadLocalNeighborListTest, DISABLED_Performance_ScalingWithParticleCount) {
    /**
     * Measure performance scaling for N = 1000, 5000, 10000, 50000
     * Verify linear or better scaling with thread-local optimization
     */

    std::cout << "[BENCHMARK] Scaling with particle count:" << std::endl;
    std::cout << "N\tBaseline(ms)\tOptimized(ms)\tSpeedup" << std::endl;

    for (int N : {1000, 5000, 10000, 50000}) {
        auto particles = CreateRandomParticles(N);

        // Benchmark baseline
        auto start = std::chrono::high_resolution_clock::now();
        fluid_force_->calculation_baseline(sim_);
        auto baseline_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - start).count();

        // Benchmark optimized
        ResetParticles(particles);
        start = std::chrono::high_resolution_clock::now();
        fluid_force_->calculation(sim_);
        auto optimized_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - start).count();

        std::cout << N << "\t" << baseline_ms << "\t\t" << optimized_ms << "\t\t"
                  << (optimized_ms > 0 ? static_cast<double>(baseline_ms) / optimized_ms : 0)
                  << std::endl;
    }
}
```

---

## Optimization 2: Distance Caching Before Sorting

### Problem Statement

Current `BHTree::neighbor_search()` sorts neighbors by distance using a lambda that recomputes distance for every comparison:

```cpp
// bhtree.cpp:141-145
std::sort(neighbor_list.begin(), neighbor_list.begin() + n_neighbor, [&](const int a, const int b) {
    const vec_t r_ia = m_periodic->calc_r_ij(pos_i, particles[a].pos);  // Computed ~O(N log N) times!
    const vec_t r_ib = m_periodic->calc_r_ij(pos_i, particles[b].pos);
    return abs2(r_ia) < abs2(r_ib);
});
```

For N neighbors, `std::sort` performs O(N log N) comparisons, so distance is computed O(N log N) times instead of just N times.

### Solution

Cache distances before sorting:

```cpp
struct NeighborWithDist {
    int id;
    real dist2;
};
std::vector<NeighborWithDist> neighbors_with_dist(n_neighbor);
for (int i = 0; i < n_neighbor; ++i) {
    neighbors_with_dist[i].id = neighbor_list[i];
    neighbors_with_dist[i].dist2 = abs2(m_periodic->calc_r_ij(pos_i, particles[neighbor_list[i]].pos));
}
std::sort(neighbors_with_dist.begin(), neighbors_with_dist.end(),
          [](const auto& a, const auto& b) { return a.dist2 < b.dist2; });
```

### Test Cases

#### 2.1 Correctness Tests

```cpp
// File: tests/test_distance_caching.cpp

class DistanceCachingTest : public ::testing::Test {
protected:
    // ... setup ...
};

TEST_F(DistanceCachingTest, SortedOrderIdentical) {
    /**
     * GIVEN: Random particle distribution
     * WHEN: Sort neighbors using cached vs recomputed distances
     * THEN: Sorted order is identical
     */

    const int N = 1000;
    auto particles = CreateRandomParticles(N);
    BuildTree(particles);

    for (auto& p : particles) {
        std::vector<int> neighbors_baseline(5000);
        std::vector<int> neighbors_cached(5000);

        int count_baseline = tree_baseline_->neighbor_search(p, neighbors_baseline, particles, false);
        int count_cached = tree_cached_->neighbor_search(p, neighbors_cached, particles, false);

        ASSERT_EQ(count_baseline, count_cached);

        for (int i = 0; i < count_baseline; ++i) {
            EXPECT_EQ(neighbors_baseline[i], neighbors_cached[i])
                << "Sort order differs at position " << i << " for particle " << p.id;
        }
    }
}

TEST_F(DistanceCachingTest, DistancesAreSorted) {
    /**
     * GIVEN: Neighbors returned from search
     * WHEN: Check distance ordering
     * THEN: Distances are monotonically non-decreasing
     */

    const int N = 500;
    auto particles = CreateRandomParticles(N);
    BuildTree(particles);

    for (auto& p : particles) {
        std::vector<int> neighbors(5000);
        int count = tree_->neighbor_search(p, neighbors, particles, false);

        for (int i = 1; i < count; ++i) {
            real dist_prev = abs2(periodic_->calc_r_ij(p.pos, particles[neighbors[i-1]].pos));
            real dist_curr = abs2(periodic_->calc_r_ij(p.pos, particles[neighbors[i]].pos));

            EXPECT_LE(dist_prev, dist_curr + 1e-14)
                << "Neighbors not sorted at position " << i << " for particle " << p.id;
        }
    }
}

TEST_F(DistanceCachingTest, PeriodicBoundaryHandledCorrectly) {
    /**
     * GIVEN: Particles near periodic boundaries
     * WHEN: Distance caching with periodic corrections
     * THEN: Distances are computed correctly across boundaries
     */

    SetUpPeriodic(0.0, 1.0);
    std::vector<SPHParticle> particles;

    // Particle at x=0.05
    particles.push_back(CreateParticle(0, vec_t(0.05, 0.5, 0.5)));
    // Particle at x=0.95 (should be closer via periodic wrap)
    particles.push_back(CreateParticle(1, vec_t(0.95, 0.5, 0.5)));
    // Particle at x=0.15 (should be farther)
    particles.push_back(CreateParticle(2, vec_t(0.15, 0.5, 0.5)));

    BuildTree(particles);

    std::vector<int> neighbors(100);
    int count = tree_->neighbor_search(particles[0], neighbors, particles, false);

    // Verify order: self (0), then particle 1 (0.1 away via wrap), then particle 2 (0.1 away direct)
    ASSERT_GE(count, 3);
    EXPECT_EQ(neighbors[0], 0);  // Self
    // particle 1 is 0.1 away via periodic wrap: |0.05 - 0.95 + 1.0| = 0.1
    // particle 2 is 0.1 away direct: |0.05 - 0.15| = 0.1
    // Both are equidistant, so order may vary
}
```

#### 2.2 Performance Tests

```cpp
TEST_F(DistanceCachingTest, DISABLED_Performance_CachingFasterThanRecompute) {
    /**
     * GIVEN: Large particle count with many neighbors per particle
     * WHEN: Run neighbor search with sorting
     * THEN: Cached distance version is faster
     *
     * Expected speedup: 2-5x depending on average neighbor count
     */

    const int N = 10000;
    auto particles = CreateRandomParticles(N);
    BuildTree(particles);

    // Baseline: recompute distances during sort
    auto start = std::chrono::high_resolution_clock::now();
    for (auto& p : particles) {
        std::vector<int> neighbors(5000);
        tree_baseline_->neighbor_search(p, neighbors, particles, false);
    }
    auto baseline_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - start).count();

    // Optimized: cache distances before sort
    start = std::chrono::high_resolution_clock::now();
    for (auto& p : particles) {
        std::vector<int> neighbors(5000);
        tree_cached_->neighbor_search(p, neighbors, particles, false);
    }
    auto cached_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - start).count();

    std::cout << "[BENCHMARK] Distance caching:" << std::endl;
    std::cout << "  Recompute: " << baseline_ms << " ms" << std::endl;
    std::cout << "  Cached:    " << cached_ms << " ms" << std::endl;
    std::cout << "  Speedup:   " << static_cast<double>(baseline_ms) / cached_ms << "x" << std::endl;

    EXPECT_LT(cached_ms, baseline_ms * 0.8)
        << "Caching should be at least 20% faster";
}

TEST_F(DistanceCachingTest, DISABLED_Performance_ScalingWithNeighborCount) {
    /**
     * Verify speedup increases with average neighbor count
     * (more comparisons = more benefit from caching)
     */

    std::cout << "[BENCHMARK] Speedup vs neighbor count:" << std::endl;
    std::cout << "SML\tAvg_Neighbors\tBaseline(ms)\tCached(ms)\tSpeedup" << std::endl;

    for (double sml : {0.05, 0.10, 0.15, 0.20}) {
        const int N = 5000;
        auto particles = CreateRandomParticlesWithSml(N, sml);
        BuildTree(particles);

        // Count average neighbors
        int total_neighbors = 0;
        for (auto& p : particles) {
            std::vector<int> neighbors(5000);
            total_neighbors += tree_->neighbor_search(p, neighbors, particles, false);
        }
        double avg_neighbors = static_cast<double>(total_neighbors) / N;

        // Benchmark
        // ... timing code ...

        std::cout << sml << "\t" << avg_neighbors << "\t\t" << baseline_ms
                  << "\t\t" << cached_ms << "\t\t" << speedup << std::endl;
    }
}
```

---

## Optimization 3: OpenMP Dynamic Scheduling

### Problem Statement

Current code uses static scheduling:

```cpp
#pragma omp parallel for  // Default is static scheduling
for(int i = 0; i < num; ++i) {
    // Neighbor search - variable cost depending on local density
}
```

With variable particle density (e.g., shock tubes, collapse simulations), some particles have many more neighbors than others, causing load imbalance.

### Solution

Use dynamic scheduling with appropriate chunk size:

```cpp
#pragma omp parallel for schedule(dynamic, 64)
for(int i = 0; i < num; ++i) {
    // ...
}
```

Or guided scheduling for better cache behavior:

```cpp
#pragma omp parallel for schedule(guided)
```

### Test Cases

#### 3.1 Correctness Tests

```cpp
// File: tests/test_openmp_scheduling.cpp

class OpenMPSchedulingTest : public ::testing::Test {
protected:
    // ... setup ...
};

TEST_F(OpenMPSchedulingTest, StaticAndDynamicProduceSameResults) {
    /**
     * GIVEN: Non-uniform particle distribution
     * WHEN: Run with static vs dynamic scheduling
     * THEN: Results are bitwise identical
     */

    const int N = 5000;
    auto particles = CreateNonUniformParticles(N);  // Clustered distribution

    // Run with static scheduling
    std::vector<vec_t> acc_static(N);
    std::vector<real> dene_static(N);
    fluid_force_static_->calculation(sim_);
    CaptureState(particles, acc_static, dene_static);

    // Reset and run with dynamic scheduling
    ResetParticles(particles);
    fluid_force_dynamic_->calculation(sim_);

    // Compare
    for (int i = 0; i < N; ++i) {
        for (int d = 0; d < DIM; ++d) {
            EXPECT_EQ(particles[i].acc[d], acc_static[i][d])
                << "acc mismatch at particle " << i;
        }
        EXPECT_EQ(particles[i].dene, dene_static[i])
            << "dene mismatch at particle " << i;
    }
}

TEST_F(OpenMPSchedulingTest, GuidedSchedulingMatchesStatic) {
    /**
     * Same test for guided scheduling
     */
    // ... similar to above ...
}

TEST_F(OpenMPSchedulingTest, OrderIndependentReduction) {
    /**
     * GIVEN: Reduction operation (e.g., h_per_v_sig minimum)
     * WHEN: Run with different scheduling
     * THEN: Reduction result is identical
     *
     * Rationale: omp_real reduction must work correctly with any scheduling
     */

    const int N = 5000;
    auto particles = CreateNonUniformParticles(N);

    pre_interaction_static_->calculation(sim_);
    real h_per_v_sig_static = sim_->get_h_per_v_sig();

    ResetParticles(particles);
    pre_interaction_dynamic_->calculation(sim_);
    real h_per_v_sig_dynamic = sim_->get_h_per_v_sig();

    EXPECT_EQ(h_per_v_sig_static, h_per_v_sig_dynamic);
}
```

#### 3.2 Thread Safety Tests

```cpp
TEST_F(OpenMPSchedulingTest, ThreadSafety_DynamicSchedulingNoRaces) {
    /**
     * GIVEN: Dynamic scheduling with small chunk size
     * WHEN: Frequent work stealing between threads
     * THEN: No data races occur
     *
     * Build with: -fsanitize=thread
     */

    const int N = 10000;
    auto particles = CreateNonUniformParticles(N);

    for (int trial = 0; trial < 100; ++trial) {
        fluid_force_dynamic_->calculation(sim_);
    }

    SUCCEED();  // If we reach here without TSAN error
}

TEST_F(OpenMPSchedulingTest, ThreadSafety_NoFalseSharing) {
    /**
     * GIVEN: Particles stored contiguously in memory
     * WHEN: Adjacent particles processed by different threads
     * THEN: No false sharing degradation
     *
     * Test strategy: Compare performance with adjacent vs strided access
     */

    // This is more of a performance observation than a correctness test
    // False sharing would show up as slower dynamic scheduling
}
```

#### 3.3 Performance Tests

```cpp
TEST_F(OpenMPSchedulingTest, DISABLED_Performance_NonUniformDistribution) {
    /**
     * GIVEN: Non-uniform particle distribution (clustered)
     * WHEN: Run with static vs dynamic scheduling
     * THEN: Dynamic scheduling is faster due to better load balance
     */

    const int N = 50000;
    auto particles = CreateClusteredParticles(N);  // 90% in center, 10% sparse

    // Static scheduling
    auto start = std::chrono::high_resolution_clock::now();
    for (int t = 0; t < 10; ++t) {
        fluid_force_static_->calculation(sim_);
    }
    auto static_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - start).count();

    // Dynamic scheduling
    ResetParticles(particles);
    start = std::chrono::high_resolution_clock::now();
    for (int t = 0; t < 10; ++t) {
        fluid_force_dynamic_->calculation(sim_);
    }
    auto dynamic_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - start).count();

    std::cout << "[BENCHMARK] OpenMP scheduling (clustered distribution):" << std::endl;
    std::cout << "  Static:  " << static_ms << " ms" << std::endl;
    std::cout << "  Dynamic: " << dynamic_ms << " ms" << std::endl;
    std::cout << "  Speedup: " << static_cast<double>(static_ms) / dynamic_ms << "x" << std::endl;

    // For clustered distribution, dynamic should be faster
    EXPECT_LT(dynamic_ms, static_ms)
        << "Dynamic scheduling should be faster for non-uniform distribution";
}

TEST_F(OpenMPSchedulingTest, DISABLED_Performance_UniformDistribution) {
    /**
     * GIVEN: Uniform particle distribution
     * WHEN: Run with static vs dynamic scheduling
     * THEN: Performance is similar (static may be slightly faster due to no overhead)
     */

    const int N = 50000;
    auto particles = CreateGridParticles(N);  // Uniform grid

    // ... benchmark code ...

    // For uniform distribution, static should be at least as fast
    EXPECT_LE(static_ms, dynamic_ms * 1.1)
        << "Static should be at least as fast for uniform distribution";
}

TEST_F(OpenMPSchedulingTest, DISABLED_Performance_ChunkSizeTuning) {
    /**
     * Find optimal chunk size for dynamic scheduling
     */

    std::cout << "[BENCHMARK] Chunk size tuning:" << std::endl;
    std::cout << "Chunk\tTime(ms)" << std::endl;

    for (int chunk : {1, 4, 16, 64, 256, 1024}) {
        // Configure chunk size
        SetDynamicChunkSize(chunk);

        // Benchmark
        // ...

        std::cout << chunk << "\t" << time_ms << std::endl;
    }
}
```

---

## Integration Tests

### Regression Tests

```cpp
// File: tests/test_optimization_regression.cpp

class OptimizationRegressionTest : public ::testing::Test {
    /**
     * End-to-end regression tests that verify all optimizations
     * together do not change simulation results
     */
};

TEST_F(OptimizationRegressionTest, ShockTubeOutputUnchanged) {
    /**
     * GIVEN: Standard Sod shock tube initial conditions
     * WHEN: Run simulation with all optimizations enabled
     * THEN: Final state matches reference data within L2 tolerance
     *
     * Reference: tests/regression/reference_data/shock_tube_3d_*.json
     */

    // Load reference data
    auto reference = LoadReferenceData("shock_tube_3d_wendland_n120.json");

    // Run simulation with optimizations
    auto sim = SetupShockTube(/*with_optimizations=*/true);
    sim->run_to_time(reference.time);

    // Compare L2 errors
    double l2_rho = ComputeL2Error(sim->particles(), reference.particles(), &ParticleData::rho);
    double l2_pres = ComputeL2Error(sim->particles(), reference.particles(), &ParticleData::pres);
    double l2_ene = ComputeL2Error(sim->particles(), reference.particles(), &ParticleData::ene);

    // Use same threshold as existing regression tests
    constexpr double L2_THRESHOLD = 1e-10;
    EXPECT_LT(l2_rho, L2_THRESHOLD);
    EXPECT_LT(l2_pres, L2_THRESHOLD);
    EXPECT_LT(l2_ene, L2_THRESHOLD);
}

TEST_F(OptimizationRegressionTest, GravityCollapseUnchanged) {
    /**
     * GIVEN: Uniform sphere with self-gravity
     * WHEN: Run collapse simulation with optimizations
     * THEN: Collapse time and final density profile match reference
     */

    // ... similar structure ...
}
```

### Combined Performance Tests

```cpp
TEST_F(OptimizationRegressionTest, DISABLED_CombinedSpeedup) {
    /**
     * Measure cumulative speedup from all optimizations
     */

    const int N = 50000;
    const int TIMESTEPS = 100;

    // Baseline: no optimizations
    auto sim_baseline = SetupSimulation(/*optimizations=*/false);
    auto start = std::chrono::high_resolution_clock::now();
    for (int t = 0; t < TIMESTEPS; ++t) {
        sim_baseline->step();
    }
    auto baseline_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - start).count();

    // Optimized: all optimizations
    auto sim_optimized = SetupSimulation(/*optimizations=*/true);
    start = std::chrono::high_resolution_clock::now();
    for (int t = 0; t < TIMESTEPS; ++t) {
        sim_optimized->step();
    }
    auto optimized_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - start).count();

    std::cout << "[BENCHMARK] Combined optimizations:" << std::endl;
    std::cout << "  Baseline:  " << baseline_ms << " ms" << std::endl;
    std::cout << "  Optimized: " << optimized_ms << " ms" << std::endl;
    std::cout << "  Speedup:   " << static_cast<double>(baseline_ms) / optimized_ms << "x" << std::endl;

    // Expect meaningful combined speedup
    EXPECT_LT(optimized_ms, baseline_ms * 0.7)
        << "Combined optimizations should give at least 30% speedup";
}
```

---

## CMake Configuration

Add to `CMakeLists.txt`:

```cmake
# Performance optimization TDD tests
add_executable(test_performance_optimizations
    tests/test_thread_local_neighbor_list.cpp
    tests/test_distance_caching.cpp
    tests/test_openmp_scheduling.cpp
    tests/test_optimization_regression.cpp
    src/fluid_force.cpp
    src/pre_interaction.cpp
    src/bhtree.cpp
    src/exhaustive_search.cpp
    src/hernquist_katz_lookup_table.cpp
    tests/softening_lookup_table.cpp
)
target_compile_features(test_performance_optimizations PUBLIC cxx_std_14)
target_compile_definitions(test_performance_optimizations PRIVATE DIM=${SPH_DIM})
target_include_directories(test_performance_optimizations PUBLIC include ${HDF5_INCLUDE_DIRS})
target_link_libraries(test_performance_optimizations
    GTest::gtest_main
    GTest::gmock
    OpenMP::OpenMP_CXX
    Boost::boost
)
target_compile_options(test_performance_optimizations PRIVATE -Wall -Wno-sign-compare)

# Thread sanitizer build for race detection
option(ENABLE_THREAD_SANITIZER "Enable ThreadSanitizer for race detection tests" OFF)
if(ENABLE_THREAD_SANITIZER)
    target_compile_options(test_performance_optimizations PRIVATE -fsanitize=thread -g)
    target_link_options(test_performance_optimizations PRIVATE -fsanitize=thread)
endif()

gtest_discover_tests(test_performance_optimizations)
```

---

## Test Execution

```bash
# Build tests
cd build && cmake -DBUILD_TESTS=ON .. && make test_performance_optimizations

# Run all tests
./test_performance_optimizations

# Run correctness tests only (fast)
./test_performance_optimizations --gtest_filter="*Correctness*:*Matches*:*Identical*"

# Run thread safety tests with TSAN
cmake -DENABLE_THREAD_SANITIZER=ON .. && make test_performance_optimizations
TSAN_OPTIONS=detect_deadlocks=0 ./test_performance_optimizations --gtest_filter="*ThreadSafety*"

# Run performance benchmarks (disabled by default)
./test_performance_optimizations --gtest_also_run_disabled_tests --gtest_filter="*Performance*"
```

---

## Implementation Order

1. **Thread-local neighbor list** (highest impact)
   - Write tests first (RED)
   - Implement in `fluid_force.cpp` and `pre_interaction.cpp`
   - Run tests (GREEN)
   - Refactor if needed

2. **Distance caching** (medium impact)
   - Write tests first (RED)
   - Implement in `bhtree.cpp`
   - Run tests (GREEN)

3. **Dynamic scheduling** (variable impact, depends on workload)
   - Write tests first (RED)
   - Implement in `fluid_force.cpp` and `pre_interaction.cpp`
   - Tune chunk size based on benchmarks
   - Run tests (GREEN)

4. **Regression tests** (run after each optimization)
   - Verify shock tube results unchanged
   - Verify gravity collapse results unchanged

---

## Summary

| Optimization | Correctness Tests | Thread Safety Tests | Performance Tests |
|--------------|-------------------|--------------------|--------------------|
| Thread-local neighbor list | 3 | 3 | 2 |
| Distance caching | 3 | 0 | 2 |
| OpenMP dynamic scheduling | 3 | 2 | 3 |
| Integration/Regression | 2 | 0 | 1 |

**Total: 24 test cases**
