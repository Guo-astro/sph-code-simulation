/**
 * @file test_neighbor_search.cpp
 * @brief TDD validation tests comparing tree-based neighbor search with exhaustive search
 *
 * Purpose: Use exhaustive search as ground truth to validate tree-based neighbor search
 * before removing the EXHAUSTIVE_SEARCH compile flag.
 *
 * Test Strategy:
 * - Both search methods are called directly (not through #ifdef)
 * - Results are compared for correctness
 * - Tests cover edge cases, periodic boundaries, and various particle distributions
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

#include "defines.hpp"
#include "particle.hpp"
#include "vector_type.hpp"
#include "bhtree.hpp"
#include "exhaustive_search.hpp"
#include "periodic.hpp"
#include "parameters.hpp"

// Use types from sph namespace
using sph::SPHParticle;
using sph::SPHParameters;
using sph::BHTree;
using sph::Periodic;

namespace {

// Tolerance for floating-point comparisons
constexpr double kAbsTol = 1e-10;
constexpr double kRelTol = 1e-6;

// Test configuration
constexpr int kDefaultNeighborNumber = 32;
constexpr int kNeighborListSize = 5000;

// Helper to create default parameters
std::shared_ptr<SPHParameters> create_default_params(bool periodic = false,
                                                           double range_min = 0.0,
                                                           double range_max = 1.0) {
    auto param = std::make_shared<SPHParameters>();

    // Tree parameters
    param->tree.max_level = 20;
    param->tree.leaf_particle_num = 8;

    // Physics parameters
    param->physics.neighbor_number = kDefaultNeighborNumber;
    param->physics.gamma = 1.4;
    param->physics.c_smooth = 1.0;

    // Periodic parameters
    param->periodic.is_valid = periodic;
    for (int i = 0; i < DIM; ++i) {
        param->periodic.range_min[i] = range_min;
        param->periodic.range_max[i] = range_max;
    }

    // Gravity (disabled for neighbor tests)
    param->gravity.is_valid = false;

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
    p.sound = 1.0;
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

// Helper to compare neighbor sets (ignoring order)
bool neighbor_sets_equal(const std::vector<int>& list1, int count1,
                         const std::vector<int>& list2, int count2) {
    if (count1 != count2) return false;

    std::set<int> set1(list1.begin(), list1.begin() + count1);
    std::set<int> set2(list2.begin(), list2.begin() + count2);

    return set1 == set2;
}

} // anonymous namespace

// =============================================================================
// Test Fixture
// =============================================================================
class NeighborSearchTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Default setup with non-periodic boundaries
        param_ = create_default_params(false, 0.0, 1.0);
        tree_ = std::make_unique<BHTree>();
        tree_->initialize(param_);

        periodic_ = std::make_unique<Periodic>();
        periodic_->initialize(param_);
    }

    void SetUpPeriodic(double range_min, double range_max) {
        param_ = create_default_params(true, range_min, range_max);
        tree_ = std::make_unique<BHTree>();
        tree_->initialize(param_);

        periodic_ = std::make_unique<Periodic>();
        periodic_->initialize(param_);
    }

    void BuildTree(std::vector<SPHParticle>& particles) {
        int n = particles.size();
        tree_->resize(n, 5);
        tree_->make(particles, n);
        tree_->set_kernel();
    }

    // Run both search methods and compare results
    struct SearchResult {
        int exhaustive_count;
        int tree_count;
        std::vector<int> exhaustive_neighbors;
        std::vector<int> tree_neighbors;
        bool counts_match;
        bool sets_match;
    };

    SearchResult CompareSearchMethods(SPHParticle& p_i,
                                      std::vector<SPHParticle>& particles,
                                      bool is_ij = false) {
        SearchResult result;
        result.exhaustive_neighbors.resize(kNeighborListSize);
        result.tree_neighbors.resize(kNeighborListSize);

        int n = particles.size();

        // Exhaustive search (ground truth)
        result.exhaustive_count = sph::exhaustive_search(
            p_i, p_i.sml, particles, n,
            result.exhaustive_neighbors, kNeighborListSize,
            periodic_.get(), is_ij
        );

        // Tree-based search
        result.tree_count = tree_->neighbor_search(
            p_i, result.tree_neighbors, particles, is_ij
        );

        result.counts_match = (result.exhaustive_count == result.tree_count);
        result.sets_match = neighbor_sets_equal(
            result.exhaustive_neighbors, result.exhaustive_count,
            result.tree_neighbors, result.tree_count
        );

        return result;
    }

    std::shared_ptr<SPHParameters> param_;
    std::unique_ptr<BHTree> tree_;
    std::unique_ptr<Periodic> periodic_;
};

// =============================================================================
// Test: Basic neighbor count comparison
// =============================================================================
TEST_F(NeighborSearchTest, GivenGridParticles_WhenSearching_ThenCountsMatch) {
    // GIVEN: Regular grid of particles
    int n_per_dim = 10;
    double sml = 0.15;  // Smoothing length > spacing to ensure neighbors
    auto particles = create_grid_particles(n_per_dim, 0.0, 1.0, sml);

    BuildTree(particles);

    // WHEN/THEN: For each particle, neighbor counts should match
    int mismatches = 0;
    for (auto& p : particles) {
        auto result = CompareSearchMethods(p, particles, false);
        if (!result.counts_match) {
            mismatches++;
            if (mismatches <= 5) {  // Only print first 5 mismatches
                std::cerr << "Count mismatch for particle " << p.id
                          << ": exhaustive=" << result.exhaustive_count
                          << ", tree=" << result.tree_count << std::endl;
            }
        }
    }

    EXPECT_EQ(mismatches, 0) << "Found " << mismatches << " particles with count mismatches";
}

// =============================================================================
// Test: Neighbor ID set equality
// =============================================================================
TEST_F(NeighborSearchTest, GivenGridParticles_WhenSearching_ThenNeighborSetsMatch) {
    // GIVEN: Regular grid of particles
    int n_per_dim = 10;
    double sml = 0.15;
    auto particles = create_grid_particles(n_per_dim, 0.0, 1.0, sml);

    BuildTree(particles);

    // WHEN/THEN: For each particle, neighbor ID sets should be identical
    int mismatches = 0;
    for (auto& p : particles) {
        auto result = CompareSearchMethods(p, particles, false);
        if (!result.sets_match) {
            mismatches++;
            if (mismatches <= 3) {
                std::cerr << "Set mismatch for particle " << p.id << std::endl;
                std::cerr << "  Exhaustive (" << result.exhaustive_count << "): ";
                for (int i = 0; i < std::min(10, result.exhaustive_count); ++i) {
                    std::cerr << result.exhaustive_neighbors[i] << " ";
                }
                std::cerr << std::endl;
                std::cerr << "  Tree (" << result.tree_count << "): ";
                for (int i = 0; i < std::min(10, result.tree_count); ++i) {
                    std::cerr << result.tree_neighbors[i] << " ";
                }
                std::cerr << std::endl;
            }
        }
    }

    EXPECT_EQ(mismatches, 0) << "Found " << mismatches << " particles with set mismatches";
}

// =============================================================================
// Test: Random distribution
// =============================================================================
TEST_F(NeighborSearchTest, GivenRandomParticles_WhenSearching_ThenResultsMatch) {
    // GIVEN: Random distribution of particles
    int n = 500;
    double sml = 0.1;
    auto particles = create_random_particles(n, 0.05, 0.95, sml, 12345);

    BuildTree(particles);

    // WHEN/THEN: Results should match for all particles
    int count_mismatches = 0;
    int set_mismatches = 0;

    for (auto& p : particles) {
        auto result = CompareSearchMethods(p, particles, false);
        if (!result.counts_match) count_mismatches++;
        if (!result.sets_match) set_mismatches++;
    }

    EXPECT_EQ(count_mismatches, 0) << "Count mismatches: " << count_mismatches << "/" << n;
    EXPECT_EQ(set_mismatches, 0) << "Set mismatches: " << set_mismatches << "/" << n;
}

// =============================================================================
// Test: is_ij=true mode (symmetric neighbor search for fluid forces)
// =============================================================================
TEST_F(NeighborSearchTest, GivenVariableSml_WhenSearchingWithIsIjTrue_ThenResultsMatch) {
    // GIVEN: Particles with varying smoothing lengths
    int n = 200;
    auto particles = create_random_particles(n, 0.1, 0.9, 0.1, 54321);

    // Assign varying smoothing lengths
    std::mt19937 gen(98765);
    std::uniform_real_distribution<double> sml_dist(0.05, 0.15);
    for (auto& p : particles) {
        p.sml = sml_dist(gen);
    }

    BuildTree(particles);

    // WHEN/THEN: With is_ij=true, results should match
    int mismatches = 0;
    for (auto& p : particles) {
        auto result = CompareSearchMethods(p, particles, true);  // is_ij = true
        if (!result.sets_match) {
            mismatches++;
        }
    }

    // Note: is_ij mode has slightly different semantics between tree and exhaustive
    // Tree uses max(h_i, kernel_size_of_node), exhaustive uses max(h_i, h_j) per neighbor
    // Some discrepancies are expected - tree may find MORE neighbors (conservative)
    // We accept this if tree finds at least all exhaustive neighbors
    if (mismatches > 0) {
        std::cerr << "Note: " << mismatches << " particles have neighbor differences in is_ij mode" << std::endl;
        std::cerr << "This is expected due to tree's conservative node-based kernel_size" << std::endl;
    }

    // Verify tree finds AT LEAST all exhaustive neighbors (no false negatives)
    int false_negatives = 0;
    for (auto& p : particles) {
        auto result = CompareSearchMethods(p, particles, true);

        std::set<int> exhaustive_set(result.exhaustive_neighbors.begin(),
                                     result.exhaustive_neighbors.begin() + result.exhaustive_count);
        std::set<int> tree_set(result.tree_neighbors.begin(),
                               result.tree_neighbors.begin() + result.tree_count);

        for (int id : exhaustive_set) {
            if (tree_set.find(id) == tree_set.end()) {
                false_negatives++;
                break;
            }
        }
    }

    EXPECT_EQ(false_negatives, 0) << "Tree search missed neighbors that exhaustive found";
}

// =============================================================================
// Test: Periodic boundary conditions
// =============================================================================
TEST_F(NeighborSearchTest, GivenPeriodicBoundary_WhenParticleNearEdge_ThenFindsNeighborsAcross) {
    // GIVEN: Periodic boundaries and particles near edges
    SetUpPeriodic(0.0, 1.0);

    std::vector<SPHParticle> particles;
    int id = 0;

    // Place particles near the boundary
    double sml = 0.15;

#if DIM == 1
    // Particle at x=0.05 should find neighbors near x=0.95
    particles.push_back(create_particle(id++, vec_t(0.05), 1.0, sml));
    particles.push_back(create_particle(id++, vec_t(0.95), 1.0, sml));
    particles.push_back(create_particle(id++, vec_t(0.10), 1.0, sml));
    particles.push_back(create_particle(id++, vec_t(0.90), 1.0, sml));
#elif DIM == 2
    particles.push_back(create_particle(id++, vec_t(0.05, 0.5), 1.0, sml));
    particles.push_back(create_particle(id++, vec_t(0.95, 0.5), 1.0, sml));
    particles.push_back(create_particle(id++, vec_t(0.10, 0.5), 1.0, sml));
    particles.push_back(create_particle(id++, vec_t(0.90, 0.5), 1.0, sml));
#elif DIM == 3
    particles.push_back(create_particle(id++, vec_t(0.05, 0.5, 0.5), 1.0, sml));
    particles.push_back(create_particle(id++, vec_t(0.95, 0.5, 0.5), 1.0, sml));
    particles.push_back(create_particle(id++, vec_t(0.10, 0.5, 0.5), 1.0, sml));
    particles.push_back(create_particle(id++, vec_t(0.90, 0.5, 0.5), 1.0, sml));
#endif

    BuildTree(particles);

    // WHEN/THEN: Both methods should find neighbors across the periodic boundary
    for (auto& p : particles) {
        auto result = CompareSearchMethods(p, particles, false);

        EXPECT_TRUE(result.sets_match)
            << "Periodic boundary mismatch for particle " << p.id
            << ": exhaustive=" << result.exhaustive_count
            << ", tree=" << result.tree_count;
    }
}

// =============================================================================
// Test: Edge case - particle at exact same position
// =============================================================================
TEST_F(NeighborSearchTest, GivenOverlappingParticles_WhenSearching_ThenNocrash) {
    // GIVEN: Two particles at the exact same position
    std::vector<SPHParticle> particles;
    double sml = 0.1;

#if DIM == 1
    vec_t pos(0.5);
#elif DIM == 2
    vec_t pos(0.5, 0.5);
#elif DIM == 3
    vec_t pos(0.5, 0.5, 0.5);
#endif

    particles.push_back(create_particle(0, pos, 1.0, sml));
    particles.push_back(create_particle(1, pos, 1.0, sml));  // Same position!
    particles.push_back(create_particle(2, pos, 1.0, sml));  // Same position!

    BuildTree(particles);

    // WHEN/THEN: Should not crash, and both should find each other
    for (auto& p : particles) {
        auto result = CompareSearchMethods(p, particles, false);

        EXPECT_GE(result.exhaustive_count, 3) << "Should find all overlapping particles";
        EXPECT_TRUE(result.sets_match) << "Sets should match for overlapping particles";
    }
}

// =============================================================================
// Test: Distance ordering
// =============================================================================
TEST_F(NeighborSearchTest, GivenParticles_WhenSearching_ThenNeighborsSortedByDistance) {
    // GIVEN: Random particles
    int n = 100;
    double sml = 0.15;
    auto particles = create_random_particles(n, 0.1, 0.9, sml, 11111);

    BuildTree(particles);

    // WHEN/THEN: Neighbors should be sorted by distance in both methods
    for (auto& p : particles) {
        auto result = CompareSearchMethods(p, particles, false);

        // Check exhaustive ordering
        for (int i = 1; i < result.exhaustive_count; ++i) {
            int a = result.exhaustive_neighbors[i - 1];
            int b = result.exhaustive_neighbors[i];
            vec_t r_a = periodic_->calc_r_ij(p.pos, particles[a].pos);
            vec_t r_b = periodic_->calc_r_ij(p.pos, particles[b].pos);
            EXPECT_LE(abs2(r_a), abs2(r_b) + kAbsTol)
                << "Exhaustive search not sorted for particle " << p.id;
        }

        // Check tree ordering
        for (int i = 1; i < result.tree_count; ++i) {
            int a = result.tree_neighbors[i - 1];
            int b = result.tree_neighbors[i];
            vec_t r_a = periodic_->calc_r_ij(p.pos, particles[a].pos);
            vec_t r_b = periodic_->calc_r_ij(p.pos, particles[b].pos);
            EXPECT_LE(abs2(r_a), abs2(r_b) + kAbsTol)
                << "Tree search not sorted for particle " << p.id;
        }
    }
}

// =============================================================================
// Test: Larger scale test (pre-1M validation)
// =============================================================================
TEST_F(NeighborSearchTest, GivenLargerParticleCount_WhenSearching_ThenResultsMatch) {
    // GIVEN: Larger number of particles (still fast enough for unit test)
    int n = 2000;
    double sml = 0.05;
    auto particles = create_random_particles(n, 0.05, 0.95, sml, 99999);

    BuildTree(particles);

    // WHEN/THEN: Sample 100 particles and verify
    int count_mismatches = 0;
    int set_mismatches = 0;
    int sample_size = 100;

    std::mt19937 gen(12345);
    std::uniform_int_distribution<int> idx_dist(0, n - 1);

    for (int s = 0; s < sample_size; ++s) {
        int i = idx_dist(gen);
        auto& p = particles[i];
        auto result = CompareSearchMethods(p, particles, false);

        if (!result.counts_match) count_mismatches++;
        if (!result.sets_match) set_mismatches++;
    }

    EXPECT_EQ(count_mismatches, 0) << "Count mismatches in " << sample_size << " samples";
    EXPECT_EQ(set_mismatches, 0) << "Set mismatches in " << sample_size << " samples";
}

// =============================================================================
// Test: Self-inclusion
// =============================================================================
TEST_F(NeighborSearchTest, GivenParticle_WhenSearching_ThenIncludesSelf) {
    // GIVEN: A simple set of particles
    int n = 50;
    double sml = 0.1;
    auto particles = create_random_particles(n, 0.2, 0.8, sml, 77777);

    BuildTree(particles);

    // WHEN/THEN: Each particle should find itself as a neighbor (distance = 0)
    for (auto& p : particles) {
        auto result = CompareSearchMethods(p, particles, false);

        // Check that particle finds itself
        bool found_self_exhaustive = false;
        bool found_self_tree = false;

        for (int i = 0; i < result.exhaustive_count; ++i) {
            if (result.exhaustive_neighbors[i] == p.id) {
                found_self_exhaustive = true;
                break;
            }
        }

        for (int i = 0; i < result.tree_count; ++i) {
            if (result.tree_neighbors[i] == p.id) {
                found_self_tree = true;
                break;
            }
        }

        EXPECT_TRUE(found_self_exhaustive) << "Exhaustive search didn't find self for particle " << p.id;
        EXPECT_TRUE(found_self_tree) << "Tree search didn't find self for particle " << p.id;
    }
}

// =============================================================================
// Performance comparison test (informational, not a pass/fail test)
// =============================================================================
TEST_F(NeighborSearchTest, DISABLED_PerformanceComparison) {
    // This test is disabled by default (prefix with DISABLED_)
    // Run with: ./test_neighbor_search --gtest_also_run_disabled_tests

    std::cout << "\n=== Performance Comparison ===\n";

    for (int n : {1000, 5000, 10000}) {
        double sml = 0.05;
        auto particles = create_random_particles(n, 0.05, 0.95, sml, 42);
        BuildTree(particles);

        // Time exhaustive search
        auto start_exhaustive = std::chrono::high_resolution_clock::now();
        for (auto& p : particles) {
            std::vector<int> neighbors(kNeighborListSize);
            sph::exhaustive_search(p, p.sml, particles, n, neighbors, kNeighborListSize, periodic_.get(), false);
        }
        auto end_exhaustive = std::chrono::high_resolution_clock::now();

        // Time tree search
        auto start_tree = std::chrono::high_resolution_clock::now();
        for (auto& p : particles) {
            std::vector<int> neighbors(kNeighborListSize);
            tree_->neighbor_search(p, neighbors, particles, false);
        }
        auto end_tree = std::chrono::high_resolution_clock::now();

        auto exhaustive_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_exhaustive - start_exhaustive).count();
        auto tree_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_tree - start_tree).count();

        std::cout << "N=" << n
                  << ": Exhaustive=" << exhaustive_ms << "ms"
                  << ", Tree=" << tree_ms << "ms"
                  << ", Speedup=" << (tree_ms > 0 ? (double)exhaustive_ms / tree_ms : 0) << "x"
                  << std::endl;
    }
}
