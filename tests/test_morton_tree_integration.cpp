// Integration tests for Morton code with BHTree
//
// Tests that verify Morton reordering:
// 1. Is called during tree construction
// 2. Does not break neighbor search
// 3. Maintains simulation correctness
// 4. Improves cache locality (particles stay spatially coherent)

#include <gtest/gtest.h>
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>
#include <numeric>

#include "defines.hpp"
#include "particle.hpp"
#include "vector_type.hpp"
#include "bhtree.hpp"
#include "morton.hpp"
#include "parameters.hpp"

using namespace sph;

class MortonTreeIntegrationTest : public ::testing::Test {
protected:
    std::shared_ptr<SPHParameters> param;
    std::unique_ptr<BHTree> tree;
    std::vector<SPHParticle> particles;

    void SetUp() override {
        param = std::make_shared<SPHParameters>();
        param->physics.gamma = 1.6667;
        param->physics.neighbor_number = 50;
        param->tree.max_level = 15;
        param->tree.leaf_particle_num = 16;
        param->gravity.constant = 0.00430091;
        param->gravity.theta = 0.5;

        tree = std::make_unique<BHTree>();
        tree->initialize(param);
    }

    SPHParticle create_particle(int id, real x, real y, real z) {
        SPHParticle p;
        p.id = id;
        p.pos[0] = x;
        p.pos[1] = y;
        p.pos[2] = z;
        p.vel[0] = p.vel[1] = p.vel[2] = 0.0;
        p.mass = 1.0;
        p.sml = 0.1;
        p.dens = 1.0;
        p.pres = 1.0;
        p.ene = 1.0;
        p.next = nullptr;
        return p;
    }

    // Create particles in a random distribution
    void create_random_particles(int n, real box_size = 1.0, unsigned seed = 42) {
        std::mt19937 rng(seed);
        std::uniform_real_distribution<real> dist(0.0, box_size);

        particles.clear();
        particles.reserve(n);
        for (int i = 0; i < n; ++i) {
            particles.push_back(create_particle(i, dist(rng), dist(rng), dist(rng)));
        }
    }

    // Create particles on a regular grid
    void create_grid_particles(int n_per_side, real box_size = 1.0) {
        particles.clear();
        int id = 0;
        real spacing = box_size / n_per_side;
        for (int i = 0; i < n_per_side; ++i) {
            for (int j = 0; j < n_per_side; ++j) {
                for (int k = 0; k < n_per_side; ++k) {
                    particles.push_back(create_particle(
                        id++,
                        (i + 0.5) * spacing,
                        (j + 0.5) * spacing,
                        (k + 0.5) * spacing
                    ));
                }
            }
        }
    }

    // Compute total distance between consecutive particles in memory
    real compute_memory_locality_score() {
        real total_dist = 0.0;
        for (size_t i = 1; i < particles.size(); ++i) {
            real dx = particles[i].pos[0] - particles[i-1].pos[0];
            real dy = particles[i].pos[1] - particles[i-1].pos[1];
            real dz = particles[i].pos[2] - particles[i-1].pos[2];
            total_dist += std::sqrt(dx*dx + dy*dy + dz*dz);
        }
        return total_dist;
    }

    // Store original particle IDs for comparison
    std::vector<int> get_particle_ids() {
        std::vector<int> ids;
        ids.reserve(particles.size());
        for (const auto& p : particles) {
            ids.push_back(p.id);
        }
        return ids;
    }
};

// Test that particles are reordered after tree construction
TEST_F(MortonTreeIntegrationTest, ParticlesReorderedAfterTreeMake) {
    create_random_particles(1000);

    // Store original order
    auto original_ids = get_particle_ids();

    // Build tree (should trigger Morton reordering)
    tree->resize(particles.size());
    tree->make(particles, particles.size());

    // Get new order
    auto new_ids = get_particle_ids();

    // Particles should be reordered (order should differ)
    // Note: With Morton ordering, particles should be sorted by Morton code
    bool order_changed = (original_ids != new_ids);

    // All original IDs should still be present
    std::vector<int> sorted_original = original_ids;
    std::vector<int> sorted_new = new_ids;
    std::sort(sorted_original.begin(), sorted_original.end());
    std::sort(sorted_new.begin(), sorted_new.end());
    EXPECT_EQ(sorted_original, sorted_new) << "Particle IDs changed during reordering";

    // Morton ordering should change the order for random particles
    EXPECT_TRUE(order_changed) << "Morton reordering didn't change particle order";
}

// Test that neighbor search still works correctly after Morton reordering
TEST_F(MortonTreeIntegrationTest, NeighborSearchCorrectAfterReordering) {
    create_grid_particles(10);  // 1000 particles

    tree->resize(particles.size());
    tree->make(particles, particles.size());
    tree->set_kernel();

    // Test neighbor search for several particles
    std::vector<int> neighbor_list(200);

    for (int test_idx = 0; test_idx < 10; ++test_idx) {
        int idx = test_idx * 100;  // Sample particles at different positions
        if (idx >= static_cast<int>(particles.size())) break;

        int n_neighbors = tree->neighbor_search(particles[idx], neighbor_list, particles);

        // Should find some neighbors (grid has particles nearby)
        EXPECT_GT(n_neighbors, 0) << "No neighbors found for particle " << idx;

        // Verify neighbors are actually within kernel radius
        real h_i = particles[idx].sml;
        for (int i = 0; i < n_neighbors; ++i) {
            int j = neighbor_list[i];
            real dx = particles[idx].pos[0] - particles[j].pos[0];
            real dy = particles[idx].pos[1] - particles[j].pos[1];
            real dz = particles[idx].pos[2] - particles[j].pos[2];
            real r = std::sqrt(dx*dx + dy*dy + dz*dz);
            real h_max = std::max(h_i, particles[j].sml);

            // Neighbor should be within 2*h (kernel support)
            EXPECT_LE(r, 2.0 * h_max * 1.1)  // 10% tolerance
                << "Neighbor " << j << " too far from particle " << idx;
        }
    }
}

// Test that Morton ordering improves memory locality
TEST_F(MortonTreeIntegrationTest, MortonImprovesMemoryLocality) {
    create_random_particles(1000);

    // Compute locality before Morton ordering
    real locality_before = compute_memory_locality_score();

    // Apply Morton reordering
    real domain_min[DIM] = {0.0, 0.0, 0.0};
    real domain_size[DIM] = {1.0, 1.0, 1.0};
    morton::sort_particles_by_morton(particles, domain_min, domain_size);

    // Compute locality after Morton ordering
    real locality_after = compute_memory_locality_score();

    // Morton ordering should improve locality (lower total distance)
    EXPECT_LT(locality_after, locality_before)
        << "Morton ordering did not improve memory locality. "
        << "Before: " << locality_before << ", After: " << locality_after;
}

// Test that particle data is preserved during reordering
TEST_F(MortonTreeIntegrationTest, ParticleDataPreservedDuringReordering) {
    create_random_particles(100);

    // Store original data
    std::map<int, std::tuple<real, real, real>> original_positions;
    std::map<int, real> original_masses;

    for (const auto& p : particles) {
        original_positions[p.id] = std::make_tuple(p.pos[0], p.pos[1], p.pos[2]);
        original_masses[p.id] = p.mass;
    }

    // Apply Morton reordering via tree make
    tree->resize(particles.size());
    tree->make(particles, particles.size());

    // Verify all data is preserved
    for (const auto& p : particles) {
        auto orig_pos = original_positions[p.id];
        EXPECT_DOUBLE_EQ(p.pos[0], std::get<0>(orig_pos)) << "Position X changed for particle " << p.id;
        EXPECT_DOUBLE_EQ(p.pos[1], std::get<1>(orig_pos)) << "Position Y changed for particle " << p.id;
        EXPECT_DOUBLE_EQ(p.pos[2], std::get<2>(orig_pos)) << "Position Z changed for particle " << p.id;
        EXPECT_DOUBLE_EQ(p.mass, original_masses[p.id]) << "Mass changed for particle " << p.id;
    }
}

// Test with large particle count (performance-oriented)
TEST_F(MortonTreeIntegrationTest, LargeParticleCountWorks) {
    const int N = 10000;
    create_random_particles(N);

    tree->resize(particles.size());

    // Should not crash or hang
    ASSERT_NO_THROW({
        tree->make(particles, particles.size());
        tree->set_kernel();
    });

    // Verify tree was built
    std::vector<int> neighbor_list(200);
    int n = tree->neighbor_search(particles[0], neighbor_list, particles);
    EXPECT_GT(n, 0) << "Tree search failed after building with " << N << " particles";
}
