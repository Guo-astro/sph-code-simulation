// Unit tests for Morton code particle reordering
//
// These tests verify that the Morton code implementation correctly:
// 1. Encodes positions to Morton codes (Z-order curve)
// 2. Sorts particles by Morton code
// 3. Preserves particle data during reordering

#include <gtest/gtest.h>
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>

#include "morton.hpp"
#include "particle.hpp"
#include "vector_type.hpp"  // for abs2()

using namespace sph;

class MortonTest : public ::testing::Test {
protected:
    // Domain bounds for testing
    real domain_min[DIM];
    real domain_size[DIM];

    void SetUp() override {
        for (int d = 0; d < DIM; ++d) {
            domain_min[d] = 0.0;
            domain_size[d] = 1.0;
        }
    }

    // Helper to create a particle at a given position
    SPHParticle create_particle(int id, real x, real y = 0.0, real z = 0.0) {
        SPHParticle p;
        p.id = id;
        p.pos[0] = x;
#if DIM >= 2
        p.pos[1] = y;
#endif
#if DIM == 3
        p.pos[2] = z;
#endif
        p.mass = 1.0;
        p.sml = 0.1;
        return p;
    }
};

// Test that Morton codes are computed correctly for corner positions
TEST_F(MortonTest, CornerPositions) {
    // Origin should have Morton code 0
    real origin[DIM] = {0.0};
    uint64_t code_origin = morton::encode(origin, domain_min, domain_size);
    EXPECT_EQ(code_origin, 0ULL);

    // Maximum position should have maximum Morton code
    real max_pos[DIM];
    for (int d = 0; d < DIM; ++d) {
        max_pos[d] = 1.0 - 1e-10;  // Just under 1.0
    }
    uint64_t code_max = morton::encode(max_pos, domain_min, domain_size);
    EXPECT_GT(code_max, 0ULL);
}

// Test that positions along X axis are ordered correctly
TEST_F(MortonTest, XAxisOrdering) {
    std::vector<uint64_t> codes;
    for (int i = 0; i < 10; ++i) {
        real pos[DIM] = {0.0};
        pos[0] = i * 0.1;
        codes.push_back(morton::encode(pos, domain_min, domain_size));
    }

    // Codes should be monotonically increasing along X axis
    for (size_t i = 1; i < codes.size(); ++i) {
        EXPECT_GE(codes[i], codes[i-1])
            << "Morton codes not monotonic at index " << i;
    }
}

// Test that nearby positions have similar Morton codes
TEST_F(MortonTest, SpatialLocality) {
    // Two nearby points should have more similar codes than distant points
    real pos_a[DIM], pos_b[DIM], pos_c[DIM];

    for (int d = 0; d < DIM; ++d) {
        pos_a[d] = 0.1;
        pos_b[d] = 0.11;  // Very close to A
        pos_c[d] = 0.9;   // Far from A
    }

    uint64_t code_a = morton::encode(pos_a, domain_min, domain_size);
    uint64_t code_b = morton::encode(pos_b, domain_min, domain_size);
    uint64_t code_c = morton::encode(pos_c, domain_min, domain_size);

    // The difference between nearby points should be smaller
    uint64_t diff_ab = (code_a > code_b) ? (code_a - code_b) : (code_b - code_a);
    uint64_t diff_ac = (code_a > code_c) ? (code_a - code_c) : (code_c - code_a);

    EXPECT_LT(diff_ab, diff_ac)
        << "Nearby points should have more similar Morton codes";
}

// Test that sorting particles by Morton code preserves particle data
TEST_F(MortonTest, SortPreservesData) {
    std::vector<SPHParticle> particles;

    // Create particles with distinct positions and data
    particles.push_back(create_particle(0, 0.9, 0.9, 0.9));
    particles.push_back(create_particle(1, 0.1, 0.1, 0.1));
    particles.push_back(create_particle(2, 0.5, 0.5, 0.5));
    particles.push_back(create_particle(3, 0.3, 0.7, 0.2));

    // Store original positions
    std::vector<vec_t> original_positions;
    for (const auto& p : particles) {
        original_positions.push_back(p.pos);
    }

    // Sort by Morton code
    auto permutation = morton::compute_sort_permutation(particles, domain_min, domain_size);
    morton::apply_permutation(particles, permutation);

    // Verify that all original positions are still present
    for (const auto& orig_pos : original_positions) {
        bool found = false;
        for (const auto& p : particles) {
            // Use abs2 (squared magnitude) for comparison
            if (abs2(p.pos - orig_pos) < 1e-20) {
                found = true;
                break;
            }
        }
        EXPECT_TRUE(found) << "Position lost during Morton sort";
    }
}

// Test that inverse permutation correctly maps indices
TEST_F(MortonTest, InversePermutation) {
    std::vector<int> permutation = {2, 0, 3, 1};
    auto inverse = morton::compute_inverse_permutation(permutation);

    // Verify: permutation[inverse[i]] == i for all i
    for (size_t i = 0; i < permutation.size(); ++i) {
        EXPECT_EQ(permutation[inverse[i]], static_cast<int>(i))
            << "Inverse permutation incorrect at index " << i;
    }

    // Verify: inverse[permutation[i]] == i for all i
    for (size_t i = 0; i < permutation.size(); ++i) {
        EXPECT_EQ(inverse[permutation[i]], static_cast<int>(i))
            << "Inverse permutation incorrect at index " << i;
    }
}

// Test sorting a large number of random particles
TEST_F(MortonTest, RandomParticlesSorted) {
    const int N = 1000;
    std::vector<SPHParticle> particles;

    std::mt19937 rng(42);  // Fixed seed for reproducibility
    std::uniform_real_distribution<real> dist(0.0, 1.0);

    for (int i = 0; i < N; ++i) {
        SPHParticle p;
        p.id = i;
        for (int d = 0; d < DIM; ++d) {
            p.pos[d] = dist(rng);
        }
        p.mass = 1.0;
        p.sml = 0.1;
        particles.push_back(p);
    }

    // Sort by Morton code
    morton::sort_particles_by_morton(particles, domain_min, domain_size);

    // Verify particles are sorted by Morton code
    std::vector<uint64_t> codes;
    for (const auto& p : particles) {
        codes.push_back(morton::encode(p.pos.get_array(), domain_min, domain_size));
    }

    for (size_t i = 1; i < codes.size(); ++i) {
        EXPECT_LE(codes[i-1], codes[i])
            << "Particles not sorted by Morton code at index " << i;
    }
}

// Test that particle IDs are updated after reordering
TEST_F(MortonTest, ParticleIDsPreserved) {
    // Create particles with specific IDs
    std::vector<SPHParticle> particles;
    particles.push_back(create_particle(100, 0.9, 0.9, 0.9));  // Original ID = 100
    particles.push_back(create_particle(200, 0.1, 0.1, 0.1));  // Original ID = 200
    particles.push_back(create_particle(300, 0.5, 0.5, 0.5));  // Original ID = 300

    // Store original IDs
    std::set<int> original_ids;
    for (const auto& p : particles) {
        original_ids.insert(p.id);
    }

    morton::sort_particles_by_morton(particles, domain_min, domain_size);

    // After reordering, original IDs should be PRESERVED (not overwritten)
    std::set<int> final_ids;
    for (const auto& p : particles) {
        final_ids.insert(p.id);
    }

    EXPECT_EQ(original_ids, final_ids) << "Particle IDs were corrupted during Morton reordering";

    // Verify specific IDs still exist
    EXPECT_TRUE(final_ids.count(100)) << "ID 100 was lost";
    EXPECT_TRUE(final_ids.count(200)) << "ID 200 was lost";
    EXPECT_TRUE(final_ids.count(300)) << "ID 300 was lost";
}

// Test edge case: empty particle list
TEST_F(MortonTest, EmptyParticleList) {
    std::vector<SPHParticle> particles;

    // Should not crash
    auto permutation = morton::compute_sort_permutation(particles, domain_min, domain_size);
    EXPECT_TRUE(permutation.empty());

    morton::sort_particles_by_morton(particles, domain_min, domain_size);
    EXPECT_TRUE(particles.empty());
}

// Test edge case: single particle
TEST_F(MortonTest, SingleParticle) {
    std::vector<SPHParticle> particles;
    particles.push_back(create_particle(0, 0.5, 0.5, 0.5));

    morton::sort_particles_by_morton(particles, domain_min, domain_size);

    EXPECT_EQ(particles.size(), 1u);
    EXPECT_EQ(particles[0].id, 0);
}

// Test with non-unit domain
TEST_F(MortonTest, NonUnitDomain) {
    // Set domain to [-10, 10]
    for (int d = 0; d < DIM; ++d) {
        domain_min[d] = -10.0;
        domain_size[d] = 20.0;
    }

    std::vector<SPHParticle> particles;

    SPHParticle p1;
    p1.id = 0;
    for (int d = 0; d < DIM; ++d) p1.pos[d] = -5.0;
    particles.push_back(p1);

    SPHParticle p2;
    p2.id = 1;
    for (int d = 0; d < DIM; ++d) p2.pos[d] = 5.0;
    particles.push_back(p2);

    // p1 at -5 should have smaller Morton code than p2 at +5
    uint64_t code1 = morton::encode(p1.pos.get_array(), domain_min, domain_size);
    uint64_t code2 = morton::encode(p2.pos.get_array(), domain_min, domain_size);

    EXPECT_LT(code1, code2);
}

#if DIM == 2
// Test 2D-specific Z-order pattern
TEST_F(MortonTest, ZOrderPattern2D) {
    // In 2D, the Z-order pattern visits cells in this order:
    // (0,0) -> (1,0) -> (0,1) -> (1,1) -> (2,0) -> ...
    // Lower-left should come before upper-right

    real pos_ll[2] = {0.1, 0.1};  // Lower-left
    real pos_ur[2] = {0.9, 0.9};  // Upper-right

    uint64_t code_ll = morton::encode(pos_ll, domain_min, domain_size);
    uint64_t code_ur = morton::encode(pos_ur, domain_min, domain_size);

    EXPECT_LT(code_ll, code_ur);
}
#endif

#if DIM == 3
// Test 3D-specific bit interleaving
TEST_F(MortonTest, BitInterleaving3D) {
    // Test that 3D Morton code correctly interleaves bits
    // Position (1,0,0) should give code with only x-bits set
    // Position (0,1,0) should give code with only y-bits set
    // Position (0,0,1) should give code with only z-bits set

    // These positions are at normalized coordinate 0.5
    real pos_x[3] = {0.5, 0.0, 0.0};
    real pos_y[3] = {0.0, 0.5, 0.0};
    real pos_z[3] = {0.0, 0.0, 0.5};

    uint64_t code_x = morton::encode(pos_x, domain_min, domain_size);
    uint64_t code_y = morton::encode(pos_y, domain_min, domain_size);
    uint64_t code_z = morton::encode(pos_z, domain_min, domain_size);

    // All three codes should be different
    EXPECT_NE(code_x, code_y);
    EXPECT_NE(code_y, code_z);
    EXPECT_NE(code_x, code_z);

    // Combined position should have higher code than any single dimension
    real pos_all[3] = {0.5, 0.5, 0.5};
    uint64_t code_all = morton::encode(pos_all, domain_min, domain_size);

    EXPECT_GT(code_all, code_x);
    EXPECT_GT(code_all, code_y);
    EXPECT_GT(code_all, code_z);
}
#endif
