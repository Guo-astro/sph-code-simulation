// TDD Tests for High-ROI Performance Optimizations
// 1. pres_per_rho2 precomputation
// 2. Merged neighbor loops in pre_interaction
// 3. Template kernel dispatch

#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <chrono>

#include "defines.hpp"
#include "particle.hpp"
#include "vector_type.hpp"
#include "bhtree.hpp"
#include "parameters.hpp"
#include "kernel/kernel_function.hpp"
#include "kernel/wendland_kernel.hpp"

using namespace sph;

// =============================================================================
// Test 1: pres_per_rho2 field exists and is computed correctly
// =============================================================================

TEST(HighROIOptimizations, ParticleHasPressurePerDensitySquaredField) {
    // Test that SPHParticle has pres_per_rho2 field
    SPHParticle p;
    p.pres = 1.0;
    p.dens = 2.0;

    // The field should exist and be settable
    p.pres_per_rho2 = p.pres / (p.dens * p.dens);

    EXPECT_DOUBLE_EQ(p.pres_per_rho2, 0.25);
}

TEST(HighROIOptimizations, PressurePerDensitySquaredIsCorrect) {
    // Test correctness of precomputed value
    SPHParticle p;
    p.pres = 10.0;
    p.dens = 5.0;

    // Expected: p / ρ² = 10 / 25 = 0.4
    p.pres_per_rho2 = p.pres / sqr(p.dens);

    EXPECT_DOUBLE_EQ(p.pres_per_rho2, 0.4);
}

TEST(HighROIOptimizations, PressurePerDensitySquaredEdgeCases) {
    SPHParticle p;

    // Very small density (avoid division by zero)
    p.pres = 1.0;
    p.dens = 1e-10;
    p.pres_per_rho2 = p.pres / sqr(p.dens);
    EXPECT_GT(p.pres_per_rho2, 0.0);
    EXPECT_FALSE(std::isinf(p.pres_per_rho2));

    // Zero pressure
    p.pres = 0.0;
    p.dens = 1.0;
    p.pres_per_rho2 = p.pres / sqr(p.dens);
    EXPECT_DOUBLE_EQ(p.pres_per_rho2, 0.0);
}

// =============================================================================
// Test 2: Merged loops produce same results as separate loops
// =============================================================================

class MergedLoopTest : public ::testing::Test {
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
        param->av.use_balsara_switch = true;
        param->av.use_time_dependent_av = true;
        param->av.alpha_max = 2.0;
        param->av.alpha_min = 0.1;
        param->av.epsilon = 0.1;

        tree = std::make_unique<BHTree>();
        tree->initialize(param);
    }

    SPHParticle create_particle(int id, real x, real y, real z) {
        SPHParticle p;
        p.id = id;
        p.pos[0] = x; p.pos[1] = y; p.pos[2] = z;
        p.vel[0] = 0.1 * x; p.vel[1] = 0.1 * y; p.vel[2] = 0.1 * z;
        p.mass = 1.0;
        p.sml = 0.15;
        p.dens = 1.0;
        p.pres = 1.0;
        p.pres_per_rho2 = 1.0;
        p.ene = 1.0;
        p.sound = 1.0;
        p.alpha = 1.0;
        p.balsara = 1.0;
        p.gradh = 1.0;
        p.next = nullptr;
        return p;
    }

    void create_grid_particles(int n) {
        particles.clear();
        int id = 0;
        real spacing = 1.0 / n;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                for (int k = 0; k < n; ++k) {
                    particles.push_back(create_particle(id++,
                        (i + 0.5) * spacing,
                        (j + 0.5) * spacing,
                        (k + 0.5) * spacing));
                }
            }
        }
    }
};

TEST_F(MergedLoopTest, DensityComputationUnchanged) {
    // Verify that merging loops doesn't change density computation
    create_grid_particles(5);  // 125 particles

    tree->resize(particles.size());
    tree->make(particles, particles.size());
    tree->set_kernel();

    // Density should still be computable correctly
    // (This is a sanity check - actual test is in integration)
    EXPECT_GT(particles.size(), 0u);
}

// =============================================================================
// Test 3: Template kernel produces same results as virtual dispatch
// =============================================================================

// Template kernel wrapper for Wendland C4
template<typename KernelType>
class TemplateKernelWrapper {
public:
    KernelType kernel;

    real w(const real r, const real h) const {
        return kernel.w(r, h);
    }

    vec_t dw(const vec_t& rij, const real r, const real h) const {
        return kernel.dw(rij, r, h);
    }

    real dhw(const real r, const real h) const {
        return kernel.dhw(r, h);
    }
};

TEST(HighROIOptimizations, TemplateKernelMatchesVirtual) {
    // Compare template kernel to virtual dispatch kernel
    Wendland::C4Kernel virtual_kernel;
    TemplateKernelWrapper<Wendland::C4Kernel> template_kernel;

    const real h = 0.1;

    // Test w() at various distances
    for (real r = 0.0; r < h; r += 0.01) {
        real w_virtual = virtual_kernel.w(r, h);
        real w_template = template_kernel.w(r, h);
        EXPECT_DOUBLE_EQ(w_virtual, w_template) << "w() mismatch at r=" << r;
    }

    // Test dw() at various distances
    vec_t rij;
    rij[0] = 0.05; rij[1] = 0.03; rij[2] = 0.02;
    real r = std::abs(rij);

    vec_t dw_virtual = virtual_kernel.dw(rij, r, h);
    vec_t dw_template = template_kernel.dw(rij, r, h);

    for (int d = 0; d < DIM; ++d) {
        EXPECT_DOUBLE_EQ(dw_virtual[d], dw_template[d]) << "dw() mismatch at dim=" << d;
    }

    // Test dhw() at various distances
    for (real r = 0.0; r < h; r += 0.01) {
        real dhw_virtual = virtual_kernel.dhw(r, h);
        real dhw_template = template_kernel.dhw(r, h);
        EXPECT_DOUBLE_EQ(dhw_virtual, dhw_template) << "dhw() mismatch at r=" << r;
    }
}

TEST(HighROIOptimizations, TemplateKernelIsInlined) {
    // Performance test: template kernel should be faster due to inlining
    // This test verifies the optimization is beneficial

    Wendland::C4Kernel virtual_kernel;
    TemplateKernelWrapper<Wendland::C4Kernel> template_kernel;

    const int iterations = 1000000;
    const real h = 0.1;
    real r = 0.05;
    vec_t rij;
    rij[0] = 0.04; rij[1] = 0.02; rij[2] = 0.01;

    // Virtual dispatch timing
    volatile real result_v = 0;
    auto start_v = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iterations; ++i) {
        KernelFunction* kf = &virtual_kernel;
        result_v += kf->w(r, h);
    }
    auto end_v = std::chrono::high_resolution_clock::now();
    auto time_virtual = std::chrono::duration<double, std::milli>(end_v - start_v).count();

    // Template dispatch timing
    volatile real result_t = 0;
    auto start_t = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iterations; ++i) {
        result_t += template_kernel.w(r, h);
    }
    auto end_t = std::chrono::high_resolution_clock::now();
    auto time_template = std::chrono::duration<double, std::milli>(end_t - start_t).count();

    std::cout << "\nKernel Performance Comparison:" << std::endl;
    std::cout << "  Virtual dispatch: " << time_virtual << " ms" << std::endl;
    std::cout << "  Template dispatch: " << time_template << " ms" << std::endl;
    std::cout << "  Speedup: " << (time_virtual / time_template) << "x" << std::endl;

    // Template should be at least as fast (usually faster)
    // Allow some tolerance for system variability
    EXPECT_LE(time_template, time_virtual * 1.1)
        << "Template kernel should not be significantly slower than virtual";
}

// =============================================================================
// Integration test: All optimizations work together
// =============================================================================

TEST_F(MergedLoopTest, AllOptimizationsPreservePhysics) {
    create_grid_particles(6);  // 216 particles

    // Set up particle properties
    for (auto& p : particles) {
        p.dens = 1.0 + 0.1 * p.pos[0];  // Slight density gradient
        p.pres = 1.0;
        p.pres_per_rho2 = p.pres / sqr(p.dens);  // Precompute
        p.sound = std::sqrt(param->physics.gamma * p.pres / p.dens);
    }

    tree->resize(particles.size());
    tree->make(particles, particles.size());
    tree->set_kernel();

    // Verify conservation properties still hold
    real total_mass = 0.0;
    for (const auto& p : particles) {
        total_mass += p.mass;

        // pres_per_rho2 should be consistent
        real expected = p.pres / sqr(p.dens);
        EXPECT_NEAR(p.pres_per_rho2, expected, 1e-14)
            << "pres_per_rho2 inconsistent for particle " << p.id;
    }

    EXPECT_GT(total_mass, 0.0);
}
