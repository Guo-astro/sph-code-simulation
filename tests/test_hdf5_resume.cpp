// Unit tests for HDF5 binary dump and resume functionality
//
// TDD Tests to verify:
// 1. HDF5 files can be written with particle data
// 2. Metadata is written correctly (time, step, physics params)
// 3. Particles can be read back for resume
// 4. Metadata can be read back for resume
// 5. Round-trip consistency (write -> read -> compare)

#include <gtest/gtest.h>
#include <vector>
#include <memory>
#include <cmath>
#include <random>
#include <cstdlib>
#include <sys/stat.h>

#include "writers/hdf5_writer.hpp"
#include "output_metadata.hpp"
#include "particle.hpp"
#include "parameters.hpp"

using namespace sph;

// Helper function to check if file exists
static bool file_exists(const std::string& path) {
    struct stat buffer;
    return (stat(path.c_str(), &buffer) == 0);
}

// Helper function to get file size
static size_t get_file_size(const std::string& path) {
    struct stat buffer;
    if (stat(path.c_str(), &buffer) == 0) {
        return buffer.st_size;
    }
    return 0;
}

// Helper function to create directory
static void create_directory(const std::string& path) {
    mkdir(path.c_str(), 0755);
}

// Helper function to remove file
static void remove_file(const std::string& path) {
    std::remove(path.c_str());
}

class HDF5ResumeTest : public ::testing::Test {
protected:
    std::string test_dir;
    std::string test_file;

    void SetUp() override {
        // Create temporary test directory
        test_dir = "/tmp/sph_test_hdf5_" + std::to_string(std::time(nullptr));
        create_directory(test_dir);
        test_file = test_dir + "/test_snapshot.h5";
    }

    void TearDown() override {
        // Clean up test files
        remove_file(test_file);
        remove_file(test_dir + "/no_compression.h5");
        remove_file(test_dir + "/max_compression.h5");
        remove_file(test_dir + "/snapshot_0021.h5");
        rmdir(test_dir.c_str());
    }

    // Helper to create a test particle with known values
    SPHParticle* create_test_particle(int id, real x, real y, real z,
                                       real vx, real vy, real vz,
                                       real mass, real sml, real dens, real ene) {
        SPHParticle* p = new SPHParticle();
        p->id = id;
        p->pos[0] = x;
#if DIM >= 2
        p->pos[1] = y;
#endif
#if DIM == 3
        p->pos[2] = z;
#endif
        p->vel[0] = vx;
#if DIM >= 2
        p->vel[1] = vy;
#endif
#if DIM == 3
        p->vel[2] = vz;
#endif
        p->mass = mass;
        p->sml = sml;
        p->dens = dens;
        p->ene = ene;
        p->pres = (1.6666666667 - 1.0) * dens * ene;
        return p;
    }

    // Helper to create test metadata
    OutputMetadata create_test_metadata(real time, int step) {
        OutputMetadata meta;
        meta.time_code = time;
        meta.step = step;
        meta.particle_count = 100;
        meta.sph_type = SPHType::GSPH;
        meta.kernel_type = KernelType::CUBIC_SPLINE;
        meta.gamma = 1.6666666667;
        meta.use_gravity = true;
        meta.gravitational_constant = 0.00430091;
        return meta;
    }
};

// Test 1: HDF5Writer can open file and write metadata
TEST_F(HDF5ResumeTest, WriteMetadata) {
    HDF5Writer writer(6);  // compression level 6
    OutputMetadata meta = create_test_metadata(1.5, 10);

    ASSERT_TRUE(writer.open(test_file, meta))
        << "Failed to open HDF5 file for writing";
    EXPECT_TRUE(writer.is_open());

    writer.close();
    EXPECT_FALSE(writer.is_open());

    // Verify file was created
    EXPECT_TRUE(file_exists(test_file))
        << "HDF5 file was not created";
}

// Test 2: Write and read particles with round-trip consistency
TEST_F(HDF5ResumeTest, WriteReadParticlesRoundTrip) {
    // Create test particles
    std::vector<SPHParticle*> particles_in;
    const int N = 100;

    std::mt19937 rng(42);
    std::uniform_real_distribution<real> dist(-10.0, 10.0);

    for (int i = 0; i < N; ++i) {
        particles_in.push_back(create_test_particle(
            i,                          // id
            dist(rng), dist(rng), dist(rng),  // position
            dist(rng) * 0.1, dist(rng) * 0.1, dist(rng) * 0.1,  // velocity
            1.0,                        // mass
            0.1 + i * 0.001,           // smoothing length
            100.0 + i,                 // density
            1.0 + i * 0.01             // internal energy
        ));
    }

    // Write particles
    OutputMetadata meta = create_test_metadata(2.5, 25);
    meta.particle_count = N;

    HDF5Writer writer(6);
    ASSERT_TRUE(writer.open(test_file, meta));
    ASSERT_TRUE(writer.write_particles(particles_in));
    writer.close();

    // Read particles back
    std::vector<SPHParticle*> particles_out;
    ASSERT_TRUE(HDF5Writer::read_particles(test_file, particles_out))
        << "Failed to read particles from HDF5 file";

    // Verify particle count
    ASSERT_EQ(particles_out.size(), particles_in.size())
        << "Particle count mismatch after round-trip";

    // Verify particle data (allowing for floating point tolerance)
    const real tol = 1e-10;
    for (int i = 0; i < N; ++i) {
        EXPECT_NEAR(particles_out[i]->pos[0], particles_in[i]->pos[0], tol)
            << "Position X mismatch for particle " << i;
#if DIM >= 2
        EXPECT_NEAR(particles_out[i]->pos[1], particles_in[i]->pos[1], tol)
            << "Position Y mismatch for particle " << i;
#endif
#if DIM == 3
        EXPECT_NEAR(particles_out[i]->pos[2], particles_in[i]->pos[2], tol)
            << "Position Z mismatch for particle " << i;
#endif
        EXPECT_NEAR(particles_out[i]->vel[0], particles_in[i]->vel[0], tol)
            << "Velocity X mismatch for particle " << i;
        EXPECT_NEAR(particles_out[i]->mass, particles_in[i]->mass, tol)
            << "Mass mismatch for particle " << i;
        EXPECT_NEAR(particles_out[i]->sml, particles_in[i]->sml, tol)
            << "Smoothing length mismatch for particle " << i;
        EXPECT_NEAR(particles_out[i]->dens, particles_in[i]->dens, tol)
            << "Density mismatch for particle " << i;
        EXPECT_NEAR(particles_out[i]->ene, particles_in[i]->ene, tol)
            << "Internal energy mismatch for particle " << i;
    }

    // Cleanup
    for (auto* p : particles_in) delete p;
    for (auto* p : particles_out) delete p;
}

// Test 3: Read metadata for resume
TEST_F(HDF5ResumeTest, ReadMetadataForResume) {
    // Write with specific metadata
    OutputMetadata meta_in = create_test_metadata(5.0, 50);
    meta_in.particle_count = 1000;
    meta_in.gamma = 1.4;
    meta_in.gravitational_constant = 0.0043;

    std::vector<SPHParticle*> particles;
    particles.push_back(create_test_particle(0, 1.0, 2.0, 3.0, 0.1, 0.2, 0.3, 1.0, 0.1, 100.0, 1.0));

    HDF5Writer writer(6);
    ASSERT_TRUE(writer.open(test_file, meta_in));
    ASSERT_TRUE(writer.write_particles(particles));
    writer.close();

    // Read metadata back
    OutputMetadata meta_out;
    ASSERT_TRUE(HDF5Writer::read_metadata(test_file, meta_out))
        << "Failed to read metadata from HDF5 file";

    // Verify metadata fields
    EXPECT_NEAR(meta_out.time_code, meta_in.time_code, 1e-10)
        << "Time mismatch in metadata";
    EXPECT_EQ(meta_out.step, meta_in.step)
        << "Step mismatch in metadata";
    EXPECT_NEAR(meta_out.gamma, meta_in.gamma, 1e-10)
        << "Gamma mismatch in metadata";
    EXPECT_NEAR(meta_out.gravitational_constant, meta_in.gravitational_constant, 1e-10)
        << "G mismatch in metadata";
    EXPECT_EQ(meta_out.use_gravity, meta_in.use_gravity)
        << "use_gravity mismatch in metadata";

    // Cleanup
    for (auto* p : particles) delete p;
}

// Test 4: Large particle count (10K particles for unit test speed)
TEST_F(HDF5ResumeTest, LargeParticleCount) {
    const int N = 10000;

    std::vector<SPHParticle*> particles_in;
    particles_in.reserve(N);

    std::mt19937 rng(123);
    std::uniform_real_distribution<real> pos_dist(-5.0, 5.0);
    std::uniform_real_distribution<real> vel_dist(-0.1, 0.1);

    for (int i = 0; i < N; ++i) {
        particles_in.push_back(create_test_particle(
            i,
            pos_dist(rng), pos_dist(rng), pos_dist(rng),
            vel_dist(rng), vel_dist(rng), vel_dist(rng),
            1.0 / N,  // total mass = 1
            0.1,
            100.0,
            1.0
        ));
    }

    OutputMetadata meta = create_test_metadata(0.0, 0);
    meta.particle_count = N;

    // Write
    HDF5Writer writer(6);
    ASSERT_TRUE(writer.open(test_file, meta));
    ASSERT_TRUE(writer.write_particles(particles_in));
    writer.close();

    // Verify file size is reasonable (binary should be efficient)
    auto file_size = get_file_size(test_file);
    EXPECT_GT(file_size, 0u) << "File is empty";
    EXPECT_LT(file_size, static_cast<size_t>(N) * 200u)
        << "File size unexpectedly large";

    // Read back
    std::vector<SPHParticle*> particles_out;
    ASSERT_TRUE(HDF5Writer::read_particles(test_file, particles_out));
    EXPECT_EQ(particles_out.size(), static_cast<size_t>(N));

    // Cleanup
    for (auto* p : particles_in) delete p;
    for (auto* p : particles_out) delete p;
}

// Test 5: File extension helper
TEST_F(HDF5ResumeTest, FileExtension) {
    EXPECT_EQ(HDF5Writer::get_extension(), ".h5")
        << "HDF5 file extension should be .h5";
}

// Test 6: Compression levels
TEST_F(HDF5ResumeTest, CompressionLevels) {
    const int N = 1000;
    std::vector<SPHParticle*> particles;

    for (int i = 0; i < N; ++i) {
        particles.push_back(create_test_particle(
            i, i * 0.01, i * 0.02, i * 0.03,
            0.0, 0.0, 0.0,
            1.0, 0.1, 100.0, 1.0
        ));
    }

    OutputMetadata meta = create_test_metadata(0.0, 0);
    meta.particle_count = N;

    // Write with no compression
    std::string file_no_comp = test_dir + "/no_compression.h5";
    HDF5Writer writer_no_comp(0);
    ASSERT_TRUE(writer_no_comp.open(file_no_comp, meta));
    ASSERT_TRUE(writer_no_comp.write_particles(particles));
    writer_no_comp.close();

    // Write with max compression
    std::string file_max_comp = test_dir + "/max_compression.h5";
    HDF5Writer writer_max_comp(9);
    ASSERT_TRUE(writer_max_comp.open(file_max_comp, meta));
    ASSERT_TRUE(writer_max_comp.write_particles(particles));
    writer_max_comp.close();

    // Compressed file should be smaller or equal (depends on data entropy)
    auto size_no_comp = get_file_size(file_no_comp);
    auto size_max_comp = get_file_size(file_max_comp);

    EXPECT_GT(size_no_comp, 0u);
    EXPECT_GT(size_max_comp, 0u);
    // With repetitive data, compression should help
    EXPECT_LE(size_max_comp, size_no_comp)
        << "Compression should reduce or maintain file size";

    // Cleanup
    for (auto* p : particles) delete p;
}

// Test 7: Resume workflow simulation
TEST_F(HDF5ResumeTest, ResumeWorkflow) {
    // Simulate Phase 1: Write relaxation checkpoint
    const int N = 100;
    std::vector<SPHParticle*> phase1_particles;

    for (int i = 0; i < N; ++i) {
        phase1_particles.push_back(create_test_particle(
            i,
            std::cos(i * 0.1), std::sin(i * 0.1), i * 0.01,  // position
            0.0, 0.0, 0.0,  // velocity (relaxed)
            1.0, 0.1, 100.0, 1.0
        ));
    }

    OutputMetadata phase1_meta = create_test_metadata(2.0, 21);
    phase1_meta.particle_count = N;
    phase1_meta.sph_type = SPHType::GSPH;
    phase1_meta.gamma = 1.6666666667;

    std::string checkpoint_file = test_dir + "/snapshot_0021.h5";
    HDF5Writer writer(6);
    ASSERT_TRUE(writer.open(checkpoint_file, phase1_meta));
    ASSERT_TRUE(writer.write_particles(phase1_particles));
    writer.close();

    // Simulate Phase 2: Resume from checkpoint
    std::vector<SPHParticle*> resumed_particles;
    OutputMetadata resumed_meta;

    ASSERT_TRUE(HDF5Writer::read_particles(checkpoint_file, resumed_particles))
        << "Failed to read particles for resume";
    ASSERT_TRUE(HDF5Writer::read_metadata(checkpoint_file, resumed_meta))
        << "Failed to read metadata for resume";

    // Verify resume state
    EXPECT_EQ(resumed_particles.size(), static_cast<size_t>(N));
    EXPECT_NEAR(resumed_meta.time_code, 2.0, 1e-10);
    EXPECT_EQ(resumed_meta.step, 21);

    // Cleanup
    for (auto* p : phase1_particles) delete p;
    for (auto* p : resumed_particles) delete p;
}

// Test 8: Error handling - invalid file path
TEST_F(HDF5ResumeTest, ErrorHandlingInvalidPath) {
    std::vector<SPHParticle*> particles;

    // Try to read from non-existent file
    EXPECT_FALSE(HDF5Writer::read_particles("/nonexistent/path/file.h5", particles))
        << "Should fail to read from non-existent file";
    EXPECT_TRUE(particles.empty());

    OutputMetadata meta;
    EXPECT_FALSE(HDF5Writer::read_metadata("/nonexistent/path/file.h5", meta))
        << "Should fail to read metadata from non-existent file";
}

// Test 9: Empty particle list
TEST_F(HDF5ResumeTest, EmptyParticleList) {
    std::vector<SPHParticle*> particles;
    OutputMetadata meta = create_test_metadata(0.0, 0);
    meta.particle_count = 0;

    HDF5Writer writer(6);
    ASSERT_TRUE(writer.open(test_file, meta));
    ASSERT_TRUE(writer.write_particles(particles));
    writer.close();

    // Read back
    std::vector<SPHParticle*> particles_out;
    ASSERT_TRUE(HDF5Writer::read_particles(test_file, particles_out));
    EXPECT_TRUE(particles_out.empty());
}
