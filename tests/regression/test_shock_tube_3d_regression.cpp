// Regression tests for 3D Shock Tube simulations
// These tests ensure numerical results remain consistent across code changes.
//
// The tests compare RAW particle data against stored reference data.
// L2 error threshold: 1e-10 (to catch any numerical changes)
//
// Reference data generated from validated simulations:
// - Wendland kernel with N=120 neighbors
// - Cubic Spline kernel with N=50 neighbors

#include <gtest/gtest.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>
#include <filesystem>

// JSON parsing
#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace {

const double L2_ERROR_THRESHOLD = 1e-10;

// Get project root directory (handles running from build_test/ or project root)
std::string getProjectRoot() {
    std::filesystem::path cwd = std::filesystem::current_path();
    // If we're in build_test, go up one level
    if (cwd.filename() == "build_test") {
        return cwd.parent_path().string();
    }
    // Check if tests/regression exists in current directory
    if (std::filesystem::exists(cwd / "tests" / "regression")) {
        return cwd.string();
    }
    // Try parent directory
    if (std::filesystem::exists(cwd.parent_path() / "tests" / "regression")) {
        return cwd.parent_path().string();
    }
    return cwd.string();
}

const std::string PROJECT_ROOT = getProjectRoot();
const std::string REGRESSION_DATA_DIR = PROJECT_ROOT + "/tests/regression/reference_data";
const std::string RESULTS_BASE_DIR = PROJECT_ROOT + "/simulations/benchmarks/shock_tube/results";

struct ParticleData {
    int id;
    double x, y, z;
    double vx, vy, vz;
    double rho, pres, ene;
    double mass, sml;
};

struct SimulationData {
    std::vector<ParticleData> particles;
    double time;
    int num_particles;
};

// Load reference data from JSON (raw particle data)
SimulationData loadReferenceData(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open reference file: " + filename);
    }

    json j;
    file >> j;

    SimulationData data;
    data.time = j["time"];
    data.num_particles = j["num_particles"];

    auto& particles_json = j["particles"];
    data.particles.reserve(particles_json.size());

    for (const auto& p : particles_json) {
        ParticleData particle;
        particle.id = p["id"];
        particle.x = p["x"];
        particle.y = p["y"];
        particle.z = p["z"];
        particle.vx = p["vx"];
        particle.vy = p["vy"];
        particle.vz = p["vz"];
        particle.rho = p["rho"];
        particle.pres = p["pres"];
        particle.ene = p["ene"];
        particle.mass = p["mass"];
        particle.sml = p["sml"];
        data.particles.push_back(particle);
    }

    // Sort by ID for consistent comparison
    std::sort(data.particles.begin(), data.particles.end(),
              [](const ParticleData& a, const ParticleData& b) { return a.id < b.id; });

    return data;
}

// Load simulation CSV (raw particle data)
SimulationData loadSimulationData(const std::string& csv_filename) {
    std::ifstream file(csv_filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open CSV file: " + csv_filename);
    }

    SimulationData data;
    data.time = 0.0;

    std::string line;
    std::vector<std::string> header;
    std::map<std::string, int> col_idx;

    while (std::getline(file, line)) {
        if (line.empty()) continue;

        // Parse comment lines for time
        if (line[0] == '#') {
            if (line.find("Time (code):") != std::string::npos) {
                size_t pos = line.find(':');
                if (pos != std::string::npos) {
                    data.time = std::stod(line.substr(pos + 1));
                }
            }
            continue;
        }

        // Parse header
        if (header.empty()) {
            std::stringstream ss(line);
            std::string col;
            int idx = 0;
            while (std::getline(ss, col, ',')) {
                header.push_back(col);
                col_idx[col] = idx++;
            }
            continue;
        }

        // Parse data
        std::stringstream ss(line);
        std::string val;
        std::vector<std::string> values;
        while (std::getline(ss, val, ',')) {
            values.push_back(val);
        }

        if (values.size() >= header.size()) {
            ParticleData p;
            p.id = std::stoi(values[col_idx["id"]]);
            p.x = std::stod(values[col_idx["pos_x"]]);
            p.y = std::stod(values[col_idx["pos_y"]]);
            p.z = std::stod(values[col_idx["pos_z"]]);
            p.vx = std::stod(values[col_idx["vel_x"]]);
            p.vy = std::stod(values[col_idx["vel_y"]]);
            p.vz = std::stod(values[col_idx["vel_z"]]);
            p.rho = std::stod(values[col_idx["dens"]]);
            p.pres = std::stod(values[col_idx["pres"]]);
            p.ene = std::stod(values[col_idx["ene"]]);
            p.mass = std::stod(values[col_idx["mass"]]);
            p.sml = std::stod(values[col_idx["sml"]]);
            data.particles.push_back(p);
        }
    }

    data.num_particles = data.particles.size();

    // Sort by ID for consistent comparison
    std::sort(data.particles.begin(), data.particles.end(),
              [](const ParticleData& a, const ParticleData& b) { return a.id < b.id; });

    return data;
}

// Compute L2 error for a specific field across all particles
double computeL2Error(const std::vector<ParticleData>& ref,
                      const std::vector<ParticleData>& test,
                      double ParticleData::*field) {
    if (ref.size() != test.size()) {
        throw std::runtime_error("Particle count mismatch in L2 error computation");
    }

    double sum_sq = 0.0;
    double sum_ref_sq = 0.0;

    for (size_t i = 0; i < ref.size(); ++i) {
        double diff = test[i].*field - ref[i].*field;
        sum_sq += diff * diff;
        sum_ref_sq += ref[i].*field * ref[i].*field;
    }

    // Relative L2 error
    if (sum_ref_sq > 1e-30) {
        return std::sqrt(sum_sq / sum_ref_sq);
    }
    return std::sqrt(sum_sq);
}

} // anonymous namespace

class ShockTube3DRegressionTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Check if we're running from the correct directory
        if (!std::filesystem::exists(REGRESSION_DATA_DIR)) {
            GTEST_SKIP() << "Reference data directory not found. Run from project root.";
        }
    }
};

// Test: Wendland kernel with N=120 neighbors
TEST_F(ShockTube3DRegressionTest, WendlandN120_RegressionCheck) {
    const std::string ref_file = REGRESSION_DATA_DIR + "/shock_tube_3d_wendland_n120.json";
    const std::string sim_file = RESULTS_BASE_DIR + "/gsph_3d_wendland_n120/snapshot_0010.csv";

    // Skip if files don't exist
    if (!std::filesystem::exists(ref_file)) {
        GTEST_SKIP() << "Reference file not found: " << ref_file;
    }
    if (!std::filesystem::exists(sim_file)) {
        GTEST_SKIP() << "Simulation file not found: " << sim_file;
    }

    // Load reference and simulation data
    SimulationData ref = loadReferenceData(ref_file);
    SimulationData sim = loadSimulationData(sim_file);

    // Check particle count matches
    ASSERT_EQ(ref.num_particles, sim.num_particles)
        << "Particle count mismatch: ref=" << ref.num_particles << ", sim=" << sim.num_particles;

    // Check time matches
    EXPECT_NEAR(ref.time, sim.time, 1e-6) << "Simulation time mismatch";

    // Compute L2 errors for each field
    double l2_rho = computeL2Error(ref.particles, sim.particles, &ParticleData::rho);
    double l2_vx = computeL2Error(ref.particles, sim.particles, &ParticleData::vx);
    double l2_vy = computeL2Error(ref.particles, sim.particles, &ParticleData::vy);
    double l2_vz = computeL2Error(ref.particles, sim.particles, &ParticleData::vz);
    double l2_pres = computeL2Error(ref.particles, sim.particles, &ParticleData::pres);
    double l2_ene = computeL2Error(ref.particles, sim.particles, &ParticleData::ene);

    std::cout << "Wendland N=120 L2 Errors (raw particle data):" << std::endl;
    std::cout << "  Density:    " << l2_rho << std::endl;
    std::cout << "  Velocity X: " << l2_vx << std::endl;
    std::cout << "  Velocity Y: " << l2_vy << std::endl;
    std::cout << "  Velocity Z: " << l2_vz << std::endl;
    std::cout << "  Pressure:   " << l2_pres << std::endl;
    std::cout << "  Energy:     " << l2_ene << std::endl;

    // Assert L2 errors are below threshold
    EXPECT_LT(l2_rho, L2_ERROR_THRESHOLD)
        << "Density L2 error " << l2_rho << " exceeds threshold " << L2_ERROR_THRESHOLD;
    EXPECT_LT(l2_vx, L2_ERROR_THRESHOLD)
        << "Velocity X L2 error " << l2_vx << " exceeds threshold " << L2_ERROR_THRESHOLD;
    EXPECT_LT(l2_pres, L2_ERROR_THRESHOLD)
        << "Pressure L2 error " << l2_pres << " exceeds threshold " << L2_ERROR_THRESHOLD;
    EXPECT_LT(l2_ene, L2_ERROR_THRESHOLD)
        << "Energy L2 error " << l2_ene << " exceeds threshold " << L2_ERROR_THRESHOLD;
}

// Test: Cubic Spline kernel with N=50 neighbors
TEST_F(ShockTube3DRegressionTest, CubicSplineN50_RegressionCheck) {
    const std::string ref_file = REGRESSION_DATA_DIR + "/shock_tube_3d_cubic_spline_n50.json";
    const std::string sim_file = RESULTS_BASE_DIR + "/gsph_3d_hll/snapshot_0010.csv";

    // Skip if files don't exist
    if (!std::filesystem::exists(ref_file)) {
        GTEST_SKIP() << "Reference file not found: " << ref_file;
    }
    if (!std::filesystem::exists(sim_file)) {
        GTEST_SKIP() << "Simulation file not found: " << sim_file;
    }

    // Load reference and simulation data
    SimulationData ref = loadReferenceData(ref_file);
    SimulationData sim = loadSimulationData(sim_file);

    // Check particle count matches
    ASSERT_EQ(ref.num_particles, sim.num_particles)
        << "Particle count mismatch: ref=" << ref.num_particles << ", sim=" << sim.num_particles;

    // Check time matches
    EXPECT_NEAR(ref.time, sim.time, 1e-6) << "Simulation time mismatch";

    // Compute L2 errors for each field
    double l2_rho = computeL2Error(ref.particles, sim.particles, &ParticleData::rho);
    double l2_vx = computeL2Error(ref.particles, sim.particles, &ParticleData::vx);
    double l2_vy = computeL2Error(ref.particles, sim.particles, &ParticleData::vy);
    double l2_vz = computeL2Error(ref.particles, sim.particles, &ParticleData::vz);
    double l2_pres = computeL2Error(ref.particles, sim.particles, &ParticleData::pres);
    double l2_ene = computeL2Error(ref.particles, sim.particles, &ParticleData::ene);

    std::cout << "Cubic Spline N=50 L2 Errors (raw particle data):" << std::endl;
    std::cout << "  Density:    " << l2_rho << std::endl;
    std::cout << "  Velocity X: " << l2_vx << std::endl;
    std::cout << "  Velocity Y: " << l2_vy << std::endl;
    std::cout << "  Velocity Z: " << l2_vz << std::endl;
    std::cout << "  Pressure:   " << l2_pres << std::endl;
    std::cout << "  Energy:     " << l2_ene << std::endl;

    // Assert L2 errors are below threshold
    EXPECT_LT(l2_rho, L2_ERROR_THRESHOLD)
        << "Density L2 error " << l2_rho << " exceeds threshold " << L2_ERROR_THRESHOLD;
    EXPECT_LT(l2_vx, L2_ERROR_THRESHOLD)
        << "Velocity X L2 error " << l2_vx << " exceeds threshold " << L2_ERROR_THRESHOLD;
    EXPECT_LT(l2_pres, L2_ERROR_THRESHOLD)
        << "Pressure L2 error " << l2_pres << " exceeds threshold " << L2_ERROR_THRESHOLD;
    EXPECT_LT(l2_ene, L2_ERROR_THRESHOLD)
        << "Energy L2 error " << l2_ene << " exceeds threshold " << L2_ERROR_THRESHOLD;
}

// Test: Verify reference data files exist and are valid
TEST_F(ShockTube3DRegressionTest, ReferenceDataIntegrity) {
    std::vector<std::string> ref_files = {
        REGRESSION_DATA_DIR + "/shock_tube_3d_wendland_n120.json",
        REGRESSION_DATA_DIR + "/shock_tube_3d_cubic_spline_n50.json"
    };

    for (const auto& ref_file : ref_files) {
        ASSERT_TRUE(std::filesystem::exists(ref_file))
            << "Reference file missing: " << ref_file;

        // Try to load and validate
        ASSERT_NO_THROW({
            SimulationData data = loadReferenceData(ref_file);
            EXPECT_GT(data.num_particles, 0) << "No particles in " << ref_file;
            EXPECT_EQ(data.particles.size(), static_cast<size_t>(data.num_particles))
                << "Particle count mismatch in " << ref_file;
        }) << "Failed to load reference file: " << ref_file;
    }
}
