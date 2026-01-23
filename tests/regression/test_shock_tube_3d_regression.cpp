/**
 * @file test_shock_tube_3d_regression.cpp
 * @brief Regression tests for 3D Sod shock tube with HCP lattice
 *
 * This test ensures algorithm consistency by comparing simulation results
 * against known reference values with strict L2 norm tolerance.
 *
 * Test validates:
 * - Density field evolution
 * - Pressure field evolution
 * - Velocity field evolution
 * - Energy conservation
 * - Mass conservation
 *
 * Any change to the GSPH algorithm that modifies these values will cause
 * the test to fail, alerting developers to potential regressions.
 */

#include <gtest/gtest.h>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <filesystem>

namespace {

// =============================================================================
// Reference values computed from validated HCP low-res simulation at t=0.05
// Nx=25, 4425 particles, GSPH with HLL Riemann solver
// =============================================================================

// L2 norms of field values
constexpr double kRefDensL2 = 8.546236217634e-01;
constexpr double kRefPresL2 = 8.357204694215e-01;
constexpr double kRefVelxL2 = 3.280763623596e-01;
constexpr double kRefEneL2 = 2.397196504234e+00;

// Conservation quantities
constexpr double kRefTotalMass = 2.250000000000e-02;
constexpr double kRefTotalKE = 1.277899042146e-03;
constexpr double kRefTotalIE = 5.363283737885e-02;

// Expected particle count
constexpr int kExpectedParticles = 4425;

// =============================================================================
// Strict tolerance for regression detection
// Any algorithm change that modifies results beyond this threshold is flagged
// =============================================================================
constexpr double kL2RelTol = 1e-8;  // 10 ppb relative tolerance
constexpr double kConservationTol = 1e-10;  // Conservation should be very strict

// =============================================================================
// Helper functions
// =============================================================================

struct SimulationData {
    std::vector<double> pos_x, pos_y, pos_z;
    std::vector<double> vel_x, vel_y, vel_z;
    std::vector<double> dens, pres, ene, mass;
    int particle_count = 0;
    bool valid = false;
};

// Parse a CSV file, skipping comment lines starting with #
SimulationData parseCSV(const std::string& filepath) {
    SimulationData data;
    std::ifstream file(filepath);
    if (!file.is_open()) {
        return data;
    }

    std::string line;
    std::vector<std::string> headers;
    bool header_found = false;

    while (std::getline(file, line)) {
        // Skip comment lines
        if (line.empty() || line[0] == '#') continue;

        std::stringstream ss(line);
        std::string cell;

        if (!header_found) {
            // Parse header row
            while (std::getline(ss, cell, ',')) {
                headers.push_back(cell);
            }
            header_found = true;
            continue;
        }

        // Parse data row
        std::vector<double> row;
        while (std::getline(ss, cell, ',')) {
            try {
                row.push_back(std::stod(cell));
            } catch (...) {
                row.push_back(0.0);
            }
        }

        // Find column indices (assuming standard column names)
        auto findCol = [&](const std::string& name) -> int {
            for (size_t i = 0; i < headers.size(); ++i) {
                if (headers[i] == name) return static_cast<int>(i);
            }
            return -1;
        };

        int idx_pos_x = findCol("pos_x");
        int idx_pos_y = findCol("pos_y");
        int idx_pos_z = findCol("pos_z");
        int idx_vel_x = findCol("vel_x");
        int idx_vel_y = findCol("vel_y");
        int idx_vel_z = findCol("vel_z");
        int idx_dens = findCol("dens");
        int idx_pres = findCol("pres");
        int idx_ene = findCol("ene");
        int idx_mass = findCol("mass");

        if (idx_pos_x >= 0 && idx_pos_x < (int)row.size()) data.pos_x.push_back(row[idx_pos_x]);
        if (idx_pos_y >= 0 && idx_pos_y < (int)row.size()) data.pos_y.push_back(row[idx_pos_y]);
        if (idx_pos_z >= 0 && idx_pos_z < (int)row.size()) data.pos_z.push_back(row[idx_pos_z]);
        if (idx_vel_x >= 0 && idx_vel_x < (int)row.size()) data.vel_x.push_back(row[idx_vel_x]);
        if (idx_vel_y >= 0 && idx_vel_y < (int)row.size()) data.vel_y.push_back(row[idx_vel_y]);
        if (idx_vel_z >= 0 && idx_vel_z < (int)row.size()) data.vel_z.push_back(row[idx_vel_z]);
        if (idx_dens >= 0 && idx_dens < (int)row.size()) data.dens.push_back(row[idx_dens]);
        if (idx_pres >= 0 && idx_pres < (int)row.size()) data.pres.push_back(row[idx_pres]);
        if (idx_ene >= 0 && idx_ene < (int)row.size()) data.ene.push_back(row[idx_ene]);
        if (idx_mass >= 0 && idx_mass < (int)row.size()) data.mass.push_back(row[idx_mass]);

        data.particle_count++;
    }

    data.valid = (data.particle_count > 0);
    return data;
}

double computeL2Norm(const std::vector<double>& v) {
    if (v.empty()) return 0.0;
    double sum = 0.0;
    for (double x : v) sum += x * x;
    return std::sqrt(sum / v.size());
}

double computeSum(const std::vector<double>& v) {
    double sum = 0.0;
    for (double x : v) sum += x;
    return sum;
}

double computeWeightedSum(const std::vector<double>& weights, const std::vector<double>& values) {
    if (weights.size() != values.size()) return 0.0;
    double sum = 0.0;
    for (size_t i = 0; i < weights.size(); ++i) {
        sum += weights[i] * values[i];
    }
    return sum;
}

double computeKineticEnergy(const std::vector<double>& mass,
                            const std::vector<double>& vx,
                            const std::vector<double>& vy,
                            const std::vector<double>& vz) {
    if (mass.size() != vx.size()) return 0.0;
    double ke = 0.0;
    for (size_t i = 0; i < mass.size(); ++i) {
        ke += 0.5 * mass[i] * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
    }
    return ke;
}

} // anonymous namespace

// =============================================================================
// Test Fixture
// =============================================================================
class ShockTube3DRegressionTest : public ::testing::Test {
protected:
    static std::string test_output_dir;
    static std::string snapshot_path;
    static bool simulation_ran;
    static SimulationData sim_data;

    static void SetUpTestSuite() {
        // Create temporary output directory
        test_output_dir = "/tmp/sph_test_shock_tube_3d_regression";
        std::filesystem::create_directories(test_output_dir);

        // Create config file
        std::string config_path = test_output_dir + "/config.json";
        std::ofstream config(config_path);
        config << R"({
  "outputDirectory": ")" << test_output_dir << R"(",
  "startTime": 0.0,
  "endTime": 0.05,
  "outputTime": 0.05,
  "energyTime": 0.01,
  "checkpoint": {"enabled": false, "resumeFile": ""},
  "output": {"formats": [{"type": "csv", "precision": 16}], "enableEnergyFile": true},
  "sample": "shock_tube_3d",
  "Nx": 25,
  "neighborNumber": 50,
  "gamma": 1.4,
  "kernel": "cubic_spline",
  "periodic": true,
  "iterativeSmoothingLength": true,
  "rangeMax": [0.5, 0.1, 0.1],
  "rangeMin": [-0.5, -0.1, -0.1],
  "SPHType": "gsph",
  "use2ndOrderGSPH": false,
  "useBalsaraSwitch": false,
  "riemannSolver": "hllc",
  "useVolumeBased": false,
  "cflSound": 0.3,
  "cflForce": 0.25
})";
        config.close();

        // Run the simulation
        std::string cmd = "./build/sph " + config_path + " > /dev/null 2>&1";
        int ret = std::system(cmd.c_str());
        simulation_ran = (ret == 0);

        // Load results
        snapshot_path = test_output_dir + "/snapshot_0001.csv";
        if (simulation_ran) {
            sim_data = parseCSV(snapshot_path);
        }
    }

    static void TearDownTestSuite() {
        // Cleanup
        std::filesystem::remove_all(test_output_dir);
    }
};

std::string ShockTube3DRegressionTest::test_output_dir;
std::string ShockTube3DRegressionTest::snapshot_path;
bool ShockTube3DRegressionTest::simulation_ran = false;
SimulationData ShockTube3DRegressionTest::sim_data;

// =============================================================================
// Tests
// =============================================================================

TEST_F(ShockTube3DRegressionTest, SimulationRuns) {
    ASSERT_TRUE(simulation_ran) << "Simulation failed to run";
}

TEST_F(ShockTube3DRegressionTest, OutputFileExists) {
    ASSERT_TRUE(simulation_ran);
    ASSERT_TRUE(std::filesystem::exists(snapshot_path))
        << "Output file not found: " << snapshot_path;
}

TEST_F(ShockTube3DRegressionTest, ParticleCountMatches) {
    ASSERT_TRUE(sim_data.valid);
    EXPECT_EQ(sim_data.particle_count, kExpectedParticles)
        << "Particle count mismatch. Algorithm may have changed particle generation.";
}

TEST_F(ShockTube3DRegressionTest, DensityL2NormRegression) {
    ASSERT_TRUE(sim_data.valid);
    double dens_L2 = computeL2Norm(sim_data.dens);
    double rel_diff = std::abs(dens_L2 - kRefDensL2) / kRefDensL2;

    EXPECT_LT(rel_diff, kL2RelTol)
        << "DENSITY L2 REGRESSION DETECTED!\n"
        << "  Reference: " << kRefDensL2 << "\n"
        << "  Computed:  " << dens_L2 << "\n"
        << "  Rel diff:  " << rel_diff << " (tol: " << kL2RelTol << ")\n"
        << "This indicates a change in the GSPH density computation algorithm.";
}

TEST_F(ShockTube3DRegressionTest, PressureL2NormRegression) {
    ASSERT_TRUE(sim_data.valid);
    double pres_L2 = computeL2Norm(sim_data.pres);
    double rel_diff = std::abs(pres_L2 - kRefPresL2) / kRefPresL2;

    EXPECT_LT(rel_diff, kL2RelTol)
        << "PRESSURE L2 REGRESSION DETECTED!\n"
        << "  Reference: " << kRefPresL2 << "\n"
        << "  Computed:  " << pres_L2 << "\n"
        << "  Rel diff:  " << rel_diff << " (tol: " << kL2RelTol << ")\n"
        << "This indicates a change in the GSPH pressure/EOS algorithm.";
}

TEST_F(ShockTube3DRegressionTest, VelocityL2NormRegression) {
    ASSERT_TRUE(sim_data.valid);
    double velx_L2 = computeL2Norm(sim_data.vel_x);
    double rel_diff = std::abs(velx_L2 - kRefVelxL2) / kRefVelxL2;

    EXPECT_LT(rel_diff, kL2RelTol)
        << "VELOCITY L2 REGRESSION DETECTED!\n"
        << "  Reference: " << kRefVelxL2 << "\n"
        << "  Computed:  " << velx_L2 << "\n"
        << "  Rel diff:  " << rel_diff << " (tol: " << kL2RelTol << ")\n"
        << "This indicates a change in the GSPH momentum/force algorithm.";
}

TEST_F(ShockTube3DRegressionTest, EnergyL2NormRegression) {
    ASSERT_TRUE(sim_data.valid);
    double ene_L2 = computeL2Norm(sim_data.ene);
    double rel_diff = std::abs(ene_L2 - kRefEneL2) / kRefEneL2;

    EXPECT_LT(rel_diff, kL2RelTol)
        << "ENERGY L2 REGRESSION DETECTED!\n"
        << "  Reference: " << kRefEneL2 << "\n"
        << "  Computed:  " << ene_L2 << "\n"
        << "  Rel diff:  " << rel_diff << " (tol: " << kL2RelTol << ")\n"
        << "This indicates a change in the GSPH energy evolution algorithm.";
}

TEST_F(ShockTube3DRegressionTest, MassConservation) {
    ASSERT_TRUE(sim_data.valid);
    double total_mass = computeSum(sim_data.mass);
    double rel_diff = std::abs(total_mass - kRefTotalMass) / kRefTotalMass;

    EXPECT_LT(rel_diff, kConservationTol)
        << "MASS CONSERVATION VIOLATED!\n"
        << "  Reference: " << kRefTotalMass << "\n"
        << "  Computed:  " << total_mass << "\n"
        << "  Rel diff:  " << rel_diff << " (tol: " << kConservationTol << ")\n"
        << "Mass should be exactly conserved in SPH.";
}

TEST_F(ShockTube3DRegressionTest, TotalInternalEnergyRegression) {
    ASSERT_TRUE(sim_data.valid);
    double total_ie = computeWeightedSum(sim_data.mass, sim_data.ene);
    double rel_diff = std::abs(total_ie - kRefTotalIE) / kRefTotalIE;

    EXPECT_LT(rel_diff, kL2RelTol)
        << "INTERNAL ENERGY REGRESSION DETECTED!\n"
        << "  Reference: " << kRefTotalIE << "\n"
        << "  Computed:  " << total_ie << "\n"
        << "  Rel diff:  " << rel_diff << " (tol: " << kL2RelTol << ")\n"
        << "This indicates a change in energy evolution.";
}

TEST_F(ShockTube3DRegressionTest, KineticEnergyRegression) {
    ASSERT_TRUE(sim_data.valid);
    double ke = computeKineticEnergy(sim_data.mass, sim_data.vel_x, sim_data.vel_y, sim_data.vel_z);
    double rel_diff = std::abs(ke - kRefTotalKE) / kRefTotalKE;

    EXPECT_LT(rel_diff, kL2RelTol)
        << "KINETIC ENERGY REGRESSION DETECTED!\n"
        << "  Reference: " << kRefTotalKE << "\n"
        << "  Computed:  " << ke << "\n"
        << "  Rel diff:  " << rel_diff << " (tol: " << kL2RelTol << ")\n"
        << "This indicates a change in velocity evolution.";
}

// =============================================================================
// Main
// =============================================================================
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
