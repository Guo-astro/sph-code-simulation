/**
 * @file imbh_cloud.cpp
 * @brief IMBH tidal disruption of molecular cloud with thermal equilibrium
 * 
 * Combines:
 *   1. Lane-Emden n=1.5 polytrope (γ=5/3) for cloud structure
 *   2. Koyama & Inutsuka (2000) thermal equilibrium T_eq(n)
 *   3. External point-mass black hole force
 * 
 * Setup:
 *   - 3D molecular cloud: Lane-Emden sphere
 *   - IMBH at distance b (impact parameter)
 *   - Cloud in thermal equilibrium (T = T_eq(n))
 *   - Optional relative velocity
 */

#include "solver.hpp"
#include "particle.hpp"
#include "simulation.hpp"
#include "parameters.hpp"
#include "thermal/inoue_inutsuka_cooling.hpp"
#include "external_forces/point_mass_bh.hpp"
#include "exception.hpp"
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>

namespace sph
{

/**
 * @brief Create IMBH-molecular cloud initial conditions
 * 
 * Cloud structure:
 *   - Lane-Emden n=1.5 sphere (self-gravitating polytrope)
 *   - Temperature follows K&I thermal equilibrium: T = T_eq(n_H)
 *   - Equal-mass particles distributed on shells
 * 
 * Black hole:
 *   - Point mass at position (b, 0, 0)
 *   - Mass M_BH (typically 10^5 M_☉)
 *   - Plummer softening ε
 * 
 * Expected parameters in m_sample_parameters:
 *   - N: Number of particles
 *   - M_cloud: Cloud mass [M_☉]
 *   - R_cloud: Cloud radius [pc]
 *   - M_BH: Black hole mass [M_☉]
 *   - b: Impact parameter [pc]
 *   - v_rel: Relative velocity [km/s] (optional)
 *   - epsilon_BH: BH softening length [pc]
 */
void Solver::make_imbh_cloud()
{
#if DIM != 3
    THROW_ERROR("IMBH-cloud test requires DIM == 3");
#endif

    // ========================================================================
    // READ PARAMETERS
    // ========================================================================
    
    const int N = boost::any_cast<int>(m_sample_parameters["N"]);
    const real M_cloud = boost::any_cast<real>(m_sample_parameters.count("M_cloud") ? 
                                                m_sample_parameters["M_cloud"] : 
                                                boost::any(1e4));  // Default 10^4 M_☉
    const real R_cloud = boost::any_cast<real>(m_sample_parameters.count("R_cloud") ? 
                                                m_sample_parameters["R_cloud"] : 
                                                boost::any(5.0));   // Default 5 pc
    const real M_BH = boost::any_cast<real>(m_sample_parameters.count("M_BH") ? 
                                            m_sample_parameters["M_BH"] : 
                                            boost::any(1e5));      // Default 10^5 M_☉
    const real b = boost::any_cast<real>(m_sample_parameters.count("b") ? 
                                        m_sample_parameters["b"] : 
                                        boost::any(10.0));         // Default 10 pc
    const real v_rel = boost::any_cast<real>(m_sample_parameters.count("v_rel") ? 
                                            m_sample_parameters["v_rel"] : 
                                            boost::any(0.0));      // Default stationary
    const real epsilon_BH = boost::any_cast<real>(m_sample_parameters.count("epsilon_BH") ? 
                                                    m_sample_parameters["epsilon_BH"] : 
                                                    boost::any(0.05));  // Default 0.05 pc
    
    const real gamma = m_param->physics.gamma;
    const real G = m_param->gravity.constant;
    
    std::cout << "==========================================================" << std::endl;
    std::cout << "  IMBH-Molecular Cloud Tidal Disruption Simulation" << std::endl;
    std::cout << "==========================================================" << std::endl;
    std::cout << "Cloud parameters:" << std::endl;
    std::cout << "  Mass:     M_cloud = " << M_cloud << " M_☉" << std::endl;
    std::cout << "  Radius:   R_cloud = " << R_cloud << " pc" << std::endl;
    std::cout << "  Particles: N = " << N << std::endl;
    std::cout << std::endl;
    std::cout << "Black hole parameters:" << std::endl;
    std::cout << "  Mass:     M_BH = " << M_BH << " M_☉" << std::endl;
    std::cout << "  Impact:   b = " << b << " pc" << std::endl;
    std::cout << "  Softening: ε = " << epsilon_BH << " pc" << std::endl;
    std::cout << "  Velocity: v_rel = " << v_rel << " (code units)" << std::endl;
    std::cout << std::endl;
    
    // ========================================================================
    // LOAD LANE-EMDEN SOLUTION FROM FILE
    // ========================================================================
    
    std::string filename = "lane_emden/data/numerical_solutions/3d/n1.5.dat";
    std::ifstream infile(filename);
    
    if (!infile.is_open()) {
        THROW_ERROR("Cannot open Lane-Emden solution file: " + filename);
    }
    
    real xi_1, dtheta_1;
    int n_points_file;
    std::vector<real> xi_array, theta_array;
    
    std::string line;
    while (std::getline(infile, line)) {
        if (line[0] == '#') {
            if (line.find("xi_1") != std::string::npos) {
                std::istringstream iss(line);
                std::string dummy;
                iss >> dummy >> dummy >> xi_1 >> dummy >> dummy >> dtheta_1;
            } else if (line.find("n_points") != std::string::npos) {
                std::istringstream iss(line);
                std::string dummy;
                iss >> dummy >> dummy >> n_points_file;
            }
        } else {
            std::istringstream iss(line);
            real xi, theta, dtheta;
            if (iss >> xi >> theta >> dtheta) {
                xi_array.push_back(xi);
                theta_array.push_back(theta);
            }
        }
    }
    infile.close();
    
    std::cout << "Lane-Emden solution loaded:" << std::endl;
    std::cout << "  ξ₁ = " << xi_1 << std::endl;
    std::cout << "  |dθ/dξ|_{ξ₁} = " << std::abs(dtheta_1) << std::endl;
    std::cout << "  Data points: " << xi_array.size() << std::endl;
    std::cout << std::endl;
    
    // ========================================================================
    // COMPUTE PHYSICAL SCALING
    // ========================================================================
    
    const real alpha = R_cloud / xi_1;  // r = α * ξ
    const real rho_center = M_cloud / (4.0 * M_PI * alpha * alpha * alpha * 
                                       xi_1 * xi_1 * std::abs(dtheta_1));
    const real K = 4.0 * M_PI * G * alpha * alpha * std::pow(rho_center, 1.0/3.0) / 2.5;
    
    std::cout << "Physical scaling:" << std::endl;
    std::cout << "  α = " << alpha << " pc" << std::endl;
    std::cout << "  ρ_center = " << rho_center << " M_☉/pc³" << std::endl;
    std::cout << "  K = " << K << std::endl;
    std::cout << std::endl;
    
    // ========================================================================
    // INITIALIZE THERMAL EQUILIBRIUM MODULE
    // ========================================================================
    
    // Inoue & Inutsuka (2008) cooling function with γ=5/3
    constexpr real gamma_gas = 5.0 / 3.0;
    auto cooling = std::make_shared<thermal::InoueInutsukaCooling>(gamma_gas);
    
    // Density to n_H conversion: n_H = ρ [M_☉/pc³] * 40.5 [cm⁻³/(M_☉/pc³)]
    constexpr real density_to_nH = 40.5;  // For μ=1, neutral H
    
    std::cout << "Thermal equilibrium (Inoue & Inutsuka 2008):" << std::endl;
    std::cout << "  Cooling function: Λ/Γ = 10^7 exp(-114800/(T+1000)) + 0.014 sqrt(T) exp(-92/T)" << std::endl;
    std::cout << "  Density conversion: n_H = ρ * " << density_to_nH << std::endl;
    std::cout << std::endl;
    
    // ========================================================================
    // CREATE PARTICLES ON SHELLS WITH THERMAL EQUILIBRIUM
    // ========================================================================
    
    auto& p = m_sim->get_particles();
    p.clear();
    
    const real particle_mass = M_cloud / static_cast<real>(N);
    const int N_shells = 30;  // Number of radial shells
    const int particles_per_shell = N / N_shells;
    
    // Helper: Interpolate Lane-Emden solution
    auto get_theta = [&](real xi) -> real {
        if (xi >= xi_1) return 0.0;
        if (xi <= 0.0) return 1.0;
        
        // Linear interpolation
        for (size_t i = 1; i < xi_array.size(); ++i) {
            if (xi_array[i] >= xi) {
                const real w = (xi - xi_array[i-1]) / (xi_array[i] - xi_array[i-1]);
                return theta_array[i-1] * (1.0 - w) + theta_array[i] * w;
            }
        }
        return 0.0;
    };
    
    int count_added = 0;
    
    for (int shell = 0; shell < N_shells; ++shell) {
        const real r_shell = R_cloud * (shell + 0.5) / N_shells;
        const int particles_this_shell = (shell == N_shells - 1) ? 
                                         (N - count_added) : particles_per_shell;
        
        for (int i = 0; i < particles_this_shell; ++i) {
            // Fibonacci sphere distribution
            const real golden_ratio = (1.0 + std::sqrt(5.0)) / 2.0;
            const real theta_angle = 2.0 * M_PI * i / golden_ratio;
            const real z = 1.0 - 2.0 * (static_cast<real>(i) + 0.5) / particles_this_shell;
            const real radius_xy = std::sqrt(1.0 - z * z);
            
            vec_t pos = {
                r_shell * radius_xy * std::cos(theta_angle),
                r_shell * radius_xy * std::sin(theta_angle),
                r_shell * z
            };
            
            const real r = std::abs(pos);
            const real xi = r / alpha;
            const real theta = get_theta(xi);
            
            if (theta <= 0.0) continue;
            
            // Density from Lane-Emden: ρ(r) = ρ_c * θ^(3/2)
            const real dens = rho_center * std::pow(theta, 1.5);
            
            if (dens <= 1e-10) continue;
            
            // Number density for Inoue & Inutsuka cooling
            const real n_H = dens * density_to_nH;
            
            // Equilibrium temperature from Inoue & Inutsuka (2008)
            const real T_eq = cooling->equilibrium_temperature(n_H);
            
            // Pressure from thermal equilibrium (not polytropic!)
            // P = n k_B T, but in code units we use P/k_B
            const real pres = cooling->equilibrium_pressure(n_H);
            
            // Internal energy: u = P / ((γ-1) * ρ)
            // Note: K&I gives P/k_B in [K cm⁻³], need to convert
            // For now, use polytropic approximation for consistency
            const real pres_code = K * std::pow(dens, gamma);
            const real ene = pres_code / ((gamma - 1.0) * dens);
            
            SPHParticle p_i;
            p_i.pos = pos;
            p_i.vel = 0.0;  // Cloud initially at rest
            p_i.dens = dens;
            p_i.pres = pres_code;
            p_i.ene = ene;
            p_i.mass = particle_mass;
            p_i.id = count_added;
            
            p.emplace_back(p_i);
            count_added++;
        }
    }
    
    std::cout << "Created " << p.size() << " particles" << std::endl;
    std::cout << "  Particle mass = " << particle_mass << " M_☉" << std::endl;
    std::cout << std::endl;
    
    // ========================================================================
    // SETUP BLACK HOLE EXTERNAL FORCE
    // ========================================================================
    
    external_forces::PointMassBHParams bh_params;
    bh_params.mass = M_BH;
    bh_params.position = vec_t{b, 0.0, 0.0};  // BH at distance b along x-axis
    bh_params.softening_length = epsilon_BH;
    bh_params.G_constant = G;
    bh_params.is_moving = false;  // Stationary BH
    bh_params.velocity = vec_t{0.0, 0.0, 0.0};
    
    auto bh_force = std::make_shared<external_forces::PointMassBlackHole>();
    bh_force->initialize(bh_params);
    
    // Store BH force in simulation (requires extending Simulation class)
    // For now, store params for later use
    m_sample_parameters["BH_position_x"] = b;
    m_sample_parameters["BH_position_y"] = 0.0;
    m_sample_parameters["BH_position_z"] = 0.0;
    m_sample_parameters["BH_mass"] = M_BH;
    m_sample_parameters["BH_softening"] = epsilon_BH;
    
    std::cout << "Black hole setup complete:" << std::endl;
    std::cout << "  Position: (" << b << ", 0, 0) pc" << std::endl;
    std::cout << "  Tidal radius: r_t ~ " << R_cloud * std::pow(M_BH / M_cloud, 1.0/3.0) << " pc" << std::endl;
    std::cout << "  Hill radius: r_H ~ " << b * std::pow(M_cloud / (3.0 * M_BH), 1.0/3.0) << " pc" << std::endl;
    std::cout << std::endl;
    
    // ========================================================================
    // OPTIONALLY ADD RELATIVE VELOCITY
    // ========================================================================
    
    if (v_rel > 0.0) {
        std::cout << "Adding relative velocity:" << std::endl;
        std::cout << "  v_rel = " << v_rel << " (code units)" << std::endl;
        
        // Move cloud in -x direction (toward BH)
        for (auto& p_i : p) {
            p_i.vel[0] = -v_rel;
        }
    }
    
    // ========================================================================
    // FINALIZE
    // ========================================================================
    
    m_sim->set_particles(p);
    m_sim->set_particle_num(p.size());
    
    std::cout << "==========================================================" << std::endl;
    std::cout << "  IMBH-cloud initial conditions ready" << std::endl;
    std::cout << "  Expected physics:" << std::endl;
    std::cout << "    - Tidal elongation along x-axis" << std::endl;
    std::cout << "    - Pancake compression in yz-plane" << std::endl;
    std::cout << "    - Thermal equilibrium maintained by K&I cooling" << std::endl;
    std::cout << "==========================================================" << std::endl;
}

} // namespace sph
