#include <cassert>
#include <algorithm>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <chrono>

#include "solver.hpp"
#include "parameters.hpp"
#include "particle.hpp"
#include "logger.hpp"
#include "exception.hpp"
#include "output_manager.hpp"
#include "output_metadata.hpp"
#include "simulation.hpp"
#include "periodic.hpp"
#include "bhtree.hpp"

#include <nlohmann/json.hpp>

// modules
#include "timestep.hpp"
#include "pre_interaction.hpp"
#include "fluid_force.hpp"
#include "gravity_force.hpp"
#include "external_forces/point_mass_bh.hpp"
#include "disph/d_pre_interaction.hpp"
#include "disph/d_fluid_force.hpp"
#include "gsph/g_pre_interaction.hpp"
#include "gsph/g_fluid_force.hpp"
#include "gdisph/gd_pre_interaction.hpp"
#include "gdisph/gd_fluid_force.hpp"
#include "srgsph/sr_pre_interaction.hpp"
#include "srgsph/sr_fluid_force.hpp"
#include "srgsph/sr_timestep.hpp"
#include "srgsph/sr_primitive_recovery.hpp"
#include "grgsph/gr_fluid_force.hpp"
#include "grgsph/gr_pre_interaction.hpp"
#include "grgsph/gr_metric.hpp"
#include "gspmhd/mhd_fluid_force.hpp"
#include "gspmhd/mhd_pre_interaction.hpp"
#include "srmhd/srmhd_fluid_force.hpp"
#include "srmhd/srmhd_pre_interaction.hpp"
#include "srmhd/srmhd_timestep.hpp"
#include "srmhd/srmhd_primitive_recovery.hpp"

// relaxation
#include "relaxation/lane_emden_relaxation.hpp"

// ghost particle boundaries
#include "shock_tube_ghost.hpp"
#include "relaxation/polytropic_slab_relaxation.hpp"
#include "relaxation/polytropic_slab_2d_relaxation.hpp"
#include "relaxation/koyama_inutsuka_relaxation.hpp"
#include "relaxation/isothermal_relaxation.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace sph
{

Solver::Solver(int argc, char * argv[])
{
    std::cout << "--------------SPH simulation-------------\n\n";
    if(argc == 1) {
        std::cerr << "how to use\n" << std::endl;
        std::cerr << "sph <paramter.json>" << std::endl;
        std::exit(EXIT_FAILURE);
    } else {
        read_parameterfile(argv[1]);
    }

    Logger::open(m_output_dir);

#ifdef _OPENMP
    WRITE_LOG << "Open MP is valid.";
    int num_threads;
    if(argc == 3) {
        num_threads = std::atoi(argv[2]);
        omp_set_num_threads(num_threads);
    } else {
        num_threads = omp_get_max_threads();
    }
    WRITE_LOG << "the number of threads = " << num_threads << "\n";
#else
    WRITE_LOG << "OpenMP is invalid.\n";
#endif
    WRITE_LOG << "parameters";

    WRITE_LOG << "output directory     = " << m_output_dir;

    WRITE_LOG << "time";
    WRITE_LOG << "* start time         = " << m_param->time.start;
    WRITE_LOG << "* end time           = " << m_param->time.end;
    WRITE_LOG << "* output time        = " << m_param->time.output;
    WRITE_LOG << "* enerty output time = " << m_param->time.energy;

    switch(m_param->type) {
    case SPHType::SSPH:
        WRITE_LOG << "SPH type: Standard SPH";
        break;
    case SPHType::DISPH:
        WRITE_LOG << "SPH type: Density Independent SPH";
        break;
    case SPHType::GSPH:
        if(m_param->gsph.is_2nd_order) {
            WRITE_LOG << "SPH type: Godunov SPH (2nd order)";
        } else {
            WRITE_LOG << "SPH type: Godunov SPH (1st order)";
        }
        if(m_param->gsph.riemann_solver == RiemannSolverType::ITERATIVE) {
            WRITE_LOG << "* Riemann solver: Iterative (van Leer 1997)";
        } else if(m_param->gsph.riemann_solver == RiemannSolverType::KITAJIMA) {
            WRITE_LOG << "* Riemann solver: Kitajima-style (Newton-Raphson)";
        } else {
            WRITE_LOG << "* Riemann solver: HLL";
        }
        break;
    case SPHType::GDISPH:
        if(m_param->gsph.is_2nd_order) {
            WRITE_LOG << "SPH type: Godunov Density-Independent SPH (2nd order)";
        } else {
            WRITE_LOG << "SPH type: Godunov Density-Independent SPH (1st order)";
        }
        break;
    case SPHType::SRGSPH:
        if(m_param->srgsph.is_2nd_order) {
            WRITE_LOG << "SPH type: Special Relativistic Godunov SPH (2nd order)";
        } else {
            WRITE_LOG << "SPH type: Special Relativistic Godunov SPH (1st order)";
        }
        WRITE_LOG << "* Speed of light c = " << m_param->srgsph.c_speed;
        break;
    case SPHType::GRGSPH:
        if(m_param->srgsph.is_2nd_order) {
            WRITE_LOG << "SPH type: General Relativistic Godunov SPH (2nd order)";
        } else {
            WRITE_LOG << "SPH type: General Relativistic Godunov SPH (1st order)";
        }
        WRITE_LOG << "* Speed of light c = " << m_param->srgsph.c_speed;
        break;
    case SPHType::GSPMHD:
        if(m_param->mhd.is_2nd_order) {
            WRITE_LOG << "SPH type: Godunov SPMHD (2nd order) - Iwasaki & Inutsuka (2011)";
        } else {
            WRITE_LOG << "SPH type: Godunov SPMHD (1st order) - Iwasaki & Inutsuka (2011)";
        }
        WRITE_LOG << "* Powell div-B correction: " << (m_param->mhd.use_powell_correction ? "enabled" : "disabled");
        break;
    case SPHType::SRMHD:
        if(m_param->srgsph.is_2nd_order) {
            WRITE_LOG << "SPH type: Special Relativistic MHD (2nd order) - SR-GSPH + GSPMHD";
        } else {
            WRITE_LOG << "SPH type: Special Relativistic MHD (1st order) - SR-GSPH + GSPMHD";
        }
        WRITE_LOG << "* Speed of light c = " << m_param->srgsph.c_speed;
        WRITE_LOG << "* Powell div-B correction: " << (m_param->mhd.use_powell_correction ? "enabled" : "disabled");
        break;
    }

    WRITE_LOG << "CFL condition";
    WRITE_LOG << "* sound speed = " << m_param->cfl.sound;
    WRITE_LOG << "* force       = " << m_param->cfl.force;

    WRITE_LOG << "Artificial Viscosity";
    WRITE_LOG << "* alpha = " << m_param->av.alpha;
    if(m_param->av.use_balsara_switch) {
        WRITE_LOG << "* use Balsara switch";
    }
    if(m_param->av.use_time_dependent_av) {
        WRITE_LOG << "* use time dependent AV";
        WRITE_LOG << "* alpha max = " << m_param->av.alpha_max;
        WRITE_LOG << "* alpha min = " << m_param->av.alpha_min;
        WRITE_LOG << "* epsilon   = " << m_param->av.epsilon;
    }

    if(m_param->ac.is_valid) {
        WRITE_LOG << "Artificial Conductivity";
        WRITE_LOG << "* alpha = " << m_param->ac.alpha;
    }

    WRITE_LOG << "Tree";
    WRITE_LOG << "* max tree level       = " << m_param->tree.max_level;
    WRITE_LOG << "* leaf particle number = " << m_param->tree.leaf_particle_num;

    WRITE_LOG << "Physics";
    WRITE_LOG << "* Neighbor number = " << m_param->physics.neighbor_number;
    WRITE_LOG << "* gamma           = " << m_param->physics.gamma;
    WRITE_LOG << "* C_smooth        = " << m_param->physics.c_smooth;

    WRITE_LOG << "Kernel";
    if(m_param->kernel == KernelType::CUBIC_SPLINE) {
        WRITE_LOG << "* Cubic Spline";
    } else if(m_param->kernel == KernelType::WENDLAND) {
        WRITE_LOG << "* Wendland";
    } else if(m_param->kernel == KernelType::GAUSSIAN) {
        WRITE_LOG << "* Gaussian";
    } else {
        THROW_ERROR("kernel is unknown.");
    }

    if(m_param->iterative_sml) {
        WRITE_LOG << "Iterative calculation for smoothing length is valid.";
    }

    if(m_param->periodic.is_valid) {
        WRITE_LOG << "Periodic boundary condition is valid.";
    }
    
    if(m_param->gravity.is_valid) {
        WRITE_LOG << "Gravity is valid.";
        WRITE_LOG << "G     = " << m_param->gravity.constant;
        WRITE_LOG << "theta = " << m_param->gravity.theta;
    }

    switch(m_sample) {
#define WRITE_SAMPLE(a, b) case a: WRITE_LOG << "Sample: " b " test"; break
        WRITE_SAMPLE(Sample::ShockTube, "shock tube");
        WRITE_SAMPLE(Sample::GreshoChanVortex, "Gresho-Chan vortex");
        WRITE_SAMPLE(Sample::PairingInstability, "Pairing Instability");
        WRITE_SAMPLE(Sample::HydroStatic, "Hydro static");
        WRITE_SAMPLE(Sample::KHI, "Kelvin-Helmholtz Instability");
        WRITE_SAMPLE(Sample::ISMCooling1D, "ISM Cooling 1D (Koyama & Inutsuka 2000)");
        WRITE_SAMPLE(Sample::Evrard, "Evrard collapse");
        WRITE_SAMPLE(Sample::EvrardColdCollapse, "Evrard cold collapse (demonstrates shock amplification)");
        WRITE_SAMPLE(Sample::LaneEmden, "Lane-Emden hydrostatic");
        WRITE_SAMPLE(Sample::LaneEmdenCylinder, "Lane-Emden cylinder (3D, radial gravity in xy-plane)");
        WRITE_SAMPLE(Sample::PolytropicSlab2D, "Polytropic slab 2D (planar, gravity in y-direction)");
        WRITE_SAMPLE(Sample::Sedov, "Sedov blast wave");
        WRITE_SAMPLE(Sample::IsothermalBonnorEbert, "Isothermal Bonnor-Ebert (self-gravitating)");
        WRITE_SAMPLE(Sample::HVCCIsothermal10K, "HVCC 10K isothermal (Oka 2017)");
        WRITE_SAMPLE(Sample::MHDShockTube1, "MHD Shock Tube 1 (Dai-Woodward)");
        WRITE_SAMPLE(Sample::MHDShockTube2, "MHD Shock Tube 2 (Strong shock)");
        WRITE_SAMPLE(Sample::SRMHDBalsara1, "SRMHD Balsara Test 1 (Relativistic MHD shock)");
        WRITE_SAMPLE(Sample::UniformCloud, "Uniform cloud (IMBH-cloud tidal interaction)");
#undef WRITE_SAMPLE
        default:
            break;
    }

    WRITE_LOG;

    // Initialize OutputManager (will be configured in read_parameterfile)
    m_snapshot_counter = 0;
}

void Solver::read_parameterfile(const char * filename)
{
    namespace pt = boost::property_tree;

    m_param = std::make_shared<SPHParameters>();

    pt::ptree input;

    std::string name_str = filename;
    std::cout << "read_parameterfile: filename = '" << name_str << "'" << std::endl;
    std::cout.flush();
    if(name_str == "shock_tube") {
        pt::read_json("simulations/benchmarks/shock_tube/shock_tube.json", input);
        m_sample = Sample::ShockTube;
        m_sample_parameters["N"] = input.get<int>("N", 100);
    } else if(name_str == "shock_tube_2d") {
        pt::read_json("simulations/benchmarks/shock_tube_2d/shock_tube_2d.json", input);
        m_sample = Sample::ShockTube2D;
        m_sample_parameters["Nx"] = input.get<int>("Nx", 200);
        m_sample_parameters["Ny"] = input.get<int>("Ny", 40);
    } else if(name_str == "vacuum") {
        pt::read_json("simulations/benchmarks/vacuum/vacuum.json", input);
        m_sample = Sample::Vacuum;
        m_sample_parameters["N"] = input.get<int>("N", 800);
    } else if(name_str == "strong_shock") {
        pt::read_json("simulations/benchmarks/strong_shock/strong_shock.json", input);
        m_sample = Sample::StrongShock;
        m_sample_parameters["N"] = input.get<int>("N", 800);
    } else if(name_str == "gresho_chan_vortex") {
        pt::read_json("simulations/stability/gresho_chan_vortex/gresho_chan_vortex.json", input);
        m_sample = Sample::GreshoChanVortex;
        m_sample_parameters["N"] = input.get<int>("N", 64);
    } else if(name_str == "pairing_instability") {
        pt::read_json("simulations/stability/pairing_instability/pairing_instability.json", input);
        m_sample = Sample::PairingInstability;
        m_sample_parameters["N"] = input.get<int>("N", 64);
    } else if(name_str == "hydrostatic") {
        pt::read_json("simulations/stability/hydrostatic/hydrostatic.json", input);
        m_sample = Sample::HydroStatic;
        m_sample_parameters["N"] = input.get<int>("N", 32);
    } else if(name_str == "gradh_study") {
        pt::read_json("simulations/stability/gradh_study/gradh_study.json", input);
        m_sample = Sample::HydroStatic;  // Uses HydroStatic or resumeFromSnapshot
        m_sample_parameters["N"] = input.get<int>("N", 64);
    } else if(name_str == "khi") {
        pt::read_json("simulations/stability/khi/khi.json", input);
        m_sample = Sample::KHI;
        m_sample_parameters["N"] = input.get<int>("N", 128);
    } else if(name_str == "ism_cooling_1d") {
        pt::read_json("simulations/astrophysics/cooling_heating/ism_cooling_1d.json", input);
        m_sample = Sample::ISMCooling1D;
        m_sample_parameters["N"] = input.get<int>("N", 1000);
        m_sample_parameters["n_min"] = input.get<real>("n_min", 0.1);
        m_sample_parameters["n_max"] = input.get<real>("n_max", 1000.0);
        m_sample_parameters["T_init"] = input.get<real>("T_init", 8000.0);
    } else if(name_str == "evrard") {
        pt::read_json("simulations/astrophysics/evrard/evrard.json", input);
        m_sample = Sample::Evrard;
        m_sample_parameters["N"] = input.get<int>("N", 20);
    } else if(name_str == "evrard_cold_collapse") {
        pt::read_json("simulations/astrophysics/evrard_cold_collapse/evrard_cold_collapse.json", input);
        m_sample = Sample::EvrardColdCollapse;
        m_sample_parameters["N"] = input.get<int>("N", 30);
    } else if(name_str == "lane_emden") {
        pt::read_json("lane_emden/lane_emden.json", input);
        m_sample = Sample::LaneEmden;
        m_sample_parameters["N"] = input.get<int>("N", 30);
    } else if(name_str == "sedov") {
        pt::read_json("simulations/benchmarks/sedov/sedov.json", input);
        m_sample = Sample::Sedov;
        m_sample_parameters["N"] = input.get<int>("N", 100);
    } else if(name_str == "sr_sod") {
        pt::read_json("simulations/relativistic/sr_sod/sr_sod.json", input);
        m_sample = Sample::SRSod;
        m_sample_parameters["N"] = input.get<int>("N", 50);
        m_sample_parameters["different_nu"] = input.get<bool>("different_nu", false);
        m_sample_parameters["testType"] = input.get<std::string>("testType", "sod");
        m_sample_parameters["particleMode"] = input.get<std::string>("particleMode", "equal_N");
    } else {
        pt::read_json(filename, input);
        
        // Try to read explicit sample type from JSON
        std::string sample_type = input.get<std::string>("sample", "");
        if (sample_type == "shock_tube") {
            m_sample = Sample::ShockTube;
            m_sample_parameters["N"] = input.get<int>("N", 100);
        } else if (sample_type == "shock_tube_2d") {
            m_sample = Sample::ShockTube2D;
            m_sample_parameters["Nx"] = input.get<int>("Nx", 200);
            m_sample_parameters["Ny"] = input.get<int>("Ny", 40);
            m_sample_parameters["Ly"] = input.get<real>("Ly", 0.2);  // Default 2D shock tube domain
        } else if (sample_type == "shock_tube_3d") {
            m_sample = Sample::ShockTube3D;
            m_sample_parameters["Nx"] = input.get<int>("Nx", 100);
            m_sample_parameters["Ny"] = input.get<int>("Ny", 10);
            m_sample_parameters["Nz"] = input.get<int>("Nz", 10);
            m_sample_parameters["Ly"] = input.get<real>("Ly", 0.1);
            m_sample_parameters["Lz"] = input.get<real>("Lz", 0.1);
        } else if (sample_type == "shock_tube_3d_cubic") {
            m_sample = Sample::ShockTube3DCubic;
            m_sample_parameters["Nx"] = input.get<int>("Nx", 100);
        } else if (sample_type == "vacuum") {
            m_sample = Sample::Vacuum;
            m_sample_parameters["N"] = input.get<int>("N", 800);
        } else if (sample_type == "strong_shock") {
            m_sample = Sample::StrongShock;
            m_sample_parameters["N"] = input.get<int>("N", 800);
        } else if (sample_type == "sedov") {
            m_sample = Sample::Sedov;
            m_sample_parameters["N"] = input.get<int>("N", 30);
        } else if (sample_type == "lane_emden") {
            m_sample = Sample::LaneEmden;
            m_sample_parameters["N"] = input.get<int>("N", 30);
            // Read Lane-Emden relaxation parameters if provided
            if(auto opt = input.get_child_optional("alpha_scaling")) {
                m_sample_parameters["alpha"] = input.get<real>("alpha_scaling");
            }
            if(auto opt = input.get_child_optional("rho_center")) {
                m_sample_parameters["rho_center"] = input.get<real>("rho_center");
            }
            if(auto opt = input.get_child_optional("K")) {
                m_sample_parameters["K"] = input.get<real>("K");
            }
            if(auto opt = input.get_child_optional("R")) {
                m_sample_parameters["R"] = input.get<real>("R");
            }
            if(auto opt = input.get_child_optional("M_total")) {
                m_sample_parameters["M_total"] = input.get<real>("M_total");
            }
        } else if (sample_type == "sr_tangent_velocity") {
            m_sample = Sample::SRTangentVelocity;
            m_sample_parameters["N"] = input.get<int>("N", 1600);
            m_sample_parameters["vt_left"] = input.get<real>("vt_left", real(0.9));
            m_sample_parameters["vt_right"] = input.get<real>("vt_right", real(0.9));
            m_sample_parameters["useGhostParticles"] = input.get<bool>("useGhostParticles", true);
            m_sample_parameters["ghostLayers"] = input.get<int>("ghostLayers", 6);
            // Non-uniform resolution support (Kitajima hypothesis: rarefaction side matters)
            if (auto opt = input.get_child_optional("N_left")) {
                m_sample_parameters["N_left"] = input.get<int>("N_left");
            }
            if (auto opt = input.get_child_optional("N_right")) {
                m_sample_parameters["N_right"] = input.get<int>("N_right");
            }
            if (auto opt = input.get_child_optional("leftResolutionRatio")) {
                m_sample_parameters["leftResolutionRatio"] = input.get<real>("leftResolutionRatio");
            }
        } else if (sample_type == "sr_rosswog") {
            m_sample = Sample::SRRosswog;
            // Rosswog (2010) arXiv:0907.4890 benchmark tests
            // testType: perturbed_shock_tube (test3), sine_advection (test5),
            //           square_advection (test6), simple_wave (test7)
            m_sample_parameters["testType"] = input.get<std::string>("testType", "perturbed_shock_tube");
            m_sample_parameters["N"] = input.get<int>("N", 500);
            m_sample_parameters["gamma"] = input.get<real>("gamma", real(5.0 / 3.0));
            m_sample_parameters["advectionVelocity"] = input.get<real>("advectionVelocity", real(0.997));
            m_sample_parameters["vMax"] = input.get<real>("vMax", real(0.7));
        } else if (sample_type == "sr_sod") {
            m_sample = Sample::SRSod;
            // SR-GSPH / GR-GSPH shock tube tests
            // testType: sod, blast_wave, strong_blast
            m_sample_parameters["N"] = input.get<int>("N", 50);
            m_sample_parameters["different_nu"] = input.get<bool>("different_nu", false);
            m_sample_parameters["testType"] = input.get<std::string>("testType", "sod");
            m_sample_parameters["particleMode"] = input.get<std::string>("particleMode", "equal_N");
            m_sample_parameters["v_left"] = input.get<real>("v_left", real(0.0));
        } else if (sample_type == "gr_schwarzschild_shock") {
            m_sample = Sample::GRSchwarzschildShock;
            // GR-GSPH Schwarzschild radial shock tube (Liptai & Price 2019)
            m_sample_parameters["N"] = input.get<int>("N", 1000);
            m_sample_parameters["blackHoleMass"] = input.get<real>("blackHoleMass", real(1.0));
            m_sample_parameters["rMin"] = input.get<real>("rMin", real(3.0));  // In units of M
            m_sample_parameters["rMax"] = input.get<real>("rMax", real(15.0)); // In units of M
            m_sample_parameters["rDiscontinuity"] = input.get<real>("rDiscontinuity", real(6.0)); // In units of M
            m_sample_parameters["ghostLayers"] = input.get<int>("ghostLayers", 6);
        } else if (sample_type == "gr_geodesic_test") {
            m_sample = Sample::GRGeodesicTest;
            // GR geodesic test (Liptai & Price 2019 Section 6.2)
            m_sample_parameters["N"] = input.get<int>("N", 10);
            m_sample_parameters["blackHoleMass"] = input.get<real>("blackHoleMass", real(1.0));
            m_sample_parameters["testType"] = input.get<std::string>("testType", "radial_infall");
            m_sample_parameters["rMin"] = input.get<real>("rMin", real(4.0));   // For radial infall
            m_sample_parameters["rMax"] = input.get<real>("rMax", real(40.0));  // For radial infall
            m_sample_parameters["rOrbit"] = input.get<real>("rOrbit", real(10.0)); // For circular orbit
        } else if (sample_type == "gr_bondi") {
            m_sample = Sample::GRBondi;
            // GR-GSPH Bondi accretion (Liptai & Price 2019 Figs 14, 18)
            m_sample_parameters["N"] = input.get<int>("N", 500);
            m_sample_parameters["blackHoleMass"] = input.get<real>("blackHoleMass", real(1.0));
            m_sample_parameters["testType"] = input.get<std::string>("testType", "geodesic");
            m_sample_parameters["rIn"] = input.get<real>("rIn", real(2.5));   // Inner boundary
            m_sample_parameters["rOut"] = input.get<real>("rOut", real(18.0)); // Outer boundary
            m_sample_parameters["ghostLayers"] = input.get<int>("ghostLayers", 6);
        } else if (sample_type == "ns_merger_2d") {
            m_sample = Sample::NSMerger2D;
            m_sample_parameters["R_star"] = input.get<real>("ns_merger.star1.radius", real(1.2));
            m_sample_parameters["rho_c"] = input.get<real>("ns_merger.star1.central_density", real(2.8));
            m_sample_parameters["separation"] = input.get<real>("ns_merger.separation", real(6.0));
            m_sample_parameters["v_collision"] = input.get<real>("ns_merger.star1.velocity_x", real(0.15));
            m_sample_parameters["n_radial"] = input.get<int>("ns_merger.star1.n_particles_radial", 30);
        } else if (sample_type == "bns_cocoon_1d") {
            m_sample = Sample::BNSCocoon1D;
            // Ejecta parameters (Gutiérrez+2024 arXiv:2408.15973v3)
            m_sample_parameters["M_ej"] = input.get<real>("bns_cocoon.M_ej", real(0.01));
            m_sample_parameters["v_min"] = input.get<real>("bns_cocoon.v_min", real(0.05));
            m_sample_parameters["v_max"] = input.get<real>("bns_cocoon.v_max", real(0.7));
            m_sample_parameters["v_break"] = input.get<real>("bns_cocoon.v_break", real(0.5));
            m_sample_parameters["beta"] = input.get<real>("bns_cocoon.beta", real(0.5));
            m_sample_parameters["alpha_tail"] = input.get<real>("bns_cocoon.alpha_tail", real(6.0));
            m_sample_parameters["t0"] = input.get<real>("bns_cocoon.t0", real(1.0));
            m_sample_parameters["gamma"] = input.get<real>("bns_cocoon.gamma", real(1.333333333));
            m_sample_parameters["rho_floor"] = input.get<real>("bns_cocoon.rho_floor", real(1e-12));
            m_sample_parameters["profile_type"] = input.get<int>("bns_cocoon.profile_type", 1);
            // Cocoon parameters
            m_sample_parameters["E_cocoon"] = input.get<real>("bns_cocoon.E_cocoon", real(1e-4));
            m_sample_parameters["r_cocoon_frac"] = input.get<real>("bns_cocoon.r_cocoon_frac", real(0.3));
            m_sample_parameters["Gamma_cocoon"] = input.get<real>("bns_cocoon.Gamma_cocoon", real(2.0));
            m_sample_parameters["n_particles"] = input.get<int>("bns_cocoon.n_particles", 500);
        } else if (sample_type == "bns_cocoon_2d") {
            m_sample = Sample::BNSCocoon2D;
            // Ejecta parameters (Gutiérrez+2024 arXiv:2408.15973v3)
            m_sample_parameters["M_ej"] = input.get<real>("bns_cocoon.M_ej", real(0.01));
            m_sample_parameters["v_min"] = input.get<real>("bns_cocoon.v_min", real(0.05));
            m_sample_parameters["v_max"] = input.get<real>("bns_cocoon.v_max", real(0.7));
            m_sample_parameters["v_break"] = input.get<real>("bns_cocoon.v_break", real(0.5));
            m_sample_parameters["beta"] = input.get<real>("bns_cocoon.beta", real(0.5));
            m_sample_parameters["alpha_tail"] = input.get<real>("bns_cocoon.alpha_tail", real(6.0));
            m_sample_parameters["t0"] = input.get<real>("bns_cocoon.t0", real(1.0));
            m_sample_parameters["gamma"] = input.get<real>("bns_cocoon.gamma", real(1.333333333));
            m_sample_parameters["rho_floor"] = input.get<real>("bns_cocoon.rho_floor", real(1e-12));
            m_sample_parameters["profile_type"] = input.get<int>("bns_cocoon.profile_type", 1);
            // Angular distribution
            m_sample_parameters["theta_polar"] = input.get<real>("bns_cocoon.theta_polar", real(0.2618));
            m_sample_parameters["polar_density_factor"] = input.get<real>("bns_cocoon.polar_density_factor", real(0.1));
            m_sample_parameters["angular_power"] = input.get<real>("bns_cocoon.angular_power", real(2.0));
            // Cocoon parameters
            m_sample_parameters["E_cocoon"] = input.get<real>("bns_cocoon.E_cocoon", real(1e-4));
            m_sample_parameters["r_cocoon_frac"] = input.get<real>("bns_cocoon.r_cocoon_frac", real(0.3));
            m_sample_parameters["theta_jet"] = input.get<real>("bns_cocoon.theta_jet", real(0.2094));
            m_sample_parameters["Gamma_cocoon"] = input.get<real>("bns_cocoon.Gamma_cocoon", real(2.0));
            m_sample_parameters["jet_injection"] = input.get<bool>("bns_cocoon.jet_injection", false);
            // Resolution
            m_sample_parameters["n_radial"] = input.get<int>("bns_cocoon.n_radial", 200);
            m_sample_parameters["n_angular"] = input.get<int>("bns_cocoon.n_angular", 160);
        } else if (sample_type == "isothermal_slab") {
            m_sample = Sample::IsothermalSlab;
            // 1D isothermal self-gravitating slab for diffusive instability study
            m_sample_parameters["N"] = input.get<int>("N", 200);
            m_sample_parameters["rho_center"] = input.get<real>("rho_center", real(1.0));
            m_sample_parameters["sound_speed"] = input.get<real>("sound_speed", real(1.0));
            m_sample_parameters["scale_height"] = input.get<real>("scale_height", real(0.5));
        } else if (sample_type == "polytropic_slab") {
            m_sample = Sample::PolytropicSlab;
            // 1D polytropic self-gravitating slab using planar Lane-Emden solution
            m_sample_parameters["N"] = input.get<int>("N", 200);
            m_sample_parameters["rho_center"] = input.get<real>("rho_center", real(1.0));
            m_sample_parameters["K"] = input.get<real>("K", real(1.0));  // Polytropic constant
        } else if (sample_type == "polytropic_slab_2d") {
            m_sample = Sample::PolytropicSlab2D;
            // 2D planar Lane-Emden slab (gravity in y-direction)
            m_sample_parameters["N"] = input.get<int>("N", 50);  // N² particles
            m_sample_parameters["rho_center"] = input.get<real>("rho_center", real(1.0));
            m_sample_parameters["K"] = input.get<real>("K", real(1.0));  // Polytropic constant
            m_sample_parameters["L_x"] = input.get<real>("L_x", real(1.0));  // Domain width in x
        } else if (sample_type == "lane_emden_cylinder") {
            m_sample = Sample::LaneEmdenCylinder;
            // 3D cylindrical Lane-Emden (radial gravity in xy-plane)
            m_sample_parameters["N"] = input.get<int>("N", 30);  // N³ particles (nominal)
            m_sample_parameters["R"] = input.get<real>("R", real(1.0));  // Cylinder radius
            m_sample_parameters["L_z"] = input.get<real>("L_z", real(1.0));  // Cylinder length
            m_sample_parameters["M_total"] = input.get<real>("M_total", real(1.0));  // Total mass
        } else if (sample_type == "sinusoidal_perturbation") {
            m_sample = Sample::SinusoidalPerturbation;
            // Periodic box with sinusoidal density perturbation for diffusion study
            m_sample_parameters["N"] = input.get<int>("N", 200);
            m_sample_parameters["wavelength"] = input.get<real>("wavelength", real(1.0));
            m_sample_parameters["amplitude"] = input.get<real>("amplitude", real(0.1));
            m_sample_parameters["rho_mean"] = input.get<real>("rho_mean", real(1.0));
            m_sample_parameters["pressure"] = input.get<real>("pressure", real(1.0));
        } else if (sample_type == "jeans_instability") {
            m_sample = Sample::JeansInstability;
            // Periodic box with sinusoidal density perturbation for Jeans instability
            m_sample_parameters["N"] = input.get<int>("N", 200);
            m_sample_parameters["wavelength"] = input.get<real>("wavelength", real(1.0));
            m_sample_parameters["amplitude"] = input.get<real>("amplitude", real(0.01));
            m_sample_parameters["rho_0"] = input.get<real>("rho_0", real(1.0));
            m_sample_parameters["c_s"] = input.get<real>("c_s", real(1.0));
        } else if (sample_type == "bonnor_ebert_ki2000") {
            m_sample = Sample::BonnorEbertKI2000;
            // K&I 2000 pressure-truncated Bonnor-Ebert sphere
            // Parameters match bonnor_ebert_ki2000.cpp sample generator
            m_sample_parameters["N"] = input.get<int>("N", 40);
            m_sample_parameters["R_cloud_pc"] = input.get<real>("R_cloud_pc", real(1.0));
            m_sample_parameters["M_cloud_Msun"] = input.get<real>("M_cloud_Msun", real(1000.0));
            m_sample_parameters["rho_center_nH"] = input.get<real>("rho_center_nH", real(100.0));
            m_sample_parameters["N_H_cm2"] = input.get<real>("N_H_cm2", real(1.0e19));
            m_sample_parameters["P_ext_K_cm3"] = input.get<real>("P_ext_K_cm3", real(1000.0));
            // Ghost envelope for pressure confinement during relaxation
            m_sample_parameters["useEnvelope"] = input.get<bool>("useEnvelope", false);
            m_sample_parameters["envelopeLayers"] = input.get<int>("envelopeLayers", 4);
        } else if (sample_type == "lane_emden_ki2000") {
            m_sample = Sample::LaneEmdenKI2000;
            // Lane-Emden density structure + K&I 2000 temperatures (NOT hydrostatic)
            m_sample_parameters["N"] = input.get<int>("N", 40);
            m_sample_parameters["R_cloud_pc"] = input.get<real>("R_cloud_pc", real(1.0));
            m_sample_parameters["M_cloud_Msun"] = input.get<real>("M_cloud_Msun", real(1000.0));
            m_sample_parameters["N_H_cm2"] = input.get<real>("N_H_cm2", real(1.0e19));
            m_sample_parameters["use_glass"] = input.get<bool>("use_glass", true);
        } else if (sample_type == "isothermal_bonnor_ebert") {
            m_sample = Sample::IsothermalBonnorEbert;
            // True isothermal Bonnor-Ebert sphere (self-gravitating hydrostatic equilibrium)
            // Solves isothermal Lane-Emden equation for stable configuration
            m_sample_parameters["N"] = input.get<int>("N", 22);
            m_sample_parameters["M_cloud"] = input.get<real>("M_cloud", real(40.0));
            m_sample_parameters["T_cloud"] = input.get<real>("T_cloud", real(10.0));
            m_sample_parameters["xi_s"] = input.get<real>("xi_s", real(6.0));  // Dimensionless truncation (critical ~6.45)
            m_sample_parameters["useEnvelope"] = input.get<bool>("useEnvelope", true);
            m_sample_parameters["envelopeLayers"] = input.get<int>("envelopeLayers", 4);
            m_sample_parameters["mu"] = input.get<real>("mu", real(1.27));  // Mean molecular weight
            // IC method: "healpix" (best), "shell" (good), or "random" (slow)
            m_sample_parameters["ic_method"] = input.get<std::string>("ic_method", std::string("healpix"));
            m_sample_parameters["use_random_ic"] = input.get<bool>("use_random_ic", false);  // Backward compat
            // Optional: specify central number density instead of mass
            // Mode 1: n_center specified -> compute M_cloud from equilibrium
            // Mode 2: M_cloud specified -> compute n_center (may be very low)
            if (input.count("n_center") > 0) {
                m_sample_parameters["n_center"] = input.get<real>("n_center");
            }
        } else if (sample_type == "hvcc_isothermal_10k") {
            m_sample = Sample::HVCCIsothermal10K;
            // 10K isothermal HVCC for IMBH-cloud interaction (Oka 2017)
            // K&I 2000 equilibrium at T=10K:
            //   N_H = 10^19 cm^-2: n_eq ~ 34 cm^-3, P/k_B ~ 920 K cm^-3
            //   N_H = 10^20 cm^-2: n_eq ~ 26 cm^-3, P/k_B ~ 710 K cm^-3
            m_sample_parameters["N"] = input.get<int>("N", 10000);
            m_sample_parameters["M_cloud"] = input.get<real>("M_cloud", real(40.0));
            m_sample_parameters["R_cloud"] = input.get<real>("R_cloud", real(0.14));
            m_sample_parameters["T_cloud"] = input.get<real>("T_cloud", real(10.0));
            m_sample_parameters["n_center"] = input.get<real>("n_center", real(1000.0));
            m_sample_parameters["N_H_cm2"] = input.get<real>("N_H_cm2", real(1.0e20));
            m_sample_parameters["M_BH"] = input.get<real>("M_BH", real(5.0e4));
            m_sample_parameters["b_impact"] = input.get<real>("b_impact", real(5.0));
            m_sample_parameters["v_rel"] = input.get<real>("v_rel", real(0.0));
            m_sample_parameters["epsilon_BH"] = input.get<real>("epsilon_BH", real(0.001));
            m_sample_parameters["useEnvelope"] = input.get<bool>("useEnvelope", true);
            m_sample_parameters["envelopeLayers"] = input.get<int>("envelopeLayers", 4);
            m_sample_parameters["P_ext"] = input.get<real>("P_ext", real(710.0));
        } else if (sample_type == "uniform_cloud") {
            m_sample = Sample::UniformCloud;
            // Uniform density cloud for IMBH-cloud tidal interaction
            // Based on observational constraints from Oka et al. (2017)
            // and first-principles analysis:
            //   M = 40 M_sun, R = 1.0-1.5 pc, T = 15 K (KI2000 equilibrium)
            //   IMBH at distance 5-10 pc creates tidal perturbations
            m_sample_parameters["N"] = input.get<int>("N", 22);  // Grid resolution (N^3 particles)
            m_sample_parameters["M_cloud"] = input.get<real>("M_cloud", real(40.0));
            m_sample_parameters["R_cloud"] = input.get<real>("R_cloud", real(1.0));
            m_sample_parameters["T_cloud"] = input.get<real>("T_cloud", real(15.0));
            m_sample_parameters["M_BH"] = input.get<real>("M_BH", real(1.0e5));
            m_sample_parameters["d_BH"] = input.get<real>("d_BH", real(5.0));
            m_sample_parameters["epsilon_BH"] = input.get<real>("epsilon_BH", real(0.01));
            m_sample_parameters["useEnvelope"] = input.get<bool>("useEnvelope", true);
            m_sample_parameters["envelopeLayers"] = input.get<int>("envelopeLayers", 4);
        } else if (sample_type == "mhd_shock_tube_1" || sample_type == "mhd_dai_woodward") {
            m_sample = Sample::MHDShockTube1;
            m_sample_parameters["N"] = input.get<int>("N", 400);
            m_sample_parameters["useGhostParticles"] = input.get<bool>("useGhostParticles", true);
            m_sample_parameters["ghostLayers"] = input.get<int>("ghostLayers", 6);
        } else if (sample_type == "mhd_shock_tube_2" || sample_type == "mhd_strong_shock") {
            m_sample = Sample::MHDShockTube2;
            m_sample_parameters["N"] = input.get<int>("N", 400);
            m_sample_parameters["useGhostParticles"] = input.get<bool>("useGhostParticles", true);
            m_sample_parameters["ghostLayers"] = input.get<int>("ghostLayers", 6);
        } else if (sample_type == "mhd_orszag_tang" || sample_type == "orszag_tang") {
            m_sample = Sample::MHDOrszagTang;
            m_sample_parameters["N"] = input.get<int>("N", 10000);
        } else if (sample_type == "srmhd_balsara_1" || sample_type == "srmhd_balsara") {
            m_sample = Sample::SRMHDBalsara1;
            m_sample_parameters["N"] = input.get<int>("N", 400);
            m_sample_parameters["useGhostParticles"] = input.get<bool>("useGhostParticles", true);
            m_sample_parameters["ghostLayers"] = input.get<int>("ghostLayers", 6);
            m_sample_parameters["useMHD"] = input.get<bool>("useMHD", true);  // Default: full MHD
        } else {
            // Try to infer sample type from SPH type and JSON content
            std::string sph_type_check = input.get<std::string>("SPHType", "");
            std::string test_type = input.get<std::string>("testType", "");
            
            if(sph_type_check == "srgsph" || sph_type_check == "grgsph") {
                // Check for BNS cocoon test types
                if(test_type == "bns_cocoon_1d" || name_str.find("bns_cocoon_1d") != std::string::npos) {
                    m_sample = Sample::BNSCocoon1D;
                    m_sample_parameters["M_ej"] = input.get<real>("bns_cocoon.M_ej", real(0.01));
                    m_sample_parameters["v_min"] = input.get<real>("bns_cocoon.v_min", real(0.05));
                    m_sample_parameters["v_max"] = input.get<real>("bns_cocoon.v_max", real(0.7));
                    m_sample_parameters["v_break"] = input.get<real>("bns_cocoon.v_break", real(0.5));
                    m_sample_parameters["beta"] = input.get<real>("bns_cocoon.beta", real(0.5));
                    m_sample_parameters["alpha_tail"] = input.get<real>("bns_cocoon.alpha_tail", real(6.0));
                    m_sample_parameters["t0"] = input.get<real>("bns_cocoon.t0", real(1.0));
                    m_sample_parameters["gamma"] = input.get<real>("bns_cocoon.gamma", real(1.333333333));
                    m_sample_parameters["rho_floor"] = input.get<real>("bns_cocoon.rho_floor", real(1e-12));
                    m_sample_parameters["profile_type"] = input.get<int>("bns_cocoon.profile_type", 1);
                    m_sample_parameters["E_cocoon"] = input.get<real>("bns_cocoon.E_cocoon", real(1e-4));
                    m_sample_parameters["r_cocoon_frac"] = input.get<real>("bns_cocoon.r_cocoon_frac", real(0.3));
                    m_sample_parameters["Gamma_cocoon"] = input.get<real>("bns_cocoon.Gamma_cocoon", real(2.0));
                    m_sample_parameters["n_particles"] = input.get<int>("bns_cocoon.n_particles", 500);
                }
                else if(test_type == "bns_cocoon_2d" || name_str.find("bns_cocoon_2d") != std::string::npos || name_str.find("bns_cocoon") != std::string::npos) {
                    m_sample = Sample::BNSCocoon2D;
                    m_sample_parameters["M_ej"] = input.get<real>("bns_cocoon.M_ej", real(0.01));
                    m_sample_parameters["v_min"] = input.get<real>("bns_cocoon.v_min", real(0.05));
                    m_sample_parameters["v_max"] = input.get<real>("bns_cocoon.v_max", real(0.7));
                    m_sample_parameters["v_break"] = input.get<real>("bns_cocoon.v_break", real(0.5));
                    m_sample_parameters["beta"] = input.get<real>("bns_cocoon.beta", real(0.5));
                    m_sample_parameters["alpha_tail"] = input.get<real>("bns_cocoon.alpha_tail", real(6.0));
                    m_sample_parameters["t0"] = input.get<real>("bns_cocoon.t0", real(1.0));
                    m_sample_parameters["gamma"] = input.get<real>("bns_cocoon.gamma", real(1.333333333));
                    m_sample_parameters["rho_floor"] = input.get<real>("bns_cocoon.rho_floor", real(1e-12));
                    m_sample_parameters["profile_type"] = input.get<int>("bns_cocoon.profile_type", 1);
                    m_sample_parameters["theta_polar"] = input.get<real>("bns_cocoon.theta_polar", real(0.2618));
                    m_sample_parameters["polar_density_factor"] = input.get<real>("bns_cocoon.polar_density_factor", real(0.1));
                    m_sample_parameters["angular_power"] = input.get<real>("bns_cocoon.angular_power", real(2.0));
                    m_sample_parameters["E_cocoon"] = input.get<real>("bns_cocoon.E_cocoon", real(1e-4));
                    m_sample_parameters["r_cocoon_frac"] = input.get<real>("bns_cocoon.r_cocoon_frac", real(0.3));
                    m_sample_parameters["theta_jet"] = input.get<real>("bns_cocoon.theta_jet", real(0.2094));
                    m_sample_parameters["Gamma_cocoon"] = input.get<real>("bns_cocoon.Gamma_cocoon", real(2.0));
                    m_sample_parameters["jet_injection"] = input.get<bool>("bns_cocoon.jet_injection", false);
                    m_sample_parameters["n_radial"] = input.get<int>("bns_cocoon.n_radial", 200);
                    m_sample_parameters["n_angular"] = input.get<int>("bns_cocoon.n_angular", 160);
                }
                // Check for NS merger test type
                else if(test_type == "ns_merger_2d" || name_str.find("ns_merger") != std::string::npos) {
                    m_sample = Sample::NSMerger2D;
                    m_sample_parameters["R_star"] = input.get<real>("ns_merger.star1.radius", real(1.2));
                    m_sample_parameters["rho_c"] = input.get<real>("ns_merger.star1.central_density", real(2.8));
                    m_sample_parameters["separation"] = input.get<real>("ns_merger.separation", real(6.0));
                    m_sample_parameters["v_collision"] = input.get<real>("ns_merger.star1.velocity_x", real(0.15));
                    m_sample_parameters["n_radial"] = input.get<int>("ns_merger.star1.n_particles_radial", 30);
                }
                // Check for SR tangent velocity test (must come before sr_sod to match more specific)
                else if(name_str.find("tangent") != std::string::npos ||
                   name_str.find("sr_tangent") != std::string::npos) {
                    m_sample = Sample::SRTangentVelocity;
                    m_sample_parameters["N"] = input.get<int>("N", 1600);
                    m_sample_parameters["vt_left"] = input.get<real>("vt_left", real(0.9));
                    m_sample_parameters["vt_right"] = input.get<real>("vt_right", real(0.9));
                    m_sample_parameters["useGhostParticles"] = input.get<bool>("useGhostParticles", true);
                    m_sample_parameters["ghostLayers"] = input.get<int>("ghostLayers", 6);
                    // Non-uniform resolution support
                    if (auto opt = input.get_child_optional("N_left")) {
                        m_sample_parameters["N_left"] = input.get<int>("N_left");
                    }
                    if (auto opt = input.get_child_optional("N_right")) {
                        m_sample_parameters["N_right"] = input.get<int>("N_right");
                    }
                    if (auto opt = input.get_child_optional("leftResolutionRatio")) {
                        m_sample_parameters["leftResolutionRatio"] = input.get<real>("leftResolutionRatio");
                    }
                }
                // Check for SR-specific test names in the path
                else if(name_str.find("sr_sod") != std::string::npos || 
                   name_str.find("sod") != std::string::npos ||
                   name_str.find("ultra") != std::string::npos ||
                   name_str.find("blast") != std::string::npos) {
                    m_sample = Sample::SRSod;
                    m_sample_parameters["N"] = input.get<int>("N", 50);
                    m_sample_parameters["different_nu"] = input.get<bool>("different_nu", false);
                    m_sample_parameters["testType"] = input.get<std::string>("testType", "sod");
                    m_sample_parameters["particleMode"] = input.get<std::string>("particleMode", "equal_N");
                    // Ultra-relativistic test: left velocity (default 0.9 = 0.9c)
                    m_sample_parameters["v_left"] = input.get<real>("v_left", real(0.0));
                } else {
                    m_sample = Sample::DoNotUse;
                }
            } else {
                m_sample = Sample::DoNotUse;
            }
        }
    }

    m_output_dir = input.get<std::string>("outputDirectory");

    // Check if we're in resume mode - if so, some parameters can have defaults
    // since they'll be loaded from the checkpoint file
    bool is_resume_mode = input.get<bool>("checkpoint.autoResume", false) || 
                          !input.get<std::string>("checkpoint.resumeFile", "").empty();
    
    if (is_resume_mode) {
        std::cout << "Resume mode detected - physics parameters will be loaded from checkpoint" << std::endl;
    }

    // time
    m_param->time.start = input.get<real>("startTime", real(0));
    m_param->time.end   = input.get<real>("endTime", real(1.0));  // Default if resuming
    if(m_param->time.end < m_param->time.start) {
        THROW_ERROR("endTime < startTime");
    }
    m_param->time.output = input.get<real>("outputTime", (m_param->time.end - m_param->time.start) / 100);
    m_param->time.energy = input.get<real>("energyTime", m_param->time.output);

    // type - use defaults if resuming (will be overridden by checkpoint)
    std::string sph_type = input.get<std::string>("SPHType", "ssph");
    if(sph_type == "ssph") {
        m_param->type = SPHType::SSPH;
    } else if(sph_type == "disph") {
        m_param->type = SPHType::DISPH;
    } else if(sph_type == "gsph") {
        m_param->type = SPHType::GSPH;
    } else if(sph_type == "gdisph") {
        m_param->type = SPHType::GDISPH;
    } else if(sph_type == "srgsph") {
        m_param->type = SPHType::SRGSPH;
    } else if(sph_type == "grgsph") {
        m_param->type = SPHType::GRGSPH;
    } else if(sph_type == "gspmhd" || sph_type == "mhd") {
        m_param->type = SPHType::GSPMHD;
    } else if(sph_type == "srmhd") {
        m_param->type = SPHType::SRMHD;
    } else {
        THROW_ERROR("Unknown SPH type");
    }

    // CFL
    m_param->cfl.sound = input.get<real>("cflSound", 0.3);
    m_param->cfl.force = input.get<real>("cflForce", 0.125);

    // Artificial Viscosity
    m_param->av.alpha = input.get<real>("avAlpha", 1.0);
    m_param->av.use_balsara_switch = input.get<bool>("useBalsaraSwitch", true);
    m_param->av.use_time_dependent_av = input.get<bool>("useTimeDependentAV", false);
    if(m_param->av.use_time_dependent_av) {
        m_param->av.alpha_max = input.get<real>("alphaMax", 2.0);
        m_param->av.alpha_min = input.get<real>("alphaMin", 0.1);
        if(m_param->av.alpha_max < m_param->av.alpha_min) {
            THROW_ERROR("alphaMax < alphaMin");
        }
        m_param->av.epsilon = input.get<real>("epsilonAV", 0.2);
    }

    // Artificial Conductivity
    m_param->ac.is_valid = input.get<bool>("useArtificialConductivity", false);
    if(m_param->ac.is_valid) {
        m_param->ac.alpha = input.get<real>("alphaAC", 1.0);
    }

    // Tree
    m_param->tree.max_level = input.get<int>("maxTreeLevel", 20);
    m_param->tree.leaf_particle_num = input.get<int>("leafParticleNumber", 1);

    // Physics - use defaults if resuming (will be overridden by checkpoint)
    m_param->physics.neighbor_number = input.get<int>("neighborNumber", 50);
    m_param->physics.gamma = input.get<real>("gamma", 1.6666666666666667);  // Default for resume
    m_param->physics.c_smooth = input.get<real>("cSmooth", 1.0);  // Smoothing length expansion factor (default=1, typical=2)

    // Kernel
    const std::string kernel_default =
        (m_param->type == SPHType::SRGSPH || m_param->type == SPHType::GRGSPH) ? "gaussian" : "cubic_spline";
    std::string kernel_name = input.get<std::string>("kernel", kernel_default);
    if(kernel_name == "cubic_spline") {
        m_param->kernel = KernelType::CUBIC_SPLINE;
    } else if(kernel_name == "wendland") {
        m_param->kernel = KernelType::WENDLAND;
    } else if(kernel_name == "gaussian") {
        m_param->kernel = KernelType::GAUSSIAN;
    } else {
        THROW_ERROR("kernel is unknown.");
    }

    if((m_param->type == SPHType::SRGSPH || m_param->type == SPHType::GRGSPH) && m_param->kernel != KernelType::GAUSSIAN) {
        THROW_ERROR("SR/GR-GSPH currently supports only the Gaussian kernel. Set kernel=\"gaussian\".");
    }

    // smoothing length
    m_param->iterative_sml = input.get<bool>("iterativeSmoothingLength", true);

    // Preserve initial density (for shock tubes - skip density recalculation in initial_smoothing)
    m_param->preserve_initial_density = input.get<bool>("preserveInitialDensity", false);

    // periodic
    m_param->periodic.is_valid = input.get<bool>("periodic", false);

    // Per-dimension periodic flags (default: all dimensions periodic if periodic=true)
    // JSON format: "periodicDimensions": [false, true, true] for periodic in y,z only (3D shock tube)
    {
        auto periodic_dims_opt = input.get_child_optional("periodicDimensions");
        if(periodic_dims_opt) {
            auto & periodic_dims = *periodic_dims_opt;
            if(periodic_dims.size() != DIM) {
                THROW_ERROR("periodicDimensions array size != DIM");
            }
            int i = 0;
            for(auto & v : periodic_dims) {
                std::string val = v.second.data();
                m_param->periodic.per_dimension[i] = (val == "true" || val == "1");
                ++i;
            }
        } else {
            // Default: all dimensions periodic if periodic=true
            for(int i = 0; i < DIM; ++i) {
                m_param->periodic.per_dimension[i] = m_param->periodic.is_valid;
            }
        }
    }

    // Always read domain boundaries (rangeMin/rangeMax) - used for ghost particles even without periodic BC
    {
        auto range_max_opt = input.get_child_optional("rangeMax");
        if(range_max_opt) {
            auto & range_max = *range_max_opt;
            if(range_max.size() != DIM) {
                THROW_ERROR("rangeMax != DIM");
            }
            int i = 0;
            for(auto & v : range_max) {
                m_param->periodic.range_max[i] = std::stod(v.second.data());
                ++i;
            }
        }
    }
    {
        auto range_min_opt = input.get_child_optional("rangeMin");
        if(range_min_opt) {
            auto & range_min = *range_min_opt;
            if(range_min.size() != DIM) {
                THROW_ERROR("rangeMin != DIM");
            }
            int i = 0;
            for(auto & v : range_min) {
                m_param->periodic.range_min[i] = std::stod(v.second.data());
                ++i;
            }
        }
    }
    
    // Boundary type for ghost particles (when periodic is false)
    // Options: "reflecting" (wall), "inflow" (fixed boundary conditions)
    std::string boundary_type_str = input.get<std::string>("boundaryType", "inflow");
    if(boundary_type_str == "reflecting" || boundary_type_str == "wall") {
        m_param->periodic.boundary_type = BoundaryType::REFLECTING;
    } else {
        m_param->periodic.boundary_type = BoundaryType::INFLOW;  // Default: inflow (fixed BC)
    }

    // gravity
    m_param->gravity.is_valid = input.get<bool>("useGravity", false);
    // Always read G for use in initial conditions (e.g., Lane-Emden)
    m_param->gravity.constant = input.get<real>("G", 1.0);
    if(m_param->gravity.is_valid) {
        m_param->gravity.theta = input.get<real>("theta", 0.5);
        // Fixed vs h-dependent gravity softening
        m_param->gravity.use_fixed_softening = input.get<bool>("useFixedGravitySoftening", false);
        m_param->gravity.fixed_softening = input.get<real>("gravitySoftening", 0.1);
        
        // Gravity softening kernel type: "hernquist_katz" (default) or "wendland_c4"
        std::string softening_type_str = input.get<std::string>("gravitySofteningType", "hernquist_katz");
        if (softening_type_str == "wendland_c4" || softening_type_str == "wendland") {
            m_param->gravity.softening_type = GravitySofteningType::WENDLAND_C4;
        } else {
            m_param->gravity.softening_type = GravitySofteningType::HERNQUIST_KATZ;
        }
        
        // 1D kernel-convolved gravity (default: true for consistency with SPH pressure)
        m_param->gravity.use_kernel_gravity_1d = input.get<bool>("useKernelGravity1D", true);
        
        // 2D kernel-convolved gravity (default: true for consistency with SPH pressure)
        m_param->gravity.use_kernel_gravity_2d = input.get<bool>("useKernelGravity2D", true);
        
        // 2D planar slab gravity: 1D gravity in y-direction (default: false)
        m_param->gravity.use_kernel_gravity_planar_2d = input.get<bool>("useKernelGravityPlanar2D", false);
        
        // 3D cylindrical gravity: 2D radial gravity in xy-plane (default: false)
        m_param->gravity.use_kernel_gravity_cylinder_3d = input.get<bool>("useKernelGravityCylinder3D", false);
    }

    // Grad-h correction: applies to ALL SPH types (SSPH, GSPH, GDISPH, DISPH)
    // Enabled by default. Disabling causes core collapse in hydrostatic tests.
    m_param->gsph.use_gradh = input.get<bool>("useGradH", true);

    // GSPH and GDISPH specific settings
    if(m_param->type == SPHType::GSPH || m_param->type == SPHType::GDISPH) {
        m_param->gsph.is_2nd_order = input.get<bool>("use2ndOrderGSPH", true);

        // Riemann solver selection: "hll", "iterative", or "kitajima"
        std::string riemann_solver_str = input.get<std::string>("riemannSolver", "hll");
        if(riemann_solver_str == "iterative") {
            m_param->gsph.riemann_solver = RiemannSolverType::ITERATIVE;
        } else if(riemann_solver_str == "kitajima") {
            m_param->gsph.riemann_solver = RiemannSolverType::KITAJIMA;
        } else {
            m_param->gsph.riemann_solver = RiemannSolverType::HLL;
        }

        // Volume-based approach (Kitajima et al.) vs mass-based (standard GSPH)
        // When enabled, uses V_p = [Σ W]^(-1) for volume and V² in force formula
        // When disabled, uses ρ = Σ m W for density and 1/ρ² in force formula
        m_param->gsph.use_volume_based = input.get<bool>("useVolumeBased", false);
        m_param->gsph.eta = input.get<real>("etaSmoothingLength", 1.0);
        m_param->gsph.c_smooth = input.get<real>("cSmooth", 2.0);
    }

    // GSPMHD specific settings (Iwasaki & Inutsuka 2011)
    if(m_param->type == SPHType::GSPMHD) {
        m_param->mhd.is_2nd_order = input.get<bool>("use2ndOrderMHD", true);
        m_param->mhd.use_powell_correction = input.get<bool>("usePowellCorrection", true);
        m_param->mhd.c_h = input.get<real>("cH", 1.2);

        WRITE_LOG << "GSPMHD (Iwasaki & Inutsuka 2011) configuration:";
        WRITE_LOG << "  2nd order: " << (m_param->mhd.is_2nd_order ? "yes" : "no");
        WRITE_LOG << "  Powell div-B correction: " << (m_param->mhd.use_powell_correction ? "yes" : "no");
    }

    // SRMHD specific settings (SR-GSPH + GSPMHD)
    if(m_param->type == SPHType::SRMHD) {
        // SR parameters
        m_param->srgsph.is_2nd_order = input.get<bool>("use2ndOrderSRGSPH", false);
        m_param->srgsph.c_speed = input.get<real>("cSpeed", 1.0);
        m_param->srgsph.c_shock = input.get<real>("cShock", 3.0);
        m_param->srgsph.c_cd = input.get<real>("cContactDiscontinuity", 0.2);
        m_param->srgsph.eta = input.get<real>("etaSmoothingLength", 1.0);

        // Riemann solver type: "hll" (default) or "hllc"
        std::string riemann_type = input.get<std::string>("riemannSolver", "hll");
        if (riemann_type == "hllc" || riemann_type == "HLLC") {
            m_param->srgsph.riemann_solver = RiemannSolverType::HLLC;
        } else {
            m_param->srgsph.riemann_solver = RiemannSolverType::HLL;
        }

        // MHD parameters
        m_param->mhd.is_2nd_order = input.get<bool>("use2ndOrderMHD", false);
        m_param->mhd.use_powell_correction = input.get<bool>("usePowellCorrection", true);
        m_param->mhd.use_mhd = input.get<bool>("useMHD", true);  // default: full MHD

        WRITE_LOG << "SRMHD (SR-GSPH + GSPMHD) configuration:";
        WRITE_LOG << "  c_speed: " << m_param->srgsph.c_speed;
        WRITE_LOG << "  Riemann solver: " << (m_param->srgsph.riemann_solver == RiemannSolverType::HLLC ? "HLLC" : "HLL");
        WRITE_LOG << "  MHD enabled: " << (m_param->mhd.use_mhd ? "yes" : "no (hydro-only, uses SR-GSPH exact solver)");
        WRITE_LOG << "  Powell div-B correction: " << (m_param->mhd.use_powell_correction ? "yes" : "no");
    }

    // SRGSPH/GRGSPH (Fixed-h formulation, §2.2)
    if(m_param->type == SPHType::SRGSPH || m_param->type == SPHType::GRGSPH) {
        m_param->srgsph.is_2nd_order = input.get<bool>("use2ndOrderSRGSPH", true);
        m_param->srgsph.c_speed = input.get<real>("cSpeed", 1.0);
        m_param->srgsph.c_shock = input.get<real>("cShock", 3.0);
        m_param->srgsph.c_cd = input.get<real>("cContactDiscontinuity", 0.2);
        m_param->srgsph.eta = input.get<real>("etaSmoothingLength", 1.0);
        
        // Riemann solver type: "exact" (default) or "hllc"
        std::string riemann_type = input.get<std::string>("riemannSolver", "exact");
        if (riemann_type == "hllc" || riemann_type == "HLLC") {
            m_param->srgsph.riemann_solver = RiemannSolverType::HLLC;
            WRITE_LOG << "Using HLLC Riemann solver for SR/GR-GSPH";
        } else {
            m_param->srgsph.riemann_solver = RiemannSolverType::EXACT;
            WRITE_LOG << "Using Exact Riemann solver for SR/GR-GSPH (with HLLC fallback)";
        }
        
        // Fixed smoothing length (§2.2): auto-compute from eta if not specified
        // h is constant for all particles and all times
        m_param->srgsph.smoothing_length = input.get<real>("fixedSmoothingLength", -1.0);
    }

    // GRGSPH-specific parameters (metric configuration)
    if(m_param->type == SPHType::GRGSPH) {
        m_param->grgsph.metric_type = input.get<std::string>("metricType", "minkowski");
        m_param->grgsph.bh_mass = input.get<real>("blackHoleMass", 1.0);
        m_param->grgsph.bh_spin = input.get<real>("blackHoleSpin", 0.0);
        m_param->grgsph.geodesic_mode = input.get<bool>("geodesicMode", false);

        WRITE_LOG << "GR-GSPH metric configuration:";
        WRITE_LOG << "  Metric type: " << m_param->grgsph.metric_type;
        if(m_param->grgsph.metric_type == "schwarzschild" || m_param->grgsph.metric_type == "kerr") {
            WRITE_LOG << "  Black hole mass M: " << m_param->grgsph.bh_mass;
        }
        if(m_param->grgsph.metric_type == "kerr") {
            WRITE_LOG << "  Black hole spin a: " << m_param->grgsph.bh_spin;
        }
        if(m_param->grgsph.geodesic_mode) {
            WRITE_LOG << "  Geodesic mode: enabled (no SPH pressure forces)";
        }
    }

    // Relaxation (for Lane-Emden)
    
    // Thermal (ISM Cooling/Heating) - Koyama & Inutsuka (2000)
    m_param->thermal.enable_cooling = input.get<bool>("enableCooling", false);
    if (m_param->thermal.enable_cooling) {
        m_param->thermal.N_H_column = input.get<real>("columnDensity", 1.0e19);
        m_param->thermal.relaxation_time = input.get<real>("thermalRelaxationTime", 0.1);

        // Unit conversion factors for cooling physics
        // Default values for Galactic unit system:
        //   Length = 1 pc, Mass = 1000 M_sun, Velocity = 1 km/s
        //   Time = L/V = 0.978 Myr = 3.086e13 s
        //   Density unit = 1000 M_sun / pc^3 = 6.77e-20 g/cm^3
        //   n_H = rho / m_n where m_n = 1.27 * m_p = 2.12e-24 g
        //   => density_to_n_H = 6.77e-20 / 2.12e-24 = 31.9 cm^-3 per code density
        //   u_to_cgs = v^2 = (1 km/s)^2 = 1e10 erg/g
        //   t_to_cgs = 0.978 Myr = 3.086e13 s
        m_param->thermal.density_to_n_H = input.get<real>("densityToNumberDensity", 31.9);
        m_param->thermal.u_to_cgs = input.get<real>("energyToCGS", 1.0e10);
        m_param->thermal.t_to_cgs = input.get<real>("timeToCGS", 3.086e13);

        WRITE_LOG << "ISM Cooling enabled (Koyama & Inutsuka 2000)";
        WRITE_LOG << "* Column density N_H = " << m_param->thermal.N_H_column << " cm^-2";
        WRITE_LOG << "* Unit conversions:";
        WRITE_LOG << "    density_to_n_H = " << m_param->thermal.density_to_n_H << " cm^-3 per code density";
        WRITE_LOG << "    u_to_cgs = " << m_param->thermal.u_to_cgs << " erg/g per code energy";
        WRITE_LOG << "    t_to_cgs = " << m_param->thermal.t_to_cgs << " s per code time";
    }

    m_use_relaxation = input.get<bool>("useRelaxation", false);
    m_relaxation_steps = input.get<int>("relaxationSteps", 0);
    m_relaxation_output_freq = input.get<int>("relaxationOutputFreq", 10);
    m_relaxation_only = input.get<bool>("relaxationOnly", false);

    // GLASS pre-relaxation to uniformize particle spacing before main relaxation
    m_use_glass_relaxation = input.get<bool>("useGlassRelaxation", false);
    m_glass_relaxation_steps = input.get<int>("glassRelaxationSteps", 5000);
    m_glass_target_neighbors = input.get<int>("glassTargetNeighbors", 50);

    if(m_use_relaxation) {
        std::cout << "Relaxation enabled: " << m_relaxation_steps << " steps" << std::endl;
        std::cout << "Relaxation output frequency: every " << m_relaxation_output_freq << " steps" << std::endl;
        if(m_use_glass_relaxation) {
            std::cout << "GLASS pre-relaxation: " << m_glass_relaxation_steps << " steps (target " << m_glass_target_neighbors << " neighbors)" << std::endl;
        }
        if(m_relaxation_only) {
            std::cout << "Relaxation-only mode: Will exit after relaxation without running simulation" << std::endl;
        }
    }
    
    // Resume configuration
    m_checkpoint_file = input.get<std::string>("checkpoint.resumeFile", "");
    
    // Also support simple "resumeFromSnapshot" parameter
    if(m_checkpoint_file.empty()) {
        m_checkpoint_file = input.get<std::string>("resumeFromSnapshot", "");
    }
    
    m_resume_from_checkpoint = !m_checkpoint_file.empty();
    m_reset_time_on_resume = input.get<bool>("resetTimeOnResume", false);
    m_ghost_density_factor = input.get<real>("ghostDensityFactor", 1.0);
    m_spherical_boundary_radius = input.get<real>("sphericalBoundaryRadius", 0.0);
    
    if(m_resume_from_checkpoint) {
        std::cout << "Resume mode enabled:" << std::endl;
        std::cout << "  - Will resume from: " << m_checkpoint_file << std::endl;
        std::cout << "  - Reset time to 0: " << (m_reset_time_on_resume ? "yes" : "no (continue from snapshot time)") << std::endl;
        std::cout << "  - JSON config as SSOT: all physics parameters from config" << std::endl;
        if(m_ghost_density_factor != 1.0) {
            std::cout << "  - Ghost density factor: " << m_ghost_density_factor << std::endl;
        }
    }
    
    // Unit system configuration
    std::string unit_type_str = input.get<std::string>("units.type", "CODE");
    if(unit_type_str == "CODE") {
        m_units = UnitSystem(); // Default CODE units
    } else if(unit_type_str == "GALACTIC") {
        real length_kpc = input.get<real>("units.length_kpc", 1.0);
        real mass_msun = input.get<real>("units.mass_msun", 1.0e10);
        real velocity_kms = input.get<real>("units.velocity_kms", 1.0);
        m_units = UnitSystem::create_galactic(length_kpc, mass_msun, velocity_kms);
    } else if(unit_type_str == "SI") {
        real length_m = input.get<real>("units.length_m", 1.0);
        real mass_kg = input.get<real>("units.mass_kg", 1.0);
        real time_s = input.get<real>("units.time_s", 1.0);
        m_units = UnitSystem::create_si(length_m, mass_kg, time_s);
    } else if(unit_type_str == "CGS") {
        real length_cm = input.get<real>("units.length_cm", 1.0);
        real mass_g = input.get<real>("units.mass_g", 1.0);
        real time_s = input.get<real>("units.time_s", 1.0);
        m_units = UnitSystem::create_cgs(length_cm, mass_g, time_s);
    } else if(unit_type_str == "IMBH_ENCOUNTER") {
        // IMBH-cloud encounter units: pc, 10^3 M_sun, km/s
        real length_pc = input.get<real>("units.length_pc", 1.0);
        real mass_1e3Msun = input.get<real>("units.mass_1e3Msun", 1.0);
        real velocity_kms = input.get<real>("units.velocity_kms", 1.0);
        m_units = UnitSystem::create_imbh_encounter(length_pc, mass_1e3Msun, velocity_kms);
    } else if(unit_type_str == "RELATIVISTIC" || unit_type_str == "SR_TEST") {
        // Relativistic natural units with c=1
        std::string preset = input.get<std::string>("units.preset", "");
        if(preset == "neutron_star") {
            real length_km = input.get<real>("units.length_km", 10.0);
            real density = input.get<real>("units.density_g_cm3", 1.0e14);
            m_units = UnitSystem::create_neutron_star(length_km, density);
        } else if(preset == "relativistic_jet") {
            real length_pc = input.get<real>("units.length_pc", 1.0);
            real density = input.get<real>("units.density_g_cm3", UnitSystem::PROTON_MASS_G);
            m_units = UnitSystem::create_relativistic_jet(length_pc, density);
        } else {
            // Default: dimensionless SR test (c=1, all scales = 1)
            m_units = UnitSystem::create_sr_test();
        }
        
        // Custom labels can override defaults
        std::string length_label = input.get<std::string>("units.labels.length", "");
        std::string density_label = input.get<std::string>("units.labels.density", "");
        std::string time_label = input.get<std::string>("units.labels.time", "");
        std::string velocity_label = input.get<std::string>("units.labels.velocity", "");
        std::string pressure_label = input.get<std::string>("units.labels.pressure", "");
        
        if(!length_label.empty()) m_units.set_length_label(length_label);
        if(!density_label.empty()) m_units.set_density_label(density_label);
        if(!time_label.empty()) m_units.set_time_label(time_label);
        if(!velocity_label.empty()) m_units.set_velocity_label(velocity_label);
        if(!pressure_label.empty()) m_units.set_pressure_label(pressure_label);
    } else {
        std::cerr << "Warning: Unknown unit type '" << unit_type_str << "', using CODE units" << std::endl;
        m_units = UnitSystem();
    }
    
    // For SR/GR-GSPH, automatically use relativistic units if not explicitly specified
    if((m_param->type == SPHType::SRGSPH || m_param->type == SPHType::GRGSPH) && unit_type_str == "CODE") {
        std::cout << "Note: SR/GR-GSPH detected, auto-selecting dimensionless relativistic units (c=1)" << std::endl;
        m_units = UnitSystem::create_sr_test();
    }
    
    std::cout << "Unit system: " << m_units.get_type_name() << std::endl;
    if(m_units.is_relativistic()) {
        std::cout << "  Speed of light (code units): c = " << m_units.get_c_code() << std::endl;
        std::cout << "  Length unit: " << m_units.get_length_unit_name() << std::endl;
        std::cout << "  Time unit: " << m_units.get_time_unit_name() << std::endl;
        std::cout << "  Velocity unit: " << m_units.get_velocity_unit_name() << std::endl;
        std::cout << "  Density unit: " << m_units.get_density_unit_name() << std::endl;
    }
    
    // IMBH external force configuration
    m_use_external_bh = input.get<bool>("imbh_parameters.enabled", false);
    m_param->external_bh.enabled = false;  // Default to disabled
    if(m_use_external_bh) {
        std::cout << "\n=== IMBH External Force Configuration ===" << std::endl;
        
        // Check if values are already in code units (skip unit conversion)
        // Use this when particle data is in code units (e.g., resuming from snapshot)
        bool use_code_units = input.get<bool>("imbh_parameters.use_code_units", false);
        if(use_code_units) {
            std::cout << "  [use_code_units=true] BH parameters already in code units, skipping conversion" << std::endl;
        }
        
        // Read BH mass
        real M_BH_input = input.get<real>("imbh_parameters.M_BH", 1.0e5);
        real M_BH_code;
        if(use_code_units) {
            // M_BH is already in code mass units
            M_BH_code = M_BH_input;
        } else {
            // M_BH is in M_sun, convert to code units
            real M_BH_cgs = M_BH_input * UnitSystem::MSUN_TO_G;
            M_BH_code = m_units.from_physical_mass(M_BH_cgs);
        }
        
        // Read initial position array
        std::vector<real> bh_pos_input(DIM);
        bh_pos_input[0] = input.get<real>("imbh_parameters.BH_initial_position.0", 0.0);
        bh_pos_input[1] = input.get<real>("imbh_parameters.BH_initial_position.1", 0.0);
#if DIM == 3
        bh_pos_input[2] = input.get<real>("imbh_parameters.BH_initial_position.2", 0.0);
#endif
        // Convert position to code units
        std::vector<real> bh_pos(DIM);
        for(int d = 0; d < DIM; ++d) {
            if(use_code_units) {
                // Position already in code units
                bh_pos[d] = bh_pos_input[d];
            } else {
                // Position in pc, convert to code units
                real pos_cgs = bh_pos_input[d] * UnitSystem::PC_TO_CM;
                bh_pos[d] = m_units.from_physical_length(pos_cgs);
            }
        }
        
        // Read initial velocity array
        std::vector<real> bh_vel_input(DIM);
        bh_vel_input[0] = input.get<real>("imbh_parameters.BH_initial_velocity.0", 0.0);
        bh_vel_input[1] = input.get<real>("imbh_parameters.BH_initial_velocity.1", 0.0);
#if DIM == 3
        bh_vel_input[2] = input.get<real>("imbh_parameters.BH_initial_velocity.2", 0.0);
#endif
        // Convert velocity to code units
        std::vector<real> bh_vel(DIM);
        for(int d = 0; d < DIM; ++d) {
            if(use_code_units) {
                // Velocity already in code units
                bh_vel[d] = bh_vel_input[d];
            } else {
                // Velocity in km/s, convert to code units
                real vel_cgs = bh_vel_input[d] * UnitSystem::KM_TO_CM;
                bh_vel[d] = m_units.from_physical_velocity(vel_cgs);
            }
        }
        
        // Check if BH is moving
        bool is_moving = false;
        for(int d = 0; d < DIM; ++d) {
            if(std::abs(bh_vel[d]) > 1e-10) {
                is_moving = true;
                break;
            }
        }
        
        // Read softening epsilon
        // Default: 0.001 pc (~200 AU) or code units if use_code_units=true
        real softening_input = input.get<real>("imbh_parameters.softening_epsilon", 0.001);
        real softening_code;
        if(use_code_units) {
            softening_code = softening_input;
        } else {
            real softening_cgs = softening_input * UnitSystem::PC_TO_CM;
            softening_code = m_units.from_physical_length(softening_cgs);
        }

        // Read sink radius
        // Default: 0.005 pc (~1000 AU) or code units if use_code_units=true
        real sink_radius_input = input.get<real>("imbh_parameters.sink_radius", 0.005);
        real sink_radius_code;
        if(use_code_units) {
            sink_radius_code = sink_radius_input;
        } else {
            real sink_radius_cgs = sink_radius_input * UnitSystem::PC_TO_CM;
            sink_radius_code = m_units.from_physical_length(sink_radius_cgs);
        }

        // Read sink enable flag (default: true for fixed BH simulations)
        bool enable_sink = input.get<bool>("imbh_parameters.enable_sink", true);

        // Construct position and velocity vectors
        vec_t pos_vec, vel_vec;
#if DIM == 2
        pos_vec = vec_t(bh_pos[0], bh_pos[1]);
        vel_vec = vec_t(bh_vel[0], bh_vel[1]);
#elif DIM == 3
        pos_vec = vec_t(bh_pos[0], bh_pos[1], bh_pos[2]);
        vel_vec = vec_t(bh_vel[0], bh_vel[1], bh_vel[2]);
#endif

        // Create and initialize external BH module
        m_external_bh = std::make_shared<external_forces::PointMassBlackHole>();

        external_forces::PointMassBHParams bh_params;
        bh_params.mass = M_BH_code;
        bh_params.position = pos_vec;
        bh_params.velocity = vel_vec;
        bh_params.softening_length = softening_code;
        bh_params.sink_radius = sink_radius_code;
        bh_params.enable_sink = enable_sink;
        bh_params.G_constant = m_param->gravity.constant;
        bh_params.is_moving = is_moving;

        m_external_bh->initialize(bh_params);
        
        // Store external BH parameters in SPHParameters for energy calculations
        m_param->external_bh.enabled = true;
        m_param->external_bh.mass = M_BH_code;
        m_param->external_bh.softening = softening_code;
        m_param->external_bh.position = pos_vec;
        m_param->external_bh.G_constant = m_param->gravity.constant;
        
        // Print BH parameters with appropriate unit labels
        std::string mass_unit = use_code_units ? "code" : "M_☉";
        std::string len_unit = use_code_units ? "code" : "pc";
        std::string vel_unit = use_code_units ? "code" : "km/s";
        
        std::cout << "  BH mass: " << M_BH_input << " " << mass_unit << " (" << M_BH_code << " code units)" << std::endl;
        std::cout << "  Initial position: [" << bh_pos_input[0];
        for(size_t i = 1; i < bh_pos_input.size(); ++i) std::cout << ", " << bh_pos_input[i];
        std::cout << "] " << len_unit << " = [" << bh_pos[0];
        for(size_t i = 1; i < bh_pos.size(); ++i) std::cout << ", " << bh_pos[i];
        std::cout << "] code" << std::endl;
        std::cout << "  Initial velocity: [" << bh_vel_input[0];
        for(size_t i = 1; i < bh_vel_input.size(); ++i) std::cout << ", " << bh_vel_input[i];
        std::cout << "] " << vel_unit << " = [" << bh_vel[0];
        for(size_t i = 1; i < bh_vel.size(); ++i) std::cout << ", " << bh_vel[i];
        std::cout << "] code" << std::endl;
        std::cout << "  Softening epsilon: " << softening_input << " " << len_unit << " = " << softening_code << " code" << std::endl;
        std::cout << "  Sink radius: " << sink_radius_input << " " << len_unit << " = " << sink_radius_code << " code" << std::endl;
        std::cout << "  Sink accretion: " << (enable_sink ? "enabled" : "disabled") << std::endl;
        std::cout << "  Is moving: " << (is_moving ? "yes (WARNING: fixed BH recommended)" : "no (fixed at origin)") << std::endl;
        std::cout << "=========================================\n" << std::endl;
        
        // Read cloud initial conditions (for shifting particles after snapshot load)
        m_has_cloud_initial_conditions = false;
        try {
            // Try to read cloud_initial_position array by iterating (boost ptree uses empty keys for arrays)
            auto& cloud_pos_node = input.get_child("imbh_parameters.cloud_initial_position");
            
            std::cout << "\n=== Cloud Initial Conditions ===" << std::endl;
            
            // Read cloud initial position (in pc, convert to code units)
            real cloud_pos_pc[DIM] = {0.0};
            int idx = 0;
            for(auto& item : cloud_pos_node) {
                if(idx < DIM) {
                    cloud_pos_pc[idx] = item.second.get_value<real>();
                }
                idx++;
            }
            
            // Read cloud initial velocity array
            auto& cloud_vel_node = input.get_child("imbh_parameters.cloud_initial_velocity");
            real cloud_vel_kms[DIM] = {0.0};
            idx = 0;
            for(auto& item : cloud_vel_node) {
                if(idx < DIM) {
                    cloud_vel_kms[idx] = item.second.get_value<real>();
                }
                idx++;
            }
            
            m_has_cloud_initial_conditions = true;
            
            // Convert from pc to code units (or use directly if use_code_units is true)
            if(use_code_units) {
                // Direct: values already in code units (pc, km/s)
#if DIM == 2
                m_cloud_initial_position = vec_t(cloud_pos_pc[0], cloud_pos_pc[1]);
                m_cloud_initial_velocity = vec_t(cloud_vel_kms[0], cloud_vel_kms[1]);
#elif DIM == 3
                m_cloud_initial_position = vec_t(cloud_pos_pc[0], cloud_pos_pc[1], cloud_pos_pc[2]);
                m_cloud_initial_velocity = vec_t(cloud_vel_kms[0], cloud_vel_kms[1], cloud_vel_kms[2]);
#endif
                std::cout << "  Using code units directly (no conversion)" << std::endl;
            } else {
                // Convert from physical units (pc, km/s) to code units via CGS
#if DIM == 2
                m_cloud_initial_position = vec_t(
                    m_units.from_physical_length(cloud_pos_pc[0] * UnitSystem::PC_TO_CM),
                    m_units.from_physical_length(cloud_pos_pc[1] * UnitSystem::PC_TO_CM)
                );
#elif DIM == 3
                m_cloud_initial_position = vec_t(
                    m_units.from_physical_length(cloud_pos_pc[0] * UnitSystem::PC_TO_CM),
                    m_units.from_physical_length(cloud_pos_pc[1] * UnitSystem::PC_TO_CM),
                    m_units.from_physical_length(cloud_pos_pc[2] * UnitSystem::PC_TO_CM)
                );
#endif
            
                // Convert velocity from km/s to code units
#if DIM == 2
                m_cloud_initial_velocity = vec_t(
                    m_units.from_physical_velocity(cloud_vel_kms[0] * UnitSystem::KM_TO_CM),
                    m_units.from_physical_velocity(cloud_vel_kms[1] * UnitSystem::KM_TO_CM)
                );
#elif DIM == 3
                m_cloud_initial_velocity = vec_t(
                    m_units.from_physical_velocity(cloud_vel_kms[0] * UnitSystem::KM_TO_CM),
                    m_units.from_physical_velocity(cloud_vel_kms[1] * UnitSystem::KM_TO_CM),
                    m_units.from_physical_velocity(cloud_vel_kms[2] * UnitSystem::KM_TO_CM)
                );
#endif
            }
            
            std::cout << "  Target cloud position: [" << cloud_pos_pc[0];
            for(int d = 1; d < DIM; ++d) std::cout << ", " << cloud_pos_pc[d];
            std::cout << "] pc = [" << m_cloud_initial_position[0];
            for(int d = 1; d < DIM; ++d) std::cout << ", " << m_cloud_initial_position[d];
            std::cout << "] code" << std::endl;
            
            std::cout << "  Target cloud velocity: [" << cloud_vel_kms[0];
            for(int d = 1; d < DIM; ++d) std::cout << ", " << cloud_vel_kms[d];
            std::cout << "] km/s = [" << m_cloud_initial_velocity[0];
            for(int d = 1; d < DIM; ++d) std::cout << ", " << m_cloud_initial_velocity[d];
            std::cout << "] code" << std::endl;
            std::cout << "  Will shift cloud particles after snapshot load" << std::endl;
            std::cout << "==================================\n" << std::endl;
        } catch(const std::exception& e) {
            std::cout << "DEBUG: cloud_initial_position not found or parse error: " << e.what() << std::endl;
            m_has_cloud_initial_conditions = false;
        } catch(...) {
            std::cout << "DEBUG: cloud_initial_position unknown exception" << std::endl;
            m_has_cloud_initial_conditions = false;
        }
    } else {
        m_has_cloud_initial_conditions = false;
    }
    
    // Output configuration - convert boost ptree to nlohmann::json for the output section
    nlohmann::json output_json;
    
    // Try to read new array-of-objects format first
    try {
        nlohmann::json formats_array = nlohmann::json::array();
        for (const auto& format : input.get_child("output.formats")) {
            nlohmann::json format_obj;
            
            // If it's an object with properties
            if (format.second.size() > 0) {
                for (const auto& prop : format.second) {
                    const std::string& key = prop.first;
                    const auto& value = prop.second;
                    
                    // Try to determine type and convert appropriately
                    if (key == "type") {
                        format_obj[key] = value.get_value<std::string>();
                    } else if (key == "precision" || key == "compression") {
                        format_obj[key] = value.get_value<int>();
                    } else if (key == "binary") {
                        format_obj[key] = value.get_value<bool>();
                    } else {
                        // Default to string for unknown properties
                        format_obj[key] = value.get_value<std::string>();
                    }
                }
            } else {
                // If it's just a string, wrap it as {"type": "string"}
                format_obj["type"] = format.second.get_value<std::string>();
            }
            
            formats_array.push_back(format_obj);
        }
        output_json["formats"] = formats_array;
    } catch (const boost::property_tree::ptree_bad_path&) {
        // Formats array not found, fall back to legacy boolean flags
        output_json["enableCSV"] = input.get<bool>("output.enableCSV", true);
        output_json["enableHDF5"] = input.get<bool>("output.enableHDF5", false);
        output_json["enableVTK"] = input.get<bool>("output.enableVTK", false);
        
        // Legacy separate options
        output_json["csvPrecision"] = input.get<int>("output.csvPrecision", 16);
        output_json["hdf5Compression"] = input.get<int>("output.hdf5Compression", 6);
        output_json["vtkBinary"] = input.get<bool>("output.vtkBinary", true);
    }
    
    output_json["enableEnergyFile"] = input.get<bool>("output.enableEnergyFile", true);
    
    // Parse config using from_json
    OutputConfig output_config = OutputConfig::from_json(output_json);
    output_config.output_unit_type = m_units.get_type();
    
    // Display enabled formats
    std::cout << "Output configuration:" << std::endl;
    std::cout << "  - Formats: ";
    for (size_t i = 0; i < output_config.formats.size(); ++i) {
        std::cout << output_config.formats[i];
        if (i < output_config.formats.size() - 1) std::cout << ", ";
    }
    std::cout << std::endl;
    std::cout << "  - Energy file: " << (output_config.enable_energy_file ? "enabled" : "disabled") << std::endl;
    
    // Create OutputManager
    m_output_manager = std::make_shared<OutputManager>(output_config, m_units, m_output_dir);
    if(!m_output_manager->initialize()) {
        THROW_ERROR("Failed to initialize OutputManager");
    }
}

void Solver::run()
{
    initialize();
    
    // If relaxation-only mode, exit immediately after initialization
    if(m_relaxation_only) {
        std::cout << "\n=== Relaxation-Only Mode Complete ===" << std::endl;
        return;
    }
    
    assert(m_sim->get_particles().size() == m_sim->get_particle_num());

    const real t_end = m_param->time.end;
    const real t_start = m_sim->get_time();  // Current time (handles resume correctly)
    real t_out = t_start + m_param->time.output;
    real t_ene = t_start + m_param->time.energy;

    // For SR/GR-GSPH and SRMHD, compute initial N from kernel sum, then update ghosts
    // This ensures consistent initial conditions for force calculation
    if(m_param->type == SPHType::SRGSPH || m_param->type == SPHType::GRGSPH || m_param->type == SPHType::SRMHD) {
        m_sim->make_tree();
        m_pre->calculation(m_sim);  // Compute N for real particles
        update_ghost_particles();   // Mirror N to ghost particles

        // For conserved variable formulation, need to compute initial forces
        // so that dS and de are initialized (not garbage) for the first predict step
        m_fforce->calculation(m_sim);
    }

    // For GSPH with gravity, we need to compute initial grav_acc
    // so the gravity-aware Riemann solver has correct gravity information
    // on the very first timestep
    if(m_param->gravity.is_valid && 
       (m_param->type == SPHType::GSPH || m_param->type == SPHType::GDISPH)) {
        m_sim->make_tree();
        m_gforce->calculation(m_sim);  // Initialize grav_acc for all particles
    }

    // Write initial snapshot
    // DEBUG: Check density before writing
    {
        auto& particles = m_sim->get_particles();
        WRITE_LOG << "DEBUG: Before write snapshot_0000, particle[0].dens = " << particles[0].dens 
                 << ", pres = " << particles[0].pres << ", mass = " << particles[0].mass;
    }
    m_output_manager->write_snapshot(m_sim, m_param, m_snapshot_counter++);
    
    // Write initial energy
    real kinetic, thermal, potential;
    compute_total_energies(kinetic, thermal, potential);
    m_output_manager->write_energy(m_sim->get_time(), kinetic, thermal, potential);

    const auto start = std::chrono::system_clock::now();
    auto t_cout_i = start;
    int loop = 0;

    real t = m_sim->get_time();
    while(t < t_end) {
        integrate();
        const real dt = m_sim->get_dt();
        const int num = m_sim->get_particle_num();
        ++loop;

        m_sim->update_time();
        t = m_sim->get_time();
        
        // 1秒毎に画面出力
        const auto t_cout_f = std::chrono::system_clock::now();
        const real t_cout_s = std::chrono::duration_cast<std::chrono::seconds>(t_cout_f - t_cout_i).count();
        if(t_cout_s >= 1.0) {
            WRITE_LOG << "loop: " << loop << ", time: " << t << ", dt: " << dt << ", num: " << num;
            t_cout_i = std::chrono::system_clock::now();
        } else {
            WRITE_LOG_ONLY << "loop: " << loop << ", time: " << t << ", dt: " << dt << ", num: " << num;
        }

        if(t > t_out) {
            m_output_manager->write_snapshot(m_sim, m_param, m_snapshot_counter++);
            t_out += m_param->time.output;
        }

        if(t > t_ene) {
            compute_total_energies(kinetic, thermal, potential);
            m_output_manager->write_energy(t, kinetic, thermal, potential);
            t_ene += m_param->time.energy;
        }
    }
    
    const auto end = std::chrono::system_clock::now();
    const real calctime = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    WRITE_LOG << "\ncalculation is finished";
    WRITE_LOG << "calclation time: " << calctime << " ms";
}

void Solver::initialize()
{
    std::cout << "Solver::initialize() starting..." << std::endl;
    m_sim = std::make_shared<Simulation>(m_param);
    std::cout << "Simulation created" << std::endl;
    
    // Check if we should resume from snapshot
    OutputMetadata snapshot_data;  // Will contain physics parameters from snapshot (SSOT)
    bool resumed = false;
    
    if(m_resume_from_checkpoint && !m_checkpoint_file.empty()) {
        std::cout << "\n=== Attempting to resume from snapshot ===" << std::endl;
        if(m_output_manager->load_for_resume(m_checkpoint_file, m_sim, &snapshot_data)) {
            resumed = true;
            std::cout << "=== Successfully resumed from snapshot ===" << std::endl;
            
            // Reset time if requested
            if(m_reset_time_on_resume) {
                m_sim->set_time(0.0);
                std::cout << "✓ Simulation time reset to 0.0" << std::endl;
            } else {
                std::cout << "  Continuing from snapshot time: " << m_sim->get_time() << std::endl;
            }
            
            // Apply cloud initial position/velocity transform for IC files
            // If resetTimeOnResume is true, we're loading a fresh IC (relaxed cloud at origin)
            // and need to translate it to the correct position with the correct velocity.
            // If resetTimeOnResume is false, we're resuming an evolved simulation and
            // particles already have their correct positions.
            if(m_has_cloud_initial_conditions && m_reset_time_on_resume) {
                std::cout << "\n=== Applying Cloud Initial Conditions ===" << std::endl;
                std::cout << "  Translating particles from IC to cloud_initial_position" << std::endl;
                std::cout << "  Adding velocity boost from cloud_initial_velocity" << std::endl;
                
                auto& particles = m_sim->get_particles();
                const int num_p = m_sim->get_particle_num();
                
                for(int i = 0; i < num_p; ++i) {
                    // Apply position offset
                    for(int d = 0; d < DIM; ++d) {
                        particles[i].pos[d] += m_cloud_initial_position[d];
                    }
                    // Apply velocity boost
                    for(int d = 0; d < DIM; ++d) {
                        particles[i].vel[d] += m_cloud_initial_velocity[d];
                    }
                }
                
                std::cout << "  ✓ Applied transform to " << num_p << " particles" << std::endl;
                std::cout << "  Cloud center: [" << m_cloud_initial_position[0];
                for(int d = 1; d < DIM; ++d) std::cout << ", " << m_cloud_initial_position[d];
                std::cout << "] pc" << std::endl;
                std::cout << "  Cloud velocity: [" << m_cloud_initial_velocity[0];
                for(int d = 1; d < DIM; ++d) std::cout << ", " << m_cloud_initial_velocity[d];
                std::cout << "] km/s" << std::endl;
                std::cout << "==========================================\n" << std::endl;
            } else if(m_has_cloud_initial_conditions) {
                std::cout << "\n=== Cloud Initial Conditions (NOT applied) ===" << std::endl;
                std::cout << "  resetTimeOnResume=false: continuing evolved simulation" << std::endl;
                std::cout << "  Particle positions/velocities taken from snapshot as-is" << std::endl;
                std::cout << "================================================\n" << std::endl;
            }
            
            // JSON Config is SSOT: All parameters come from JSON, not snapshot metadata
            // Snapshot provides ONLY particle data: positions, velocities, masses, densities
            std::cout << "\nJSON Config as SSOT: Using ALL parameters from config file" << std::endl;
            std::cout << "Snapshot provides: particle positions, velocities, masses, densities" << std::endl;
            std::cout << "\nPhysics parameters from JSON config:" << std::endl;
            std::cout << "  - gamma: " << m_param->physics.gamma << std::endl;
            std::cout << "  - G: " << m_param->gravity.constant << std::endl;
            std::cout << "  - SPH type: " << (int)m_param->type << std::endl;
            std::cout << "  - Kernel: " << (int)m_param->kernel << std::endl;
            std::cout << "  - Neighbor number: " << m_param->physics.neighbor_number << std::endl;
            std::cout << "  - Alpha (AV): " << m_param->av.alpha << std::endl;
            std::cout << "  - Balsara switch: " << (m_param->av.use_balsara_switch ? "enabled" : "disabled") << std::endl;
            std::cout << "  - Time-dependent AV: " << (m_param->av.use_time_dependent_av ? "enabled" : "disabled") << std::endl;
            std::cout << "  - Gravity: " << (m_param->gravity.is_valid ? "enabled" : "disabled") << std::endl;
            std::cout << "  - End time: " << m_param->time.end << std::endl;
            std::cout << "  - Output interval: " << m_param->time.output << std::endl;
            
            // Restore Lane-Emden relaxation parameters if this is a relaxation snapshot
            if(snapshot_data.is_relaxation && m_sample == Sample::LaneEmden) {
                m_sample_parameters["alpha"] = snapshot_data.alpha_scaling;
                m_sample_parameters["rho_center"] = snapshot_data.rho_center;
                m_sample_parameters["K"] = snapshot_data.K;
                m_sample_parameters["R"] = snapshot_data.R;
                m_sample_parameters["M_total"] = snapshot_data.M_total;
                std::cout << "  - Restored Lane-Emden parameters from snapshot" << std::endl;
            }
            
            // Apply ghost density factor if specified
            // This strengthens/weakens the ghost envelope pressure confinement
            if(m_ghost_density_factor != 1.0) {
                std::cout << "\n=== Applying Ghost Density Factor ===" << std::endl;
                auto& particles = m_sim->get_particles();
                const int num_p = m_sim->get_particle_num();
                const real gamma = m_param->physics.gamma;
                int ghost_count = 0;
                real old_rho_sum = 0.0, new_rho_sum = 0.0;
                
                for(int i = 0; i < num_p; ++i) {
                    if(particles[i].is_ghost) {
                        old_rho_sum += particles[i].dens;
                        // Scale density
                        particles[i].dens *= m_ghost_density_factor;
                        new_rho_sum += particles[i].dens;
                        // Recalculate pressure for isothermal/polytropic EOS: P = ρ × u × (γ-1)
                        // Since u stays constant, P scales with ρ
                        particles[i].pres *= m_ghost_density_factor;
                        ghost_count++;
                    }
                }
                
                if(ghost_count > 0) {
                    std::cout << "  Modified " << ghost_count << " ghost particles" << std::endl;
                    std::cout << "  Avg density: " << (old_rho_sum/ghost_count) << " → " 
                              << (new_rho_sum/ghost_count) << " M☉/pc³" << std::endl;
                    std::cout << "  Pressure also scaled by " << m_ghost_density_factor << std::endl;
                    std::cout << "======================================\n" << std::endl;
                }
            }
            
            // For ALL resumes, find the highest existing snapshot number and continue from there
            // This handles cases where output frequency changed between runs and avoids overwriting
            int max_existing = 0;
            for(int i = 0; i < 10000; ++i) {
                std::ostringstream ss;
                ss << m_output_dir << "/snapshot_" << std::setfill('0') << std::setw(4) << i << ".csv";
                std::ifstream test(ss.str());
                if(test.good()) {
                    max_existing = i;
                } else {
                    break;  // Stop at first missing snapshot
                }
            }
            m_snapshot_counter = max_existing + 1;
            std::cout << "  - Snapshot counter set to: " << m_snapshot_counter << " (next: snapshot_" 
                      << std::setfill('0') << std::setw(4) << m_snapshot_counter << ".csv)" << std::endl;
        } else {
            std::cerr << "Warning: Failed to load snapshot, starting from scratch" << std::endl;
        }
    }
    
    // If not resumed, create initial conditions
    if(!resumed) {
        make_initial_condition();
        std::cout << "Initial condition made, particle_num = " << m_sim->get_particle_num() << std::endl;

        // Create ghost particles for shock tube samples (non-periodic boundary)
        if(!m_param->periodic.is_valid &&
           (m_sample == Sample::ShockTube || m_sample == Sample::ShockTube2D || m_sample == Sample::ShockTube3D)) {
            std::cout << "Creating ghost particles for shock tube boundary..." << std::endl;

            // Get domain parameters from sample_parameters
            const int Nx = boost::any_cast<int>(m_sample_parameters["Nx"]);
            const real Lx = 1.0;  // Shock tube domain is always [0, 1] in x
            const real dx = Lx / Nx;

#if DIM == 1
            const real Ly = 1.0;  // Not used in 1D
            const real Lz = 1.0;  // Not used in 1D
            const real dy = 1.0;
            const real dz = 1.0;
#elif DIM == 2
            const int Ny = boost::any_cast<int>(m_sample_parameters["Ny"]);
            const real Ly = boost::any_cast<real>(m_sample_parameters["Ly"]);
            const real Lz = 1.0;  // Not used in 2D
            const real dy = Ly / Ny;
            const real dz = 1.0;
#else  // DIM == 3
            const int Ny = boost::any_cast<int>(m_sample_parameters["Ny"]);
            const int Nz = boost::any_cast<int>(m_sample_parameters["Nz"]);
            const real Ly = boost::any_cast<real>(m_sample_parameters["Ly"]);
            const real Lz = boost::any_cast<real>(m_sample_parameters["Lz"]);
            const real dy = Ly / Ny;
            const real dz = Lz / Nz;
#endif

            const real gamma = m_param->physics.gamma;
            // Kernel support for Wendland C4 is ~2h, where h ≈ 2.3*dx for N_neighbor=50
            // So we need ~5 ghost layers to fully cover kernel support
            const int n_ghost_layers = 5;

            auto& particles = m_sim->get_particles();
            create_shock_tube_ghost_particles(particles, n_ghost_layers, dx, dy, dz, Ly, Lz, gamma);

            // Update particle count to include ghosts
            m_sim->set_particle_num(particles.size());
            std::cout << "Ghost particles created, total particle_num = " << particles.size() << std::endl;
        }
    }

    m_timestep = std::make_shared<TimeStep>();
    if(m_param->type == SPHType::SSPH) {
        m_pre = std::make_shared<PreInteraction>();
        m_fforce = std::make_shared<FluidForce>();
    } else if(m_param->type == SPHType::DISPH) {
        m_pre = std::make_shared<disph::PreInteraction>();
        m_fforce = std::make_shared<disph::FluidForce>();
    } else if(m_param->type == SPHType::GSPH) {
        m_pre = std::make_shared<gsph::PreInteraction>();
        m_fforce = std::make_shared<gsph::FluidForce>();
    } else if(m_param->type == SPHType::GDISPH) {
        m_pre = std::make_shared<gdisph::PreInteraction>();
        m_fforce = std::make_shared<gdisph::FluidForce>();
    } else if(m_param->type == SPHType::SRGSPH) {
        m_timestep = std::make_shared<srgsph::TimeStep>();  // Use SR timestep
        m_pre = std::make_shared<srgsph::PreInteraction>();
        m_fforce = std::make_shared<srgsph::FluidForce>();
    } else if(m_param->type == SPHType::GRGSPH) {
        m_timestep = std::make_shared<srgsph::TimeStep>();  // Use SR timestep (same for GR)

        // Create GR pre-interaction with metric-aware Lorentz factor
        auto gr_pre = std::make_shared<grgsph::GRPreInteraction>();
        auto gr_fforce = std::make_shared<grgsph::GRFluidForce>();

        // Set up the metric based on configuration
        // Create two metrics: one for pre-interaction, one for force
        auto metric_pre = grgsph::create_metric(
            m_param->grgsph.metric_type,
            m_param->grgsph.bh_mass,
            m_param->grgsph.bh_spin
        );
        auto metric_force = grgsph::create_metric(
            m_param->grgsph.metric_type,
            m_param->grgsph.bh_mass,
            m_param->grgsph.bh_spin
        );

        gr_pre->set_metric(std::move(metric_pre));
        gr_fforce->set_metric(std::move(metric_force));

        m_pre = gr_pre;
        m_fforce = gr_fforce;
    } else if(m_param->type == SPHType::GSPMHD) {
        m_pre = std::make_shared<gspmhd::PreInteraction>();
        m_fforce = std::make_shared<gspmhd::FluidForce>();
    } else if(m_param->type == SPHType::SRMHD) {
        m_timestep = std::make_shared<srmhd::TimeStep>();
        // For hydro-only mode, use SR-GSPH modules (known to work)
        if(!m_param->mhd.use_mhd) {
            std::cout << "[SRMHD] Hydro-only mode: using SR-GSPH pre-interaction and force" << std::endl;
            m_pre = std::make_shared<srgsph::PreInteraction>();
            m_fforce = std::make_shared<srgsph::FluidForce>();
        } else {
            m_pre = std::make_shared<srmhd::PreInteraction>();
            m_fforce = std::make_shared<srmhd::FluidForce>();
        }
    }
    m_gforce = std::make_shared<GravityForce>();

    // GSPH, GDISPH, SRGSPH, and GRGSPH require gradient arrays for MUSCL
    if(m_param->type == SPHType::GSPH || m_param->type == SPHType::GDISPH ||
       m_param->type == SPHType::SRGSPH || m_param->type == SPHType::GRGSPH) {
        std::vector<std::string> names;
        names.push_back("grad_density");
        names.push_back("grad_pressure");
        names.push_back("grad_velocity_0");
#if DIM == 2
        names.push_back("grad_velocity_1");
#elif DIM == 3
        names.push_back("grad_velocity_1");
        names.push_back("grad_velocity_2");
#endif
        // SRGSPH/GRGSPH also needs number density gradient
        if(m_param->type == SPHType::SRGSPH || m_param->type == SPHType::GRGSPH) {
            names.push_back("grad_number_density");
        }
        m_sim->add_vector_array(names);
    }

    m_timestep->initialize(m_param);
    m_pre->initialize(m_param);
    m_fforce->initialize(m_param);
    m_gforce->initialize(m_param);
    
    // Initialize relaxation for Lane-Emden if enabled
    if(m_use_relaxation && m_sample == Sample::LaneEmden) {
        std::cout << "\n=== Initializing Lane-Emden Relaxation ===" << std::endl;
        m_lane_emden_relax = std::make_shared<LaneEmdenRelaxation>();
        
        // Get parameters from sample_parameters (set in make_lane_emden)
        LaneEmdenRelaxationParams relax_params;
        relax_params.alpha_scaling = boost::any_cast<real>(m_sample_parameters["alpha"]);
        relax_params.rho_center = boost::any_cast<real>(m_sample_parameters["rho_center"]);
        relax_params.K = boost::any_cast<real>(m_sample_parameters["K"]);
        relax_params.R = boost::any_cast<real>(m_sample_parameters["R"]);
        relax_params.M_total = boost::any_cast<real>(m_sample_parameters["M_total"]);
        relax_params.G = m_param->gravity.constant;
        relax_params.gamma = m_param->physics.gamma;
        
        m_lane_emden_relax->initialize(relax_params);
        std::cout << "=== Relaxation Initialized ===" << std::endl;
        
        // Initialize sound speed and tree for relaxation calculations
        {
            auto& particles = m_sim->get_particles();
            const int num_p = m_sim->get_particle_num();
            const real gamma = m_param->physics.gamma;
            const real c_sound_factor = gamma * (gamma - 1.0);
            for(int i = 0; i < num_p; ++i) {
                particles[i].sound = std::sqrt(c_sound_factor * particles[i].ene);
            }
            
            auto tree = m_sim->get_tree();
            tree->resize(num_p);
            tree->make(particles, num_p);
        }
        
        // Run relaxation phase
        std::cout << "\n=== Starting Relaxation Phase (" << m_relaxation_steps << " steps) ===" << std::endl;
        
        // Check if resuming from checkpoint
        int start_step = 0;
        real accumulated_time = 0.0;
        
        if(resumed && snapshot_data.is_relaxation) {
            start_step = snapshot_data.relaxation_step;
            accumulated_time = snapshot_data.accumulated_time;
            std::cout << "Resuming from step " << start_step << " (time=" << accumulated_time << ")" << std::endl;
        }
        
        // Calculate target step: for fresh runs it's m_relaxation_steps,
        // for resumed runs it's start_step + m_relaxation_steps (additional steps)
        int target_step = start_step + m_relaxation_steps;
        std::cout << "Will run from step " << start_step << " to " << target_step 
                  << " (" << m_relaxation_steps << " steps)" << std::endl;
        
        // Prepare relaxation metadata for snapshot output
        OutputMetadata relax_meta;
        relax_meta.is_relaxation = true;
        relax_meta.relaxation_total_steps = target_step;  // Use target_step for total
        relax_meta.alpha_scaling = relax_params.alpha_scaling;
        relax_meta.rho_center = relax_params.rho_center;
        relax_meta.K = relax_params.K;
        relax_meta.R = relax_params.R;
        relax_meta.M_total = relax_params.M_total;
        
        std::cout << "Progress: [" << std::string(50, ' ') << "] 0%" << std::flush;
        
        int output_counter = 0;
        int last_percent = -1;
        int last_sub_percent = -1;
        
        // Timing for ETA calculation
        auto start_wall_time = std::chrono::steady_clock::now();
        double avg_step_time = 0.0;
        int timing_samples = 0;
        
        for(int step = start_step; step < target_step; ++step) {
            auto step_start = std::chrono::steady_clock::now();
            // Rebuild tree for accurate neighbor search
            // Particles move during relaxation, so tree must be updated
            auto& particles = m_sim->get_particles();
            const int num_p = m_sim->get_particle_num();
            
            auto tree = m_sim->get_tree();
            tree->resize(num_p);
            tree->make(particles, num_p);
            
            // Calculate SPH forces (pressure, gravity, etc.)
            m_pre->calculation(m_sim);
            m_fforce->calculation(m_sim);
            if(m_param->gravity.is_valid) {
                m_gforce->calculation(m_sim);
            }
            
            // Apply relaxation: subtract equilibrium forces to get net force
            // Net force = SPH forces - Lane-Emden analytical pressure gradient
            m_lane_emden_relax->apply_relaxation(m_sim, 0.0);  // damping_factor unused now
            
            // Calculate timestep dynamically based on CFL condition
            m_timestep->calculation(m_sim);
            real dt_relax = m_sim->get_dt();

            // Integrate positions with net acceleration
            // STEEPEST DESCENT: Zero velocities, move in direction of force
            // Use Δx = a·dt (NOT ½at² which is extremely slow)
            auto * periodic = m_sim->get_periodic().get();
            
#pragma omp parallel for
            for(int i = 0; i < num_p; ++i) {
                // Skip ghost particles
                if(particles[i].is_ghost) continue;
                
                // Zero velocities (quasi-static relaxation)
                particles[i].vel[0] = 0.0;
                particles[i].vel[1] = 0.0;
#if DIM == 3
                particles[i].vel[2] = 0.0;
#endif
                
                // Kinematic integration with v₀=0: Δx = ½at²
                const real half_dt2 = 0.5 * dt_relax * dt_relax;
                particles[i].pos[0] += particles[i].acc[0] * half_dt2;
                particles[i].pos[1] += particles[i].acc[1] * half_dt2;
#if DIM == 3
                particles[i].pos[2] += particles[i].acc[2] * half_dt2;
#endif
                
                periodic->apply(particles[i].pos);
            }

            // Remove particles that escaped beyond cloud boundary
            m_lane_emden_relax->remove_escaping_particles(m_sim, 1.1);

            // Update accumulated time
            accumulated_time += dt_relax;
            
            // Track step timing for ETA
            auto step_end = std::chrono::steady_clock::now();
            double step_duration = std::chrono::duration<double>(step_end - step_start).count();
            
            // Exponential moving average for step time (more weight on recent)
            if(timing_samples == 0) {
                avg_step_time = step_duration;
            } else {
                avg_step_time = 0.9 * avg_step_time + 0.1 * step_duration;
            }
            timing_samples++;
            
            // Update progress bar (calculate based on steps completed from start)
            int steps_completed = step - start_step + 1;  // +1 because step 0 is first completed step
            int percent = (steps_completed * 100) / m_relaxation_steps;
            
            // Calculate sub-progress within current output interval
            int steps_in_interval = step % m_relaxation_output_freq;
            int sub_percent = (steps_in_interval * 100) / m_relaxation_output_freq;
            
            // Update display every step for real-time feedback
            bool should_update = true;
            
            if(should_update) {
                int bar_width = 30;
                int filled = (steps_completed * bar_width) / m_relaxation_steps;
                
                // Calculate max acceleration for display
                real max_acc = 0.0;
                for(const auto& pi : particles) {
                    real a = std::sqrt(pi.acc[0]*pi.acc[0] + pi.acc[1]*pi.acc[1] + pi.acc[2]*pi.acc[2]);
                    max_acc = std::max(max_acc, a);
                }
                
                // Calculate ETA
                int steps_remaining = m_relaxation_steps - steps_completed;
                double eta_seconds = steps_remaining * avg_step_time;
                int eta_hours = static_cast<int>(eta_seconds / 3600);
                int eta_mins = static_cast<int>((eta_seconds - eta_hours * 3600) / 60);
                int eta_secs = static_cast<int>(eta_seconds) % 60;
                
                std::ostringstream eta_str;
                if(eta_hours > 0) {
                    eta_str << eta_hours << "h" << std::setw(2) << std::setfill('0') << eta_mins << "m";
                } else if(eta_mins > 0) {
                    eta_str << eta_mins << "m" << std::setw(2) << std::setfill('0') << eta_secs << "s";
                } else {
                    eta_str << eta_secs << "s";
                }
                
                // Use \33[2K to clear line, then \r to return to start
                std::cout << "\33[2K\r[" << std::string(filled, '=') << std::string(bar_width - filled, ' ') 
                          << "] " << std::setw(3) << std::setfill(' ') << percent << "% "
                          << step << "/" << target_step
                          << " ETA:" << eta_str.str()
                          << " a=" << std::fixed << std::setprecision(0) << max_acc
                          << std::flush;
                last_percent = percent;
                last_sub_percent = sub_percent;
            }
            
            // Output snapshots at specified frequency
            if(step % m_relaxation_output_freq == 0 || step == m_relaxation_steps - 1) {
                // Print newline to clear progress bar before snapshot message
                std::cout << std::endl;
                
                // Set simulation time to accumulated physical time
                m_sim->set_time(accumulated_time);
                
                // Prepare particles for output
                auto& p = m_sim->get_particles();
                const int num = m_sim->get_particle_num();
                const real gamma = m_param->physics.gamma;
                const real c_sound_factor = gamma * (gamma - 1.0);
                const real alpha = m_param->av.alpha;
                
#pragma omp parallel for
                for(int i = 0; i < num; ++i) {
                    p[i].alpha = alpha;
                    p[i].balsara = 1.0;
                    p[i].sound = std::sqrt(c_sound_factor * p[i].ene);
                }
                
                // Update relaxation metadata for current state
                relax_meta.relaxation_step = step;
                relax_meta.accumulated_time = accumulated_time;
                
                // Output snapshot with Lane-Emden metadata
                m_output_manager->write_snapshot(m_sim, m_param, m_snapshot_counter++, &relax_meta);
                output_counter++;
            }
        }
        
        // Final newline after progress bar
        std::cout << std::endl;
        std::cout << "=== Relaxation Complete ===" << std::endl;
        
        // If relaxation-only mode, output results and return early
        if(m_relaxation_only) {
            std::cout << "\n=== Relaxation-Only Mode: Outputting Results ===\n" << std::endl;
            
            // Reset time to 0 for output
            m_sim->set_time(0.0);
            
            // Calculate final state for output
            auto & p = m_sim->get_particles();
            const int num = m_sim->get_particle_num();
            const real gamma = m_param->physics.gamma;
            const real c_sound = gamma * (gamma - 1.0);
            
            const real alpha = m_param->av.alpha;
#pragma omp parallel for
            for(int i = 0; i < num; ++i) {
                p[i].alpha = alpha;
                p[i].balsara = 1.0;
                p[i].sound = std::sqrt(c_sound * p[i].ene);
            }
            
            auto tree = m_sim->get_tree();
            tree->resize(num);
            tree->make(p, num);
            
            // Calculate final forces and derivatives
            m_pre->calculation(m_sim);
            m_fforce->calculation(m_sim);
            if(m_param->gravity.is_valid) {
                m_gforce->calculation(m_sim);
            }
            
            // Update relaxation metadata for final state
            relax_meta.relaxation_step = m_relaxation_steps;
            relax_meta.accumulated_time = accumulated_time;
            
            // Output relaxed configuration with Lane-Emden metadata
            m_output_manager->write_snapshot(m_sim, m_param, m_snapshot_counter++, &relax_meta);
            
            // Output energy
            real kinetic, thermal, potential;
            compute_total_energies(kinetic, thermal, potential);
            m_output_manager->write_energy(0.0, kinetic, thermal, potential);
            
            std::cout << "=== Relaxed Configuration Saved ===" << std::endl;
            std::cout << "Check output directory for results" << std::endl;
            std::cout << "=== Exiting (No Simulation Run) ===\n" << std::endl;
            
            return;  // Exit initialize() early
        }
        
        std::cout << "=== Starting Main Simulation ===\n" << std::endl;
        
        // Reset time to 0 after relaxation
        m_sim->set_time(0.0);
    }

    // Initialize relaxation for PolytropicSlab if enabled
    if(m_use_relaxation && m_sample == Sample::PolytropicSlab) {
        std::cout << "\n=== Initializing Polytropic Slab Relaxation ===" << std::endl;
        m_polytropic_slab_relax = std::make_shared<PolytropicSlabRelaxation>();
        
        // Get base parameters
        const real rho_center = boost::any_cast<real>(m_sample_parameters["rho_center"]);
        const real K = boost::any_cast<real>(m_sample_parameters["K"]);
        const real G = m_param->gravity.constant;
        const real gamma = m_param->physics.gamma;
        
        // Derive other parameters
        // n = 1/(γ-1), polytropic index
        const real n = 1.0 / (gamma - 1.0);
        
        // α² = K(n+1)ρ_c^(1-n) / (2πG) for planar geometry
        const real alpha_sq = K * (n + 1.0) * std::pow(rho_center, 1.0 - n) / (2.0 * M_PI * G);
        const real alpha = std::sqrt(alpha_sq);
        
        // Need to solve Lane-Emden to get ξ₁, but let relaxation class do that
        // For now, estimate x_surface ≈ 2 (will be refined by the relaxation module)
        // Actually, the PolytropicSlabRelaxation will compute this internally
        
        PolytropicSlabRelaxationParams relax_params;
        relax_params.rho_center = rho_center;
        relax_params.K = K;
        relax_params.gamma = gamma;
        relax_params.n = n;
        relax_params.alpha_scaling = alpha;
        relax_params.x_surface = 0.0;  // Will be computed by relaxation module
        
        m_polytropic_slab_relax->initialize(relax_params);
        std::cout << "=== Relaxation Initialized ===" << std::endl;
        
        // Initialize sound speed and tree for relaxation calculations
        {
            auto& particles = m_sim->get_particles();
            const int num_p = m_sim->get_particle_num();
            const real gamma = m_param->physics.gamma;
            const real c_sound_factor = gamma * (gamma - 1.0);
            for(int i = 0; i < num_p; ++i) {
                particles[i].sound = std::sqrt(c_sound_factor * particles[i].ene);
            }
            
            auto tree = m_sim->get_tree();
            tree->resize(num_p);
            tree->make(particles, num_p);
        }
        
        // Run relaxation phase
        std::cout << "\n=== Starting Relaxation Phase (" << m_relaxation_steps << " steps) ===" << std::endl;
        
        int start_step = 0;
        real accumulated_time = 0.0;
        
        if(resumed && snapshot_data.is_relaxation) {
            start_step = snapshot_data.relaxation_step;
            accumulated_time = snapshot_data.accumulated_time;
            std::cout << "Resuming from step " << start_step << " (time=" << accumulated_time << ")" << std::endl;
        }
        
        int target_step = start_step + m_relaxation_steps;
        std::cout << "Will run from step " << start_step << " to " << target_step 
                  << " (" << m_relaxation_steps << " steps)" << std::endl;
        
        // Prepare relaxation metadata
        OutputMetadata relax_meta;
        relax_meta.is_relaxation = true;
        relax_meta.relaxation_total_steps = target_step;
        relax_meta.rho_center = relax_params.rho_center;
        relax_meta.K = relax_params.K;
        
        std::cout << "Progress: [" << std::string(50, ' ') << "] 0%" << std::flush;
        
        int output_counter = 0;
        int last_percent = -1;
        int last_sub_percent = -1;
        
        auto start_wall_time = std::chrono::steady_clock::now();
        double avg_step_time = 0.0;
        int timing_samples = 0;
        
        for(int step = start_step; step < target_step; ++step) {
            auto step_start = std::chrono::steady_clock::now();
            auto& particles = m_sim->get_particles();
            const int num_p = m_sim->get_particle_num();
            
            auto tree = m_sim->get_tree();
            tree->resize(num_p);
            tree->make(particles, num_p);
            
            // Calculate SPH forces
            m_pre->calculation(m_sim);
            m_fforce->calculation(m_sim);
            if(m_param->gravity.is_valid) {
                m_gforce->calculation(m_sim);
            }
            
            // Apply relaxation: subtract equilibrium forces
            m_polytropic_slab_relax->apply_relaxation(m_sim, 0.0);
            
            // Calculate timestep
            m_timestep->calculation(m_sim);
            real dt_relax = m_sim->get_dt();

            // Integrate positions with net acceleration (zero velocity constraint)
            auto * periodic = m_sim->get_periodic().get();
            real max_acc = 0.0;

#pragma omp parallel for reduction(max:max_acc)
            for(int i = 0; i < num_p; ++i) {
                // Zero velocities (constraint)
                particles[i].vel[0] = 0.0;
#if DIM >= 2
                particles[i].vel[1] = 0.0;
#endif
#if DIM == 3
                particles[i].vel[2] = 0.0;
#endif
                
                // Integrate position: Δx = ½at²
                particles[i].pos[0] += 0.5 * particles[i].acc[0] * dt_relax * dt_relax;
#if DIM >= 2
                particles[i].pos[1] += 0.5 * particles[i].acc[1] * dt_relax * dt_relax;
#endif
#if DIM == 3
                particles[i].pos[2] += 0.5 * particles[i].acc[2] * dt_relax * dt_relax;
#endif
                
                periodic->apply(particles[i].pos);
                
                real acc_mag = std::abs(particles[i].acc[0]);
#if DIM >= 2
                acc_mag = std::max(acc_mag, std::abs(particles[i].acc[1]));
#endif
#if DIM == 3
                acc_mag = std::max(acc_mag, std::abs(particles[i].acc[2]));
#endif
                max_acc = std::max(max_acc, acc_mag);
            }
            
            accumulated_time += dt_relax;
            
            // Progress bar update
            auto step_end = std::chrono::steady_clock::now();
            double step_time = std::chrono::duration<double>(step_end - step_start).count();
            timing_samples++;
            avg_step_time = ((timing_samples - 1) * avg_step_time + step_time) / timing_samples;
            
            int percent = (step - start_step + 1) * 100 / m_relaxation_steps;
            int sub_percent = ((step - start_step + 1) * 1000 / m_relaxation_steps) % 10;
            
            if(percent != last_percent || sub_percent != last_sub_percent) {
                int bar_width = 50;
                int filled = percent * bar_width / 100;
                
                int remaining_steps = target_step - step - 1;
                double eta_seconds = remaining_steps * avg_step_time;
                int eta_hours = static_cast<int>(eta_seconds / 3600);
                int eta_mins = static_cast<int>((eta_seconds - eta_hours * 3600) / 60);
                int eta_secs = static_cast<int>(eta_seconds - eta_hours * 3600 - eta_mins * 60);
                
                std::ostringstream eta_str;
                if(eta_hours > 0) {
                    eta_str << eta_hours << "h" << std::setw(2) << std::setfill('0') << eta_mins << "m";
                } else if(eta_mins > 0) {
                    eta_str << eta_mins << "m" << std::setw(2) << std::setfill('0') << eta_secs << "s";
                } else {
                    eta_str << eta_secs << "s";
                }
                
                std::cout << "\33[2K\r[" << std::string(filled, '=') << std::string(bar_width - filled, ' ') 
                          << "] " << std::setw(3) << std::setfill(' ') << percent << "% "
                          << step << "/" << target_step
                          << " ETA:" << eta_str.str()
                          << " a=" << std::fixed << std::setprecision(0) << max_acc
                          << std::flush;
                last_percent = percent;
                last_sub_percent = sub_percent;
            }
            
            // Output snapshots
            if(step % m_relaxation_output_freq == 0 || step == m_relaxation_steps - 1) {
                std::cout << std::endl;
                m_sim->set_time(accumulated_time);
                
                auto& p = m_sim->get_particles();
                const int num = m_sim->get_particle_num();
                const real gamma = m_param->physics.gamma;
                const real c_sound_factor = gamma * (gamma - 1.0);
                const real alpha = m_param->av.alpha;
                
#pragma omp parallel for
                for(int i = 0; i < num; ++i) {
                    p[i].alpha = alpha;
                    p[i].balsara = 1.0;
                    p[i].sound = std::sqrt(c_sound_factor * p[i].ene);
                }
                
                relax_meta.relaxation_step = step;
                relax_meta.accumulated_time = accumulated_time;
                
                m_output_manager->write_snapshot(m_sim, m_param, m_snapshot_counter++, &relax_meta);
                output_counter++;
            }
        }
        
        std::cout << std::endl;
        std::cout << "=== Relaxation Complete ===" << std::endl;
        
        // If relaxation-only mode, output and exit
        if(m_relaxation_only) {
            std::cout << "\n=== Relaxation-Only Mode: Outputting Results ===\n" << std::endl;
            
            m_sim->set_time(0.0);
            
            auto & p = m_sim->get_particles();
            const int num = m_sim->get_particle_num();
            const real gamma = m_param->physics.gamma;
            const real c_sound = gamma * (gamma - 1.0);
            const real alpha = m_param->av.alpha;
            
#pragma omp parallel for
            for(int i = 0; i < num; ++i) {
                p[i].alpha = alpha;
                p[i].balsara = 1.0;
                p[i].sound = std::sqrt(c_sound * p[i].ene);
            }
            
            auto tree = m_sim->get_tree();
            tree->resize(num);
            tree->make(p, num);
            
            m_pre->calculation(m_sim);
            m_fforce->calculation(m_sim);
            if(m_param->gravity.is_valid) {
                m_gforce->calculation(m_sim);
            }
            
            relax_meta.relaxation_step = m_relaxation_steps;
            relax_meta.accumulated_time = accumulated_time;
            
            m_output_manager->write_snapshot(m_sim, m_param, m_snapshot_counter++, &relax_meta);
            
            real kinetic, thermal, potential;
            compute_total_energies(kinetic, thermal, potential);
            m_output_manager->write_energy(0.0, kinetic, thermal, potential);
            
            std::cout << "=== Relaxed Configuration Saved ===" << std::endl;
            std::cout << "Check output directory for results" << std::endl;
            std::cout << "=== Exiting (No Simulation Run) ===\n" << std::endl;
            
            return;
        }
        
        std::cout << "=== Starting Main Simulation ===\n" << std::endl;
        m_sim->set_time(0.0);
    }

    // Initialize relaxation for PolytropicSlab2D if enabled
    if(m_use_relaxation && m_sample == Sample::PolytropicSlab2D) {
        std::cout << "\n=== Initializing 2D Polytropic Slab Relaxation ===" << std::endl;
        m_polytropic_slab_2d_relax = std::make_shared<PolytropicSlab2DRelaxation>();
        
        // Get base parameters
        const real rho_center = boost::any_cast<double>(m_sample_parameters["rho_center"]);
        const real K = boost::any_cast<double>(m_sample_parameters["K"]);
        const real L_x = boost::any_cast<double>(m_sample_parameters["L_x"]);
        const real G = m_param->gravity.constant;
        const real gamma = m_param->physics.gamma;
        
        // Derive other parameters
        // n = 1/(γ-1), polytropic index
        const real n = 1.0 / (gamma - 1.0);
        
        // α² = K(n+1)ρ_c^(1-n) / (2πG) for planar geometry
        const real alpha_sq = K * (n + 1.0) * std::pow(rho_center, 1.0 - n) / (2.0 * M_PI * G);
        const real alpha = std::sqrt(alpha_sq);
        
        PolytropicSlab2DRelaxationParams relax_params;
        relax_params.rho_center = rho_center;
        relax_params.K = K;
        relax_params.gamma = gamma;
        relax_params.n = n;
        relax_params.alpha_scaling = alpha;
        relax_params.y_surface = 0.0;  // Will be computed by relaxation module
        relax_params.L_x = L_x;
        
        m_polytropic_slab_2d_relax->initialize(relax_params);
        std::cout << "=== 2D Relaxation Initialized ===" << std::endl;
        
        // Initialize sound speed and tree for relaxation calculations
        {
            auto& particles = m_sim->get_particles();
            const int num_p = m_sim->get_particle_num();
            const real gamma = m_param->physics.gamma;
            const real c_sound_factor = gamma * (gamma - 1.0);
            for(int i = 0; i < num_p; ++i) {
                particles[i].sound = std::sqrt(c_sound_factor * particles[i].ene);
            }
            
            auto tree = m_sim->get_tree();
            tree->resize(num_p);
            tree->make(particles, num_p);
        }
        
        // Run relaxation phase
        std::cout << "\n=== Starting 2D Relaxation Phase (" << m_relaxation_steps << " steps) ===" << std::endl;
        
        int start_step = 0;
        real accumulated_time = 0.0;
        
        // Metadata for relaxation snapshots
        OutputMetadata relax_meta;
        relax_meta.is_relaxation = true;
        relax_meta.relaxation_total_steps = m_relaxation_steps;
        relax_meta.accumulated_time = 0.0;
        
        // Progress tracking
        int output_counter = 0;
        int last_percent = -1;
        int last_sub_percent = -1;
        int target_step = start_step + m_relaxation_steps;
        
        // Timing for ETA
        auto start_time = std::chrono::steady_clock::now();
        double avg_step_time = 0.0;
        int timing_samples = 0;
        
        for(int step = start_step; step < target_step; ++step) {
            auto step_start = std::chrono::steady_clock::now();
            
            // Update particle properties
            auto& particles = m_sim->get_particles();
            const int num_p = m_sim->get_particle_num();
            
            auto tree = m_sim->get_tree();
            tree->resize(num_p);
            tree->make(particles, num_p);
            
            // Calculate SPH forces
            m_pre->calculation(m_sim);
            m_fforce->calculation(m_sim);
            if(m_param->gravity.is_valid) {
                m_gforce->calculation(m_sim);
            }
            
            // Apply relaxation: velocity damping to let particles settle
            // We don't subtract analytical forces - instead let SPH find its own equilibrium
            m_polytropic_slab_2d_relax->apply_relaxation(m_sim, 0.0); 
            
            // Calculate timestep
            m_timestep->calculation(m_sim);
            real dt_relax = m_sim->get_dt();

            // Get polytropic constant K and gamma for isentropic EOS reset
            const real K_poly = m_polytropic_slab_2d_relax->get_K();
            const real gamma_poly = m_polytropic_slab_2d_relax->get_gamma();
            
            // Integrate positions with net acceleration (zero velocity constraint)
            auto * periodic = m_sim->get_periodic().get();
            real max_acc = 0.0;
            
#pragma omp parallel for reduction(max:max_acc)
            for(int i = 0; i < num_p; ++i) {
                // Zero velocities (constraint for quasi-static relaxation)
                particles[i].vel[0] = 0.0;
                particles[i].vel[1] = 0.0;
                
                // Integrate position: Δx = ½at² (no velocity term)
                particles[i].pos[0] += 0.5 * particles[i].acc[0] * dt_relax * dt_relax;
                particles[i].pos[1] += 0.5 * particles[i].acc[1] * dt_relax * dt_relax;
                
                periodic->apply(particles[i].pos);
                
                // Reset internal energy to isentropic value: P = Kρ^γ, u = P/((γ-1)ρ)
                // This maintains the isentropic EOS during relaxation
                const real rho = particles[i].dens;
                if(rho > 0) {
                    const real pres_isentropic = K_poly * std::pow(rho, gamma_poly);
                    particles[i].pres = pres_isentropic;
                    particles[i].ene = pres_isentropic / ((gamma_poly - 1.0) * rho);
                    particles[i].sound = std::sqrt(gamma_poly * pres_isentropic / rho);
                }
                
                real acc_mag = std::max(std::abs(particles[i].acc[0]), std::abs(particles[i].acc[1]));
                max_acc = std::max(max_acc, acc_mag);
            }
            
            accumulated_time += dt_relax;
            
            // Progress bar update
            auto step_end = std::chrono::steady_clock::now();
            double step_time = std::chrono::duration<double>(step_end - step_start).count();
            timing_samples++;
            avg_step_time = ((timing_samples - 1) * avg_step_time + step_time) / timing_samples;
            
            int percent = (step - start_step + 1) * 100 / m_relaxation_steps;
            int sub_percent = ((step - start_step + 1) * 1000 / m_relaxation_steps) % 10;
            
            if(percent != last_percent || sub_percent != last_sub_percent) {
                int bar_width = 50;
                int filled = percent * bar_width / 100;
                
                int remaining_steps = target_step - step - 1;
                double eta_seconds = remaining_steps * avg_step_time;
                int eta_mins = static_cast<int>(eta_seconds / 60);
                int eta_secs = static_cast<int>(eta_seconds - eta_mins * 60);
                
                std::cout << "\33[2K\r[" << std::string(filled, '=') << std::string(bar_width - filled, ' ') 
                          << "] " << std::setw(3) << std::setfill(' ') << percent << "% "
                          << step << "/" << target_step
                          << " ETA:" << eta_mins << "m" << eta_secs << "s"
                          << " a=" << std::fixed << std::setprecision(2) << max_acc
                          << std::flush;
                last_percent = percent;
                last_sub_percent = sub_percent;
            }
            
            // Output snapshots
            if(step % m_relaxation_output_freq == 0 || step == target_step - 1) {
                std::cout << std::endl;
                m_sim->set_time(accumulated_time);
                
                auto& p = m_sim->get_particles();
                const int num = m_sim->get_particle_num();
                const real gamma = m_param->physics.gamma;
                const real c_sound_factor = gamma * (gamma - 1.0);
                const real alpha = m_param->av.alpha;
                
#pragma omp parallel for
                for(int i = 0; i < num; ++i) {
                    p[i].alpha = alpha;
                    p[i].balsara = 1.0;
                    p[i].sound = std::sqrt(c_sound_factor * p[i].ene);
                }
                
                relax_meta.relaxation_step = step;
                relax_meta.accumulated_time = accumulated_time;
                
                m_output_manager->write_snapshot(m_sim, m_param, m_snapshot_counter++, &relax_meta);
                output_counter++;
            }
        }
        
        std::cout << std::endl;
        std::cout << "=== 2D Relaxation Complete ===" << std::endl;
        
        // If relaxation-only mode, output and exit
        if(m_relaxation_only) {
            std::cout << "\n=== Relaxation-Only Mode: Outputting Results ===\n" << std::endl;
            
            m_sim->set_time(0.0);
            
            auto & p = m_sim->get_particles();
            const int num = m_sim->get_particle_num();
            const real gamma = m_param->physics.gamma;
            const real c_sound = gamma * (gamma - 1.0);
            const real alpha = m_param->av.alpha;
            
#pragma omp parallel for
            for(int i = 0; i < num; ++i) {
                p[i].alpha = alpha;
                p[i].balsara = 1.0;
                p[i].sound = std::sqrt(c_sound * p[i].ene);
            }
            
            auto tree = m_sim->get_tree();
            tree->resize(num);
            tree->make(p, num);
            
            m_pre->calculation(m_sim);
            m_fforce->calculation(m_sim);
            if(m_param->gravity.is_valid) {
                m_gforce->calculation(m_sim);
            }
            
            relax_meta.relaxation_step = m_relaxation_steps;
            relax_meta.accumulated_time = accumulated_time;
            
            m_output_manager->write_snapshot(m_sim, m_param, m_snapshot_counter++, &relax_meta);
            
            real kinetic, thermal, potential;
            compute_total_energies(kinetic, thermal, potential);
            m_output_manager->write_energy(0.0, kinetic, thermal, potential);
            
            std::cout << "=== 2D Relaxed Configuration Saved ===" << std::endl;
            std::cout << "=== Exiting (No Simulation Run) ===\n" << std::endl;
            
            return;
        }
        
        std::cout << "=== Starting Main Simulation ===\n" << std::endl;
        m_sim->set_time(0.0);
    }

    // Initialize relaxation for BonnorEbertKI2000 if enabled
    if(m_use_relaxation && m_sample == Sample::BonnorEbertKI2000) {
        std::cout << "\n=== Initializing K&I 2000 Bonnor-Ebert Relaxation ===" << std::endl;
        m_ki2000_relax = std::make_shared<KoyamaInutsukaRelaxation>();

        // Try to load profile from file saved by IC generator (SSOT pattern)
        std::string profile_file = m_output_dir + "/ki2000_profile.dat";
        bool loaded_from_file = false;

        if (m_sample_parameters.count("ki2000_profile_file")) {
            profile_file = boost::any_cast<std::string>(m_sample_parameters["ki2000_profile_file"]);
        }

        // Try to load from file first (ensures same profile as IC)
        loaded_from_file = m_ki2000_relax->load_profile_from_file(profile_file);

        if (!loaded_from_file) {
            // Fallback: compute profile (may differ from IC!)
            std::cout << "WARNING: Could not load profile from " << profile_file << std::endl;
            std::cout << "WARNING: Computing profile independently - may differ from IC!" << std::endl;

            // Get parameters from sample_parameters
            const real R_cloud = boost::any_cast<double>(m_sample_parameters["R_cloud"]);
            const real M_cloud = boost::any_cast<double>(m_sample_parameters["M_cloud"]);
            const real P_ext_K_cm3 = boost::any_cast<double>(m_sample_parameters["P_ext_K_cm3"]);
            const real rho_center_code = boost::any_cast<double>(m_sample_parameters["rho_center_code"]);
            const real N_H_cm2 = boost::any_cast<double>(m_sample_parameters["N_H_cm2"]);
            const real density_to_n = boost::any_cast<double>(m_sample_parameters["density_to_n"]);
            const real G = m_param->gravity.constant;
            const real gamma = m_param->physics.gamma;

            KIRelaxationParams relax_params;
            relax_params.R_cloud = R_cloud;
            relax_params.M_cloud = M_cloud;
            relax_params.P_ext = P_ext_K_cm3;  // P/k_B in K cm^-3
            relax_params.rho_center = rho_center_code;
            relax_params.N_H = N_H_cm2;
            relax_params.G = G;
            relax_params.gamma = gamma;
            relax_params.density_to_n = density_to_n;
            relax_params.pressure_to_cgs = 1.0;  // Will be set properly

            m_ki2000_relax->initialize(relax_params);
        }
        std::cout << "=== K&I 2000 Relaxation Initialized ===" << std::endl;
        std::cout << "  Profile R_cloud = " << m_ki2000_relax->get_R_cloud() << " [code]" << std::endl;

        // Initialize sound speed and tree for relaxation calculations
        {
            auto& particles = m_sim->get_particles();
            const int num_p = m_sim->get_particle_num();
            const real gamma = m_param->physics.gamma;
            const real c_sound_factor = gamma * (gamma - 1.0);
            for(int i = 0; i < num_p; ++i) {
                particles[i].sound = std::sqrt(c_sound_factor * particles[i].ene);
            }

            auto tree = m_sim->get_tree();
            tree->resize(num_p);
            tree->make(particles, num_p);
        }

        // Run relaxation phase
        std::cout << "\n=== Starting K&I 2000 Relaxation Phase (" << m_relaxation_steps << " steps) ===" << std::endl;

        int start_step = 0;
        real accumulated_time = 0.0;

        // Metadata for relaxation snapshots
        OutputMetadata relax_meta;
        relax_meta.is_relaxation = true;
        relax_meta.relaxation_total_steps = m_relaxation_steps;
        relax_meta.accumulated_time = 0.0;

        // Progress tracking
        int output_counter = 0;
        int last_percent = -1;
        int last_sub_percent = -1;
        int target_step = start_step + m_relaxation_steps;

        // Timing for ETA
        auto start_time = std::chrono::steady_clock::now();
        double avg_step_time = 0.0;
        int timing_samples = 0;

        for(int step = start_step; step < target_step; ++step) {
            auto step_start = std::chrono::steady_clock::now();

            // Update particle properties
            auto& particles = m_sim->get_particles();
            const int num_p = m_sim->get_particle_num();

            auto tree = m_sim->get_tree();
            tree->resize(num_p);
            tree->make(particles, num_p);

            // Calculate SPH forces
            m_pre->calculation(m_sim);
            m_fforce->calculation(m_sim);
            // Note: gravity is OFF during relaxation

            // Apply K&I 2000 relaxation: subtract analytical pressure gradient
            m_ki2000_relax->apply_relaxation(m_sim, 0.0);

            // Calculate timestep
            m_timestep->calculation(m_sim);
            real dt_relax = m_sim->get_dt();

            // Integrate positions with net acceleration (zero velocity constraint)
            auto * periodic = m_sim->get_periodic().get();
            real max_acc = 0.0;

#pragma omp parallel for reduction(max:max_acc)
            for(int i = 0; i < num_p; ++i) {
                // Skip ghost/envelope particles - they provide fixed pressure boundary
                if(particles[i].is_ghost) {
                    particles[i].vel = 0.0;
                    particles[i].acc = 0.0;
                    continue;
                }

                // Zero velocities (constraint for quasi-static relaxation)
                particles[i].vel[0] = 0.0;
                particles[i].vel[1] = 0.0;
#if DIM == 3
                particles[i].vel[2] = 0.0;
#endif

                // Integrate position: Δx = ½at² (no velocity term)
                particles[i].pos[0] += 0.5 * particles[i].acc[0] * dt_relax * dt_relax;
                particles[i].pos[1] += 0.5 * particles[i].acc[1] * dt_relax * dt_relax;
#if DIM == 3
                particles[i].pos[2] += 0.5 * particles[i].acc[2] * dt_relax * dt_relax;
#endif

                if(periodic) periodic->apply(particles[i].pos);

                real acc_mag = std::abs(particles[i].acc[0]);
                acc_mag = std::max(acc_mag, std::abs(particles[i].acc[1]));
#if DIM == 3
                acc_mag = std::max(acc_mag, std::abs(particles[i].acc[2]));
#endif
                max_acc = std::max(max_acc, acc_mag);
            }

            accumulated_time += dt_relax;

            // Progress bar update
            auto step_end = std::chrono::steady_clock::now();
            double step_time = std::chrono::duration<double>(step_end - step_start).count();
            timing_samples++;
            avg_step_time = ((timing_samples - 1) * avg_step_time + step_time) / timing_samples;

            int percent = (step - start_step + 1) * 100 / m_relaxation_steps;
            int sub_percent = ((step - start_step + 1) * 1000 / m_relaxation_steps) % 10;

            if(percent != last_percent || sub_percent != last_sub_percent) {
                int bar_width = 50;
                int filled = percent * bar_width / 100;

                int remaining_steps = target_step - step - 1;
                double eta_seconds = remaining_steps * avg_step_time;
                int eta_mins = static_cast<int>(eta_seconds / 60);
                int eta_secs = static_cast<int>(eta_seconds - eta_mins * 60);

                std::cout << "\33[2K\r[" << std::string(filled, '=') << std::string(bar_width - filled, ' ')
                          << "] " << std::setw(3) << std::setfill(' ') << percent << "% "
                          << step << "/" << target_step
                          << " ETA:" << eta_mins << "m" << eta_secs << "s"
                          << " a=" << std::fixed << std::setprecision(2) << max_acc
                          << std::flush;
                last_percent = percent;
                last_sub_percent = sub_percent;
            }

            // Output snapshots
            if(step % m_relaxation_output_freq == 0 || step == target_step - 1) {
                std::cout << std::endl;
                m_sim->set_time(accumulated_time);

                auto& p = m_sim->get_particles();
                const int num = m_sim->get_particle_num();
                const real gamma = m_param->physics.gamma;
                const real c_sound_factor = gamma * (gamma - 1.0);
                const real alpha = m_param->av.alpha;

#pragma omp parallel for
                for(int i = 0; i < num; ++i) {
                    p[i].alpha = alpha;
                    p[i].balsara = 1.0;
                    p[i].sound = std::sqrt(c_sound_factor * p[i].ene);
                }

                relax_meta.relaxation_step = step;
                relax_meta.accumulated_time = accumulated_time;

                m_output_manager->write_snapshot(m_sim, m_param, m_snapshot_counter++, &relax_meta);
                output_counter++;
            }
        }

        std::cout << std::endl;
        std::cout << "=== K&I 2000 Relaxation Complete ===" << std::endl;

        // If relaxation-only mode, output and exit
        if(m_relaxation_only) {
            std::cout << "\n=== Relaxation-Only Mode: Outputting Results ===\n" << std::endl;

            m_sim->set_time(0.0);

            auto & p = m_sim->get_particles();
            const int num = m_sim->get_particle_num();
            const real gamma = m_param->physics.gamma;
            const real c_sound = gamma * (gamma - 1.0);
            const real alpha = m_param->av.alpha;

#pragma omp parallel for
            for(int i = 0; i < num; ++i) {
                p[i].alpha = alpha;
                p[i].balsara = 1.0;
                p[i].sound = std::sqrt(c_sound * p[i].ene);
            }

            auto tree = m_sim->get_tree();
            tree->resize(num);
            tree->make(p, num);

            m_pre->calculation(m_sim);
            m_fforce->calculation(m_sim);

            relax_meta.relaxation_step = m_relaxation_steps;
            relax_meta.accumulated_time = accumulated_time;

            m_output_manager->write_snapshot(m_sim, m_param, m_snapshot_counter++, &relax_meta);

            real kinetic, thermal, potential;
            compute_total_energies(kinetic, thermal, potential);
            m_output_manager->write_energy(0.0, kinetic, thermal, potential);

            std::cout << "=== K&I 2000 Relaxed Configuration Saved ===" << std::endl;
            std::cout << "=== Exiting (No Simulation Run) ===\n" << std::endl;

            return;
        }

        std::cout << "=== Starting Main Simulation ===\n" << std::endl;
        m_sim->set_time(0.0);
    }

    // Initialize relaxation for isothermal sphere samples (HVCCIsothermal10K and IsothermalBonnorEbert)
    if(m_use_relaxation && (m_sample == Sample::HVCCIsothermal10K || m_sample == Sample::IsothermalBonnorEbert)) {
        std::cout << "\n=== Initializing TRUE Bonnor-Ebert Relaxation ===" << std::endl;
        m_isothermal_relax = std::make_shared<IsothermalRelaxation>();

        // Get parameters from sample_parameters
        const real T_cloud = boost::any_cast<double>(m_sample_parameters["T_cloud"]);
        const real rho_center = boost::any_cast<double>(m_sample_parameters["rho_center_code"]);
        const real r_0 = boost::any_cast<double>(m_sample_parameters["r_c_code"]);  // r_0 stored as r_c_code
        const real R_cloud = boost::any_cast<double>(m_sample_parameters["R_cloud_code"]);
        const real P_ext = boost::any_cast<double>(m_sample_parameters["P_ext"]);
        const real density_to_n = boost::any_cast<double>(m_sample_parameters["density_to_n"]);
        const real xi_s = boost::any_cast<double>(m_sample_parameters["xi_s"]);
        const real mu = m_sample_parameters.count("mu") ?
                        boost::any_cast<double>(m_sample_parameters["mu"]) : 1.27;
        const real G = m_param->gravity.constant;

        IsothermalRelaxationParams relax_params;
        relax_params.T_cloud = T_cloud;
        relax_params.rho_center = rho_center;
        relax_params.r_0 = r_0;
        relax_params.R_cloud = R_cloud;
        relax_params.xi_s = xi_s;
        relax_params.P_ext = P_ext;
        relax_params.G = G;
        relax_params.density_to_n = density_to_n;
        relax_params.mu = mu;

        m_isothermal_relax->initialize(relax_params);
        std::cout << "=== TRUE Bonnor-Ebert Relaxation Initialized ===" << std::endl;

        // Initialize sound speed and tree for relaxation calculations
        {
            auto& particles = m_sim->get_particles();
            const int num_p = m_sim->get_particle_num();
            const real gamma = m_param->physics.gamma;
            const real c_sound_factor = gamma * (gamma - 1.0);
            for(int i = 0; i < num_p; ++i) {
                particles[i].sound = std::sqrt(c_sound_factor * particles[i].ene);
            }

            auto tree = m_sim->get_tree();
            tree->resize(num_p);
            tree->make(particles, num_p);
        }

        // GLASS pre-relaxation phase to uniformize particle spacing
        if(m_use_glass_relaxation) {
            std::cout << "\n=== Starting GLASS Pre-Relaxation (" << m_glass_relaxation_steps << " steps) ===" << std::endl;
            std::cout << "Target neighbor count: " << m_glass_target_neighbors << std::endl;

            auto& particles = m_sim->get_particles();
            const int num_p = m_sim->get_particle_num();
            auto* kernel = m_sim->get_kernel().get();
            auto* periodic = m_sim->get_periodic().get();

            // Get R_cloud for boundary conditions
            const real R_cloud = boost::any_cast<double>(m_sample_parameters["R_cloud_code"]);

            // GLASS uses repulsive forces between nearby particles
            // F_ij = k * (1 - r/h)^2 * r_hat for r < h
            // This pushes particles apart until uniform spacing achieved
            const real glass_strength = 0.1;  // Strength of repulsive force

            auto glass_start = std::chrono::steady_clock::now();
            double avg_step_time = 0.0;
            int timing_samples = 0;

            for(int step = 0; step < m_glass_relaxation_steps; ++step) {
                auto step_start = std::chrono::steady_clock::now();

                // Build tree for neighbor search
                auto tree = m_sim->get_tree();
                tree->resize(num_p);
                tree->make(particles, num_p);
                // Compute density and smoothing length
                m_pre->calculation(m_sim);

                // Compute GLASS repulsive forces
                const int max_neighbors = 200;

#pragma omp parallel for
                for(int i = 0; i < num_p; ++i) {
                    auto& p_i = particles[i];
                    if(p_i.is_ghost) continue;

                    std::vector<int> neighbor_list(max_neighbors);
                    int n_neighbor = tree->neighbor_search(p_i, neighbor_list, particles, true);
                    // GLASS repulsive force: push apart particles that are too close
                    vec_t f_glass(0.0);
                    for(int n = 0; n < n_neighbor; ++n) {
                        int j = neighbor_list[n];
                        auto& p_j = particles[j];

                        vec_t r_ij = periodic->calc_r_ij(p_i.pos, p_j.pos);
                        real r = std::abs(r_ij);
                        real h_avg = 0.5 * (p_i.sml + p_j.sml);

                        if(r < h_avg && r > 1e-10) {
                            // Repulsive force: stronger when particles closer
                            real q = r / h_avg;
                            real f_mag = glass_strength * (1.0 - q) * (1.0 - q) / r;
                            f_glass -= r_ij * f_mag;  // Repel (opposite direction of r_ij)
                        }
                    }
                    p_i.acc = f_glass;
                }

                // Calculate timestep based on maximum acceleration
                real max_acc = 0.0;
                for(int i = 0; i < num_p; ++i) {
                    if(!particles[i].is_ghost) {
                        max_acc = std::max(max_acc, std::abs(particles[i].acc));
                    }
                }
                real dt = (max_acc > 0) ? 0.01 * std::sqrt(particles[0].sml / max_acc) : 1e-6;
                dt = std::min(dt, real(0.001));

                // Update positions: Δx = 0.5 * a * dt²
#pragma omp parallel for
                for(int i = 0; i < num_p; ++i) {
                    if(particles[i].is_ghost) continue;

                    particles[i].pos[0] += 0.5 * particles[i].acc[0] * dt * dt;
                    particles[i].pos[1] += 0.5 * particles[i].acc[1] * dt * dt;
#if DIM == 3
                    particles[i].pos[2] += 0.5 * particles[i].acc[2] * dt * dt;
#endif
                    // Keep particles inside cloud radius
                    real r = std::abs(particles[i].pos);
                    if(r > R_cloud * 0.98) {
                        particles[i].pos *= (R_cloud * 0.98 / r);
                    }
                }

                // Timing and progress
                auto step_end = std::chrono::steady_clock::now();
                double step_time = std::chrono::duration<double>(step_end - step_start).count();
                timing_samples++;
                avg_step_time = ((timing_samples - 1) * avg_step_time + step_time) / timing_samples;

                // Progress bar every 1%
                if((step + 1) % (m_glass_relaxation_steps / 100 + 1) == 0 || step == m_glass_relaxation_steps - 1) {
                    int percent = (step + 1) * 100 / m_glass_relaxation_steps;
                    int bar_width = 40;
                    int filled = (step + 1) * bar_width / m_glass_relaxation_steps;

                    // Compute neighbor statistics
                    real avg_neighbors = 0.0;
                    real min_neighbors = 1e10, max_neighbors_stat = 0.0;
                    int cloud_count = 0;
                    for(int i = 0; i < num_p; ++i) {
                        if(!particles[i].is_ghost) {
                            avg_neighbors += particles[i].neighbor;
                            min_neighbors = std::min(min_neighbors, (real)particles[i].neighbor);
                            max_neighbors_stat = std::max(max_neighbors_stat, (real)particles[i].neighbor);
                            cloud_count++;
                        }
                    }
                    if(cloud_count > 0) avg_neighbors /= cloud_count;

                    int eta_secs = static_cast<int>((m_glass_relaxation_steps - step - 1) * avg_step_time);

                    std::cout << "\r[GLASS] [";
                    for(int i = 0; i < bar_width; ++i) std::cout << (i < filled ? "#" : " ");
                    std::cout << "] " << std::setw(3) << percent << "% "
                              << "N=" << std::fixed << std::setprecision(1) << avg_neighbors
                              << " [" << (int)min_neighbors << "-" << (int)max_neighbors_stat << "]"
                              << " a=" << std::scientific << std::setprecision(1) << max_acc
                              << " ETA:" << eta_secs << "s" << std::fixed << std::flush;
                }
            }
            std::cout << std::endl;
            std::cout << "=== GLASS Pre-Relaxation Complete ===" << std::endl;

            // Re-compute density after GLASS
            m_pre->calculation(m_sim);
        }

        // Run relaxation phase
        std::cout << "\n=== Starting Isothermal Relaxation Phase (" << m_relaxation_steps << " steps) ===" << std::endl;

        int start_step = 0;
        real accumulated_time = 0.0;

        // Metadata for relaxation snapshots
        OutputMetadata relax_meta;
        relax_meta.is_relaxation = true;
        relax_meta.relaxation_total_steps = m_relaxation_steps;
        relax_meta.accumulated_time = 0.0;

        // Progress tracking
        int output_counter = 0;
        int last_percent = -1;
        int last_sub_percent = -1;
        int target_step = start_step + m_relaxation_steps;

        // Timing for ETA
        auto start_time = std::chrono::steady_clock::now();
        double avg_step_time = 0.0;
        int timing_samples = 0;

        for(int step = start_step; step < target_step; ++step) {
            auto step_start = std::chrono::steady_clock::now();

            auto& particles = m_sim->get_particles();
            const int num_p = m_sim->get_particle_num();

            // Build tree
            auto tree = m_sim->get_tree();
            tree->resize(num_p);
            tree->make(particles, num_p);

            // ================================================================
            // HYBRID RELAXATION: Analytical force + velocity damping
            // ================================================================
            // - Subtract analytical equilibrium acceleration (prevents expansion)
            // - Apply velocity damping (removes oscillations)
            // ================================================================

            // Step 1: Compute SPH quantities
            m_pre->calculation(m_sim);
            m_fforce->calculation(m_sim);

            // Step 2: Apply isothermal relaxation (subtracts equilibrium forces)
            m_isothermal_relax->apply_relaxation(m_sim);

            // Step 3: Calculate timestep
            m_timestep->calculation(m_sim);
            real dt = m_sim->get_dt();
            dt = std::max(dt, real(1.0e-10));

            // ================================================================
            // QUASI-STATIC RELAXATION (same as Lane-Emden):
            // - Zero velocities every step (quasi-static)
            // - Kinematic integration: Δx = ½at² (with v₀=0)
            // - External pressure confinement via envelope particles
            // - Remove escaping particles beyond R_cloud
            // ================================================================

            // Step 4: Update positions with acceleration (velocities zeroed)
            auto * periodic = m_sim->get_periodic().get();
#pragma omp parallel for
            for(int i = 0; i < num_p; ++i) {
                if(particles[i].is_ghost) continue;

                // Zero velocities (quasi-static relaxation)
                particles[i].vel[0] = 0.0;
                particles[i].vel[1] = 0.0;
#if DIM == 3
                particles[i].vel[2] = 0.0;
#endif

                // Kinematic integration with v₀=0: Δx = ½at²
                const real half_dt2 = 0.5 * dt * dt;
                particles[i].pos[0] += particles[i].acc[0] * half_dt2;
                particles[i].pos[1] += particles[i].acc[1] * half_dt2;
#if DIM == 3
                particles[i].pos[2] += particles[i].acc[2] * half_dt2;
#endif

                periodic->apply(particles[i].pos);
            }

            // Remove particles that escaped beyond cloud boundary
            // External pressure confinement via envelope particles keeps most particles inside
            m_isothermal_relax->remove_escaping_particles(m_sim, 1.1);

            accumulated_time += dt;

            // Timing
            auto step_end = std::chrono::steady_clock::now();
            double step_time = std::chrono::duration<double>(step_end - step_start).count();
            timing_samples++;
            avg_step_time = ((timing_samples - 1) * avg_step_time + step_time) / timing_samples;

            int percent = (step - start_step + 1) * 100 / m_relaxation_steps;
            int sub_percent = ((step - start_step + 1) * 1000 / m_relaxation_steps) % 10;

            if(percent != last_percent || sub_percent != last_sub_percent) {
                int bar_width = 50;
                int filled = (step - start_step + 1) * bar_width / m_relaxation_steps;

                // Calculate max acceleration (Lane-Emden convergence criterion)
                // With velocities zeroed, we monitor |a_net| = |a_SPH - a_eq| -> 0
                real max_acc = 0.0;
                for(int i = 0; i < num_p; ++i) {
                    if(!particles[i].is_ghost) {
                        real a2 = particles[i].acc[0] * particles[i].acc[0]
                                + particles[i].acc[1] * particles[i].acc[1];
#if DIM == 3
                        a2 += particles[i].acc[2] * particles[i].acc[2];
#endif
                        max_acc = std::max(max_acc, std::sqrt(a2));
                    }
                }

                // Calculate ETA
                int steps_remaining = m_relaxation_steps - (step - start_step + 1);
                double eta_seconds = steps_remaining * avg_step_time;
                int eta_mins = static_cast<int>(eta_seconds / 60);
                int eta_secs = static_cast<int>(eta_seconds) % 60;

                std::cout << "\r[";
                for(int i = 0; i < bar_width; ++i) {
                    std::cout << (i < filled ? "#" : " ");
                }
                std::cout << "] " << std::setw(3) << percent << "." << sub_percent << "% "
                          << step - start_step + 1 << "/" << m_relaxation_steps
                          << " ETA:" << eta_mins << "m" << std::setw(2) << std::setfill('0') << eta_secs << "s"
                          << std::setfill(' ')
                          << " |a|=" << std::scientific << std::setprecision(2) << max_acc
                          << std::fixed
                          << std::flush;

                last_percent = percent;
                last_sub_percent = sub_percent;
            }

            // Output snapshots
            if(step % m_relaxation_output_freq == 0 || step == m_relaxation_steps - 1) {
                std::cout << std::endl;
                m_sim->set_time(accumulated_time);

                // Update sound speed before output
                const real gamma = m_param->physics.gamma;
                const real c_sound_factor = gamma * (gamma - 1.0);
                for(int i = 0; i < num_p; ++i) {
                    particles[i].sound = std::sqrt(c_sound_factor * particles[i].ene);
                }

                relax_meta.relaxation_step = step;
                relax_meta.accumulated_time = accumulated_time;

                m_output_manager->write_snapshot(m_sim, m_param, m_snapshot_counter++, &relax_meta);
                output_counter++;

                std::cout << "Wrote snapshot " << m_snapshot_counter - 1
                          << " at step " << step << " (t=" << accumulated_time << ")" << std::endl;
            }
        }

        std::cout << "\n=== Isothermal Relaxation Complete ===" << std::endl;
        std::cout << "Total relaxation time: " << accumulated_time << " [code units]" << std::endl;

        if(m_relaxation_only) {
            // Output final relaxed configuration
            auto& p = m_sim->get_particles();
            const int num = m_sim->get_particle_num();
            const real gamma = m_param->physics.gamma;
            const real c_sound_factor = gamma * (gamma - 1.0);

            auto tree = m_sim->get_tree();
            tree->resize(num);
            tree->make(p, num);

            m_pre->calculation(m_sim);
            m_fforce->calculation(m_sim);

            relax_meta.relaxation_step = m_relaxation_steps;
            relax_meta.accumulated_time = accumulated_time;

            m_output_manager->write_snapshot(m_sim, m_param, m_snapshot_counter++, &relax_meta);

            real kinetic, thermal, potential;
            compute_total_energies(kinetic, thermal, potential);
            m_output_manager->write_energy(0.0, kinetic, thermal, potential);

            std::cout << "=== Isothermal Relaxed Configuration Saved ===" << std::endl;
            std::cout << "=== Exiting (No Simulation Run) ===\n" << std::endl;

            return;
        }

        // ================================================================
        // REMOVE GHOST PARTICLES AFTER RELAXATION
        // ================================================================
        // Ghost particles provided external pressure confinement during
        // relaxation. For the actual simulation (e.g., BH flyby), we remove
        // them so the cloud can freely evolve without artificial boundaries.
        // ================================================================
        {
            auto& particles = m_sim->get_particles();
            const int num_before = m_sim->get_particle_num();

            // Count ghost particles before removal
            int ghost_count = 0;
            for(int i = 0; i < num_before; ++i) {
                if(particles[i].is_ghost) ghost_count++;
            }

            if(ghost_count > 0) {
                std::cout << "\n=== Removing Ghost Particles ===" << std::endl;
                std::cout << "Particles before: " << num_before << std::endl;
                std::cout << "Ghost particles: " << ghost_count << std::endl;

                // Remove ghost particles using erase-remove idiom
                auto new_end = std::remove_if(particles.begin(), particles.end(),
                    [](const SPHParticle& p) { return p.is_ghost; });
                particles.erase(new_end, particles.end());

                const int num_after = static_cast<int>(particles.size());
                m_sim->get_particle_num() = num_after;

                std::cout << "Particles after: " << num_after << std::endl;
                std::cout << "Removed: " << (num_before - num_after) << " ghost particles" << std::endl;
                std::cout << "=================================\n" << std::endl;
            }
        }

        std::cout << "=== Starting Main Simulation ===\n" << std::endl;
        m_sim->set_time(0.0);
    }

    auto & p = m_sim->get_particles();
    const int num = m_sim->get_particle_num();
    const real gamma = m_param->physics.gamma;
    const real c_sound = gamma * (gamma - 1.0);

    assert(p.size() == num);
    const real alpha = m_param->av.alpha;
#pragma omp parallel for
    for(int i = 0; i < num; ++i) {
        p[i].alpha = alpha;
        p[i].balsara = 1.0;
        p[i].sound = std::sqrt(c_sound * p[i].ene);
    }

    auto tree = m_sim->get_tree();
    tree->resize(num);
    tree->make(p, num);

    m_pre->calculation(m_sim);
    
    // For SR/GR-GSPH and SRMHD hydro-only, update ghost particles after computing N but before forces
    if(m_param->type == SPHType::SRGSPH || m_param->type == SPHType::GRGSPH ||
       (m_param->type == SPHType::SRMHD && !m_param->mhd.use_mhd)) {
        update_ghost_particles();
    }
    
    // Compute gravity BEFORE fluid force so grav_acc is available for
    // gravity-aware Riemann solver in GSPH
    if(m_param->gravity.is_valid) {
        m_gforce->calculation(m_sim);
    }
    
    m_fforce->calculation(m_sim);
    
    // Add gravity acceleration to total acceleration AFTER fluid force
    if(m_param->gravity.is_valid) {
        auto & p = m_sim->get_particles();
        const int num = m_sim->get_particle_num();
#pragma omp parallel for
        for(int i = 0; i < num; ++i) {
            p[i].acc += p[i].grav_acc;
        }
    }
    
    // Apply external BH force if enabled (must be called in initialize too!)
    if(m_use_external_bh) {
        m_external_bh->calculation(m_sim);
    }
}

void Solver::integrate()
{
    m_timestep->calculation(m_sim);

    predict();
    m_sim->make_tree();
    // First compute densities and smoothing lengths for real particles
    m_pre->calculation(m_sim);
    
    // Then update ghost particles to mirror the freshly computed values
    // This ensures ghosts have current N, h, etc. for force calculation
    update_ghost_particles();
    
    // Compute gravity BEFORE fluid force so that grav_acc is available
    // for the gravity-aware Riemann solver in GSPH
    if(m_param->gravity.is_valid) {
        m_gforce->calculation(m_sim);
    }
    
    m_fforce->calculation(m_sim);
    
    // Add gravity acceleration to total acceleration AFTER fluid force
    // (fluid force overwrites acc, so gravity must be added after)
    if(m_param->gravity.is_valid) {
        auto & p = m_sim->get_particles();
        const int num = m_sim->get_particle_num();
#pragma omp parallel for
        for(int i = 0; i < num; ++i) {
            p[i].acc += p[i].grav_acc;
        }
    }
    
    // Apply external BH force if enabled
    if(m_use_external_bh) {
        m_external_bh->calculation(m_sim);

        // Apply sink accretion (remove particles that fall into the BH)
        // This is called AFTER force calculation to use updated velocities
        // Particles inside r_sink that are bound and moving inward are removed.
        // BH mass remains fixed (fixed-potential approximation).
        m_external_bh->apply_sink_accretion(m_sim);

        // Update BH position if it's moving (typically fixed at origin)
        const real dt = m_sim->get_dt();
        m_external_bh->update_position(dt);
    }

    correct();
}

void Solver::predict()
{
    auto & p = m_sim->get_particles();
    const int num = m_sim->get_particle_num();
    auto * periodic = m_sim->get_periodic().get();
    const real dt = m_sim->get_dt();
    const real gamma = m_param->physics.gamma;
    const real c_sound = gamma * (gamma - 1.0);

    assert(p.size() == num);

    // SR/GR-GSPH: Integrate CONSERVED variables (S, e, N)
    // Standard SPH: Integrate PRIMITIVE variables (v, u, ρ)
    // SRMHD hydro-only: Use SR-GSPH integration (proven to work)
    const bool use_srgsph_integration = (m_param->type == SPHType::SRGSPH ||
                                         m_param->type == SPHType::GRGSPH ||
                                         (m_param->type == SPHType::SRMHD && !m_param->mhd.use_mhd));
    if(use_srgsph_integration) {
        const real c_speed = m_param->srgsph.c_speed;

#pragma omp parallel for
        for(int i = 0; i < num; ++i) {
            // Skip ghost particles - they are fixed and don't evolve
            if(p[i].is_ghost) continue;
            
            // === SR-GSPH PREDICTOR-CORRECTOR TIME INTEGRATION ===
            // Implements 2nd-order Heun's method (predictor-corrector)
            //
            // PREDICTOR: Predict S, e at t^(n+1) using derivatives at t^n
            // Then force calculation computes NEW derivatives at t^(n+1)
            // Then CORRECTOR uses average of old and new derivatives

            // STEP 1: Save old derivatives for corrector step
            p[i].dS_old = p[i].dS;
            p[i].de_old = p[i].de;
            p[i].dS_t_old = p[i].dS_t;  // Tangent momentum derivative

            // STEP 2: Predict with Euler step (half-step for position)
            vec_t S_half = p[i].S + p[i].dS * (0.5 * dt);
            real e_half = p[i].e + p[i].de * (0.5 * dt);
            real S_t_half = p[i].S_t + p[i].dS_t * (0.5 * dt);  // Tangent momentum

            // CRITICAL: Floor energy BEFORE primitive recovery to prevent instability
            // If e < 1, primitive recovery gives H→1, P→0, and v = S/(γH) can explode!
            if (e_half < 1.0) {
                e_half = 1.0 + 1e-10;
            }

            // Limit half-step S/e ratio to prevent superluminal velocities
            // Tangent momentum is passive and handled separately
            const real S_normal_mag = std::sqrt(inner_product(S_half, S_half));
            const real S_half_ratio = S_normal_mag / std::max(e_half, 1.0e-10);
            if(S_half_ratio > 0.9999) {  // Only intervene if approaching c
                const real scale = 0.9999 / S_half_ratio;
                S_half = S_half * scale;
                // Don't scale S_t_half - it's passive
            }

            // STEP 3: Full Euler prediction
            // dS and de are already PER-BARYON (divided by nu in normalize_sr_derivatives)
            // Python reference: p.S += (p.dSdt / p.nu) * dt
            // Since normalize_sr_derivatives already divided by nu, we just multiply by dt
            p[i].S += p[i].dS * dt;      // Predicted S using OLD derivative
            p[i].e += p[i].de * dt;      // Predicted e using OLD derivative
            p[i].S_t += p[i].dS_t * dt;  // Predicted S_t using OLD derivative

            // CRITICAL: Floor energy BEFORE primitive recovery (full step)
            // This prevents the runaway instability where low e → low H → high v → lower e
            if (p[i].e < 1.0) {
                p[i].e = 1.0 + 1e-10;
            }

            // Limit full-step S/e ratio to prevent superluminal velocities
            {
                const real S_mag = std::sqrt(inner_product(p[i].S, p[i].S));
                const real S_ratio = S_mag / std::max(p[i].e, 1.0e-10);
                if (S_ratio > 0.9999) {
                    const real scale = 0.9999 / S_ratio;
                    p[i].S = p[i].S * scale;
                }
            }

            // Recover primitive variables at half-step for position update
            // Use conserved S_t to recover v_t (S_t = γHv_t is conserved, not v_t)
            auto prim_half = srgsph::PrimitiveRecovery::conserved_to_primitive_with_tangent(
                S_half, S_t_half, e_half, p[i].N, gamma, c_speed
            );
            
            // Update position using half-step velocity (drift)
            p[i].pos += prim_half.vel * dt;
            periodic->apply(p[i].pos);
            
            // Recover primitive variables at full step
            // Use conserved S_t to recover v_t (S_t = γHv_t is conserved, not v_t)
            auto prim_full = srgsph::PrimitiveRecovery::conserved_to_primitive_with_tangent(
                p[i].S, p[i].S_t, p[i].e, p[i].N, gamma, c_speed
            );
            
            // Store PRIMITIVE variables for output and next step's Riemann solver
            // NOTE: These fields have DIFFERENT meanings in SRGSPH vs standard SPH!
            //   - vel: primitive velocity v (NOT time-integrated, just recovered)
            //   - ene: primitive internal energy u (NOT conserved energy e!)
            //   - dens: rest-frame density n (recovered from N/γ in primitive recovery)
            p[i].vel = prim_full.vel;       // Primitive velocity v
            
            // Update vel_t from primitive recovery (v_t is recovered from conserved S_t)
            p[i].vel_t = prim_full.vel_t;

            // Safety clamp: ensure |v|² + v_t² < 1 (subluminal)
            // In DIM=1: v² = v_x² + v_t² (tangent velocity)
            // In DIM≥2: v² = |vel|² (no tangent velocity)
            {
                real v2 = inner_product(p[i].vel, p[i].vel);
#if DIM == 1
                const real v_t2 = p[i].vel_t * p[i].vel_t;
                v2 += v_t2;
#endif
                if (v2 > 0.9999) {
                    static int clamp_count = 0;
                    const real scale = std::sqrt(0.99 / v2);
                    if (clamp_count < 10) {
                        WRITE_LOG << "[PREDICT CLAMP] particle " << i
                                  << " |v|²=" << v2
                                  << " clamped by factor " << scale;
                        ++clamp_count;
                    }
                    p[i].vel = p[i].vel * scale;
#if DIM == 1
                    p[i].vel_t *= scale;
#endif
                }
            }

            p[i].vel_p = prim_half.vel;     // Half-step velocity (for position update)
            p[i].ene = prim_full.pressure / ((gamma - 1.0) * prim_full.density);  // u = P/[(γ-1)n]
            p[i].ene_p = prim_half.pressure / ((gamma - 1.0) * prim_half.density);
            p[i].pres = prim_full.pressure;
            p[i].dens = prim_full.density;  // Rest-frame density n (for output)
            p[i].sound = prim_full.sound_speed;
            p[i].gamma_lor = prim_full.gamma_lor;
            p[i].enthalpy = prim_full.enthalpy;

            // NOTE: S_t is CONSERVED (dS_t/dt = 0), do NOT update it here!
            // v_t is recovered from S_t/(γH) in primitive recovery above.
        }
    } else if(m_param->type == SPHType::GSPMHD) {
        // === GSPMHD TIME INTEGRATION (Iwasaki & Inutsuka 2011) ===
        // Standard SPH integration plus B field and MHD velocity components
#pragma omp parallel for
        for(int i = 0; i < num; ++i) {
            if(p[i].is_ghost) continue;  // Ghost particles are fixed for stability

            // Save old derivatives for predictor-corrector
            p[i].dB_old = p[i].dB;

            // k -> k+1/2 (half-step for position update)
            p[i].vel_p = p[i].vel + p[i].acc * (0.5 * dt);
            p[i].ene_p = p[i].ene + p[i].dene * (0.5 * dt);

            // k -> k+1 (full step for velocities, energy, B)
            p[i].pos += p[i].vel_p * dt;
            p[i].vel += p[i].acc * dt;
            p[i].vy_mhd += p[i].acc_y_mhd * dt;  // MHD transverse velocity
            p[i].vz_mhd += p[i].acc_z_mhd * dt;
            p[i].ene += p[i].dene * dt;
            p[i].B += p[i].dB * dt;  // B field evolution

            // Floor internal energy
            if(p[i].ene < 1e-30) p[i].ene = 1e-30;
            p[i].sound = std::sqrt(c_sound * p[i].ene);

            periodic->apply(p[i].pos);
        }
    } else if(m_param->type == SPHType::SRMHD) {
        // === SRMHD TIME INTEGRATION (matching SR-GSPH structure) ===
        // CRITICAL: Must integrate CONSERVED variables (S, e, S_t) like SR-GSPH!
        // Also: do primitive recovery IN predictor to match SR-GSPH behavior.
        const real c_speed = m_param->srgsph.c_speed;

#pragma omp parallel for
        for(int i = 0; i < num; ++i) {
            if(p[i].is_ghost) continue;

            // STEP 1: Save old derivatives for corrector step
            p[i].dS_old = p[i].dS;
            p[i].de_old = p[i].de;
            p[i].dS_t_old = p[i].dS_t;
            p[i].dB_old = p[i].dB;

            // STEP 2: Compute half-step conserved variables
            vec_t S_half = p[i].S + p[i].dS * (0.5 * dt);
            real e_half = p[i].e + p[i].de * (0.5 * dt);
            real S_t_half = p[i].S_t + p[i].dS_t * (0.5 * dt);

            // CRITICAL: Floor e_half BEFORE primitive recovery to prevent instability
            if (e_half < 1.0) {
                e_half = 1.0 + 1e-10;
            }

            // Limit half-step S/e ratio to prevent superluminal velocities
            const real S_normal_mag = std::sqrt(inner_product(S_half, S_half));
            const real S_half_ratio = S_normal_mag / std::max(e_half, 1.0e-10);
            if (S_half_ratio > 0.9999) {
                const real scale = 0.9999 / S_half_ratio;
                S_half = S_half * scale;
            }

            // STEP 3: Full step prediction
            p[i].S += p[i].dS * dt;
            p[i].e += p[i].de * dt;
            p[i].S_t += p[i].dS_t * dt;
            p[i].B += p[i].dB * dt;

            // MHD transverse velocity evolution (1D only)
            // In relativistic MHD, acc_y_mhd is dS_y/dt (momentum rate), not dv_y/dt
            // Convert: dv_y/dt = dS_y/dt / (γH)
            // Also limit transverse velocity to keep total v² < 1
#if DIM == 1
            if (m_param->mhd.use_mhd) {
                const real gamma_H = p[i].gamma_lor * p[i].enthalpy;
                const real inv_gamma_H = 1.0 / std::max(gamma_H, 1.0);
                p[i].vy_mhd += p[i].acc_y_mhd * inv_gamma_H * dt;
                p[i].vz_mhd += p[i].acc_z_mhd * inv_gamma_H * dt;

                // Limit transverse velocity to keep total velocity subluminal
                // v² = vx² + vy² + vz² < 1
                const real vx2 = p[i].vel[0] * p[i].vel[0];
                const real vt2 = p[i].vy_mhd * p[i].vy_mhd + p[i].vz_mhd * p[i].vz_mhd;
                const real v2_max = 0.9999 * 0.9999 - vx2;  // Max allowed vt²
                if (vt2 > v2_max && vt2 > 1e-30) {
                    const real scale = std::sqrt(std::max(v2_max, 0.0) / vt2);
                    p[i].vy_mhd *= scale;
                    p[i].vz_mhd *= scale;
                }
            }
#endif

            // CRITICAL: Floor canonical energy BEFORE primitive recovery
            if (p[i].e < 1.0) p[i].e = 1.0 + 1e-10;

            // Limit full-step S/e ratio to prevent superluminal velocities
            {
                const real S_mag = std::sqrt(inner_product(p[i].S, p[i].S));
                const real S_ratio = S_mag / std::max(p[i].e, 1.0e-10);
                if (S_ratio > 0.9999) {
                    const real scale = 0.9999 / S_ratio;
                    p[i].S = p[i].S * scale;
                }
            }

            // STEP 4: Recover primitives at half-step for position update
            auto prim_half = srmhd::PrimitiveRecovery::conserved_to_primitive_with_tangent(
                S_half, S_t_half, e_half, p[i].N, gamma, c_speed
            );

            // STEP 5: Position update using half-step velocity
            p[i].pos += prim_half.vel * dt;
            periodic->apply(p[i].pos);

            // STEP 6: Recover primitives at full step
            auto prim_full = srmhd::PrimitiveRecovery::conserved_to_primitive_with_tangent(
                p[i].S, p[i].S_t, p[i].e, p[i].N, gamma, c_speed
            );

            // STEP 7: Store primitives
            p[i].vel = prim_full.vel;
            p[i].vel_t = prim_full.vel_t;
            p[i].vel_p = prim_half.vel;
            p[i].dens = prim_full.density;
            p[i].pres = prim_full.pressure;
            p[i].gamma_lor = prim_full.gamma_lor;
            p[i].enthalpy = prim_full.enthalpy;
            p[i].sound = prim_full.sound_speed;
            p[i].ene = prim_full.pressure / ((gamma - 1.0) * prim_full.density);

            // STEP 8: Velocity clamp (matching SR-GSPH)
            {
                real v2 = inner_product(p[i].vel, p[i].vel);
#if DIM == 1
                v2 += p[i].vel_t * p[i].vel_t;
#endif
                if (v2 > 0.9999) {
                    const real scale = std::sqrt(0.99 / v2);
                    p[i].vel = p[i].vel * scale;
#if DIM == 1
                    p[i].vel_t *= scale;
#endif
                }
            }
        }
    } else {
        // === STANDARD SPH TIME INTEGRATION ===
        // Integrate PRIMITIVE variables: v (velocity), u (internal energy)
#pragma omp parallel for
        for(int i = 0; i < num; ++i) {
            // Skip ghost particles - they are fixed and don't evolve
            if(p[i].is_ghost) continue;

            // k -> k+1/2
            p[i].vel_p = p[i].vel + p[i].acc * (0.5 * dt);
            p[i].ene_p = p[i].ene + p[i].dene * (0.5 * dt);

            // k -> k+1
            p[i].pos += p[i].vel_p * dt;
            p[i].vel += p[i].acc * dt;
            p[i].ene += p[i].dene * dt;
            p[i].sound = std::sqrt(c_sound * p[i].ene);

            periodic->apply(p[i].pos);
            
            // Apply spherical reflecting boundary if enabled
            if(m_spherical_boundary_radius > 0.0) {
                real r2 = p[i].pos[0] * p[i].pos[0] + p[i].pos[1] * p[i].pos[1];
#if DIM == 3
                r2 += p[i].pos[2] * p[i].pos[2];
#endif
                real r = std::sqrt(r2);
                if(r > m_spherical_boundary_radius) {
                    // Reflecting boundary: push particle back and reflect velocity
                    real scale = m_spherical_boundary_radius * 0.99 / r;
                    p[i].pos[0] *= scale;
                    p[i].pos[1] *= scale;
#if DIM == 3
                    p[i].pos[2] *= scale;
#endif
                    // Also reduce radial velocity component
                    real v_r = (p[i].vel[0] * p[i].pos[0] + p[i].vel[1] * p[i].pos[1]) / r;
#if DIM == 3
                    v_r += p[i].vel[2] * p[i].pos[2] / r;
#endif
                    if(v_r > 0) {  // Only if moving outward
                        // Reflect radial velocity component
                        real n_x = p[i].pos[0] / r;
                        real n_y = p[i].pos[1] / r;
                        p[i].vel[0] -= 2.0 * v_r * n_x;
                        p[i].vel[1] -= 2.0 * v_r * n_y;
#if DIM == 3
                        real n_z = p[i].pos[2] / r;
                        p[i].vel[2] -= 2.0 * v_r * n_z;
#endif
                    }
                }
            }
        }
    }
}

void Solver::correct()
{
    auto & p = m_sim->get_particles();
    const int num = m_sim->get_particle_num();
    const real dt = m_sim->get_dt();
    const real gamma = m_param->physics.gamma;
    const real c_sound = gamma * (gamma - 1.0);

    assert(p.size() == num);

    // SR/GR-GSPH uses different correction (conserved variables)
    // SRMHD hydro-only: Use SR-GSPH integration (proven to work)
    const bool use_srgsph_correction = (m_param->type == SPHType::SRGSPH ||
                                        m_param->type == SPHType::GRGSPH ||
                                        (m_param->type == SPHType::SRMHD && !m_param->mhd.use_mhd));
    if(use_srgsph_correction) {
        const real c_speed = m_param->srgsph.c_speed;
        
#pragma omp parallel for
        for(int i = 0; i < num; ++i) {
            // Skip ghost particles - they are fixed and don't evolve
            if(p[i].is_ghost) continue;
            
            // === SR-GSPH CORRECTOR STEP ===
            // Now we have NEW derivatives (dS, de) from force calculation at predicted state
            // Apply corrector: use AVERAGE of old and new derivatives
            //
            // Current state: S_pred = S^n + dS_old × dt  (from predict step)
            // Corrector:     S_corr = S^n + 0.5 × (dS_old + dS_new) × dt
            //
            // Algebra: S_corr = S_pred - dS_old×dt + 0.5×(dS_old + dS_new)×dt
            //                 = S_pred + 0.5×(dS_new - dS_old)×dt

            // Check if dS_old is zero (first timestep)
            const real dS_old_mag = std::sqrt(inner_product(p[i].dS_old, p[i].dS_old));
            if(dS_old_mag < 1.0e-20 && std::abs(p[i].de_old) < 1.0e-20) {
                // First timestep: dS_old not valid, keep Euler prediction
                // No correction needed - already have S = S^n + dS*dt from predict
            } else {
                // Normal corrector step
                // dS and de are already PER-BARYON (divided by nu in normalize_sr_derivatives)
                p[i].S = p[i].S + (p[i].dS - p[i].dS_old) * (0.5 * dt);
                p[i].e = p[i].e + (p[i].de - p[i].de_old) * (0.5 * dt);
                p[i].S_t = p[i].S_t + (p[i].dS_t - p[i].dS_t_old) * (0.5 * dt);  // Tangent momentum
            }

            // CRITICAL: Floor energy to prevent runaway instability
            // For relativistic canonical energy, e = γH - P/(Nc²) ≥ 1 for physical states
            // If e < 1, primitive recovery gives H → 1, P → 0, and v = S/(γH) can explode!
            if(p[i].e < 1.0) {
                p[i].e = 1.0 + 1e-10;
            }

            // Safety: Limit S/e ratio to prevent superluminal velocities
            // For tangent velocity tests, only consider normal momentum
            // (tangent momentum is passive and handled separately in primitive recovery)
            const real S_normal_mag = std::sqrt(inner_product(p[i].S, p[i].S));
            const real S_over_e = S_normal_mag / std::max(p[i].e, 1.0e-10);
            if(S_over_e > 0.9999) {  // Only intervene if approaching c
                const real scale = 0.9999 / S_over_e;
                p[i].S = p[i].S * scale;
                // Don't scale S_t - it's a passive scalar conserved independently
            }

            // Floor on N
            p[i].N = std::max(p[i].N, 1.0e-6);

            // Recover primitive variables from CORRECTED conserved variables
            // Use conserved S_t to recover v_t (S_t = γHv_t is conserved, not v_t)
            auto prim = sph::srgsph::PrimitiveRecovery::conserved_to_primitive_with_tangent(
                p[i].S, p[i].S_t, p[i].e, p[i].N, gamma, c_speed
            );
            
            // Update primitive variables from recovered state
            p[i].vel = prim.vel;

            // Update vel_t from primitive recovery (v_t is recovered from conserved S_t)
            p[i].vel_t = prim.vel_t;

            // Safety clamp: ensure |v_x| + v_t^2 < 1 (subluminal with tangent velocity)
            {
                const real v_t2 = p[i].vel_t * p[i].vel_t;
                const real max_vx = std::sqrt(std::max(0.9999 - v_t2, 0.01));
                if (std::abs(p[i].vel[0]) > max_vx) {
                    static int clamp_count = 0;
                    if (clamp_count < 10) {
                        WRITE_LOG << "[CORRECT CLAMP] particle " << i
                                  << " vel[0]=" << p[i].vel[0]
                                  << " clamped to " << std::copysign(max_vx, p[i].vel[0])
                                  << " (v_t=" << p[i].vel_t << ")";
                        ++clamp_count;
                    }
                    p[i].vel[0] = std::copysign(max_vx, p[i].vel[0]);
                }
            }

            p[i].ene = prim.pressure / ((gamma - 1.0) * prim.density);
            p[i].pres = prim.pressure;
            p[i].dens = prim.density;
            p[i].sound = prim.sound_speed;
            p[i].gamma_lor = prim.gamma_lor;
            p[i].enthalpy = prim.enthalpy;

            // NOTE: S_t is CONSERVED (dS_t/dt = 0), do NOT update it here!
            // v_t is recovered from S_t/(γH) in primitive recovery above.
        }
    } else if(m_param->type == SPHType::GSPMHD) {
        // GSPMHD correction with B field
#pragma omp parallel for
        for(int i = 0; i < num; ++i) {
            if(p[i].is_ghost) continue;

            // Corrector: v^{n+1} = v^{n+1/2} + 0.5 * a^{n+1} * dt
            p[i].vel = p[i].vel_p + p[i].acc * (0.5 * dt);
            p[i].ene = p[i].ene_p + p[i].dene * (0.5 * dt);

            // B field corrector (use average of old and new derivatives)
            // B^{n+1}_corr = B^{n+1}_pred + 0.5 * (dB^{n+1} - dB^n) * dt
            p[i].B += (p[i].dB - p[i].dB_old) * (0.5 * dt);

            // Floor internal energy
            if(p[i].ene < 1e-30) p[i].ene = 1e-30;
            p[i].sound = std::sqrt(c_sound * p[i].ene);
        }
    } else if(m_param->type == SPHType::SRMHD) {
        // === SRMHD CORRECTION: Conserved variables (S, e, S_t) ===
        // CRITICAL: Must correct CONSERVED variables, not primitives!
        // The new derivatives (dS, de, dS_t) are computed from the predicted state.
        // Correction: apply 0.5*(dS_new - dS_old)*dt to the predicted S, e, S_t
        // Primitives (vel, ene) will be recovered in next pre_interaction call.

#pragma omp parallel for
        for(int i = 0; i < num; ++i) {
            if(p[i].is_ghost) continue;

            // Corrector for CONSERVED variables:
            // X^{n+1} = X^{predicted} + 0.5 * (dX^{new} - dX^{old}) * dt
            p[i].S = p[i].S + (p[i].dS - p[i].dS_old) * (0.5 * dt);
            p[i].e = p[i].e + (p[i].de - p[i].de_old) * (0.5 * dt);
            p[i].S_t = p[i].S_t + (p[i].dS_t - p[i].dS_t_old) * (0.5 * dt);

            // B field correction (for MHD)
            p[i].B += (p[i].dB - p[i].dB_old) * (0.5 * dt);

            // Safety: floor canonical energy
            // e = γH - P/(Nc²) ≈ 1 for non-relativistic cold gas
            if (p[i].e < 1.0) p[i].e = 1.0 + 1e-10;
        }
        // NOTE: Primitive variables (vel, pres, dens) will be recovered
        // from (S, e, S_t, N) in the next pre_interaction call.
    } else {
        // Standard SPH correction
#pragma omp parallel for
        for(int i = 0; i < num; ++i) {
            // Skip ghost particles - they are fixed and don't evolve
            if(p[i].is_ghost) continue;

            p[i].vel = p[i].vel_p + p[i].acc * (0.5 * dt);
            p[i].ene = p[i].ene_p + p[i].dene * (0.5 * dt);
            p[i].sound = std::sqrt(c_sound * p[i].ene);
        }
    }

    // Update ghost particles by mirroring nearest real particles
    update_ghost_particles();
}

void Solver::update_ghost_particles()
{
    // For INFLOW boundary: ghost particles maintain their initial state
    // (position, velocity, density, pressure, B field, etc.)
    // This provides fixed boundary conditions.
    //
    // Ghost particles do NOT move and do NOT get updated.
    // They provide SPH kernel support at the boundaries.
    //
    // For shock tube tests: The gap that forms between ghosts and real
    // particles when bulk flow moves them away is physical (rarefaction wave).

    // Update shock tube ghost particles if this is a shock tube sample
    if(m_sample == Sample::ShockTube || m_sample == Sample::ShockTube2D || m_sample == Sample::ShockTube3D) {
        const real gamma = m_param->physics.gamma;
        auto& particles = m_sim->get_particles();
        update_shock_tube_ghost_particles(particles, gamma);
    }
}

void Solver::make_initial_condition()
{
    std::cout << "make_initial_condition() called, m_sample = " << static_cast<int>(m_sample) << std::endl;
    std::cout.flush();
    switch(m_sample) {
#define MAKE_SAMPLE(a, b) case a: std::cout << "Calling make_" #b "()" << std::endl; std::cout.flush(); make_##b(); break
        MAKE_SAMPLE(Sample::ShockTube, shock_tube);
        MAKE_SAMPLE(Sample::ShockTube2D, shock_tube_2d);
        MAKE_SAMPLE(Sample::ShockTube3D, shock_tube_3d);
        MAKE_SAMPLE(Sample::ShockTube3DCubic, shock_tube_3d_cubic);
        MAKE_SAMPLE(Sample::Vacuum, vacuum);
        MAKE_SAMPLE(Sample::StrongShock, strong_shock);
        MAKE_SAMPLE(Sample::PressureEquilibrium, pressure_equilibrium);
        MAKE_SAMPLE(Sample::GreshoChanVortex, gresho_chan_vortex);
        MAKE_SAMPLE(Sample::PairingInstability, pairing_instability);
        MAKE_SAMPLE(Sample::HydroStatic, hydrostatic);
        MAKE_SAMPLE(Sample::KHI, khi);
        MAKE_SAMPLE(Sample::ISMCooling1D, ism_cooling_1d);
        MAKE_SAMPLE(Sample::Evrard, evrard);
        MAKE_SAMPLE(Sample::EvrardColdCollapse, evrard_cold_collapse);
        MAKE_SAMPLE(Sample::LaneEmden, lane_emden);
        MAKE_SAMPLE(Sample::LaneEmdenCylinder, lane_emden_cylinder);
        MAKE_SAMPLE(Sample::PolytropicSlab2D, polytropic_slab_2d);
        MAKE_SAMPLE(Sample::Sedov, sedov);
        MAKE_SAMPLE(Sample::SRSod, sr_sod);
        MAKE_SAMPLE(Sample::SRTangentVelocity, sr_tangent_velocity);
        MAKE_SAMPLE(Sample::SRRosswog, sr_rosswog);
        MAKE_SAMPLE(Sample::GRSchwarzschildShock, gr_schwarzschild_shock);
        MAKE_SAMPLE(Sample::GRGeodesicTest, gr_geodesic_test);
        MAKE_SAMPLE(Sample::GRBondi, gr_bondi);
        MAKE_SAMPLE(Sample::NSMerger2D, ns_merger_2d);
        MAKE_SAMPLE(Sample::BNSCocoon1D, bns_cocoon_1d);
        MAKE_SAMPLE(Sample::BNSCocoon2D, bns_cocoon_2d);
        MAKE_SAMPLE(Sample::IsothermalSlab, isothermal_slab);
        MAKE_SAMPLE(Sample::PolytropicSlab, polytropic_slab);
        MAKE_SAMPLE(Sample::SinusoidalPerturbation, sinusoidal_perturbation);
        MAKE_SAMPLE(Sample::JeansInstability, jeans_instability);
        MAKE_SAMPLE(Sample::BonnorEbertKI2000, bonnor_ebert_ki2000);
        MAKE_SAMPLE(Sample::LaneEmdenKI2000, lane_emden_ki2000);
        MAKE_SAMPLE(Sample::IsothermalBonnorEbert, isothermal_bonnor_ebert);
        MAKE_SAMPLE(Sample::HVCCIsothermal10K, hvcc_isothermal_10k);
        MAKE_SAMPLE(Sample::MHDShockTube1, mhd_shock_tube_1);
        MAKE_SAMPLE(Sample::MHDShockTube2, mhd_shock_tube_2);
        MAKE_SAMPLE(Sample::MHDOrszagTang, mhd_orszag_tang);
        MAKE_SAMPLE(Sample::SRMHDBalsara1, srmhd_balsara_1);
        MAKE_SAMPLE(Sample::UniformCloud, uniform_cloud);
        case Sample::DoNotUse:

            // サンプルを使わない場合はここを実装
            std::cout << "Sample::DoNotUse" << std::endl;
            std::cout.flush();
            break;
        default:
            THROW_ERROR("unknown sample type.");
#undef MAKE_SAMPLE
    }
    std::cout << "make_initial_condition() completed" << std::endl;
    std::cout.flush();
}

void Solver::compute_total_energies(real& kinetic, real& thermal, real& potential) const
{
    kinetic = 0.0;
    thermal = 0.0;
    potential = 0.0;
    
    const auto& particles = m_sim->get_particles();
    const int num = m_sim->get_particle_num();
    const bool use_gravity = m_param->gravity.is_valid;
    
    // For SR/GR-GSPH, compute relativistic energies from canonical variables
    if(m_param->type == SPHType::SRGSPH || m_param->type == SPHType::GRGSPH) {
        // In special relativistic Godunov SPH (Kitajima et al. 2025):
        // - Total energy = sum of canonical energy contributions
        // - Canonical energy per baryon: e = γH - P/(Nc²) where H = 1 + ε + P/(ρc²)
        // - Kinetic-like: related to momentum S = γHv
        // - Thermal-like: rest mass + internal energy contribution
        //
        // For output, we decompose the relativistic energy into:
        // - kinetic:   sum_i m_i * (γ_i - 1) * c²  (relativistic kinetic energy)
        // - thermal:   sum_i m_i * ε_i             (internal energy)
        // - potential: gravitational (if enabled)
        //
        // Total relativistic energy: E_total = sum_i m_i * e_i * c²
        // where e_i is the canonical energy per baryon stored in particle
        
        const real c_speed = m_param->srgsph.c_speed;
        const real c2 = c_speed * c_speed;
        
        real rel_kinetic = 0.0;
        real rel_thermal = 0.0;
        real rel_potential = 0.0;
        
        #pragma omp parallel for reduction(+:rel_kinetic,rel_thermal,rel_potential)
        for(int i = 0; i < num; ++i) {
            const auto& p = particles[i];
            
            // Relativistic kinetic energy: m * (γ - 1) * c²
            // γ = Lorentz factor = 1/sqrt(1 - v²/c²)
            real vsq = inner_product(p.vel, p.vel);
            real gamma_lor = 1.0 / std::sqrt(1.0 - vsq / c2);
            rel_kinetic += p.mass * (gamma_lor - 1.0) * c2;
            
            // Internal (thermal) energy: m * ε (specific internal energy)
            rel_thermal += p.mass * p.ene;
            
            // Gravitational potential energy
            if(use_gravity) {
                rel_potential += 0.5 * p.mass * p.phi;
            }
        }
        
        // Add external BH potential if enabled
        if(m_use_external_bh) {
            #pragma omp parallel for reduction(+:rel_potential)
            for(int i = 0; i < num; ++i) {
                const auto& p = particles[i];
                real phi_bh = m_external_bh->potential(p.pos);
                rel_potential += p.mass * phi_bh;  // No 0.5 factor for external potential
            }
        }
        
        kinetic = rel_kinetic;
        thermal = rel_thermal;
        potential = rel_potential;
        
    } else {
        // Standard Newtonian SPH energy computation
        #pragma omp parallel for reduction(+:kinetic,thermal,potential)
        for(int i = 0; i < num; ++i) {
            const auto& p = particles[i];
            
            // Kinetic energy: 0.5 * m * v^2
            real vsq = inner_product(p.vel, p.vel);
            kinetic += 0.5 * p.mass * vsq;
            
            // Thermal energy: m * u
            thermal += p.mass * p.ene;
            
            // Gravitational potential energy: 0.5 * m * phi
            if(use_gravity) {
                potential += 0.5 * p.mass * p.phi;
            }
        }
        
        // Add external BH potential if enabled
        if(m_use_external_bh) {
            #pragma omp parallel for reduction(+:potential)
            for(int i = 0; i < num; ++i) {
                const auto& p = particles[i];
                real phi_bh = m_external_bh->potential(p.pos);
                potential += p.mass * phi_bh;  // No 0.5 factor for external potential
            }
        }
    }
}

}