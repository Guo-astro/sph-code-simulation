#include <cassert>

#include <iostream>
#include <iomanip>
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

// relaxation
#include "relaxation/lane_emden_relaxation.hpp"

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
        WRITE_SAMPLE(Sample::Evrard, "Evrard collapse");
        WRITE_SAMPLE(Sample::EvrardColdCollapse, "Evrard cold collapse (demonstrates shock amplification)");
        WRITE_SAMPLE(Sample::LaneEmden, "Lane-Emden hydrostatic");
        WRITE_SAMPLE(Sample::Sedov, "Sedov blast wave");
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
        pt::read_json("sample/shock_tube/shock_tube.json", input);
        m_sample = Sample::ShockTube;
        m_sample_parameters["N"] = input.get<int>("N", 100);
    } else if(name_str == "shock_tube_2d") {
        pt::read_json("sample/shock_tube_2d/shock_tube_2d.json", input);
        m_sample = Sample::ShockTube2D;
        m_sample_parameters["Nx"] = input.get<int>("Nx", 200);
        m_sample_parameters["Ny"] = input.get<int>("Ny", 40);
    } else if(name_str == "vacuum") {
        pt::read_json("sample/vacuum/vacuum.json", input);
        m_sample = Sample::Vacuum;
        m_sample_parameters["N"] = input.get<int>("N", 800);
    } else if(name_str == "strong_shock") {
        pt::read_json("sample/strong_shock/strong_shock.json", input);
        m_sample = Sample::StrongShock;
        m_sample_parameters["N"] = input.get<int>("N", 800);
    } else if(name_str == "gresho_chan_vortex") {
        pt::read_json("sample/gresho_chan_vortex/gresho_chan_vortex.json", input);
        m_sample = Sample::GreshoChanVortex;
        m_sample_parameters["N"] = input.get<int>("N", 64);
    } else if(name_str == "pairing_instability") {
        pt::read_json("sample/pairing_instability/pairing_instability.json", input);
        m_sample = Sample::PairingInstability;
        m_sample_parameters["N"] = input.get<int>("N", 64);
    } else if(name_str == "hydrostatic") {
        pt::read_json("sample/hydrostatic/hydrostatic.json", input);
        m_sample = Sample::HydroStatic;
        m_sample_parameters["N"] = input.get<int>("N", 32);
    } else if(name_str == "khi") {
        pt::read_json("sample/khi/khi.json", input);
        m_sample = Sample::KHI;
        m_sample_parameters["N"] = input.get<int>("N", 128);
    } else if(name_str == "evrard") {
        pt::read_json("sample/evrard/evrard.json", input);
        m_sample = Sample::Evrard;
        m_sample_parameters["N"] = input.get<int>("N", 20);
    } else if(name_str == "evrard_cold_collapse") {
        pt::read_json("sample/evrard_cold_collapse/evrard_cold_collapse.json", input);
        m_sample = Sample::EvrardColdCollapse;
        m_sample_parameters["N"] = input.get<int>("N", 30);
    } else if(name_str == "lane_emden") {
        pt::read_json("sample/lane_emden/lane_emden.json", input);
        m_sample = Sample::LaneEmden;
        m_sample_parameters["N"] = input.get<int>("N", 30);
    } else if(name_str == "sedov") {
        pt::read_json("sample/sedov/sedov.json", input);
        m_sample = Sample::Sedov;
        m_sample_parameters["N"] = input.get<int>("N", 100);
    } else if(name_str == "sr_sod") {
        pt::read_json("sample/sr_sod/sr_sod.json", input);
        m_sample = Sample::SRSod;
        m_sample_parameters["N"] = input.get<int>("N", 50);
        m_sample_parameters["different_nu"] = input.get<bool>("different_nu", false);
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
        } else if (sample_type == "vacuum") {
            m_sample = Sample::Vacuum;
            m_sample_parameters["N"] = input.get<int>("N", 800);
        } else if (sample_type == "strong_shock") {
            m_sample = Sample::StrongShock;
            m_sample_parameters["N"] = input.get<int>("N", 800);
        } else if (sample_type == "sedov") {
            m_sample = Sample::Sedov;
            m_sample_parameters["N"] = input.get<int>("N", 30);
        } else {
            // Try to infer sample type from SPH type and JSON content
            std::string sph_type_check = input.get<std::string>("SPHType", "");
            if(sph_type_check == "srgsph") {
                // Check for SR-specific test names in the path
                if(name_str.find("sr_sod") != std::string::npos || 
                   name_str.find("sod") != std::string::npos) {
                    m_sample = Sample::SRSod;
                    m_sample_parameters["N"] = input.get<int>("N", 50);
                    m_sample_parameters["different_nu"] = input.get<bool>("different_nu", false);
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

    // Kernel
    std::string kernel_name = input.get<std::string>("kernel", "cubic_spline");
    if(kernel_name == "cubic_spline") {
        m_param->kernel = KernelType::CUBIC_SPLINE;
    } else if(kernel_name == "wendland") {
        m_param->kernel = KernelType::WENDLAND;
    } else if(kernel_name == "gaussian") {
        m_param->kernel = KernelType::GAUSSIAN;
    } else {
        THROW_ERROR("kernel is unknown.");
    }

    // smoothing length
    m_param->iterative_sml = input.get<bool>("iterativeSmoothingLength", true);

    // periodic
    m_param->periodic.is_valid = input.get<bool>("periodic", false);
    if(m_param->periodic.is_valid) {
        {
            auto & range_max = input.get_child("rangeMax");
            if(range_max.size() != DIM) {
                THROW_ERROR("rangeMax != DIM");
            }
            int i = 0;
            for(auto & v : range_max) {
                m_param->periodic.range_max[i] = std::stod(v.second.data());
                ++i;
            }
        }

        {
            auto & range_min = input.get_child("rangeMin");
            if(range_min.size() != DIM) {
                THROW_ERROR("rangeMax != DIM");
            }
            int i = 0;
            for(auto & v : range_min) {
                m_param->periodic.range_min[i] = std::stod(v.second.data());
                ++i;
            }
        }
    }

    // gravity
    m_param->gravity.is_valid = input.get<bool>("useGravity", false);
    // Always read G for use in initial conditions (e.g., Lane-Emden)
    m_param->gravity.constant = input.get<real>("G", 1.0);
    if(m_param->gravity.is_valid) {
        m_param->gravity.theta = input.get<real>("theta", 0.5);
    }

    // GSPH
    if(m_param->type == SPHType::GSPH) {
        m_param->gsph.is_2nd_order = input.get<bool>("use2ndOrderGSPH", true);
        
        // Riemann solver selection: "hll" or "iterative"
        std::string riemann_solver_str = input.get<std::string>("riemannSolver", "hll");
        if(riemann_solver_str == "iterative") {
            m_param->gsph.riemann_solver = RiemannSolverType::ITERATIVE;
        } else {
            m_param->gsph.riemann_solver = RiemannSolverType::HLL;
        }
    }
    
    // SRGSPH
    if(m_param->type == SPHType::SRGSPH) {
        m_param->srgsph.is_2nd_order = input.get<bool>("use2ndOrderSRGSPH", true);
        m_param->srgsph.c_speed = input.get<real>("cSpeed", 1.0);
        m_param->srgsph.c_shock = input.get<real>("cShock", 3.0);
        m_param->srgsph.c_cd = input.get<real>("cContactDiscontinuity", 1.0);
        m_param->srgsph.eta = input.get<real>("etaSmoothingLength", 1.0);
        m_param->srgsph.c_smooth = input.get<real>("cSmoothGradient", 2.0);
    }
    
    // Relaxation (for Lane-Emden)
    m_use_relaxation = input.get<bool>("useRelaxation", false);
    m_relaxation_steps = input.get<int>("relaxationSteps", 0);
    m_relaxation_output_freq = input.get<int>("relaxationOutputFreq", 10);
    m_relaxation_only = input.get<bool>("relaxationOnly", false);
    
    if(m_use_relaxation) {
        std::cout << "Relaxation enabled: " << m_relaxation_steps << " steps" << std::endl;
        std::cout << "Relaxation output frequency: every " << m_relaxation_output_freq << " steps" << std::endl;
        if(m_relaxation_only) {
            std::cout << "Relaxation-only mode: Will exit after relaxation without running simulation" << std::endl;
        }
    }
    
    // Resume configuration (SSOT mode)
    m_checkpoint_file = input.get<std::string>("checkpoint.resumeFile", "");
    m_resume_from_checkpoint = !m_checkpoint_file.empty();
    
    if(m_resume_from_checkpoint) {
        std::cout << "Resume mode enabled (SSOT):" << std::endl;
        std::cout << "  - Will resume from: " << m_checkpoint_file << std::endl;
        std::cout << "  - All physics parameters will be loaded from snapshot metadata" << std::endl;
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
    } else {
        std::cerr << "Warning: Unknown unit type '" << unit_type_str << "', using CODE units" << std::endl;
        m_units = UnitSystem();
    }
    
    std::cout << "Unit system: " << m_units.get_type_name() << std::endl;
    
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
    real t_out = m_param->time.output;
    real t_ene = m_param->time.energy;

    // Write initial snapshot
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
            
            // SSOT: Override physics parameters from snapshot metadata
            std::cout << "Applying physics parameters from snapshot (SSOT):" << std::endl;
            m_param->physics.gamma = snapshot_data.gamma;
            std::cout << "  - gamma: " << m_param->physics.gamma << std::endl;
            
            m_param->gravity.constant = snapshot_data.gravitational_constant;
            m_param->gravity.is_valid = snapshot_data.use_gravity;
            std::cout << "  - G: " << m_param->gravity.constant << " (gravity " 
                      << (m_param->gravity.is_valid ? "enabled" : "disabled") << ")" << std::endl;
            
            m_param->type = snapshot_data.sph_type;
            std::cout << "  - SPH type: " << (int)m_param->type << std::endl;
            
            m_param->kernel = snapshot_data.kernel_type;
            std::cout << "  - Kernel: " << (int)m_param->kernel << std::endl;
            
            m_param->physics.neighbor_number = snapshot_data.neighbor_number;
            std::cout << "  - Neighbor number: " << m_param->physics.neighbor_number << std::endl;
            
            m_param->av.use_balsara_switch = snapshot_data.use_balsara;
            m_param->av.use_time_dependent_av = snapshot_data.use_time_dependent_av;
            std::cout << "  - Balsara switch: " << (m_param->av.use_balsara_switch ? "enabled" : "disabled") << std::endl;
            std::cout << "  - Time-dependent AV: " << (m_param->av.use_time_dependent_av ? "enabled" : "disabled") << std::endl;
            
            // Restore Lane-Emden relaxation parameters if this is a relaxation snapshot
            if(snapshot_data.is_relaxation && m_sample == Sample::LaneEmden) {
                m_sample_parameters["alpha"] = snapshot_data.alpha_scaling;
                m_sample_parameters["rho_center"] = snapshot_data.rho_center;
                m_sample_parameters["K"] = snapshot_data.K;
                m_sample_parameters["R"] = snapshot_data.R;
                m_sample_parameters["M_total"] = snapshot_data.M_total;
                std::cout << "  - Restored Lane-Emden parameters from snapshot" << std::endl;
            }
        } else {
            std::cerr << "Warning: Failed to load snapshot, starting from scratch" << std::endl;
        }
    }
    
    // If not resumed, create initial conditions
    if(!resumed) {
        make_initial_condition();
        std::cout << "Initial condition made, particle_num = " << m_sim->get_particle_num() << std::endl;
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
    }
    m_gforce = std::make_shared<GravityForce>();

    // GSPH, GDISPH, and SRGSPH require gradient arrays for MUSCL
    if(m_param->type == SPHType::GSPH || m_param->type == SPHType::GDISPH || m_param->type == SPHType::SRGSPH) {
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
            
#ifndef EXHAUSTIVE_SEARCH
            auto tree = m_sim->get_tree();
            tree->resize(num_p);
            tree->make(particles, num_p);
#endif
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
        
        std::cout << "Progress: [" << std::string(50, ' ') << "] 0%" << std::flush;
        
        int output_counter = 0;
        int last_percent = -1;
        
        for(int step = start_step; step < m_relaxation_steps; ++step) {
            // Rebuild tree for accurate neighbor search
            // Particles move during relaxation, so tree must be updated
            auto& particles = m_sim->get_particles();
            const int num_p = m_sim->get_particle_num();
            
#ifndef EXHAUSTIVE_SEARCH
            auto tree = m_sim->get_tree();
            tree->resize(num_p);
            tree->make(particles, num_p);
#endif
            
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
            
            // For stability during relaxation, use a smaller fraction of CFL timestep
            dt_relax *= 0.1;  // Safety factor for relaxation
            
            // Integrate positions with net acceleration
            // Zero velocities and integrate position directly from acceleration
            auto * periodic = m_sim->get_periodic().get();
            
#pragma omp parallel for
            for(int i = 0; i < num_p; ++i) {
                // Set velocity to zero (constraint)
                particles[i].vel[0] = 0.0;
                particles[i].vel[1] = 0.0;
#if DIM == 3
                particles[i].vel[2] = 0.0;
#endif
                
                // Integrate position using net acceleration: Δx = ½at²
                particles[i].pos[0] += 0.5 * particles[i].acc[0] * dt_relax * dt_relax;
                particles[i].pos[1] += 0.5 * particles[i].acc[1] * dt_relax * dt_relax;
#if DIM == 3
                particles[i].pos[2] += 0.5 * particles[i].acc[2] * dt_relax * dt_relax;
#endif
                
                periodic->apply(particles[i].pos);
            }
            
            // Update accumulated time
            accumulated_time += dt_relax;
            
            // Update progress bar
            int percent = (step * 100) / m_relaxation_steps;
            if(percent != last_percent) {
                int bar_width = 50;
                int filled = (step * bar_width) / m_relaxation_steps;
                
                // Calculate max acceleration for display
                real max_acc = 0.0;
                for(const auto& pi : particles) {
                    real a = std::sqrt(pi.acc[0]*pi.acc[0] + pi.acc[1]*pi.acc[1] + pi.acc[2]*pi.acc[2]);
                    max_acc = std::max(max_acc, a);
                }
                
                std::cout << "\rProgress: [" << std::string(filled, '=') << std::string(bar_width - filled, ' ') 
                          << "] " << percent << "% | Step " << step << "/" << m_relaxation_steps
                          << " | a_max=" << std::fixed << std::setprecision(2) << max_acc
                          << " | dt=" << std::scientific << std::setprecision(2) << dt_relax
                          << " | t=" << std::fixed << std::setprecision(3) << accumulated_time
                          << std::flush;
                last_percent = percent;
            }
            
            // Output snapshots at specified frequency
            if(step % m_relaxation_output_freq == 0 || step == m_relaxation_steps - 1) {
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
                
                // Output snapshot (silently)
                m_output_manager->write_snapshot(m_sim, m_param, m_snapshot_counter++);
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
            
#ifndef EXHAUSTIVE_SEARCH
            auto tree = m_sim->get_tree();
            tree->resize(num);
            tree->make(p, num);
#endif
            
            // Calculate final forces and derivatives
            m_pre->calculation(m_sim);
            m_fforce->calculation(m_sim);
            if(m_param->gravity.is_valid) {
                m_gforce->calculation(m_sim);
            }
            
            // Output relaxed configuration
            m_output_manager->write_snapshot(m_sim, m_param, m_snapshot_counter++);
            
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

#ifndef EXHAUSTIVE_SEARCH
    auto tree = m_sim->get_tree();
    tree->resize(num);
    tree->make(p, num);
#endif

    m_pre->calculation(m_sim);
    m_fforce->calculation(m_sim);
    m_gforce->calculation(m_sim);
}

void Solver::integrate()
{
    m_timestep->calculation(m_sim);

    predict();
#ifndef EXHAUSTIVE_SEARCH
    m_sim->make_tree();
#endif
    m_pre->calculation(m_sim);
    m_fforce->calculation(m_sim);
    m_gforce->calculation(m_sim);
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

    // SR-GSPH uses different integration (conserved variables)
    if(m_param->type == SPHType::SRGSPH) {
        const real c_speed = m_param->srgsph.c_speed;
        
        // DIAGNOSTIC: Print force values at first few steps
        static int step_count = 0;
        if(step_count < 3) {
            std::cerr << "\n=== STEP " << step_count << " FORCE DIAGNOSTIC ===" << std::endl;
            for(int i = 0; i < std::min(5, num); ++i) {
                std::cerr << "Particle " << i << ": S=" << std::abs(p[i].S) 
                          << " dS/dt=" << std::abs(p[i].dS) 
                          << " N=" << p[i].N << " nu=" << p[i].nu
                          << " |dS|*dt=" << std::abs(p[i].dS * dt) << std::endl;
            }
            step_count++;
        }
        
#pragma omp parallel for
        for(int i = 0; i < num; ++i) {
            // For SR-GSPH: Integrate conserved variables S and e
            // Half-step for predictor (k -> k+1/2)
            vec_t S_half = p[i].S + p[i].dS * (0.5 * dt);
            real e_half = p[i].e + p[i].de * (0.5 * dt);
            
            // Full step for position and conserved variables (k -> k+1)
            p[i].S += p[i].dS * dt;
            p[i].e += p[i].de * dt;
            
            // Recover primitive variables at half-step for position update
            auto prim_half = srgsph::PrimitiveRecovery::conserved_to_primitive(
                S_half, e_half, p[i].N, gamma, c_speed
            );
            
            // Update position using half-step velocity (drift)
            p[i].pos += prim_half.vel * dt;
            periodic->apply(p[i].pos);
            
            // Recover primitive variables at full step
            auto prim_full = srgsph::PrimitiveRecovery::conserved_to_primitive(
                p[i].S, p[i].e, p[i].N, gamma, c_speed
            );
            
            // Store primitive variables for output and next step
            p[i].vel = prim_full.vel;
            p[i].vel_p = prim_half.vel;  // Store half-step velocity
            // Internal energy: u = P/[(γ-1)ρ] for ideal gas EOS
            p[i].ene = prim_full.pressure / ((gamma - 1.0) * prim_full.density);
            p[i].ene_p = prim_half.pressure / ((gamma - 1.0) * prim_half.density);
            p[i].pres = prim_full.pressure;
            p[i].dens = prim_full.density;
            p[i].sound = prim_full.sound_speed;
            p[i].gamma_lor = prim_full.gamma_lor;
            p[i].enthalpy = prim_full.enthalpy;
        }
    } else {
        // Standard SPH integration
#pragma omp parallel for
        for(int i = 0; i < num; ++i) {
            // k -> k+1/2
            p[i].vel_p = p[i].vel + p[i].acc * (0.5 * dt);
            p[i].ene_p = p[i].ene + p[i].dene * (0.5 * dt);

            // k -> k+1
            p[i].pos += p[i].vel_p * dt;
            p[i].vel += p[i].acc * dt;
            p[i].ene += p[i].dene * dt;
            p[i].sound = std::sqrt(c_sound * p[i].ene);

            periodic->apply(p[i].pos);
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

    // SR-GSPH uses different correction (conserved variables)
    if(m_param->type == SPHType::SRGSPH) {
        const real c_speed = m_param->srgsph.c_speed;
        
#pragma omp parallel for
        for(int i = 0; i < num; ++i) {
            // For SR-GSPH: Correct using half-step conserved variables
            // Compute S and e at half-step from predictor
            vec_t S_half = p[i].S - p[i].dS * (0.5 * dt);  // Reverse half-step
            real e_half = p[i].e - p[i].de * (0.5 * dt);
            
            // Now apply corrector using new derivatives
            p[i].S = S_half + p[i].dS * (0.5 * dt);
            p[i].e = e_half + p[i].de * (0.5 * dt);
            
            // Recover primitive variables at corrected state
            auto prim = sph::srgsph::PrimitiveRecovery::conserved_to_primitive(
                p[i].S, p[i].e, p[i].N, gamma, c_speed
            );
            
            // Update primitive variables
            p[i].vel = prim.vel;
            // Internal energy: u = P/[(γ-1)ρ] for ideal gas EOS
            p[i].ene = prim.pressure / ((gamma - 1.0) * prim.density);
            p[i].pres = prim.pressure;
            p[i].dens = prim.density;
            p[i].sound = prim.sound_speed;
            p[i].gamma_lor = prim.gamma_lor;
            p[i].enthalpy = prim.enthalpy;
        }
    } else {
        // Standard SPH correction
#pragma omp parallel for
        for(int i = 0; i < num; ++i) {
            p[i].vel = p[i].vel_p + p[i].acc * (0.5 * dt);
            p[i].ene = p[i].ene_p + p[i].dene * (0.5 * dt);
            p[i].sound = std::sqrt(c_sound * p[i].ene);
        }
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
        MAKE_SAMPLE(Sample::Vacuum, vacuum);
        MAKE_SAMPLE(Sample::StrongShock, strong_shock);
        MAKE_SAMPLE(Sample::PressureEquilibrium, pressure_equilibrium);
        MAKE_SAMPLE(Sample::GreshoChanVortex, gresho_chan_vortex);
        MAKE_SAMPLE(Sample::PairingInstability, pairing_instability);
        MAKE_SAMPLE(Sample::HydroStatic, hydrostatic);
        MAKE_SAMPLE(Sample::KHI, khi);
        MAKE_SAMPLE(Sample::Evrard, evrard);
        MAKE_SAMPLE(Sample::EvrardColdCollapse, evrard_cold_collapse);
        MAKE_SAMPLE(Sample::LaneEmden, lane_emden);
        MAKE_SAMPLE(Sample::Sedov, sedov);
        MAKE_SAMPLE(Sample::SRSod, sr_sod);
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
}

}