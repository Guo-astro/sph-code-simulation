#include "output_manager.hpp"
#include "logger.hpp"
#include <cmath>
#include <sys/stat.h>
#include <errno.h>
#include <iomanip>

namespace sph {

OutputManager::OutputManager(const OutputConfig& config,
                           const UnitSystem& units,
                           const std::string& output_dir)
    : m_config(config)
    , m_units(units)
    , m_output_dir(output_dir)
    , m_csv_writer(nullptr)
    , m_hdf5_writer(nullptr)
    , m_vtk_writer(nullptr)
{
    if (m_config.is_format_enabled("csv")) {
        m_csv_writer = std::make_unique<CSVWriter>(m_config.csv_precision);
    }
    
    if (m_config.is_format_enabled("hdf5")) {
        m_hdf5_writer = std::make_unique<HDF5Writer>(m_config.hdf5_compression);
    }
    
    if (m_config.is_format_enabled("vtk")) {
        m_vtk_writer = std::make_unique<VTKWriter>(m_config.vtk_binary);
    }
}

OutputManager::~OutputManager() {
    if (m_energy_file.is_open()) {
        m_energy_file.close();
    }
}

bool OutputManager::initialize() {
    // Create output directory if it doesn't exist
    struct stat st;
    if (stat(m_output_dir.c_str(), &st) != 0) {
        // Directory doesn't exist, create it
        if (mkdir(m_output_dir.c_str(), 0755) != 0) {
            WRITE_LOG << "Failed to create output directory: " << m_output_dir;
            return false;
        }
        WRITE_LOG << "Created output directory: " << m_output_dir;
    } else if (!S_ISDIR(st.st_mode)) {
        WRITE_LOG << "Output path exists but is not a directory: " << m_output_dir;
        return false;
    }
    
    // Open energy file if energy output is enabled
    if (m_config.enable_energy_file) {
        std::string energy_path = m_output_dir + "/energy.dat";
        m_energy_file.open(energy_path);
        if (!m_energy_file.is_open()) {
            WRITE_LOG << "Failed to open energy file: " << energy_path;
            return false;
        }
        
        // Write energy file header
        m_energy_file << std::scientific << std::setprecision(m_config.csv_precision);
        m_energy_file << "# Energy log for SPH simulation\n";
        m_energy_file << "# Units: " << m_units.get_type_name() << "\n";
        m_energy_file << "# Time unit: " << m_units.get_time_unit_name() << "\n";
        m_energy_file << "# Energy unit: " << m_units.get_energy_unit_name() << "\n";
        m_energy_file << "# Columns: time  kinetic  thermal  potential  total\n";
        m_energy_file << "#\n";
        
        WRITE_LOG << "Opened energy file: " << energy_path;
    }
    
    return true;
}

bool OutputManager::write_snapshot(std::shared_ptr<Simulation> sim,
                                   std::shared_ptr<SPHParameters> params,
                                   int count) {
    // Build metadata
    OutputMetadata metadata = build_metadata(sim, params, count);
    
    // Get particles as pointers
    std::vector<SPHParticle*> particle_ptrs;
    particle_ptrs.reserve(sim->get_particle_num());
    for (size_t i = 0; i < sim->get_particles().size(); ++i) {
        particle_ptrs.push_back(&sim->get_particles()[i]);
    }
    
    bool success = true;
    
    // Write CSV if enabled
    if (m_config.is_format_enabled("csv") && m_csv_writer) {
        std::string csv_path = generate_filename("snapshot", count, ".csv");
        
        if (!m_csv_writer->open(csv_path, metadata)) {
            WRITE_LOG << "Failed to open CSV file: " << csv_path;
            success = false;
        } else {
            if (!m_csv_writer->write_particles(particle_ptrs)) {
                WRITE_LOG << "Failed to write CSV data: " << csv_path;
                success = false;
            } else {
                WRITE_LOG << "Wrote snapshot CSV: " << csv_path;
            }
            m_csv_writer->close();
        }
    }
    
    // Write HDF5 if enabled
    if (m_config.is_format_enabled("hdf5") && m_hdf5_writer) {
        std::string hdf5_path = generate_filename("snapshot", count, ".h5");
        
        if (!m_hdf5_writer->open(hdf5_path, metadata)) {
            WRITE_LOG << "Failed to open HDF5 file: " << hdf5_path;
            success = false;
        } else {
            if (!m_hdf5_writer->write_particles(particle_ptrs)) {
                WRITE_LOG << "Failed to write HDF5 data: " << hdf5_path;
                success = false;
            } else {
                WRITE_LOG << "Wrote snapshot HDF5: " << hdf5_path;
            }
            m_hdf5_writer->close();
        }
    }
    
    // Write VTK if enabled
    if (m_config.is_format_enabled("vtk") && m_vtk_writer) {
        std::string vtk_path = generate_filename("snapshot", count, ".vtk");
        
        if (!m_vtk_writer->open(vtk_path, metadata)) {
            WRITE_LOG << "Failed to open VTK file: " << vtk_path;
            success = false;
        } else {
            if (!m_vtk_writer->write_particles(particle_ptrs)) {
                WRITE_LOG << "Failed to write VTK data: " << vtk_path;
                success = false;
            } else {
                WRITE_LOG << "Wrote snapshot VTK: " << vtk_path;
            }
            m_vtk_writer->close();
        }
    }
    
    return success;
}

bool OutputManager::load_for_resume(const std::string& filepath,
                                    std::shared_ptr<Simulation> sim,
                                    OutputMetadata* output_meta) {
    WRITE_LOG << "Loading snapshot for resume from: " << filepath;
    
    // Determine file format from extension
    std::string ext = filepath.substr(filepath.find_last_of('.'));
    bool is_csv = (ext == ".csv");
    bool is_hdf5 = (ext == ".h5" || ext == ".hdf5");
    bool is_vtk = (ext == ".vtk");
    
    if (!is_csv && !is_hdf5 && !is_vtk) {
        WRITE_LOG << "Unsupported checkpoint file format: " << filepath;
        WRITE_LOG << "Supported formats: .csv, .h5, .hdf5, .vtk";
        return false;
    }
    
    // Read particles from file (reading as vector of pointers)
    std::vector<SPHParticle*> particle_ptrs;
    OutputMetadata metadata;
    
    if (is_hdf5) {
        // HDF5 format
        if (!HDF5Writer::read_particles(filepath, particle_ptrs)) {
            WRITE_LOG << "Failed to read particles from HDF5 checkpoint: " << filepath;
            return false;
        }
        
        if (!HDF5Writer::read_metadata(filepath, metadata)) {
            WRITE_LOG << "Failed to read metadata from HDF5 checkpoint: " << filepath;
            for (auto p : particle_ptrs) delete p;
            return false;
        }
    } else if (is_csv) {
        // CSV format
        if (!CSVWriter::read_particles(filepath, particle_ptrs)) {
            WRITE_LOG << "Failed to read particles from CSV checkpoint: " << filepath;
            return false;
        }
        
        if (!CSVWriter::read_metadata(filepath, metadata)) {
            WRITE_LOG << "Failed to read metadata from CSV checkpoint: " << filepath;
            for (auto p : particle_ptrs) delete p;
            return false;
        }
    } else if (is_vtk) {
        // VTK format - note: VTK typically doesn't store full checkpoint metadata
        // We'll read particles but may need to get metadata from a companion file
        WRITE_LOG << "Warning: VTK format may not contain complete checkpoint metadata";
        WRITE_LOG << "Consider using HDF5 or CSV for checkpoints if you need full metadata";
        
        // For now, VTK resume is not fully supported - suggest using HDF5 or CSV
        WRITE_LOG << "ERROR: VTK checkpoint resume not yet implemented";
        WRITE_LOG << "Please use HDF5 (.h5) or CSV (.csv) checkpoint files for resume";
        return false;
    }
    
    // Return output metadata if requested (for physics parameters - SSOT)
    if (output_meta != nullptr) {
        *output_meta = metadata;
    }
    
    // Convert pointers to vector of values for Simulation
    std::vector<SPHParticle> particles;
    particles.reserve(particle_ptrs.size());
    for (auto p : particle_ptrs) {
        particles.push_back(*p);
        delete p;  // Clean up allocated memory
    }
    
    // Set particles in simulation (now uses vector directly)
    sim->get_particles() = particles;
    sim->get_particle_num() = particles.size();
    sim->get_time() = metadata.time_code;
    
    WRITE_LOG << "Successfully loaded " << particles.size() << " particles from snapshot";
    WRITE_LOG << "Resume from step " << metadata.step << ", time " << metadata.time_code;
    WRITE_LOG << "Snapshot physics (SSOT): gamma=" << metadata.gamma 
              << ", G=" << metadata.gravitational_constant
              << ", SPH type=" << static_cast<int>(metadata.sph_type)
              << ", kernel=" << static_cast<int>(metadata.kernel_type);
    
    return true;
}

bool OutputManager::write_energy(real time, real kinetic, real thermal, real potential) {
    if (!m_config.enable_energy_file || !m_energy_file.is_open()) {
        return false;
    }
    
    real total = kinetic + thermal + potential;
    
    m_energy_file << time << "  "
                 << kinetic << "  "
                 << thermal << "  "
                 << potential << "  "
                 << total << "\n";
    
    // Flush to ensure data is written
    m_energy_file.flush();
    
    return true;
}

OutputMetadata OutputManager::build_metadata(std::shared_ptr<Simulation> sim,
                                             std::shared_ptr<SPHParameters> params,
                                             int step) const {
    OutputMetadata metadata;
    
    // Provenance
    metadata.simulation_name = "SPH Simulation";
    metadata.timestamp = OutputMetadata::generate_timestamp();
    metadata.format_version = 1;
    
    // Simulation state
    metadata.step = step;
    metadata.time_code = sim->get_time();
    metadata.time_physical = m_units.to_physical_time(sim->get_time());
    metadata.particle_count = sim->get_particle_num();
    
    // Physics parameters from SPHParameters
    metadata.sph_type = params->type;
    metadata.kernel_type = params->kernel;
    metadata.gamma = params->physics.gamma;
    metadata.gravitational_constant = params->gravity.constant;
    metadata.neighbor_number = params->physics.neighbor_number;
    metadata.use_gravity = params->gravity.is_valid;
    metadata.use_balsara = params->av.use_balsara_switch;
    metadata.use_time_dependent_av = params->av.use_time_dependent_av;
    
    // Units
    metadata.units = m_units;
    
    // Compute energies - get particles as pointers
    std::vector<SPHParticle*> particle_ptrs;
    particle_ptrs.reserve(sim->get_particle_num());
    for (size_t i = 0; i < sim->get_particles().size(); ++i) {
        particle_ptrs.push_back(&sim->get_particles()[i]);
    }
    
    compute_energies(particle_ptrs, params->gravity.is_valid,
                    metadata.kinetic_energy,
                    metadata.thermal_energy,
                    metadata.potential_energy);
    
    metadata.total_energy = metadata.kinetic_energy +
                           metadata.thermal_energy +
                           metadata.potential_energy;
    
    return metadata;
}

void OutputManager::compute_energies(const std::vector<SPHParticle*>& particles,
                                     bool use_gravity,
                                     real& kinetic,
                                     real& thermal,
                                     real& potential) const {
    kinetic = 0.0;
    thermal = 0.0;
    potential = 0.0;
    
    #pragma omp parallel for reduction(+:kinetic,thermal,potential)
    for (size_t i = 0; i < particles.size(); ++i) {
        const SPHParticle* p = particles[i];
        
        // Kinetic energy: 0.5 * m * v^2
        real vsq = inner_product(p->vel, p->vel);
        kinetic += 0.5 * p->mass * vsq;
        
        // Thermal energy: m * u
        thermal += p->mass * p->ene;
        
        // Gravitational potential energy: 0.5 * m * phi
        if (use_gravity) {
            potential += 0.5 * p->mass * p->phi;
        }
    }
}

std::string OutputManager::generate_filename(const std::string& prefix,
                                            int count,
                                            const std::string& extension) const {
    // Format: output_dir/prefix_0000.ext
    std::ostringstream oss;
    oss << m_output_dir << "/" << prefix << "_"
        << std::setfill('0') << std::setw(4) << count
        << extension;
    return oss.str();
}

} // namespace sph
