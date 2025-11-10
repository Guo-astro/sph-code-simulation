#include "output_manager.hpp"
#include "logger.hpp"
#include <cmath>
#include <sys/stat.h>
#include <errno.h>
#include <iomanip>

namespace sph {

OutputManager::OutputManager(const OutputConfig& config,
                           const UnitSystem& units,
                           const std::string& output_dir,
                           const std::string& checkpoint_dir)
    : m_config(config)
    , m_units(units)
    , m_output_dir(output_dir)
    , m_checkpoint_dir(checkpoint_dir.empty() ? output_dir + "/checkpoints" : checkpoint_dir)
    , m_csv_writer(nullptr)
    , m_hdf5_writer(nullptr)
{
    if (m_config.enable_csv) {
        m_csv_writer = std::make_unique<CSVWriter>(m_config.csv_precision);
    }
    
    if (m_config.enable_hdf5) {
        m_hdf5_writer = std::make_unique<HDF5Writer>(m_config.hdf5_compression);
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
    
    // Create checkpoint directory if different from output directory
    if (m_checkpoint_dir != m_output_dir) {
        if (stat(m_checkpoint_dir.c_str(), &st) != 0) {
            // Directory doesn't exist, create it
            if (mkdir(m_checkpoint_dir.c_str(), 0755) != 0) {
                WRITE_LOG << "Failed to create checkpoint directory: " << m_checkpoint_dir;
                return false;
            }
            WRITE_LOG << "Created checkpoint directory: " << m_checkpoint_dir;
        } else if (!S_ISDIR(st.st_mode)) {
            WRITE_LOG << "Checkpoint path exists but is not a directory: " << m_checkpoint_dir;
            return false;
        }
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
    OutputMetadata metadata = build_metadata(sim, params, count, false);
    
    // Get particles as pointers
    std::vector<SPHParticle*> particle_ptrs;
    particle_ptrs.reserve(sim->get_particle_num());
    for (size_t i = 0; i < sim->get_particles().size(); ++i) {
        particle_ptrs.push_back(&sim->get_particles()[i]);
    }
    
    bool success = true;
    
    // Write CSV if enabled
    if (m_config.enable_csv && m_csv_writer) {
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
    if (m_config.enable_hdf5 && m_hdf5_writer) {
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
    
    return success;
}

bool OutputManager::write_checkpoint(std::shared_ptr<Simulation> sim,
                                     std::shared_ptr<SPHParameters> params,
                                     const CheckpointMetadata& checkpoint_meta,
                                     int step) {
    // Build metadata
    OutputMetadata metadata = build_metadata(sim, params, step, true);
    
    // Add checkpoint-specific data
    metadata.checkpoint_data = checkpoint_meta;
    metadata.is_checkpoint = true;
    
    // Get particles as pointers
    std::vector<SPHParticle*> particle_ptrs;
    particle_ptrs.reserve(sim->get_particle_num());
    for (size_t i = 0; i < sim->get_particles().size(); ++i) {
        particle_ptrs.push_back(&sim->get_particles()[i]);
    }
    
    bool success = true;
    
    // Write CSV if enabled (use checkpoint directory)
    if (m_config.enable_csv && m_csv_writer) {
        std::string csv_path = generate_checkpoint_filename("checkpoint", step, ".csv");
        
        if (!m_csv_writer->open(csv_path, metadata)) {
            WRITE_LOG << "Failed to open checkpoint CSV: " << csv_path;
            success = false;
        } else {
            if (!m_csv_writer->write_particles(particle_ptrs)) {
                WRITE_LOG << "Failed to write checkpoint CSV data: " << csv_path;
                success = false;
            } else {
                WRITE_LOG << "Wrote checkpoint CSV: " << csv_path;
            }
            m_csv_writer->close();
        }
    }
    
    // Write HDF5 if enabled (HDF5 is primary for checkpoints, use checkpoint directory)
    if (m_config.enable_hdf5 && m_hdf5_writer) {
        std::string hdf5_path = generate_checkpoint_filename("checkpoint", step, ".h5");
        
        if (!m_hdf5_writer->open(hdf5_path, metadata)) {
            WRITE_LOG << "Failed to open checkpoint HDF5: " << hdf5_path;
            success = false;
        } else {
            if (!m_hdf5_writer->write_particles(particle_ptrs)) {
                WRITE_LOG << "Failed to write checkpoint HDF5 data: " << hdf5_path;
                success = false;
            } else {
                WRITE_LOG << "Wrote checkpoint HDF5: " << hdf5_path;
            }
            m_hdf5_writer->close();
        }
    }
    
    return success;
}

bool OutputManager::load_for_resume(const std::string& filepath,
                                    std::shared_ptr<Simulation> sim,
                                    CheckpointMetadata& checkpoint_meta) {
    WRITE_LOG << "Loading checkpoint from: " << filepath;
    
    // Read particles from HDF5 (reading as vector of pointers)
    std::vector<SPHParticle*> particle_ptrs;
    
    if (!HDF5Writer::read_particles(filepath, particle_ptrs)) {
        WRITE_LOG << "Failed to read particles from checkpoint: " << filepath;
        return false;
    }
    
    // Read metadata from HDF5
    OutputMetadata metadata;
    if (!HDF5Writer::read_metadata(filepath, metadata)) {
        WRITE_LOG << "Failed to read metadata from checkpoint: " << filepath;
        // Clean up allocated particles
        for (auto p : particle_ptrs) delete p;
        return false;
    }
    
    // Verify this is a checkpoint file
    if (!metadata.is_checkpoint) {
        WRITE_LOG << "File is not a checkpoint: " << filepath;
        for (auto p : particle_ptrs) delete p;
        return false;
    }
    
    // Populate checkpoint metadata
    checkpoint_meta = metadata.checkpoint_data;
    
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
    
    WRITE_LOG << "Successfully loaded " << particles.size() << " particles from checkpoint";
    WRITE_LOG << "Resume from step " << checkpoint_meta.step << ", time " << checkpoint_meta.time;
    
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
                                             int step,
                                             bool is_checkpoint) const {
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
    
    // Checkpoint flag
    metadata.is_checkpoint = is_checkpoint;
    
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

std::string OutputManager::generate_checkpoint_filename(const std::string& prefix,
                                                       int count,
                                                       const std::string& extension) const {
    // Format: checkpoint_dir/prefix_0000.ext
    std::ostringstream oss;
    oss << m_checkpoint_dir << "/" << prefix << "_"
        << std::setfill('0') << std::setw(4) << count
        << extension;
    return oss.str();
}

} // namespace sph
