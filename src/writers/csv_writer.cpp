#include "writers/csv_writer.hpp"
#include <iomanip>
#include <sstream>
#include <cmath>

namespace sph {

CSVWriter::CSVWriter(int precision)
    : m_file()
    , m_precision(precision)
{
}

CSVWriter::~CSVWriter() {
    close();
}

bool CSVWriter::open(const std::string& filepath, const OutputMetadata& metadata) {
    // Close existing file if open
    if (m_file.is_open()) {
        close();
    }
    
    // Open new file
    m_file.open(filepath);
    if (!m_file.is_open()) {
        return false;
    }
    
    // Set precision
    m_file << std::scientific << std::setprecision(m_precision);
    
    // Write header
    write_header(metadata);
    
    // Write column names
    write_column_names();
    
    return true;
}

void CSVWriter::write_header(const OutputMetadata& metadata) {
    m_file << "# SPH Simulation Output - CSV Format v" << metadata.format_version << "\n";
    m_file << "# Timestamp: " << metadata.timestamp << "\n";
    m_file << "# Simulation: " << metadata.simulation_name << "\n";
    m_file << "#\n";
    
    // Unit system information
    m_file << "# === Unit System ===\n";
    m_file << "# Type: " << metadata.units.get_type_name() << "\n";
    m_file << "# Length: " << metadata.units.get_length_unit_name() 
           << " (" << std::scientific << metadata.units.get_length_to_cgs() << " cm)\n";
    m_file << "# Mass: " << metadata.units.get_mass_unit_name() 
           << " (" << metadata.units.get_mass_to_cgs() << " g)\n";
    m_file << "# Time: " << metadata.units.get_time_unit_name() 
           << " (" << metadata.units.get_time_to_cgs() << " s)\n";
    m_file << "# Velocity: " << metadata.units.get_velocity_unit_name() 
           << " (" << metadata.units.get_velocity_to_cgs() << " cm/s)\n";
    m_file << std::defaultfloat;
    
    // Relativistic units additional info
    if (metadata.units.is_relativistic()) {
        m_file << "# c (code units): " << metadata.units.get_c_code() << "\n";
        m_file << "# Density: " << metadata.units.get_density_unit_name()
               << " (" << std::scientific << metadata.units.get_density_to_cgs() << " g/cm³)\n";
        m_file << "# Pressure: " << metadata.units.get_pressure_unit_name()
               << " (" << metadata.units.get_pressure_to_cgs() << " dyne/cm²)\n";
        m_file << std::defaultfloat;
    }
    m_file << "#\n";
    
    // Simulation state
    m_file << "# === Simulation State ===\n";
    m_file << "# Step: " << metadata.step << "\n";
    m_file << "# Time (code): " << metadata.time_code << "\n";
    m_file << "# Time (physical): " << metadata.time_physical << " " 
           << metadata.units.get_time_unit_name() << "\n";
    m_file << "# Particle Count: " << metadata.particle_count << "\n";
    m_file << "#\n";
    
    // Physics parameters
    m_file << "# === Physics Parameters ===\n";
    m_file << "# Gamma: " << metadata.gamma << "\n";
    m_file << "# G: " << metadata.gravitational_constant << "\n";
    m_file << "# Neighbor Number: " << metadata.neighbor_number << "\n";
    m_file << "# SPH Type: " << metadata.get_sph_type_name() << "\n";
    m_file << "# Kernel: " << metadata.get_kernel_type_name() << "\n";
    m_file << "# Gravity: " << (metadata.use_gravity ? "enabled" : "disabled") << "\n";
    m_file << "# Balsara Switch: " << (metadata.use_balsara ? "enabled" : "disabled") << "\n";
    m_file << "# Time-Dependent AV: " << (metadata.use_time_dependent_av ? "enabled" : "disabled") << "\n";
    m_file << "#\n";
    
    // Energy diagnostics
    m_file << "# === Energy ===\n";
    m_file << "# Kinetic: " << metadata.kinetic_energy << " [code units]\n";
    m_file << "# Thermal: " << metadata.thermal_energy << " [code units]\n";
    m_file << "# Potential: " << metadata.potential_energy << " [code units]\n";
    m_file << "# Total: " << metadata.total_energy << " [code units]\n";
    m_file << "#\n";
    
    // Relaxation information (for Lane-Emden)
    if (metadata.is_relaxation) {
        m_file << "# === Lane-Emden Relaxation ===\n";
        m_file << "# Relaxation Step: " << metadata.relaxation_step << " / " 
               << metadata.relaxation_total_steps << "\n";
        m_file << "# Accumulated Time: " << metadata.accumulated_time << "\n";
        m_file << "# Alpha Scaling: " << metadata.alpha_scaling << "\n";
        m_file << "# Central Density: " << metadata.rho_center << "\n";
        m_file << "# Polytropic K: " << metadata.K << "\n";
        m_file << "# Radius: " << metadata.R << "\n";
        m_file << "# Total Mass: " << metadata.M_total << "\n";
        m_file << "#\n";
    }
    
    // Column description
    m_file << "# === Columns ===\n";
    m_file << "# id: Particle ID\n";
#if DIM == 1
    m_file << "# pos_x: Position [code units]\n";
    m_file << "# vel_x: Normal velocity [code units]\n";
    m_file << "# vel_t: Tangent velocity [code units] (for SR tests)\n";
    m_file << "# acc_x: Acceleration [code units]\n";
#elif DIM == 2
    m_file << "# pos_x, pos_y: Position [code units]\n";
    m_file << "# vel_x, vel_y: Velocity [code units]\n";
    m_file << "# acc_x, acc_y: Acceleration [code units]\n";
#else
    m_file << "# pos_x, pos_y, pos_z: Position [code units]\n";
    m_file << "# vel_x, vel_y, vel_z: Velocity [code units]\n";
    m_file << "# acc_x, acc_y, acc_z: Acceleration [code units]\n";
#endif
    m_file << "# mass: Particle mass [code units]\n";
    m_file << "# dens: Density [code units]\n";
    m_file << "# pres: Pressure [code units]\n";
    m_file << "# ene: Specific internal energy [code units]\n";
    m_file << "# sml: Smoothing length [code units]\n";
    m_file << "# sound: Sound speed [code units]\n";
    m_file << "# alpha: Artificial viscosity coefficient\n";
    m_file << "# balsara: Balsara factor\n";
    m_file << "# gradh: grad-h term factor\n";
    m_file << "# phi: Gravitational potential [code units]\n";
#if DIM == 1
    m_file << "# grav_acc_x: Gravitational acceleration [code units]\n";
#elif DIM == 2
    m_file << "# grav_acc_x, grav_acc_y: Gravitational acceleration [code units]\n";
#else
    m_file << "# grav_acc_x, grav_acc_y, grav_acc_z: Gravitational acceleration [code units]\n";
#endif
    m_file << "# neighbor: Number of neighbors\n";
    m_file << "# is_ghost: Ghost particle flag (0=real, 1=ghost)\n";
    m_file << "#\n";
}

void CSVWriter::write_column_names() {
    m_file << "id,";
#if DIM == 1
    m_file << "pos_x,";
    m_file << "vel_x,";
    m_file << "vel_t,";  // Tangent velocity for SR tests
    m_file << "acc_x,";
#elif DIM == 2
    m_file << "pos_x,pos_y,";
    m_file << "vel_x,vel_y,";
    m_file << "acc_x,acc_y,";
#else
    m_file << "pos_x,pos_y,pos_z,";
    m_file << "vel_x,vel_y,vel_z,";
    m_file << "acc_x,acc_y,acc_z,";
#endif
    m_file << "mass,dens,pres,ene,sml,sound,";
#if DIM == 1
    m_file << "alpha,balsara,gradh,phi,grav_acc_x,";
#elif DIM == 2
    m_file << "alpha,balsara,gradh,phi,grav_acc_x,grav_acc_y,";
#else
    m_file << "alpha,balsara,gradh,phi,grav_acc_x,grav_acc_y,grav_acc_z,";
#endif
    // MHD fields (always 3D for B, regardless of DIM)
    m_file << "B_x,B_y,B_z,vy_mhd,vz_mhd,";
    m_file << "neighbor,is_ghost\n";
}

bool CSVWriter::write_particles(const std::vector<SPHParticle*>& particles) {
    if (!m_file.is_open()) {
        return false;
    }
    
    for (const auto* p : particles) {
        // ID
        m_file << p->id << ",";
        
        // Position (dimension-aware to avoid UB)
#if DIM == 1
        m_file << p->pos[0] << ",";
#elif DIM == 2
        m_file << p->pos[0] << "," << p->pos[1] << ",";
#else
        m_file << p->pos[0] << "," << p->pos[1] << "," << p->pos[2] << ",";
#endif
        
        // Velocity (dimension-aware to avoid UB)
#if DIM == 1
        m_file << p->vel[0] << ",";
        m_file << p->vel_t << ",";  // Tangent velocity for SR tests
#elif DIM == 2
        m_file << p->vel[0] << "," << p->vel[1] << ",";
#else
        m_file << p->vel[0] << "," << p->vel[1] << "," << p->vel[2] << ",";
#endif
        
        // Acceleration (dimension-aware to avoid UB)
#if DIM == 1
        m_file << p->acc[0] << ",";
#elif DIM == 2
        m_file << p->acc[0] << "," << p->acc[1] << ",";
#else
        m_file << p->acc[0] << "," << p->acc[1] << "," << p->acc[2] << ",";
#endif
        
        // Scalar fields
        m_file << p->mass << ",";
        m_file << p->dens << ",";
        m_file << p->pres << ",";
        m_file << p->ene << ",";
        m_file << p->sml << ",";
        m_file << p->sound << ",";
        m_file << p->alpha << ",";
        m_file << p->balsara << ",";
        m_file << p->gradh << ",";
        m_file << p->phi << ",";
        
        // Gravitational acceleration (dimension-aware)
#if DIM == 1
        m_file << p->grav_acc[0] << ",";
#elif DIM == 2
        m_file << p->grav_acc[0] << "," << p->grav_acc[1] << ",";
#else
        m_file << p->grav_acc[0] << "," << p->grav_acc[1] << "," << p->grav_acc[2] << ",";
#endif
        
        // MHD fields (B is always vec3_t, regardless of DIM)
        m_file << p->B[0] << "," << p->B[1] << "," << p->B[2] << ",";
        m_file << p->vy_mhd << "," << p->vz_mhd << ",";

        m_file << p->neighbor << ",";
        m_file << (p->is_ghost ? 1 : 0) << "\n";
    }
    
    m_file.flush();
    return true;
}

void CSVWriter::close() {
    if (m_file.is_open()) {
        m_file.close();
    }
}

bool CSVWriter::read_metadata(const std::string& filepath, OutputMetadata& metadata) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        return false;
    }
    
    std::string line;
    int line_count = 0;
    
    // Parse header lines (lines 1-50)
    while (line_count < 50 && std::getline(file, line)) {
        line_count++;
        
        // Skip empty lines
        if (line.empty() || line[0] != '#') {
            continue;
        }
        
        // Extract key-value from "# Key: Value" format
        auto colon_pos = line.find(':');
        if (colon_pos == std::string::npos) {
            continue;
        }
        
        std::string key = line.substr(2, colon_pos - 2); // Skip "# "
        std::string value = line.substr(colon_pos + 2); // Skip ": "
        
        // Trim whitespace
        key.erase(0, key.find_first_not_of(" \t"));
        key.erase(key.find_last_not_of(" \t") + 1);
        value.erase(0, value.find_first_not_of(" \t"));
        value.erase(value.find_last_not_of(" \t") + 1);
        
        // Parse known fields
        if (key == "Timestamp") {
            metadata.timestamp = value;
        } else if (key == "Simulation") {
            metadata.simulation_name = value;
        } else if (key == "Step") {
            metadata.step = std::stoi(value);
        } else if (key == "Time (code)") {
            metadata.time_code = std::stod(value);
        } else if (key == "Particle Count") {
            metadata.particle_count = std::stoi(value);
        } else if (key == "Gamma") {
            metadata.gamma = std::stod(value);
        } else if (key == "G") {
            metadata.gravitational_constant = std::stod(value);
        } else if (key == "Neighbor Number") {
            metadata.neighbor_number = std::stoi(value);
        } else if (key == "Kinetic") {
            // Extract value before " [code units]"
            auto bracket_pos = value.find(" [");
            if (bracket_pos != std::string::npos) {
                metadata.kinetic_energy = std::stod(value.substr(0, bracket_pos));
            }
        } else if (key == "Thermal") {
            auto bracket_pos = value.find(" [");
            if (bracket_pos != std::string::npos) {
                metadata.thermal_energy = std::stod(value.substr(0, bracket_pos));
            }
        } else if (key == "Potential") {
            auto bracket_pos = value.find(" [");
            if (bracket_pos != std::string::npos) {
                metadata.potential_energy = std::stod(value.substr(0, bracket_pos));
            }
        } else if (key == "Total") {
            auto bracket_pos = value.find(" [");
            if (bracket_pos != std::string::npos) {
                metadata.total_energy = std::stod(value.substr(0, bracket_pos));
            }
        } else if (key == "Relaxation Step") {
            // Format: "step / total"
            metadata.is_relaxation = true;
            auto slash_pos = value.find('/');
            if (slash_pos != std::string::npos) {
                metadata.relaxation_step = std::stoi(value.substr(0, slash_pos));
                metadata.relaxation_total_steps = std::stoi(value.substr(slash_pos + 2));
            }
        } else if (key == "Accumulated Time") {
            metadata.accumulated_time = std::stod(value);
        } else if (key == "Alpha Scaling") {
            metadata.alpha_scaling = std::stod(value);
        } else if (key == "Central Density") {
            metadata.rho_center = std::stod(value);
        } else if (key == "Polytropic K") {
            metadata.K = std::stod(value);
        } else if (key == "Radius") {
            metadata.R = std::stod(value);
        } else if (key == "Total Mass") {
            metadata.M_total = std::stod(value);
        }
    }
    
    file.close();
    return true;
}

// Helper to parse double, treating out-of-range denormals as zero
static double safe_stod(const std::string& s, double default_val = 0.0) {
    try {
        double val = std::stod(s);
        // Check for denormals (very small values that may be uninitialized memory)
        if (std::abs(val) < 1e-300 && val != 0.0) {
            return default_val;
        }
        return val;
    } catch (const std::out_of_range&) {
        // Denormal or out of range - treat as default
        return default_val;
    } catch (const std::invalid_argument&) {
        // Invalid input - treat as default
        return default_val;
    }
}

bool CSVWriter::read_particles(const std::string& filepath, std::vector<SPHParticle*>& particles) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        return false;
    }
    
    // Skip header lines until we find the column names line (starts with "id,")
    std::string line;
    while (std::getline(file, line)) {
        // Check if this is the column names line
        if (line.rfind("id,", 0) == 0) {
            // Found the column header, next line is the first data line
            break;
        }
    }
    
    // Check if we found the header
    if (!file.good()) {
        file.close();
        return false;
    }
    
    // Clear existing particles
    particles.clear();
    
    // Read particle data line by line
    while (std::getline(file, line)) {
        if (line.empty()) {
            continue;
        }
        
        std::istringstream ss(line);
        std::string field;
        std::vector<std::string> fields;
        
        // Split by comma
        while (std::getline(ss, field, ',')) {
            fields.push_back(field);
        }
        
        // Verify we have expected number of columns
        // DIM=3: 25 columns (id, pos*3, vel*3, acc*3, mass, dens, pres, ene, sml, sound, alpha, balsara, gradh, phi, grav_acc*3, neighbor, is_ghost)
        // DIM=2: 21 columns, DIM=1: 18 columns (with vel_t)
        // Also accept old format without grav_acc (22 columns for DIM=3)
#if DIM == 3
        const size_t expected_cols_new = 25;
        const size_t expected_cols_old = 22;
#elif DIM == 2
        const size_t expected_cols_new = 21;
        const size_t expected_cols_old = 18;
#else
        const size_t expected_cols_new = 18;
        const size_t expected_cols_old = 15;
#endif
        const bool has_grav_acc = (fields.size() == expected_cols_new);
        if (fields.size() != expected_cols_new && fields.size() != expected_cols_old) {
            // Clean up allocated particles
            for (auto* p : particles) {
                delete p;
            }
            particles.clear();
            file.close();
            return false;
        }
        
        // Allocate new particle
        auto* p = new SPHParticle();
        
        try {
            // Parse fields in order matching write_column_names()
            int idx = 0;
            p->id = std::stoi(fields[idx++]);
            
            // Position (dimension-aware)
#if DIM == 1
            p->pos[0] = std::stod(fields[idx++]);
#elif DIM == 2
            p->pos[0] = std::stod(fields[idx++]);
            p->pos[1] = std::stod(fields[idx++]);
#else
            p->pos[0] = std::stod(fields[idx++]);
            p->pos[1] = std::stod(fields[idx++]);
            p->pos[2] = std::stod(fields[idx++]);
#endif
            
            // Velocity (dimension-aware)
#if DIM == 1
            p->vel[0] = std::stod(fields[idx++]);
            p->vel_t = std::stod(fields[idx++]);  // Tangent velocity for 1D SR tests
#elif DIM == 2
            p->vel[0] = std::stod(fields[idx++]);
            p->vel[1] = std::stod(fields[idx++]);
#else
            p->vel[0] = std::stod(fields[idx++]);
            p->vel[1] = std::stod(fields[idx++]);
            p->vel[2] = std::stod(fields[idx++]);
#endif
            
            // Acceleration (dimension-aware)
#if DIM == 1
            p->acc[0] = std::stod(fields[idx++]);
#elif DIM == 2
            p->acc[0] = std::stod(fields[idx++]);
            p->acc[1] = std::stod(fields[idx++]);
#else
            p->acc[0] = std::stod(fields[idx++]);
            p->acc[1] = std::stod(fields[idx++]);
            p->acc[2] = std::stod(fields[idx++]);
#endif
            
            // Scalar fields
            p->mass = std::stod(fields[idx++]);
            p->dens = std::stod(fields[idx++]);
            p->pres = std::stod(fields[idx++]);
            p->ene = std::stod(fields[idx++]);
            p->sml = std::stod(fields[idx++]);
            p->sound = std::stod(fields[idx++]);
            p->alpha = safe_stod(fields[idx++], 1.0);    // May be uninitialized
            p->balsara = safe_stod(fields[idx++], 1.0);  // May be uninitialized
            p->gradh = safe_stod(fields[idx++], 1.0);    // May be uninitialized - default to 1.0
            p->phi = safe_stod(fields[idx++], 0.0);      // May be uninitialized
            
            // Gravitational acceleration (optional - new format)
            if (has_grav_acc) {
#if DIM == 1
                p->grav_acc[0] = safe_stod(fields[idx++], 0.0);
#elif DIM == 2
                p->grav_acc[0] = safe_stod(fields[idx++], 0.0);
                p->grav_acc[1] = safe_stod(fields[idx++], 0.0);
#else
                p->grav_acc[0] = safe_stod(fields[idx++], 0.0);
                p->grav_acc[1] = safe_stod(fields[idx++], 0.0);
                p->grav_acc[2] = safe_stod(fields[idx++], 0.0);
#endif
            } else {
                // Old format without grav_acc - initialize to zero
                p->grav_acc = vec_t(0.0);
            }
            
            p->neighbor = std::stoi(fields[idx++]);
            p->is_ghost = (std::stoi(fields[idx++]) != 0);
            
            particles.push_back(p);
            
        } catch (const std::exception&) {
            // Failed to parse this line
            delete p;
            // Clean up all particles
            for (auto* particle : particles) {
                delete particle;
            }
            particles.clear();
            file.close();
            return false;
        }
    }
    
    file.close();
    return !particles.empty();
}

} // namespace sph
