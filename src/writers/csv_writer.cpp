#include "writers/csv_writer.hpp"
#include <iomanip>
#include <sstream>

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
           << " (" << metadata.units.get_length_to_cgs() << " cm)\n";
    m_file << "# Mass: " << metadata.units.get_mass_unit_name() 
           << " (" << metadata.units.get_mass_to_cgs() << " g)\n";
    m_file << "# Time: " << metadata.units.get_time_unit_name() 
           << " (" << metadata.units.get_time_to_cgs() << " s)\n";
    m_file << "# Velocity: " << metadata.units.get_velocity_unit_name() 
           << " (" << metadata.units.get_velocity_to_cgs() << " cm/s)\n";
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
    
    // Checkpoint information (if applicable)
    if (metadata.is_checkpoint && metadata.checkpoint_data.is_relaxation) {
        m_file << "# === Lane-Emden Relaxation ===\n";
        m_file << "# Relaxation Step: " << metadata.checkpoint_data.relaxation_step << " / " 
               << metadata.checkpoint_data.relaxation_total_steps << "\n";
        m_file << "# Accumulated Time: " << metadata.checkpoint_data.accumulated_time << "\n";
        m_file << "# Alpha Scaling: " << metadata.checkpoint_data.alpha_scaling << "\n";
        m_file << "# Central Density: " << metadata.checkpoint_data.rho_center << "\n";
        m_file << "# Polytropic K: " << metadata.checkpoint_data.K << "\n";
        m_file << "# Radius: " << metadata.checkpoint_data.R << "\n";
        m_file << "# Total Mass: " << metadata.checkpoint_data.M_total << "\n";
        if (!metadata.checkpoint_data.preset_name.empty()) {
            m_file << "# Preset: " << metadata.checkpoint_data.preset_name << "\n";
        }
        m_file << "#\n";
    }
    
    // Column description
    m_file << "# === Columns ===\n";
    m_file << "# id: Particle ID\n";
    m_file << "# pos_x, pos_y, pos_z: Position [code units]\n";
    m_file << "# vel_x, vel_y, vel_z: Velocity [code units]\n";
    m_file << "# acc_x, acc_y, acc_z: Acceleration [code units]\n";
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
    m_file << "# neighbor: Number of neighbors\n";
    m_file << "#\n";
}

void CSVWriter::write_column_names() {
    m_file << "id,";
    m_file << "pos_x,pos_y,pos_z,";
    m_file << "vel_x,vel_y,vel_z,";
    m_file << "acc_x,acc_y,acc_z,";
    m_file << "mass,dens,pres,ene,sml,sound,";
    m_file << "alpha,balsara,gradh,phi,neighbor\n";
}

bool CSVWriter::write_particles(const std::vector<SPHParticle*>& particles) {
    if (!m_file.is_open()) {
        return false;
    }
    
    for (const auto* p : particles) {
        // ID
        m_file << p->id << ",";
        
        // Position
        m_file << p->pos[0] << "," << p->pos[1] << "," << p->pos[2] << ",";
        
        // Velocity
        m_file << p->vel[0] << "," << p->vel[1] << "," << p->vel[2] << ",";
        
        // Acceleration
        m_file << p->acc[0] << "," << p->acc[1] << "," << p->acc[2] << ",";
        
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
        m_file << p->neighbor << "\n";
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
            metadata.is_checkpoint = true;
            metadata.checkpoint_data.is_relaxation = true;
            auto slash_pos = value.find('/');
            if (slash_pos != std::string::npos) {
                metadata.checkpoint_data.relaxation_step = std::stoi(value.substr(0, slash_pos));
                metadata.checkpoint_data.relaxation_total_steps = std::stoi(value.substr(slash_pos + 2));
            }
        } else if (key == "Accumulated Time") {
            metadata.checkpoint_data.accumulated_time = std::stod(value);
        } else if (key == "Alpha Scaling") {
            metadata.checkpoint_data.alpha_scaling = std::stod(value);
        } else if (key == "Central Density") {
            metadata.checkpoint_data.rho_center = std::stod(value);
        } else if (key == "Polytropic K") {
            metadata.checkpoint_data.K = std::stod(value);
        } else if (key == "Radius") {
            metadata.checkpoint_data.R = std::stod(value);
        } else if (key == "Total Mass") {
            metadata.checkpoint_data.M_total = std::stod(value);
        } else if (key == "Preset") {
            metadata.checkpoint_data.preset_name = value;
        }
    }
    
    file.close();
    return true;
}

bool CSVWriter::read_particles(const std::string& filepath, std::vector<SPHParticle*>& particles) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        return false;
    }
    
    // Skip 51 header lines (50 metadata + 1 column names)
    std::string line;
    for (int i = 0; i < 51; ++i) {
        if (!std::getline(file, line)) {
            file.close();
            return false;
        }
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
        
        // Verify we have all 21 columns
        if (fields.size() != 21) {
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
            
            // Position
            p->pos[0] = std::stod(fields[idx++]);
            p->pos[1] = std::stod(fields[idx++]);
            p->pos[2] = std::stod(fields[idx++]);
            
            // Velocity
            p->vel[0] = std::stod(fields[idx++]);
            p->vel[1] = std::stod(fields[idx++]);
            p->vel[2] = std::stod(fields[idx++]);
            
            // Acceleration
            p->acc[0] = std::stod(fields[idx++]);
            p->acc[1] = std::stod(fields[idx++]);
            p->acc[2] = std::stod(fields[idx++]);
            
            // Scalar fields
            p->mass = std::stod(fields[idx++]);
            p->dens = std::stod(fields[idx++]);
            p->pres = std::stod(fields[idx++]);
            p->ene = std::stod(fields[idx++]);
            p->sml = std::stod(fields[idx++]);
            p->sound = std::stod(fields[idx++]);
            p->alpha = std::stod(fields[idx++]);
            p->balsara = std::stod(fields[idx++]);
            p->gradh = std::stod(fields[idx++]);
            p->phi = std::stod(fields[idx++]);
            p->neighbor = std::stoi(fields[idx++]);
            
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
