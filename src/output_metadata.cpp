#include "output_metadata.hpp"
#include <iomanip>
#include <sstream>

namespace sph {

// OutputMetadata implementation
OutputMetadata::OutputMetadata()
    : simulation_name("SPH Simulation")
    , timestamp(generate_timestamp())
    , format_version(1)
    , step(0)
    , time_code(0.0)
    , time_physical(0.0)
    , particle_count(0)
    , sph_type(SPHType::SSPH)
    , kernel_type(KernelType::CUBIC_SPLINE)
    , gamma(1.4)
    , gravitational_constant(1.0)
    , neighbor_number(32)
    , use_gravity(false)
    , use_balsara(true)
    , use_time_dependent_av(false)
    , units()
    , is_checkpoint(false)
    , checkpoint_data()
    , kinetic_energy(0.0)
    , thermal_energy(0.0)
    , potential_energy(0.0)
    , total_energy(0.0)
{
}

std::string OutputMetadata::generate_timestamp() {
    std::time_t now = std::time(nullptr);
    std::tm* utc_time = std::gmtime(&now);
    
    std::ostringstream oss;
    oss << std::put_time(utc_time, "%Y-%m-%dT%H:%M:%SZ");
    return oss.str();
}

std::string OutputMetadata::get_sph_type_name() const {
    switch (sph_type) {
        case SPHType::SSPH:  return "Standard SPH";
        case SPHType::DISPH: return "Density Independent SPH";
        case SPHType::GSPH:  return "Godunov SPH";
        default:             return "Unknown";
    }
}

std::string OutputMetadata::get_kernel_type_name() const {
    switch (kernel_type) {
        case KernelType::CUBIC_SPLINE: return "Cubic Spline";
        case KernelType::WENDLAND:     return "Wendland C4";
        default:                       return "Unknown";
    }
}

nlohmann::json OutputMetadata::to_json() const {
    nlohmann::json j;
    
    // Provenance
    j["simulation_name"] = simulation_name;
    j["timestamp"] = timestamp;
    j["format_version"] = format_version;
    
    // Simulation state
    j["step"] = step;
    j["time_code"] = time_code;
    j["time_physical"] = time_physical;
    j["particle_count"] = particle_count;
    
    // Physics parameters
    j["sph_type"] = get_sph_type_name();
    j["sph_type_enum"] = static_cast<int>(sph_type);
    j["kernel_type"] = get_kernel_type_name();
    j["kernel_type_enum"] = static_cast<int>(kernel_type);
    j["gamma"] = gamma;
    j["gravitational_constant"] = gravitational_constant;
    j["neighbor_number"] = neighbor_number;
    j["use_gravity"] = use_gravity;
    j["use_balsara"] = use_balsara;
    j["use_time_dependent_av"] = use_time_dependent_av;
    
    // Unit system
    j["units"] = units.to_json();
    
    // Checkpoint data (if applicable)
    j["is_checkpoint"] = is_checkpoint;
    if (is_checkpoint) {
        nlohmann::json chk;
        chk["version"] = checkpoint_data.version;
        chk["time"] = checkpoint_data.time;
        chk["step"] = checkpoint_data.step;
        chk["particle_num"] = checkpoint_data.particle_num;
        chk["is_relaxation"] = checkpoint_data.is_relaxation;
        
        if (checkpoint_data.is_relaxation) {
            chk["relaxation_step"] = checkpoint_data.relaxation_step;
            chk["relaxation_total_steps"] = checkpoint_data.relaxation_total_steps;
            chk["accumulated_time"] = checkpoint_data.accumulated_time;
            chk["alpha_scaling"] = checkpoint_data.alpha_scaling;
            chk["rho_center"] = checkpoint_data.rho_center;
            chk["K"] = checkpoint_data.K;
            chk["R"] = checkpoint_data.R;
            chk["M_total"] = checkpoint_data.M_total;
            chk["config_hash"] = checkpoint_data.config_hash;
            chk["preset_name"] = checkpoint_data.preset_name;
        }
        
        j["checkpoint"] = chk;
    }
    
    // Energy diagnostics
    j["kinetic_energy"] = kinetic_energy;
    j["thermal_energy"] = thermal_energy;
    j["potential_energy"] = potential_energy;
    j["total_energy"] = total_energy;
    
    return j;
}

OutputMetadata OutputMetadata::from_json(const nlohmann::json& j) {
    OutputMetadata meta;
    
    // Provenance
    meta.simulation_name = j.value("simulation_name", "SPH Simulation");
    meta.timestamp = j.value("timestamp", generate_timestamp());
    meta.format_version = j.value("format_version", 1);
    
    // Simulation state
    meta.step = j.value("step", 0);
    meta.time_code = j.value("time_code", 0.0);
    meta.time_physical = j.value("time_physical", 0.0);
    meta.particle_count = j.value("particle_count", 0);
    
    // Physics parameters
    meta.sph_type = static_cast<SPHType>(j.value("sph_type_enum", 0));
    meta.kernel_type = static_cast<KernelType>(j.value("kernel_type_enum", 0));
    meta.gamma = j.value("gamma", 1.4);
    meta.gravitational_constant = j.value("gravitational_constant", 1.0);
    meta.neighbor_number = j.value("neighbor_number", 32);
    meta.use_gravity = j.value("use_gravity", false);
    meta.use_balsara = j.value("use_balsara", true);
    meta.use_time_dependent_av = j.value("use_time_dependent_av", false);
    
    // Unit system
    if (j.contains("units")) {
        meta.units = UnitSystem::from_json(j["units"]);
    }
    
    // Checkpoint data
    meta.is_checkpoint = j.value("is_checkpoint", false);
    if (meta.is_checkpoint && j.contains("checkpoint")) {
        const auto& chk = j["checkpoint"];
        meta.checkpoint_data.version = chk.value("version", 1);
        meta.checkpoint_data.time = chk.value("time", 0.0);
        meta.checkpoint_data.step = chk.value("step", 0);
        meta.checkpoint_data.particle_num = chk.value("particle_num", 0);
        meta.checkpoint_data.is_relaxation = chk.value("is_relaxation", false);
        
        if (meta.checkpoint_data.is_relaxation) {
            meta.checkpoint_data.relaxation_step = chk.value("relaxation_step", 0);
            meta.checkpoint_data.relaxation_total_steps = chk.value("relaxation_total_steps", 0);
            meta.checkpoint_data.accumulated_time = chk.value("accumulated_time", 0.0);
            meta.checkpoint_data.alpha_scaling = chk.value("alpha_scaling", 1.0);
            meta.checkpoint_data.rho_center = chk.value("rho_center", 0.0);
            meta.checkpoint_data.K = chk.value("K", 0.0);
            meta.checkpoint_data.R = chk.value("R", 0.0);
            meta.checkpoint_data.M_total = chk.value("M_total", 0.0);
            meta.checkpoint_data.config_hash = chk.value("config_hash", "");
            meta.checkpoint_data.preset_name = chk.value("preset_name", "");
        }
    }
    
    // Energy diagnostics
    meta.kinetic_energy = j.value("kinetic_energy", 0.0);
    meta.thermal_energy = j.value("thermal_energy", 0.0);
    meta.potential_energy = j.value("potential_energy", 0.0);
    meta.total_energy = j.value("total_energy", 0.0);
    
    return meta;
}

// OutputConfig implementation
OutputConfig::OutputConfig()
    : formats({"csv", "hdf5"})  // Default to CSV and HDF5
    , enable_energy_file(true)
    , csv_precision(16)
    , hdf5_compression(6)
    , hdf5_single_file_series(false)
    , vtk_binary(true)
    , output_unit_type(UnitSystem::Type::CODE)
{
}

bool OutputConfig::is_format_enabled(const std::string& format) const {
    return std::find(formats.begin(), formats.end(), format) != formats.end();
}

OutputConfig OutputConfig::from_json(const nlohmann::json& j) {
    OutputConfig config;
    
    // Check for "output" section in JSON
    const nlohmann::json& output_section = j.contains("output") ? j["output"] : j;
    
    // New format: array of objects [{"type": "csv", "precision": 16}, ...]
    if (output_section.contains("formats") && output_section["formats"].is_array() 
        && !output_section["formats"].empty() && output_section["formats"][0].is_object()) {
        
        config.formats.clear();  // Clear defaults
        
        for (const auto& format_obj : output_section["formats"]) {
            std::string type = format_obj.value("type", "");
            if (type.empty()) continue;
            
            config.formats.push_back(type);
            
            // Parse format-specific options
            if (type == "csv") {
                config.csv_precision = format_obj.value("precision", config.csv_precision);
            } else if (type == "hdf5") {
                config.hdf5_compression = format_obj.value("compression", config.hdf5_compression);
            } else if (type == "vtk") {
                config.vtk_binary = format_obj.value("binary", config.vtk_binary);
            }
        }
    }
    // Legacy format 1: simple string array ["csv", "hdf5"]
    else if (output_section.contains("formats") && output_section["formats"].is_array()) {
        config.formats.clear();
        const auto& formats_array = output_section["formats"];
        for (const auto& fmt : formats_array) {
            config.formats.push_back(fmt.get<std::string>());
        }
        
        // Read separate option fields
        config.csv_precision = output_section.value("csvPrecision", config.csv_precision);
        config.hdf5_compression = output_section.value("hdf5Compression", config.hdf5_compression);
        config.vtk_binary = output_section.value("vtkBinary", config.vtk_binary);
    }
    // Legacy format 2: individual enable flags
    else {
        config.formats.clear();
        if (output_section.value("enableCSV", true)) {
            config.formats.push_back("csv");
        }
        if (output_section.value("enableHDF5", true)) {
            config.formats.push_back("hdf5");
        }
        if (output_section.value("enableVTK", false)) {
            config.formats.push_back("vtk");
        }
        
        // Read separate option fields
        config.csv_precision = output_section.value("csvPrecision", config.csv_precision);
        config.hdf5_compression = output_section.value("hdf5Compression", config.hdf5_compression);
        config.vtk_binary = output_section.value("vtkBinary", config.vtk_binary);
    }
    
    // Energy file
    config.enable_energy_file = output_section.value("enableEnergyFile", config.enable_energy_file);
    
    // Unit system for output
    if (j.contains("units")) {
        const auto& units = j["units"];
        std::string system = units.value("type", "CODE");
        
        if (system == "GALACTIC" || system == "galactic") {
            config.output_unit_type = UnitSystem::Type::GALACTIC;
        } else if (system == "SI" || system == "si") {
            config.output_unit_type = UnitSystem::Type::SI;
        } else if (system == "CGS" || system == "cgs") {
            config.output_unit_type = UnitSystem::Type::CGS;
        } else {
            config.output_unit_type = UnitSystem::Type::CODE;
        }
    }
    
    return config;
}

nlohmann::json OutputConfig::to_json() const {
    nlohmann::json j;
    
    // Output formats array with embedded configuration
    nlohmann::json formats_array = nlohmann::json::array();
    for (const auto& fmt : formats) {
        nlohmann::json format_obj;
        format_obj["type"] = fmt;
        
        if (fmt == "csv") {
            format_obj["precision"] = csv_precision;
        } else if (fmt == "hdf5") {
            format_obj["compression"] = hdf5_compression;
        } else if (fmt == "vtk") {
            format_obj["binary"] = vtk_binary;
        }
        
        formats_array.push_back(format_obj);
    }
    j["formats"] = formats_array;
    
    // Energy file
    j["enableEnergyFile"] = enable_energy_file;
    
    // Unit system
    std::string unit_type_str;
    switch (output_unit_type) {
        case UnitSystem::Type::CODE: unit_type_str = "CODE"; break;
        case UnitSystem::Type::GALACTIC: unit_type_str = "GALACTIC"; break;
        case UnitSystem::Type::SI: unit_type_str = "SI"; break;
        case UnitSystem::Type::CGS: unit_type_str = "CGS"; break;
    }
    j["outputUnitSystem"] = unit_type_str;
    
    return j;
}

} // namespace sph
