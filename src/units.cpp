#include "units.hpp"
#include <stdexcept>
#include <sstream>

namespace sph {

// Static constexpr member definitions (required for C++14 ODR-use)
constexpr real UnitSystem::SPEED_OF_LIGHT_CGS;
constexpr real UnitSystem::GRAVITATIONAL_CONSTANT_CGS;
constexpr real UnitSystem::KPC_TO_CM;
constexpr real UnitSystem::MSUN_TO_G;
constexpr real UnitSystem::KM_TO_CM;
constexpr real UnitSystem::MYR_TO_S;
constexpr real UnitSystem::PC_TO_CM;
constexpr real UnitSystem::FM_TO_CM;
constexpr real UnitSystem::PROTON_MASS_G;

// Default constructor - CODE units (no conversion)
UnitSystem::UnitSystem() : m_type(Type::CODE) {
    m_length_to_cgs = 1.0;
    m_mass_to_cgs = 1.0;
    m_time_to_cgs = 1.0;
    m_c_code = SPEED_OF_LIGHT_CGS;  // c in code units (no dimensionless conversion)
    
    // Default labels
    m_length_label = "code_length";
    m_mass_label = "code_mass";
    m_time_label = "code_time";
    m_velocity_label = "code_velocity";
    m_energy_label = "code_energy";
    m_density_label = "code_density";
    m_pressure_label = "code_pressure";
    
    compute_derived_factors();
}

// Constructor with explicit unit specification
UnitSystem::UnitSystem(Type type, real code_length, real code_mass, real code_velocity)
    : m_type(type) {
    
    if (type == Type::CODE) {
        // Code units: no conversion
        m_length_to_cgs = 1.0;
        m_mass_to_cgs = 1.0;
        m_time_to_cgs = 1.0;
        m_c_code = SPEED_OF_LIGHT_CGS;
        
        m_length_label = "code_length";
        m_mass_label = "code_mass";
        m_time_label = "code_time";
        m_velocity_label = "code_velocity";
        m_energy_label = "code_energy";
        m_density_label = "code_density";
        m_pressure_label = "code_pressure";
    } else if (type == Type::GALACTIC) {
        // Galactic units: convert to CGS
        m_length_to_cgs = code_length * KPC_TO_CM;
        m_mass_to_cgs = code_mass * MSUN_TO_G;
        // Time derived from length/velocity
        m_time_to_cgs = (code_length * KPC_TO_CM) / (code_velocity * KM_TO_CM);
        m_c_code = SPEED_OF_LIGHT_CGS / (code_velocity * KM_TO_CM);
        
        m_length_label = "kpc";
        m_mass_label = "M_sun";
        m_time_label = "Myr";
        m_velocity_label = "km/s";
        m_energy_label = "erg";
        m_density_label = "M_sun/kpc³";
        m_pressure_label = "erg/cm³";
    } else if (type == Type::SI) {
        // SI units: convert to CGS
        m_length_to_cgs = code_length * 100.0;  // m to cm
        m_mass_to_cgs = code_mass * 1000.0;      // kg to g
        m_time_to_cgs = code_length / code_velocity;  // Consistency: t = L/v
        m_c_code = SPEED_OF_LIGHT_CGS / (code_velocity * 100.0);  // c in m/s units
        
        m_length_label = "m";
        m_mass_label = "kg";
        m_time_label = "s";
        m_velocity_label = "m/s";
        m_energy_label = "J";
        m_density_label = "kg/m³";
        m_pressure_label = "Pa";
    } else if (type == Type::CGS) {
        // Already in CGS
        m_length_to_cgs = code_length;
        m_mass_to_cgs = code_mass;
        m_time_to_cgs = code_length / code_velocity;  // Consistency
        m_c_code = SPEED_OF_LIGHT_CGS / code_velocity;
        
        m_length_label = "cm";
        m_mass_label = "g";
        m_time_label = "s";
        m_velocity_label = "cm/s";
        m_energy_label = "erg";
        m_density_label = "g/cm³";
        m_pressure_label = "dyne/cm²";
    } else if (type == Type::RELATIVISTIC) {
        // Relativistic natural units: c = 1
        // Will be set by factory method
        m_length_to_cgs = code_length;
        m_mass_to_cgs = code_mass;
        m_time_to_cgs = code_length / SPEED_OF_LIGHT_CGS;  // t = L/c
        m_c_code = 1.0;  // c = 1 in natural units
        
        m_length_label = "L";
        m_mass_label = "M";
        m_time_label = "L/c";
        m_velocity_label = "c";
        m_energy_label = "Mc²";
        m_density_label = "n₀";
        m_pressure_label = "P₀";
    }
    
    compute_derived_factors();
}

// Factory methods
UnitSystem UnitSystem::create_galactic(real code_length_kpc, real code_mass_msun, real code_velocity_kms) {
    return UnitSystem(Type::GALACTIC, code_length_kpc, code_mass_msun, code_velocity_kms);
}

UnitSystem UnitSystem::create_si(real code_length_m, real code_mass_kg, real code_velocity_ms) {
    return UnitSystem(Type::SI, code_length_m, code_mass_kg, code_velocity_ms);
}

UnitSystem UnitSystem::create_cgs(real code_length_cm, real code_mass_g, real code_velocity_cms) {
    return UnitSystem(Type::CGS, code_length_cm, code_mass_g, code_velocity_cms);
}

// Relativistic unit system factory methods
UnitSystem UnitSystem::create_relativistic(
    real code_length_cm,
    real code_density_g_cm3,
    const std::string& length_label,
    const std::string& density_label
) {
    UnitSystem units;
    units.m_type = Type::RELATIVISTIC;
    
    // In relativistic units: c = 1 (code velocity = c)
    // Length is specified directly in cm
    units.m_length_to_cgs = code_length_cm;
    
    // Time = Length / c (since c=1 in code, physical time = L_phys / c_phys)
    units.m_time_to_cgs = code_length_cm / SPEED_OF_LIGHT_CGS;
    
    // Velocity is c (code velocity 1 = c physical)
    units.m_velocity_to_cgs = SPEED_OF_LIGHT_CGS;
    
    // Density: specified directly
    units.m_density_to_cgs = code_density_g_cm3;
    
    // Mass derived from density and length: M = ρ * L³
    units.m_mass_to_cgs = code_density_g_cm3 * code_length_cm * code_length_cm * code_length_cm;
    
    // Energy = M * c²
    units.m_energy_to_cgs = units.m_mass_to_cgs * SPEED_OF_LIGHT_CGS * SPEED_OF_LIGHT_CGS;
    
    // Pressure = Energy / Volume = ρ * c²
    units.m_pressure_to_cgs = code_density_g_cm3 * SPEED_OF_LIGHT_CGS * SPEED_OF_LIGHT_CGS;
    
    // Speed of light is 1 in code units
    units.m_c_code = 1.0;
    
    // Set display labels
    units.m_length_label = length_label;
    units.m_density_label = density_label;
    units.m_time_label = length_label + "/c";
    units.m_velocity_label = "c";
    units.m_mass_label = "M";
    units.m_energy_label = "Mc²";
    units.m_pressure_label = "P₀";
    
    return units;
}

UnitSystem UnitSystem::create_sr_test() {
    // Dimensionless units for SR test problems
    // All quantities are dimensionless: c=1, and all conversions are 1
    UnitSystem units;
    units.m_type = Type::RELATIVISTIC;
    
    // No physical scaling - pure dimensionless
    units.m_length_to_cgs = 1.0;
    units.m_mass_to_cgs = 1.0;
    units.m_time_to_cgs = 1.0;
    units.m_velocity_to_cgs = 1.0;
    units.m_energy_to_cgs = 1.0;
    units.m_density_to_cgs = 1.0;
    units.m_pressure_to_cgs = 1.0;
    units.m_c_code = 1.0;
    
    // Display labels for test problems
    units.m_length_label = "L";
    units.m_mass_label = "M";
    units.m_time_label = "L/c";
    units.m_velocity_label = "c";
    units.m_energy_label = "E";
    units.m_density_label = "n₀";
    units.m_pressure_label = "P₀";
    
    return units;
}

UnitSystem UnitSystem::create_neutron_star(real length_km, real density_scale) {
    // Neutron star merger units
    // Typical: L ~ 10 km, ρ ~ 10^14 g/cm³
    real length_cm = length_km * KM_TO_CM;
    
    UnitSystem units = create_relativistic(length_cm, density_scale, "km", "ρ_nuc");
    
    // Override time label to be more intuitive
    units.m_time_label = "ms";  // milliseconds typical for NS mergers
    
    return units;
}

UnitSystem UnitSystem::create_relativistic_jet(real length_pc, real density_scale) {
    // Relativistic jet units
    // Typical: L ~ 1 pc, ρ ~ proton/cm³
    real length_cm = length_pc * PC_TO_CM;
    
    UnitSystem units = create_relativistic(length_cm, density_scale, "pc", "n_ISM");
    
    return units;
}

UnitSystem UnitSystem::create_imbh_encounter(real length_pc, real mass_1e3Msun, real velocity_kms) {
    // IMBH-cloud encounter units
    // Optimized for 10^4-10^6 M_sun IMBH interacting with 10^2-10^5 M_sun molecular clouds
    // Typical scales: L ~ pc, M ~ 10^3 M_sun, V ~ km/s
    //
    // Example: With defaults (1 pc, 1000 M_sun, 1 km/s):
    //   - IMBH M = 10^5 M_sun = 100 code_mass
    //   - Cloud M = 10^4 M_sun = 10 code_mass
    //   - Cloud R = 5 pc = 5 code_length
    //   - Velocity = 10 km/s = 10 code_velocity
    //   - Time = (1 pc) / (1 km/s) = 0.978 kyr
    
    real length_cm = length_pc * PC_TO_CM;
    real mass_g = mass_1e3Msun * 1.0e3 * MSUN_TO_G;  // 10^3 M_sun to grams
    real velocity_cms = velocity_kms * KM_TO_CM;
    
    // Convert pc to kpc for GALACTIC unit system
    real length_kpc = length_pc / 1000.0;  // 1 pc = 0.001 kpc
    
    UnitSystem units(Type::GALACTIC, length_kpc, mass_1e3Msun * 1.0e3, velocity_kms);
    
    // Override labels for clarity
    units.m_length_label = "pc";
    units.m_mass_label = "1000 M_sun";
    units.m_time_label = "kyr";
    units.m_velocity_label = "km/s";
    units.m_energy_label = "erg";
    units.m_density_label = "M_sun/pc³";
    units.m_pressure_label = "erg/cm³";
    
    return units;
}

// Compute derived conversion factors
void UnitSystem::compute_derived_factors() {
    m_velocity_to_cgs = m_length_to_cgs / m_time_to_cgs;
    m_energy_to_cgs = m_mass_to_cgs * m_velocity_to_cgs * m_velocity_to_cgs;
    m_density_to_cgs = m_mass_to_cgs / (m_length_to_cgs * m_length_to_cgs * m_length_to_cgs);
    m_pressure_to_cgs = m_energy_to_cgs / (m_length_to_cgs * m_length_to_cgs * m_length_to_cgs);
}

// Conversion to physical units (CGS)
real UnitSystem::to_physical_length(real code_val) const {
    return code_val * m_length_to_cgs;
}

real UnitSystem::to_physical_mass(real code_val) const {
    return code_val * m_mass_to_cgs;
}

real UnitSystem::to_physical_time(real code_val) const {
    return code_val * m_time_to_cgs;
}

real UnitSystem::to_physical_velocity(real code_val) const {
    return code_val * m_velocity_to_cgs;
}

real UnitSystem::to_physical_energy(real code_val) const {
    return code_val * m_energy_to_cgs;
}

real UnitSystem::to_physical_density(real code_val) const {
    return code_val * m_density_to_cgs;
}

real UnitSystem::to_physical_pressure(real code_val) const {
    return code_val * m_pressure_to_cgs;
}

// Conversion from physical units (CGS)
real UnitSystem::from_physical_length(real phys_val) const {
    return phys_val / m_length_to_cgs;
}

real UnitSystem::from_physical_mass(real phys_val) const {
    return phys_val / m_mass_to_cgs;
}

real UnitSystem::from_physical_time(real phys_val) const {
    return phys_val / m_time_to_cgs;
}

real UnitSystem::from_physical_velocity(real phys_val) const {
    return phys_val / m_velocity_to_cgs;
}

real UnitSystem::from_physical_energy(real phys_val) const {
    return phys_val / m_energy_to_cgs;
}

real UnitSystem::from_physical_density(real phys_val) const {
    return phys_val / m_density_to_cgs;
}

real UnitSystem::from_physical_pressure(real phys_val) const {
    return phys_val / m_pressure_to_cgs;
}

// Type name accessors
std::string UnitSystem::get_type_name() const {
    switch (m_type) {
        case Type::CODE:        return "Code";
        case Type::RELATIVISTIC: return "Relativistic";
        case Type::GALACTIC:    return "Galactic";
        case Type::SI:          return "SI";
        case Type::CGS:         return "CGS";
        default:                return "Unknown";
    }
}

std::string UnitSystem::get_length_unit_name() const {
    if (!m_length_label.empty()) return m_length_label;
    switch (m_type) {
        case Type::CODE:        return "code_length";
        case Type::RELATIVISTIC: return "L";
        case Type::GALACTIC:    return "kpc";
        case Type::SI:          return "m";
        case Type::CGS:         return "cm";
        default:                return "unknown";
    }
}

std::string UnitSystem::get_mass_unit_name() const {
    if (!m_mass_label.empty()) return m_mass_label;
    switch (m_type) {
        case Type::CODE:        return "code_mass";
        case Type::RELATIVISTIC: return "M";
        case Type::GALACTIC:    return "M_sun";
        case Type::SI:          return "kg";
        case Type::CGS:         return "g";
        default:                return "unknown";
    }
}

std::string UnitSystem::get_time_unit_name() const {
    if (!m_time_label.empty()) return m_time_label;
    switch (m_type) {
        case Type::CODE:        return "code_time";
        case Type::RELATIVISTIC: return "L/c";
        case Type::GALACTIC:    return "Myr";
        case Type::SI:          return "s";
        case Type::CGS:         return "s";
        default:                return "unknown";
    }
}

std::string UnitSystem::get_velocity_unit_name() const {
    if (!m_velocity_label.empty()) return m_velocity_label;
    switch (m_type) {
        case Type::CODE:        return "code_velocity";
        case Type::RELATIVISTIC: return "c";
        case Type::GALACTIC:    return "km/s";
        case Type::SI:          return "m/s";
        case Type::CGS:         return "cm/s";
        default:                return "unknown";
    }
}

std::string UnitSystem::get_energy_unit_name() const {
    if (!m_energy_label.empty()) return m_energy_label;
    switch (m_type) {
        case Type::CODE:        return "code_energy";
        case Type::RELATIVISTIC: return "Mc²";
        case Type::GALACTIC:    return "erg";
        case Type::SI:          return "J";
        case Type::CGS:         return "erg";
        default:                return "unknown";
    }
}

std::string UnitSystem::get_density_unit_name() const {
    if (!m_density_label.empty()) return m_density_label;
    switch (m_type) {
        case Type::CODE:        return "code_density";
        case Type::RELATIVISTIC: return "n₀";
        case Type::GALACTIC:    return "M_sun/kpc^3";
        case Type::SI:          return "kg/m^3";
        case Type::CGS:         return "g/cm^3";
        default:                return "unknown";
    }
}

std::string UnitSystem::get_pressure_unit_name() const {
    if (!m_pressure_label.empty()) return m_pressure_label;
    switch (m_type) {
        case Type::CODE:        return "code_pressure";
        case Type::RELATIVISTIC: return "P₀";
        case Type::GALACTIC:    return "erg/cm^3";
        case Type::SI:          return "Pa";
        case Type::CGS:         return "dyne/cm^2";
        default:                return "unknown";
    }
}

// JSON serialization
nlohmann::json UnitSystem::to_json() const {
    nlohmann::json j;
    
    j["type"] = get_type_name();
    j["type_enum"] = static_cast<int>(m_type);
    
    // Speed of light in code units
    j["c_code"] = m_c_code;
    j["is_relativistic"] = (m_type == Type::RELATIVISTIC);
    
    // Unit names (using custom labels if set)
    j["length_unit"] = get_length_unit_name();
    j["mass_unit"] = get_mass_unit_name();
    j["time_unit"] = get_time_unit_name();
    j["velocity_unit"] = get_velocity_unit_name();
    j["energy_unit"] = get_energy_unit_name();
    j["density_unit"] = get_density_unit_name();
    j["pressure_unit"] = get_pressure_unit_name();
    
    // Display labels (custom labels for visualization)
    j["labels"] = {
        {"length", m_length_label},
        {"mass", m_mass_label},
        {"time", m_time_label},
        {"velocity", m_velocity_label},
        {"energy", m_energy_label},
        {"density", m_density_label},
        {"pressure", m_pressure_label}
    };
    
    // Conversion factors to CGS
    j["length_to_cgs"] = m_length_to_cgs;
    j["mass_to_cgs"] = m_mass_to_cgs;
    j["time_to_cgs"] = m_time_to_cgs;
    j["velocity_to_cgs"] = m_velocity_to_cgs;
    j["energy_to_cgs"] = m_energy_to_cgs;
    j["density_to_cgs"] = m_density_to_cgs;
    j["pressure_to_cgs"] = m_pressure_to_cgs;
    
    return j;
}

UnitSystem UnitSystem::from_json(const nlohmann::json& j) {
    UnitSystem units;
    
    // Read type
    int type_enum = j.value("type_enum", 0);
    units.m_type = static_cast<Type>(type_enum);
    
    // Read c_code
    units.m_c_code = j.value("c_code", 1.0);
    
    // Read conversion factors
    units.m_length_to_cgs = j.value("length_to_cgs", 1.0);
    units.m_mass_to_cgs = j.value("mass_to_cgs", 1.0);
    units.m_time_to_cgs = j.value("time_to_cgs", 1.0);
    units.m_velocity_to_cgs = j.value("velocity_to_cgs", 1.0);
    units.m_energy_to_cgs = j.value("energy_to_cgs", 1.0);
    units.m_density_to_cgs = j.value("density_to_cgs", 1.0);
    units.m_pressure_to_cgs = j.value("pressure_to_cgs", 1.0);
    
    // Read display labels
    if (j.contains("labels")) {
        const auto& labels = j["labels"];
        units.m_length_label = labels.value("length", "");
        units.m_mass_label = labels.value("mass", "");
        units.m_time_label = labels.value("time", "");
        units.m_velocity_label = labels.value("velocity", "");
        units.m_energy_label = labels.value("energy", "");
        units.m_density_label = labels.value("density", "");
        units.m_pressure_label = labels.value("pressure", "");
    }
    
    return units;
}

} // namespace sph
