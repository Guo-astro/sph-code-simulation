#include "units.hpp"
#include <stdexcept>
#include <sstream>

namespace sph {

// Default constructor - CODE units (no conversion)
UnitSystem::UnitSystem() : m_type(Type::CODE) {
    m_length_to_cgs = 1.0;
    m_mass_to_cgs = 1.0;
    m_time_to_cgs = 1.0;
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
    } else if (type == Type::GALACTIC) {
        // Galactic units: convert to CGS
        m_length_to_cgs = code_length * KPC_TO_CM;
        m_mass_to_cgs = code_mass * MSUN_TO_G;
        // Time derived from length/velocity
        m_time_to_cgs = (code_length * KPC_TO_CM) / (code_velocity * KM_TO_CM);
    } else if (type == Type::SI) {
        // SI units: convert to CGS
        m_length_to_cgs = code_length * 100.0;  // m to cm
        m_mass_to_cgs = code_mass * 1000.0;      // kg to g
        m_time_to_cgs = code_length / code_velocity;  // Consistency: t = L/v
    } else if (type == Type::CGS) {
        // Already in CGS
        m_length_to_cgs = code_length;
        m_mass_to_cgs = code_mass;
        m_time_to_cgs = code_length / code_velocity;  // Consistency
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
        case Type::CODE:     return "Code";
        case Type::GALACTIC: return "Galactic";
        case Type::SI:       return "SI";
        case Type::CGS:      return "CGS";
        default:             return "Unknown";
    }
}

std::string UnitSystem::get_length_unit_name() const {
    switch (m_type) {
        case Type::CODE:     return "code_length";
        case Type::GALACTIC: return "kpc";
        case Type::SI:       return "m";
        case Type::CGS:      return "cm";
        default:             return "unknown";
    }
}

std::string UnitSystem::get_mass_unit_name() const {
    switch (m_type) {
        case Type::CODE:     return "code_mass";
        case Type::GALACTIC: return "M_sun";
        case Type::SI:       return "kg";
        case Type::CGS:      return "g";
        default:             return "unknown";
    }
}

std::string UnitSystem::get_time_unit_name() const {
    switch (m_type) {
        case Type::CODE:     return "code_time";
        case Type::GALACTIC: return "Myr";
        case Type::SI:       return "s";
        case Type::CGS:      return "s";
        default:             return "unknown";
    }
}

std::string UnitSystem::get_velocity_unit_name() const {
    switch (m_type) {
        case Type::CODE:     return "code_velocity";
        case Type::GALACTIC: return "km/s";
        case Type::SI:       return "m/s";
        case Type::CGS:      return "cm/s";
        default:             return "unknown";
    }
}

std::string UnitSystem::get_energy_unit_name() const {
    switch (m_type) {
        case Type::CODE:     return "code_energy";
        case Type::GALACTIC: return "erg";
        case Type::SI:       return "J";
        case Type::CGS:      return "erg";
        default:             return "unknown";
    }
}

std::string UnitSystem::get_density_unit_name() const {
    switch (m_type) {
        case Type::CODE:     return "code_density";
        case Type::GALACTIC: return "M_sun/kpc^3";
        case Type::SI:       return "kg/m^3";
        case Type::CGS:      return "g/cm^3";
        default:             return "unknown";
    }
}

std::string UnitSystem::get_pressure_unit_name() const {
    switch (m_type) {
        case Type::CODE:     return "code_pressure";
        case Type::GALACTIC: return "erg/cm^3";
        case Type::SI:       return "Pa";
        case Type::CGS:      return "dyne/cm^2";
        default:             return "unknown";
    }
}

// JSON serialization
nlohmann::json UnitSystem::to_json() const {
    nlohmann::json j;
    
    j["type"] = get_type_name();
    j["type_enum"] = static_cast<int>(m_type);
    
    // Unit names
    j["length_unit"] = get_length_unit_name();
    j["mass_unit"] = get_mass_unit_name();
    j["time_unit"] = get_time_unit_name();
    j["velocity_unit"] = get_velocity_unit_name();
    j["energy_unit"] = get_energy_unit_name();
    j["density_unit"] = get_density_unit_name();
    j["pressure_unit"] = get_pressure_unit_name();
    
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
    
    // Read conversion factors
    units.m_length_to_cgs = j.value("length_to_cgs", 1.0);
    units.m_mass_to_cgs = j.value("mass_to_cgs", 1.0);
    units.m_time_to_cgs = j.value("time_to_cgs", 1.0);
    units.m_velocity_to_cgs = j.value("velocity_to_cgs", 1.0);
    units.m_energy_to_cgs = j.value("energy_to_cgs", 1.0);
    units.m_density_to_cgs = j.value("density_to_cgs", 1.0);
    units.m_pressure_to_cgs = j.value("pressure_to_cgs", 1.0);
    
    return units;
}

} // namespace sph
