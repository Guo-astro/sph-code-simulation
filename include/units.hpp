#ifndef SPH_UNITS_HPP
#define SPH_UNITS_HPP

#include "defines.hpp"
#include <nlohmann/json.hpp>
#include <string>
#include <cmath>

namespace sph {

/**
 * @brief Unit system management for SPH simulations
 * 
 * Provides conversion between code units and physical units (Galactic, SI, CGS).
 * Follows RAII principles and immutable design - unit system defined at construction.
 */
class UnitSystem {
public:
    enum class Type {
        CODE,           ///< Dimensionless code units (no conversion)
        RELATIVISTIC,   ///< Natural units with c=1, for SR hydrodynamics
        GALACTIC,       ///< kpc, M_sun, km/s (astrophysics)
        SI,             ///< m, kg, s (engineering)
        CGS             ///< cm, g, s (physics)
    };

    // Physical constants in CGS units
    static constexpr real GRAVITATIONAL_CONSTANT_CGS = 6.67430e-8;  // cm^3 g^-1 s^-2
    static constexpr real SPEED_OF_LIGHT_CGS = 2.99792458e10;       // cm/s
    static constexpr real KPC_TO_CM = 3.085677581e21;               // 1 kpc in cm
    static constexpr real MSUN_TO_G = 1.98847e33;                   // 1 M_sun in g
    static constexpr real KM_TO_CM = 1.0e5;                         // 1 km in cm
    static constexpr real MYR_TO_S = 3.15576e13;                    // 1 Myr in s
    static constexpr real PC_TO_CM = 3.085677581e18;                // 1 pc in cm
    static constexpr real FM_TO_CM = 1.0e-13;                       // 1 fm in cm
    static constexpr real PROTON_MASS_G = 1.6726219e-24;            // proton mass in g

    /**
     * @brief Default constructor - creates CODE unit system (no conversion)
     */
    UnitSystem();

    /**
     * @brief Create unit system with specified type and code unit definitions
     * 
     * @param type Unit system type (GALACTIC, SI, CGS, CODE)
     * @param code_length Code unit length in specified system (e.g., 1.0 kpc for Galactic)
     * @param code_mass Code unit mass in specified system (e.g., 1.0 M_sun for Galactic)
     * @param code_velocity Code unit velocity in specified system (e.g., 1.0 km/s for Galactic)
     */
    UnitSystem(Type type, real code_length, real code_mass, real code_velocity);

    /**
     * @brief Create Galactic unit system (kpc, M_sun, km/s)
     * 
     * @param code_length_kpc Code unit length in kpc
     * @param code_mass_msun Code unit mass in M_sun
     * @param code_velocity_kms Code unit velocity in km/s
     * @return UnitSystem configured for Galactic units
     */
    static UnitSystem create_galactic(real code_length_kpc, real code_mass_msun, real code_velocity_kms);

    /**
     * @brief Create SI unit system (m, kg, s)
     * 
     * @param code_length_m Code unit length in meters
     * @param code_mass_kg Code unit mass in kilograms
     * @param code_velocity_ms Code unit velocity in m/s
     * @return UnitSystem configured for SI units
     */
    static UnitSystem create_si(real code_length_m, real code_mass_kg, real code_velocity_ms);

    /**
     * @brief Create CGS unit system (cm, g, s)
     * 
     * @param code_length_cm Code unit length in cm
     * @param code_mass_g Code unit mass in g
     * @param code_velocity_cms Code unit velocity in cm/s
     * @return UnitSystem configured for CGS units
     */
    static UnitSystem create_cgs(real code_length_cm, real code_mass_g, real code_velocity_cms);

    /**
     * @brief Create relativistic unit system with c=1 (natural units)
     * 
     * For special relativistic hydrodynamics where the speed of light c=1 in code units.
     * All velocities are measured in units of c, time in units of length/c.
     * 
     * @param code_length_cm Physical length scale (1 code_length = X cm)
     * @param code_density_g_cm3 Physical density scale (1 code_density = X g/cm³)
     * @param length_label Display label for length (default: "L")
     * @param density_label Display label for density (default: "n₀")
     * @return UnitSystem configured for relativistic natural units
     */
    static UnitSystem create_relativistic(
        real code_length_cm = 1.0,
        real code_density_g_cm3 = 1.0,
        const std::string& length_label = "L",
        const std::string& density_label = "n₀"
    );

    /**
     * @brief Create dimensionless relativistic unit system for test problems
     * 
     * Convenience method for SR shock tube tests where c=1 and all
     * quantities are dimensionless (code units = physical units).
     * 
     * @return UnitSystem configured for dimensionless SR tests
     */
    static UnitSystem create_sr_test();

    /**
     * @brief Create neutron star merger unit system
     * 
     * Typical scales: L ~ 10 km, ρ ~ 10^14 g/cm³, c = 1
     * 
     * @param length_km Code unit length in km (default: 10 km)
     * @param density_scale Density in g/cm³ (default: 10^14 g/cm³)
     * @return UnitSystem configured for NS merger simulations
     */
    static UnitSystem create_neutron_star(
        real length_km = 10.0,
        real density_scale = 1.0e14
    );

    /**
     * @brief Create relativistic jet unit system
     * 
     * Typical scales: L ~ 1 pc, ρ ~ ISM density, c = 1
     * 
     * @param length_pc Code unit length in parsecs (default: 1 pc)
     * @param density_scale Density in g/cm³ (default: proton/cm³ ~ 1.67e-24)
     * @return UnitSystem configured for relativistic jet simulations
     */
    static UnitSystem create_relativistic_jet(
        real length_pc = 1.0,
        real density_scale = PROTON_MASS_G
    );

    /**
     * @brief Create IMBH-cloud encounter unit system
     * 
     * Optimized for IMBH (10^4-10^6 M_sun) - molecular cloud (10^2-10^5 M_sun) interactions.
     * Typical scales: L ~ pc, M ~ 10^3 M_sun, V ~ km/s
     * 
     * @param length_pc Code unit length in parsecs (default: 1 pc)
     * @param mass_1e3Msun Code unit mass in 10^3 M_sun (default: 1.0 → 1000 M_sun)
     * @param velocity_kms Code unit velocity in km/s (default: 1.0)
     * @return UnitSystem configured for IMBH-cloud encounter simulations
     */
    static UnitSystem create_imbh_encounter(
        real length_pc = 1.0,
        real mass_1e3Msun = 1.0,
        real velocity_kms = 1.0
    );

    // Conversion methods: code units → physical units (CGS)
    real to_physical_length(real code_val) const;
    real to_physical_mass(real code_val) const;
    real to_physical_time(real code_val) const;
    real to_physical_velocity(real code_val) const;
    real to_physical_energy(real code_val) const;
    real to_physical_density(real code_val) const;
    real to_physical_pressure(real code_val) const;

    // Conversion methods: physical units (CGS) → code units
    real from_physical_length(real phys_val) const;
    real from_physical_mass(real phys_val) const;
    real from_physical_time(real phys_val) const;
    real from_physical_velocity(real phys_val) const;
    real from_physical_energy(real phys_val) const;
    real from_physical_density(real phys_val) const;
    real from_physical_pressure(real phys_val) const;

    // Accessors
    Type get_type() const { return m_type; }
    std::string get_type_name() const;
    
    std::string get_length_unit_name() const;
    std::string get_mass_unit_name() const;
    std::string get_time_unit_name() const;
    std::string get_velocity_unit_name() const;
    std::string get_energy_unit_name() const;
    std::string get_density_unit_name() const;
    std::string get_pressure_unit_name() const;

    real get_length_to_cgs() const { return m_length_to_cgs; }
    real get_mass_to_cgs() const { return m_mass_to_cgs; }
    real get_time_to_cgs() const { return m_time_to_cgs; }
    real get_velocity_to_cgs() const { return m_velocity_to_cgs; }
    real get_energy_to_cgs() const { return m_energy_to_cgs; }
    real get_density_to_cgs() const { return m_density_to_cgs; }
    real get_pressure_to_cgs() const { return m_pressure_to_cgs; }

    /**
     * @brief Get speed of light in code units
     * 
     * For relativistic simulations, this is 1.0 (natural units).
     * For non-relativistic simulations, this is c in the code velocity units.
     * 
     * @return Speed of light in code units
     */
    real get_c_code() const { return m_c_code; }

    /**
     * @brief Check if this is a relativistic unit system
     * @return true if Type::RELATIVISTIC
     */
    bool is_relativistic() const { return m_type == Type::RELATIVISTIC; }

    // Custom display labels (for relativistic units)
    void set_length_label(const std::string& label) { m_length_label = label; }
    void set_density_label(const std::string& label) { m_density_label = label; }
    void set_time_label(const std::string& label) { m_time_label = label; }
    void set_velocity_label(const std::string& label) { m_velocity_label = label; }
    void set_pressure_label(const std::string& label) { m_pressure_label = label; }

    // JSON serialization
    nlohmann::json to_json() const;
    static UnitSystem from_json(const nlohmann::json& j);

private:
    Type m_type;

    // Conversion factors: code_unit * factor = CGS_unit
    real m_length_to_cgs;
    real m_mass_to_cgs;
    real m_time_to_cgs;
    real m_velocity_to_cgs;   // derived
    real m_energy_to_cgs;     // derived: mass * velocity^2
    real m_density_to_cgs;    // derived: mass / length^3
    real m_pressure_to_cgs;   // derived: energy / length^3

    // Speed of light in code units (1.0 for relativistic, c_cgs/velocity_to_cgs for others)
    real m_c_code;

    // Custom display labels for units
    std::string m_length_label;
    std::string m_mass_label;
    std::string m_time_label;
    std::string m_velocity_label;
    std::string m_energy_label;
    std::string m_density_label;
    std::string m_pressure_label;

    /**
     * @brief Compute derived conversion factors from base units
     * 
     * Must be called after setting m_length_to_cgs, m_mass_to_cgs, m_time_to_cgs
     */
    void compute_derived_factors();
};

} // namespace sph

#endif // SPH_UNITS_HPP
