"""
SPH Unit System - Python Bridge

This module provides Python classes that mirror the C++ UnitSystem class,
allowing consistent unit handling between simulation and visualization.

Physical Constants (CGS):
- Speed of light: c = 2.99792458e10 cm/s
- Gravitational constant: G = 6.67430e-8 cm³ g⁻¹ s⁻²
- 1 kpc = 3.085677581e21 cm
- 1 M_sun = 1.98847e33 g
- 1 pc = 3.085677581e18 cm
- Proton mass = 1.6726219e-24 g

Usage:
    # For SR test problems (dimensionless)
    units = RelativisticUnits.create_sr_test()
    
    # For neutron star mergers
    units = RelativisticUnits.create_neutron_star(length_km=10, density=1e14)
    
    # From simulation config
    units = UnitSystem.from_config(config_dict)
"""

from dataclasses import dataclass, field
from typing import Optional, Dict, Any
from enum import Enum
import json


# Physical constants in CGS
SPEED_OF_LIGHT_CGS = 2.99792458e10      # cm/s
GRAVITATIONAL_CONSTANT_CGS = 6.67430e-8 # cm³ g⁻¹ s⁻²
KPC_TO_CM = 3.085677581e21              # 1 kpc in cm
MSUN_TO_G = 1.98847e33                  # 1 M_sun in g
KM_TO_CM = 1.0e5                        # 1 km in cm
MYR_TO_S = 3.15576e13                   # 1 Myr in s
PC_TO_CM = 3.085677581e18               # 1 pc in cm
FM_TO_CM = 1.0e-13                      # 1 fm in cm
PROTON_MASS_G = 1.6726219e-24           # proton mass in g


class UnitType(Enum):
    """Unit system types matching C++ UnitSystem::Type"""
    CODE = 0
    RELATIVISTIC = 1
    GALACTIC = 2
    SI = 3
    CGS = 4


@dataclass
class UnitSystem:
    """
    Python bridge for C++ UnitSystem class.
    
    Provides conversion between code units and physical units (CGS).
    """
    
    # Unit type
    unit_type: UnitType = UnitType.CODE
    
    # Conversion factors: code_unit * factor = CGS_unit
    length_to_cgs: float = 1.0
    mass_to_cgs: float = 1.0
    time_to_cgs: float = 1.0
    velocity_to_cgs: float = 1.0
    energy_to_cgs: float = 1.0
    density_to_cgs: float = 1.0
    pressure_to_cgs: float = 1.0
    
    # Speed of light in code units
    c_code: float = SPEED_OF_LIGHT_CGS
    
    # Display labels
    length_label: str = "code_length"
    mass_label: str = "code_mass"
    time_label: str = "code_time"
    velocity_label: str = "code_velocity"
    energy_label: str = "code_energy"
    density_label: str = "code_density"
    pressure_label: str = "code_pressure"
    
    def __post_init__(self):
        """Compute derived units if not set"""
        if self.velocity_to_cgs == 1.0 and self.unit_type != UnitType.CODE:
            self.velocity_to_cgs = self.length_to_cgs / self.time_to_cgs
        if self.energy_to_cgs == 1.0 and self.unit_type != UnitType.CODE:
            self.energy_to_cgs = self.mass_to_cgs * self.velocity_to_cgs ** 2
        if self.density_to_cgs == 1.0 and self.unit_type != UnitType.CODE:
            self.density_to_cgs = self.mass_to_cgs / self.length_to_cgs ** 3
        if self.pressure_to_cgs == 1.0 and self.unit_type != UnitType.CODE:
            self.pressure_to_cgs = self.energy_to_cgs / self.length_to_cgs ** 3
    
    # === Conversion to physical (CGS) ===
    def to_physical_length(self, code_val: float) -> float:
        return code_val * self.length_to_cgs
    
    def to_physical_mass(self, code_val: float) -> float:
        return code_val * self.mass_to_cgs
    
    def to_physical_time(self, code_val: float) -> float:
        return code_val * self.time_to_cgs
    
    def to_physical_velocity(self, code_val: float) -> float:
        return code_val * self.velocity_to_cgs
    
    def to_physical_density(self, code_val: float) -> float:
        return code_val * self.density_to_cgs
    
    def to_physical_pressure(self, code_val: float) -> float:
        return code_val * self.pressure_to_cgs
    
    # === Conversion from physical (CGS) ===
    def from_physical_length(self, phys_val: float) -> float:
        return phys_val / self.length_to_cgs
    
    def from_physical_velocity(self, phys_val: float) -> float:
        return phys_val / self.velocity_to_cgs
    
    def from_physical_density(self, phys_val: float) -> float:
        return phys_val / self.density_to_cgs
    
    # === Type checking ===
    def is_relativistic(self) -> bool:
        return self.unit_type == UnitType.RELATIVISTIC
    
    def get_type_name(self) -> str:
        return self.unit_type.name.capitalize()
    
    # === Factory methods ===
    @classmethod
    def create_code(cls) -> 'UnitSystem':
        """Create default CODE unit system (no conversion)"""
        return cls(unit_type=UnitType.CODE)
    
    @classmethod
    def create_galactic(cls, length_kpc: float = 1.0, mass_msun: float = 1e10, 
                       velocity_kms: float = 1.0) -> 'UnitSystem':
        """Create Galactic unit system (kpc, M_sun, km/s)"""
        length_cgs = length_kpc * KPC_TO_CM
        mass_cgs = mass_msun * MSUN_TO_G
        velocity_cgs = velocity_kms * KM_TO_CM
        time_cgs = length_cgs / velocity_cgs
        
        return cls(
            unit_type=UnitType.GALACTIC,
            length_to_cgs=length_cgs,
            mass_to_cgs=mass_cgs,
            time_to_cgs=time_cgs,
            velocity_to_cgs=velocity_cgs,
            c_code=SPEED_OF_LIGHT_CGS / velocity_cgs,
            length_label="kpc",
            mass_label="M_sun",
            time_label="Myr",
            velocity_label="km/s",
            density_label="M_sun/kpc³",
            pressure_label="erg/cm³"
        )
    
    @classmethod
    def from_json(cls, j: Dict[str, Any]) -> 'UnitSystem':
        """Load from JSON (matches C++ UnitSystem::from_json)"""
        type_enum = j.get('type_enum', 0)
        
        return cls(
            unit_type=UnitType(type_enum),
            length_to_cgs=j.get('length_to_cgs', 1.0),
            mass_to_cgs=j.get('mass_to_cgs', 1.0),
            time_to_cgs=j.get('time_to_cgs', 1.0),
            velocity_to_cgs=j.get('velocity_to_cgs', 1.0),
            energy_to_cgs=j.get('energy_to_cgs', 1.0),
            density_to_cgs=j.get('density_to_cgs', 1.0),
            pressure_to_cgs=j.get('pressure_to_cgs', 1.0),
            c_code=j.get('c_code', 1.0),
            length_label=j.get('labels', {}).get('length', ''),
            mass_label=j.get('labels', {}).get('mass', ''),
            time_label=j.get('labels', {}).get('time', ''),
            velocity_label=j.get('labels', {}).get('velocity', ''),
            energy_label=j.get('labels', {}).get('energy', ''),
            density_label=j.get('labels', {}).get('density', ''),
            pressure_label=j.get('labels', {}).get('pressure', '')
        )
    
    @classmethod
    def from_config(cls, config: Dict[str, Any]) -> 'UnitSystem':
        """
        Load unit system from simulation config JSON.
        
        Supports:
        - units.type: "CODE", "RELATIVISTIC", "GALACTIC", etc.
        - units.preset: "neutron_star", "relativistic_jet"
        - units.labels.*: Custom display labels
        
        For SR-GSPH configs (SPHType = "srgsph"), defaults to relativistic units.
        """
        units_config = config.get('units', {})
        unit_type_str = units_config.get('type', 'CODE').upper()
        
        # Auto-detect SR-GSPH
        sph_type = config.get('SPHType', '').lower()
        if sph_type == 'srgsph' and unit_type_str == 'CODE':
            unit_type_str = 'RELATIVISTIC'
        
        if unit_type_str in ('RELATIVISTIC', 'SR_TEST'):
            preset = units_config.get('preset', '')
            if preset == 'neutron_star':
                length_km = units_config.get('length_km', 10.0)
                density = units_config.get('density_g_cm3', 1e14)
                return RelativisticUnits.create_neutron_star(length_km, density)
            elif preset == 'relativistic_jet':
                length_pc = units_config.get('length_pc', 1.0)
                density = units_config.get('density_g_cm3', PROTON_MASS_G)
                return RelativisticUnits.create_relativistic_jet(length_pc, density)
            elif units_config.get('length_cm') or units_config.get('density_gcm3'):
                # Custom scales specified directly
                length_cm = units_config.get('length_cm', 1.0)
                density = units_config.get('density_gcm3', 1.0)
                length_label = units_config.get('labels', {}).get('length', 'L')
                density_label = units_config.get('labels', {}).get('density', 'n₀')
                return RelativisticUnits.create_relativistic(
                    length_cm, density, length_label, density_label
                )
            else:
                units = RelativisticUnits.create_sr_test()
                # Apply custom labels
                labels = units_config.get('labels', {})
                if labels.get('length'): units.length_label = labels['length']
                if labels.get('density'): units.density_label = labels['density']
                if labels.get('time'): units.time_label = labels['time']
                if labels.get('velocity'): units.velocity_label = labels['velocity']
                if labels.get('pressure'): units.pressure_label = labels['pressure']
                return units
        elif unit_type_str == 'GALACTIC':
            return cls.create_galactic(
                units_config.get('length_kpc', 1.0),
                units_config.get('mass_msun', 1e10),
                units_config.get('velocity_kms', 1.0)
            )
        else:
            return cls.create_code()
    
    def to_json(self) -> Dict[str, Any]:
        """Convert to JSON (matches C++ UnitSystem::to_json)"""
        return {
            'type': self.get_type_name(),
            'type_enum': self.unit_type.value,
            'c_code': self.c_code,
            'is_relativistic': self.is_relativistic(),
            'length_unit': self.length_label,
            'mass_unit': self.mass_label,
            'time_unit': self.time_label,
            'velocity_unit': self.velocity_label,
            'energy_unit': self.energy_label,
            'density_unit': self.density_label,
            'pressure_unit': self.pressure_label,
            'labels': {
                'length': self.length_label,
                'mass': self.mass_label,
                'time': self.time_label,
                'velocity': self.velocity_label,
                'energy': self.energy_label,
                'density': self.density_label,
                'pressure': self.pressure_label,
            },
            'length_to_cgs': self.length_to_cgs,
            'mass_to_cgs': self.mass_to_cgs,
            'time_to_cgs': self.time_to_cgs,
            'velocity_to_cgs': self.velocity_to_cgs,
            'energy_to_cgs': self.energy_to_cgs,
            'density_to_cgs': self.density_to_cgs,
            'pressure_to_cgs': self.pressure_to_cgs,
        }


@dataclass
class RelativisticUnits(UnitSystem):
    """
    Relativistic unit system with c=1 (natural units).
    
    For special relativistic hydrodynamics simulations.
    All velocities are measured in units of c, time in units of length/c.
    """
    
    def __post_init__(self):
        self.unit_type = UnitType.RELATIVISTIC
        self.c_code = 1.0
        super().__post_init__()
    
    @classmethod
    def create_sr_test(cls) -> 'RelativisticUnits':
        """
        Create dimensionless relativistic unit system for test problems.
        
        All quantities are dimensionless: c=1, and all conversions are 1.
        Used for SR shock tube tests (Sod, blast wave, etc.)
        """
        return cls(
            unit_type=UnitType.RELATIVISTIC,
            length_to_cgs=1.0,
            mass_to_cgs=1.0,
            time_to_cgs=1.0,
            velocity_to_cgs=1.0,
            energy_to_cgs=1.0,
            density_to_cgs=1.0,
            pressure_to_cgs=1.0,
            c_code=1.0,
            length_label="L",
            mass_label="M",
            time_label="L/c",
            velocity_label="c",
            energy_label="E",
            density_label="n₀",
            pressure_label="P₀"
        )
    
    @classmethod
    def create_relativistic(cls, code_length_cm: float = 1.0, 
                           code_density_g_cm3: float = 1.0,
                           length_label: str = "L",
                           density_label: str = "n₀") -> 'RelativisticUnits':
        """
        Create relativistic unit system with specified physical scales.
        
        Args:
            code_length_cm: Physical length scale (1 code_length = X cm)
            code_density_g_cm3: Physical density scale (1 code_density = X g/cm³)
            length_label: Display label for length
            density_label: Display label for density
        """
        # In relativistic units: c = 1 (code velocity = c)
        time_cgs = code_length_cm / SPEED_OF_LIGHT_CGS
        velocity_cgs = SPEED_OF_LIGHT_CGS
        mass_cgs = code_density_g_cm3 * code_length_cm ** 3
        energy_cgs = mass_cgs * SPEED_OF_LIGHT_CGS ** 2
        pressure_cgs = code_density_g_cm3 * SPEED_OF_LIGHT_CGS ** 2
        
        return cls(
            unit_type=UnitType.RELATIVISTIC,
            length_to_cgs=code_length_cm,
            mass_to_cgs=mass_cgs,
            time_to_cgs=time_cgs,
            velocity_to_cgs=velocity_cgs,
            energy_to_cgs=energy_cgs,
            density_to_cgs=code_density_g_cm3,
            pressure_to_cgs=pressure_cgs,
            c_code=1.0,
            length_label=length_label,
            mass_label="M",
            time_label=f"{length_label}/c",
            velocity_label="c",
            energy_label="Mc²",
            density_label=density_label,
            pressure_label="P₀"
        )
    
    @classmethod
    def create_neutron_star(cls, length_km: float = 10.0, 
                           density_scale: float = 1e14) -> 'RelativisticUnits':
        """
        Create neutron star merger unit system.
        
        Typical scales: L ~ 10 km, ρ ~ 10^14 g/cm³, c = 1
        
        Args:
            length_km: Code unit length in km (default: 10 km)
            density_scale: Density in g/cm³ (default: 10^14 g/cm³)
        """
        length_cm = length_km * KM_TO_CM
        units = cls.create_relativistic(length_cm, density_scale, "km", "ρ_nuc")
        units.time_label = "ms"  # Typical for NS mergers
        return units
    
    @classmethod
    def create_relativistic_jet(cls, length_pc: float = 1.0,
                               density_scale: float = PROTON_MASS_G) -> 'RelativisticUnits':
        """
        Create relativistic jet unit system.
        
        Typical scales: L ~ 1 pc, ρ ~ ISM density (proton/cm³), c = 1
        
        Args:
            length_pc: Code unit length in parsecs (default: 1 pc)
            density_scale: Density in g/cm³ (default: proton mass ~ 1.67e-24)
        """
        length_cm = length_pc * PC_TO_CM
        return cls.create_relativistic(length_cm, density_scale, "pc", "n_ISM")
    
    @classmethod
    def create_grb_afterglow(cls, length_cm: float = 1e16,
                            density_scale: float = PROTON_MASS_G) -> 'RelativisticUnits':
        """
        Create GRB afterglow unit system.
        
        Typical scales: L ~ 10^16 cm, ρ ~ ISM/wind density, c = 1
        
        Args:
            length_cm: Code unit length in cm (default: 10^16 cm)
            density_scale: Density in g/cm³ (default: proton mass)
        """
        return cls.create_relativistic(length_cm, density_scale, "10¹⁶ cm", "n_ext")
    
    # === Convenience methods for relativistic quantities ===
    def velocity_in_c(self, code_velocity: float) -> float:
        """Convert code velocity to fraction of c"""
        return code_velocity  # In relativistic units, v is already in units of c
    
    def lorentz_factor(self, code_velocity: float) -> float:
        """Compute Lorentz factor from code velocity"""
        v = self.velocity_in_c(code_velocity)
        return 1.0 / (1.0 - v ** 2) ** 0.5
    
    def format_time(self, t: float) -> str:
        """Format time with proper units"""
        return f"t = {t:.4f} [{self.time_label}]"
    
    def format_velocity(self, v: float) -> str:
        """Format velocity with units (as fraction of c)"""
        return f"{v:.4f} c"


# === Convenience functions ===
def parse_units_from_snapshot(header_lines: list) -> Optional[UnitSystem]:
    """
    Parse unit system from CSV snapshot header comments.
    
    Args:
        header_lines: List of header comment lines starting with '#'
    
    Returns:
        UnitSystem if found, None otherwise
    """
    units_data = {}
    labels = {}
    
    for line in header_lines:
        if not line.startswith('#'):
            continue
        line = line[1:].strip()  # Remove '#' and whitespace
        
        if line.startswith('Type:'):
            type_str = line.split(':', 1)[1].strip()
            if type_str.lower() == 'relativistic':
                units_data['type_enum'] = UnitType.RELATIVISTIC.value
            elif type_str.lower() == 'galactic':
                units_data['type_enum'] = UnitType.GALACTIC.value
            elif type_str.lower() == 'cgs':
                units_data['type_enum'] = UnitType.CGS.value
            else:
                units_data['type_enum'] = UnitType.CODE.value
        
        elif line.startswith('c (code units):'):
            units_data['c_code'] = float(line.split(':', 1)[1].strip())
        
        elif line.startswith('Length:'):
            parts = line.split('(')
            if len(parts) > 1:
                labels['length'] = parts[0].split(':', 1)[1].strip()
                # Extract CGS value
                cgs_str = parts[1].split()[0]
                try:
                    units_data['length_to_cgs'] = float(cgs_str)
                except ValueError:
                    pass
        
        elif line.startswith('Density:'):
            parts = line.split('(')
            if len(parts) > 1:
                labels['density'] = parts[0].split(':', 1)[1].strip()
                cgs_str = parts[1].split()[0]
                try:
                    units_data['density_to_cgs'] = float(cgs_str)
                except ValueError:
                    pass
    
    if not units_data:
        return None
    
    units_data['labels'] = labels
    return UnitSystem.from_json(units_data)


if __name__ == '__main__':
    # Test unit system creation
    print("=== UnitSystem Tests ===\n")
    
    # Test 1: SR test units (dimensionless)
    sr_test = RelativisticUnits.create_sr_test()
    print("SR Test Units:")
    print(f"  Type: {sr_test.get_type_name()}")
    print(f"  c (code): {sr_test.c_code}")
    print(f"  Length: {sr_test.length_label}")
    print(f"  Time: {sr_test.time_label}")
    print(f"  Velocity: {sr_test.velocity_label}")
    print(f"  Density: {sr_test.density_label}")
    print()
    
    # Test 2: Neutron star units
    ns_units = RelativisticUnits.create_neutron_star(10.0, 1e14)
    print("Neutron Star Units:")
    print(f"  Type: {ns_units.get_type_name()}")
    print(f"  Length scale: {ns_units.length_to_cgs / KM_TO_CM:.1f} km")
    print(f"  Density scale: {ns_units.density_to_cgs:.2e} g/cm³")
    print(f"  Time scale: {ns_units.time_to_cgs * 1000:.4f} ms")
    print()
    
    # Test 3: From config
    config = {
        'SPHType': 'srgsph',
        'units': {
            'type': 'relativistic',
            'preset': 'neutron_star',
            'length_km': 10.0,
            'density_g_cm3': 1e14
        }
    }
    from_cfg = UnitSystem.from_config(config)
    print("From Config (NS preset):")
    print(f"  Type: {from_cfg.get_type_name()}")
    print(f"  Is Relativistic: {from_cfg.is_relativistic()}")
    print()
    
    # Test 4: JSON round-trip
    json_data = sr_test.to_json()
    print("JSON output (keys):", list(json_data.keys()))
    restored = UnitSystem.from_json(json_data)
    print(f"Restored type: {restored.get_type_name()}")
    print(f"Restored c_code: {restored.c_code}")
    print()
    
    print("=== All Tests Passed ===")
