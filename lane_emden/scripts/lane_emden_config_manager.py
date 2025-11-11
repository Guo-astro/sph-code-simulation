#!/usr/bin/env python3
"""
Lane-Emden Configuration Manager

Unified tool for managing Lane-Emden preset configurations.
Supports creating, updating, listing, and validating presets.

Usage:
    python3 config_manager.py list
    python3 config_manager.py create --name custom_n3 --polytrope 3.0 --dimension 3
    python3 config_manager.py update --preset polytrope_n1.5_3d --relax-steps 200
    python3 config_manager.py validate --preset polytrope_n1.5_3d
"""

import json
import argparse
import sys
from pathlib import Path
from typing import Dict, Any, Optional
import shutil

# Directory paths relative to script location
SCRIPT_DIR = Path(__file__).parent
LANE_EMDEN_ROOT = SCRIPT_DIR / "lane_emden"
PRESETS_DIR = LANE_EMDEN_ROOT / "config" / "presets"
TEMPLATES_DIR = LANE_EMDEN_ROOT / "config" / "templates"
DATA_DIR = LANE_EMDEN_ROOT / "data" / "numerical_solutions"
RESULTS_DIR = LANE_EMDEN_ROOT / "results"

# Physics constants
POLYTROPE_TABLE = {
    0.0: {"gamma": float('inf'), "description": "Incompressible fluid"},
    0.5: {"gamma": 3.0, "description": "Very stiff EOS (white dwarfs)"},
    1.0: {"gamma": 2.0, "description": "Isothermal limit (molecular clouds)"},
    1.5: {"gamma": 5.0/3.0, "description": "Monatomic ideal gas (main sequence stars)"},
    2.0: {"gamma": 1.5, "description": "Intermediate polytrope"},
    3.0: {"gamma": 4.0/3.0, "description": "Radiation-dominated (massive stars)"},
    4.0: {"gamma": 1.25, "description": "Degenerate intermediate"},
    5.0: {"gamma": 1.2, "description": "Degenerate gas (low-mass stars)"}
}


class ConfigManager:
    """Manage Lane-Emden preset configurations"""
    
    def __init__(self):
        self.ensure_directories()
    
    def ensure_directories(self):
        """Ensure required directories exist"""
        PRESETS_DIR.mkdir(parents=True, exist_ok=True)
        TEMPLATES_DIR.mkdir(parents=True, exist_ok=True)
        RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    
    def list_presets(self, verbose: bool = False):
        """List all available presets"""
        if not PRESETS_DIR.exists():
            print("No presets directory found.")
            return
        
        presets = sorted(PRESETS_DIR.glob("*.json"))
        
        if not presets:
            print("No presets found.")
            return
        
        print(f"\n{'='*70}")
        print(f"Available Lane-Emden Presets ({len(presets)} total)")
        print(f"{'='*70}\n")
        
        for preset_file in presets:
            with open(preset_file, 'r') as f:
                config = json.load(f)
            
            name = config.get("name", preset_file.stem)
            desc = config.get("description", "No description")
            dim = config.get("dimension", "?")
            n = config["physics"].get("polytropic_index", "?")
            gamma = config["physics"].get("adiabatic_index", "?")
            
            print(f"📋 {name}")
            print(f"   Dimension: {dim}D | n={n} | γ={gamma:.3f}")
            print(f"   {desc}")
            
            if verbose:
                relax_steps = config["relaxation"].get("steps", "?")
                particles = config["initial_conditions"].get("particle_count", "?")
                print(f"   Relaxation: {relax_steps} steps | Particles: {particles}")
                use_case = config.get("use_case", "")
                if use_case:
                    print(f"   Use case: {use_case}")
            
            print()
    
    def create_preset(self, name: str, dimension: int, polytrope: float,
                     gamma: Optional[float] = None, particles: Optional[int] = None,
                     relax_steps: int = 100):
        """Create a new preset from template"""
        
        # Auto-calculate gamma if not provided
        if gamma is None:
            if polytrope in POLYTROPE_TABLE:
                gamma = POLYTROPE_TABLE[polytrope]["gamma"]
            elif polytrope > 0:
                gamma = 1.0 + 1.0 / polytrope
            else:
                print(f"Error: Cannot auto-calculate γ for n={polytrope}")
                return False
        
        # Default particle counts
        if particles is None:
            particles = 5400 if dimension == 3 else 3600
        
        # Load template
        template_file = TEMPLATES_DIR / f"base_{dimension}d.json"
        if not template_file.exists():
            print(f"Error: Template {template_file} not found.")
            return False
        
        with open(template_file, 'r') as f:
            config = json.load(f)
        
        # Update configuration
        config["name"] = name
        config["dimension"] = dimension
        config["physics"]["polytropic_index"] = polytrope
        config["physics"]["adiabatic_index"] = gamma
        config["initial_conditions"]["particle_count"] = particles
        config["relaxation"]["steps"] = relax_steps
        
        # Set data file path
        data_file = f"lane_emden/data/numerical_solutions/{dimension}d/n{polytrope}.dat"
        config["initial_conditions"]["data_file"] = data_file
        
        # Set output directory
        config["simulation"]["output_directory"] = f"lane_emden/results/{name}"
        
        # Add description
        if polytrope in POLYTROPE_TABLE:
            physics_desc = POLYTROPE_TABLE[polytrope]["description"]
            config["description"] = f"{dimension}D n={polytrope} polytrope (γ={gamma:.3f}) - {physics_desc}"
        else:
            config["description"] = f"{dimension}D n={polytrope} polytrope (γ={gamma:.3f})"
        
        # Save preset
        preset_file = PRESETS_DIR / f"{name}.json"
        with open(preset_file, 'w') as f:
            json.dump(config, f, indent=2)
        
        print(f"✓ Created preset: {preset_file}")
        print(f"  n={polytrope}, γ={gamma:.6f}, {dimension}D, {particles} particles")
        
        # Check if data file exists
        data_path = DATA_DIR / f"{dimension}d" / f"n{polytrope}.dat"
        if not data_path.exists():
            print(f"\n⚠️  Warning: Data file does not exist: {data_path}")
            print(f"   Generate it with: python3 scripts/generators/solve_lane_emden.py --n {polytrope} --dim {dimension}")
        
        return True
    
    def update_preset(self, preset_name: str, **kwargs):
        """Update an existing preset"""
        preset_file = PRESETS_DIR / f"{preset_name}.json"
        
        if not preset_file.exists():
            print(f"Error: Preset '{preset_name}' not found.")
            return False
        
        with open(preset_file, 'r') as f:
            config = json.load(f)
        
        # Update fields
        updates = []
        if kwargs.get("relax_steps") is not None:
            config["relaxation"]["steps"] = kwargs["relax_steps"]
            updates.append(f"relaxation.steps = {kwargs['relax_steps']}")
        
        if kwargs.get("end_time") is not None:
            config["simulation"]["end_time"] = kwargs["end_time"]
            updates.append(f"simulation.end_time = {kwargs['end_time']}")
        
        if kwargs.get("particles") is not None:
            config["initial_conditions"]["particle_count"] = kwargs["particles"]
            updates.append(f"particle_count = {kwargs['particles']}")
        
        if kwargs.get("gamma") is not None:
            config["physics"]["adiabatic_index"] = kwargs["gamma"]
            updates.append(f"adiabatic_index = {kwargs['gamma']}")
        
        if not updates:
            print("No updates specified.")
            return False
        
        # Save updated config
        with open(preset_file, 'w') as f:
            json.dump(config, f, indent=2)
        
        print(f"✓ Updated preset: {preset_name}")
        for update in updates:
            print(f"  - {update}")
        
        return True
    
    def validate_preset(self, preset_name: str):
        """Validate a preset configuration"""
        preset_file = PRESETS_DIR / f"{preset_name}.json"
        
        if not preset_file.exists():
            print(f"❌ Error: Preset '{preset_name}' not found.")
            return False
        
        try:
            with open(preset_file, 'r') as f:
                config = json.load(f)
        except json.JSONDecodeError as e:
            print(f"❌ Error: Invalid JSON in {preset_file}")
            print(f"   {e}")
            return False
        
        errors = []
        warnings = []
        
        # Required fields
        required_sections = ["physics", "initial_conditions", "relaxation", "simulation", "numerical"]
        for section in required_sections:
            if section not in config:
                errors.append(f"Missing required section: {section}")
        
        # Check data file exists
        data_file = config.get("initial_conditions", {}).get("data_file")
        if data_file:
            # Convert to absolute path
            data_path = Path(data_file)
            if not data_path.is_absolute():
                data_path = LANE_EMDEN_ROOT.parent / data_file
            
            if not data_path.exists():
                warnings.append(f"Data file not found: {data_file}")
        
        # Check physics consistency
        n = config.get("physics", {}).get("polytropic_index")
        gamma = config.get("physics", {}).get("adiabatic_index")
        if n is not None and gamma is not None and n > 0:
            expected_gamma = 1.0 + 1.0 / n
            if abs(gamma - expected_gamma) > 0.01:
                warnings.append(f"γ={gamma:.3f} inconsistent with n={n} (expected γ={expected_gamma:.3f})")
        
        # Check dimension vs data file
        dim = config.get("dimension")
        if dim and data_file:
            if f"/{dim}d/" not in data_file and f"_{dim}d" not in data_file:
                warnings.append(f"Dimension {dim}D may not match data file: {data_file}")
        
        # Print results
        print(f"\nValidation results for: {preset_name}")
        print(f"{'='*70}")
        
        if errors:
            print(f"\n❌ Errors ({len(errors)}):")
            for error in errors:
                print(f"   - {error}")
        
        if warnings:
            print(f"\n⚠️  Warnings ({len(warnings)}):")
            for warning in warnings:
                print(f"   - {warning}")
        
        if not errors and not warnings:
            print("\n✓ Configuration is valid!")
            
            # Print summary
            print(f"\nConfiguration summary:")
            print(f"  Dimension: {dim}D")
            print(f"  Polytrope: n={n}, γ={gamma:.6f}")
            print(f"  Particles: {config['initial_conditions'].get('particle_count', '?')}")
            print(f"  Relaxation: {config['relaxation'].get('steps', '?')} steps")
        
        return len(errors) == 0
    
    def generate_sample_config(self, preset_name: str, sample_dir: Path):
        """Generate old-style sample config from preset (for backward compatibility)"""
        preset_file = PRESETS_DIR / f"{preset_name}.json"
        
        if not preset_file.exists():
            print(f"Error: Preset '{preset_name}' not found.")
            return False
        
        with open(preset_file, 'r') as f:
            preset = json.load(f)
        
        # Convert to old format
        sample_config = {
            "outputDirectory": preset["simulation"]["output_directory"],
            "endTime": preset["simulation"]["end_time"],
            "outputTime": preset["simulation"]["output_time"],
            "avAlpha": preset["artificial_viscosity"]["alpha"],
            "neighborNumber": preset["numerical"]["neighbor_number"],
            "useBalsaraSwitch": preset["artificial_viscosity"]["use_balsara_switch"],
            "useTimeDependentAV": preset["artificial_viscosity"]["use_time_dependent_av"],
            "useArtificialConductivity": preset["artificial_viscosity"].get("use_artificial_conductivity", False),
            "leafParticleNumber": preset["numerical"]["leaf_particle_number"],
            "gamma": preset["physics"]["adiabatic_index"],
            "kernel": preset["numerical"]["kernel"],
            "N": preset["initial_conditions"].get("shells", 30),
            "periodic": preset["numerical"]["periodic"],
            "useGravity": preset["numerical"]["use_gravity"],
            "G": preset["physics"]["gravitational_constant"],
            "SPHType": preset["numerical"]["sph_type"],
            "iterativeSmoothingLength": preset["numerical"]["iterative_smoothing_length"],
            "useRelaxation": preset["relaxation"]["enabled"],
            "relaxationSteps": preset["relaxation"]["steps"],
            "relaxationOutputFreq": preset["relaxation"]["output_frequency"],
            "relaxationOnly": False,
            "output": {
                "formats": [
                    {
                        "type": "csv",
                        "precision": 16
                    },
                    {
                        "type": "hdf5",
                        "compression": 6
                    }
                ],
                "enableEnergyFile": True
            }
        }
        
        output_file = sample_dir / "lane_emden.json"
        with open(output_file, 'w') as f:
            json.dump(sample_config, f, indent=2)
        
        print(f"✓ Generated sample config: {output_file}")
        return True


def main():
    parser = argparse.ArgumentParser(
        description="Lane-Emden Configuration Manager",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s list
  %(prog)s list --verbose
  %(prog)s create --name custom_n3_3d --polytrope 3.0 --dimension 3
  %(prog)s update --preset polytrope_n1.5_3d --relax-steps 200
  %(prog)s validate --preset polytrope_n1.5_3d
        """
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Command to execute')
    
    # List command
    list_parser = subparsers.add_parser('list', help='List all presets')
    list_parser.add_argument('--verbose', '-v', action='store_true', help='Show detailed information')
    
    # Create command
    create_parser = subparsers.add_parser('create', help='Create new preset')
    create_parser.add_argument('--name', required=True, help='Preset name (e.g., custom_n3_3d)')
    create_parser.add_argument('--polytrope', '-n', type=float, required=True, help='Polytropic index')
    create_parser.add_argument('--dimension', '-d', type=int, choices=[2, 3], required=True, help='Dimension (2 or 3)')
    create_parser.add_argument('--gamma', '-g', type=float, help='Adiabatic index (auto-calculated if not provided)')
    create_parser.add_argument('--particles', '-p', type=int, help='Particle count')
    create_parser.add_argument('--relax-steps', '-r', type=int, default=100, help='Relaxation steps (default: 100)')
    
    # Update command
    update_parser = subparsers.add_parser('update', help='Update existing preset')
    update_parser.add_argument('--preset', required=True, help='Preset name to update')
    update_parser.add_argument('--relax-steps', type=int, help='Update relaxation steps')
    update_parser.add_argument('--end-time', type=float, help='Update simulation end time')
    update_parser.add_argument('--particles', type=int, help='Update particle count')
    update_parser.add_argument('--gamma', type=float, help='Update adiabatic index')
    
    # Validate command
    validate_parser = subparsers.add_parser('validate', help='Validate preset configuration')
    validate_parser.add_argument('--preset', required=True, help='Preset name to validate')
    
    # Apply command (generate C++ config from preset)
    apply_parser = subparsers.add_parser('apply', help='Apply preset to sample directory (generate C++ config)')
    apply_parser.add_argument('--preset', required=True, help='Preset name to apply')
    apply_parser.add_argument('--sample-dir', default='sample/lane_emden', help='Sample directory path (default: sample/lane_emden)')
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        return 1
    
    manager = ConfigManager()
    
    if args.command == 'list':
        manager.list_presets(verbose=args.verbose)
    
    elif args.command == 'create':
        success = manager.create_preset(
            name=args.name,
            dimension=args.dimension,
            polytrope=args.polytrope,
            gamma=args.gamma,
            particles=args.particles,
            relax_steps=args.relax_steps
        )
        return 0 if success else 1
    
    elif args.command == 'update':
        success = manager.update_preset(
            preset_name=args.preset,
            relax_steps=args.relax_steps,
            end_time=args.end_time,
            particles=args.particles,
            gamma=args.gamma
        )
        return 0 if success else 1
    
    elif args.command == 'validate':
        success = manager.validate_preset(args.preset)
        return 0 if success else 1
    
    elif args.command == 'apply':
        sample_dir = Path(args.sample_dir)
        sample_dir.mkdir(parents=True, exist_ok=True)
        success = manager.generate_sample_config(args.preset, sample_dir)
        return 0 if success else 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
