#!/usr/bin/env python3
"""
Generate Initial Conditions for BNS Cocoon Simulations

Creates SRGSPH initial conditions from:
1. Radice+2018 GR simulation data
2. Analytic fits (Hotokezaka, Bauswein, etc.)

Output: JSON config file for SPH simulation

Usage:
    # From Radice+2018 data
    python generate_ic.py --source radice2018 --model SFHo_M140140 \
        --t0 1.0 --output config/presets/ejecta_radice.json
    
    # From analytic profile
    python generate_ic.py --source analytic --profile hotokezaka \
        --total-mass 0.01 --avg-velocity 0.2 \
        --output config/presets/ejecta_hotokezaka.json
    
    # With cocoon injection
    python generate_ic.py --source analytic --profile hotokezaka \
        --cocoon-energy 1e50 --cocoon-angle 30 \
        --output config/presets/cocoon_fiducial.json
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from radice_loader import RadiceDataLoader, EjectaProfile as RadiceEjectaProfile
from homologous_expansion import HomologousExpander, SpatialProfile
from analytic_profiles import (
    EjectaParameters, 
    create_profile,
    HotokezakaProfile,
    BausweinProfile,
)
from particle_sampler import ParticleSampler, ParticleDistribution


def generate_from_radice(args) -> dict:
    """
    Generate IC from Radice+2018 data.
    
    Args:
        args: Parsed command-line arguments
        
    Returns:
        Config dictionary for SPH simulation
    """
    # Load data
    loader = RadiceDataLoader(args.radice_dir)
    profile = loader.load_model(args.model)
    loader.print_summary(profile)
    
    # Homologous expansion
    expander = HomologousExpander(
        t0=args.t0,
        length_unit=args.length_unit,
    )
    
    if profile.dM_dv_dOmega is not None and args.dimension == 2:
        # 2D profile
        spatial = expander.expand_2d(
            profile.v_inf, 
            profile.cos_theta,
            profile.dM_dv_dOmega
        )
    else:
        # 1D profile
        spatial = expander.expand_1d(profile.v_inf, profile.dM_dv)
    
    # Add ISM if requested
    if args.n_ism > 0:
        spatial = expander.add_external_medium(spatial, n_ism=args.n_ism)
    
    # Sample particles
    sampler = ParticleSampler(
        n_particles=args.n_particles,
        gamma=args.gamma,
        random_seed=args.seed,
    )
    
    if args.dimension == 2 and spatial.rho_2d is not None:
        particles = sampler.sample_2d(
            spatial.r, spatial.theta,
            spatial.rho_2d, spatial.v_r_2d,
            polytropic_K=args.polytropic_k,
        )
    else:
        particles = sampler.sample_1d(
            spatial.r, spatial.rho, spatial.v_r,
            polytropic_K=args.polytropic_k,
        )
    
    # Add cocoon if requested
    if args.cocoon_energy > 0:
        particles = sampler.add_cocoon(
            particles,
            cocoon_energy=args.cocoon_energy,
            cocoon_opening_angle=args.cocoon_angle,
            cocoon_velocity=args.cocoon_velocity,
            n_cocoon_particles=args.n_cocoon_particles,
        )
    
    # Build config
    config = build_sph_config(args, particles, spatial, source='radice2018', model=args.model)
    
    return config


def generate_from_analytic(args) -> dict:
    """
    Generate IC from analytic profile.
    
    Args:
        args: Parsed command-line arguments
        
    Returns:
        Config dictionary for SPH simulation
    """
    # Create profile parameters
    params = EjectaParameters(
        total_mass=args.total_mass,
        avg_velocity=args.avg_velocity,
        max_velocity=args.max_velocity,
        min_velocity=args.min_velocity,
        velocity_power=args.velocity_power,
        angular_conc=args.angular_conc,
        polar_opening=args.polar_opening,
    )
    
    # Create profile
    profile = create_profile(args.profile, params=params)
    
    # Generate velocity and angle grids
    v = np.linspace(params.min_velocity, params.max_velocity, 100)
    theta = np.linspace(0.01, np.pi - 0.01, 50)
    
    # Get distributions
    dM_dv = profile.velocity_distribution(v)
    
    # Homologous expansion
    expander = HomologousExpander(
        t0=args.t0,
        length_unit=args.length_unit,
    )
    
    if args.dimension == 2:
        dM_dv_dOmega = profile.velocity_angle_distribution(v, theta)
        spatial = expander.expand_2d(v, np.cos(theta), dM_dv_dOmega)
    else:
        spatial = expander.expand_1d(v, dM_dv)
    
    # Add ISM
    if args.n_ism > 0:
        spatial = expander.add_external_medium(spatial, n_ism=args.n_ism)
    
    # Sample particles
    sampler = ParticleSampler(
        n_particles=args.n_particles,
        gamma=args.gamma,
        random_seed=args.seed,
    )
    
    if args.dimension == 2 and spatial.rho_2d is not None:
        particles = sampler.sample_2d(
            spatial.r, spatial.theta,
            spatial.rho_2d, spatial.v_r_2d,
            polytropic_K=args.polytropic_k,
        )
    else:
        particles = sampler.sample_1d(
            spatial.r, spatial.rho, spatial.v_r,
            polytropic_K=args.polytropic_k,
        )
    
    # Add cocoon
    if args.cocoon_energy > 0:
        particles = sampler.add_cocoon(
            particles,
            cocoon_energy=args.cocoon_energy,
            cocoon_opening_angle=args.cocoon_angle,
            cocoon_velocity=args.cocoon_velocity,
            n_cocoon_particles=args.n_cocoon_particles,
        )
    
    # Build config
    config = build_sph_config(args, particles, spatial, 
                              source='analytic', model=args.profile)
    
    return config


def build_sph_config(args, particles: ParticleDistribution, 
                     spatial: SpatialProfile,
                     source: str, model: str) -> dict:
    """
    Build SPH simulation config dictionary.
    
    Args:
        args: Command-line arguments
        particles: Sampled particle distribution
        spatial: Spatial profile
        source: Data source type
        model: Model name
        
    Returns:
        Config dictionary
    """
    # Compute domain size
    r_max = np.max(np.linalg.norm(particles.positions, axis=1)) * 1.5
    
    config = {
        "simulation": {
            "name": f"BNS Cocoon - {model}",
            "description": f"Cocoon shock breakout simulation from {source} ejecta model"
        },
        
        "output": {
            "directory": f"sample/bns_cocoon/results/{args.output.stem}",
            "formats": ["csv"],
            "energyFile": True
        },
        
        "time": {
            "start": 0.0,
            "end": args.end_time,
            "outputInterval": args.output_interval,
            "energyInterval": args.energy_interval
        },
        
        "domain": {
            "dimension": args.dimension,
            "size": r_max,
            "boundary": "open"
        },
        
        "physics": {
            "gamma": args.gamma,
            "G": 0.0,  # No gravity for SR simulations
            "c": 1.0,  # Relativistic
            "neighborNumber": args.neighbor_number
        },
        
        "sph": {
            "type": "srgsph",
            "kernel": "gaussian",
            "balsara": True,
            "alpha": 1.0,
            "order": 2
        },
        
        "timestep": {
            "CFL_sound": 0.3,
            "CFL_force": 0.25
        },
        
        "sample": {
            "type": "bns_cocoon",
            "source": source,
            "model": model,
        },
        
        "bns_cocoon": {
            "t0": args.t0,
            "length_unit_cm": spatial.length_unit,
            "density_unit_cgs": spatial.density_unit,
            "total_mass_msun": args.total_mass if hasattr(args, 'total_mass') else 0.01,
            "avg_velocity_c": args.avg_velocity if hasattr(args, 'avg_velocity') else 0.2,
            "n_particles": particles.n_particles,
            "polytropic_K": args.polytropic_k,
            "n_ism_per_cc": args.n_ism,
        },
        
        "cocoon": {
            "enabled": args.cocoon_energy > 0,
            "energy_erg": args.cocoon_energy,
            "opening_angle_deg": args.cocoon_angle,
            "velocity_c": args.cocoon_velocity,
            "n_particles": args.n_cocoon_particles,
        },
        
        "unit_system": {
            "type": "relativistic",
            "length_cm": spatial.length_unit,
            "mass_g": spatial.density_unit * spatial.length_unit**3,
            "c": 1.0,
        },
    }
    
    # Add initial particle data
    config["initial_particles"] = particles.to_dict()
    
    return config


def main():
    parser = argparse.ArgumentParser(
        description='Generate BNS cocoon simulation initial conditions',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # From Radice+2018 data
  python generate_ic.py --source radice2018 --model SFHo_M140140 \\
      --radice-dir ../data/radice2018 --output config.json
  
  # From analytic Hotokezaka profile
  python generate_ic.py --source analytic --profile hotokezaka \\
      --total-mass 0.01 --avg-velocity 0.2 --output config.json
  
  # With cocoon energy injection
  python generate_ic.py --source analytic --profile hotokezaka \\
      --cocoon-energy 1e50 --cocoon-angle 30 --output config.json
"""
    )
    
    # Source selection
    parser.add_argument('--source', choices=['radice2018', 'analytic'],
                       required=True, help='Data source for ejecta model')
    
    # Radice-specific options
    parser.add_argument('--radice-dir', default='../data/radice2018',
                       help='Path to Radice+2018 data directory')
    parser.add_argument('--model', default='SFHo_M140140',
                       help='Model name from Radice+2018 dataset')
    
    # Analytic profile options
    parser.add_argument('--profile', 
                       choices=['hotokezaka', 'bauswein', 'powerlaw', 'fasttail'],
                       default='hotokezaka',
                       help='Analytic profile type')
    parser.add_argument('--total-mass', type=float, default=0.01,
                       help='Total ejecta mass [M_sun]')
    parser.add_argument('--avg-velocity', type=float, default=0.2,
                       help='Average ejecta velocity [c]')
    parser.add_argument('--max-velocity', type=float, default=0.5,
                       help='Maximum velocity [c]')
    parser.add_argument('--min-velocity', type=float, default=0.05,
                       help='Minimum velocity [c]')
    parser.add_argument('--velocity-power', type=float, default=-2.0,
                       help='Velocity power-law index')
    parser.add_argument('--angular-conc', type=float, default=2.0,
                       help='Angular concentration parameter')
    parser.add_argument('--polar-opening', type=float, default=60.0,
                       help='Polar opening angle [degrees]')
    
    # Expansion parameters
    parser.add_argument('--t0', type=float, default=1.0,
                       help='Reference time after merger [seconds]')
    parser.add_argument('--length-unit', type=float, default=1e10,
                       help='Length unit [cm]')
    
    # External medium
    parser.add_argument('--n-ism', type=float, default=1e-3,
                       help='ISM number density [cm^-3]')
    
    # Simulation parameters
    parser.add_argument('--dimension', '-d', type=int, choices=[1, 2], default=1,
                       help='Simulation dimension')
    parser.add_argument('--n-particles', type=int, default=10000,
                       help='Number of ejecta particles')
    parser.add_argument('--gamma', type=float, default=4/3,
                       help='Adiabatic index')
    parser.add_argument('--polytropic-k', type=float, default=1.0,
                       help='Polytropic constant K')
    parser.add_argument('--neighbor-number', type=int, default=30,
                       help='SPH neighbor number')
    parser.add_argument('--seed', type=int, default=42,
                       help='Random seed')
    
    # Cocoon parameters
    parser.add_argument('--cocoon-energy', type=float, default=0.0,
                       help='Cocoon energy [erg]')
    parser.add_argument('--cocoon-angle', type=float, default=30.0,
                       help='Cocoon half-opening angle [degrees]')
    parser.add_argument('--cocoon-velocity', type=float, default=0.9,
                       help='Initial cocoon velocity [c]')
    parser.add_argument('--n-cocoon-particles', type=int, default=1000,
                       help='Number of cocoon particles')
    
    # Time parameters
    parser.add_argument('--end-time', type=float, default=100.0,
                       help='Simulation end time [code units]')
    parser.add_argument('--output-interval', type=float, default=1.0,
                       help='Snapshot output interval')
    parser.add_argument('--energy-interval', type=float, default=0.1,
                       help='Energy log interval')
    
    # Output
    parser.add_argument('--output', '-o', type=Path, required=True,
                       help='Output JSON config file')
    
    args = parser.parse_args()
    
    # Generate IC based on source
    if args.source == 'radice2018':
        config = generate_from_radice(args)
    else:
        config = generate_from_analytic(args)
    
    # Save config
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, 'w') as f:
        json.dump(config, f, indent=2)
    
    print(f"\n✓ Generated config: {args.output}")
    print(f"  Particles: {config['bns_cocoon']['n_particles']}")
    print(f"  Dimension: {config['domain']['dimension']}D")
    print(f"  End time: {config['time']['end']}")
    if config['cocoon']['enabled']:
        print(f"  Cocoon energy: {config['cocoon']['energy_erg']:.2e} erg")


if __name__ == '__main__':
    main()
