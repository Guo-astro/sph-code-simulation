#!/usr/bin/env python3
"""
Generate all IMBH-cloud simulation category configurations.

This script generates JSON config files for the systematic category system:
- CAT1: Adiabatic (no cooling)
- CAT2: Radiative (Inoue-Inutsuka cooling)  
- CAT3: Full physics (Radiative + enhanced self-gravity)
- CAT_OKA: Exact Oka et al. (2017) orbit parameters from Nature Astronomy paper

Resolution levels:
- A: 61k particles (quick tests)
- B: 200k particles (production)
- C: 1M particles (high resolution)

Impact parameters (standard scenarios):
- b1p5: b=1.5 pc (deep disruption)
- b3: b=3.0 pc (partial disruption)
- b6: b=6.0 pc (weak interaction)

Oka et al. (2017) scenario:
- Uses EXACT orbit from their N-body model: (X,Y)=(9.8,-0.65)pc, (vX,vY)=(-8.19,0.4)km/s
- Pericentre ~1 pc (from tidal mass calculation r_peri=1 pc)
- Original cloud σ_r=0.2 pc; our SPH cloud R=1.13 pc is ~5.6× larger
- This is a "scaled deep encounter" - expect strong tidal disruption

Usage:
    python generate_category_configs.py [--output-dir DIR] [--dry-run]
    python generate_category_configs.py --category CAT_OKA  # Generate only Oka et al. configs
"""

import json
import argparse
from pathlib import Path


# ============================================================================
# CATEGORY DEFINITIONS
# ============================================================================
# ALL categories include cloud self-gravity (useGravity: true).
# The difference is in thermal physics and gravity accuracy (theta parameter).

CATEGORIES = {
    "CAT1": {
        "name": "Adiabatic",
        "description": "Pure gas dynamics with γ=5/3, self-gravity enabled",
        "thermal": {
            "cooling": False,
            "heating": False,
            "comment": "CAT1: Adiabatic evolution with gamma=5/3, self-gravity ON"
        },
        "theta": 0.5,  # Standard Barnes-Hut opening angle
        "self_gravity_note": "enabled"
    },
    "CAT2": {
        "name": "Radiative (Inoue-Inutsuka)",
        "description": "ISM cooling/heating function, self-gravity enabled",
        "thermal": {
            "cooling": True,
            "heating": True,
            "coolingType": "inoue_inutsuka",
            "comment": "CAT2: Inoue & Inutsuka (2008) cooling, self-gravity ON"
        },
        "theta": 0.5,
        "self_gravity_note": "enabled"
    },
    "CAT3": {
        "name": "Full Physics",
        "description": "Radiative + enhanced self-gravity accuracy for fragmentation",
        "thermal": {
            "cooling": True,
            "heating": True,
            "coolingType": "inoue_inutsuka",
            "comment": "CAT3: Full physics - cooling + enhanced self-gravity accuracy"
        },
        "theta": 0.3,  # Lower theta for better fragmentation capture
        "self_gravity_note": "enhanced_accuracy"
    },
    # Special category: Exact Oka et al. (2017) Nature Astronomy parameters
    "CAT_OKA": {
        "name": "Oka et al. (2017) Exact",
        "description": "Exact parameters from Nature Astronomy paper for CO-0.40-0.22",
        "thermal": {
            "cooling": True,
            "heating": True,
            "coolingType": "inoue_inutsuka",
            "comment": "CAT_OKA: Oka et al. exact parameters with radiative cooling"
        },
        "theta": 0.3,  # High accuracy for comparison with observations
        "self_gravity_note": "enhanced_accuracy",
        "is_special": True  # Flag for special handling
    }
}

RESOLUTIONS = {
    "A": {
        "name": "61k (Low)",
        "particles": 61000,
        "folder": "A_61k",
        "ic_path": "sample/imbh_cloud/results/IC/61k/snapshot_0000.csv"
    },
    "B": {
        "name": "200k (Mid)",
        "particles": 200000,
        "folder": "B_200k",
        "ic_path": "sample/imbh_cloud/results/IC/200k/snapshot_0000.csv"
    },
    "C": {
        "name": "1M (High)",
        "particles": 1000000,
        "folder": "C_1M",
        "ic_path": "sample/imbh_cloud/results/IC/1M/snapshot_0000.csv"
    }
}

SCENARIOS = {
    "b1p5": {
        "impact_parameter_pc": 1.5,
        "b_over_r_tidal": 0.41,
        "regime": "deep_disruption",
        "y_offset": 1.5
    },
    "b3": {
        "impact_parameter_pc": 3.0,
        "b_over_r_tidal": 0.83,
        "regime": "partial_disruption",
        "y_offset": 3.0
    },
    "b6": {
        "impact_parameter_pc": 6.0,
        "b_over_r_tidal": 1.65,
        "regime": "weak_interaction",
        "y_offset": 6.0
    },
    # Special scenario: Oka et al. (2017) orbit scaled for R=1.13 pc cloud
    # Original Oka: test particles with σ=0.2 pc, pericentre ~0 pc (head-on)
    # Scaled: pericentre = 1.5 × R_cloud = 1.69 pc (grazing encounter)
    # Orbit designed to produce similar strong tidal effects
    "oka": {
        "impact_parameter_pc": 5.17,  # Calculated from hyperbolic orbit
        "b_over_r_tidal": 0.32,       # r_peri/r_tidal = 1.69/5.24
        "regime": "oka_grazing_encounter",
        "y_offset": 5.17,             # Impact parameter as y-offset
        "x_start": 20.0,              # Start at +x (Oka-style geometry)
        "vx_start": -9.35,            # Approaching BH
        "vy_start": 4.48,             # Small +y for orbit curvature
        "pericentre_pc": 1.69,        # 1.5 × cloud radius
        "eccentricity": 1.241,        # Hyperbolic orbit
        "is_oka_orbit": True          # Flag for special handling
    }
}

# Fixed parameters (standard categories)
CLOUD_MASS_MSUN = 1000
CLOUD_RADIUS_PC = 1.13
BH_MASS_MSUN = 100000
APPROACH_VELOCITY_KMS = 10.0
TIDAL_RADIUS_PC = 5.24  # R_cloud × (M_BH/M_cloud)^(1/3) = 1.13 × 4.64


def generate_config(cat_key: str, res_key: str, scen_key: str) -> dict:
    """Generate a single configuration dictionary."""
    
    cat = CATEGORIES[cat_key]
    res = RESOLUTIONS[res_key]
    scen = SCENARIOS[scen_key]
    
    config_name = f"{cat_key}_{res_key}_{scen_key}"
    
    # Check if this is an Oka-style orbit (different initial conditions)
    is_oka_orbit = scen.get("is_oka_orbit", False)
    
    if is_oka_orbit:
        # Oka-style: start from +x with -y offset
        x_start = scen["x_start"]
        y_start = -scen["y_offset"]  # Negative y-offset
        vx_start = scen["vx_start"]
        vy_start = scen["vy_start"]
        approach_vel = abs(vx_start)
    else:
        # Standard: start from -x approaching along +x
        x_start = -20.0
        y_start = scen["y_offset"]
        vx_start = APPROACH_VELOCITY_KMS
        vy_start = 0.0
        approach_vel = APPROACH_VELOCITY_KMS
    
    return {
        "name": config_name,
        "description": f"Category {cat_key[-1]} ({cat['name']}), Resolution {res_key} ({res['name']}), Impact b={scen['impact_parameter_pc']}pc",
        
        "category_info": {
            "category": cat_key,
            "category_name": cat["name"],
            "resolution": res_key,
            "resolution_name": res["name"],
            "scenario": scen_key,
            "impact_parameter_pc": scen["impact_parameter_pc"]
        },
        
        "physics_summary": {
            "cloud_mass_Msun": CLOUD_MASS_MSUN,
            "cloud_radius_pc": CLOUD_RADIUS_PC,
            "bh_mass_Msun": BH_MASS_MSUN,
            "impact_parameter_pc": scen["impact_parameter_pc"],
            "approach_velocity_kms": approach_vel,
            "tidal_radius_pc": TIDAL_RADIUS_PC,
            "b_over_r_tidal": scen["b_over_r_tidal"],
            "regime": scen["regime"],
            "thermal_physics": "adiabatic" if not cat["thermal"]["cooling"] else "radiative_inoue_inutsuka",
            "self_gravity": cat["self_gravity_note"],
            "particles": res["particles"],
            **({"pericentre_pc": scen["pericentre_pc"], "eccentricity": scen["eccentricity"]} if is_oka_orbit else {})
        },
        
        "units": {
            "type": "IMBH_ENCOUNTER",
            "length_pc": 1.0,
            "mass_1e3Msun": 1.0,
            "velocity_kms": 1.0,
            "comment": "1 code length = 1 pc, 1 code mass = 1000 Msun, 1 code velocity = 1 km/s"
        },
        
        "outputDirectory": f"sample/imbh_cloud/results/{cat_key}/{res['folder']}/{scen_key}",
        "outputFormat": "csv",
        "enableEnergyFile": True,
        
        "SPHType": "gsph",
        "use2ndOrderGSPH": False,
        "kernel": "wendland",
        "neighborNumber": 50,
        "riemannSolver": "hll",
        "gamma": 5.0/3.0,
        
        "cflSound": 0.3,
        "cflForce": 0.125,
        
        "maxTreeLevel": 20,
        "leafParticleNumber": 1,
        "iterativeSmoothingLength": True,
        "periodic": False,
        
        "startTime": 0.0,
        "endTime": 4.0,
        "outputTime": 0.02,
        
        "checkpointInterval": 1000,
        "checkpointPath": f"sample/imbh_cloud/results/{cat_key}/{res['folder']}/{scen_key}/checkpoints",
        
        "useGravity": True,
        "G": 1.0,
        "theta": cat["theta"],
        "gravitySofteningType": "wendland_c4",
        
        "resumeFromSnapshot": res["ic_path"],
        "resetTimeOnResume": True,
        
        "initialCondition": {
            "transform": {
                "translate": [x_start, y_start, 0.0],
                "velocity_boost": [vx_start, vy_start, 0.0]
            },
            "comment": f"Cloud starts at ({x_start}, {y_start}, 0) pc with velocity ({vx_start}, {vy_start}, 0) km/s" + 
                       (f" [Oka-style hyperbolic orbit, pericentre={scen['pericentre_pc']:.2f} pc]" if is_oka_orbit else "")
        },
        
        "imbh_parameters": {
            "enabled": True,
            "M_BH": float(BH_MASS_MSUN),
            "BH_initial_position": [0.0, 0.0, 0.0],
            "BH_initial_velocity": [0.0, 0.0, 0.0],
            "softening_epsilon": 0.05,
            "cloud_initial_position": [x_start, y_start, 0.0],
            "cloud_initial_velocity": [vx_start, vy_start, 0.0]
        },
        
        "externalForces": {
            "pointMassBH": {
                "enabled": True,
                "mass": 100.0,  # In code units (100 × 1000 M☉ = 10^5 M☉)
                "position": [0.0, 0.0, 0.0],
                "softening": 0.05,
                "comment": "IMBH at origin: M_BH = 100 code = 1e5 Msun"
            }
        },
        
        "thermal": cat["thermal"]
    }


def write_config(config: dict, output_dir: Path, dry_run: bool = False):
    """Write configuration to JSON file."""
    
    cat_key = config["category_info"]["category"]
    res_folder = RESOLUTIONS[config["category_info"]["resolution"]]["folder"]
    scen_key = config["category_info"]["scenario"]
    
    # Create directory structure: categories/CAT1/A_60k/
    config_dir = output_dir / cat_key / res_folder
    config_file = config_dir / f"{scen_key}.json"
    
    if dry_run:
        print(f"  [DRY-RUN] Would create: {config_file}")
        return
    
    config_dir.mkdir(parents=True, exist_ok=True)
    
    with open(config_file, 'w') as f:
        json.dump(config, f, indent=2)
    
    print(f"  ✓ Created: {config_file}")


def main():
    parser = argparse.ArgumentParser(description="Generate IMBH-cloud category configurations")
    parser.add_argument("--output-dir", type=str, 
                        default="sample/imbh_cloud/config/presets/categories",
                        help="Output directory for config files")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show what would be created without writing files")
    parser.add_argument("--category", type=str, choices=["CAT1", "CAT2", "CAT3", "CAT_OKA"],
                        help="Generate only specific category")
    parser.add_argument("--resolution", type=str, choices=["A", "B", "C"],
                        help="Generate only specific resolution")
    parser.add_argument("--scenario", type=str, choices=["b1p5", "b3", "b6", "oka"],
                        help="Generate only specific scenario")
    parser.add_argument("--oka-only", action="store_true",
                        help="Generate only CAT_OKA configurations with 'oka' scenario")
    
    args = parser.parse_args()
    output_dir = Path(args.output_dir)
    
    print("=" * 70)
    print("  IMBH-Cloud Category Configuration Generator")
    print("=" * 70)
    print()
    
    # Handle --oka-only shortcut
    if args.oka_only:
        categories = ["CAT_OKA"]
        resolutions = list(RESOLUTIONS.keys())
        scenarios = ["oka"]
    else:
        # Determine which configs to generate
        # CAT_OKA only makes sense with "oka" scenario
        if args.category == "CAT_OKA":
            categories = ["CAT_OKA"]
            scenarios = ["oka"]
        elif args.scenario == "oka":
            # "oka" scenario only available for CAT_OKA
            categories = ["CAT_OKA"]
            scenarios = ["oka"]
        else:
            categories = [args.category] if args.category else [k for k in CATEGORIES.keys() if k != "CAT_OKA"]
            scenarios = [args.scenario] if args.scenario else [k for k in SCENARIOS.keys() if k != "oka"]
        
        resolutions = [args.resolution] if args.resolution else list(RESOLUTIONS.keys())
    
    total = len(categories) * len(resolutions) * len(scenarios)
    print(f"Generating {total} configuration files...")
    print(f"Output directory: {output_dir}")
    print()
    
    count = 0
    for cat_key in categories:
        print(f"\n[{cat_key}] {CATEGORIES[cat_key]['name']}")
        for res_key in resolutions:
            for scen_key in scenarios:
                config = generate_config(cat_key, res_key, scen_key)
                write_config(config, output_dir, dry_run=args.dry_run)
                count += 1
    
    print()
    print("=" * 70)
    print(f"  Generated {count} configuration files")
    print("=" * 70)
    
    if not args.dry_run:
        print()
        print("Next steps:")
        print("  1. Generate initial conditions for each resolution:")
        print("     make -f sample/imbh_cloud/Makefile.hydrostatic hydrostatic_61k")
        print("     make -f sample/imbh_cloud/Makefile.hydrostatic hydrostatic_200k")
        print("     make -f sample/imbh_cloud/Makefile.hydrostatic hydrostatic_1M")
        print()
        print("  2. Run simulations:")
        print("     make -f sample/imbh_cloud/Makefile.categories CAT1_A_b3")
        print("     make -f sample/imbh_cloud/Makefile.categories CAT_OKA_A_oka")
        print("     make -f sample/imbh_cloud/Makefile.categories CAT3_B_b3")


if __name__ == "__main__":
    main()
