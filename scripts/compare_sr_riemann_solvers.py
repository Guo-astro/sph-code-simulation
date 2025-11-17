#!/usr/bin/env python3
"""
Compare EXACT vs ITERATIVE Riemann solvers for SR-GSPH
Runs both solvers on SR Sod shock tube and plots results
"""

import subprocess
import json
import os
import shutil
from pathlib import Path

def run_simulation(config_file, output_name):
    """Run simulation and return output directory"""
    print(f"\n{'='*60}")
    print(f"Running {output_name}...")
    print(f"{'='*60}")
    
    result = subprocess.run(
        ["./build/sph", config_file],
        capture_output=True,
        text=True
    )
    
    if result.returncode != 0:
        print(f"ERROR in {output_name}:")
        print(result.stderr)
        return None
    
    # Extract output directory from config
    with open(config_file) as f:
        config = json.load(f)
    return config["outputDirectory"]

def main():
    # Base directory
    base_dir = Path(__file__).parent
    
    # Create test configs
    configs = {
        "EXACT": {
            "outputDirectory": "sample/sr_sod/results/comparison_exact",
            "SPHType": "srgsph",
            "riemannSolverSRGSPH": "EXACT",
            "endTime": 0.35,
            "outputTime": 0.01,
        },
        "ITERATIVE": {
            "outputDirectory": "sample/sr_sod/results/comparison_iterative",
            "SPHType": "srgsph",
            "riemannSolverSRGSPH": "ITERATIVE",
            "endTime": 0.35,
            "outputTime": 0.01,
        }
    }
    
    # Load base config
    base_config_file = "sr_sod.json"
    with open(base_config_file) as f:
        base_config = json.load(f)
    
    results = {}
    
    for solver_name, overrides in configs.items():
        # Create config
        config = base_config.copy()
        config.update(overrides)
        
        # Write config file
        config_file = f"sr_sod_{solver_name.lower()}_test.json"
        with open(config_file, 'w') as f:
            json.dump(config, f, indent=2)
        
        # Run simulation
        output_dir = run_simulation(config_file, solver_name)
        if output_dir:
            results[solver_name] = output_dir
        
        # Clean up config
        os.remove(config_file)
    
    print(f"\n{'='*60}")
    print("Comparison Complete!")
    print(f"{'='*60}")
    print("\nResults:")
    for solver, output in results.items():
        print(f"  {solver:12s}: {output}")
    
    print("\nTo visualize:")
    print("  python3 sample/sr_sod/scripts/compare_solvers.py")

if __name__ == "__main__":
    main()
