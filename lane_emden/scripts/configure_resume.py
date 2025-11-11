#!/usr/bin/env python3
"""
Configure minimal resume setup for Lane-Emden checkpoint resume.
This implements SSOT - physics parameters come from checkpoint, not config.
"""
import json
import sys
import os

def main():
    if len(sys.argv) < 2:
        print("Usage: configure_resume.py <checkpoint_file>", file=sys.stderr)
        sys.exit(1)
    
    checkpoint_file = sys.argv[1]
    
    # Make absolute path if relative
    if not os.path.isabs(checkpoint_file):
        checkpoint_file = os.path.abspath(checkpoint_file)
    
    # Verify checkpoint exists
    if not os.path.exists(checkpoint_file):
        print(f"ERROR: Checkpoint file not found: {checkpoint_file}", file=sys.stderr)
        sys.exit(1)
    
    # Create minimal resume configuration (SSOT compliant)
    # Only runtime control parameters - no physics!
    resume_config = {
        "outputDirectory": "lane_emden/results/polytrope_n1.5_3d",
        "startTime": 0.0,
        "endTime": 5.0,
        "outputTime": 0.05,
        "energyTime": 0.05,
        "checkpoint": {
            "enabled": True,
            "saveInterval": 1000,
            "directory": "lane_emden/results/polytrope_n1.5_3d/checkpoints",
            "saveOnExit": True,
            "maxCheckpoints": 5,
            "autoResume": False,
            "resumeFile": checkpoint_file
        },
        "output": {
            "formats": [{"type": "csv", "precision": 16}],
            "enableEnergyFile": True
        }
    }
    
    # Write to lane_emden.json
    output_path = "sample/lane_emden/lane_emden.json"
    with open(output_path, 'w') as f:
        json.dump(resume_config, f, indent=2)
    
    print(f"✓ Configured to resume from: {checkpoint_file}")
    print("✓ Using SSOT mode - physics parameters will be loaded from checkpoint")
    print(f"✓ Config written to: {output_path}")

if __name__ == "__main__":
    main()
