#!/usr/bin/env python3
"""
Generate custom Lane-Emden configuration for 50k particle relaxation test.
"""
import json
import sys

if len(sys.argv) < 3:
    print("Usage: generate_config.py <template_json> <output_json>")
    sys.exit(1)

template_file = sys.argv[1]
output_file = sys.argv[2]

# Load template
with open(template_file, 'r') as f:
    config = json.load(f)

# Modify parameters for 50k test
config['N'] = 75  # ~50,000 particles
config['relaxationSteps'] = 1000
config['relaxationOutputFreq'] = 100
config['SPHType'] = 'gdisph'
config['kernel'] = 'wendland'
config['outputDirectory'] = 'sample/imbh_cloud/results/lane_emden_50k_relax'
config['use2ndOrderGSPH'] = True
config['useBalsaraSwitch'] = True

# Save
with open(output_file, 'w') as f:
    json.dump(config, f, indent=2)

print(f"✓ Generated config: {output_file}")
print(f"  N_shells = {config['N']} (→ ~{config['N']**3} particles)")
print(f"  Relaxation steps = {config['relaxationSteps']}")
