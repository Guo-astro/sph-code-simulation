"""
Quick test of SRGSPH with minimal particles
"""

import numpy as np
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from srgsph.simulator import SRGSPH1D

print("Testing SRGSPH with minimal setup...")

# Very small test - 20 particles
n_left = 15
n_right = 5
positions = np.concatenate([
    np.linspace(0, 0.5, n_left, endpoint=False),
    np.linspace(0.5, 1.0, n_right)
])

velocities = np.zeros(20)
densities = np.concatenate([np.ones(n_left), 0.125*np.ones(n_right)])
pressures = np.concatenate([np.ones(n_left), 0.1*np.ones(n_right)])

print(f"Particles: {len(positions)}")

# Create simulator
sim = SRGSPH1D(gamma_c=5.0/3.0, c=1.0, use_variable_h=True, use_muscl=False)

print("Setting up particles...")
sim.setup_particles(positions, velocities, densities, pressures)

print(f"Initial time: {sim.time}")
print(f"Computing timestep...")
dt = sim.compute_timestep()
print(f"Timestep: {dt:.6e}")

print("Taking one step...")
sim.step_forward(dt)
print(f"Time after step: {sim.time}")

print("Getting solution...")
sol = sim.get_solution()
print(f"Solution arrays: {sol.keys()}")
print(f"Position range: [{sol['x'].min():.4f}, {sol['x'].max():.4f}]")
print(f"Velocity range: [{sol['v'].min():.4f}, {sol['v'].max():.4f}]")

print("\n✓ Basic SRGSPH test passed!")
