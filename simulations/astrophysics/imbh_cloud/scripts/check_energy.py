#!/usr/bin/env python3
"""Quick check of energy.dat values vs expectations"""

import numpy as np

# Load energy data
energy_file = "simulations/astrophysics/imbh_cloud/results/Mc1e3_Mbh1e5_b3_v10/adiabatic_61k_gsph/energy.dat"

try:
    data = np.genfromtxt(energy_file, comments='#')
    if data.ndim == 1:
        data = data.reshape(1, -1)
    
    print("Energy data from energy.dat:")
    print("-" * 60)
    
    # First few rows
    print(f"{'Time':>10} {'Kinetic':>12} {'Thermal':>10} {'Potential':>12} {'Total':>12}")
    for i in range(min(5, len(data))):
        print(f"{data[i,0]:>10.4f} {data[i,1]:>12.4f} {data[i,2]:>10.4f} {data[i,3]:>12.4f} {data[i,4]:>12.4f}")
    
    print("...")
    
    # Around perihelion (t~2 code units = 100 snapshots)
    idx_100 = min(100, len(data)-1)
    print(f"\nAt index ~{idx_100} (should be t~2):")
    print(f"{data[idx_100,0]:>10.4f} {data[idx_100,1]:>12.4f} {data[idx_100,2]:>10.4f} {data[idx_100,3]:>12.4f} {data[idx_100,4]:>12.4f}")
    
    # Check conservation
    print(f"\n{'Energy Conservation Check':^60}")
    print("-" * 60)
    print(f"Initial total energy: {data[0,4]:.6f}")
    print(f"Final total energy: {data[-1,4]:.6f}")
    print(f"Difference: {abs(data[-1,4] - data[0,4]):.6f}")
    print(f"Relative error: {abs(data[-1,4] - data[0,4]) / abs(data[0,4]) * 100:.4f}%")
    
    # Expected BH contribution at t=0
    print(f"\n{'Expected BH Contribution at t=0':^60}")
    print("-" * 60)
    r_0 = np.sqrt(20**2 + 3**2)  # Initial distance
    M_BH = 100  # code units
    U_BH_0 = -M_BH / r_0
    print(f"Initial distance to BH: {r_0:.2f} pc")
    print(f"Expected BH potential: {U_BH_0:.2f} code units")
    print(f"Actual initial potential: {data[0,3]:.2f} code units")
    print(f"Difference: {data[0,3] - U_BH_0:.2f} (self-gravity contribution)")
    
    if data[0,3] > -8:  # If PE is less negative than expected with BH
        print("\n⚠️  WARNING: Initial potential may be MISSING the BH contribution!")
        print(f"   Expected ~{-5.8 + U_BH_0:.1f}, got {data[0,3]:.1f}")
    else:
        print("\n✓ Initial potential includes BH contribution")
        
except Exception as e:
    print(f"Error: {e}")
