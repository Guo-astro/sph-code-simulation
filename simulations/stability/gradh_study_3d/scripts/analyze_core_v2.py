#!/usr/bin/env python3
"""
Analyze the CORE behavior - always look at particles nearest to center at each snapshot.
This captures what happens at r~0 regardless of which particle is there.
"""
import pandas as pd
import numpy as np
import os

GRADH_DIR = "/Users/guo/Downloads/sphcode/simulations/stability/gradh_study_3d/results/kernel_test/gradh_hk"
NOGRADH_DIR = "/Users/guo/Downloads/sphcode/simulations/stability/gradh_study_3d/results/kernel_test/no_gradh_hk"

def read_snapshot(directory, snap_num):
    """Read a snapshot file in CSV format."""
    fname = os.path.join(directory, f"snapshot_{snap_num:04d}.csv")
    if not os.path.exists(fname):
        return None
    
    # Read CSV, skipping comment lines
    df = pd.read_csv(fname, comment='#')
    df['r'] = np.sqrt(df['pos_x']**2 + df['pos_y']**2 + df['pos_z']**2)
    # Rename columns for consistency
    df = df.rename(columns={'dens': 'rho', 'pres': 'P', 'sml': 'h', 'gradh': 'omega'})
    return df

def get_core_particles(df, n=10):
    """Get the n particles closest to center."""
    return df.nsmallest(n, 'r')

print("=" * 80)
print("CORE ANALYSIS: What happens at r~0 over time")
print("=" * 80)

# Check how many snapshots exist
import glob
gradh_files = sorted(glob.glob(os.path.join(GRADH_DIR, "snapshot_*.csv")))
nogradh_files = sorted(glob.glob(os.path.join(NOGRADH_DIR, "snapshot_*.csv")))
print(f"\nFound {len(gradh_files)} grad-h snapshots, {len(nogradh_files)} no-gradh snapshots")

# Analyze key snapshots
snapshots = [0, 10, 20, 50, 100, 150, 200, 250, 300]

print("\n" + "=" * 80)
print("CORE DENSITY EVOLUTION (average of 10 innermost particles)")
print("=" * 80)

print(f"\n{'snap':>4}  {'r_max_gradh':>12}  {'ρ_gradh':>10}  {'Ω_gradh':>10}  {'r_max_nogr':>12}  {'ρ_nogradh':>10}")
print("-" * 80)

data_gradh = []
data_nogradh = []

for snap in snapshots:
    df_gradh = read_snapshot(GRADH_DIR, snap)
    df_nogradh = read_snapshot(NOGRADH_DIR, snap)
    
    if df_gradh is None:
        print(f"  Missing gradh snapshot {snap}")
        continue
    if df_nogradh is None:
        print(f"  Missing no-gradh snapshot {snap}")
        continue
    
    core_gradh = get_core_particles(df_gradh, 10)
    core_nogradh = get_core_particles(df_nogradh, 10)
    
    rho_g = core_gradh['rho'].mean()
    omega_g = core_gradh['omega'].mean()
    r_max_g = core_gradh['r'].max()
    h_g = core_gradh['h'].mean()
    
    rho_ng = core_nogradh['rho'].mean()
    r_max_ng = core_nogradh['r'].max()
    h_ng = core_nogradh['h'].mean()
    
    print(f"{snap:4d}  {r_max_g:12.6f}  {rho_g:10.4f}  {omega_g:10.4f}  {r_max_ng:12.6f}  {rho_ng:10.4f}")
    
    data_gradh.append({'snap': snap, 'rho': rho_g, 'omega': omega_g, 'r_max': r_max_g, 'h': h_g})
    data_nogradh.append({'snap': snap, 'rho': rho_ng, 'r_max': r_max_ng, 'h': h_ng})

if len(data_gradh) == 0:
    print("ERROR: No data found!")
    exit(1)

print("\n" + "=" * 80)
print("KEY FINDING: COMPARING DENSITY COMPRESSION RATIO")
print("=" * 80)

# Initial values
rho0_g = data_gradh[0]['rho']
rho0_ng = data_nogradh[0]['rho']
omega0 = data_gradh[0]['omega']
h0_g = data_gradh[0]['h']

print(f"\nInitial core density: ρ₀ = {rho0_g:.4f}")
print(f"Initial Ω at core: Ω₀ = {omega0:.4f}")
print(f"Initial h at core: h₀ = {h0_g:.6f}")

print(f"\n{'snap':>4}  {'(ρ/ρ₀)_gradh':>12}  {'(ρ/ρ₀)_nogradh':>14}  {'Ratio':>8}  {'Ω':>8}  {'h/h₀':>8}")
print("-" * 70)

for dg, dng in zip(data_gradh, data_nogradh):
    ratio_g = dg['rho'] / rho0_g
    ratio_ng = dng['rho'] / rho0_ng
    ratio = ratio_ng / ratio_g if ratio_g > 0 else np.inf
    h_ratio = dg['h'] / h0_g
    print(f"{dg['snap']:4d}  {ratio_g:12.4f}  {ratio_ng:14.4f}  {ratio:8.2f}x  {dg['omega']:8.4f}  {h_ratio:8.4f}")

print("\n" + "=" * 80)
print("Ω BEHAVIOR ANALYSIS")
print("=" * 80)

print("""
Ω is defined as: Ω = 1 / (1 + (h/Dρ) dρ/dh)

For adiabatic change: ρ ∝ h^(-D) → dρ/dh = -D·ρ/h
So: (h/Dρ) dρ/dh = -1
And: Ω = 1 / (1 - 1) → undefined!

But this assumes UNIFORM compression. Reality is different:
- Core compresses more than outer regions
- h-ρ relationship deviates from ρ ∝ h^(-3)
""")

print("\n--- Checking h vs ρ relationship ---")
print(f"\n{'snap':>4}  {'ρ/ρ₀':>10}  {'h/h₀':>10}  {'theory h/h₀':>12}  {'deviation':>10}")
print("-" * 60)

for dg in data_gradh:
    rho_ratio = dg['rho'] / rho0_g
    h_ratio = dg['h'] / h0_g
    h_theory = rho_ratio ** (-1/3)  # From ρ ∝ h^(-3)
    deviation = (h_ratio - h_theory) / h_theory * 100
    print(f"{dg['snap']:4d}  {rho_ratio:10.4f}  {h_ratio:10.4f}  {h_theory:12.4f}  {deviation:+10.2f}%")

print("\n" + "=" * 80)
print("FORCE ANALYSIS: WHY Ω STABILIZES")
print("=" * 80)

print("""
The Standard SPH force is:
  F_std = m·(P_i/ρ_i²·∇W_i + P_j/ρ_j²·∇W_j)

For uniform perturbation δρ:
  δF_std/F_std = δρ/ρ  (coefficient = 1)

The GSPH force with grad-h is:
  F_gsph = m·p*·(Ω_i/ρ_i²·∇W_i + Ω_j/ρ_j²·∇W_j)

The shared p* changes response, but Ω provides additional correction.
""")

# Track force-like quantity F ∝ P·Ω·|∇W|/ρ² = P·Ω·h^(-4)/ρ²
print("\n--- Force-like quantity F ∝ P·Ω/(ρ²·h⁴) ---")
print(f"\n{'snap':>4}  {'P':>10}  {'ρ':>10}  {'h':>10}  {'Ω':>8}  {'F/F₀':>10}  {'P/P₀':>10}  {'ρ/ρ₀':>10}")
print("-" * 90)

P0_g = None
for snap in snapshots:
    df = read_snapshot(GRADH_DIR, snap)
    if df is None:
        continue
    core = get_core_particles(df, 10)
    
    P = core['P'].mean()
    rho = core['rho'].mean()
    h = core['h'].mean()
    omega = core['omega'].mean()
    
    if P0_g is None:
        P0_g = P
        rho0 = rho
        h0 = h
        F0 = P * omega / (rho**2 * h**4)
    
    F = P * omega / (rho**2 * h**4)
    F_ratio = F / F0
    P_ratio = P / P0_g
    rho_ratio = rho / rho0
    
    print(f"{snap:4d}  {P:10.4f}  {rho:10.4f}  {h:10.6f}  {omega:8.4f}  {F_ratio:10.4f}  {P_ratio:10.4f}  {rho_ratio:10.4f}")

print("\n" + "=" * 80)
print("CRITICAL COMPARISON: FORCE RESPONSE TO DENSITY")
print("=" * 80)

print("""
For stability against gravity, we need:
  δF_pressure / F_pressure ≥ δF_gravity / F_gravity

Gravity gives: δF_g/F_g = δρ/ρ (linear response)

So we need pressure force to respond at least linearly to density.
Let's calculate the ACTUAL response coefficient:
""")

print("\n--- Calculating δF/F vs δρ/ρ ---")
print(f"\n{'snap':>4}  {'δρ/ρ₀':>10}  {'δF/F₀':>10}  {'δF/F / δρ/ρ':>12}  Comment")
print("-" * 70)

for i, snap in enumerate(snapshots[1:], 1):  # Skip first
    df = read_snapshot(GRADH_DIR, snap)
    if df is None:
        continue
    core = get_core_particles(df, 10)
    
    P = core['P'].mean()
    rho = core['rho'].mean()
    h = core['h'].mean()
    omega = core['omega'].mean()
    
    F = P * omega / (rho**2 * h**4)
    F_ratio = F / F0
    delta_F = F_ratio - 1
    delta_rho = rho / rho0 - 1
    
    if abs(delta_rho) > 0.01:
        coeff = delta_F / delta_rho
        if coeff >= 0.9:
            comment = "STABLE ✓"
        else:
            comment = f"weak response"
    else:
        coeff = float('nan')
        comment = "small δρ"
    
    print(f"{snap:4d}  {delta_rho:+10.4f}  {delta_F:+10.4f}  {coeff:12.4f}  {comment}")

print("\n" + "=" * 80)
print("NO GRAD-H COMPARISON")
print("=" * 80)

P0_ng = None
print(f"\n{'snap':>4}  {'P':>10}  {'ρ':>10}  {'h':>10}  {'F/F₀':>10}  {'P/P₀':>10}  {'ρ/ρ₀':>10}  {'δF/δρ coeff':>12}")
print("-" * 100)

for snap in snapshots:
    df = read_snapshot(NOGRADH_DIR, snap)
    if df is None:
        continue
    core = get_core_particles(df, 10)
    
    P = core['P'].mean()
    rho = core['rho'].mean()
    h = core['h'].mean()
    omega = 1.0  # No grad-h means Ω = 1
    
    if P0_ng is None:
        P0_ng = P
        rho0_ng = rho
        h0_ng = h
        F0_ng = P * omega / (rho**2 * h**4)
    
    F = P * omega / (rho**2 * h**4)
    F_ratio = F / F0_ng
    P_ratio = P / P0_ng
    rho_ratio = rho / rho0_ng
    
    delta_F = F_ratio - 1
    delta_rho = rho_ratio - 1
    if abs(delta_rho) > 0.01:
        coeff = delta_F / delta_rho
    else:
        coeff = float('nan')
    
    print(f"{snap:4d}  {P:10.4f}  {rho:10.4f}  {h:10.6f}  {F_ratio:10.4f}  {P_ratio:10.4f}  {rho_ratio:10.4f}  {coeff:12.4f}")

print("\n" + "=" * 80)
print("CONCLUSION")
print("=" * 80)

print("""
The Ω factor adjusts the force response to density changes:

WITHOUT Ω (no grad-h):
- Force F ∝ P/(ρ²h⁴) with Ω = 1
- When density increases, h decreases (h ∝ ρ^(-1/3))
- |∇W| ∝ h^(-4) increases rapidly (∝ ρ^(4/3))
- But in GSPH with shared p*, the response is mismatched
- Result: Positive feedback → COLLAPSE

WITH Ω (grad-h enabled):
- Ω adjusts based on local gradient ∂ρ/∂h
- This modification restores proper coupling
- Derived from variational principle → energy conserving
- Result: Stable oscillation around equilibrium
""")
