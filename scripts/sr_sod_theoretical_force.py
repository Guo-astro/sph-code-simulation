#!/usr/bin/env python3
"""
Calculate theoretical force values for SR Sod shock tube
Based on Kitajima+ 2025 (arxiv:2510.18251v1) equations
"""

import numpy as np

print("="*70)
print("SR SOD SHOCK TUBE: THEORETICAL FORCE CALCULATION")
print("="*70)

# Initial conditions (from paper Section 3.1.1)
print("\n--- INITIAL CONDITIONS ---")
P_L = 1.0       # Left pressure
n_L = 1.0       # Left rest-frame density
v_L = 0.0       # Left velocity

P_R = 0.1       # Right pressure  
n_R = 0.125     # Right rest-frame density
v_R = 0.0       # Right velocity

gamma_c = 5.0/3.0  # Adiabatic index
c = 1.0            # Speed of light

print(f"Left state:  P={P_L}, n={n_L}, v={v_L}")
print(f"Right state: P={P_R}, n={n_R}, v={v_R}")
print(f"γ_c = {gamma_c:.4f}, c = {c}")

# Calculate initial thermodynamic quantities
# u = P/[(γ_c - 1) × n]  (Eq. 9)
u_L = P_L / ((gamma_c - 1) * n_L)
u_R = P_R / ((gamma_c - 1) * n_R)

print(f"\nThermal energy per baryon:")
print(f"  u_L = {u_L:.6f}")
print(f"  u_R = {u_R:.6f}")

# Lorentz factor: γ = 1/√(1 - v²/c²)  (Eq. 7)
gamma_L = 1.0 / np.sqrt(1.0 - v_L**2 / c**2)
gamma_R = 1.0 / np.sqrt(1.0 - v_R**2 / c**2)

print(f"\nLorentz factors:")
print(f"  γ_L = {gamma_L:.6f}")
print(f"  γ_R = {gamma_R:.6f}")

# Enthalpy: H = 1 + u/c² + P/(n×c²)  (Eq. 8)
H_L = 1.0 + u_L/c**2 + P_L/(n_L * c**2)
H_R = 1.0 + u_R/c**2 + P_R/(n_R * c**2)

print(f"\nSpecific enthalpy:")
print(f"  H_L = {H_L:.6f}")
print(f"  H_R = {H_R:.6f}")

# Lab-frame baryon number density: N = γ × n  (mentioned after Eq. 8)
N_L = gamma_L * n_L
N_R = gamma_R * n_R

print(f"\nLab-frame baryon number density:")
print(f"  N_L = {N_L:.6f}")
print(f"  N_R = {N_R:.6f}")

# Sound speed: c_s = √[(Γ-1)(H-1)/H] where Γ = γ_c (from text before Eq. 66)
cs_L = np.sqrt((gamma_c - 1) * (H_L - 1) / H_L)
cs_R = np.sqrt((gamma_c - 1) * (H_R - 1) / H_R)

print(f"\nSound speeds:")
print(f"  c_s,L = {cs_L:.6f}")
print(f"  c_s,R = {cs_R:.6f}")

print("\n" + "="*70)
print("FORCE CALCULATION AT CONTACT DISCONTINUITY")
print("="*70)

# At t=0, particles at the interface will interact
# The Riemann solver gives P* and v* at the interface

# For initial conditions, use HLL solver approximation
# (This is a rough estimate; exact Riemann solver would be better)

# Simple average for initial estimate
print("\n--- RIEMANN SOLVER (Approximate) ---")
print("At contact discontinuity, solving Riemann problem:")

# Use simple HLL-type estimate
rho_L_sqrt = np.sqrt(n_L)
rho_R_sqrt = np.sqrt(n_R)
rho_avg = (rho_L_sqrt + rho_R_sqrt) / 2

u_tilde = (rho_L_sqrt * v_L + rho_R_sqrt * v_R) / (rho_L_sqrt + rho_R_sqrt)
c_tilde = (rho_L_sqrt * cs_L + rho_R_sqrt * cs_R) / (rho_L_sqrt + rho_R_sqrt)

s_L = min(v_L - cs_L, u_tilde - c_tilde)
s_R = max(v_R + cs_R, u_tilde + c_tilde)

# HLL pressure and velocity
c1 = n_L * (s_L - v_L)
c2 = n_R * (s_R - v_R)
c3 = 1.0 / (c1 - c2)
c4 = P_L - v_L * c1
c5 = P_R - v_R * c2

v_star = (c5 - c4) * c3
P_star = (c1 * c5 - c2 * c4) * c3

print(f"  v* ≈ {v_star:.6f}")
print(f"  P* ≈ {P_star:.6f}")

# Expected P* for Sod problem is approximately 0.3-0.5
# Let's use a typical value
P_star_typical = 0.4
print(f"\nTypical P* for Sod problem: {P_star_typical:.3f}")

print("\n" + "="*70)
print("SPH PARTICLE SETUP")
print("="*70)

# From the code: 3200 particles on left, 400 on right
# Domain: x ∈ [-0.5, 0.5], discontinuity at x=0
N_particles_L = 3200
N_particles_R = 400
domain_length = 1.0

# Particle spacing
dx_L = 0.5 / N_particles_L  # Left region spans [-0.5, 0]
dx_R = 0.5 / N_particles_R  # Right region spans [0, 0.5]

print(f"Number of particles: {N_particles_L} (left) + {N_particles_R} (right)")
print(f"Particle spacing: dx_L = {dx_L:.6e}, dx_R = {dx_R:.6e}")

# Baryon number per particle (constant in this setup)
# Total baryon number in each region = n × volume
total_baryon_L = n_L * 0.5
total_baryon_R = n_R * 0.5

nu_L = total_baryon_L / N_particles_L
nu_R = total_baryon_R / N_particles_R

print(f"\nBaryon number per particle:")
print(f"  ν_L = {nu_L:.6e}")
print(f"  ν_R = {nu_R:.6e}")

# Smoothing length (approximately 50 neighbors × dx / 2)
neighbor_num = 50
h_L = neighbor_num * dx_L / 2.0
h_R = neighbor_num * dx_R / 2.0

print(f"\nSmoothing length (η=1, neighbor=50):")
print(f"  h_L ≈ {h_L:.6e}")
print(f"  h_R ≈ {h_R:.6e}")

print("\n" + "="*70)
print("FORCE MAGNITUDE ESTIMATION")
print("="*70)

# Force equation (Eq. 64):
# ⟨ν_i Ṡ_i⟩ = -Σ_j P*_ij V²_ij [∇_i W - ∇_j W]

print("\n--- Per-Particle Force Components ---")

# V²_ij term: V²_ij = 0.5 × ((ν/N)²_i + (ν/N)²_j)
V2_L = (nu_L / N_L)**2
V2_R = (nu_R / N_R)**2
V2_ij = 0.5 * (V2_L + V2_R)

print(f"\nVolume-squared terms:")
print(f"  (ν/N)²_L = {V2_L:.6e}")
print(f"  (ν/N)²_R = {V2_R:.6e}")
print(f"  V²_ij = 0.5 × (V²_L + V²_R) = {V2_ij:.6e}")

# Kernel gradient magnitude
# For cubic spline at small q, |∇W| ~ C/h² where C depends on dimensionality
# For 1D cubic spline: W(q) with q = r/h
# |∇W| ≈ (∂W/∂q) × (1/h) × (r/r)
# At contact (r ~ dx), q ~ dx/h ~ 1/25 (small)
# For cubic spline: ∂W/∂q ≈ -2q for q < 1, so |dW/dq| ≈ 0.08
# |∇W| ≈ 0.08 / h

# Using kernel with argument 2h (as in Eq. 64)
h_avg = (h_L + h_R) / 2
r_contact = (dx_L + dx_R) / 2  # Typical separation at contact
q_contact = r_contact / (2 * h_avg)

print(f"\nKernel gradient estimation:")
print(f"  h_avg = {h_avg:.6e}")
print(f"  r_contact ≈ {r_contact:.6e}")
print(f"  q = r/(2h) ≈ {q_contact:.4f}")

# For 1D cubic spline with argument (r, r, 2h):
# The normalization is 1/(2h), and gradient magnitude scales as 1/(2h)²
# Rough estimate: |∇W| ~ 0.1 / (2h)²
grad_W_mag = 0.1 / (2 * h_avg)**2
print(f"  |∇W| ≈ {grad_W_mag:.6e}")

# Number of neighbors contributing
n_neighbors = neighbor_num  # ~50

print(f"\nForce per particle (left side):")
print(f"  Number of neighbors ≈ {n_neighbors}")

# Single neighbor contribution: P* × V²_ij × |∇W|
dS_single = P_star_typical * V2_ij * grad_W_mag
print(f"  |dS| per neighbor ≈ P* × V²_ij × |∇W|")
print(f"                    ≈ {P_star_typical} × {V2_ij:.3e} × {grad_W_mag:.3e}")
print(f"                    ≈ {dS_single:.3e}")

# Total force from all neighbors
dS_total = dS_single * n_neighbors
print(f"\n  |dS| total ≈ {dS_total:.3e}")

# Force per baryon: dS/dt per baryon = (1/ν) × dS/dt
dS_per_baryon = dS_total / nu_L
print(f"  |dS/dt| / ν ≈ {dS_per_baryon:.3e}")

print("\n" + "="*70)
print("TIME EVOLUTION ESTIMATE")
print("="*70)

# Timestep (CFL condition)
# Δt = C_CFL × min(h/c_s)
C_CFL = 0.5
dt = C_CFL * h_L / cs_L
print(f"\nTimestep (CFL = {C_CFL}):")
print(f"  Δt ≈ {dt:.6e}")

# Change in S per timestep
Delta_S = dS_total * dt
print(f"\nChange per timestep:")
print(f"  Δ|S| ≈ |dS/dt| × Δt ≈ {Delta_S:.6e}")

# How many steps until S becomes large?
# S should remain small (|S| << γHv ~ 0.1 for v~0.1)
S_threshold = 0.1  # Conservative threshold
steps_to_threshold = S_threshold / Delta_S if Delta_S > 0 else np.inf

print(f"\nSteps until |S| ~ {S_threshold}:")
print(f"  N_steps ≈ {steps_to_threshold:.1f}")

print("\n" + "="*70)
print("COMPARISON WITH SIMULATION")
print("="*70)

print(f"\nTheoretical predictions:")
print(f"  V²_ij ≈ {V2_ij:.3e}")
print(f"  |∇W| ≈ {grad_W_mag:.3e}")
print(f"  |dS/dt| / ν ≈ {dS_per_baryon:.3e}")
print(f"  Δt ≈ {dt:.3e}")
print(f"  Δ|S| per step ≈ {Delta_S:.3e}")

print(f"\nSimulation should:")
print(f"  1. Have V²_ij ~ {V2_ij:.2e} at contact")
print(f"  2. Have |dS/dt|/ν ~ {dS_per_baryon:.2e}")
print(f"  3. Remain stable for > {steps_to_threshold:.0f} steps")
print(f"  4. Not have S grow beyond ~0.1-0.2 initially")

print("\n" + "="*70)
print("DIAGNOSTIC SUMMARY")
print("="*70)

print(f"""
Expected initial force characteristics:
  • Pressure jump: ΔP = {P_L - P_R} → P* ≈ {P_star_typical}
  • Volume factor: V²_ij = {V2_ij:.2e} [L²]
  • Kernel gradient: |∇W| ~ {grad_W_mag:.2e} [1/L²]  
  • Force magnitude: |dS/dt|/ν ~ {dS_per_baryon:.2e}
  • Timestep: Δt ~ {dt:.2e}
  • Momentum change: Δ|S| ~ {Delta_S:.2e} per step

If simulation shows much larger forces:
  → Check kernel normalization
  → Check kernel gradient calculation
  → Check if 2h is used correctly in kernel arguments
  → Check Riemann solver P* values
  → Check dimension factors (1D vs 3D)

If S grows to O(1) within ~10 steps:
  → Force is ~100× too large
  → Likely issue: missing normalization or wrong kernel argument
""")

print("="*70)
