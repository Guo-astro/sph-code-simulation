#!/usr/bin/env python3
"""
Comprehensive SR-GSPH diagnostic comparing ALL physics values with theory.
Checks: sound speed, enthalpy, Lorentz factor, thermodynamic quantities, forces.
"""

import numpy as np
import sys

# Physical constants
c = 1.0  # Speed of light (in code units)
gamma_c = 5.0/3.0  # Adiabatic index

print("="*80)
print("COMPREHENSIVE SR-SOD DIAGNOSTIC: THEORY vs SIMULATION")
print("="*80)

# ============================================================================
# Initial Conditions (from paper Section 3.1.1)
# ============================================================================
print("\n" + "="*80)
print("INITIAL CONDITIONS (from paper)")
print("="*80)

# Left state
P_L = 1.0
n_L = 1.0  # Rest-frame density
v_L = 0.0

# Right state
P_R = 0.1
n_R = 0.125
v_R = 0.0

print(f"Left state:  P = {P_L:.3f}, n = {n_L:.3f}, v = {v_L:.3f}")
print(f"Right state: P = {P_R:.3f}, n = {n_R:.3f}, v = {v_R:.3f}")
print(f"Adiabatic index: γ_c = {gamma_c:.4f}")

# ============================================================================
# Thermodynamic Quantities (Equations 8, 9)
# ============================================================================
print("\n" + "="*80)
print("THERMODYNAMIC QUANTITIES (Eqs. 8, 9)")
print("="*80)

# Internal energy per baryon: P = (γ_c - 1) * n * u
u_L = P_L / ((gamma_c - 1.0) * n_L)
u_R = P_R / ((gamma_c - 1.0) * n_R)

print(f"\nInternal energy per baryon u (Eq. 9):")
print(f"  Left:  u = {u_L:.6f}")
print(f"  Right: u = {u_R:.6f}")

# Enthalpy per baryon: H = 1 + u/c² + P/(n*c²)
H_L = 1.0 + u_L/(c**2) + P_L/(n_L * c**2)
H_R = 1.0 + u_R/(c**2) + P_R/(n_R * c**2)

print(f"\nEnthalpy per baryon H (Eq. 8):")
print(f"  Left:  H = {H_L:.6f}")
print(f"  Right: H = {H_R:.6f}")

# Lorentz factor (Eq. 7) - should be 1.0 for v=0
gamma_L = 1.0 / np.sqrt(1.0 - v_L**2/c**2)
gamma_R = 1.0 / np.sqrt(1.0 - v_R**2/c**2)

print(f"\nLorentz factor γ (Eq. 7):")
print(f"  Left:  γ = {gamma_L:.6f}")
print(f"  Right: γ = {gamma_R:.6f}")

# Lab-frame density N = γ * n
N_L = gamma_L * n_L
N_R = gamma_R * n_R

print(f"\nLab-frame density N = γ*n:")
print(f"  Left:  N = {N_L:.6f}")
print(f"  Right: N = {N_R:.6f}")

# ============================================================================
# SOUND SPEED (Equation 66 - CRITICAL!)
# ============================================================================
print("\n" + "="*80)
print("SOUND SPEED (Eq. 66 - CRITICAL CHECK)")
print("="*80)

# From paper Eq. 66: c_s = sqrt[(γ_c - 1)(H - 1) / H]
# This is the CORRECT formula!

X = gamma_c / (gamma_c - 1.0)
print(f"\nX = γ_c/(γ_c-1) = {X:.6f}")

cs_L_correct = np.sqrt((gamma_c - 1.0) * (H_L - 1.0) / H_L)
cs_R_correct = np.sqrt((gamma_c - 1.0) * (H_R - 1.0) / H_R)

print(f"\nCORRECT sound speed c_s (Eq. 66):")
print(f"  c_s = sqrt[(γ_c - 1)(H - 1) / H]")
print(f"  Left:  c_s = {cs_L_correct:.6f}")
print(f"  Right: c_s = {cs_R_correct:.6f}")

# Common WRONG formula (for comparison)
cs_L_wrong = np.sqrt(gamma_c * P_L / (n_L * H_L))
cs_R_wrong = np.sqrt(gamma_c * P_R / (n_R * H_R))

print(f"\nWRONG sound speed (common mistake):")
print(f"  c_s = sqrt[γ_c * P / (n * H)]")
print(f"  Left:  c_s = {cs_L_wrong:.6f}  (WRONG!)")
print(f"  Right: c_s = {cs_R_wrong:.6f}  (WRONG!)")

print(f"\nRatio (correct/wrong):")
print(f"  Left:  {cs_L_correct/cs_L_wrong:.6f}")
print(f"  Right: {cs_R_correct/cs_R_wrong:.6f}")

# ============================================================================
# RIEMANN SOLVER (approximate)
# ============================================================================
print("\n" + "="*80)
print("RIEMANN SOLVER (approximate)")
print("="*80)

# Use simple average for P* and v* (exact solver would be better)
P_star = 0.5 * (P_L + P_R)
v_star = 0.5 * (v_L + v_R)

print(f"Approximate Riemann solution:")
print(f"  P* ≈ {P_star:.6f}")
print(f"  v* ≈ {v_star:.6f}")

# ============================================================================
# SPH PARTICLE SETUP
# ============================================================================
print("\n" + "="*80)
print("SPH PARTICLE SETUP")
print("="*80)

N_left = 3200
N_right = 400
L = 1.0  # Domain length

# Particle spacing
dx_left = L/2.0 / N_left
dx_right = L/2.0 / N_right

print(f"Number of particles:")
print(f"  Left:  {N_left}")
print(f"  Right: {N_right}")
print(f"Particle spacing:")
print(f"  Left:  dx = {dx_left:.6e}")
print(f"  Right: dx = {dx_right:.6e}")

# Baryon number per particle (constant)
nu = n_L * dx_left  # Left particles
print(f"\nBaryon number per particle:")
print(f"  ν = n_L * dx_L = {nu:.6e}")

# Particle volume (for GSPH)
Vp_L = dx_left  # In 1D, volume = spacing
Vp_R = dx_right

print(f"\nParticle volume V_p:")
print(f"  Left:  V_p = {Vp_L:.6e}")
print(f"  Right: V_p = {Vp_R:.6e}")

# Verify N = ν / V_p
N_L_check = nu / Vp_L
N_R_check = nu / Vp_R

print(f"\nVerify N = ν/V_p:")
print(f"  Left:  N = {N_L_check:.6f} (should be {N_L:.6f})")
print(f"  Right: N = {N_R_check:.6f} (should be {N_R:.6f})")

# ============================================================================
# KERNEL QUANTITIES
# ============================================================================
print("\n" + "="*80)
print("KERNEL QUANTITIES (1D Gaussian)")
print("="*80)

# Smoothing length
h = 2.0 * dx_left
print(f"Smoothing length h = 2*dx_L = {h:.6e}")

# Kernel normalization (1D Gaussian)
# W(x, h) = (1/(h*sqrt(π))) * exp(-x²/h²)
norm_1D = 1.0 / (h * np.sqrt(np.pi))
print(f"Kernel normalization (1D): 1/(h√π) = {norm_1D:.6e}")

# At contact (x = 0), typical particle separation ~ dx_left
r_contact = dx_left
r_scaled = r_contact / h

print(f"\nAt contact discontinuity:")
print(f"  Separation r ≈ {r_contact:.6e}")
print(f"  Scaled r/h ≈ {r_scaled:.6f}")

# Kernel gradient magnitude (approximate)
# dW/dx ≈ -(2x/h²) * W(x,h)
# At r ~ dx: |dW/dx| ~ (2*dx/h²) * norm
kernel_val = norm_1D * np.exp(-(r_contact/h)**2)
dW_dx = (2.0 * r_contact / h**2) * kernel_val

print(f"  Kernel W(r, h) ≈ {kernel_val:.6e}")
print(f"  |dW/dx| ≈ {dW_dx:.6e}")

# For argument 2h (as used in code)
norm_2h = 1.0 / (2.0*h * np.sqrt(np.pi))
kernel_val_2h = norm_2h * np.exp(-(r_contact/(2.0*h))**2)
dW_dx_2h = (2.0 * r_contact / (2.0*h)**2) * kernel_val_2h

print(f"\nWith 2h argument (used in force calculation):")
print(f"  Kernel W(r, 2h) ≈ {kernel_val_2h:.6e}")
print(f"  |dW/dx|(2h) ≈ {dW_dx_2h:.6e}")

# ============================================================================
# V²_ij CALCULATION (Eq. 29, 63)
# ============================================================================
print("\n" + "="*80)
print("V²_ij CALCULATION (Eqs. 29, 63)")
print("="*80)

# At contact, we have one left particle and one right particle
# V²_ij = 0.5 * ((ν/N_i)² + (ν/N_j)²)

Vsq_i = (nu / N_L)**2
Vsq_j = (nu / N_R)**2
Vsq_ij = 0.5 * (Vsq_i + Vsq_j)

print(f"V² components:")
print(f"  (ν/N_L)² = {Vsq_i:.6e}")
print(f"  (ν/N_R)² = {Vsq_j:.6e}")
print(f"  V²_ij = 0.5*(V²_i + V²_j) = {Vsq_ij:.6e}")

print(f"\nDimensional check:")
print(f"  ν has units: [baryon number]")
print(f"  N has units: [baryon number/length]")
print(f"  ν/N has units: [length]")
print(f"  (ν/N)² has units: [length²]  ✓")

# ============================================================================
# FORCE CALCULATION (Eq. 64)
# ============================================================================
print("\n" + "="*80)
print("FORCE CALCULATION (Eq. 64)")
print("="*80)

# ⟨ν_i * dS/dt⟩ = -Σ_j P*_ij * V²_ij * [∇_i W - ∇_j W]
# For single pair at contact: dS/dt ≈ -(P_star * V²_ij * 2 * |dW/dx|) / ν

# Using 2h kernel gradient (as in code)
force_per_nu = P_star * Vsq_ij * 2.0 * dW_dx_2h

print(f"Force calculation (Eq. 64):")
print(f"  P* = {P_star:.6f}")
print(f"  V²_ij = {Vsq_ij:.6e}")
print(f"  2*|∇W|(2h) = {2.0*dW_dx_2h:.6e}")
print(f"  |dS/dt|/ν ≈ P* * V²_ij * 2|∇W| = {force_per_nu:.6e}")

# ============================================================================
# TIMESTEP AND EVOLUTION
# ============================================================================
print("\n" + "="*80)
print("TIMESTEP AND EVOLUTION")
print("="*80)

# CFL condition: Δt = C_CFL * h / c_s
C_CFL = 0.3
dt = C_CFL * h / cs_L_correct

print(f"CFL condition (Eq. 73):")
print(f"  C_CFL = {C_CFL}")
print(f"  c_s (left) = {cs_L_correct:.6f}")
print(f"  Δt = C_CFL * h / c_s = {dt:.6e}")

# Expected change in S per timestep
dS_per_step = force_per_nu * dt

print(f"\nExpected momentum change per timestep:")
print(f"  Δ|S| = (|dS/dt|/ν) * Δt = {dS_per_step:.6e}")

# How many steps to reach S ~ 0.1?
S_threshold = 0.1
steps_to_threshold = S_threshold / dS_per_step

print(f"\nSteps to reach |S| ~ {S_threshold}:")
print(f"  N_steps ≈ {steps_to_threshold:.1f}")

# ============================================================================
# SUMMARY OF EXPECTED VALUES
# ============================================================================
print("\n" + "="*80)
print("SUMMARY: EXPECTED VALUES IN SIMULATION")
print("="*80)

print(f"\nThermodynamic quantities:")
print(f"  H_L = {H_L:.6f},  H_R = {H_R:.6f}")
print(f"  c_s(L) = {cs_L_correct:.6f},  c_s(R) = {cs_R_correct:.6f}  (Eq. 66!)")
print(f"  N_L = {N_L:.6f},  N_R = {N_R:.6f}")

print(f"\nSPH quantities:")
print(f"  ν = {nu:.6e}")
print(f"  V²_ij ≈ {Vsq_ij:.6e} [L²]")
print(f"  |∇W|(2h) ≈ {dW_dx_2h:.6e} [1/L²]")

print(f"\nForce and evolution:")
print(f"  |dS/dt|/ν ≈ {force_per_nu:.6e}")
print(f"  Δt ≈ {dt:.6e}")
print(f"  Δ|S| per step ≈ {dS_per_step:.6e}")

print(f"\nStability prediction:")
print(f"  Simulation should remain stable for > {int(steps_to_threshold)} steps")
print(f"  S should grow VERY SLOWLY (Δ|S| ~ {dS_per_step:.2e} per step)")

# ============================================================================
# DIAGNOSTIC CHECKS TO PERFORM
# ============================================================================
print("\n" + "="*80)
print("DIAGNOSTIC CHECKS FOR SIMULATION CODE")
print("="*80)

print(f"\n1. CHECK SOUND SPEED (CRITICAL!):")
print(f"   Expected: c_s = sqrt[(γ_c-1)(H-1)/H] = {cs_L_correct:.6f}")
print(f"   If code uses: c_s = sqrt[γ_c*P/(n*H)] = {cs_L_wrong:.6f}  ← WRONG!")
print(f"   Ratio: {cs_L_correct/cs_L_wrong:.6f}")

print(f"\n2. CHECK ENTHALPY:")
print(f"   Expected: H = 1 + u/c² + P/(n*c²) = {H_L:.6f}")

print(f"\n3. CHECK V²_ij:")
print(f"   Expected: V²_ij = 0.5*((ν/N)²_i + (ν/N)²_j) ≈ {Vsq_ij:.6e}")

print(f"\n4. CHECK KERNEL GRADIENT:")
print(f"   Expected: |∇W|(2h) ≈ {dW_dx_2h:.6e} [1/L²]")

print(f"\n5. CHECK FORCE MAGNITUDE:")
print(f"   Expected: |dS/dt|/ν ≈ {force_per_nu:.6e}")
print(f"   If observed is ~100× larger → likely sound speed error!")
print(f"   If observed is ~70,000× larger → multiple errors!")

print(f"\n6. CHECK TIMESTEP:")
print(f"   Expected: Δt ≈ {dt:.6e}")
print(f"   If Δt is much larger → CFL condition violated")

print(f"\n7. CHECK S EVOLUTION:")
print(f"   Expected: Δ|S| ≈ {dS_per_step:.6e} per step")
print(f"   If S grows to O(1) in <10 steps → force is WAY too large!")

print("\n" + "="*80)
print("END OF DIAGNOSTIC")
print("="*80)
