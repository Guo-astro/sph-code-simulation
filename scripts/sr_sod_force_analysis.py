#!/usr/bin/env python3
"""
SR Sod Shock Tube Force Analysis

Calculate expected force magnitudes for SR-GSPH and compare to simulation output.
Based on Kitajima et al. 2025 (arxiv:2510.18251v1) Equations 64-65.
"""

import numpy as np
import sys

print("=" * 80)
print("SR-GSPH Force Analysis for Sod Shock Tube")
print("=" * 80)

# ============================================================================
# Initial Conditions (from sr_sod.json)
# ============================================================================
print("\n--- Initial Conditions ---")
# Left state (x < 0)
P_L = 1.0
n_L = 1.0  # rest-frame baryon density
v_L = 0.0

# Right state (x > 0)
P_R = 0.1
n_R = 0.125
v_R = 0.0

# Particle properties
nu = 0.00015625  # baryon number per particle (constant)
N_particles_left = 3200
N_particles_right = 400
total_particles = 3600

# SPH parameters
neighbor_num = 50
gamma_eos = 5.0/3.0
c_speed = 1.0

print(f"Left region:  P = {P_L}, n = {n_L}, v = {v_L}")
print(f"Right region: P = {P_R}, n = {n_R}, v = {v_R}")
print(f"Baryon number per particle: ν = {nu}")
print(f"EOS: γ = {gamma_eos}")
print(f"Speed of light: c = {c_speed}")
print(f"Target neighbors: {neighbor_num}")

# ============================================================================
# Derived Quantities
# ============================================================================
print("\n--- Derived Quantities (Initial State) ---")

# For non-relativistic initial conditions (v ≈ 0), γ ≈ 1
gamma_L = 1.0 / np.sqrt(1 - v_L**2 / c_speed**2)
gamma_R = 1.0 / np.sqrt(1 - v_R**2 / c_speed**2)

print(f"Lorentz factors: γ_L = {gamma_L:.6f}, γ_R = {gamma_R:.6f}")

# Specific enthalpy: H = 1 + u + P/(nc²) where u = P/[(γ-1)n] for ideal gas
# For non-relativistic: u = P/[(γ-1)n]
u_L = P_L / ((gamma_eos - 1.0) * n_L)
u_R = P_R / ((gamma_eos - 1.0) * n_R)

H_L = 1.0 + u_L + P_L / (n_L * c_speed**2)
H_R = 1.0 + u_R + P_R / (n_R * c_speed**2)

print(f"Specific internal energy: u_L = {u_L:.6f}, u_R = {u_R:.6f}")
print(f"Specific enthalpy: H_L = {H_L:.6f}, H_R = {H_R:.6f}")

# Lab-frame baryon number density: N = γn
N_L = gamma_L * n_L
N_R = gamma_R * n_R

print(f"Lab-frame density: N_L = {N_L:.6f}, N_R = {N_R:.6f}")

# Particle volume: V = ν/N
V_L = nu / N_L
V_R = nu / N_R

print(f"Particle volume: V_L = {V_L:.6e}, V_R = {V_R:.6e}")
print(f"Volume squared: V²_L = {V_L**2:.6e}, V²_R = {V_R**2:.6e}")

# Sound speed: c_s = √[(γ-1)(H-1)/H]
cs_L = np.sqrt((gamma_eos - 1.0) * (H_L - 1.0) / H_L)
cs_R = np.sqrt((gamma_eos - 1.0) * (H_R - 1.0) / H_R)

print(f"Sound speed: c_s,L = {cs_L:.6f}, c_s,R = {cs_R:.6f}")

# ============================================================================
# Particle Spacing and Smoothing Length
# ============================================================================
print("\n--- Particle Spacing and Smoothing Length ---")

# Domain: x ∈ [-0.5, 0.5], discontinuity at x = 0
domain_length = 1.0
x_left = -0.5
x_right = 0.5

# Particle spacing
dx_L = nu / N_L  # spacing = ν/N for uniform 1D
dx_R = nu / N_R

print(f"Left region particle spacing: dx_L = {dx_L:.6e}")
print(f"Right region particle spacing: dx_R = {dx_R:.6e}")

# Smoothing length for neighbor_num neighbors in 1D
# h ≈ (neighbor_num/2) × dx (since kernel support is 2h)
h_L = (neighbor_num / 2.0) * dx_L
h_R = (neighbor_num / 2.0) * dx_R

print(f"Expected smoothing length: h_L ≈ {h_L:.6e}, h_R ≈ {h_R:.6e}")

# From simulation output (observed values)
h_L_obs = 0.00398  # from diagnostic
h_R_obs = 0.0316   # from diagnostic

print(f"Observed smoothing length: h_L = {h_L_obs:.6e}, h_R = {h_R_obs:.6e}")
print(f"Ratio (obs/expected): h_L = {h_L_obs/h_L:.3f}, h_R = {h_R_obs/h_R:.3f}")

# ============================================================================
# Kernel Gradient Magnitude
# ============================================================================
print("\n--- Kernel Gradient Estimate ---")

# For cubic spline kernel in 1D with support 2h:
# W(q,h) = σ/h × f(q) where q = r/h, σ = 2/3 for 1D
# ∇W = (dW/dq) × (dq/dr) × (r̂) = (σ/h²) × f'(q) × (r̂)
# At small q (near neighbor): |f'(q)| ~ q for cubic spline
# Maximum |f'(q)| ≈ 1 at q ≈ 1

sigma_1d = 2.0/3.0  # normalization constant for 1D cubic spline

# Typical neighbor at r ≈ dx
r_typical_L = dx_L
r_typical_R = dx_R

# q = r / (2h) for our kernel call
q_L = r_typical_L / (2.0 * h_L_obs)
q_R = r_typical_R / (2.0 * h_R_obs)

print(f"Typical neighbor separation: r_L ≈ {r_typical_L:.6e}, r_R ≈ {r_typical_R:.6e}")
print(f"Dimensionless distance: q_L ≈ {q_L:.4f}, q_R ≈ {q_R:.4f}")

# Gradient magnitude estimate: |∇W| ~ σ/(2h)² for our kernel call with 2h
dW_mag_L = sigma_1d / (2.0 * h_L_obs)**2
dW_mag_R = sigma_1d / (2.0 * h_R_obs)**2

print(f"Estimated |∇W|: left ≈ {dW_mag_L:.2e}, right ≈ {dW_mag_R:.2e}")

# From simulation (observed)
dW_obs_L = 4820  # from diagnostic
dW_obs_R = 46305  # from diagnostic

print(f"Observed |∇W|: left = {dW_obs_L:.2e}, right = {dW_obs_R:.2e}")
print(f"Ratio (obs/estimate): left = {dW_obs_L/dW_mag_L:.1f}, right = {dW_obs_R/dW_mag_R:.1f}")

# ============================================================================
# Riemann Solver Estimate
# ============================================================================
print("\n--- Riemann Solver (Contact Discontinuity) ---")

# For contact in Sod problem, P* is between P_L and P_R
# Exact solution gives P* ≈ 0.3 for standard Sod
# For our P_L=1.0, P_R=0.1, expect similar

P_star_estimate = 0.3  # rough estimate
v_star_estimate = 0.2  # contact velocity

print(f"Expected P* ≈ {P_star_estimate} (between {P_R} and {P_L})")
print(f"Expected v* ≈ {v_star_estimate}")

# From simulation (observed)
P_star_obs = 0.668  # from diagnostic
print(f"Observed P* = {P_star_obs:.3f}")

# ============================================================================
# Force Calculation (Equation 64)
# ============================================================================
print("\n--- Force Calculation (Eq. 64) ---")
print("⟨ν_i Ṡ_i⟩ = -Σ_j P*_{ij} V²_{ij} [∇_i W - ∇_j W]")

# Consider left-side particle (i) interacting with another left-side particle (j)
print("\nCase 1: Left particle i with left neighbor j")
V_i = V_L
V_j = V_L
V2_ij = 0.5 * (V_i**2 + V_j**2)

print(f"  V_i = V_j = {V_L:.6e}")
print(f"  V²_ij = 0.5(V_i² + V_j²) = {V2_ij:.6e}")

# For similar h, ∇W_ij ≈ 0, but let's assume 10% difference
h_i = h_L_obs
h_j = h_L_obs * 1.1  # 10% larger
dW_i = sigma_1d / (2.0 * h_i)**2
dW_j = sigma_1d / (2.0 * h_j)**2
dW_ij = abs(dW_i - dW_j)

print(f"  Assuming h_j = 1.1×h_i: |∇W_ij| ≈ {dW_ij:.2e}")

# Force contribution from one neighbor
dS_contrib = P_star_obs * V2_ij * dW_ij

print(f"  dS_contrib = P* × V²_ij × |∇W_ij| = {dS_contrib:.6e}")

# Total force from ~50 neighbors
dS_total = neighbor_num * dS_contrib

print(f"  Total |dS| ≈ {neighbor_num} × {dS_contrib:.6e} = {dS_total:.6e}")

# Divide by ν to get Ṡ
dS_over_nu = dS_total / nu

print(f"  |dS|/ν = {dS_over_nu:.2e}")

print("\nCase 2: Using observed kernel gradient from diagnostic")
# Use actual observed values
dW_ij_obs = 1523  # from diagnostic where h_i ≠ h_j
V2_ij_obs = 7.29e-08  # from V² diagnostic

dS_contrib_obs = P_star_obs * V2_ij_obs * dW_ij_obs
dS_total_obs = neighbor_num * dS_contrib_obs
dS_over_nu_obs = dS_total_obs / nu

print(f"  Observed: P* = {P_star_obs}, V²_ij = {V2_ij_obs:.2e}, |∇W_ij| = {dW_ij_obs}")
print(f"  dS_contrib = {dS_contrib_obs:.6e}")
print(f"  Total |dS| ≈ {dS_total_obs:.6e}")
print(f"  |dS|/ν = {dS_over_nu_obs:.2e}")

# ============================================================================
# Comparison with Simulation
# ============================================================================
print("\n" + "=" * 80)
print("COMPARISON WITH SIMULATION")
print("=" * 80)

# From simulation diagnostics
dS_over_nu_sim = 807238  # from Predict call 1

print(f"\nSimulation reports: |dS|/ν = {dS_over_nu_sim:.2e}")
print(f"Expected (Case 1):  |dS|/ν = {dS_over_nu_obs:.2e}")
print(f"Discrepancy factor: {dS_over_nu_sim / dS_over_nu_obs:.1f}×")

# ============================================================================
# Root Cause Investigation
# ============================================================================
print("\n" + "=" * 80)
print("ROOT CAUSE INVESTIGATION")
print("=" * 80)

print("\n1. Check V²_{ij} formula:")
print(f"   Current: V²_ij = 0.5(V_i² + V_j²) = {V2_ij:.6e}")
print(f"   Units: [V] = [ν]/[N] = [baryon]/[baryon/L] = [L]")
print(f"   So [V²] = [L²]")

print("\n2. Check force equation units:")
print("   LHS: ⟨ν Ṡ⟩ has units [baryon][momentum/time] = [Force]")
print("   RHS: P × V² × ∇W")
print(f"   [P] = [Force/L²]")
print(f"   [∇W] = [1/L^(d+1)] = [1/L²] for d=1")
print(f"   [V²] should give [L⁴] for dimensional consistency")
print("   BUT: (ν/N)² gives [L²] not [L⁴]!")

print("\n3. Possible issue: V²_{ij} definition")
print("   Paper Eq. 29-30 defines V²_{ij}(h) as an INTEGRAL:")
print("   V²_{ij}(h) ≡ ∫ [1/N²(x)] W(x-x_i, h) dx")
print("   This is NOT the same as (V_i)² = (ν_i/N_i)²!")

print("\n4. Alternative interpretation:")
print("   If V²_{ij} means the square of volume FACTOR, not volume squared:")
print(f"   Try: V²_ij = ν² × 0.5(1/N_i² + 1/N_j²)")
alt_V2_ij = nu**2 * 0.5 * (1.0/N_L**2 + 1.0/N_L**2)
print(f"   V²_ij = {nu}² × 0.5(1/{N_L}² + 1/{N_L}²) = {alt_V2_ij:.6e}")

dS_contrib_alt = P_star_obs * alt_V2_ij * dW_ij_obs
dS_total_alt = neighbor_num * dS_contrib_alt
dS_over_nu_alt = dS_total_alt / nu

print(f"   Then: dS_contrib = {dS_contrib_alt:.6e}")
print(f"   Total |dS|/ν = {dS_over_nu_alt:.2e}")
print(f"   Discrepancy: {dS_over_nu_sim / dS_over_nu_alt:.1f}×")

print("\n5. Check if ν² factor is already included:")
print(f"   Current V²_ij already has ν² in numerator")
print(f"   V_i = ν/N = {nu}/{N_L} = {V_L:.6e}")
print(f"   V_i² = ν²/N² = {nu**2}/{N_L**2} = {V_L**2:.6e}")
print(f"   This ALREADY includes ν²!")

print("\n6. Dimensional analysis of actual equation:")
print("   If we use V²_ij = (ν/N)²:")
print(f"   Force ~ P × (ν²/N²) × (1/h²)")
print(f"   Force ~ [F/L²] × [ν²·L²/ν²] × [1/L²] = [F/L²]")
print("   This gives [F/L²] not [F] - MISSING ν FACTOR?")

print("\n7. Compare to non-relativistic GSPH:")
print("   Non-rel: ma = Σ m_j P* (1/ρ_i² + 1/ρ_j²) ∇W")
print("   With m = ρV: a = Σ (m_j/ρ_i²) P* ∇W")
print(f"   SR analog: Ṡ should be ~ Σ (ν_j P* / N_i²) ∇W")
print("   But our equation has V²_ij = (ν/N)², giving extra ν factor!")

print("\n" + "=" * 80)
print("CONCLUSION")
print("=" * 80)
print("\n*** SUSPECTED ROOT CAUSE ***")
print("The V²_{ij} term in Eq. 64 likely does NOT mean (V_i)² = (ν_i/N_i)²")
print("Instead, it may be a DENSITY-WEIGHTED INTEGRAL as in Eq. 29-30,")
print("or the notation V²_{ij} has different units than [volume²].")
print("\nThe force is systematically too large by factor ~1000×")
print("This suggests V²_{ij} should be ~1000× smaller than current calculation.")
print("\nRECOMMENDATION: Check paper Eq. 29-30 carefully for V²_{ij} definition.")
print("=" * 80)
