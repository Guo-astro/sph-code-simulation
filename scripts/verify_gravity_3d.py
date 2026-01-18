#!/usr/bin/env python3
"""
Verify 3D gravity calculation against analytic solution for a uniform sphere.

This script:
1. Creates a uniform density sphere particle distribution
2. Computes direct N-body gravity with softening (using same kernels as SPH code)
3. Compares with analytic solution
4. Plots both for visual verification

Analytic solution for uniform sphere:
- Inside (r < R):  g = (4/3) * π * G * ρ * r  (linear with radius)
- Outside (r ≥ R): g = G * M / r²            (point mass)
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# =============================================================================
# Physical Constants and Parameters
# =============================================================================
G = 1.0           # Gravitational constant (code units)
M_total = 1.0     # Total mass
R_sphere = 1.0    # Sphere radius
N_particles = 2000  # Number of particles
h_softening = 0.1   # Softening length (fixed)

# Derived quantities
rho = M_total / (4/3 * np.pi * R_sphere**3)  # Uniform density
m_particle = M_total / N_particles            # Mass per particle

# =============================================================================
# Hernquist-Katz Softening Kernels (same as in C++ code)
# =============================================================================
def hernquist_katz_f(r, h):
    """Potential softening kernel f(r,h) - Hernquist & Katz 1989"""
    e = h * 0.5
    u = r / e

    if u < 1.0:
        return (-0.5 * u * u * (1.0/3.0 - 3.0/20.0 * u*u + u*u*u/20.0) + 1.4) / e
    elif u < 2.0:
        return -1.0 / (15.0 * r) + (-u*u * (4.0/3.0 - u + 0.3*u*u - u*u*u/30.0) + 1.6) / e
    else:
        return 1.0 / r

def hernquist_katz_g(r, h):
    """Force softening kernel g(r,h) - Hernquist & Katz 1989"""
    e = h * 0.5
    u = r / e

    if u < 1.0:
        return (4.0/3.0 - 1.2*u*u + 0.5*u*u*u) / (e**3)
    elif u < 2.0:
        return (-1.0/15.0 + 8.0/3.0*u**3 - 3.0*u**4 + 1.2*u**5 - u**6/6.0) / (r**3)
    else:
        return 1.0 / (r**3)

# Vectorized versions
hernquist_katz_f_vec = np.vectorize(hernquist_katz_f)
hernquist_katz_g_vec = np.vectorize(hernquist_katz_g)

# =============================================================================
# Particle Distribution
# =============================================================================
def create_uniform_sphere(n_particles, radius, seed=42):
    """Create uniformly distributed particles in a sphere using rejection sampling"""
    rng = np.random.default_rng(seed)

    positions = []
    while len(positions) < n_particles:
        # Generate random points in cube
        batch_size = n_particles * 2
        x = rng.uniform(-radius, radius, batch_size)
        y = rng.uniform(-radius, radius, batch_size)
        z = rng.uniform(-radius, radius, batch_size)

        # Keep only points inside sphere
        r2 = x**2 + y**2 + z**2
        mask = r2 <= radius**2

        for i in range(len(x)):
            if mask[i] and len(positions) < n_particles:
                positions.append([x[i], y[i], z[i]])

    return np.array(positions)

# =============================================================================
# Gravity Computation
# =============================================================================
def compute_gravity_direct(positions, masses, h):
    """
    Compute gravitational acceleration using direct N-body summation
    with Hernquist-Katz softening (same as SPH code with lookup tables)
    """
    n = len(positions)
    accelerations = np.zeros((n, 3))

    for i in range(n):
        for j in range(n):
            if i == j:
                continue

            r_ij = positions[i] - positions[j]
            r = np.linalg.norm(r_ij)

            if r < 1e-10:
                continue

            # Use Hernquist-Katz force kernel
            g_val = hernquist_katz_g(r, h)

            # F_i = -G * m_i * m_j * g(r,h) * r_ij
            # a_i = F_i / m_i = -G * m_j * g(r,h) * r_ij
            accelerations[i] -= G * masses[j] * g_val * r_ij

    return accelerations

def compute_gravity_analytic(positions, M_total, R_sphere):
    """
    Compute analytic gravitational acceleration for uniform sphere

    Inside (r < R):  g = (4/3) * π * G * ρ * r = G * M * r / R³
    Outside (r ≥ R): g = G * M / r²
    """
    n = len(positions)
    accelerations = np.zeros((n, 3))

    for i in range(n):
        r_vec = positions[i]
        r = np.linalg.norm(r_vec)

        if r < 1e-10:
            continue

        r_hat = r_vec / r

        if r < R_sphere:
            # Inside: linear with r
            g_mag = G * M_total * r / (R_sphere**3)
        else:
            # Outside: inverse square
            g_mag = G * M_total / (r**2)

        # Acceleration points toward center (negative r_hat)
        accelerations[i] = -g_mag * r_hat

    return accelerations

# =============================================================================
# Main Computation
# =============================================================================
print("=" * 60)
print("3D Gravity Verification: Uniform Sphere")
print("=" * 60)
print(f"Parameters:")
print(f"  N_particles = {N_particles}")
print(f"  M_total     = {M_total}")
print(f"  R_sphere    = {R_sphere}")
print(f"  h_softening = {h_softening}")
print(f"  ρ (density) = {rho:.4f}")
print()

# Create particle distribution
print("Creating uniform sphere distribution...")
positions = create_uniform_sphere(N_particles, R_sphere)
masses = np.full(N_particles, m_particle)

# Compute radii
radii = np.linalg.norm(positions, axis=1)

# Compute gravity using direct summation with softening
print("Computing gravity (direct N-body with Hernquist-Katz softening)...")
acc_numerical = compute_gravity_direct(positions, masses, h_softening)

# Compute analytic gravity
print("Computing analytic gravity...")
acc_analytic = compute_gravity_analytic(positions, M_total, R_sphere)

# Compute magnitudes (acceleration points toward center, so magnitude is |a|)
g_numerical = np.linalg.norm(acc_numerical, axis=1)
g_analytic = np.linalg.norm(acc_analytic, axis=1)

# =============================================================================
# Error Analysis
# =============================================================================
# Compute relative error (excluding very small accelerations to avoid division issues)
mask = g_analytic > 1e-10
relative_error = np.zeros_like(g_numerical)
relative_error[mask] = np.abs(g_numerical[mask] - g_analytic[mask]) / g_analytic[mask]

print()
print("Error Analysis:")
print(f"  Mean relative error:   {np.mean(relative_error[mask]):.4f}")
print(f"  Max relative error:    {np.max(relative_error[mask]):.4f}")
print(f"  Median relative error: {np.median(relative_error[mask]):.4f}")

# Separate interior vs exterior
interior_mask = radii < R_sphere * 0.9
exterior_mask = radii > R_sphere * 0.1  # Actually all are inside for this test
print(f"  Interior error (r < 0.9R): {np.mean(relative_error[interior_mask & mask]):.4f}")

# =============================================================================
# Plotting
# =============================================================================
print()
print("Generating plots...")

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Plot 1: Gravity magnitude vs radius
ax1 = axes[0, 0]
# Sort by radius for cleaner lines
sort_idx = np.argsort(radii)
r_sorted = radii[sort_idx]
g_num_sorted = g_numerical[sort_idx]
g_ana_sorted = g_analytic[sort_idx]

ax1.scatter(radii, g_numerical, s=5, alpha=0.5, label='Numerical (H-K softening)', c='blue')
ax1.scatter(radii, g_analytic, s=5, alpha=0.5, label='Analytic (uniform sphere)', c='red')

# Add analytic curve
r_theory = np.linspace(0.01, R_sphere, 100)
g_theory = G * M_total * r_theory / (R_sphere**3)  # Inside uniform sphere
ax1.plot(r_theory, g_theory, 'k--', linewidth=2, label='Theory: g = GMr/R³')

ax1.set_xlabel('Radius r', fontsize=12)
ax1.set_ylabel('|g| (acceleration magnitude)', fontsize=12)
ax1.set_title('Gravitational Acceleration vs Radius', fontsize=14)
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3)

# Plot 2: Relative error vs radius
ax2 = axes[0, 1]
ax2.scatter(radii[mask], relative_error[mask] * 100, s=5, alpha=0.5, c='green')
ax2.axhline(y=10, color='r', linestyle='--', label='10% error threshold')
ax2.axhline(y=20, color='orange', linestyle='--', label='20% error threshold')
ax2.set_xlabel('Radius r', fontsize=12)
ax2.set_ylabel('Relative Error (%)', fontsize=12)
ax2.set_title('Relative Error vs Radius', fontsize=14)
ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.3)
ax2.set_ylim(0, min(50, np.max(relative_error[mask]) * 100 * 1.1))

# Plot 3: Numerical vs Analytic scatter
ax3 = axes[1, 0]
ax3.scatter(g_analytic, g_numerical, s=5, alpha=0.5, c='purple')
max_g = max(np.max(g_analytic), np.max(g_numerical))
ax3.plot([0, max_g], [0, max_g], 'k--', linewidth=2, label='Perfect agreement')
ax3.set_xlabel('Analytic |g|', fontsize=12)
ax3.set_ylabel('Numerical |g|', fontsize=12)
ax3.set_title('Numerical vs Analytic Comparison', fontsize=14)
ax3.legend(fontsize=10)
ax3.grid(True, alpha=0.3)
ax3.set_aspect('equal')

# Plot 4: 3D particle distribution colored by error
ax4 = axes[1, 1]
ax4 = fig.add_subplot(2, 2, 4, projection='3d')
# Subsample for 3D plot
subsample = 500
idx = np.random.choice(N_particles, min(subsample, N_particles), replace=False)
scatter = ax4.scatter(positions[idx, 0], positions[idx, 1], positions[idx, 2],
                      c=relative_error[idx] * 100, cmap='RdYlGn_r', s=10, alpha=0.7)
ax4.set_xlabel('X')
ax4.set_ylabel('Y')
ax4.set_zlabel('Z')
ax4.set_title('Particle Distribution (color = error %)', fontsize=14)
plt.colorbar(scatter, ax=ax4, label='Relative Error (%)', shrink=0.6)

plt.tight_layout()
plt.savefig('/Users/yansongguo/personal/research/shamrock-wrapper/sph-code-simulation/scripts/gravity_3d_verification.png',
            dpi=150, bbox_inches='tight')
print(f"Saved: gravity_3d_verification.png")
plt.show()

# =============================================================================
# Summary
# =============================================================================
print()
print("=" * 60)
print("VERIFICATION SUMMARY")
print("=" * 60)
mean_err = np.mean(relative_error[mask]) * 100
if mean_err < 20:
    print(f"✓ PASSED: Mean error {mean_err:.1f}% < 20% threshold")
    print("  The softened gravity matches analytic solution within expected")
    print("  particle discretization error.")
else:
    print(f"✗ WARNING: Mean error {mean_err:.1f}% >= 20% threshold")
    print("  Error higher than expected. Check softening parameters.")

print()
print("Note: Error is expected due to:")
print("  1. Particle discretization (finite N)")
print("  2. Gravitational softening (smooths forces at small r)")
print("  3. Random particle placement (Poisson noise)")
print()
print("The lookup table integration is VERIFIED if numerical results")
print("follow the same trend as the analytic solution.")
