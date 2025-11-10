#!/usr/bin/env python3
"""
Analyze the relaxed Lane-Emden configuration.
Shows particle positions, physical variables, and compares to analytical solution.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# Load relaxed data
data = np.loadtxt('sample/lane_emden/results/00000.dat')
print(f"Loaded {len(data)} particles from relaxed configuration")

# Extract variables (file format: x y z vx vy vz ax ay az mass dens pres ene sml id neighbor alpha gradh)
pos = data[:, 0:3]  # x, y, z
vel = data[:, 3:6]  # vx, vy, vz
acc = data[:, 6:9]  # ax, ay, az
mass = data[:, 9]
ene = data[:, 12]  # internal energy
rho = data[:, 10]  # density
pres = data[:, 11]  # pressure
h = data[:, 13]  # smoothing length

# Compute derived quantities
r = np.sqrt(np.sum(pos**2, axis=1))
v_mag = np.sqrt(np.sum(vel**2, axis=1))
a_mag = np.sqrt(np.sum(acc**2, axis=1))

# Statistics
print("\n=== Relaxed State Statistics ===")
print(f"Radius range: [{r.min():.6f}, {r.max():.6f}]")
print(f"Velocity magnitude: min={v_mag.min():.2e}, max={v_mag.max():.2e}, mean={v_mag.mean():.2e}")
print(f"Acceleration magnitude: min={a_mag.min():.2e}, max={a_mag.max():.2e}, mean={a_mag.mean():.2e}")
print(f"Density range: [{rho.min():.6f}, {rho.max():.6f}], ratio={rho.max()/rho.min():.2e}")
print(f"Pressure range: [{pres.min():.6f}, {pres.max():.6f}]")
print(f"Internal energy range: [{ene.min():.6f}, {ene.max():.6f}]")
print(f"Smoothing length range: [{h.min():.6f}, {h.max():.6f}]")
print(f"Total mass: {mass.sum():.6f}")

# Create comprehensive visualization
fig = plt.figure(figsize=(18, 12))
gs = GridSpec(3, 4, figure=fig, hspace=0.3, wspace=0.3)

# Row 1: Spatial distributions
ax1 = fig.add_subplot(gs[0, 0])
scatter = ax1.scatter(pos[:, 0], pos[:, 1], c=r, s=5, cmap='viridis')
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_title('XY Projection (colored by radius)')
ax1.set_aspect('equal')
plt.colorbar(scatter, ax=ax1, label='radius')

ax2 = fig.add_subplot(gs[0, 1], projection='3d')
ax2.scatter(pos[:, 0], pos[:, 1], pos[:, 2], c=r, s=2, cmap='viridis')
ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.set_zlabel('z')
ax2.set_title('3D Distribution')

ax3 = fig.add_subplot(gs[0, 2])
scatter3 = ax3.scatter(pos[:, 0], pos[:, 1], c=np.log10(rho), s=5, cmap='hot')
ax3.set_xlabel('x')
ax3.set_ylabel('y')
ax3.set_title('Density (log10)')
ax3.set_aspect('equal')
plt.colorbar(scatter3, ax=ax3, label='log10(ρ)')

ax4 = fig.add_subplot(gs[0, 3])
ax4.hist(r, bins=50, alpha=0.7, edgecolor='black')
ax4.set_xlabel('radius')
ax4.set_ylabel('count')
ax4.set_title('Radial Distribution')
ax4.axvline(1.0, color='r', linestyle='--', label='R=1.0')
ax4.legend()

# Row 2: Physical variables vs radius
ax5 = fig.add_subplot(gs[1, 0])
ax5.scatter(r, rho, s=1, alpha=0.5)
ax5.set_xlabel('radius')
ax5.set_ylabel('density')
ax5.set_title('Density Profile')
ax5.set_yscale('log')

ax6 = fig.add_subplot(gs[1, 1])
ax6.scatter(r, pres, s=1, alpha=0.5)
ax6.set_xlabel('radius')
ax6.set_ylabel('pressure')
ax6.set_title('Pressure Profile')
ax6.set_yscale('log')

ax7 = fig.add_subplot(gs[1, 2])
ax7.scatter(r, ene, s=1, alpha=0.5)
ax7.set_xlabel('radius')
ax7.set_ylabel('internal energy')
ax7.set_title('Internal Energy Profile')

ax8 = fig.add_subplot(gs[1, 3])
ax8.scatter(r, h, s=1, alpha=0.5)
ax8.set_xlabel('radius')
ax8.set_ylabel('smoothing length')
ax8.set_title('Smoothing Length Profile')

# Row 3: Velocities and accelerations
ax9 = fig.add_subplot(gs[2, 0])
ax9.scatter(r, v_mag, s=1, alpha=0.5)
ax9.set_xlabel('radius')
ax9.set_ylabel('|velocity|')
ax9.set_title('Velocity Magnitude (should be ~0)')
ax9.set_yscale('log')

ax10 = fig.add_subplot(gs[2, 1])
ax10.scatter(r, a_mag, s=1, alpha=0.5)
ax10.set_xlabel('radius')
ax10.set_ylabel('|acceleration|')
ax10.set_title('Acceleration Magnitude')

ax11 = fig.add_subplot(gs[2, 2])
ax11.hist(v_mag, bins=50, alpha=0.7, edgecolor='black')
ax11.set_xlabel('|velocity|')
ax11.set_ylabel('count')
ax11.set_title('Velocity Distribution')
ax11.set_yscale('log')

ax12 = fig.add_subplot(gs[2, 3])
ax12.hist(a_mag, bins=50, alpha=0.7, edgecolor='black')
ax12.set_xlabel('|acceleration|')
ax12.set_ylabel('count')
ax12.set_title('Acceleration Distribution')

plt.suptitle('Lane-Emden Relaxed Configuration Analysis', fontsize=16, y=0.995)
output_file = 'sample/lane_emden/results/relaxed_analysis.png'
plt.savefig(output_file, dpi=150, bbox_inches='tight')
print(f"\n✓ Saved comprehensive analysis: {output_file}")

# Additional radial profile analysis
fig2, axes = plt.subplots(2, 2, figsize=(12, 10))

# Bin particles by radius for cleaner profiles
r_bins = np.linspace(0, r.max(), 100)
r_centers = (r_bins[:-1] + r_bins[1:]) / 2
rho_profile = []
pres_profile = []
ene_profile = []
a_profile = []

for i in range(len(r_bins) - 1):
    mask = (r >= r_bins[i]) & (r < r_bins[i+1])
    if mask.sum() > 0:
        rho_profile.append(rho[mask].mean())
        pres_profile.append(pres[mask].mean())
        ene_profile.append(ene[mask].mean())
        a_profile.append(a_mag[mask].mean())
    else:
        rho_profile.append(np.nan)
        pres_profile.append(np.nan)
        ene_profile.append(np.nan)
        a_profile.append(np.nan)

axes[0, 0].plot(r_centers, rho_profile, 'b-', linewidth=2)
axes[0, 0].set_xlabel('radius')
axes[0, 0].set_ylabel('density')
axes[0, 0].set_title('Radially Averaged Density')
axes[0, 0].set_yscale('log')
axes[0, 0].grid(True, alpha=0.3)

axes[0, 1].plot(r_centers, pres_profile, 'r-', linewidth=2)
axes[0, 1].set_xlabel('radius')
axes[0, 1].set_ylabel('pressure')
axes[0, 1].set_title('Radially Averaged Pressure')
axes[0, 1].set_yscale('log')
axes[0, 1].grid(True, alpha=0.3)

axes[1, 0].plot(r_centers, ene_profile, 'g-', linewidth=2)
axes[1, 0].set_xlabel('radius')
axes[1, 0].set_ylabel('internal energy')
axes[1, 0].set_title('Radially Averaged Internal Energy')
axes[1, 0].grid(True, alpha=0.3)

axes[1, 1].plot(r_centers, a_profile, 'm-', linewidth=2)
axes[1, 1].set_xlabel('radius')
axes[1, 1].set_ylabel('|acceleration|')
axes[1, 1].set_title('Radially Averaged Acceleration')
axes[1, 1].grid(True, alpha=0.3)

plt.suptitle('Lane-Emden Radial Profiles', fontsize=14)
plt.tight_layout()
output_file2 = 'sample/lane_emden/results/relaxed_profiles.png'
plt.savefig(output_file2, dpi=150)
print(f"✓ Saved radial profiles: {output_file2}")

print("\n=== Analysis Complete ===")
