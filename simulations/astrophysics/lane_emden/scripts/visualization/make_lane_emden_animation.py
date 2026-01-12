#!/usr/bin/env python3
"""
Generate animation and final snapshot for Lane-Emden simulation results.
Usage: python3 scripts/make_lane_emden_animation.py [data_dir] [output_prefix]
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import sys
import os
from pathlib import Path
from matplotlib.animation import FuncAnimation, PillowWriter

# Try to import tqdm for progress bar
try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False
    print("Note: Install 'tqdm' for progress bars: pip install tqdm")

# Get data directory and output prefix from command line or use defaults
data_dir = sys.argv[1] if len(sys.argv) > 1 else "lane_emden/results/polytrope_n1.5_3d"
output_prefix = sys.argv[2] if len(sys.argv) > 2 else "lane_emden"

print('Loading data files...')
# Try CSV first (new format), then fall back to .dat (old format)
files = sorted(glob.glob(f'{data_dir}/snapshot_*.csv'))
if len(files) == 0:
    files = sorted(glob.glob(f'{data_dir}/[0-9]*.dat'))  # Old format fallback

if len(files) == 0:
    print('ERROR: No data files found!')
    print(f'Searched in: {data_dir}')
    print('Looking for: snapshot_*.csv or [0-9]*.dat')
    sys.exit(1)

file_format = 'csv' if files[0].endswith('.csv') else 'dat'
print(f'Found {len(files)} {file_format.upper()} files')

# Create animation with 5 subplots
fig = plt.figure(figsize=(20, 8))
ax1 = fig.add_subplot(2, 3, 1)  # XY projection
ax2 = fig.add_subplot(2, 3, 2)  # Density profile
ax3 = fig.add_subplot(2, 3, 3)  # Force profile
ax4 = fig.add_subplot(2, 3, 4, projection='3d')  # 3D view
ax5 = fig.add_subplot(2, 3, 5)  # Radial distribution
ax6 = fig.add_subplot(2, 3, 6)  # Acceleration magnitude profile

# Progress bar for animation frames
pbar = None
colorbars_created = False

def update(frame):
    global pbar, colorbars_created
    if pbar is None and HAS_TQDM:
        pbar = tqdm(total=len(files), desc='Generating frames', unit='frame')
    
    # Only clear axes, don't remove colorbars
    for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
        ax.clear()
    
    # Load data based on file format
    if file_format == 'csv':
        # CSV format: skip 51 header lines (metadata + column names)
        data = np.loadtxt(files[frame], delimiter=',', skiprows=51)
        # CSV column order: id,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,acc_x,acc_y,acc_z,mass,dens,pres,ene,sml,sound,alpha,balsara,gradh,phi,neighbor
        pos = data[:, 1:4]      # pos_x, pos_y, pos_z (columns 1, 2, 3)
        acc = data[:, 7:10]     # acc_x, acc_y, acc_z (columns 7, 8, 9)
        dens = data[:, 11]      # dens (column 11)
    else:
        # Old .dat format: space-delimited, no header
        data = np.loadtxt(files[frame])
        # File format: x y z vx vy vz ax ay az mass dens pres ene sml id neighbor alpha gradh
        pos = data[:, 0:3]      # x, y, z (columns 0, 1, 2)
        acc = data[:, 6:9]      # ax, ay, az (columns 6, 7, 8)
        dens = data[:, 10]      # density (column 10)
    
    r = np.sqrt(np.sum(pos**2, axis=1))
    acc_mag = np.sqrt(np.sum(acc**2, axis=1))
    
    # 1. XY projection colored by density
    ax1.scatter(pos[:, 0], pos[:, 1], c=dens, s=2, cmap='hot', vmin=0, vmax=2.0)
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_title(f'XY Projection (Frame {frame})')
    ax1.set_xlim(-1.5, 1.5)
    ax1.set_ylim(-1.5, 1.5)
    ax1.set_aspect('equal')
    
    # 2. Density profile (radial bins)
    r_bins = np.linspace(0, 1.5, 30)
    r_centers = 0.5 * (r_bins[1:] + r_bins[:-1])
    dens_profile = []
    for i in range(len(r_bins)-1):
        mask = (r >= r_bins[i]) & (r < r_bins[i+1])
        if np.sum(mask) > 0:
            dens_profile.append(np.mean(dens[mask]))
        else:
            dens_profile.append(0)
    
    ax2.plot(r_centers, dens_profile, 'b-', linewidth=2)
    ax2.scatter(r, dens, c='lightblue', s=1, alpha=0.3)
    ax2.set_xlabel('Radius')
    ax2.set_ylabel('Density')
    ax2.set_title('Density Profile')
    ax2.set_xlim(0, 1.5)
    ax2.set_ylim(0, 2.0)
    ax2.grid(True, alpha=0.3)
    
    # 3. Force/Acceleration profile
    acc_profile = []
    for i in range(len(r_bins)-1):
        mask = (r >= r_bins[i]) & (r < r_bins[i+1])
        if np.sum(mask) > 0:
            acc_profile.append(np.mean(acc_mag[mask]))
        else:
            acc_profile.append(0)
    
    ax3.plot(r_centers, acc_profile, 'r-', linewidth=2)
    ax3.scatter(r, acc_mag, c='pink', s=1, alpha=0.3)
    ax3.set_xlabel('Radius')
    ax3.set_ylabel('|Acceleration|')
    ax3.set_title('Force/Acceleration Profile')
    ax3.set_xlim(0, 1.5)
    ax3.grid(True, alpha=0.3)
    
    # 4. 3D view colored by density
    ax4.scatter(pos[:, 0], pos[:, 1], pos[:, 2], c=dens, s=2, cmap='hot', vmin=0, vmax=2.0)
    ax4.set_xlabel('x')
    ax4.set_ylabel('y')
    ax4.set_zlabel('z')
    ax4.set_title('3D View')
    ax4.set_xlim(-1.5, 1.5)
    ax4.set_ylim(-1.5, 1.5)
    ax4.set_zlim(-1.5, 1.5)
    
    # 5. Radial distribution histogram
    ax5.hist(r, bins=50, range=(0, 1.5), alpha=0.7, color='green')
    ax5.set_xlabel('Radius')
    ax5.set_ylabel('Particle Count')
    ax5.set_title('Radial Distribution')
    ax5.grid(True, alpha=0.3)
    
    # 6. Acceleration magnitude distribution
    ax6.hist(acc_mag, bins=50, alpha=0.7, color='orange')
    ax6.set_xlabel('|Acceleration|')
    ax6.set_ylabel('Count')
    ax6.set_title('Acceleration Distribution')
    ax6.set_yscale('log')
    ax6.grid(True, alpha=0.3)
    
    if pbar and HAS_TQDM:
        pbar.update(1)
    elif not HAS_TQDM and frame % 10 == 0:
        print(f'  Processing frame {frame}/{len(files)}...')
    
    return []

print('Creating animation...')
anim = FuncAnimation(fig, update, frames=len(files), interval=100, blit=False)

# Create output directory for animations
animation_dir = f'{data_dir}/animations'
Path(animation_dir).mkdir(parents=True, exist_ok=True)

gif_file = f'{animation_dir}/{output_prefix}.gif'
print('\nSaving animation (this may take a while)...')
anim.save(gif_file, writer=PillowWriter(fps=10), progress_callback=lambda i, n: None)
if pbar and HAS_TQDM:
    pbar.close()
print(f'✓ Saved: {gif_file}')
plt.close()

# Create final snapshot with all profiles
print('Creating final snapshot...')
fig2 = plt.figure(figsize=(20, 8))
ax1 = fig2.add_subplot(2, 3, 1)
ax2 = fig2.add_subplot(2, 3, 2)
ax3 = fig2.add_subplot(2, 3, 3)
ax4 = fig2.add_subplot(2, 3, 4, projection='3d')
ax5 = fig2.add_subplot(2, 3, 5)
ax6 = fig2.add_subplot(2, 3, 6)

# Load final snapshot with correct format
if file_format == 'csv':
    data = np.loadtxt(files[-1], delimiter=',', skiprows=51)
    # CSV column order: id,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,acc_x,acc_y,acc_z,mass,dens,pres,ene,sml,sound,alpha,balsara,gradh,phi,neighbor
    pos = data[:, 1:4]      # pos_x, pos_y, pos_z (columns 1, 2, 3)
    acc = data[:, 7:10]     # acc_x, acc_y, acc_z (columns 7, 8, 9)
    dens = data[:, 11]      # dens (column 11)
else:
    data = np.loadtxt(files[-1])
    # File format: x y z vx vy vz ax ay az mass dens pres ene sml id neighbor alpha gradh
    pos = data[:, 0:3]      # x, y, z (columns 0, 1, 2)
    acc = data[:, 6:9]      # ax, ay, az (columns 6, 7, 8)
    dens = data[:, 10]      # density (column 10)

r = np.sqrt(np.sum(pos**2, axis=1))
acc_mag = np.sqrt(np.sum(acc**2, axis=1))

# 1. XY projection colored by density
sc1 = ax1.scatter(pos[:, 0], pos[:, 1], c=dens, s=5, cmap='hot', vmin=0, vmax=2.0)
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_title('XY Projection (Final State)')
ax1.set_xlim(-1.5, 1.5)
ax1.set_ylim(-1.5, 1.5)
ax1.set_aspect('equal')
plt.colorbar(sc1, ax=ax1, label='Density')

# 2. Density profile
r_bins = np.linspace(0, 1.5, 30)
r_centers = 0.5 * (r_bins[1:] + r_bins[:-1])
dens_profile = []
for i in range(len(r_bins)-1):
    mask = (r >= r_bins[i]) & (r < r_bins[i+1])
    if np.sum(mask) > 0:
        dens_profile.append(np.mean(dens[mask]))
    else:
        dens_profile.append(0)

ax2.plot(r_centers, dens_profile, 'b-', linewidth=3, label='Mean')
ax2.scatter(r, dens, c='lightblue', s=2, alpha=0.3, label='Particles')
ax2.set_xlabel('Radius')
ax2.set_ylabel('Density')
ax2.set_title('Final Density Profile')
ax2.set_xlim(0, 1.5)
ax2.set_ylim(0, 2.0)
ax2.grid(True, alpha=0.3)
ax2.legend()

# 3. Force/Acceleration profile
acc_profile = []
for i in range(len(r_bins)-1):
    mask = (r >= r_bins[i]) & (r < r_bins[i+1])
    if np.sum(mask) > 0:
        acc_profile.append(np.mean(acc_mag[mask]))
    else:
        acc_profile.append(0)

ax3.plot(r_centers, acc_profile, 'r-', linewidth=3, label='Mean')
ax3.scatter(r, acc_mag, c='pink', s=2, alpha=0.3, label='Particles')
ax3.set_xlabel('Radius')
ax3.set_ylabel('|Acceleration|')
ax3.set_title('Final Force/Acceleration Profile')
ax3.set_xlim(0, 1.5)
ax3.grid(True, alpha=0.3)
ax3.legend()

# 4. 3D view colored by density
sc4 = ax4.scatter(pos[:, 0], pos[:, 1], pos[:, 2], c=dens, s=5, cmap='hot', vmin=0, vmax=2.0)
ax4.set_xlabel('x')
ax4.set_ylabel('y')
ax4.set_zlabel('z')
ax4.set_title('3D View (Final State)')
ax4.set_xlim(-1.5, 1.5)
ax4.set_ylim(-1.5, 1.5)
ax4.set_zlim(-1.5, 1.5)

# 5. Radial distribution
ax5.hist(r, bins=50, range=(0, 1.5), alpha=0.7, color='green')
ax5.set_xlabel('Radius')
ax5.set_ylabel('Particle Count')
ax5.set_title('Final Radial Distribution')
ax5.grid(True, alpha=0.3)

# 6. Acceleration magnitude distribution
ax6.hist(acc_mag, bins=50, alpha=0.7, color='orange')
ax6.set_xlabel('|Acceleration|')
ax6.set_ylabel('Count')
ax6.set_title('Final Acceleration Distribution')
ax6.set_yscale('log')
ax6.grid(True, alpha=0.3)

plt.tight_layout()
png_file = f'{animation_dir}/{output_prefix}_final.png'
plt.savefig(png_file, dpi=150, bbox_inches='tight')
print(f'✓ Saved: {png_file}')
print('\n✓ Done!')
