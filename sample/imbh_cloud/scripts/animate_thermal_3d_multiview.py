#!/usr/bin/env python3
"""
IMBH Tidal Disruption: Multi-View 3D Animation
===============================================
Creates comprehensive 3D animations with multiple viewing angles:
- Edge-on view (XY plane, looking along Z-axis)
- Face-on view (XZ plane, looking along Y-axis)
- Perspective view (3D rotational)

Shows thermal comparison (WITH vs WITHOUT cooling) for each view

Usage:
    python animate_thermal_3d_multiview.py <results_dir_cool> <results_dir_nocool> <output_dir>
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from mpl_toolkits.mplot3d import Axes3D
import sys
from pathlib import Path

# GALACTIC UNITS
TIME_UNIT = 0.978  # Myr per code unit

# IMBH Parameters (CORRECTED PHYSICS)
M_BH = 1e5  # M_☉
BH_POS = np.array([0.0, 0.0, 0.0])  # pc - IMBH is STATIONARY at origin
BH_VEL = np.array([0.0, 0.0, 0.0])  # km/s - IMBH does not move

# Cloud initial conditions (from configuration)
CLOUD_INITIAL_POS = np.array([-20.0, 3.0, 0.0])  # pc - Cloud starts here
CLOUD_INITIAL_VEL = np.array([10.0, 0.0, 0.0])   # km/s - Cloud moves toward IMBH
IMPACT_PARAM = 3.0  # pc

def load_snapshot(filepath):
    """Load a single CSV snapshot"""
    try:
        with open(filepath, 'r') as f:
            skiprows = 0
            for line in f:
                if line.startswith('id,'):
                    break
                skiprows += 1
        
        data = np.loadtxt(filepath, delimiter=',', skiprows=skiprows+1)
        
        # Positions are in kpc, convert to pc
        # Particles are at ~origin, need to shift to match cloud initial position
        return {
            'id': data[:, 0].astype(int),
            'pos': data[:, 1:4] * 1000.0,  # kpc -> pc
            'vel': data[:, 4:7],           # km/s
            'mass': data[:, 10],           # 10^10 M_☉
            'rho': data[:, 11],            # code density units
            'P': data[:, 12],              # code pressure units
            'u': data[:, 13],              # code internal energy units
            'h': data[:, 14] * 1000.0      # kpc -> pc
        }
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return None

def shift_cloud_to_initial_position(data, time):
    """
    Shift cloud particles from loaded position (~ origin) to correct initial position
    The cloud starts at (-20, 3, 0) pc and moves with v=(10, 0, 0) km/s
    """
    if data is None:
        return None
    
    # Calculate cloud center from loaded positions
    cloud_center_loaded = np.array([
        data['pos'][:, 0].mean(),
        data['pos'][:, 1].mean(),
        data['pos'][:, 2].mean()
    ])
    
    # Cloud position at given time (moving with constant velocity initially)
    cloud_center_at_time = CLOUD_INITIAL_POS + CLOUD_INITIAL_VEL * time
    
    # Shift all particles
    shift = cloud_center_at_time - cloud_center_loaded
    data['pos'] = data['pos'] + shift
    
    return data

def create_edge_on_animation(snapshots_cool, snapshots_nocool, output_dir):
    """Create edge-on view animation (XY plane, along tidal tail direction)"""
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    print("Creating edge-on view animation (XY projection)...")
    
    # Load all snapshots
    n_frames = min(len(snapshots_cool), len(snapshots_nocool))
    data_cool = []
    data_nocool = []
    times = []
    
    for i in range(n_frames):
        d_cool = load_snapshot(snapshots_cool[i])
        d_nocool = load_snapshot(snapshots_nocool[i])
        
        if d_cool is None or d_nocool is None:
            continue
        
        time_idx = int(snapshots_cool[i].stem.split('_')[1])
        time_code = time_idx * 0.02
        times.append(time_code)
        
        # Shift cloud to correct position
        d_cool = shift_cloud_to_initial_position(d_cool, time_code)
        d_nocool = shift_cloud_to_initial_position(d_nocool, time_code)
        
        data_cool.append(d_cool)
        data_nocool.append(d_nocool)
    
    if len(data_cool) == 0:
        print("No valid snapshots found!")
        return
    
    print(f"Loaded {len(data_cool)} valid snapshot pairs")
    
    # Determine plot limits (use final snapshot for max extent)
    final_pos = np.vstack([data_cool[-1]['pos'], data_nocool[-1]['pos']])
    xlim = [final_pos[:, 0].min() - 10, final_pos[:, 0].max() + 10]
    ylim = [final_pos[:, 1].min() - 10, final_pos[:, 1].max() + 10]
    
    # Density colormap range
    rho_all = np.concatenate([d['rho'] for d in data_cool + data_nocool])
    vmin_rho = np.log10(np.percentile(rho_all[rho_all > 0], 1))
    vmax_rho = np.log10(np.percentile(rho_all[rho_all > 0], 99))
    
    # Create figure
    fig, (ax_cool, ax_nocool) = plt.subplots(1, 2, figsize=(18, 7))
    
    def init():
        ax_cool.clear()
        ax_nocool.clear()
        return []
    
    def update(frame):
        ax_cool.clear()
        ax_nocool.clear()
        
        t = times[frame]
        t_myr = t * TIME_UNIT
        
        # WITH cooling
        pos_cool = data_cool[frame]['pos']
        rho_cool = data_cool[frame]['rho']
        log_rho_cool = np.log10(rho_cool + 1e-20)
        
        ax_cool.scatter(pos_cool[:, 0], pos_cool[:, 1],
                       c=log_rho_cool, cmap='viridis', s=3, alpha=0.7,
                       vmin=vmin_rho, vmax=vmax_rho)
        
        # Plot IMBH at origin
        ax_cool.scatter([BH_POS[0]], [BH_POS[1]], c='red', s=400, marker='*',
                       edgecolors='black', linewidth=2, zorder=10,
                       label=f'IMBH ({M_BH:.0e} M☉)')
        
        # Mark impact parameter line
        ax_cool.axhline(y=IMPACT_PARAM, color='red', linestyle=':', alpha=0.3, linewidth=1.5)
        ax_cool.text(xlim[1]-5, IMPACT_PARAM+0.5, f'b = {IMPACT_PARAM} pc', 
                    color='red', alpha=0.7, fontsize=10)
        
        ax_cool.set_xlabel('x [pc]', fontsize=13)
        ax_cool.set_ylabel('y [pc]', fontsize=13)
        ax_cool.set_xlim(xlim)
        ax_cool.set_ylim(ylim)
        ax_cool.set_aspect('equal')
        ax_cool.set_title(f'EDGE-ON VIEW: WITH Cooling\nt = {t_myr:.3f} Myr | N = {len(pos_cool)}', 
                         fontsize=13, fontweight='bold')
        ax_cool.legend(loc='upper left', fontsize=11)
        ax_cool.grid(True, alpha=0.3)
        
        # WITHOUT cooling
        pos_nocool = data_nocool[frame]['pos']
        rho_nocool = data_nocool[frame]['rho']
        log_rho_nocool = np.log10(rho_nocool + 1e-20)
        
        ax_nocool.scatter(pos_nocool[:, 0], pos_nocool[:, 1],
                         c=log_rho_nocool, cmap='viridis', s=3, alpha=0.7,
                         vmin=vmin_rho, vmax=vmax_rho)
        
        # Plot IMBH at origin
        ax_nocool.scatter([BH_POS[0]], [BH_POS[1]], c='red', s=400, marker='*',
                         edgecolors='black', linewidth=2, zorder=10,
                         label=f'IMBH ({M_BH:.0e} M☉)')
        
        # Mark impact parameter line
        ax_nocool.axhline(y=IMPACT_PARAM, color='red', linestyle=':', alpha=0.3, linewidth=1.5)
        ax_nocool.text(xlim[1]-5, IMPACT_PARAM+0.5, f'b = {IMPACT_PARAM} pc', 
                      color='red', alpha=0.7, fontsize=10)
        
        ax_nocool.set_xlabel('x [pc]', fontsize=13)
        ax_nocool.set_ylabel('y [pc]', fontsize=13)
        ax_nocool.set_xlim(xlim)
        ax_nocool.set_ylim(ylim)
        ax_nocool.set_aspect('equal')
        ax_nocool.set_title(f'EDGE-ON VIEW: WITHOUT Cooling\nt = {t_myr:.3f} Myr | N = {len(pos_nocool)}', 
                           fontsize=13, fontweight='bold')
        ax_nocool.legend(loc='upper left', fontsize=11)
        ax_nocool.grid(True, alpha=0.3)
        
        plt.tight_layout()
        return [ax_cool, ax_nocool]
    
    # Create animation
    anim = FuncAnimation(fig, update, frames=len(data_cool), init_func=init,
                        interval=100, blit=False)
    
    # Save
    output_file = output_path / 'thermal_comparison_edge_on.gif'
    print(f"Saving edge-on animation to {output_file}...")
    writer = PillowWriter(fps=10)
    anim.save(output_file, writer=writer, dpi=120)
    print(f"✓ Edge-on animation saved: {output_file}")
    
    plt.close()

def create_face_on_animation(snapshots_cool, snapshots_nocool, output_dir):
    """Create face-on view animation (XZ plane, perpendicular to impact parameter)"""
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    print("Creating face-on view animation (XZ projection)...")
    
    # Load all snapshots
    n_frames = min(len(snapshots_cool), len(snapshots_nocool))
    data_cool = []
    data_nocool = []
    times = []
    
    for i in range(n_frames):
        d_cool = load_snapshot(snapshots_cool[i])
        d_nocool = load_snapshot(snapshots_nocool[i])
        
        if d_cool is None or d_nocool is None:
            continue
        
        time_idx = int(snapshots_cool[i].stem.split('_')[1])
        time_code = time_idx * 0.02
        times.append(time_code)
        
        # Shift cloud to correct position
        d_cool = shift_cloud_to_initial_position(d_cool, time_code)
        d_nocool = shift_cloud_to_initial_position(d_nocool, time_code)
        
        data_cool.append(d_cool)
        data_nocool.append(d_nocool)
    
    if len(data_cool) == 0:
        print("No valid snapshots found!")
        return
    
    # Determine plot limits
    final_pos = np.vstack([data_cool[-1]['pos'], data_nocool[-1]['pos']])
    xlim = [final_pos[:, 0].min() - 10, final_pos[:, 0].max() + 10]
    zlim = [final_pos[:, 2].min() - 10, final_pos[:, 2].max() + 10]
    
    # Density colormap range
    rho_all = np.concatenate([d['rho'] for d in data_cool + data_nocool])
    vmin_rho = np.log10(np.percentile(rho_all[rho_all > 0], 1))
    vmax_rho = np.log10(np.percentile(rho_all[rho_all > 0], 99))
    
    # Create figure
    fig, (ax_cool, ax_nocool) = plt.subplots(1, 2, figsize=(18, 7))
    
    def init():
        ax_cool.clear()
        ax_nocool.clear()
        return []
    
    def update(frame):
        ax_cool.clear()
        ax_nocool.clear()
        
        t = times[frame]
        t_myr = t * TIME_UNIT
        
        # WITH cooling
        pos_cool = data_cool[frame]['pos']
        rho_cool = data_cool[frame]['rho']
        log_rho_cool = np.log10(rho_cool + 1e-20)
        
        ax_cool.scatter(pos_cool[:, 0], pos_cool[:, 2],
                       c=log_rho_cool, cmap='viridis', s=3, alpha=0.7,
                       vmin=vmin_rho, vmax=vmax_rho)
        
        # Plot IMBH at origin
        ax_cool.scatter([BH_POS[0]], [BH_POS[2]], c='red', s=400, marker='*',
                       edgecolors='black', linewidth=2, zorder=10,
                       label=f'IMBH ({M_BH:.0e} M☉)')
        
        ax_cool.set_xlabel('x [pc]', fontsize=13)
        ax_cool.set_ylabel('z [pc]', fontsize=13)
        ax_cool.set_xlim(xlim)
        ax_cool.set_ylim(zlim)
        ax_cool.set_aspect('equal')
        ax_cool.set_title(f'FACE-ON VIEW: WITH Cooling\nt = {t_myr:.3f} Myr | N = {len(pos_cool)}', 
                         fontsize=13, fontweight='bold')
        ax_cool.legend(loc='upper left', fontsize=11)
        ax_cool.grid(True, alpha=0.3)
        
        # WITHOUT cooling
        pos_nocool = data_nocool[frame]['pos']
        rho_nocool = data_nocool[frame]['rho']
        log_rho_nocool = np.log10(rho_nocool + 1e-20)
        
        ax_nocool.scatter(pos_nocool[:, 0], pos_nocool[:, 2],
                         c=log_rho_nocool, cmap='viridis', s=3, alpha=0.7,
                         vmin=vmin_rho, vmax=vmax_rho)
        
        # Plot IMBH at origin
        ax_nocool.scatter([BH_POS[0]], [BH_POS[2]], c='red', s=400, marker='*',
                         edgecolors='black', linewidth=2, zorder=10,
                         label=f'IMBH ({M_BH:.0e} M☉)')
        
        ax_nocool.set_xlabel('x [pc]', fontsize=13)
        ax_nocool.set_ylabel('z [pc]', fontsize=13)
        ax_nocool.set_xlim(xlim)
        ax_nocool.set_ylim(zlim)
        ax_nocool.set_aspect('equal')
        ax_nocool.set_title(f'FACE-ON VIEW: WITHOUT Cooling\nt = {t_myr:.3f} Myr | N = {len(pos_nocool)}', 
                           fontsize=13, fontweight='bold')
        ax_nocool.legend(loc='upper left', fontsize=11)
        ax_nocool.grid(True, alpha=0.3)
        
        plt.tight_layout()
        return [ax_cool, ax_nocool]
    
    # Create animation
    anim = FuncAnimation(fig, update, frames=len(data_cool), init_func=init,
                        interval=100, blit=False)
    
    # Save
    output_file = output_path / 'thermal_comparison_face_on.gif'
    print(f"Saving face-on animation to {output_file}...")
    writer = PillowWriter(fps=10)
    anim.save(output_file, writer=writer, dpi=120)
    print(f"✓ Face-on animation saved: {output_file}")
    
    plt.close()

def create_3d_rotating_animation(snapshots_cool, snapshots_nocool, output_dir):
    """Create 3D perspective animation with slow rotation"""
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    print("Creating 3D rotating animation...")
    
    # Load all snapshots
    n_frames = min(len(snapshots_cool), len(snapshots_nocool))
    data_cool = []
    data_nocool = []
    times = []
    
    for i in range(n_frames):
        d_cool = load_snapshot(snapshots_cool[i])
        d_nocool = load_snapshot(snapshots_nocool[i])
        
        if d_cool is None or d_nocool is None:
            continue
        
        time_idx = int(snapshots_cool[i].stem.split('_')[1])
        time_code = time_idx * 0.02
        times.append(time_code)
        
        # Shift cloud to correct position
        d_cool = shift_cloud_to_initial_position(d_cool, time_code)
        d_nocool = shift_cloud_to_initial_position(d_nocool, time_code)
        
        data_cool.append(d_cool)
        data_nocool.append(d_nocool)
    
    if len(data_cool) == 0:
        print("No valid snapshots found!")
        return
    
    # Determine plot limits
    final_pos = np.vstack([data_cool[-1]['pos'], data_nocool[-1]['pos']])
    xlim = [final_pos[:, 0].min() - 10, final_pos[:, 0].max() + 10]
    ylim = [final_pos[:, 1].min() - 10, final_pos[:, 1].max() + 10]
    zlim = [final_pos[:, 2].min() - 10, final_pos[:, 2].max() + 10]
    
    # Density colormap range
    rho_all = np.concatenate([d['rho'] for d in data_cool + data_nocool])
    vmin_rho = np.percentile(rho_all, 1)
    vmax_rho = np.percentile(rho_all, 99)
    
    # Create figure
    fig = plt.figure(figsize=(18, 8))
    ax_cool = fig.add_subplot(121, projection='3d')
    ax_nocool = fig.add_subplot(122, projection='3d')
    
    def init():
        ax_cool.clear()
        ax_nocool.clear()
        return []
    
    def update(frame):
        ax_cool.clear()
        ax_nocool.clear()
        
        t = times[frame]
        t_myr = t * TIME_UNIT
        
        # Rotation angle
        azim = 45 + frame * 1.5  # Slow rotation
        
        # WITH cooling
        pos_cool = data_cool[frame]['pos']
        rho_cool = data_cool[frame]['rho']
        
        ax_cool.scatter(pos_cool[:, 0], pos_cool[:, 1], pos_cool[:, 2],
                       c=rho_cool, cmap='viridis', s=1, alpha=0.6,
                       vmin=vmin_rho, vmax=vmax_rho)
        
        # Plot IMBH at origin
        ax_cool.scatter([BH_POS[0]], [BH_POS[1]], [BH_POS[2]], 
                       c='red', s=300, marker='*', edgecolors='black', linewidth=2,
                       label=f'IMBH ({M_BH:.0e} M☉)')
        
        ax_cool.set_xlabel('x [pc]')
        ax_cool.set_ylabel('y [pc]')
        ax_cool.set_zlabel('z [pc]')
        ax_cool.set_xlim(xlim)
        ax_cool.set_ylim(ylim)
        ax_cool.set_zlim(zlim)
        ax_cool.set_title(f'3D VIEW: WITH Cooling\nt = {t_myr:.3f} Myr', fontsize=12, fontweight='bold')
        ax_cool.legend(loc='upper left')
        ax_cool.view_init(elev=20, azim=azim)
        
        # WITHOUT cooling
        pos_nocool = data_nocool[frame]['pos']
        rho_nocool = data_nocool[frame]['rho']
        
        ax_nocool.scatter(pos_nocool[:, 0], pos_nocool[:, 1], pos_nocool[:, 2],
                         c=rho_nocool, cmap='viridis', s=1, alpha=0.6,
                         vmin=vmin_rho, vmax=vmax_rho)
        
        # Plot IMBH at origin
        ax_nocool.scatter([BH_POS[0]], [BH_POS[1]], [BH_POS[2]], 
                         c='red', s=300, marker='*', edgecolors='black', linewidth=2,
                         label=f'IMBH ({M_BH:.0e} M☉)')
        
        ax_nocool.set_xlabel('x [pc]')
        ax_nocool.set_ylabel('y [pc]')
        ax_nocool.set_zlabel('z [pc]')
        ax_nocool.set_xlim(xlim)
        ax_nocool.set_ylim(ylim)
        ax_nocool.set_zlim(zlim)
        ax_nocool.set_title(f'3D VIEW: WITHOUT Cooling\nt = {t_myr:.3f} Myr', fontsize=12, fontweight='bold')
        ax_nocool.legend(loc='upper left')
        ax_nocool.view_init(elev=20, azim=azim)
        
        plt.tight_layout()
        return [ax_cool, ax_nocool]
    
    # Create animation
    anim = FuncAnimation(fig, update, frames=len(data_cool), init_func=init,
                        interval=100, blit=False)
    
    # Save
    output_file = output_path / 'thermal_comparison_3d_rotating.gif'
    print(f"Saving 3D rotating animation to {output_file}...")
    writer = PillowWriter(fps=10)
    anim.save(output_file, writer=writer, dpi=100)
    print(f"✓ 3D rotating animation saved: {output_file}")
    
    plt.close()

def main():
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)
    
    results_cool = sys.argv[1]
    results_nocool = sys.argv[2]
    output_dir = sys.argv[3]
    
    # Find snapshots
    snapshots_cool = sorted(Path(results_cool).glob('snapshot_*.csv'))
    snapshots_nocool = sorted(Path(results_nocool).glob('snapshot_*.csv'))
    
    if not snapshots_cool or not snapshots_nocool:
        print("Error: No snapshots found!")
        print(f"  WITH cooling: {len(snapshots_cool)} snapshots in {results_cool}")
        print(f"  WITHOUT cooling: {len(snapshots_nocool)} snapshots in {results_nocool}")
        sys.exit(1)
    
    print("=" * 70)
    print("IMBH TIDAL DISRUPTION: MULTI-VIEW 3D ANIMATION")
    print("=" * 70)
    print(f"WITH cooling:    {results_cool} ({len(snapshots_cool)} snapshots)")
    print(f"WITHOUT cooling: {results_nocool} ({len(snapshots_nocool)} snapshots)")
    print(f"Output:          {output_dir}")
    print("=" * 70)
    print("\nPhysics Setup (CORRECTED):")
    print(f"  • IMBH: Stationary at origin [0, 0, 0] pc")
    print(f"  • Cloud: Initial position [{CLOUD_INITIAL_POS[0]}, {CLOUD_INITIAL_POS[1]}, {CLOUD_INITIAL_POS[2]}] pc")
    print(f"  • Cloud: Initial velocity [{CLOUD_INITIAL_VEL[0]}, {CLOUD_INITIAL_VEL[1]}, {CLOUD_INITIAL_VEL[2]}] km/s")
    print(f"  • Impact parameter b = {IMPACT_PARAM} pc")
    print("=" * 70)
    print()
    
    # Create all three views
    create_edge_on_animation(snapshots_cool, snapshots_nocool, output_dir)
    create_face_on_animation(snapshots_cool, snapshots_nocool, output_dir)
    create_3d_rotating_animation(snapshots_cool, snapshots_nocool, output_dir)
    
    print("\n" + "=" * 70)
    print("✓ ALL MULTI-VIEW ANIMATIONS COMPLETE!")
    print("=" * 70)
    print(f"\nOutput files in {output_dir}:")
    print("  - thermal_comparison_edge_on.gif      : XY projection (tidal tail view)")
    print("  - thermal_comparison_face_on.gif      : XZ projection (perpendicular view)")
    print("  - thermal_comparison_3d_rotating.gif  : 3D perspective with rotation")

if __name__ == '__main__':
    main()
