#!/usr/bin/env python3
"""
IMBH Tidal Disruption - 3D Animation with Black Hole
=====================================================

Creates a rotating 3D animation showing the full spatial arrangement:
- SPH particles colored by density/temperature/velocity
- IMBH position marked clearly at origin
- Rotating view to show 3D structure
- Optional: trajectory path, tidal radius indicator

Units: IMBH_ENCOUNTER ([L]=pc, [M]=1000 M_sun, [V]=km/s, [T]=0.978 Myr)

Usage:
    python animate_3d_view.py <results_dir> [-o output.gif] [--mode density]
    python animate_3d_view.py <results_dir> --mode velocity --rotate
    python animate_3d_view.py <results_dir> --mode temperature --static-view
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.colors import LogNorm, Normalize, LinearSegmentedColormap
import matplotlib.patheffects as path_effects
import argparse
from pathlib import Path
import glob
import warnings
import json
warnings.filterwarnings('ignore', category=RuntimeWarning)

# =============================================================================
# HIGH CONTRAST DARK THEME
# =============================================================================

DARK_BG = '#0d0d12'
DARK_PANEL = '#15151f'
GRID_COLOR = '#3a3a4a'
TEXT_COLOR = '#f0f0f0'
ACCENT_COLOR = '#00ffff'
BH_COLOR = '#ff4444'
TRAJECTORY_COLOR = '#ffaa00'

# Physical constants
TIME_UNIT = 0.978  # Myr
GAMMA = 5.0 / 3.0
MU = 2.3
M_H = 1.67e-24
K_B = 1.38e-16
VELOCITY_TO_CGS = 1e5
TEMP_CONVERSION = MU * M_H / K_B * VELOCITY_TO_CGS**2

# =============================================================================
# COLORMAPS
# =============================================================================

def create_high_contrast_colormaps():
    """Create high-contrast, perceptually uniform colormaps for SPH visualization."""
    
    # Density: Dark blue -> Cyan -> Green -> Yellow -> Orange -> White
    # Good for log-scale density, emphasizes high-density regions
    density_colors = [
        (0.0, '#0d0221'),   # Deep purple-black (low density void)
        (0.15, '#1a237e'),  # Deep blue
        (0.30, '#0277bd'),  # Ocean blue
        (0.45, '#00acc1'),  # Cyan
        (0.55, '#26a69a'),  # Teal
        (0.65, '#66bb6a'),  # Green
        (0.75, '#ffeb3b'),  # Yellow
        (0.85, '#ff9800'),  # Orange
        (0.92, '#f44336'),  # Red
        (1.0, '#ffffff'),   # White (highest density)
    ]
    
    # Temperature: Deep blue -> Blue -> Cyan -> Green -> Yellow -> Orange -> Red -> Pink
    # Classic "cool to hot" that's intuitive for temperature
    temp_colors = [
        (0.0, '#000033'),   # Deep navy (coldest)
        (0.12, '#0d47a1'),  # Dark blue
        (0.25, '#1976d2'),  # Blue
        (0.35, '#00bcd4'),  # Cyan
        (0.45, '#4caf50'),  # Green
        (0.55, '#8bc34a'),  # Light green
        (0.65, '#ffeb3b'),  # Yellow
        (0.75, '#ff9800'),  # Orange
        (0.85, '#f44336'),  # Red
        (0.92, '#e91e63'),  # Pink
        (1.0, '#ffffff'),   # White (hottest)
    ]
    
    # Velocity: Purple -> Blue -> Cyan -> Green -> Yellow -> Red
    # Diverging-style for velocity magnitude
    velocity_colors = [
        (0.0, '#1a237e'),   # Deep indigo (stationary)
        (0.2, '#3949ab'),   # Indigo
        (0.35, '#5c6bc0'),  # Light indigo
        (0.5, '#26c6da'),   # Cyan (moderate)
        (0.65, '#66bb6a'),  # Green
        (0.8, '#ffca28'),   # Amber
        (0.9, '#ff7043'),   # Deep orange
        (1.0, '#d32f2f'),   # Red (fast)
    ]
    
    def make_cmap(colors, name):
        positions = [c[0] for c in colors]
        hex_colors = [c[1] for c in colors]
        rgb_colors = []
        for h in hex_colors:
            h = h.lstrip('#')
            rgb_colors.append(tuple(int(h[i:i+2], 16)/255 for i in (0, 2, 4)))
        
        cdict = {'red': [], 'green': [], 'blue': []}
        for pos, (r, g, b) in zip(positions, rgb_colors):
            cdict['red'].append((pos, r, r))
            cdict['green'].append((pos, g, g))
            cdict['blue'].append((pos, b, b))
        
        return LinearSegmentedColormap(name, cdict, N=256)
    
    return {
        'density': make_cmap(density_colors, 'density_3d'),
        'temperature': make_cmap(temp_colors, 'temp_3d'),
        'velocity': make_cmap(velocity_colors, 'vel_3d'),
    }

CUSTOM_CMAPS = create_high_contrast_colormaps()


# =============================================================================
# DATA LOADING
# =============================================================================

def load_snapshot(filepath):
    """Load a single CSV snapshot."""
    try:
        with open(filepath, 'r') as f:
            skip_lines = 0
            for line in f:
                if line.startswith('#'):
                    skip_lines += 1
                else:
                    break
        
        data = np.genfromtxt(filepath, delimiter=',', skip_header=skip_lines, names=True)
        
        snapshot = {
            'pos_x': data['pos_x'],
            'pos_y': data['pos_y'],
            'pos_z': data['pos_z'],
            'vel_x': data['vel_x'],
            'vel_y': data['vel_y'],
            'vel_z': data['vel_z'],
            'dens': data['dens'],
            'pres': data['pres'],
            'mass': data['mass'],
            'sml': data['sml'],
            'sound': data['sound'],
        }
        
        return snapshot
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return None


def load_config(results_dir):
    """Try to load config from results directory or parent."""
    # Try to find config.json in results or nearby
    config_paths = [
        Path(results_dir) / 'config.json',
        Path(results_dir).parent / 'config.json',
    ]
    
    for config_path in config_paths:
        if config_path.exists():
            try:
                with open(config_path) as f:
                    return json.load(f)
            except Exception:
                pass
    
    return None


# =============================================================================
# PHYSICS
# =============================================================================

def compute_temperature(pres, dens):
    """Compute temperature from pressure and density."""
    specific_energy = pres / (dens + 1e-30)
    temperature = specific_energy * TEMP_CONVERSION / GAMMA
    return temperature


def compute_velocity_magnitude(vel_x, vel_y, vel_z, mass=None):
    """Compute velocity magnitude relative to center of mass."""
    if mass is None:
        mass = np.ones_like(vel_x)
    
    total_mass = np.sum(mass)
    v_com_x = np.sum(mass * vel_x) / total_mass
    v_com_y = np.sum(mass * vel_y) / total_mass
    v_com_z = np.sum(mass * vel_z) / total_mass
    
    dvel_x = vel_x - v_com_x
    dvel_y = vel_y - v_com_y
    dvel_z = vel_z - v_com_z
    
    return np.sqrt(dvel_x**2 + dvel_y**2 + dvel_z**2)


def compute_center_of_mass(pos_x, pos_y, pos_z, mass):
    """Compute center of mass position."""
    total_mass = np.sum(mass)
    com_x = np.sum(mass * pos_x) / total_mass
    com_y = np.sum(mass * pos_y) / total_mass
    com_z = np.sum(mass * pos_z) / total_mass
    return com_x, com_y, com_z


# =============================================================================
# 3D ANIMATION
# =============================================================================

def create_3d_animation(results_dir, output_file, mode='density', fps=15,
                        rotate=False, elevation=25, show_trajectory=True,
                        bh_position=(0, 0, 0), bh_mass=100, downsample=1):
    """
    Create 3D animation of IMBH-cloud encounter.
    
    Parameters:
    -----------
    results_dir : str
        Directory containing snapshot CSV files
    output_file : str
        Output GIF filename
    mode : str
        Color mode: 'density', 'temperature', 'velocity'
    fps : int
        Frames per second
    rotate : bool
        Enable rotating view (azimuth changes with time)
    elevation : float
        Camera elevation angle in degrees
    show_trajectory : bool
        Show cloud trajectory path
    bh_position : tuple
        Black hole position (x, y, z) in code units
    bh_mass : float
        Black hole mass in code units (for tidal radius)
    downsample : int
        Downsample particles by this factor (for speed)
    """
    
    snapshot_files = sorted(glob.glob(f"{results_dir}/snapshot_*.csv"))
    if not snapshot_files:
        print(f"Error: No snapshots found in {results_dir}")
        return False
    
    print(f"Found {len(snapshot_files)} snapshots")
    
    # Load first snapshot to determine bounds
    first = load_snapshot(snapshot_files[0])
    if first is None:
        return False
    
    # Calculate bounds including BH position
    all_x = [first['pos_x'].min(), first['pos_x'].max(), bh_position[0]]
    all_y = [first['pos_y'].min(), first['pos_y'].max(), bh_position[1]]
    all_z = [first['pos_z'].min(), first['pos_z'].max(), bh_position[2]]
    
    # Load last snapshot for extended bounds
    last = load_snapshot(snapshot_files[-1])
    if last is not None:
        all_x.extend([last['pos_x'].min(), last['pos_x'].max()])
        all_y.extend([last['pos_y'].min(), last['pos_y'].max()])
        all_z.extend([last['pos_z'].min(), last['pos_z'].max()])
    
    pad = 2.0
    xlim = (min(all_x) - pad, max(all_x) + pad)
    ylim = (min(all_y) - pad, max(all_y) + pad)
    zlim = (min(all_z) - pad, max(all_z) + pad)
    
    # Make symmetric around origin for better view
    max_range = max(abs(xlim[0]), abs(xlim[1]), abs(ylim[0]), abs(ylim[1]), 
                    abs(zlim[0]), abs(zlim[1]))
    xlim = (-max_range, max_range)
    ylim = (-max_range, max_range)
    zlim = (-max_range * 0.3, max_range * 0.3)  # Z is typically smaller
    
    # Pre-compute GLOBAL color ranges from all snapshots for fixed colorbar
    print("Computing global color ranges from all snapshots...")
    global_dens_min, global_dens_max = np.inf, -np.inf
    global_temp_min, global_temp_max = np.inf, -np.inf
    global_vel_min, global_vel_max = np.inf, -np.inf
    
    for i, snap_file in enumerate(snapshot_files):
        snap = load_snapshot(snap_file)
        if snap is not None:
            # Density range
            dens = snap['dens']
            global_dens_min = min(global_dens_min, np.percentile(dens, 1))
            global_dens_max = max(global_dens_max, np.percentile(dens, 99))
            
            # Temperature range
            temp = compute_temperature(snap['pres'], snap['dens'])
            global_temp_min = min(global_temp_min, np.percentile(temp, 1))
            global_temp_max = max(global_temp_max, np.percentile(temp, 99))
            
            # Velocity range
            vel_mag = compute_velocity_magnitude(
                snap['vel_x'], snap['vel_y'], snap['vel_z'], snap['mass'])
            global_vel_min = min(global_vel_min, np.percentile(vel_mag, 1))
            global_vel_max = max(global_vel_max, np.percentile(vel_mag, 99))
        
        if (i + 1) % 10 == 0:
            print(f"  Scanned {i + 1}/{len(snapshot_files)} snapshots...")
    
    # Add small padding to ranges
    global_dens_min = max(global_dens_min * 0.5, 1e-3)
    global_dens_max = global_dens_max * 2.0
    global_temp_min = max(global_temp_min * 0.5, 1.0)
    global_temp_max = global_temp_max * 2.0
    global_vel_max = max(global_vel_max * 1.2, 1.0)
    
    print(f"  Density range: [{global_dens_min:.2e}, {global_dens_max:.2e}] code")
    print(f"  Temperature range: [{global_temp_min:.1f}, {global_temp_max:.1f}] K")
    print(f"  Velocity range: [0, {global_vel_max:.2f}] km/s")
    
    # Compute tidal radius for reference
    # r_tidal = (M_cloud / M_BH)^(1/3) * r_peri
    # With M_cloud ~ 1 (1000 Msun), M_BH = 100 (1e5 Msun)
    r_tidal = (1.0 / bh_mass) ** (1/3) * 3.0  # ~3.6 pc for b=3pc
    
    # Setup figure
    fig = plt.figure(figsize=(14, 12), facecolor=DARK_BG)
    ax = fig.add_subplot(111, projection='3d', facecolor=DARK_BG)
    
    # Configure 3D axes
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_zlim(zlim)
    ax.set_xlabel('X (pc)', fontsize=12, color=TEXT_COLOR, labelpad=10)
    ax.set_ylabel('Y (pc)', fontsize=12, color=TEXT_COLOR, labelpad=10)
    ax.set_zlabel('Z (pc)', fontsize=12, color=TEXT_COLOR, labelpad=10)
    
    # Dark theme for 3D
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor(GRID_COLOR)
    ax.yaxis.pane.set_edgecolor(GRID_COLOR)
    ax.zaxis.pane.set_edgecolor(GRID_COLOR)
    ax.tick_params(colors=TEXT_COLOR, labelsize=9)
    ax.xaxis.line.set_color(GRID_COLOR)
    ax.yaxis.line.set_color(GRID_COLOR)
    ax.zaxis.line.set_color(GRID_COLOR)
    
    # Initial view
    ax.view_init(elev=elevation, azim=-60)
    
    # Setup colormap and normalization using GLOBAL ranges (fixed across all frames)
    if mode == 'density':
        cmap = CUSTOM_CMAPS['density']
        norm = LogNorm(vmin=global_dens_min, vmax=global_dens_max)
        cbar_label = 'Density (code units)'
        title_mode = 'Density'
    elif mode == 'temperature':
        cmap = CUSTOM_CMAPS['temperature']
        norm = LogNorm(vmin=global_temp_min, vmax=global_temp_max)
        cbar_label = 'Temperature (K)'
        title_mode = 'Temperature'
    elif mode == 'velocity':
        cmap = CUSTOM_CMAPS['velocity']
        norm = Normalize(vmin=0, vmax=global_vel_max)
        cbar_label = 'Internal Velocity (km/s)'
        title_mode = 'Velocity'
    else:
        cmap = CUSTOM_CMAPS['density']
        norm = LogNorm(vmin=global_dens_min, vmax=global_dens_max)
        cbar_label = 'Density'
        title_mode = 'Density'
    
    # Create scatter for particles (will update)
    scatter = ax.scatter([], [], [], s=5, c=[], cmap=cmap, norm=norm, 
                         alpha=0.7, edgecolors='none', depthshade=True)
    
    # Add colorbar
    cbar = fig.colorbar(scatter, ax=ax, shrink=0.6, pad=0.1, aspect=20)
    cbar.set_label(cbar_label, fontsize=12, color=TEXT_COLOR)
    cbar.ax.tick_params(colors=TEXT_COLOR, labelsize=10)
    cbar.outline.set_edgecolor(GRID_COLOR)
    
    # Draw IMBH at origin with glowing effect
    bh_x, bh_y, bh_z = bh_position
    
    # Glow rings around BH
    for size, alpha in [(800, 0.05), (600, 0.1), (400, 0.15), (250, 0.25)]:
        ax.scatter([bh_x], [bh_y], [bh_z], s=size, c=BH_COLOR, alpha=alpha, 
                   edgecolors='none', depthshade=False)
    
    # Central BH marker
    ax.scatter([bh_x], [bh_y], [bh_z], s=200, c='white', 
               marker='*', edgecolors=BH_COLOR, linewidths=2,
               depthshade=False, zorder=100)
    
    # Draw tidal radius circle (in XY plane at z=0)
    theta_circle = np.linspace(0, 2*np.pi, 100)
    tidal_x = bh_x + r_tidal * np.cos(theta_circle)
    tidal_y = bh_y + r_tidal * np.sin(theta_circle)
    tidal_z = np.zeros_like(theta_circle) + bh_z
    ax.plot(tidal_x, tidal_y, tidal_z, '--', color=ACCENT_COLOR, alpha=0.5, 
            linewidth=1.5, label=f'Tidal radius ({r_tidal:.1f} pc)')
    
    # Pre-compute full COM orbit trajectory for all snapshots
    print("Pre-computing COM orbit trajectory...")
    full_trajectory_x = []
    full_trajectory_y = []
    full_trajectory_z = []
    for snap_file in snapshot_files:
        snap = load_snapshot(snap_file)
        if snap is not None:
            com_x, com_y, com_z = compute_center_of_mass(
                snap['pos_x'], snap['pos_y'], snap['pos_z'], snap['mass'])
            full_trajectory_x.append(com_x)
            full_trajectory_y.append(com_y)
            full_trajectory_z.append(com_z)
    
    # Draw the full COM orbit line (static, shows complete path)
    if show_trajectory and len(full_trajectory_x) > 1:
        ax.plot(full_trajectory_x, full_trajectory_y, full_trajectory_z, 
                '-', color=TRAJECTORY_COLOR, alpha=0.4, linewidth=1.5, 
                label='COM orbit (full)')
    
    # Trajectory marker for current position (will be updated)
    trajectory_marker, = ax.plot([], [], [], 'o', color=TRAJECTORY_COLOR, 
                                  markersize=8, markeredgecolor='white',
                                  markeredgewidth=1.5, label='Cloud COM')
    
    # Partial trajectory line showing path up to current time
    trajectory_line, = ax.plot([], [], [], '-', color=TRAJECTORY_COLOR, 
                                alpha=0.9, linewidth=2.5)
    
    # Text annotations
    ax.set_title(f'IMBH-Cloud Encounter - {title_mode} (3D View)', 
                 fontsize=16, fontweight='bold', color=TEXT_COLOR, pad=20)
    
    time_text = ax.text2D(0.02, 0.98, '', transform=ax.transAxes, fontsize=14,
                          fontweight='bold', va='top', color=TEXT_COLOR,
                          path_effects=[path_effects.withStroke(linewidth=3, foreground=DARK_BG)])
    
    info_text = ax.text2D(0.02, 0.02, '', transform=ax.transAxes, fontsize=10,
                          va='bottom', color=TEXT_COLOR, family='monospace',
                          path_effects=[path_effects.withStroke(linewidth=2, foreground=DARK_BG)])
    
    # Legend
    ax.legend(loc='upper right', fontsize=10, facecolor=DARK_PANEL, 
              edgecolor=GRID_COLOR, labelcolor=TEXT_COLOR)
    
    # Add BH label
    ax.text(bh_x, bh_y, bh_z + max_range * 0.15, 'IMBH\n(10⁵ M☉)', 
            fontsize=11, color=BH_COLOR, ha='center', va='bottom',
            fontweight='bold')
    
    plt.tight_layout()
    
    print(f"Creating 3D {title_mode.lower()} animation...")
    print(f"  Frames: {len(snapshot_files)}")
    print(f"  FPS: {fps}")
    print(f"  Output: {output_file}")
    if rotate:
        print("  Rotation: enabled (360° over animation)")
    
    def update(frame_idx):
        snapshot = load_snapshot(snapshot_files[frame_idx])
        if snapshot is None:
            return scatter, trajectory_line, trajectory_marker
        
        # Downsample if needed
        n_particles = len(snapshot['pos_x'])
        if downsample > 1:
            idx = np.arange(0, n_particles, downsample)
        else:
            idx = np.arange(n_particles)
        
        # Get positions
        x = snapshot['pos_x'][idx]
        y = snapshot['pos_y'][idx]
        z = snapshot['pos_z'][idx]
        
        # Compute color data
        if mode == 'density':
            color_data = snapshot['dens'][idx]
        elif mode == 'temperature':
            color_data = compute_temperature(snapshot['pres'][idx], snapshot['dens'][idx])
        elif mode == 'velocity':
            color_data = compute_velocity_magnitude(
                snapshot['vel_x'][idx], snapshot['vel_y'][idx], snapshot['vel_z'][idx],
                snapshot['mass'][idx])
        else:
            color_data = snapshot['dens'][idx]
        
        # Update scatter
        scatter._offsets3d = (x, y, z)
        scatter.set_array(color_data)
        
        # Get current COM from pre-computed trajectory
        com_x = full_trajectory_x[frame_idx]
        com_y = full_trajectory_y[frame_idx]
        com_z = full_trajectory_z[frame_idx]
        
        # Update trajectory line (path up to current time)
        if show_trajectory and frame_idx > 0:
            trajectory_line.set_data(
                np.array(full_trajectory_x[:frame_idx+1]), 
                np.array(full_trajectory_y[:frame_idx+1]))
            trajectory_line.set_3d_properties(np.array(full_trajectory_z[:frame_idx+1]))
            # Update COM marker position
            trajectory_marker.set_data([com_x], [com_y])
            trajectory_marker.set_3d_properties([com_z])
        
        # Compute distance to BH
        r_to_bh = np.sqrt((com_x - bh_x)**2 + (com_y - bh_y)**2 + (com_z - bh_z)**2)
        
        # Update time text
        time_code = frame_idx * 0.02  # Assuming outputTime=0.02
        time_myr = time_code * TIME_UNIT
        time_text.set_text(f't = {time_myr:.3f} Myr')
        
        # Update info text
        info_str = f'Cloud CoM: ({com_x:.1f}, {com_y:.1f}, {com_z:.1f}) pc\n'
        info_str += f'Distance to BH: {r_to_bh:.2f} pc\n'
        info_str += f'Particles: {n_particles:,}'
        info_text.set_text(info_str)
        
        # Rotate view if enabled
        if rotate:
            azim = -60 + (frame_idx / len(snapshot_files)) * 120  # Rotate 120° total
            ax.view_init(elev=elevation, azim=azim)
        
        if (frame_idx + 1) % 5 == 0 or frame_idx == 0:
            print(f"  Frame {frame_idx + 1}/{len(snapshot_files)}: t={time_myr:.3f} Myr, r_BH={r_to_bh:.2f} pc")
        
        return scatter, trajectory_line, trajectory_marker
    
    # Create animation
    anim = FuncAnimation(fig, update, frames=len(snapshot_files), 
                         interval=1000/fps, blit=False)
    
    # Save
    print("Saving animation...")
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    writer = PillowWriter(fps=fps)
    anim.save(output_file, writer=writer, dpi=100)
    
    plt.close(fig)
    print(f"✓ 3D animation saved: {output_file}")
    
    return True


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='IMBH Tidal Disruption - 3D Animation with Black Hole',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python animate_3d_view.py results/Mc1e3_Mbh1e5_b3_v10/adiabatic_61k_gsph
  python animate_3d_view.py results/ -o animations/3d_density.gif --mode density
  python animate_3d_view.py results/ --mode temperature --no-rotate --elev 45
  python animate_3d_view.py results/ --mode velocity --downsample 2 --fps 20
        """
    )
    
    parser.add_argument('results_dir', type=str,
                        help='Directory containing snapshot_XXXX.csv files')
    parser.add_argument('-o', '--output', type=str, default=None,
                        help='Output GIF file (default: <results_dir>/3d_<mode>.gif)')
    parser.add_argument('--mode', type=str, default='density',
                        choices=['density', 'temperature', 'velocity'],
                        help='Coloring mode (default: density)')
    parser.add_argument('--fps', type=int, default=15,
                        help='Frames per second (default: 15)')
    parser.add_argument('--rotate', action='store_true',
                        help='Enable view rotation (default: no rotation)')
    parser.add_argument('--elev', type=float, default=25,
                        help='Camera elevation angle in degrees (default: 25)')
    parser.add_argument('--no-trajectory', action='store_true',
                        help='Disable trajectory path')
    parser.add_argument('--downsample', type=int, default=1,
                        help='Downsample particles by this factor (default: 1)')
    parser.add_argument('--bh-mass', type=float, default=100.0,
                        help='Black hole mass in code units (default: 100 = 1e5 Msun)')
    parser.add_argument('--bh-x', type=float, default=0.0,
                        help='Black hole X position (default: 0)')
    parser.add_argument('--bh-y', type=float, default=0.0,
                        help='Black hole Y position (default: 0)')
    parser.add_argument('--bh-z', type=float, default=0.0,
                        help='Black hole Z position (default: 0)')
    
    args = parser.parse_args()
    
    # Determine output file
    if args.output:
        output_file = args.output
    else:
        output_file = str(Path(args.results_dir) / f'3d_{args.mode}.gif')
    
    print("=================================================================")
    print("IMBH Tidal Disruption - 3D Animation with Black Hole")
    print("=================================================================")
    print(f"Results: {args.results_dir}")
    print(f"Output:  {output_file}")
    print(f"Mode:    {args.mode}")
    print(f"BH pos:  ({args.bh_x}, {args.bh_y}, {args.bh_z})")
    print(f"BH mass: {args.bh_mass} (code) = {args.bh_mass * 1000:.0e} Msun")
    print("=================================================================")
    
    success = create_3d_animation(
        results_dir=args.results_dir,
        output_file=output_file,
        mode=args.mode,
        fps=args.fps,
        rotate=args.rotate,
        elevation=args.elev,
        show_trajectory=not args.no_trajectory,
        bh_position=(args.bh_x, args.bh_y, args.bh_z),
        bh_mass=args.bh_mass,
        downsample=args.downsample,
    )
    
    if success:
        print("")
        print("=================================================================")
        print("✓ 3D animation complete!")
        print("=================================================================")
    else:
        print("")
        print("=================================================================")
        print("✗ Animation failed!")
        print("=================================================================")
        return 1
    
    return 0


if __name__ == '__main__':
    exit(main())
