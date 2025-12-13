#!/usr/bin/env python3
"""
IMBH Tidal Disruption - Shock Diagnostics Animation
====================================================

Explicitly visualizes shock wave propagation during tidal disruption:

1. VERTICAL COMPRESSION SHOCK (Z-direction):
   - Tidal force compresses cloud vertically
   - Forms pancake structure
   - 1D plot: P(z), T(z), ρ(z) for particles in central column
   
2. IN-PLANE TIDAL STRETCHING (X-Y plane):
   - Differential acceleration stretches cloud
   - Creates tidal tails and potential in-plane shocks
   - 1D plot: P(x), T(x), ρ(x) along tidal axis

3. SHOCK TEMPERATURE MAP:
   - Highlights regions where T > T_initial (shock-heated)
   - Color encodes shock strength: T/T_initial

Physics:
- Sound speed c_s ~ 1 km/s, initial T ~ 50 K
- Shock heating: T_shock / T_0 = (P_shock / P_0)^((γ-1)/γ) for adiabatic shock
- Strong shock limit: T_shock / T_0 ≈ (γ-1)/(γ+1) * (M^2) for Mach M

Units: IMBH_ENCOUNTER ([L]=pc, [M]=1000 M_sun, [V]=km/s, [T]=0.978 Myr)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.colors import LogNorm, Normalize, LinearSegmentedColormap
from matplotlib.gridspec import GridSpec
import argparse
from pathlib import Path
import glob
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)

# =============================================================================
# DARK THEME COLORS
# =============================================================================

DARK_BG = '#0d0d12'
DARK_PANEL = '#15151f'
GRID_COLOR = '#3a3a4a'
TEXT_COLOR = '#f0f0f0'
ACCENT_COLOR = '#00ffff'

# Line colors for 1D plots - HIGH CONTRAST, COLORBLIND-FRIENDLY
COLOR_DENSITY = '#00ff00'    # Pure bright green
COLOR_PRESSURE = '#ff8800'   # Bright orange  
COLOR_TEMP = '#ff00ff'       # Magenta (distinct from others)
COLOR_MACH = '#00ccff'       # Cyan
COLOR_INITIAL = '#888888'    # Gray (for reference)

# =============================================================================
# PHYSICAL CONSTANTS & UNITS
# =============================================================================

TIME_UNIT = 0.978  # Myr
VELOCITY_UNIT = 1.0  # km/s
LENGTH_UNIT = 1.0  # pc
MASS_UNIT = 1000.0  # M_sun

# No display scaling - use simulation coordinates directly
# Simulation cloud: R=0.56 pc, M=1000 M☉, n=25000 cm⁻³, T~50 K
DISPLAY_SCALE_FACTOR = 1.0

GAMMA = 5.0 / 3.0
MU = 2.3
M_H = 1.67e-24
K_B = 1.38e-16
VELOCITY_TO_CGS = 1e5
TEMP_CONVERSION = MU * M_H / K_B * VELOCITY_TO_CGS**2

# =============================================================================
# COLORMAPS
# =============================================================================

def create_shock_colormaps():
    """Create colormaps optimized for shock visualization."""
    
    def make_cmap(colors, name):
        return LinearSegmentedColormap.from_list(name, colors, N=256)
    
    # SHOCK TEMPERATURE: Blue (cold/unshocked) -> White (initial) -> Yellow -> Red (hot/shocked)
    # Centered at T/T0 = 1 (white)
    shock_temp_colors = [
        (0.0, '#0044ff'),    # Very cold (T << T0) - blue
        (0.2, '#0088ff'),    # Cold - light blue
        (0.4, '#00cccc'),    # Cool - cyan
        (0.5, '#ffffff'),    # T = T0 (initial) - WHITE
        (0.6, '#ffff00'),    # Slightly heated - yellow
        (0.75, '#ff8800'),   # Moderately shocked - orange
        (0.9, '#ff0000'),    # Strong shock - red
        (1.0, '#ff00ff'),    # Extreme shock - magenta
    ]
    
    # DENSITY: Dark blue -> Cyan -> Green -> Yellow -> Red
    density_colors = [
        (0.0, '#000066'),
        (0.2, '#0066cc'),
        (0.4, '#00cc99'),
        (0.6, '#66ff33'),
        (0.8, '#ffcc00'),
        (1.0, '#ff3300'),
    ]
    
    # PRESSURE: Purple -> Blue -> Cyan -> Green -> Yellow -> Red
    pressure_colors = [
        (0.0, '#220066'),
        (0.2, '#0044cc'),
        (0.4, '#00aacc'),
        (0.6, '#44cc44'),
        (0.8, '#ffcc00'),
        (1.0, '#ff2200'),
    ]
    
    return {
        'shock_temp': make_cmap(shock_temp_colors, 'shock_temp'),
        'density': make_cmap(density_colors, 'density'),
        'pressure': make_cmap(pressure_colors, 'pressure'),
    }

CMAPS = create_shock_colormaps()

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
        
        return {
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
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return None


def compute_temperature(pres, dens):
    """Compute temperature from pressure and density."""
    specific_energy = pres / (dens + 1e-30)
    temperature = specific_energy * TEMP_CONVERSION / GAMMA
    return temperature


def compute_internal_mach(vel_x, vel_y, vel_z, sound, mass):
    """Compute internal Mach number relative to CoM."""
    total_mass = np.sum(mass)
    v_com_x = np.sum(mass * vel_x) / total_mass
    v_com_y = np.sum(mass * vel_y) / total_mass
    v_com_z = np.sum(mass * vel_z) / total_mass
    
    dvel_x = vel_x - v_com_x
    dvel_y = vel_y - v_com_y
    dvel_z = vel_z - v_com_z
    
    v_internal = np.sqrt(dvel_x**2 + dvel_y**2 + dvel_z**2)
    mach = v_internal / (sound + 1e-30)
    
    return mach, v_internal


def compute_com(pos_x, pos_y, pos_z, mass):
    """Compute center of mass."""
    total_mass = np.sum(mass)
    com_x = np.sum(mass * pos_x) / total_mass
    com_y = np.sum(mass * pos_y) / total_mass
    com_z = np.sum(mass * pos_z) / total_mass
    return com_x, com_y, com_z


# =============================================================================
# 1D PROFILE EXTRACTION
# =============================================================================

def extract_vertical_profile(snapshot, com, column_radius=0.1, n_bins=50):
    """
    Extract 1D vertical (Z) profile from a column of particles near the CoM.
    
    This shows the VERTICAL COMPRESSION SHOCK structure.
    
    Parameters:
    -----------
    snapshot : dict - particle data
    com : tuple - (x, y, z) center of mass
    column_radius : float - radius of column in X-Y plane (pc)
    n_bins : int - number of bins in Z direction
    
    Returns:
    --------
    z_bins : array - bin centers
    profiles : dict - binned profiles of density, pressure, temperature, mach
    """
    com_x, com_y, com_z = com
    
    # Select particles in a vertical column around CoM
    dx = snapshot['pos_x'] - com_x
    dy = snapshot['pos_y'] - com_y
    r_xy = np.sqrt(dx**2 + dy**2)
    
    in_column = r_xy < column_radius
    
    if np.sum(in_column) < 10:
        # Too few particles, expand column
        column_radius *= 2
        in_column = r_xy < column_radius
    
    # Get particles in column
    z_col = snapshot['pos_z'][in_column] - com_z
    dens_col = snapshot['dens'][in_column]
    pres_col = snapshot['pres'][in_column]
    mass_col = snapshot['mass'][in_column]
    sound_col = snapshot['sound'][in_column]
    
    temp_col = compute_temperature(pres_col, dens_col)
    
    # Velocity relative to CoM for Mach calculation
    total_mass = np.sum(snapshot['mass'])
    v_com_z = np.sum(snapshot['mass'] * snapshot['vel_z']) / total_mass
    vel_z_rel = snapshot['vel_z'][in_column] - v_com_z
    mach_z = np.abs(vel_z_rel) / (sound_col + 1e-30)
    
    # Bin the data
    z_range = np.max(np.abs(z_col)) * 1.1
    z_edges = np.linspace(-z_range, z_range, n_bins + 1)
    z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])
    
    # Mass-weighted binning
    profiles = {
        'z': z_centers,
        'density': np.zeros(n_bins),
        'pressure': np.zeros(n_bins),
        'temperature': np.zeros(n_bins),
        'mach_z': np.zeros(n_bins),
        'n_particles': np.zeros(n_bins),
    }
    
    bin_indices = np.digitize(z_col, z_edges) - 1
    bin_indices = np.clip(bin_indices, 0, n_bins - 1)
    
    for i in range(n_bins):
        mask = bin_indices == i
        if np.sum(mask) > 0:
            mass_bin = mass_col[mask]
            total_bin_mass = np.sum(mass_bin)
            
            profiles['density'][i] = np.sum(mass_bin * dens_col[mask]) / total_bin_mass
            profiles['pressure'][i] = np.sum(mass_bin * pres_col[mask]) / total_bin_mass
            profiles['temperature'][i] = np.sum(mass_bin * temp_col[mask]) / total_bin_mass
            profiles['mach_z'][i] = np.sum(mass_bin * mach_z[mask]) / total_bin_mass
            profiles['n_particles'][i] = np.sum(mask)
    
    profiles['column_radius'] = column_radius
    profiles['n_in_column'] = np.sum(in_column)
    
    return profiles


def extract_horizontal_profile(snapshot, com, slice_thickness=0.1, n_bins=80):
    """
    Extract 1D horizontal (X) profile from a slice through the midplane.
    
    This shows the IN-PLANE TIDAL STRETCHING and potential shocks.
    
    Parameters:
    -----------
    snapshot : dict - particle data
    com : tuple - (x, y, z) center of mass
    slice_thickness : float - thickness of slice in Z and Y directions (pc)
    n_bins : int - number of bins in X direction
    
    Returns:
    --------
    profiles : dict - binned profiles along X
    """
    com_x, com_y, com_z = com
    
    # Select particles in a horizontal slice along X-axis
    dz = np.abs(snapshot['pos_z'] - com_z)
    dy = np.abs(snapshot['pos_y'] - com_y)
    
    in_slice = (dz < slice_thickness) & (dy < slice_thickness)
    
    if np.sum(in_slice) < 10:
        slice_thickness *= 2
        in_slice = (dz < slice_thickness) & (dy < slice_thickness)
    
    # Get particles in slice
    x_slice = snapshot['pos_x'][in_slice] - com_x
    dens_slice = snapshot['dens'][in_slice]
    pres_slice = snapshot['pres'][in_slice]
    mass_slice = snapshot['mass'][in_slice]
    sound_slice = snapshot['sound'][in_slice]
    
    temp_slice = compute_temperature(pres_slice, dens_slice)
    
    # Velocity for Mach
    total_mass = np.sum(snapshot['mass'])
    v_com_x = np.sum(snapshot['mass'] * snapshot['vel_x']) / total_mass
    vel_x_rel = snapshot['vel_x'][in_slice] - v_com_x
    mach_x = np.abs(vel_x_rel) / (sound_slice + 1e-30)
    
    # Bin the data
    x_range = np.max(np.abs(x_slice)) * 1.1
    if x_range < 0.1:
        x_range = 1.0  # Minimum range
    x_edges = np.linspace(-x_range, x_range, n_bins + 1)
    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
    
    profiles = {
        'x': x_centers,
        'density': np.zeros(n_bins),
        'pressure': np.zeros(n_bins),
        'temperature': np.zeros(n_bins),
        'mach_x': np.zeros(n_bins),
        'n_particles': np.zeros(n_bins),
    }
    
    bin_indices = np.digitize(x_slice, x_edges) - 1
    bin_indices = np.clip(bin_indices, 0, n_bins - 1)
    
    for i in range(n_bins):
        mask = bin_indices == i
        if np.sum(mask) > 0:
            mass_bin = mass_slice[mask]
            total_bin_mass = np.sum(mass_bin)
            
            profiles['density'][i] = np.sum(mass_bin * dens_slice[mask]) / total_bin_mass
            profiles['pressure'][i] = np.sum(mass_bin * pres_slice[mask]) / total_bin_mass
            profiles['temperature'][i] = np.sum(mass_bin * temp_slice[mask]) / total_bin_mass
            profiles['mach_x'][i] = np.sum(mass_bin * mach_x[mask]) / total_bin_mass
            profiles['n_particles'][i] = np.sum(mask)
    
    profiles['slice_thickness'] = slice_thickness
    profiles['n_in_slice'] = np.sum(in_slice)
    
    return profiles


# =============================================================================
# ANIMATION
# =============================================================================

def create_shock_diagnostic_animation(results_dir, output_file, fps=15):
    """
    Create animation with shock diagnostic panels:
    
    Layout (3 rows x 3 cols):
    ┌─────────────────────┬─────────────────────┬─────────────────────┐
    │   DENSITY (XY)      │   DENSITY (XZ)      │  SHOCK TEMP (XY)    │
    │   with column box   │   with slice box    │  T/T_initial        │
    ├─────────────────────┼─────────────────────┼─────────────────────┤
    │   VERTICAL SHOCK    │   IN-PLANE SHOCK    │   PHASE INFO        │
    │   ρ,P,T vs Z        │   ρ,P,T vs X        │   + Legend          │
    └─────────────────────┴─────────────────────┴─────────────────────┘
    """
    
    snapshot_files = sorted(glob.glob(f"{results_dir}/snapshot_*.csv"))
    if not snapshot_files:
        print(f"Error: No snapshots found in {results_dir}")
        return False
    
    print(f"Found {len(snapshot_files)} snapshots")
    
    # Load first snapshot for initial values
    first = load_snapshot(snapshot_files[0])
    if first is None:
        return False
    
    # Store initial values for normalization
    T_initial = np.median(compute_temperature(first['pres'], first['dens']))
    rho_initial = np.median(first['dens'])
    P_initial = np.median(first['pres'])
    
    # Calculate actual cloud parameters from first snapshot for display
    total_mass = first['mass'].sum()  # code units (1 code_mass = 1000 M_sun)
    com_init = compute_com(first['pos_x'], first['pos_y'], first['pos_z'], first['mass'])
    r_from_com = np.sqrt((first['pos_x'] - com_init[0])**2 + 
                         (first['pos_y'] - com_init[1])**2 + 
                         (first['pos_z'] - com_init[2])**2)
    r_cloud_rms = np.sqrt(np.average(r_from_com**2, weights=first['mass']))
    
    # Physical parameters (measured from simulation)
    M_CLOUD_MSUN = total_mass * 1000  # Convert code units to M_sun
    R_CLOUD_PC = r_cloud_rms          # Already in pc
    T_CLOUD_K = T_initial             # Already in K
    
    # Fixed parameters from config (IMBH setup)
    M_BH_MSUN = 1e5      # Black hole mass [M_sun] - from config
    V_INIT = 10.0        # Initial velocity [km/s] - from config
    B_PARAM = 3.0        # Impact parameter [pc] - from config  
    X_INIT = com_init[0] # Initial X position [pc] - measured
    Y_INIT = com_init[1] # Initial Y position [pc] - measured
    
    # Display scale factor (1.0 = no scaling, use simulation coordinates)
    s = DISPLAY_SCALE_FACTOR
    
    print(f"Initial values: T0={T_initial:.1e} K, ρ0={rho_initial:.2e}, P0={P_initial:.2e}")
    print(f"Cloud parameters: M={M_CLOUD_MSUN:.0f} M☉, R_rms={R_CLOUD_PC:.3f} pc, T0={T_CLOUD_K:.1f} K")
    
    # =========================================================================
    # FIXED AXIS RANGES - ALL PLOTS USE FIXED RANGES FOR CLARITY
    # =========================================================================
    
    # 1D profile Y-axis ranges (normalized values: quantity / initial_value)
    Y_RANGE_DENSITY = (0, 5.0)      # ρ/ρ₀: expect compression up to ~5x
    Y_RANGE_PRESSURE = (0, 10.0)    # P/P₀: pressure can spike more due to P ∝ ρ^γ
    Y_RANGE_TEMPERATURE = (0, 5.0)  # T/T₀: temperature from shock heating
    
    # 1D profile X-axis ranges (position relative to CoM, in pc)
    Z_RANGE_FIXED = (-0.4, 0.4)     # Z - Z_CoM range for vertical column
    X_RANGE_FIXED = (-1.5, 1.5)     # X - X_CoM range for horizontal slice
    
    # Fixed viewport in CoM frame: cloud stays centered, IMBH moves through frame
    # Symmetric ranges centered on cloud CoM (which is at origin in this frame)
    XLIM_FIXED = (-3.0, 3.0)        # X range in CoM frame
    YLIM_FIXED = (-3.0, 3.0)        # Y range in CoM frame
    ZLIM_FIXED = (-1.0, 1.0)        # Z range in CoM frame (smaller, vertical)
    
    print(f"Fixed viewport (CoM frame): X=[{XLIM_FIXED[0]}, {XLIM_FIXED[1]}], Y=[{YLIM_FIXED[0]}, {YLIM_FIXED[1]}], Z=[{ZLIM_FIXED[0]}, {ZLIM_FIXED[1]}]")
    
    # =========================================================================
    # CREATE FIGURE
    # =========================================================================
    
    fig = plt.figure(figsize=(20, 14), facecolor=DARK_BG)
    
    # GridSpec layout
    gs = GridSpec(2, 3, figure=fig, height_ratios=[1.2, 1],
                  hspace=0.25, wspace=0.25,
                  left=0.05, right=0.95, top=0.92, bottom=0.08)
    
    # Top row: 2D scatter plots
    ax_dens_xy = fig.add_subplot(gs[0, 0])
    ax_dens_xz = fig.add_subplot(gs[0, 1])
    ax_shock_temp = fig.add_subplot(gs[0, 2])
    
    # Bottom row: 1D profiles
    ax_vertical = fig.add_subplot(gs[1, 0])
    ax_horizontal = fig.add_subplot(gs[1, 1])
    ax_info = fig.add_subplot(gs[1, 2])
    
    for ax in [ax_dens_xy, ax_dens_xz, ax_shock_temp, ax_vertical, ax_horizontal, ax_info]:
        ax.set_facecolor(DARK_PANEL)
    
    # =========================================================================
    # TOP ROW: 2D SCATTER PLOTS
    # =========================================================================
    
    # Density XY view
    scatter_dens_xy = ax_dens_xy.scatter([], [], s=6, c=[], cmap=CMAPS['density'],
                                          norm=LogNorm(vmin=1e-2, vmax=1e2),
                                          alpha=0.85, edgecolors='none')
    cbar1 = fig.colorbar(scatter_dens_xy, ax=ax_dens_xy, shrink=0.8, pad=0.02)
    cbar1.set_label('Density', fontsize=10, color=TEXT_COLOR)
    cbar1.ax.tick_params(colors=TEXT_COLOR, labelsize=8)
    
    # Box showing vertical column selection
    column_box_xy, = ax_dens_xy.plot([], [], 'w-', linewidth=2, alpha=0.8)
    
    ax_dens_xy.set_xlim(XLIM_FIXED)
    ax_dens_xy.set_ylim(YLIM_FIXED)
    ax_dens_xy.set_xlabel('X - X_CoM (pc)', fontsize=10, color=TEXT_COLOR)
    ax_dens_xy.set_ylabel('Y - Y_CoM (pc)', fontsize=10, color=TEXT_COLOR)
    ax_dens_xy.set_title('Density (XY, CoM frame) + Vertical Column', fontsize=12, color=TEXT_COLOR, fontweight='bold')
    ax_dens_xy.tick_params(colors=TEXT_COLOR, labelsize=8)
    ax_dens_xy.set_aspect('equal')
    for spine in ax_dens_xy.spines.values():
        spine.set_color(GRID_COLOR)
    
    # Density XZ view
    scatter_dens_xz = ax_dens_xz.scatter([], [], s=6, c=[], cmap=CMAPS['density'],
                                          norm=LogNorm(vmin=1e-2, vmax=1e2),
                                          alpha=0.85, edgecolors='none')
    cbar2 = fig.colorbar(scatter_dens_xz, ax=ax_dens_xz, shrink=0.8, pad=0.02)
    cbar2.set_label('Density', fontsize=10, color=TEXT_COLOR)
    cbar2.ax.tick_params(colors=TEXT_COLOR, labelsize=8)
    
    # Box showing horizontal slice selection
    slice_box_xz, = ax_dens_xz.plot([], [], 'c-', linewidth=2, alpha=0.8)
    
    ax_dens_xz.set_xlim(XLIM_FIXED)
    ax_dens_xz.set_ylim(ZLIM_FIXED)
    ax_dens_xz.set_xlabel('X - X_CoM (pc)', fontsize=10, color=TEXT_COLOR)
    ax_dens_xz.set_ylabel('Z - Z_CoM (pc)', fontsize=10, color=TEXT_COLOR)
    ax_dens_xz.set_title('Density (XZ, CoM frame) + Horizontal Slice', fontsize=12, color=TEXT_COLOR, fontweight='bold')
    ax_dens_xz.tick_params(colors=TEXT_COLOR, labelsize=8)
    ax_dens_xz.set_aspect('equal')
    for spine in ax_dens_xz.spines.values():
        spine.set_color(GRID_COLOR)
    
    # Shock Temperature (T / T_initial)
    # Using diverging colormap centered at T/T0 = 1
    scatter_shock = ax_shock_temp.scatter([], [], s=6, c=[], cmap=CMAPS['shock_temp'],
                                           norm=LogNorm(vmin=0.3, vmax=30),
                                           alpha=0.85, edgecolors='none')
    cbar3 = fig.colorbar(scatter_shock, ax=ax_shock_temp, shrink=0.8, pad=0.02)
    cbar3.set_label('T / T₀ (shock heating)', fontsize=10, color=TEXT_COLOR)
    cbar3.ax.tick_params(colors=TEXT_COLOR, labelsize=8)
    
    ax_shock_temp.set_xlim(XLIM_FIXED)
    ax_shock_temp.set_ylim(YLIM_FIXED)
    ax_shock_temp.set_xlabel('X - X_CoM (pc)', fontsize=10, color=TEXT_COLOR)
    ax_shock_temp.set_ylabel('Y - Y_CoM (pc)', fontsize=10, color=TEXT_COLOR)
    ax_shock_temp.set_title('Shock Temperature (CoM, T/T₀=1: White)', fontsize=12, color=TEXT_COLOR, fontweight='bold')
    ax_shock_temp.tick_params(colors=TEXT_COLOR, labelsize=8)
    ax_shock_temp.set_aspect('equal')
    for spine in ax_shock_temp.spines.values():
        spine.set_color(GRID_COLOR)
    
    # IMBH markers - will be updated dynamically since we're in CoM frame
    # IMBH position in CoM frame = -CoM (since IMBH is at origin in lab frame)
    imbh_glow_xy = []
    imbh_glow_xz = []
    imbh_glow_shock = []
    
    for ax, glow_list in [(ax_dens_xy, imbh_glow_xy), (ax_shock_temp, imbh_glow_shock)]:
        for s, alpha in [(300, 0.03), (240, 0.06), (160, 0.12), (80, 0.25)]:
            sc = ax.scatter([0], [0], s=s, c='#00ff88', alpha=alpha, edgecolors='none', zorder=100)
            glow_list.append(sc)
        sc = ax.scatter([0], [0], s=120, c='white', marker='*',
                       edgecolors='#00ff88', linewidths=2, zorder=101)
        glow_list.append(sc)
    
    for s, alpha in [(300, 0.03), (240, 0.06), (160, 0.12), (80, 0.25)]:
        sc = ax_dens_xz.scatter([0], [0], s=s, c='#00ff88', alpha=alpha, edgecolors='none', zorder=100)
        imbh_glow_xz.append(sc)
    sc = ax_dens_xz.scatter([0], [0], s=120, c='white', marker='*',
                       edgecolors='#00ff88', linewidths=2, zorder=101)
    imbh_glow_xz.append(sc)
    
    # =========================================================================
    # BOTTOM ROW: 1D PROFILES
    # =========================================================================
    
    # VERTICAL PROFILE (Z direction) - Shows vertical compression shock
    ax_vertical.set_title('VERTICAL COMPRESSION SHOCK', fontsize=12, color='#ff6666', fontweight='bold')
    ax_vertical.set_xlabel('Z - Z_CoM (pc)', fontsize=10, color=TEXT_COLOR)
    ax_vertical.tick_params(colors=TEXT_COLOR, labelsize=8)
    ax_vertical.grid(True, alpha=0.2, color=GRID_COLOR)
    for spine in ax_vertical.spines.values():
        spine.set_color(GRID_COLOR)
    
    # Create twin axes for different quantities
    ax_v_dens = ax_vertical
    ax_v_pres = ax_vertical.twinx()
    ax_v_temp = ax_vertical.twinx()
    
    # HIDE all Y-axes on vertical profile (left plot) - labels will be on right plot
    ax_v_dens.tick_params(left=False, labelleft=False)
    ax_v_dens.spines['left'].set_visible(False)
    ax_v_pres.spines['right'].set_visible(False)
    ax_v_temp.spines['right'].set_visible(False)
    ax_v_pres.tick_params(right=False, labelright=False)
    ax_v_temp.tick_params(right=False, labelright=False)
    
    line_v_dens, = ax_v_dens.plot([], [], '-', color=COLOR_DENSITY, linewidth=2.5, label='ρ/ρ₀')
    line_v_pres, = ax_v_pres.plot([], [], '-', color=COLOR_PRESSURE, linewidth=2.5, label='P/P₀')
    line_v_temp, = ax_v_temp.plot([], [], '-', color=COLOR_TEMP, linewidth=2.5, label='T/T₀ (actual)')
    
    # =========================================================================
    # ANALYTICAL ADIABATIC REFERENCES (dashed lines)
    # For polytropic/isentropic gas with γ=5/3 (Lane-Emden n=1.5):
    #   P/P₀ = (ρ/ρ₀)^γ           (γ=5/3)
    #   T/T₀ = (ρ/ρ₀)^(γ-1)       (γ-1=2/3)
    # Initial cloud follows this relation (hydrostatic equilibrium)
    # Deviation above = shock heating (entropy increase)
    # =========================================================================
    # Use bright, high-contrast colors for dashed reference lines
    COLOR_P_ADIABATIC = '#00ffff'  # Bright cyan - distinct from orange pressure
    COLOR_T_ADIABATIC = '#ffff00'  # Bright yellow - distinct from pink temperature
    
    line_v_pres_adiabatic, = ax_v_pres.plot([], [], '--', color=COLOR_P_ADIABATIC, linewidth=2, 
                                             label='P_adiabatic', alpha=0.9)
    line_v_temp_adiabatic, = ax_v_temp.plot([], [], '--', color=COLOR_T_ADIABATIC, linewidth=2, 
                                             label='T_adiabatic', alpha=0.9)
    
    # Initial reference lines
    ax_v_dens.axhline(y=1, color=COLOR_DENSITY, linestyle=':', alpha=0.4, linewidth=1)
    ax_v_pres.axhline(y=1, color=COLOR_PRESSURE, linestyle=':', alpha=0.4, linewidth=1)
    ax_v_temp.axhline(y=1, color=COLOR_TEMP, linestyle=':', alpha=0.4, linewidth=1)
    
    # No Y-axis labels on left plot - they will be on the right plot
    
    # HORIZONTAL PROFILE (X direction) - Shows in-plane tidal stretching
    ax_horizontal.set_title('IN-PLANE TIDAL STRETCHING', fontsize=12, color='#66ff66', fontweight='bold')
    ax_horizontal.set_xlabel('X - X_CoM (pc)', fontsize=10, color=TEXT_COLOR)
    ax_horizontal.tick_params(colors=TEXT_COLOR, labelsize=8)
    ax_horizontal.grid(True, alpha=0.2, color=GRID_COLOR)
    for spine in ax_horizontal.spines.values():
        spine.set_color(GRID_COLOR)
    
    ax_h_dens = ax_horizontal
    ax_h_pres = ax_horizontal.twinx()
    ax_h_temp = ax_horizontal.twinx()
    
    # SHOW Y-axes on horizontal profile (right plot) with labels
    # Move pressure and temperature axes outward to avoid overlap
    ax_h_pres.spines['right'].set_position(('outward', 0))
    ax_h_temp.spines['right'].set_position(('outward', 70))
    
    line_h_dens, = ax_h_dens.plot([], [], '-', color=COLOR_DENSITY, linewidth=2.5, label='ρ/ρ₀')
    line_h_pres, = ax_h_pres.plot([], [], '-', color=COLOR_PRESSURE, linewidth=2.5, label='P/P₀')
    line_h_temp, = ax_h_temp.plot([], [], '-', color=COLOR_TEMP, linewidth=2.5, label='T/T₀ (actual)')
    
    # ANALYTICAL ADIABATIC REFERENCES for horizontal profile
    line_h_pres_adiabatic, = ax_h_pres.plot([], [], '--', color=COLOR_P_ADIABATIC, linewidth=2,
                                             label='P_adiabatic', alpha=0.9)
    line_h_temp_adiabatic, = ax_h_temp.plot([], [], '--', color=COLOR_T_ADIABATIC, linewidth=2,
                                             label='T_adiabatic', alpha=0.9)
    
    ax_h_dens.axhline(y=1, color=COLOR_DENSITY, linestyle=':', alpha=0.4, linewidth=1)
    ax_h_pres.axhline(y=1, color=COLOR_PRESSURE, linestyle=':', alpha=0.4, linewidth=1)
    ax_h_temp.axhline(y=1, color=COLOR_TEMP, linestyle=':', alpha=0.4, linewidth=1)
    
    # Y-axis labels with high contrast and larger font - ALL ON RIGHT PLOT
    ax_h_dens.set_ylabel('ρ/ρ₀', fontsize=12, color=COLOR_DENSITY, fontweight='bold')
    ax_h_pres.set_ylabel('P/P₀', fontsize=12, color=COLOR_PRESSURE, fontweight='bold')
    ax_h_temp.set_ylabel('T/T₀', fontsize=12, color=COLOR_TEMP, fontweight='bold')
    
    ax_h_dens.tick_params(axis='y', colors=COLOR_DENSITY, labelsize=9)
    ax_h_pres.tick_params(axis='y', colors=COLOR_PRESSURE, labelsize=9)
    ax_h_temp.tick_params(axis='y', colors=COLOR_TEMP, labelsize=9)
    
    ax_h_pres.spines['right'].set_color(COLOR_PRESSURE)
    ax_h_temp.spines['right'].set_color(COLOR_TEMP)
    
    # INFO PANEL
    ax_info.axis('off')
    
    # Title
    title_text = fig.suptitle('IMBH Tidal Disruption - Shock Diagnostics (CoM Frame)\nt = 0.000 Myr',
                              fontsize=16, color=TEXT_COLOR, fontweight='bold', y=0.97)
    
    # =========================================================================
    # SIMULATION PARAMETERS (from first snapshot - physics info panel)
    # =========================================================================
    # Format cloud mass for display (use scientific notation if >= 10000)
    if M_CLOUD_MSUN >= 10000:
        mass_str = f"{M_CLOUD_MSUN:.1e}"
    else:
        mass_str = f"{M_CLOUD_MSUN:.0f}"
    
    # Static simulation parameters - TOP of info panel
    params_text = f"""SIMULATION SETUP:
━━━━━━━━━━━━━━━━━━━━
IMBH:  M = {M_BH_MSUN:.0e} M☉
Cloud: M = {mass_str} M☉
       R = {R_CLOUD_PC:.2f} pc
       T₀= {T_CLOUD_K:.0f} K
       γ = 5/3 polytrope
Initial: x₀={X_INIT:.1f}pc
         y₀={Y_INIT:.1f}pc (b={B_PARAM}pc)
         v₀={V_INIT}km/s (+x)"""
    ax_info.text(0.05, 0.98, params_text, fontsize=7, color='#88ccff',
                 transform=ax_info.transAxes, family='monospace',
                 verticalalignment='top')
    
    # Dynamic info text - MIDDLE of panel (updated each frame)
    info_text = ax_info.text(0.05, 0.58, '', fontsize=8, color=TEXT_COLOR,
                             transform=ax_info.transAxes, family='monospace',
                             verticalalignment='top')
    
    # Analytical references - LOWER MIDDLE
    physics_text = """REFERENCES (dashed):
━━━━━━━━━━━━━━━━━━━━
P/P₀ = (ρ/ρ₀)^(5/3)
T/T₀ = (ρ/ρ₀)^(2/3)
solid > dashed → SHOCK"""
    ax_info.text(0.05, 0.32, physics_text, fontsize=7, color='#aaaaaa',
                 transform=ax_info.transAxes, family='monospace',
                 verticalalignment='top')
    
    # Legend - BOTTOM of panel
    legend_items = [
        (COLOR_DENSITY, '━', 'ρ/ρ₀'),
        (COLOR_PRESSURE, '━', 'P/P₀'),
        (COLOR_P_ADIABATIC, '--', 'P_ad'),
        (COLOR_TEMP, '━', 'T/T₀'),
        (COLOR_T_ADIABATIC, '--', 'T_ad'),
    ]
    
    for i, (color, style, label) in enumerate(legend_items):
        y_pos = 0.12 - i * 0.025
        ax_info.text(0.05, y_pos, style, fontsize=9, color=color,
                     transform=ax_info.transAxes, fontweight='bold')
        ax_info.text(0.14, y_pos, label, fontsize=8, color=TEXT_COLOR,
                     transform=ax_info.transAxes)
    
    # =========================================================================
    # ANIMATION UPDATE
    # =========================================================================
    
    def update(frame_idx):
        snapshot = load_snapshot(snapshot_files[frame_idx])
        if snapshot is None:
            return []
        
        time_code = frame_idx * 0.025
        time_myr = time_code * TIME_UNIT
        
        # Compute center of mass
        com = compute_com(snapshot['pos_x'], snapshot['pos_y'], snapshot['pos_z'], snapshot['mass'])
        com_x, com_y, com_z = com
        
        # =====================================================================
        # TRANSFORM TO CoM FRAME: cloud stays centered, IMBH moves
        # =====================================================================
        pos_x_com = snapshot['pos_x'] - com_x  # X relative to CoM
        pos_y_com = snapshot['pos_y'] - com_y  # Y relative to CoM
        pos_z_com = snapshot['pos_z'] - com_z  # Z relative to CoM
        
        # IMBH position in CoM frame (IMBH is at origin in lab frame)
        imbh_x_com = -com_x
        imbh_y_com = -com_y
        imbh_z_com = -com_z

        # Compute temperature and shock ratio
        temperature = compute_temperature(snapshot['pres'], snapshot['dens'])
        T_ratio = temperature / T_initial
        
        # Update 2D scatter plots - NOW IN CoM FRAME!
        scatter_dens_xy.set_offsets(np.column_stack([pos_x_com, pos_y_com]))
        scatter_dens_xy.set_array(snapshot['dens'])
        
        scatter_dens_xz.set_offsets(np.column_stack([pos_x_com, pos_z_com]))
        scatter_dens_xz.set_array(snapshot['dens'])
        
        scatter_shock.set_offsets(np.column_stack([pos_x_com, pos_y_com]))
        scatter_shock.set_array(T_ratio)
        
        # Update IMBH marker positions (they move in CoM frame)
        for sc in imbh_glow_xy:
            sc.set_offsets([[imbh_x_com, imbh_y_com]])
        for sc in imbh_glow_shock:
            sc.set_offsets([[imbh_x_com, imbh_y_com]])
        for sc in imbh_glow_xz:
            sc.set_offsets([[imbh_x_com, imbh_z_com]])
        
        # Extract 1D profiles (these already use CoM-relative coordinates internally)
        v_prof = extract_vertical_profile(snapshot, com, column_radius=0.15)
        h_prof = extract_horizontal_profile(snapshot, com, slice_thickness=0.15)
        
        # Draw selection boxes on 2D plots - NOW IN CoM FRAME (centered at 0,0)
        # Vertical column box (circle in XY view, centered at origin)
        theta = np.linspace(0, 2*np.pi, 50)
        col_r = v_prof['column_radius']
        column_box_xy.set_data(col_r * np.cos(theta), col_r * np.sin(theta))
        
        # Horizontal slice box (rectangle in XZ view, centered at origin in CoM frame)
        slice_t = h_prof['slice_thickness']
        x_extent_box = np.max(np.abs(h_prof['x'])) if len(h_prof['x']) > 0 else 1.0
        box_x = [-x_extent_box, x_extent_box, x_extent_box, -x_extent_box, -x_extent_box]
        box_z = [-slice_t, -slice_t, slice_t, slice_t, -slice_t]
        slice_box_xz.set_data(box_x, box_z)
        
        # Update vertical profile
        z_data = v_prof['z']
        valid_v = v_prof['n_particles'] > 0
        
        if np.any(valid_v):
            rho_v = v_prof['density'][valid_v] / rho_initial
            P_v = v_prof['pressure'][valid_v] / P_initial
            T_v = v_prof['temperature'][valid_v] / T_initial
            z_v = z_data[valid_v]
            
            line_v_dens.set_data(z_v, rho_v)
            line_v_pres.set_data(z_v, P_v)
            line_v_temp.set_data(z_v, T_v)
            
            # ================================================================
            # ANALYTICAL ADIABATIC REFERENCES:
            #   P_ad/P₀ = (ρ/ρ₀)^γ        (γ=5/3)
            #   T_ad/T₀ = (ρ/ρ₀)^(γ-1)    (γ-1=2/3)
            # If actual > adiabatic → SHOCK (entropy increase!)
            # ================================================================
            P_adiabatic_v = rho_v ** GAMMA           # = rho^(5/3) for gamma=5/3
            T_adiabatic_v = rho_v ** (GAMMA - 1)     # = rho^(2/3) for gamma=5/3
            line_v_pres_adiabatic.set_data(z_v, P_adiabatic_v)
            line_v_temp_adiabatic.set_data(z_v, T_adiabatic_v)
            
            # FIXED X-axis range for 1D profiles
            ax_v_dens.set_xlim(Z_RANGE_FIXED)
            
            # FIXED Y-AXIS RANGES - don't change during animation!
            ax_v_dens.set_ylim(Y_RANGE_DENSITY)
            ax_v_pres.set_ylim(Y_RANGE_PRESSURE)
            ax_v_temp.set_ylim(Y_RANGE_TEMPERATURE)
        
        # Update horizontal profile
        x_data = h_prof['x']
        valid_h = h_prof['n_particles'] > 0
        
        if np.any(valid_h):
            rho_h = h_prof['density'][valid_h] / rho_initial
            P_h = h_prof['pressure'][valid_h] / P_initial
            T_h = h_prof['temperature'][valid_h] / T_initial
            x_h = x_data[valid_h]
            
            line_h_dens.set_data(x_h, rho_h)
            line_h_pres.set_data(x_h, P_h)
            line_h_temp.set_data(x_h, T_h)
            
            # ANALYTICAL ADIABATIC REFERENCES for horizontal profile
            P_adiabatic_h = rho_h ** GAMMA           # = rho^(5/3)
            T_adiabatic_h = rho_h ** (GAMMA - 1)     # = rho^(2/3)
            line_h_pres_adiabatic.set_data(x_h, P_adiabatic_h)
            line_h_temp_adiabatic.set_data(x_h, T_adiabatic_h)
            
            # FIXED X-axis range
            ax_h_dens.set_xlim(X_RANGE_FIXED)
            
            # FIXED Y-AXIS RANGES - don't change during animation!
            ax_h_dens.set_ylim(Y_RANGE_DENSITY)
            ax_h_pres.set_ylim(Y_RANGE_PRESSURE)
            ax_h_temp.set_ylim(Y_RANGE_TEMPERATURE)
        
        # Update title
        title_text.set_text(f'IMBH Tidal Disruption - Shock Diagnostics\nt = {time_myr:.3f} Myr')
        
        # Update info
        T_max_ratio = np.max(T_ratio)
        rho_max = np.max(snapshot['dens']) / rho_initial
        P_max = np.max(snapshot['pres']) / P_initial
        
        # Distance from IMBH to cloud CoM
        dist_imbh = np.sqrt(com_x**2 + com_y**2 + com_z**2)
        
        # Compact info string - removed SHAPE, focused on dynamics
        info_str = f"""t = {time_myr:.2f} Myr
━━━━━━━━━━━━━━━━━━━━
CoM: ({com_x:.1f}, {com_y:.1f}, {com_z:.1f}) pc
IMBH→CoM: {dist_imbh:.2f} pc

EXTREMA:
  T_max/T₀ = {T_max_ratio:.1f}
  ρ_max/ρ₀ = {rho_max:.1f}
  P_max/P₀ = {P_max:.1f}"""
        
        info_text.set_text(info_str)
        
        if frame_idx % 20 == 0:
            print(f"  Frame {frame_idx+1}/{len(snapshot_files)}: t={time_myr:.3f} Myr, T_max/T0={T_max_ratio:.1f}")
        
        return []
    
    # =========================================================================
    # CREATE ANIMATION
    # =========================================================================
    
    print("Creating shock diagnostics animation...")
    print(f"  Frames: {len(snapshot_files)}")
    print(f"  FPS: {fps}")
    print(f"  Output: {output_file}")
    
    anim = FuncAnimation(fig, update, frames=len(snapshot_files),
                        interval=1000/fps, blit=False)
    
    print("Saving animation...")
    writer = PillowWriter(fps=fps)
    anim.save(output_file, writer=writer, dpi=100)
    
    print(f"✓ Animation saved: {output_file}")
    plt.close(fig)
    
    return True


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description='IMBH Tidal Disruption - Shock Diagnostics Animation')
    parser.add_argument('results_dir', type=str, help='Directory containing snapshot CSV files')
    parser.add_argument('--output', '-o', type=str, default='shock_diagnostics.gif',
                        help='Output animation filename')
    parser.add_argument('--fps', type=int, default=15, help='Frames per second')
    
    args = parser.parse_args()
    
    print("=" * 65)
    print("IMBH Tidal Disruption - Shock Diagnostics")
    print("=" * 65)
    print(f"Results: {args.results_dir}")
    print(f"Output:  {args.output}")
    print("=" * 65)
    
    success = create_shock_diagnostic_animation(args.results_dir, args.output, args.fps)
    
    if success:
        print("\n" + "=" * 65)
        print("✓ Shock diagnostics animation complete!")
        print("=" * 65)
    else:
        print("\n✗ Animation failed!")
        return 1
    
    return 0


if __name__ == '__main__':
    exit(main())
