#!/usr/bin/env python3
"""
IMBH Tidal Disruption Animation - High Contrast Dark Theme
===========================================================

Creates stunning, publication-quality animations with:
- HIGH CONTRAST colormaps for easy reading
- Dark background (ParaView/VTK style)
- Energy evolution with analytical predictions
- Modern visual aesthetics

Units: IMBH_ENCOUNTER ([L]=pc, [M]=1000 M_sun, [V]=km/s, [T]=0.978 Myr)

Energy Evolution Physics:
- Kinetic Energy: Increases during IMBH flyby due to tidal acceleration
- Thermal Energy: Increases from shock heating and PdV work  
- Potential Energy: Becomes more negative as cloud approaches IMBH
- Total Energy: Should be conserved (check for numerical errors)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.colors import LogNorm, Normalize, LinearSegmentedColormap
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patheffects as path_effects
import argparse
from pathlib import Path
import glob
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)

# =============================================================================
# HIGH CONTRAST DARK THEME - Optimized for Readability
# =============================================================================

# Dark background colors
DARK_BG = '#0d0d12'
DARK_PANEL = '#15151f'
GRID_COLOR = '#3a3a4a'
TEXT_COLOR = '#f0f0f0'
ACCENT_COLOR = '#00ffff'  # Bright cyan
HIGHLIGHT_COLOR = '#ff4444'  # Bright red

# High-contrast color scheme for lines
COLORS = {
    'kinetic': '#00ff00',      # Bright green
    'thermal': '#ff6600',      # Bright orange  
    'potential': '#ff00ff',    # Bright magenta
    'total': '#ffffff',        # White
    'analytic': '#ffff00',     # Yellow (for analytical)
}

# =============================================================================
# HIGH CONTRAST COLORMAPS - Optimized for visibility
# =============================================================================

def create_high_contrast_colormaps():
    """Create high-contrast colormaps for easy visualization."""
    
    # DENSITY: Black -> Blue -> Cyan -> Yellow -> White (high contrast)
    density_colors = [
        (0.0, '#000010'),    # Near black
        (0.1, '#000080'),    # Deep blue
        (0.25, '#0066ff'),   # Bright blue
        (0.4, '#00ffff'),    # Cyan
        (0.55, '#00ff66'),   # Green-cyan
        (0.7, '#ffff00'),    # Yellow
        (0.85, '#ff8800'),   # Orange
        (1.0, '#ffffff'),    # White
    ]
    
    # TEMPERATURE: BRIGHT visible colors - Blue -> Cyan -> Green -> Yellow -> Orange -> Red -> White
    # Starts with visible blue instead of black!
    temp_colors = [
        (0.0, '#0044ff'),    # Bright blue (COLD - visible!)
        (0.12, '#0088ff'),   # Light blue
        (0.25, '#00cccc'),   # Cyan
        (0.38, '#00ff88'),   # Green-cyan
        (0.5, '#88ff00'),    # Yellow-green
        (0.62, '#ffff00'),   # Yellow
        (0.75, '#ff8800'),   # Orange
        (0.88, '#ff0000'),   # Red (HOT)
        (1.0, '#ffffff'),    # White (HOTTEST)
    ]
    
    # MACH NUMBER: Blue -> Cyan -> White -> Yellow -> Orange -> Red (diverging at M=1)
    mach_colors = [
        (0.0, '#0000aa'),    # Deep blue (M=0)
        (0.15, '#0066ff'),   # Blue
        (0.3, '#00ccff'),    # Light blue
        (0.45, '#66ffff'),   # Cyan
        (0.5, '#ffffff'),    # White (M=1 sonic)
        (0.55, '#ffff66'),   # Light yellow
        (0.7, '#ffcc00'),    # Yellow
        (0.85, '#ff6600'),   # Orange
        (1.0, '#ff0000'),    # Red (M=5 hypersonic)
    ]
    
    # SHOCK: Purple -> Magenta -> Red -> Orange -> Yellow (dramatic)
    shock_colors = [
        (0.0, '#100020'),    # Dark purple
        (0.2, '#440066'),    # Purple
        (0.4, '#880088'),    # Magenta
        (0.6, '#cc0044'),    # Red-magenta
        (0.8, '#ff6600'),    # Orange
        (1.0, '#ffff00'),    # Yellow
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
        'density': make_cmap(density_colors, 'density_hc'),
        'temperature': make_cmap(temp_colors, 'temp_hc'),
        'mach': make_cmap(mach_colors, 'mach_hc'),
        'shock': make_cmap(shock_colors, 'shock_hc'),
    }

CUSTOM_CMAPS = create_high_contrast_colormaps()

# =============================================================================
# PHYSICAL CONSTANTS
# =============================================================================

TIME_UNIT = 0.978  # Myr
VELOCITY_UNIT = 1.0  # km/s
LENGTH_UNIT = 1.0  # pc
MASS_UNIT = 1000.0  # M_sun

GAMMA = 5.0 / 3.0
MU = 2.3
M_H = 1.67e-24
K_B = 1.38e-16
VELOCITY_TO_CGS = 1e5
TEMP_CONVERSION = MU * M_H / K_B * VELOCITY_TO_CGS**2

# IMBH parameters for analytical estimates
M_IMBH = 1000.0  # M_sun in code units
R_PERI = 3.0  # pc (perihelion distance)
G_CODE = 1.0  # Gravitational constant in code units


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
        
        if 'ene' in data.dtype.names:
            snapshot['ene'] = data['ene']
        
        return snapshot
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return None


def load_energy_file(filepath):
    """Load energy evolution data."""
    try:
        data = np.genfromtxt(filepath, comments='#')
        if data.ndim == 1:
            data = data.reshape(1, -1)
        return {
            'time': data[:, 0],
            'kinetic': data[:, 1],
            'thermal': data[:, 2],
            'potential': data[:, 3],
            'total': data[:, 4],
        }
    except Exception as e:
        print(f"Error loading energy file {filepath}: {e}")
        return None


# =============================================================================
# ANALYTICAL ENERGY PREDICTIONS
# =============================================================================

def compute_analytical_energy(t, E_kin_0, E_therm_0, E_pot_0, E_tot_0,
                               v_rel=10.0, b=3.0, M_imbh=1000.0):
    """
    Compute analytical energy predictions for tidal encounter.
    
    During a hyperbolic encounter:
    - Kinetic energy: KE ~ KE_0 + delta_KE from tidal work
    - Thermal energy: Increases from shock heating (irreversible)
    - Potential energy: Varies with distance to IMBH
    - Total energy: Should be conserved
    
    Parameters:
    -----------
    t : array
        Time in Myr
    v_rel : float
        Relative velocity at infinity (km/s -> code units)
    b : float
        Impact parameter (pc)
    M_imbh : float
        IMBH mass (M_sun)
    """
    t_code = t / TIME_UNIT  # Convert to code units
    
    # Characteristic timescale: crossing time ~ b / v_rel
    t_cross = b / v_rel * TIME_UNIT  # in Myr
    
    # Closest approach time (assume t=0 is start, closest at t ~ some time)
    t_peri = 1.0  # Myr - approximate time of perihelion
    
    # Analytical predictions:
    
    # 1. Total energy should be CONSERVED
    E_tot_analytic = np.ones_like(t) * E_tot_0
    
    # 2. Potential energy varies with distance
    # During approach: becomes more negative
    # After perihelion: becomes less negative
    # Approximate as: E_pot ~ E_pot_0 * (1 + A * exp(-(t-t_peri)^2/sigma^2))
    sigma = 0.5  # Width of encounter in Myr
    approach_factor = 1.5  # How much deeper the potential gets
    E_pot_analytic = E_pot_0 * (1 + (approach_factor - 1) * np.exp(-((t - t_peri)**2) / (2 * sigma**2)))
    
    # 3. Kinetic energy: compensates potential + includes tidal heating
    # KE = E_tot - E_pot - E_therm
    # During encounter, kinetic energy increases due to tidal acceleration
    
    # 4. Thermal energy: increases due to shock heating
    # This is IRREVERSIBLE - thermal energy only increases
    # Model as: E_therm ~ E_therm_0 * (1 + heating_factor * (1 - exp(-t/t_heat)))
    t_heat = 0.8  # Heating timescale in Myr
    max_heating = 2.0  # Factor by which thermal energy can increase
    E_therm_analytic = E_therm_0 * (1 + (max_heating - 1) * (1 - np.exp(-t / t_heat)))
    
    # 5. Kinetic from conservation (if total is conserved)
    E_kin_analytic = E_tot_analytic - E_pot_analytic - E_therm_analytic
    
    return {
        'kinetic': E_kin_analytic,
        'thermal': E_therm_analytic,
        'potential': E_pot_analytic,
        'total': E_tot_analytic,
    }


# =============================================================================
# PHYSICS CALCULATIONS
# =============================================================================

def compute_temperature(pres, dens):
    """Compute temperature from pressure and density."""
    specific_energy = pres / (dens + 1e-30)
    temperature = specific_energy * TEMP_CONVERSION / GAMMA
    return temperature


def compute_mach_number(vel_x, vel_y, vel_z, sound, mass=None):
    """
    Compute INTERNAL Mach number relative to center-of-mass motion.
    
    The key insight: we want to measure internal turbulence/shocks, not bulk motion.
    
    For a cloud orbiting an IMBH at v_orb ~ 4 km/s with sound speed c_s ~ 1 km/s:
    - WRONG: M = |v_total| / c_s ~ 4 everywhere (all red)
    - RIGHT: M = |v - v_CoM| / c_s ~ 0-1 initially, increases during tidal disruption
    
    This correctly shows:
    - Initial hydrostatic cloud: low internal Mach (blue/green)
    - Tidal shocking: high internal Mach in shock fronts (orange/red)
    - Turbulent wake: moderate internal Mach (yellow/green)
    """
    # If mass not provided, use equal weighting
    if mass is None:
        mass = np.ones_like(vel_x)
    
    total_mass = np.sum(mass)
    
    # Compute mass-weighted center-of-mass velocity (bulk motion)
    v_com_x = np.sum(mass * vel_x) / total_mass
    v_com_y = np.sum(mass * vel_y) / total_mass
    v_com_z = np.sum(mass * vel_z) / total_mass
    
    # Velocity relative to CoM (internal motion only)
    dvel_x = vel_x - v_com_x
    dvel_y = vel_y - v_com_y
    dvel_z = vel_z - v_com_z
    
    # Internal velocity magnitude
    v_internal = np.sqrt(dvel_x**2 + dvel_y**2 + dvel_z**2)
    
    # Internal Mach number
    mach = v_internal / (sound + 1e-30)
    
    return mach, v_internal


def detect_shocks(dens, pres, rho_0=None, P_0=None, threshold=1.5):
    """Detect shock regions."""
    if rho_0 is None:
        rho_0 = np.median(dens[dens > 0])
    if P_0 is None:
        P_0 = np.median(pres[pres > 0])
    
    P_adiabatic = P_0 * (dens / rho_0)**GAMMA
    shock_param = pres / (P_adiabatic + 1e-30)
    shocked = shock_param > threshold
    
    return shocked, shock_param


def compute_tidal_deformation(pos_x, pos_y, pos_z, mass):
    """Compute tidal deformation metrics."""
    total_mass = np.sum(mass)
    
    com_x = np.sum(pos_x * mass) / total_mass
    com_y = np.sum(pos_y * mass) / total_mass
    com_z = np.sum(pos_z * mass) / total_mass
    
    dx = pos_x - com_x
    dy = pos_y - com_y
    dz = pos_z - com_z
    
    r_inplane = np.sqrt(dx**2 + dy**2)
    r_vertical = np.abs(dz)
    
    r_xy_rms = np.sqrt(np.sum(mass * r_inplane**2) / total_mass)
    r_z_rms = np.sqrt(np.sum(mass * r_vertical**2) / total_mass)
    
    deformation_ratio = r_xy_rms / (r_z_rms + 1e-10)
    
    return (com_x, com_y, com_z), r_xy_rms, r_z_rms, deformation_ratio


def style_dark_axis(ax, xlabel='', ylabel='', title='', is_3d=False):
    """Apply dark theme styling to an axis."""
    ax.set_facecolor(DARK_PANEL)
    
    if is_3d:
        ax.set_xlabel(xlabel, fontsize=10, color=TEXT_COLOR, labelpad=5)
        ax.set_ylabel(ylabel, fontsize=10, color=TEXT_COLOR, labelpad=5)
        ax.set_zlabel('Z (pc)', fontsize=10, color=TEXT_COLOR, labelpad=5)
        ax.tick_params(colors=TEXT_COLOR, labelsize=8)
        ax.grid(True, alpha=0.15, color=GRID_COLOR)
    else:
        ax.set_xlabel(xlabel, fontsize=11, color=TEXT_COLOR, fontweight='bold')
        ax.set_ylabel(ylabel, fontsize=11, color=TEXT_COLOR, fontweight='bold')
        ax.tick_params(colors=TEXT_COLOR, labelsize=9)
        for spine in ax.spines.values():
            spine.set_color(GRID_COLOR)
            spine.set_linewidth(1.5)
        ax.grid(True, alpha=0.25, color=GRID_COLOR, linestyle='-', linewidth=0.5)
    
    if title:
        ax.set_title(title, fontsize=13, fontweight='bold', color=TEXT_COLOR, pad=10)


# =============================================================================
# SINGLE MODE ANIMATION
# =============================================================================

def create_single_mode_animation(results_dir, output_file, mode='density', 
                                  xlim=None, ylim=None, fps=15):
    """Create single panel animation with high-contrast visualization."""
    
    snapshot_files = sorted(glob.glob(f"{results_dir}/snapshot_*.csv"))
    if not snapshot_files:
        print(f"Error: No snapshots found in {results_dir}")
        return False
    
    print(f"Found {len(snapshot_files)} snapshots")
    
    first = load_snapshot(snapshot_files[0])
    if first is None:
        return False
    
    if xlim is None:
        pad = 1.0
        xlim = (first['pos_x'].min() - pad, first['pos_x'].max() + pad)
    if ylim is None:
        pad = 1.0
        ylim = (first['pos_y'].min() - pad, first['pos_y'].max() + pad)
    
    fig, ax = plt.subplots(figsize=(14, 12), facecolor=DARK_BG)
    ax.set_facecolor(DARK_BG)
    
    if mode == 'density':
        cmap = CUSTOM_CMAPS['density']
        norm = LogNorm(vmin=1e-2, vmax=1e2)
        cbar_label = 'Log Density (code units)'
        title_prefix = 'Gas Density'
    elif mode == 'temperature':
        cmap = CUSTOM_CMAPS['temperature']
        norm = LogNorm(vmin=1e2, vmax=1e7)
        cbar_label = 'Temperature (K)'
        title_prefix = 'Temperature'
    elif mode == 'mach':
        cmap = CUSTOM_CMAPS['mach']
        norm = Normalize(vmin=0, vmax=5)
        cbar_label = 'Mach Number'
        title_prefix = 'Mach Number'
    elif mode == 'shock':
        cmap = CUSTOM_CMAPS['shock']
        norm = LogNorm(vmin=1, vmax=100)
        cbar_label = 'Shock Parameter (P/P_adiabatic)'
        title_prefix = 'Shock Structure'
    else:
        cmap = CUSTOM_CMAPS['density']
        norm = LogNorm(vmin=1e-2, vmax=1e2)
        cbar_label = 'Log Density'
        title_prefix = 'Gas'
    
    scatter = ax.scatter([], [], s=10, c=[], cmap=cmap, norm=norm, 
                        alpha=0.9, edgecolors='none')
    
    cbar = fig.colorbar(scatter, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label(cbar_label, fontsize=14, color=TEXT_COLOR, fontweight='bold')
    cbar.ax.tick_params(colors=TEXT_COLOR, labelsize=11)
    cbar.outline.set_edgecolor(GRID_COLOR)
    
    # IMBH with glow
    for s, alpha in [(300, 0.03), (250, 0.06), (200, 0.1), (150, 0.2)]:
        ax.scatter([0], [0], s=s, c=ACCENT_COLOR, alpha=alpha, edgecolors='none', zorder=100)
    ax.scatter([0], [0], s=120, c='white', marker='*', 
               edgecolors=ACCENT_COLOR, linewidths=2, zorder=101)
    
    time_text = ax.text(0.02, 0.98, '', transform=ax.transAxes, fontsize=16,
                       fontweight='bold', va='top', color=TEXT_COLOR,
                       path_effects=[path_effects.withStroke(linewidth=3, foreground=DARK_BG)])
    
    physics_text = ax.text(0.98, 0.02, '', transform=ax.transAxes, fontsize=11,
                          va='bottom', ha='right', color=TEXT_COLOR,
                          family='monospace', alpha=0.9,
                          path_effects=[path_effects.withStroke(linewidth=2, foreground=DARK_BG)])
    
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel('X (pc)', fontsize=14, color=TEXT_COLOR, fontweight='bold')
    ax.set_ylabel('Y (pc)', fontsize=14, color=TEXT_COLOR, fontweight='bold')
    ax.tick_params(colors=TEXT_COLOR, labelsize=11)
    ax.set_aspect('equal')
    for spine in ax.spines.values():
        spine.set_color(GRID_COLOR)
    ax.grid(True, alpha=0.15, color=GRID_COLOR)
    
    ax.set_title(f'{title_prefix} - IMBH Tidal Disruption', 
                fontsize=18, fontweight='bold', color=TEXT_COLOR, pad=15)
    
    def update(frame_idx):
        snapshot = load_snapshot(snapshot_files[frame_idx])
        if snapshot is None:
            return scatter,
        
        if mode == 'density':
            color_data = snapshot['dens']
        elif mode == 'temperature':
            color_data = compute_temperature(snapshot['pres'], snapshot['dens'])
        elif mode == 'mach':
            color_data, _ = compute_mach_number(
                snapshot['vel_x'], snapshot['vel_y'], snapshot['vel_z'], 
                snapshot['sound'], snapshot['mass'])
        elif mode == 'shock':
            _, color_data = detect_shocks(snapshot['dens'], snapshot['pres'])
        else:
            color_data = snapshot['dens']
        
        scatter.set_offsets(np.column_stack([snapshot['pos_x'], snapshot['pos_y']]))
        scatter.set_array(color_data)
        
        time_code = frame_idx * 0.025
        time_myr = time_code * TIME_UNIT
        
        mach, v_internal = compute_mach_number(
            snapshot['vel_x'], snapshot['vel_y'], snapshot['vel_z'], 
            snapshot['sound'], snapshot['mass'])
        temp = compute_temperature(snapshot['pres'], snapshot['dens'])
        
        time_text.set_text(f't = {time_myr:.3f} Myr')
        
        physics_info = (f'N = {len(snapshot["dens"]):,}\n'
                       f'ρ_max = {snapshot["dens"].max():.1e}\n'
                       f'T_max = {temp.max():.1e} K\n'
                       f'M_int_max = {mach.max():.1f}')
        physics_text.set_text(physics_info)
        
        if frame_idx % 20 == 0:
            print(f"  Frame {frame_idx+1}/{len(snapshot_files)}: t={time_myr:.3f} Myr")
        
        return scatter, time_text, physics_text
    
    print(f"Creating {mode} animation...")
    anim = FuncAnimation(fig, update, frames=len(snapshot_files), 
                        interval=1000//fps, blit=False)
    
    print(f"Saving to {output_file}...")
    writer = PillowWriter(fps=fps)
    anim.save(output_file, writer=writer, dpi=100)
    plt.close(fig)
    
    print(f"✓ Animation saved: {output_file}")
    return True


# =============================================================================
# MAIN TIDAL ANIMATION WITH ENERGY PANEL
# =============================================================================

def create_tidal_disruption_animation(results_dir, output_file, xlim=None, ylim=None, fps=15):
    """
    Create comprehensive 6-panel animation with ENERGY EVOLUTION panel.
    
    Layout:
    ┌─────────────────────┬─────────────────────┬─────────────────────┐
    │   DENSITY (XY)      │   TEMPERATURE       │   MACH NUMBER       │
    ├─────────────────────┼─────────────────────┼─────────────────────┤
    │   DENSITY (XZ)      │   ENERGY EVOLUTION  │   PHYSICS INFO      │
    └─────────────────────┴─────────────────────┴─────────────────────┘
    """
    
    snapshot_files = sorted(glob.glob(f"{results_dir}/snapshot_*.csv"))
    energy_file = f"{results_dir}/energy.dat"
    
    if not snapshot_files:
        print(f"Error: No snapshots found in {results_dir}")
        return False
    
    print(f"Found {len(snapshot_files)} snapshots")
    
    # Load energy data
    energy_data = load_energy_file(energy_file)
    if energy_data is not None:
        print(f"Loaded energy data: {len(energy_data['time'])} points")
        # Convert time from code units to Myr
        energy_time_myr = energy_data['time'] * TIME_UNIT
    else:
        print("Warning: No energy data found")
        energy_time_myr = None
    
    first = load_snapshot(snapshot_files[0])
    if first is None:
        return False
    
    if xlim is None:
        xlim = (-3, 3)  # CoM frame: cloud stays centered
    if ylim is None:
        ylim = (-3, 3)  # CoM frame: symmetric around origin
    
    # Compute analytical energy predictions if we have initial data
    analytic_energy = None
    if energy_data is not None and len(energy_data['time']) > 0:
        t_analytic = np.linspace(0, energy_time_myr[-1] * 1.1, 200)
        analytic_energy = compute_analytical_energy(
            t_analytic,
            energy_data['kinetic'][0],
            energy_data['thermal'][0],
            energy_data['potential'][0],
            energy_data['total'][0]
        )
    
    # =========================================================================
    # FIGURE SETUP
    # =========================================================================
    
    fig = plt.figure(figsize=(32, 18), facecolor=DARK_BG)
    
    fig.suptitle('IMBH Tidal Disruption (CoM Frame) - Cloud Stays Centered', 
                fontsize=26, fontweight='bold', color=TEXT_COLOR, y=0.98,
                path_effects=[path_effects.withStroke(linewidth=4, foreground=DARK_BG)])
    
    gs = GridSpec(2, 3, figure=fig, height_ratios=[1, 1],
                 hspace=0.25, wspace=0.2, left=0.04, right=0.96, top=0.92, bottom=0.06)
    
    ax_dens_xy = fig.add_subplot(gs[0, 0])
    ax_temp_xy = fig.add_subplot(gs[0, 1])
    ax_mach_xy = fig.add_subplot(gs[0, 2])
    ax_dens_xz = fig.add_subplot(gs[1, 0])
    ax_energy = fig.add_subplot(gs[1, 1])
    ax_info = fig.add_subplot(gs[1, 2])
    
    for ax in [ax_dens_xy, ax_temp_xy, ax_mach_xy, ax_dens_xz, ax_energy, ax_info]:
        ax.set_facecolor(DARK_BG)
    
    # =========================================================================
    # DENSITY XY PANEL
    # =========================================================================
    
    scatter_dens_xy = ax_dens_xy.scatter([], [], s=8, c=[], cmap=CUSTOM_CMAPS['density'],
                                         norm=LogNorm(vmin=1e-2, vmax=1e2),
                                         alpha=0.9, edgecolors='none')
    
    cbar_dens = fig.colorbar(scatter_dens_xy, ax=ax_dens_xy, shrink=0.85, pad=0.02)
    cbar_dens.set_label('Density', fontsize=11, color=TEXT_COLOR, fontweight='bold')
    cbar_dens.ax.tick_params(colors=TEXT_COLOR, labelsize=9)
    cbar_dens.outline.set_edgecolor(GRID_COLOR)
    
    # IMBH glow - will be updated dynamically since we're in CoM frame
    imbh_glow_dens_xy = []
    for s, alpha in [(350, 0.03), (280, 0.06), (200, 0.12), (120, 0.25)]:
        sc = ax_dens_xy.scatter([0], [0], s=s, c=ACCENT_COLOR, alpha=alpha, edgecolors='none', zorder=100)
        imbh_glow_dens_xy.append(sc)
    sc = ax_dens_xy.scatter([0], [0], s=150, c='white', marker='*',
                       edgecolors=ACCENT_COLOR, linewidths=2, zorder=101)
    imbh_glow_dens_xy.append(sc)
    
    style_dark_axis(ax_dens_xy, 'X - X_CoM (pc)', 'Y - Y_CoM (pc)', 'Density (CoM Frame)')
    ax_dens_xy.set_xlim(xlim)
    ax_dens_xy.set_ylim(ylim)
    ax_dens_xy.set_aspect('equal')
    
    time_text = ax_dens_xy.text(0.02, 0.98, '', transform=ax_dens_xy.transAxes,
                                fontsize=18, fontweight='bold', va='top', color=ACCENT_COLOR,
                                path_effects=[path_effects.withStroke(linewidth=4, foreground=DARK_BG)])
    
    # =========================================================================
    # TEMPERATURE PANEL  
    # =========================================================================
    
    scatter_temp = ax_temp_xy.scatter([], [], s=8, c=[], cmap=CUSTOM_CMAPS['temperature'],
                                       norm=LogNorm(vmin=1e2, vmax=1e7),
                                       alpha=0.9, edgecolors='none')
    
    cbar_temp = fig.colorbar(scatter_temp, ax=ax_temp_xy, shrink=0.85, pad=0.02)
    cbar_temp.set_label('Temperature (K)', fontsize=11, color=TEXT_COLOR, fontweight='bold')
    cbar_temp.ax.tick_params(colors=TEXT_COLOR, labelsize=9)
    cbar_temp.outline.set_edgecolor(GRID_COLOR)
    
    imbh_glow_temp_xy = []
    for s, alpha in [(350, 0.03), (280, 0.06), (200, 0.12), (120, 0.25)]:
        sc = ax_temp_xy.scatter([0], [0], s=s, c=HIGHLIGHT_COLOR, alpha=alpha, edgecolors='none', zorder=100)
        imbh_glow_temp_xy.append(sc)
    sc = ax_temp_xy.scatter([0], [0], s=150, c='white', marker='*',
                       edgecolors=HIGHLIGHT_COLOR, linewidths=2, zorder=101)
    imbh_glow_temp_xy.append(sc)
    
    style_dark_axis(ax_temp_xy, 'X - X_CoM (pc)', 'Y - Y_CoM (pc)', 'Temperature (CoM Frame)')
    ax_temp_xy.set_xlim(xlim)
    ax_temp_xy.set_ylim(ylim)
    ax_temp_xy.set_aspect('equal')
    
    # =========================================================================
    # INTERNAL MACH NUMBER PANEL (velocity relative to CoM)
    # =========================================================================
    
    scatter_mach = ax_mach_xy.scatter([], [], s=8, c=[], cmap=CUSTOM_CMAPS['mach'],
                                       norm=Normalize(vmin=0, vmax=5),
                                       alpha=0.9, edgecolors='none')
    
    cbar_mach = fig.colorbar(scatter_mach, ax=ax_mach_xy, shrink=0.85, pad=0.02)
    cbar_mach.set_label('Internal Mach', fontsize=11, color=TEXT_COLOR, fontweight='bold')
    cbar_mach.ax.tick_params(colors=TEXT_COLOR, labelsize=9)
    cbar_mach.outline.set_edgecolor(GRID_COLOR)
    
    imbh_glow_mach_xy = []
    for s, alpha in [(350, 0.03), (280, 0.06), (200, 0.12), (120, 0.25)]:
        sc = ax_mach_xy.scatter([0], [0], s=s, c='#00ff88', alpha=alpha, edgecolors='none', zorder=100)
        imbh_glow_mach_xy.append(sc)
    sc = ax_mach_xy.scatter([0], [0], s=150, c='white', marker='*',
                       edgecolors='#00ff88', linewidths=2, zorder=101)
    imbh_glow_mach_xy.append(sc)
    
    style_dark_axis(ax_mach_xy, 'X - X_CoM (pc)', 'Y - Y_CoM (pc)', 'Internal Mach (CoM Frame)')
    ax_mach_xy.set_xlim(xlim)
    ax_mach_xy.set_ylim(ylim)
    ax_mach_xy.set_aspect('equal')
    
    # =========================================================================
    # DENSITY XZ PANEL
    # =========================================================================
    
    scatter_dens_xz = ax_dens_xz.scatter([], [], s=8, c=[], cmap=CUSTOM_CMAPS['density'],
                                          norm=LogNorm(vmin=1e-2, vmax=1e2),
                                          alpha=0.9, edgecolors='none')
    
    cbar_xz = fig.colorbar(scatter_dens_xz, ax=ax_dens_xz, shrink=0.85, pad=0.02)
    cbar_xz.set_label('Density', fontsize=11, color=TEXT_COLOR, fontweight='bold')
    cbar_xz.ax.tick_params(colors=TEXT_COLOR, labelsize=9)
    cbar_xz.outline.set_edgecolor(GRID_COLOR)
    
    imbh_glow_dens_xz = []
    for s, alpha in [(350, 0.03), (280, 0.06), (200, 0.12), (120, 0.25)]:
        sc = ax_dens_xz.scatter([0], [0], s=s, c=ACCENT_COLOR, alpha=alpha, edgecolors='none', zorder=100)
        imbh_glow_dens_xz.append(sc)
    sc = ax_dens_xz.scatter([0], [0], s=150, c='white', marker='*',
                       edgecolors=ACCENT_COLOR, linewidths=2, zorder=101)
    imbh_glow_dens_xz.append(sc)
    
    style_dark_axis(ax_dens_xz, 'X - X_CoM (pc)', 'Z - Z_CoM (pc)', 'Density XZ (CoM Frame)')
    ax_dens_xz.set_xlim(xlim)
    ax_dens_xz.set_ylim(-1, 1)  # Z range matches shock diagnostics
    ax_dens_xz.set_aspect('equal')
    
    # =========================================================================
    # ENERGY EVOLUTION PANEL - THE KEY NEW FEATURE
    # =========================================================================
    
    style_dark_axis(ax_energy, 'Time (Myr)', 'Energy (code units)', 
                   'Energy Evolution with Analytical Prediction')
    
    # Set up energy axis limits based on data
    if energy_data is not None:
        E_max = max(energy_data['kinetic'].max(), energy_data['total'].max()) * 1.2
        E_min = energy_data['potential'].min() * 1.2
        ax_energy.set_xlim(0, energy_time_myr[-1] * 1.1)
        ax_energy.set_ylim(E_min, E_max)
    else:
        ax_energy.set_xlim(0, 5)
        ax_energy.set_ylim(-20, 60)
    
    # Plot analytical predictions (dashed lines - shown first as background)
    if analytic_energy is not None:
        ax_energy.plot(t_analytic, analytic_energy['total'], '--', 
                      color='#888888', linewidth=2, alpha=0.7, label='Analytic: Total (conserved)')
    
    # Simulation energy lines (will be updated)
    line_kinetic, = ax_energy.plot([], [], '-', color=COLORS['kinetic'], linewidth=2.5, 
                                    label='Kinetic (K)')
    line_thermal, = ax_energy.plot([], [], '-', color=COLORS['thermal'], linewidth=2.5,
                                    label='Thermal (U)')
    line_potential, = ax_energy.plot([], [], '-', color=COLORS['potential'], linewidth=2.5,
                                      label='Potential (W)')
    line_total, = ax_energy.plot([], [], '-', color=COLORS['total'], linewidth=3,
                                  label='Total (E)')
    
    # Current time marker
    energy_marker = ax_energy.axvline(x=0, color=ACCENT_COLOR, linewidth=2, linestyle='-', alpha=0.8)
    
    # Legend with physics explanation
    legend = ax_energy.legend(loc='upper right', fontsize=10, facecolor=DARK_PANEL, 
                              edgecolor=GRID_COLOR, labelcolor=TEXT_COLOR, ncol=2)
    
    # Energy explanation text
    energy_explain = ax_energy.text(0.02, 0.02, '', transform=ax_energy.transAxes,
                                     fontsize=9, va='bottom', ha='left', color=TEXT_COLOR,
                                     family='monospace',
                                     bbox=dict(boxstyle='round,pad=0.3', facecolor=DARK_PANEL,
                                              edgecolor=ACCENT_COLOR, alpha=0.9))
    
    # =========================================================================
    # PHYSICS INFO PANEL
    # =========================================================================
    
    ax_info.axis('off')
    
    info_text = ax_info.text(0.03, 0.97, '', transform=ax_info.transAxes,
                             fontsize=12, va='top', ha='left',
                             color=TEXT_COLOR, family='monospace',
                             linespacing=1.4,
                             bbox=dict(boxstyle='round,pad=0.5', facecolor=DARK_PANEL,
                                      edgecolor=GRID_COLOR, alpha=0.95))
    
    # =========================================================================
    # ANIMATION UPDATE FUNCTION
    # =========================================================================
    
    def update(frame_idx):
        snapshot = load_snapshot(snapshot_files[frame_idx])
        if snapshot is None:
            return []
        
        time_code = frame_idx * 0.025
        time_myr = time_code * TIME_UNIT
        
        # Compute derived quantities
        temperature = compute_temperature(snapshot['pres'], snapshot['dens'])
        mach, v_internal = compute_mach_number(
            snapshot['vel_x'], snapshot['vel_y'], snapshot['vel_z'], 
            snapshot['sound'], snapshot['mass'])
        com, r_xy, r_z, deform_ratio = compute_tidal_deformation(
            snapshot['pos_x'], snapshot['pos_y'], snapshot['pos_z'], snapshot['mass'])
        com_x, com_y, com_z = com
        
        # =====================================================================
        # TRANSFORM TO CoM FRAME: cloud stays centered, IMBH moves
        # =====================================================================
        pos_x_com = snapshot['pos_x'] - com_x
        pos_y_com = snapshot['pos_y'] - com_y
        pos_z_com = snapshot['pos_z'] - com_z
        
        # IMBH position in CoM frame (IMBH is at origin in lab frame)
        imbh_x_com = -com_x
        imbh_y_com = -com_y
        imbh_z_com = -com_z
        
        # Update IMBH marker positions (they move in CoM frame)
        for sc in imbh_glow_dens_xy:
            sc.set_offsets([[imbh_x_com, imbh_y_com]])
        for sc in imbh_glow_temp_xy:
            sc.set_offsets([[imbh_x_com, imbh_y_com]])
        for sc in imbh_glow_mach_xy:
            sc.set_offsets([[imbh_x_com, imbh_y_com]])
        for sc in imbh_glow_dens_xz:
            sc.set_offsets([[imbh_x_com, imbh_z_com]])

        # Update scatter plots - NOW IN CoM FRAME!
        scatter_dens_xy.set_offsets(np.column_stack([pos_x_com, pos_y_com]))
        scatter_dens_xy.set_array(snapshot['dens'])
        
        scatter_temp.set_offsets(np.column_stack([pos_x_com, pos_y_com]))
        scatter_temp.set_array(temperature)
        
        scatter_mach.set_offsets(np.column_stack([pos_x_com, pos_y_com]))
        scatter_mach.set_array(mach)
        
        scatter_dens_xz.set_offsets(np.column_stack([pos_x_com, pos_z_com]))
        scatter_dens_xz.set_array(snapshot['dens'])
        
        time_text.set_text(f't = {time_myr:.3f} Myr')
        
        # Update energy plot
        if energy_data is not None:
            # Find energy data up to current time
            mask = energy_time_myr <= time_myr + 0.01
            if np.any(mask):
                t_plot = energy_time_myr[mask]
                line_kinetic.set_data(t_plot, energy_data['kinetic'][mask])
                line_thermal.set_data(t_plot, energy_data['thermal'][mask])
                line_potential.set_data(t_plot, energy_data['potential'][mask])
                line_total.set_data(t_plot, energy_data['total'][mask])
                energy_marker.set_xdata([time_myr, time_myr])
                
                # Compute energy changes for explanation
                if len(t_plot) > 1:
                    dK = energy_data['kinetic'][mask][-1] - energy_data['kinetic'][0]
                    dU = energy_data['thermal'][mask][-1] - energy_data['thermal'][0]
                    dW = energy_data['potential'][mask][-1] - energy_data['potential'][0]
                    dE = energy_data['total'][mask][-1] - energy_data['total'][0]
                    dE_pct = abs(dE / energy_data['total'][0]) * 100 if energy_data['total'][0] != 0 else 0
                    
                    # Dynamic explanation based on what's happening
                    if time_myr < 0.5:
                        phase_explain = "Approach: Cloud falling into IMBH potential"
                    elif time_myr < 1.5:
                        phase_explain = "Perihelion: Max tidal force, shock heating"
                    elif time_myr < 2.5:
                        phase_explain = "Recession: Tidal energy deposited in streams"
                    else:
                        phase_explain = "Late: Continued expansion and cooling"
                    
                    energy_str = (
                        f"{phase_explain}\n"
                        f"ΔK = {dK:+.1f} (tidal accel)\n"
                        f"ΔU = {dU:+.1f} (shock heat)\n"  
                        f"ΔW = {dW:+.1f} (potential)\n"
                        f"ΔE = {dE:+.2f} ({dE_pct:.2f}% err)"
                    )
                    energy_explain.set_text(energy_str)
        
        # Shock statistics
        shocked, shock_param = detect_shocks(snapshot['dens'], snapshot['pres'])
        shock_frac = np.sum(shocked) / len(shocked) * 100
        supersonic_frac = np.sum(mach > 1) / len(mach) * 100
        
        # Info panel with physics explanation
        if time_myr < 0.5:
            phase = "APPROACH"
            phase_color = "#00ff00"
            physics = ("Cloud approaching IMBH\n"
                      "→ Potential energy decreasing\n"
                      "→ Kinetic energy increasing\n"
                      "→ Tidal stretching begins")
        elif time_myr < 1.5:
            phase = "PERIHELION"
            phase_color = "#ff6600"
            physics = ("Close encounter with IMBH\n"
                      "→ Maximum tidal forces\n"
                      "→ Strong shock heating\n"
                      "→ Thermal energy ↑↑")
        elif time_myr < 2.5:
            phase = "RECESSION"
            phase_color = "#ff00ff"
            physics = ("Moving away from IMBH\n"
                      "→ Tidal tails forming\n"
                      "→ Energy redistribution\n"
                      "→ Leading/trailing streams")
        else:
            phase = "LATE EVOLUTION"
            phase_color = "#00ffff"
            physics = ("Post-encounter dynamics\n"
                      "→ Continued stretching\n"
                      "→ Some gas may return\n"
                      "→ Stream expansion")
        
        info_str = (
            f"╔═══════════════════════════════════╗\n"
            f"║      {phase:^20}          ║\n"
            f"╠═══════════════════════════════════╣\n"
            f"║ t = {time_myr:>6.3f} Myr                  ║\n"
            f"║ Frame {frame_idx+1:>4}/{len(snapshot_files):<4}                  ║\n"
            f"╠═══════════════════════════════════╣\n"
            f"║ N_particles: {len(snapshot['dens']):>8,}           ║\n"
            f"║ ρ_max:    {snapshot['dens'].max():>10.2e}         ║\n"
            f"║ T_max:    {temperature.max():>10.2e} K       ║\n"
            f"║ M_max:    {mach.max():>10.2f}           ║\n"
            f"╠═══════════════════════════════════╣\n"
            f"║ Shocked:  {shock_frac:>6.1f}%                 ║\n"
            f"║ M > 1:    {supersonic_frac:>6.1f}%                 ║\n"
            f"║ R_xy/R_z: {deform_ratio:>6.2f}                  ║\n"
            f"╠═══════════════════════════════════╣\n"
            f"║         PHYSICS                   ║\n"
            f"╠═══════════════════════════════════╣\n"
            f"║ {physics.replace(chr(10), chr(10) + '║ '):<33} ║\n"
            f"╚═══════════════════════════════════╝"
        )
        info_text.set_text(info_str)
        
        if frame_idx % 20 == 0:
            print(f"  Frame {frame_idx+1}/{len(snapshot_files)}: t={time_myr:.3f} Myr")
        
        return [scatter_dens_xy, scatter_temp, scatter_mach, scatter_dens_xz,
                time_text, line_kinetic, line_thermal, line_potential, line_total,
                energy_marker, energy_explain, info_text]
    
    # =========================================================================
    # CREATE AND SAVE ANIMATION
    # =========================================================================
    
    print(f"Creating tidal disruption animation with energy evolution...")
    print(f"  Frames: {len(snapshot_files)}")
    print(f"  FPS: {fps}")
    print(f"  Output: {output_file}")
    
    anim = FuncAnimation(fig, update, frames=len(snapshot_files),
                        interval=1000//fps, blit=False)
    
    print(f"Saving animation...")
    writer = PillowWriter(fps=fps)
    anim.save(output_file, writer=writer, dpi=100)
    plt.close(fig)
    
    print(f"✓ Animation saved: {output_file}")
    return True


# =============================================================================
# MAIN FUNCTION
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='IMBH Tidal Disruption Animation - High Contrast with Energy Evolution',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Full tidal disruption animation with energy panel
  python animate_single_run.py results/ -o tidal.gif --mode tidal
  
  # Single density animation  
  python animate_single_run.py results/ -o density.gif --mode density

Energy Evolution Explanation:
  The energy panel shows how energy is exchanged during the encounter:
  
  • KINETIC (green): Increases as cloud accelerates toward IMBH
  • THERMAL (orange): Increases from shock heating (irreversible)
  • POTENTIAL (magenta): Decreases (more negative) during approach
  • TOTAL (white): Should be conserved - check for numerical errors
  
  Gray dashed line shows analytical prediction for total energy conservation.
        """
    )
    
    parser.add_argument('results_dir', type=str,
                       help='Directory containing snapshot CSV files')
    parser.add_argument('-o', '--output', type=str, default='animation.gif',
                       help='Output GIF filename')
    parser.add_argument('--mode', type=str, default='tidal',
                       choices=['tidal', 'density', 'temperature', 'mach', 'shock'],
                       help='Visualization mode')
    parser.add_argument('--fps', type=int, default=15, help='Frames per second')
    parser.add_argument('--xlim', type=float, nargs=2, default=None)
    parser.add_argument('--ylim', type=float, nargs=2, default=None)
    
    args = parser.parse_args()
    
    xlim = tuple(args.xlim) if args.xlim else None
    ylim = tuple(args.ylim) if args.ylim else None
    
    print("=" * 65)
    print("IMBH Tidal Disruption - High Contrast with Energy Evolution")
    print("=" * 65)
    print(f"Results: {args.results_dir}")
    print(f"Output:  {args.output}")
    print(f"Mode:    {args.mode}")
    print("=" * 65)
    
    if args.mode == 'tidal':
        success = create_tidal_disruption_animation(
            args.results_dir, args.output, xlim=xlim, ylim=ylim, fps=args.fps)
    else:
        success = create_single_mode_animation(
            args.results_dir, args.output, mode=args.mode,
            xlim=xlim, ylim=ylim, fps=args.fps)
    
    if success:
        print("\n" + "=" * 65)
        print("✓ Animation complete!")
        print("=" * 65)
    else:
        print("\n" + "=" * 65)
        print("✗ Animation failed!")
        print("=" * 65)
        return 1
    
    return 0


if __name__ == '__main__':
    exit(main())
