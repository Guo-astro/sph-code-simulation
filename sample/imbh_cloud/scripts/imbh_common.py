#!/usr/bin/env python3
"""
IMBH Visualization Common Module
================================

Shared utilities for all IMBH tidal disruption visualization scripts.
Provides consistent colormaps, physics calculations, and data loading.

This is the Single Source of Truth (SSOT) for:
- Physical constants and unit conversions
- Colormap definitions
- Data loading functions
- Physics calculations (temperature, Mach, COM, etc.)
- Analytical energy model

Units: IMBH_ENCOUNTER ([L]=pc, [M]=1000 M_sun, [V]=km/s, [T]=0.978 Myr)
"""

import numpy as np
import glob
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)

# =============================================================================
# PHYSICAL CONSTANTS & UNIT CONVERSIONS
# =============================================================================

TIME_UNIT = 0.978  # 1 code time = 0.978 Myr
VELOCITY_UNIT = 1.0  # km/s
LENGTH_UNIT = 1.0  # pc
MASS_UNIT = 1000.0  # M_sun

GAMMA = 5.0 / 3.0
MU = 2.3  # Mean molecular weight
M_H = 1.67e-24  # Hydrogen mass (g)
K_B = 1.38e-16  # Boltzmann constant (erg/K)
VELOCITY_TO_CGS = 1e5  # km/s to cm/s
TEMP_CONVERSION = MU * M_H / K_B * VELOCITY_TO_CGS**2

# =============================================================================
# DARK THEME COLOR PALETTE
# =============================================================================

DARK_BG = '#0d0d12'
DARK_PANEL = '#15151f'
GRID_COLOR = '#3a3a4a'
TEXT_COLOR = '#f0f0f0'
ACCENT_COLOR = '#00ffff'
HIGHLIGHT_COLOR = '#ff4444'
BH_COLOR = '#ff4444'
TRAJECTORY_COLOR = '#ffaa00'

# Energy panel colors
ENERGY_COLORS = {
    'kinetic': '#00ff00',      # Bright green
    'thermal': '#ff6600',      # Bright orange  
    'potential': '#ff00ff',    # Bright magenta
    'total': '#ffffff',        # White
    'analytic': '#ffff00',     # Yellow
}

# 1D profile line colors (colorblind-friendly)
LINE_COLORS = {
    'density': '#00ff00',    # Pure bright green
    'pressure': '#ff8800',   # Bright orange  
    'temperature': '#ff00ff',  # Magenta
    'mach': '#00ccff',       # Cyan
    'initial': '#888888',    # Gray (reference)
}


# =============================================================================
# COLORMAP DEFINITIONS
# =============================================================================

def hex_to_rgb(hex_color):
    """Convert hex color to RGB tuple (0-1 range)."""
    h = hex_color.lstrip('#')
    return tuple(int(h[i:i+2], 16)/255 for i in (0, 2, 4))


def create_colormap_dict(colors, name):
    """Create matplotlib colormap dict from list of (position, hex_color) tuples."""
    from matplotlib.colors import LinearSegmentedColormap
    
    positions = [c[0] for c in colors]
    rgb_colors = [hex_to_rgb(c[1]) for c in colors]
    
    cdict = {'red': [], 'green': [], 'blue': []}
    for pos, (r, g, b) in zip(positions, rgb_colors):
        cdict['red'].append((pos, r, r))
        cdict['green'].append((pos, g, g))
        cdict['blue'].append((pos, b, b))
    
    return LinearSegmentedColormap(name, cdict, N=256)


# Colormap color definitions
DENSITY_COLORS = [
    (0.0, '#000010'),    # Near black
    (0.1, '#000080'),    # Deep blue
    (0.25, '#0066ff'),   # Bright blue
    (0.4, '#00ffff'),    # Cyan
    (0.55, '#00ff66'),   # Green-cyan
    (0.7, '#ffff00'),    # Yellow
    (0.85, '#ff8800'),   # Orange
    (1.0, '#ffffff'),    # White
]

TEMPERATURE_COLORS = [
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

MACH_COLORS = [
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

SHOCK_COLORS = [
    (0.0, '#100020'),    # Dark purple
    (0.2, '#440066'),    # Purple
    (0.4, '#880088'),    # Magenta
    (0.6, '#cc0044'),    # Red-magenta
    (0.8, '#ff6600'),    # Orange
    (1.0, '#ffff00'),    # Yellow
]

VELOCITY_COLORS = [
    (0.0, '#1a237e'),   # Deep indigo (stationary)
    (0.2, '#3949ab'),   # Indigo
    (0.35, '#5c6bc0'),  # Light indigo
    (0.5, '#26c6da'),   # Cyan (moderate)
    (0.65, '#66bb6a'),  # Green
    (0.8, '#ffca28'),   # Amber
    (0.9, '#ff7043'),   # Deep orange
    (1.0, '#d32f2f'),   # Red (fast)
]

SHOCK_TEMP_COLORS = [
    (0.0, '#0044ff'),    # Cold - blue
    (0.2, '#0088ff'),    # Light blue
    (0.4, '#00cccc'),    # Cyan
    (0.5, '#ffffff'),    # T = T0 (initial) - WHITE
    (0.6, '#ffff00'),    # Slightly heated - yellow
    (0.75, '#ff8800'),   # Moderately shocked - orange
    (0.9, '#ff0000'),    # Strong shock - red
    (1.0, '#ff00ff'),    # Extreme shock - magenta
]


def get_colormaps():
    """Get dictionary of all custom colormaps."""
    return {
        'density': create_colormap_dict(DENSITY_COLORS, 'density_hc'),
        'temperature': create_colormap_dict(TEMPERATURE_COLORS, 'temp_hc'),
        'mach': create_colormap_dict(MACH_COLORS, 'mach_hc'),
        'shock': create_colormap_dict(SHOCK_COLORS, 'shock_hc'),
        'velocity': create_colormap_dict(VELOCITY_COLORS, 'velocity_hc'),
        'shock_temp': create_colormap_dict(SHOCK_TEMP_COLORS, 'shock_temp'),
    }


def get_plotly_colorscale(mode):
    """Get Plotly-compatible colorscale for a given mode."""
    colorscales = {
        'density': [
            [0.0, '#0d0221'], [0.15, '#1a237e'], [0.30, '#0277bd'],
            [0.45, '#00acc1'], [0.55, '#26a69a'], [0.65, '#66bb6a'],
            [0.75, '#ffeb3b'], [0.85, '#ff9800'], [0.92, '#f44336'], [1.0, '#ffffff']
        ],
        'temperature': [
            [0.0, '#000033'], [0.12, '#0d47a1'], [0.25, '#1976d2'],
            [0.35, '#00bcd4'], [0.45, '#4caf50'], [0.55, '#8bc34a'],
            [0.65, '#ffeb3b'], [0.75, '#ff9800'], [0.85, '#f44336'],
            [0.92, '#e91e63'], [1.0, '#ffffff']
        ],
        'velocity': [
            [0.0, '#1a237e'], [0.2, '#3949ab'], [0.35, '#5c6bc0'],
            [0.5, '#26c6da'], [0.65, '#66bb6a'], [0.8, '#ffca28'],
            [0.9, '#ff7043'], [1.0, '#d32f2f']
        ],
        'mach': [
            [0.0, '#0000aa'], [0.15, '#0066ff'], [0.3, '#00ccff'],
            [0.45, '#66ffff'], [0.5, '#ffffff'], [0.55, '#ffff66'],
            [0.7, '#ffcc00'], [0.85, '#ff6600'], [1.0, '#ff0000']
        ],
        'shock': [
            [0.0, '#100020'], [0.2, '#440066'], [0.4, '#880088'],
            [0.6, '#cc0044'], [0.8, '#ff6600'], [1.0, '#ffff00']
        ],
    }
    return colorscales.get(mode, colorscales['density'])


# =============================================================================
# DATA LOADING
# =============================================================================

def load_snapshot(filepath):
    """
    Load a single CSV snapshot file.
    
    Returns:
    --------
    dict with keys: pos_x, pos_y, pos_z, vel_x, vel_y, vel_z, dens, pres, mass, sml, sound
    Optional: ene (specific internal energy)
    """
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
            'sound': data['sound'],
        }
        
        if 'sml' in data.dtype.names:
            snapshot['sml'] = data['sml']
        if 'ene' in data.dtype.names:
            snapshot['ene'] = data['ene']
        
        return snapshot
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return None


def load_snapshots(results_dir, verbose=True):
    """
    Load all snapshot files from a results directory.
    
    Returns:
    --------
    list of snapshot dicts, sorted by time
    """
    snapshot_files = sorted(glob.glob(f"{results_dir}/snapshot_*.csv"))
    if not snapshot_files:
        print(f"Error: No snapshots found in {results_dir}")
        return []
    
    if verbose:
        print(f"Found {len(snapshot_files)} snapshots")
    
    all_data = []
    for i, f in enumerate(snapshot_files):
        data = load_snapshot(f)
        if data:
            all_data.append(data)
            if verbose and (i + 1) % 10 == 0:
                print(f"  Loaded {i + 1}/{len(snapshot_files)} snapshots")
    
    if verbose:
        print(f"  Loaded {len(all_data)} valid snapshots")
    
    return all_data


def load_energy_file(filepath):
    """
    Load energy evolution data from energy.dat file.
    
    Returns:
    --------
    dict with keys: time, kinetic, thermal, potential, total
    """
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
# PHYSICS CALCULATIONS
# =============================================================================

def compute_temperature(pres, dens):
    """Compute temperature from pressure and density."""
    specific_energy = pres / (dens + 1e-30)
    temperature = specific_energy * TEMP_CONVERSION / GAMMA
    return temperature


def compute_velocity_magnitude(vel_x, vel_y, vel_z, mass=None):
    """Compute velocity magnitude relative to center-of-mass motion."""
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


def compute_mach_number(vel_x, vel_y, vel_z, sound, mass=None):
    """
    Compute INTERNAL Mach number relative to center-of-mass motion.
    
    This removes the bulk motion of the cloud, so we see only
    internal velocity dispersion and shock velocities.
    """
    v_internal = compute_velocity_magnitude(vel_x, vel_y, vel_z, mass)
    mach = v_internal / (sound + 1e-30)
    return mach


def compute_center_of_mass(pos_x, pos_y, pos_z, mass):
    """Compute center of mass position."""
    total_mass = np.sum(mass)
    com_x = np.sum(mass * pos_x) / total_mass
    com_y = np.sum(mass * pos_y) / total_mass
    com_z = np.sum(mass * pos_z) / total_mass
    return com_x, com_y, com_z


def compute_shock_indicator(temp, temp_initial=50.0):
    """
    Compute shock indicator as temperature ratio T/T_initial.
    
    Values > 1 indicate shock heating.
    """
    return temp / temp_initial


def compute_tidal_radius(cloud_mass, bh_mass, cloud_radius=0.56):
    """
    Compute tidal radius where IMBH tidal force equals cloud self-gravity.
    
    r_tidal = r_cloud * (M_BH / M_cloud)^(1/3)
    
    Parameters:
    -----------
    cloud_mass : float - cloud mass in code units (1 = 1000 M_sun)
    bh_mass : float - IMBH mass in code units (100 = 10^5 M_sun)
    cloud_radius : float - cloud radius in pc
    
    Returns:
    --------
    r_tidal in pc
    """
    return cloud_radius * (bh_mass / cloud_mass) ** (1/3)


# =============================================================================
# GLOBAL COLOR RANGE COMPUTATION
# =============================================================================

def compute_global_ranges(all_data, percentile_low=1, percentile_high=99):
    """
    Compute global color ranges from all snapshots for consistent colorbars.
    
    Returns:
    --------
    dict with ranges for density, temperature, velocity, mach, shock
    """
    all_dens = np.concatenate([d['dens'] for d in all_data])
    all_temp = np.concatenate([compute_temperature(d['pres'], d['dens']) for d in all_data])
    all_vel = np.concatenate([compute_velocity_magnitude(d['vel_x'], d['vel_y'], d['vel_z'], d['mass']) 
                              for d in all_data])
    all_mach = np.concatenate([compute_mach_number(d['vel_x'], d['vel_y'], d['vel_z'], d['sound'], d['mass']) 
                               for d in all_data])
    all_shock = np.concatenate([compute_shock_indicator(compute_temperature(d['pres'], d['dens'])) 
                                for d in all_data])
    
    return {
        'density': {
            'min': np.log10(max(np.percentile(all_dens, percentile_low), 1e-10)),
            'max': np.log10(np.percentile(all_dens, percentile_high)),
        },
        'temperature': {
            'min': np.log10(max(np.percentile(all_temp, percentile_low), 1.0)),
            'max': np.log10(np.percentile(all_temp, percentile_high)),
        },
        'velocity': {
            'min': 0,
            'max': np.percentile(all_vel, percentile_high) * 1.1,
        },
        'mach': {
            'min': 0,
            'max': min(np.percentile(all_mach, percentile_high) * 1.2, 10.0),
        },
        'shock': {
            'min': 0,
            'max': min(np.percentile(all_shock, percentile_high) * 1.2, 100.0),
        },
    }


# =============================================================================
# ANALYTICAL ENERGY MODEL
# =============================================================================

def compute_analytical_energy(t, E_kin_0, E_therm_0, E_pot_0, E_tot_0,
                               v_cloud=10.0, b=3.0, M_BH=100.0, x_0=-20.0):
    """
    Compute analytical energy predictions for IMBH tidal encounter.
    
    Physics:
    ========
    IMBH is stationary at origin. Cloud starts at (x_0, b, 0) and moves with
    velocity (v_cloud, 0, 0). The cloud falls into the IMBH's gravitational
    potential, converting potential energy to kinetic energy.
    
    Energy is conserved: E_total = KE + PE + U_thermal = const
    
    Parameters:
    -----------
    t : array - Time in Myr
    E_kin_0, E_therm_0, E_pot_0, E_tot_0 : float - Initial energies
    v_cloud : float - Cloud velocity (km/s)
    b : float - Impact parameter (pc)
    M_BH : float - IMBH mass in code units (100 = 10^5 M_sun)
    x_0 : float - Initial x position (pc)
    
    Returns:
    --------
    dict with 'kinetic', 'thermal', 'potential', 'total', 'r_to_BH', 't_peri_myr'
    """
    # Convert time to code units
    t_code = t / TIME_UNIT
    
    # Cloud trajectory (straight line)
    x_cloud = x_0 + v_cloud * t_code
    r_to_BH = np.sqrt(x_cloud**2 + b**2)
    
    # Time of closest approach
    t_peri_code = abs(x_0) / v_cloud
    t_peri_myr = t_peri_code * TIME_UNIT
    
    # Initial distance to BH
    r_0 = np.sqrt(x_0**2 + b**2)
    
    # IMBH potential change
    PE_BH = -M_BH / r_to_BH
    PE_BH_0 = -M_BH / r_0
    delta_PE_BH = PE_BH - PE_BH_0
    
    # Energy evolution (exact energy conservation)
    E_pot_analytic = E_pot_0 + delta_PE_BH
    E_kin_analytic = E_kin_0 - delta_PE_BH
    
    # Thermal energy decreases due to adiabatic expansion
    perihelion_factor = np.minimum(r_0 / r_to_BH, 3.0)
    expansion_cooling = 1.0 / (1.0 + 0.5 * (perihelion_factor - 1.0)**2)
    E_therm_analytic = E_therm_0 * expansion_cooling
    
    # Total energy
    E_tot_analytic = E_kin_analytic + E_therm_analytic + E_pot_analytic
    
    return {
        'kinetic': E_kin_analytic,
        'thermal': E_therm_analytic,
        'potential': E_pot_analytic,
        'total': E_tot_analytic,
        'r_to_BH': r_to_BH,
        't_peri_myr': t_peri_myr,
    }


# =============================================================================
# 1D PROFILE EXTRACTION
# =============================================================================

def extract_vertical_profile(snapshot, com, column_radius=0.1, n_bins=50):
    """
    Extract 1D vertical (Z) profile from a column of particles near the CoM.
    Shows VERTICAL COMPRESSION structure.
    """
    com_x, com_y, com_z = com
    
    dx = snapshot['pos_x'] - com_x
    dy = snapshot['pos_y'] - com_y
    r_xy = np.sqrt(dx**2 + dy**2)
    
    in_column = r_xy < column_radius
    
    if np.sum(in_column) < 10:
        column_radius *= 2
        in_column = r_xy < column_radius
    
    z_col = snapshot['pos_z'][in_column] - com_z
    dens_col = snapshot['dens'][in_column]
    pres_col = snapshot['pres'][in_column]
    mass_col = snapshot['mass'][in_column]
    sound_col = snapshot['sound'][in_column]
    temp_col = compute_temperature(pres_col, dens_col)
    
    total_mass = np.sum(snapshot['mass'])
    v_com_z = np.sum(snapshot['mass'] * snapshot['vel_z']) / total_mass
    vel_z_rel = snapshot['vel_z'][in_column] - v_com_z
    mach_z = np.abs(vel_z_rel) / (sound_col + 1e-30)
    
    z_range = np.max(np.abs(z_col)) * 1.1
    z_edges = np.linspace(-z_range, z_range, n_bins + 1)
    z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])
    
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
    Shows IN-PLANE TIDAL STRETCHING structure.
    """
    com_x, com_y, com_z = com
    
    dz = np.abs(snapshot['pos_z'] - com_z)
    dy = np.abs(snapshot['pos_y'] - com_y)
    in_slice = (dz < slice_thickness) & (dy < slice_thickness)
    
    if np.sum(in_slice) < 10:
        slice_thickness *= 2
        in_slice = (dz < slice_thickness) & (dy < slice_thickness)
    
    x_slice = snapshot['pos_x'][in_slice] - com_x
    dens_slice = snapshot['dens'][in_slice]
    pres_slice = snapshot['pres'][in_slice]
    mass_slice = snapshot['mass'][in_slice]
    sound_slice = snapshot['sound'][in_slice]
    temp_slice = compute_temperature(pres_slice, dens_slice)
    
    total_mass = np.sum(snapshot['mass'])
    v_com_x = np.sum(snapshot['mass'] * snapshot['vel_x']) / total_mass
    vel_x_rel = snapshot['vel_x'][in_slice] - v_com_x
    mach_x = np.abs(vel_x_rel) / (sound_slice + 1e-30)
    
    x_range = np.max(np.abs(x_slice)) * 1.1
    if x_range < 0.1:
        x_range = 1.0
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
# COLORBAR TITLES
# =============================================================================

COLORBAR_TITLES = {
    'density': 'Log₁₀(ρ)',
    'temperature': 'Log₁₀(T [K])',
    'velocity': 'Internal v (km/s)',
    'mach': 'Mach Number',
    'shock': 'Shock Indicator (T/T₀)',
}


if __name__ == '__main__':
    print("IMBH Common Module - Shared utilities for visualization scripts")
    print("=" * 60)
    print("Physical Constants:")
    print(f"  TIME_UNIT = {TIME_UNIT} Myr")
    print(f"  GAMMA = {GAMMA}")
    print(f"  MU = {MU}")
    print(f"  TEMP_CONVERSION = {TEMP_CONVERSION:.3e}")
    print("")
    print("Available colormaps:", list(get_colormaps().keys()))
    print("Available Plotly colorscales:", ['density', 'temperature', 'velocity', 'mach', 'shock'])
