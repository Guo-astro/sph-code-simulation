#!/usr/bin/env python3
"""
SR-GSPH Sod Shock Tube Visualization - Single Source of Truth (SSOT)

This is the authoritative visualization script for SR Sod shock tube results.
It provides:
  - Animation of time evolution
  - Static plots (initial, final, evolution)
  - Proper physical units display
  - Timestep and simulation time tracking
  - Support for all test types (Sod, blast wave, strong blast, etc.)
  - Analytic solution overlay for benchmark comparison

Reference: Kitajima et al. (2025) arXiv:2510.18251v1

Usage:
    python sr_sod_visualize.py <results_dir> [--animate] [--plot] [--format gif|mp4]
    python sr_sod_visualize.py <results_dir> --analytic  # Enable analytic solution overlay
    
Examples:
    python sr_sod_visualize.py ../results/kitajima --animate
    python sr_sod_visualize.py ../results/blast_wave --plot --animate --format mp4
    python sr_sod_visualize.py ../results/sod_kitajima --analytic --animate  # Benchmark comparison
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter, FFMpegWriter
from pathlib import Path
import json
import sys
import argparse
from dataclasses import dataclass
from typing import Optional, List, Tuple, Dict, Any

# Add parent directories to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent / 'docs'))
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent / 'scripts'))

try:
    from relativistic_riemann_solver import RelativisiticRiemannSolver
    RIEMANN_SOLVER_AVAILABLE = True
except ImportError:
    RIEMANN_SOLVER_AVAILABLE = False
    print("Warning: Relativistic Riemann solver not found. Analytic solution disabled.")

# Import unit system
try:
    from units.unit_system import RelativisticUnits
    UNIT_SYSTEM_AVAILABLE = True
except ImportError:
    UNIT_SYSTEM_AVAILABLE = False
    print("Note: Unit system module not found. Using default units.")

# =============================================================================
# Physical Units Configuration (Code Units for SR hydrodynamics)
# =============================================================================
@dataclass
class PhysicalUnits:
    """
    Physical units configuration for SR-GSPH simulations.
    
    If the unit system module is available, this class can be initialized
    from a config file or CSV header. Otherwise, it uses default code units.
    """
    # In code units: c = 1, so velocities are in units of c
    c_speed: float = 1.0           # Speed of light [code units]
    length_unit: str = "L"         # Length unit label
    time_unit: str = "L/c"         # Time unit label  
    density_unit: str = "n₀"       # Rest-frame density unit
    pressure_unit: str = "P₀"      # Pressure unit
    velocity_unit: str = "c"       # Velocity in units of c
    energy_unit: str = "E₀"        # Energy unit
    
    # Conversion factors to CGS (optional, for physical unit display)
    length_to_cgs: float = 1.0     # cm
    time_to_cgs: float = 1.0       # s
    density_to_cgs: float = 1.0    # g/cm³
    velocity_to_cgs: float = 1.0   # cm/s
    pressure_to_cgs: float = 1.0   # dyn/cm²
    
    @classmethod
    def from_config(cls, config: Dict[str, Any]) -> 'PhysicalUnits':
        """
        Create PhysicalUnits from a config dictionary.
        Uses the new unit system if available.
        """
        if UNIT_SYSTEM_AVAILABLE:
            units = RelativisticUnits.from_config(config)
            return cls(
                c_speed=units.c_code,
                length_unit=units.length_label,
                time_unit=units.time_label,
                density_unit=units.density_label,
                pressure_unit=units.pressure_label,
                velocity_unit=units.velocity_label,
                energy_unit=units.energy_label,
                length_to_cgs=units.length_to_cgs,
                time_to_cgs=units.time_to_cgs,
                density_to_cgs=units.density_to_cgs,
                velocity_to_cgs=units.velocity_to_cgs,
                pressure_to_cgs=units.pressure_to_cgs,
            )
        else:
            # Fallback to defaults
            return cls()
    
    @classmethod
    def from_csv_header(cls, metadata: Dict[str, str]) -> 'PhysicalUnits':
        """
        Create PhysicalUnits from CSV metadata.
        Uses the new unit system if available.
        """
        if UNIT_SYSTEM_AVAILABLE:
            units = RelativisticUnits.from_csv_header(metadata)
            return cls(
                c_speed=units.c_code,
                length_unit=units.length_label,
                time_unit=units.time_label,
                density_unit=units.density_label,
                pressure_unit=units.pressure_label,
                velocity_unit=units.velocity_label,
                energy_unit=units.energy_label,
                length_to_cgs=units.length_to_cgs,
                time_to_cgs=units.time_to_cgs,
                density_to_cgs=units.density_to_cgs,
                velocity_to_cgs=units.velocity_to_cgs,
                pressure_to_cgs=units.pressure_to_cgs,
            )
        else:
            # Parse from metadata manually
            c_code = float(metadata.get('c_code', '1.0'))
            return cls(c_speed=c_code)
    
    def format_time(self, t: float) -> str:
        """Format time with proper units"""
        return f"t = {t:.4f} [{self.time_unit}]"
    
    def format_velocity(self, v: float) -> str:
        """Format velocity with units"""
        return f"{v:.4f} c"
    
    def to_physical_time(self, t_code: float) -> float:
        """Convert code units time to physical (CGS) time"""
        return t_code * self.time_to_cgs
    
    def to_physical_length(self, x_code: float) -> float:
        """Convert code units length to physical (CGS) length"""
        return x_code * self.length_to_cgs
    
    def to_physical_density(self, rho_code: float) -> float:
        """Convert code units density to physical (CGS) density"""
        return rho_code * self.density_to_cgs
    
    def to_physical_velocity(self, v_code: float) -> float:
        """Convert code units velocity to physical (CGS) velocity"""
        return v_code * self.velocity_to_cgs
    
    def to_physical_pressure(self, p_code: float) -> float:
        """Convert code units pressure to physical (CGS) pressure"""
        return p_code * self.pressure_to_cgs


# =============================================================================
# Analytic Solution Provider
# =============================================================================
class AnalyticSolution:
    """
    Provides exact analytic solutions for SR shock tube problems.
    Uses the relativistic Riemann solver from Marti & Mueller (1994).
    """
    
    # Initial conditions for each test type (from Kitajima et al. 2025)
    TEST_CONDITIONS = {
        'sod': {
            'P_left': 1.0, 'n_left': 1.0, 'v_left': 0.0,
            'P_right': 0.1, 'n_right': 0.125, 'v_right': 0.0,
            'gamma': 5.0/3.0,
        },
        'blast_wave': {
            'P_left': 40.0/3.0, 'n_left': 10.0, 'v_left': 0.0,
            'P_right': 1.0e-6, 'n_right': 1.0, 'v_right': 0.0,
            'gamma': 5.0/3.0,
        },
        'strong_blast': {
            'P_left': 1000.0, 'n_left': 1.0, 'v_left': 0.0,
            'P_right': 0.01, 'n_right': 1.0, 'v_right': 0.0,
            'gamma': 5.0/3.0,
        },
        'ultra_relativistic': {
            'P_left': 1.0, 'n_left': 1.0, 'v_left': 0.9,
            'P_right': 10.0, 'n_right': 1.0, 'v_right': 0.0,
            'gamma': 5.0/3.0,
        },
    }
    
    def __init__(self, test_type: str = 'sod', gamma: float = 5.0/3.0):
        """
        Initialize analytic solution provider.
        
        Args:
            test_type: Type of test ('sod', 'blast_wave', 'strong_blast', 'ultra_relativistic')
            gamma: Adiabatic index (default 5/3 for ideal gas)
        """
        if not RIEMANN_SOLVER_AVAILABLE:
            raise ImportError("Relativistic Riemann solver not available")
        
        # Get initial conditions
        self.test_type = test_type
        conditions = self.TEST_CONDITIONS.get(test_type, self.TEST_CONDITIONS['sod'])
        
        self.P_left = conditions['P_left']
        self.n_left = conditions['n_left']
        self.v_left = conditions['v_left']
        self.P_right = conditions['P_right']
        self.n_right = conditions['n_right']
        self.v_right = conditions['v_right']
        self.gamma = gamma if gamma else conditions['gamma']
        
        # Initialize solver
        self.solver = RelativisiticRiemannSolver(self.gamma)
        self.solver.set_initial_states(
            self.P_left, self.n_left, self.v_left,
            self.P_right, self.n_right, self.v_right
        )
    
    def compute(self, t: float, n_points: int = 500) -> Dict[str, np.ndarray]:
        """
        Compute analytic solution at time t.
        
        Args:
            t: Time for the solution
            n_points: Number of spatial points
            
        Returns:
            Dictionary with 'x', 'dens', 'pres', 'vel' arrays
        """
        if t <= 0:
            # Return initial conditions
            x = np.linspace(0, 1, n_points)
            dens = np.where(x < 0.5, self.n_left, self.n_right)
            pres = np.where(x < 0.5, self.P_left, self.P_right)
            vel = np.where(x < 0.5, self.v_left, self.v_right)
            return {'x': x - 0.5, 'dens': dens, 'pres': pres, 'vel': vel}
        
        try:
            x, pressure, density, velocity, _ = self.solver.solve(t, n=n_points)
            # Shift x from [0,1] to [-0.5, 0.5] to match simulation domain
            return {
                'x': x - 0.5,
                'dens': density,
                'pres': pressure,
                'vel': velocity
            }
        except Exception as e:
            print(f"Warning: Analytic solution failed at t={t}: {e}")
            return None


# =============================================================================
# Data I/O Functions
# =============================================================================
def read_snapshot(filepath: Path) -> Optional[Dict[str, np.ndarray]]:
    """
    Read CSV snapshot file with proper header parsing.
    
    Returns dict with arrays: pos_x, vel_x, dens, pres, ene, etc.
    """
    data_rows = []
    header = None
    metadata = {}
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('#'):
                # Parse metadata comments
                if ':' in line:
                    key, _, value = line[1:].partition(':')
                    metadata[key.strip()] = value.strip()
                continue
            if line.startswith('id,'):
                header = line.split(',')
                continue
            data_rows.append(line.split(','))
    
    if not data_rows or header is None:
        return None
    
    # Convert to numpy arrays
    result = {'_metadata': metadata}
    n_particles = len(data_rows)
    
    for col_idx, col_name in enumerate(header):
        values = np.zeros(n_particles)
        for row_idx, row in enumerate(data_rows):
            try:
                values[row_idx] = float(row[col_idx])
            except (ValueError, IndexError):
                values[row_idx] = np.nan
        result[col_name] = values
    
    return result


def read_config(config_path: Path) -> Dict[str, Any]:
    """Read JSON configuration file"""
    try:
        with open(config_path, 'r') as f:
            return json.load(f)
    except Exception as e:
        print(f"Warning: Could not read config {config_path}: {e}")
        return {}


def find_config(results_dir: Path) -> Optional[Dict[str, Any]]:
    """Find and read configuration file for the simulation"""
    results_dir = results_dir.resolve()
    
    # Check for config.json in results directory first
    direct_config = results_dir / 'config.json'
    if direct_config.exists():
        config = read_config(direct_config)
        if config:
            return config
    
    # Check parent directories for config
    for parent in [results_dir.parent, results_dir.parent.parent]:
        parent_config = parent / 'config.json'
        if parent_config.exists():
            config = read_config(parent_config)
            if config:
                return config
    
    # Search preset configs by matching outputDirectory to results_dir
    preset_dir = results_dir.parent.parent / 'config' / 'presets'
    if preset_dir.exists():
        for preset_path in preset_dir.glob('*.json'):
            config = read_config(preset_path)
            if config:
                out_dir = config.get('outputDirectory', '')
                if out_dir:
                    # Normalize and check if this config's output matches results_dir
                    out_path = Path(out_dir).resolve() if Path(out_dir).is_absolute() else (preset_dir.parent.parent.parent / out_dir).resolve()
                    if out_path == results_dir or results_dir.name in str(out_dir):
                        return config
    
    # Fallback: try to match by directory name pattern in preset names
    dir_name = results_dir.name.lower()
    if preset_dir.exists():
        # Priority order for matching
        patterns = [
            ('sod_kitajima', 'sr_sod_kitajima.json'),
            ('sod_sharp', 'sr_sod_sharp.json'),
            ('sod', 'sr_sod.json'),
            ('strong_blast', 'sr_strong_blast.json'),
            ('blast', 'sr_blast_wave.json'),
            ('ultra', 'sr_ultra_relat.json'),
        ]
        for pattern, preset_name in patterns:
            if pattern in dir_name:
                preset_path = preset_dir / preset_name
                if preset_path.exists():
                    return read_config(preset_path)
    
    return None


def detect_test_info(results_dir: Path, config: Optional[Dict] = None) -> Dict[str, str]:
    """Detect test type and solver information from directory name or config"""
    info = {
        'test_type': 'SR Sod Shock Tube',
        'solver': 'EXACT',
        'solver_name': 'Kitajima et al. (2025)',
        'reference': 'arXiv:2510.18251v1'
    }
    
    dir_name = results_dir.name.lower()
    
    # Detect test type from directory name
    if 'blast' in dir_name and 'strong' in dir_name:
        info['test_type'] = 'Strong Relativistic Blast Wave'
        info['section'] = 'Section 3.1.3'
    elif 'blast' in dir_name:
        info['test_type'] = 'Relativistic Blast Wave'
        info['section'] = 'Section 3.1.2'
    elif 'ultra' in dir_name:
        info['test_type'] = 'Ultra-Relativistic Shock Tube'
    elif 'kitajima' in dir_name or 'sod' in dir_name:
        info['test_type'] = 'Sod Shock Tube'
        info['section'] = 'Section 3.1.1'
    
    # Detect solver type
    if 'iterative' in dir_name:
        info['solver'] = 'ITERATIVE'
        info['solver_name'] = 'van Leer Relativistic'
    elif 'hllc' in dir_name:
        info['solver'] = 'HLLC'
        info['solver_name'] = 'HLLC Approximate'
    
    # Override with config if available
    if config:
        if 'testType' in config:
            test_type = config['testType']
            if test_type == 'sod':
                info['test_type'] = 'Sod Shock Tube (Section 3.1.1)'
            elif test_type == 'blast_wave':
                info['test_type'] = 'Relativistic Blast Wave (Section 3.1.2)'
            elif test_type == 'strong_blast':
                info['test_type'] = 'Strong Blast Wave (Section 3.1.3)'
        
        if 'riemannSolverSRGSPH' in config:
            info['solver'] = config['riemannSolverSRGSPH']
    
    return info


# =============================================================================
# Visualization Class
# =============================================================================
class SRSodVisualizer:
    """
    Single Source of Truth visualizer for SR-GSPH Sod shock tube simulations.
    Supports analytic solution overlay for benchmark comparison.
    Now supports config-driven unit systems.
    """
    
    def __init__(self, results_dir: str | Path, enable_analytic: bool = False):
        self.results_dir = Path(results_dir).resolve()
        self.enable_analytic = enable_analytic and RIEMANN_SOLVER_AVAILABLE
        self.analytic_solver = None
        
        # Verify directory exists
        if not self.results_dir.exists():
            # Try relative to script location
            script_dir = Path(__file__).parent
            alt_path = (script_dir.parent / results_dir).resolve()
            if alt_path.exists():
                self.results_dir = alt_path
            else:
                raise FileNotFoundError(f"Results directory not found: {results_dir}\n"
                                       f"  Tried: {Path(results_dir).resolve()}\n"
                                       f"  Tried: {alt_path}")
        
        # Load configuration
        self.config = find_config(self.results_dir) or {}
        self.test_info = detect_test_info(self.results_dir, self.config)
        
        # Get output time interval from config
        self.dt_output = self.config.get('outputTime', 0.01)
        self.gamma = self.config.get('gamma', 5.0/3.0)
        
        # Find snapshots
        self.snapshot_files = sorted(self.results_dir.glob('snapshot_*.csv'))
        if not self.snapshot_files:
            raise FileNotFoundError(f"No snapshot files found in {self.results_dir}")
        
        print(f"Found {len(self.snapshot_files)} snapshots")
        print(f"Test type: {self.test_info['test_type']}")
        print(f"Riemann solver: {self.test_info['solver']} ({self.test_info['solver_name']})")
        
        # Initialize units - try from CSV first, then config, then defaults
        self.units = self._init_units()
        
        # Initialize analytic solver if requested
        if self.enable_analytic:
            self._init_analytic_solver()
        
        # Load all data
        self.snapshots = []
        self.times = []
        self._load_all_snapshots()
        
        # Compute global ranges for consistent axes
        self._compute_global_ranges()
    
    def _init_units(self) -> PhysicalUnits:
        """
        Initialize units from available sources:
        1. CSV header metadata (has actual simulation units)
        2. Config file (json)  
        3. Defaults
        """
        # Try to read units from first snapshot
        if self.snapshot_files:
            data = read_snapshot(self.snapshot_files[0])
            if data and '_metadata' in data:
                metadata = data['_metadata']
                if 'unit_type' in metadata or 'c_code' in metadata:
                    units = PhysicalUnits.from_csv_header(metadata)
                    print(f"Units: Loaded from CSV ({units.c_speed=})")
                    return units
        
        # Try config
        if self.config and 'units' in self.config:
            units = PhysicalUnits.from_config(self.config)
            print(f"Units: Loaded from config ({units.c_speed=})")
            return units
        
        # Defaults
        print("Units: Using defaults (c=1 code units)")
        return PhysicalUnits()
    
    def _init_analytic_solver(self):
        """Initialize the analytic Riemann solver based on test type"""
        # Detect test type from directory or config
        test_type = 'sod'  # default
        dir_name = self.results_dir.name.lower()
        
        if 'blast' in dir_name and 'strong' in dir_name:
            test_type = 'strong_blast'
        elif 'blast' in dir_name:
            test_type = 'blast_wave'
        elif 'ultra' in dir_name:
            test_type = 'ultra_relativistic'
        
        # Override from config
        if 'testType' in self.config:
            cfg_type = self.config['testType']
            if cfg_type == 'blast_wave':
                test_type = 'blast_wave'
            elif cfg_type == 'strong_blast':
                test_type = 'strong_blast'
        
        try:
            self.analytic_solver = AnalyticSolution(test_type, self.gamma)
            print(f"✓ Analytic solution enabled for test type: {test_type}")
        except Exception as e:
            print(f"Warning: Could not initialize analytic solver: {e}")
            self.enable_analytic = False
    
    def _load_all_snapshots(self):
        """Load all snapshot data"""
        for snap_file in self.snapshot_files:
            data = read_snapshot(snap_file)
            if data is not None:
                self.snapshots.append(data)
                # Extract snapshot number and compute time
                snap_num = int(snap_file.stem.split('_')[-1])
                self.times.append(snap_num * self.dt_output)
        
        print(f"Loaded {len(self.snapshots)} snapshots")
        print(f"Time range: t = {self.times[0]:.4f} to {self.times[-1]:.4f} [{self.units.time_unit}]")
    
    def _compute_global_ranges(self):
        """Compute global min/max for consistent axis scaling"""
        all_dens = np.concatenate([s['dens'] for s in self.snapshots])
        all_vel = np.concatenate([s['vel_x'] for s in self.snapshots])
        all_pres = np.concatenate([s['pres'] for s in self.snapshots])
        
        margin = 0.1  # 10% margin
        
        self.dens_range = (
            all_dens.min() * (1 - margin),
            all_dens.max() * (1 + margin)
        )
        self.vel_range = (
            all_vel.min() - abs(all_vel.max() - all_vel.min()) * margin,
            all_vel.max() + abs(all_vel.max() - all_vel.min()) * margin
        )
        self.pres_range = (
            max(all_pres.min() * 0.5, 1e-10),
            all_pres.max() * 2
        )
    
    def _create_figure(self) -> Tuple[plt.Figure, np.ndarray]:
        """Create figure with consistent styling"""
        fig, axes = plt.subplots(3, 1, figsize=(12, 10))
        fig.patch.set_facecolor('#fafafa')
        
        # Colors for each quantity
        self.colors = {
            'dens': '#2563eb',   # Blue
            'pres': '#16a34a',   # Green
            'vel': '#dc2626',    # Red
        }
        
        return fig, axes
    
    def _setup_axes(self, axes: np.ndarray, show_units: bool = True):
        """Configure axes with labels and styling"""
        labels = [
            ('Rest-Frame Density', 'n', self.units.density_unit, self.colors['dens']),
            ('Pressure', 'P', self.units.pressure_unit, self.colors['pres']),
            ('Velocity', 'v', self.units.velocity_unit, self.colors['vel']),
        ]
        
        for ax, (name, symbol, unit, color) in zip(axes, labels):
            if show_units:
                ylabel = f'{name} {symbol} [{unit}]'
            else:
                ylabel = f'{name} {symbol}'
            
            ax.set_ylabel(ylabel, fontsize=12, fontweight='bold', color=color)
            ax.set_facecolor('white')
            ax.grid(True, alpha=0.3, linestyle='--')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_color(color)
            ax.spines['left'].set_linewidth(1.5)
            ax.tick_params(axis='y', colors=color)
        
        # X-axis label on bottom plot
        axes[-1].set_xlabel(f'Position x [{self.units.length_unit}]', fontsize=12, fontweight='bold')
        
        # Set axis limits
        for ax in axes:
            ax.set_xlim(-0.55, 0.55)
        
        axes[0].set_ylim(*self.dens_range)
        axes[1].set_yscale('log')
        axes[1].set_ylim(*self.pres_range)
        axes[2].set_ylim(*self.vel_range)
        axes[2].axhline(y=0, color='black', linestyle='-', alpha=0.5, linewidth=1)
    
    def _make_title(self) -> str:
        """Generate plot title with test info"""
        title = f"SR-GSPH: {self.test_info['test_type']}\n"
        title += f"Riemann Solver: {self.test_info['solver']} ({self.test_info['solver_name']})"
        return title
    
    def plot_snapshot(self, frame_idx: int, axes: np.ndarray, 
                     show_particles: bool = True, show_line: bool = True):
        """Plot a single snapshot on given axes"""
        data = self.snapshots[frame_idx]
        
        # Sort by position
        sort_idx = np.argsort(data['pos_x'])
        x = data['pos_x'][sort_idx]
        dens = data['dens'][sort_idx]
        pres = data['pres'][sort_idx]
        vel = data['vel_x'][sort_idx]
        
        quantities = [
            (dens, 'dens'),
            (pres, 'pres'),
            (vel, 'vel'),
        ]
        
        artists = []
        for ax, (y, key) in zip(axes, quantities):
            if show_particles:
                scatter = ax.scatter(data['pos_x'], data[key.replace('vel', 'vel_x')], 
                                    s=8, alpha=0.6, c=self.colors[key],
                                    edgecolors='white', linewidths=0.3)
                artists.append(scatter)
            
            if show_line:
                line, = ax.plot(x, y, color=self.colors[key], linewidth=1.5, alpha=0.8)
                artists.append(line)
        
        return artists
    
    def create_animation(self, output_format: str = 'gif', fps: int = 12) -> Path:
        """
        Create animation showing time evolution with optional analytic solution overlay.
        
        Parameters:
            output_format: 'gif' or 'mp4'
            fps: Frames per second
        
        Returns:
            Path to output file
        """
        analytic_mode = self.enable_analytic and self.analytic_solver is not None
        mode_str = "with ANALYTIC overlay" if analytic_mode else ""
        print(f"\nCreating animation ({output_format.upper()}) {mode_str}...")
        
        fig, axes = self._create_figure()
        self._setup_axes(axes)
        
        # Title and time display
        title = self._make_title()
        if analytic_mode:
            title += "\n(Black line = Exact solution)"
        title_text = fig.suptitle(title, fontsize=14, fontweight='bold', y=0.98)
        
        # Time and timestep info box
        info_box = fig.text(0.5, 0.93, '', ha='center', fontsize=11,
                           bbox=dict(boxstyle='round,pad=0.5', facecolor='white', 
                                    edgecolor='#333', alpha=0.9, linewidth=1.5),
                           family='monospace')
        
        # Precompute sorted data for efficiency
        sorted_data = []
        for data in self.snapshots:
            sort_idx = np.argsort(data['pos_x'])
            sorted_data.append({
                'x': data['pos_x'][sort_idx],
                'x_raw': data['pos_x'],
                'dens': data['dens'][sort_idx],
                'dens_raw': data['dens'],
                'pres': data['pres'][sort_idx],
                'pres_raw': data['pres'],
                'vel': data['vel_x'][sort_idx],
                'vel_raw': data['vel_x'],
            })
        
        # Precompute analytic solutions if enabled
        analytic_data = []
        if analytic_mode:
            print("  Computing analytic solutions for each timestep...")
            for t in self.times:
                sol = self.analytic_solver.compute(t, n_points=500)
                analytic_data.append(sol)
        
        # Initialize plot elements - simulation data
        scatter_dens = axes[0].scatter([], [], s=8, alpha=0.6, c=self.colors['dens'],
                                       edgecolors='white', linewidths=0.3, label='SPH')
        scatter_pres = axes[1].scatter([], [], s=8, alpha=0.6, c=self.colors['pres'],
                                       edgecolors='white', linewidths=0.3)
        scatter_vel = axes[2].scatter([], [], s=8, alpha=0.6, c=self.colors['vel'],
                                      edgecolors='white', linewidths=0.3)
        
        line_dens, = axes[0].plot([], [], color=self.colors['dens'], linewidth=1.5, alpha=0.8)
        line_pres, = axes[1].plot([], [], color=self.colors['pres'], linewidth=1.5, alpha=0.8)
        line_vel, = axes[2].plot([], [], color=self.colors['vel'], linewidth=1.5, alpha=0.8)
        
        # Initialize analytic solution lines (black dashed)
        analytic_lines = []
        if analytic_mode:
            line_dens_exact, = axes[0].plot([], [], 'k-', linewidth=2, alpha=0.9, label='Exact')
            line_pres_exact, = axes[1].plot([], [], 'k-', linewidth=2, alpha=0.9)
            line_vel_exact, = axes[2].plot([], [], 'k-', linewidth=2, alpha=0.9)
            analytic_lines = [line_dens_exact, line_pres_exact, line_vel_exact]
            axes[0].legend(loc='upper right', fontsize=9)
        
        plt.tight_layout(rect=[0, 0, 1, 0.91])
        
        def animate(frame_idx):
            data = sorted_data[frame_idx]
            t = self.times[frame_idx]
            
            # Update scatter plots
            scatter_dens.set_offsets(np.c_[data['x_raw'], data['dens_raw']])
            scatter_pres.set_offsets(np.c_[data['x_raw'], data['pres_raw']])
            scatter_vel.set_offsets(np.c_[data['x_raw'], data['vel_raw']])
            
            # Update line plots
            line_dens.set_data(data['x'], data['dens'])
            line_pres.set_data(data['x'], data['pres'])
            line_vel.set_data(data['x'], data['vel'])
            
            # Update analytic solution if enabled
            if analytic_mode and analytic_data[frame_idx] is not None:
                sol = analytic_data[frame_idx]
                analytic_lines[0].set_data(sol['x'], sol['dens'])
                analytic_lines[1].set_data(sol['x'], sol['pres'])
                analytic_lines[2].set_data(sol['x'], sol['vel'])
            
            # Update time info with physical units
            dt = self.times[1] - self.times[0] if len(self.times) > 1 else self.dt_output
            info_text = (f"Time: t = {t:.4f} [{self.units.time_unit}]  |  "
                        f"Δt_out = {dt:.4f}  |  "
                        f"Frame {frame_idx + 1}/{len(self.snapshots)}")
            info_box.set_text(info_text)
            
            result = [scatter_dens, scatter_pres, scatter_vel, line_dens, line_pres, line_vel, info_box]
            result.extend(analytic_lines)
            return result
        
        # Create animation
        anim = animation.FuncAnimation(
            fig, animate, frames=len(self.snapshots),
            interval=1000 // fps, blit=True, repeat=True
        )
        
        # Save
        output_name = f"sr_sod_animation_{self.test_info['solver'].lower()}"
        if analytic_mode:
            output_name += "_vs_exact"
        
        if output_format.lower() == 'mp4':
            output_file = self.results_dir / f'{output_name}.mp4'
            try:
                writer = FFMpegWriter(fps=fps, bitrate=2500, codec='libx264',
                                     extra_args=['-pix_fmt', 'yuv420p'])
                anim.save(str(output_file), writer=writer, dpi=120)
                print(f"✓ MP4 saved: {output_file}")
            except Exception as e:
                print(f"FFmpeg error: {e}")
                print("Falling back to GIF...")
                output_format = 'gif'
        
        if output_format.lower() == 'gif':
            output_file = self.results_dir / f'{output_name}.gif'
            writer = PillowWriter(fps=fps)
            anim.save(str(output_file), writer=writer, dpi=100)
            print(f"✓ GIF saved: {output_file}")
        
        plt.close(fig)
        
        # Print summary
        file_size = output_file.stat().st_size / 1024 / 1024
        print(f"  File size: {file_size:.2f} MB")
        print(f"  Frames: {len(self.snapshots)}")
        print(f"  Time range: t = {self.times[0]:.4f} to {self.times[-1]:.4f}")
        if analytic_mode:
            print(f"  ✓ Includes exact analytic solution overlay")
        
        return output_file
    
    def plot_final_state(self) -> Path:
        """Plot the final state of the simulation"""
        print("\nPlotting final state...")
        
        fig, axes = self._create_figure()
        self._setup_axes(axes)
        
        # Plot final snapshot
        self.plot_snapshot(-1, axes)
        
        # Title with time
        t_final = self.times[-1]
        title = self._make_title() + f"\nFinal State: t = {t_final:.4f} [{self.units.time_unit}]"
        fig.suptitle(title, fontsize=13, fontweight='bold')
        
        plt.tight_layout(rect=[0, 0, 1, 0.93])
        
        output_file = self.results_dir / 'sr_sod_final_state.png'
        plt.savefig(output_file, dpi=150, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        
        print(f"✓ Final state saved: {output_file}")
        return output_file
    
    def plot_evolution(self, n_times: int = 6) -> Path:
        """Plot evolution at multiple time slices"""
        print("\nPlotting time evolution...")
        
        fig, axes = self._create_figure()
        self._setup_axes(axes, show_units=True)
        
        # Select time indices
        n_times = min(n_times, len(self.snapshots))
        indices = np.linspace(0, len(self.snapshots) - 1, n_times, dtype=int)
        
        # Color gradient for time evolution
        cmap = plt.cm.viridis
        colors_time = [cmap(i / (n_times - 1)) for i in range(n_times)]
        
        for idx, color in zip(indices, colors_time):
            data = self.snapshots[idx]
            sort_idx = np.argsort(data['pos_x'])
            x = data['pos_x'][sort_idx]
            t = self.times[idx]
            
            label = f't = {t:.3f}'
            axes[0].plot(x, data['dens'][sort_idx], color=color, linewidth=1.5, 
                        alpha=0.8, label=label)
            axes[1].plot(x, data['pres'][sort_idx], color=color, linewidth=1.5, alpha=0.8)
            axes[2].plot(x, data['vel_x'][sort_idx], color=color, linewidth=1.5, alpha=0.8)
        
        axes[0].legend(loc='upper right', fontsize=9, framealpha=0.9)
        
        fig.suptitle(self._make_title() + '\nTime Evolution', fontsize=13, fontweight='bold')
        plt.tight_layout(rect=[0, 0, 1, 0.93])
        
        output_file = self.results_dir / 'sr_sod_evolution.png'
        plt.savefig(output_file, dpi=150, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        
        print(f"✓ Evolution plot saved: {output_file}")
        return output_file


# =============================================================================
# Main Entry Point
# =============================================================================
def main():
    parser = argparse.ArgumentParser(
        description='SR-GSPH Sod Shock Tube Visualization (SSOT)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s results/kitajima --animate
  %(prog)s results/blast_wave --plot --animate --format mp4
  %(prog)s results/strong_blast --all
  %(prog)s results/sod_kitajima --analytic --animate  # Benchmark with exact solution
        """
    )
    
    parser.add_argument('results_dir', type=str,
                       help='Directory containing snapshot_*.csv files')
    parser.add_argument('--animate', '-a', action='store_true',
                       help='Create animation')
    parser.add_argument('--plot', '-p', action='store_true',
                       help='Create static plots (final state + evolution)')
    parser.add_argument('--all', action='store_true',
                       help='Create all visualizations')
    parser.add_argument('--analytic', '--exact', action='store_true',
                       help='Enable analytic (exact) solution overlay for benchmark comparison')
    parser.add_argument('--format', '-f', choices=['gif', 'mp4'], default='gif',
                       help='Animation format (default: gif)')
    parser.add_argument('--fps', type=int, default=12,
                       help='Animation frames per second (default: 12)')
    
    args = parser.parse_args()
    
    # Default to all if nothing specified
    if not (args.animate or args.plot or args.all):
        args.all = True
    
    # Check if analytic solution is available
    if args.analytic and not RIEMANN_SOLVER_AVAILABLE:
        print("⚠️  Warning: Analytic solution requested but Riemann solver not available")
        print("   Install from: docs/relativistic_riemann_solver.py")
        args.analytic = False
    
    try:
        viz = SRSodVisualizer(args.results_dir, enable_analytic=args.analytic)
        
        outputs = []
        
        if args.plot or args.all:
            outputs.append(viz.plot_final_state())
            outputs.append(viz.plot_evolution())
        
        if args.animate or args.all:
            outputs.append(viz.create_animation(output_format=args.format, fps=args.fps))
        
        print("\n" + "=" * 60)
        print("✅ Visualization complete!")
        print("=" * 60)
        print("Output files:")
        for f in outputs:
            print(f"  📁 {f}")
        
    except FileNotFoundError as e:
        print(f"❌ Error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
