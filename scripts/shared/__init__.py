# Shared physics modules for SPH visualization and analysis
# Single Source of Truth (SSOT) for common functionality

from .lane_emden import (
    # Constants
    XI_1_N15,
    DTHETA_1_N15,
    KNOWN_XI_1,
    # Data classes
    LaneEmdenSolution,
    ICParameters,
    # Core functions - data file based (preferred)
    load_lane_emden_solution,
    get_density_profile,
    get_pressure_profile,
    get_analytic_profile,
    # IC presets
    get_ic_preset,
    list_ic_presets,
    IC_PRESETS,
    # ODE solvers (fallback when data files unavailable)
    solve_lane_emden,
    solve_lane_emden_spherical,
    solve_lane_emden_cylindrical,
    solve_lane_emden_planar,
    get_density_profile_from_solver,
    # Utility functions
    lane_emden_mass,
    lane_emden_alpha,
    lane_emden_profile,  # Deprecated but kept for compatibility
)

from .snapshot import (
    # Data classes
    SnapshotMetadata,
    # Core functions
    read_snapshot,
    parse_header,
    get_time_from_snapshot,
    find_snapshots,
    load_snapshot_series,
    # Utility functions
    compute_radii,
    bin_radial_profile,
)

__all__ = [
    # Lane-Emden
    'XI_1_N15', 'DTHETA_1_N15', 'KNOWN_XI_1',
    'LaneEmdenSolution', 'ICParameters',
    'load_lane_emden_solution', 'get_density_profile', 'get_pressure_profile',
    'get_analytic_profile', 'get_ic_preset', 'list_ic_presets', 'IC_PRESETS',
    'solve_lane_emden', 'solve_lane_emden_spherical', 
    'solve_lane_emden_cylindrical', 'solve_lane_emden_planar',
    'get_density_profile_from_solver',
    'lane_emden_mass', 'lane_emden_alpha', 'lane_emden_profile',
    # Snapshot
    'SnapshotMetadata',
    'read_snapshot', 'parse_header', 'get_time_from_snapshot',
    'find_snapshots', 'load_snapshot_series',
    'compute_radii', 'bin_radial_profile',
]
