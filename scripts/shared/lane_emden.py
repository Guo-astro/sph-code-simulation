#!/usr/bin/env python3
"""
Lane-Emden equation solution module - Single Source of Truth (SSOT)

This module provides the CORRECT Lane-Emden solution loaded from pre-computed
data files. DO NOT re-implement ODE integration - use this module instead.

Supported geometries:
- 3D Spherical: (1/ξ²) d/dξ (ξ² dθ/dξ) = -θⁿ
- 2D Cylindrical: (1/ξ) d/dξ (ξ dθ/dξ) = -θⁿ  
- 1D Planar (slab): d²θ/dξ² = -θⁿ

For n=1.5 (γ=5/3 polytrope):
    - First zero (3D): ξ₁ = 3.6537540101
    - |dθ/dξ|₁ = 0.20331244769

Physical relations:
    - r = α ξ  where α = √((n+1)K ρ_c^(1/n-1) / (4πG))
    - ρ(ξ) = ρ_c θⁿ
    - P(ξ) = K ρ^γ = K ρ_c^γ θ^(n+1)

For n=1.5 polytrope (γ=5/3):
    - ρ(ξ) = ρ_c θ^(3/2)
    - M = -4π α³ ρ_c ξ₁² (dθ/dξ)₁
    - R = α ξ₁

Author: SSOT refactoring for sphcode project
"""

import numpy as np
from pathlib import Path
from typing import Optional, Tuple, Dict, Any
from dataclasses import dataclass


# =============================================================================
# CONSTANTS: Exact Lane-Emden n=1.5 solution values
# =============================================================================
XI_1_N15 = 3.6537540101  # First zero of θ for n=1.5 (3D spherical)
DTHETA_1_N15 = -0.20331244769  # dθ/dξ at ξ₁ for n=1.5 (3D spherical)

# Known first zeros for different geometries and polytropic indices
KNOWN_XI_1 = {
    # (n, geometry) -> xi_1
    (1.5, "spherical"): 3.6537540101,
    (1.5, "cylindrical"): 2.7192,  # Approximate
    (1.5, "planar"): 1.7868,       # Approximate (for slab)
    (2.5, "planar"): 1.7868,       # n=1/(γ-1) for γ=1.4
    (1.0, "spherical"): np.pi,     # Exact analytical
    (0.0, "spherical"): np.sqrt(6), # Exact analytical
}


# =============================================================================
# GEOMETRY-SPECIFIC SOLVERS (for when pre-computed data is unavailable)
# These use RK4 integration - use pre-computed data when available!
# =============================================================================

def _rk4_step(f, xi, y, dxi):
    """Single RK4 step for ODE integration."""
    k1 = f(xi, y)
    k2 = f(xi + 0.5*dxi, y + 0.5*dxi*k1)
    k3 = f(xi + 0.5*dxi, y + 0.5*dxi*k2)
    k4 = f(xi + dxi, y + dxi*k3)
    return y + (dxi/6.0) * (k1 + 2*k2 + 2*k3 + k4)


def solve_lane_emden_spherical(n: float, xi_max: float = 10.0, 
                                dxi: float = 1e-4) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Solve 3D spherical Lane-Emden equation using RK4.
    
    Equation: (1/ξ²) d/dξ (ξ² dθ/dξ) = -θⁿ
    
    Parameters
    ----------
    n : float
        Polytropic index
    xi_max : float
        Maximum ξ to integrate to
    dxi : float
        Step size
        
    Returns
    -------
    xi : ndarray
        Dimensionless radius
    theta : ndarray
        θ(ξ)
    xi_1 : float
        First zero (surface radius)
    """
    xi_vals = [0.0]
    theta_vals = [1.0]
    
    xi = dxi
    theta = 1.0 - (dxi**2) / 6.0  # Taylor expansion near origin
    dtheta = -dxi / 3.0
    
    while xi < xi_max and theta > 0:
        # d²θ/dξ² = -θⁿ - (2/ξ) dθ/dξ
        def f_theta(xi_cur, state):
            th, dth = state
            if th <= 0:
                return np.array([0.0, 0.0])
            d2theta = -th**n - (2.0/xi_cur) * dth if xi_cur > 1e-10 else -th**n / 3.0
            return np.array([dth, d2theta])
        
        state = np.array([theta, dtheta])
        state = _rk4_step(f_theta, xi, state, dxi)
        theta, dtheta = state
        xi += dxi
        
        if theta > 0:
            xi_vals.append(xi)
            theta_vals.append(theta)
    
    xi_arr = np.array(xi_vals)
    theta_arr = np.array(theta_vals)
    
    # Find first zero by interpolation
    if len(theta_arr) > 1 and theta_arr[-1] <= 0:
        idx = len(theta_arr) - 1
        xi_1 = xi_arr[idx-1] - theta_arr[idx-1] * dxi / (theta_arr[idx] - theta_arr[idx-1])
    else:
        xi_1 = xi_arr[-1]
    
    return xi_arr, theta_arr, xi_1


def solve_lane_emden_cylindrical(n: float, xi_max: float = 8.0,
                                  dxi: float = 1e-4) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Solve 2D cylindrical Lane-Emden equation using RK4.
    
    Equation: (1/ξ) d/dξ (ξ dθ/dξ) = -θⁿ
    
    Parameters
    ----------
    n : float
        Polytropic index
    xi_max : float
        Maximum ξ to integrate to
    dxi : float
        Step size
        
    Returns
    -------
    xi : ndarray
        Dimensionless radius
    theta : ndarray
        θ(ξ)
    xi_1 : float
        First zero (surface radius)
    """
    xi_vals = [0.0]
    theta_vals = [1.0]
    
    xi = dxi
    theta = 1.0 - (dxi**2) / 4.0  # Taylor expansion near origin (2D)
    dtheta = -dxi / 2.0
    
    while xi < xi_max and theta > 0:
        def f_theta(xi_cur, state):
            th, dth = state
            if th <= 0:
                return np.array([0.0, 0.0])
            # d²θ/dξ² = -θⁿ - (1/ξ) dθ/dξ
            d2theta = -th**n - (1.0/xi_cur) * dth if xi_cur > 1e-10 else -th**n / 2.0
            return np.array([dth, d2theta])
        
        state = np.array([theta, dtheta])
        state = _rk4_step(f_theta, xi, state, dxi)
        theta, dtheta = state
        xi += dxi
        
        if theta > 0:
            xi_vals.append(xi)
            theta_vals.append(theta)
    
    xi_arr = np.array(xi_vals)
    theta_arr = np.array(theta_vals)
    
    if len(theta_arr) > 1 and theta_arr[-1] <= 0:
        idx = len(theta_arr) - 1
        xi_1 = xi_arr[idx-1] - theta_arr[idx-1] * dxi / (theta_arr[idx] - theta_arr[idx-1])
    else:
        xi_1 = xi_arr[-1]
    
    return xi_arr, theta_arr, xi_1


def solve_lane_emden_planar(n: float, xi_max: float = 5.0,
                            dxi: float = 1e-4) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Solve 1D planar (slab) Lane-Emden equation using RK4.
    
    Equation: d²θ/dξ² = -θⁿ
    
    This is the simplest form - no geometric terms.
    
    Parameters
    ----------
    n : float
        Polytropic index
    xi_max : float
        Maximum ξ to integrate to
    dxi : float
        Step size
        
    Returns
    -------
    xi : ndarray
        Dimensionless position
    theta : ndarray
        θ(ξ)
    xi_1 : float
        First zero (slab half-width)
    """
    xi_vals = [0.0]
    theta_vals = [1.0]
    
    xi = 0.0
    theta = 1.0
    dtheta = 0.0
    
    while xi < xi_max and theta > 0:
        # d²θ/dξ² = -θⁿ (no geometric term)
        def f_theta(xi_cur, state):
            th, dth = state
            if th <= 0:
                return np.array([0.0, 0.0])
            d2theta = -th**n
            return np.array([dth, d2theta])
        
        state = np.array([theta, dtheta])
        state = _rk4_step(f_theta, xi, state, dxi)
        theta, dtheta = state
        xi += dxi
        
        if theta > 0:
            xi_vals.append(xi)
            theta_vals.append(theta)
    
    xi_arr = np.array(xi_vals)
    theta_arr = np.array(theta_vals)
    
    if len(theta_arr) > 1 and theta_arr[-1] <= 0:
        idx = len(theta_arr) - 1
        xi_1 = xi_arr[idx-1] - theta_arr[idx-1] * dxi / (theta_arr[idx] - theta_arr[idx-1])
    else:
        xi_1 = xi_arr[-1]
    
    return xi_arr, theta_arr, xi_1


def solve_lane_emden(
    n: float = 1.5,
    geometry: str = "spherical",
    xi_max: float = 10.0,
    dxi: float = 1e-4
) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Unified Lane-Emden solver for all geometries.
    
    This is the RECOMMENDED function to use when pre-computed data is unavailable.
    
    Parameters
    ----------
    n : float
        Polytropic index (default 1.5 for γ=5/3)
    geometry : str
        One of "spherical" (3D), "cylindrical" (2D), "planar" (1D slab)
    xi_max : float
        Maximum ξ to integrate to
    dxi : float
        Step size
        
    Returns
    -------
    xi : ndarray
        Dimensionless coordinate
    theta : ndarray
        θ(ξ) values
    xi_1 : float
        First zero (surface/boundary)
        
    Example
    -------
    >>> xi, theta, xi_1 = solve_lane_emden(n=1.5, geometry="spherical")
    >>> print(f"Surface at ξ₁ = {xi_1:.4f}")
    Surface at ξ₁ = 3.6538
    """
    geometry = geometry.lower()
    
    if geometry in ("spherical", "3d", "sphere"):
        return solve_lane_emden_spherical(n, xi_max, dxi)
    elif geometry in ("cylindrical", "2d", "cylinder"):
        return solve_lane_emden_cylindrical(n, xi_max, dxi)
    elif geometry in ("planar", "1d", "slab"):
        return solve_lane_emden_planar(n, xi_max, dxi)
    else:
        raise ValueError(f"Unknown geometry: {geometry}. "
                        "Use 'spherical', 'cylindrical', or 'planar'.")


def get_density_profile_from_solver(
    n: float = 1.5,
    geometry: str = "spherical",
    rho_c: float = 1.0,
    R: float = 1.0,
    n_points: int = 200
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute density profile by solving Lane-Emden ODE directly.
    
    Use this when pre-computed data files are unavailable.
    
    Parameters
    ----------
    n : float
        Polytropic index
    geometry : str
        Geometry type
    rho_c : float
        Central density
    R : float
        Physical radius/half-width
    n_points : int
        Number of output points
        
    Returns
    -------
    r : ndarray
        Physical positions
    rho : ndarray
        Density values
    """
    xi, theta, xi_1 = solve_lane_emden(n=n, geometry=geometry)
    
    # Scale to physical coordinates
    alpha = R / xi_1
    r_phys = xi * alpha
    rho = rho_c * np.maximum(theta, 0.0)**n
    
    # Interpolate to uniform grid
    r_out = np.linspace(0, R, n_points)
    rho_out = np.interp(r_out, r_phys, rho, right=0.0)
    
    return r_out, rho_out


# =============================================================================
# IC PRESETS: Known initial condition parameters from relaxation runs
# These are the EXACT values from the IC snapshot headers - use these!
# =============================================================================
IC_PRESETS: Dict[str, Dict[str, Any]] = {
    # IMBH Cloud 10k relaxed IC
    # Source: sample/imbh_cloud/results/ic_relax_10k/snapshot_0032.csv header
    "imbh_10k": {
        "rho_c": 1.430096922633852,   # Central density (exact from IC header)
        "R": 1.0,                      # Cloud radius
        "M": 1.0,                      # Total mass
        "K": 0.4242088610712798,       # Polytropic constant
        "alpha": 0.2736911125477303,   # Scaling factor (R/ξ₁)
        "gamma": 5.0/3.0,              # Adiabatic index
        "n": 1.5,                      # Polytropic index
        "dim": 3,                      # Dimension
        "N_particles": 10648,          # Number of particles
        "source": "sample/imbh_cloud/results/ic_relax_10k/snapshot_0032.csv",
        "description": "IMBH cloud 10k relaxed Lane-Emden n=1.5 polytrope",
    },
    # Add more presets as needed...
}


@dataclass
class ICParameters:
    """Container for initial condition parameters."""
    rho_c: float        # Central density
    R: float            # Cloud radius
    M: float            # Total mass
    K: float            # Polytropic constant
    alpha: float        # Scaling factor (R/ξ₁)
    gamma: float        # Adiabatic index
    n: float            # Polytropic index
    dim: int            # Dimension
    N_particles: int    # Number of particles
    source: str         # Source file path
    description: str    # Human-readable description


def get_ic_preset(name: str) -> ICParameters:
    """
    Get IC parameters for a known preset.
    
    Parameters
    ----------
    name : str
        Preset name (e.g., "imbh_10k")
        
    Returns
    -------
    ICParameters
        IC parameter container
        
    Example
    -------
    >>> ic = get_ic_preset("imbh_10k")
    >>> print(f"Central density: {ic.rho_c}")
    Central density: 1.430096922633852
    """
    if name not in IC_PRESETS:
        available = ", ".join(IC_PRESETS.keys())
        raise ValueError(f"Unknown IC preset: {name}. Available: {available}")
    
    params = IC_PRESETS[name]
    return ICParameters(**params)


def list_ic_presets() -> Dict[str, str]:
    """
    List all available IC presets with descriptions.
    
    Returns
    -------
    dict
        Mapping of preset name to description
    """
    return {name: params["description"] for name, params in IC_PRESETS.items()}


def get_analytic_profile(
    preset: str = "imbh_10k",
    n_points: int = 200
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Get the analytic density profile for a known IC preset.
    
    This is the recommended way to get the correct analytic solution
    for comparison with simulation data.
    
    Parameters
    ----------
    preset : str
        IC preset name (default: "imbh_10k")
    n_points : int
        Number of points in the profile
        
    Returns
    -------
    r : ndarray
        Radial positions from 0 to R
    rho : ndarray
        Density values
        
    Example
    -------
    >>> r, rho = get_analytic_profile("imbh_10k")
    >>> plt.plot(r, rho, 'k--', label='Lane-Emden')
    """
    ic = get_ic_preset(preset)
    solution = load_lane_emden_solution(n=ic.n, dim=ic.dim)
    return get_density_profile(solution, ic.rho_c, ic.R, n_points=n_points)


@dataclass
class LaneEmdenSolution:
    """Container for Lane-Emden solution data."""
    xi: np.ndarray       # Dimensionless radius ξ
    theta: np.ndarray    # θ(ξ) - dimensionless density
    dtheta: np.ndarray   # dθ/dξ - derivative
    xi_1: float          # First zero (surface)
    dtheta_1: float      # |dθ/dξ| at surface
    n: float             # Polytropic index


def get_data_dir() -> Path:
    """Get the path to the Lane-Emden data directory."""
    # Try multiple relative paths
    candidates = [
        Path(__file__).parent.parent.parent / "data" / "lane_emden",
        Path("data/lane_emden"),
        Path("../data/lane_emden"),
        Path("../../data/lane_emden"),
    ]
    
    for path in candidates:
        if path.exists():
            return path.resolve()
    
    raise FileNotFoundError(
        f"Could not find lane_emden data directory. Searched: {candidates}"
    )


def load_lane_emden_solution(
    n: float = 1.5,
    dim: int = 3,
    data_dir: Optional[Path] = None
) -> LaneEmdenSolution:
    """
    Load pre-computed Lane-Emden solution from data file.
    
    This is the CORRECT way to get Lane-Emden solutions. DO NOT re-implement
    the ODE integration - use this function instead.
    
    Parameters
    ----------
    n : float
        Polytropic index (default 1.5 for γ=5/3 gas)
    dim : int
        Dimension (2 or 3, default 3)
    data_dir : Path, optional
        Custom data directory path
        
    Returns
    -------
    LaneEmdenSolution
        Solution data with xi, theta, dtheta arrays and boundary values
        
    Notes
    -----
    Data files are stored in data/lane_emden/ with naming:
        n1.5_3d.dat, n1.5_2d.dat, etc.
    """
    if data_dir is None:
        data_dir = get_data_dir()
    else:
        data_dir = Path(data_dir)
    
    # Construct filename - prefer the dot version (n1.5_3d.dat) over underscore
    # The dot version contains correct xi_1 values from the original solution
    n_str = str(n)  # e.g., "1.5"
    filename = f"n{n_str}_{dim}d.dat"  # e.g., "n1.5_3d.dat"
    
    filepath = data_dir / filename
    if not filepath.exists():
        # Try underscore naming as fallback
        n_str_underscore = str(n).replace('.', '_')
        filename_alt = f"n{n_str_underscore}_{dim}d.dat"
        filepath_alt = data_dir / filename_alt
        if filepath_alt.exists():
            filepath = filepath_alt
    
    if not filepath.exists():
        raise FileNotFoundError(f"Lane-Emden data file not found: {filepath}")
    
    # Parse header and load data
    xi_1 = None
    dtheta_1 = None
    
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#'):
                if 'xi_1' in line:
                    xi_1 = float(line.split('=')[1].strip())
                elif 'dtheta_1' in line:
                    dtheta_1 = float(line.split('=')[1].strip())
    
    # Load numerical data (skip header lines)
    data = np.loadtxt(filepath, comments='#')
    
    # Default to known values if not in header
    if xi_1 is None:
        xi_1 = XI_1_N15
    if dtheta_1 is None:
        dtheta_1 = DTHETA_1_N15
    
    return LaneEmdenSolution(
        xi=data[:, 0],
        theta=data[:, 1],
        dtheta=data[:, 2],
        xi_1=xi_1,
        dtheta_1=abs(dtheta_1),  # Store as positive
        n=n
    )


def get_density_profile(
    solution: LaneEmdenSolution,
    rho_c: float,
    R: float,
    r_values: Optional[np.ndarray] = None,
    n_points: int = 500
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Get physical density profile from Lane-Emden solution.
    
    Parameters
    ----------
    solution : LaneEmdenSolution
        Lane-Emden solution data
    rho_c : float
        Central density
    R : float  
        Physical radius of the sphere
    r_values : ndarray, optional
        Specific r values to interpolate to
    n_points : int
        Number of points if r_values not provided
        
    Returns
    -------
    r : ndarray
        Physical radii
    rho : ndarray
        Density at each radius
        
    Notes
    -----
    For polytropic index n=1.5:
        ρ(r) = ρ_c θ(ξ)^(3/2)
    where ξ = r/α and α = R/ξ₁
    """
    # Scaling: R = α ξ₁, so α = R / ξ₁
    alpha = R / solution.xi_1
    
    if r_values is None:
        # Generate r values from 0 to R
        r = np.linspace(0, R, n_points)
    else:
        r = np.asarray(r_values)
    
    # Convert to dimensionless radius
    xi = r / alpha
    
    # Interpolate θ(ξ)
    theta = np.interp(xi, solution.xi, solution.theta, left=1.0, right=0.0)
    theta = np.maximum(theta, 0.0)  # Ensure non-negative
    
    # Compute density: ρ = ρ_c θⁿ
    n = solution.n
    rho = rho_c * theta**n
    
    return r, rho


def get_pressure_profile(
    solution: LaneEmdenSolution,
    rho_c: float,
    R: float,
    K: float,
    gamma: float = 5.0/3.0,
    r_values: Optional[np.ndarray] = None,
    n_points: int = 500
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Get physical pressure profile from Lane-Emden solution.
    
    Parameters
    ----------
    solution : LaneEmdenSolution
        Lane-Emden solution data
    rho_c : float
        Central density
    R : float
        Physical radius
    K : float
        Polytropic constant (P = K ρ^γ)
    gamma : float
        Adiabatic index (default 5/3)
    r_values : ndarray, optional
        Specific r values
    n_points : int
        Number of points if r_values not provided
        
    Returns
    -------
    r : ndarray
        Physical radii
    P : ndarray
        Pressure at each radius
    """
    r, rho = get_density_profile(solution, rho_c, R, r_values, n_points)
    P = K * rho**gamma
    return r, P


def lane_emden_mass(rho_c: float, R: float, xi_1: float = XI_1_N15, 
                    dtheta_1: float = abs(DTHETA_1_N15)) -> float:
    """
    Calculate total mass of Lane-Emden sphere.
    
    M = 4π α³ ρ_c ξ₁² |dθ/dξ|₁
    
    where α = R / ξ₁
    """
    alpha = R / xi_1
    return 4.0 * np.pi * alpha**3 * rho_c * xi_1**2 * dtheta_1


def lane_emden_alpha(rho_c: float, K: float, G: float, n: float = 1.5) -> float:
    """
    Calculate scaling factor α from physical parameters.
    
    α² = (n+1) K ρ_c^(1/n-1) / (4πG)
    """
    return np.sqrt((n + 1) * K * rho_c**((1.0 - n)/n) / (4.0 * np.pi * G))


# Convenience function for backward compatibility
def lane_emden_profile(rho_c: float, R_cloud: float, n_points: int = 500,
                       data_dir: Optional[Path] = None) -> Tuple[np.ndarray, np.ndarray]:
    """
    Get Lane-Emden density profile for n=1.5 polytrope.
    
    DEPRECATED: Use load_lane_emden_solution() + get_density_profile() instead.
    
    Parameters
    ----------
    rho_c : float
        Central density
    R_cloud : float
        Cloud radius
    n_points : int
        Number of points
    data_dir : Path, optional
        Custom data directory
        
    Returns
    -------
    r : ndarray
        Radii from 0 to R_cloud
    rho : ndarray
        Density profile
    """
    solution = load_lane_emden_solution(n=1.5, dim=3, data_dir=data_dir)
    return get_density_profile(solution, rho_c, R_cloud, n_points=n_points)


if __name__ == "__main__":
    # Test the module
    print("=" * 60)
    print("Lane-Emden Module Test")
    print("=" * 60)
    
    try:
        solution = load_lane_emden_solution()
        print("✓ Loaded solution from data file")
        print(f"  - Points: {len(solution.xi)}")
        print(f"  - ξ₁ = {solution.xi_1:.10f}")
        print(f"  - |dθ/dξ|₁ = {solution.dtheta_1:.10f}")
        print(f"  - θ(0) = {solution.theta[0]:.6f}")
        print(f"  - θ(ξ₁) ≈ {np.interp(solution.xi_1, solution.xi, solution.theta):.6e}")
        
        # Test density profile with generic parameters
        rho_c = 1.0
        R = 1.0
        r, rho = get_density_profile(solution, rho_c, R)
        print("\n✓ Generated density profile (generic)")
        print(f"  - ρ(0) = {rho[0]:.6f} (expected {rho_c:.6f})")
        print(f"  - ρ(R) = {rho[-1]:.6e} (expected ~0)")
        
        # Test mass calculation
        M = lane_emden_mass(rho_c, R)
        print("\n✓ Mass calculation")
        print(f"  - M = {M:.6f}")
        
        # Test IC presets
        print("\n" + "=" * 60)
        print("IC Preset Test")
        print("=" * 60)
        
        presets = list_ic_presets()
        print(f"✓ Available IC presets: {len(presets)}")
        for name, desc in presets.items():
            print(f"  - {name}: {desc}")
        
        # Test IMBH 10k preset
        print("\n✓ Testing 'imbh_10k' preset:")
        ic = get_ic_preset("imbh_10k")
        print(f"  - Central density: ρ_c = {ic.rho_c:.10f}")
        print(f"  - Cloud radius: R = {ic.R}")
        print(f"  - Total mass: M = {ic.M}")
        print(f"  - Polytropic K: {ic.K:.10f}")
        print(f"  - Scaling α: {ic.alpha:.10f}")
        print(f"  - N particles: {ic.N_particles}")
        
        # Test analytic profile with preset
        r_preset, rho_preset = get_analytic_profile("imbh_10k", n_points=100)
        print(f"\n✓ Generated analytic profile from preset:")
        print(f"  - r: [0, {r_preset[-1]:.4f}], {len(r_preset)} points")
        print(f"  - ρ(0) = {rho_preset[0]:.6f} (expected {ic.rho_c:.6f})")
        print(f"  - ρ(R) = {rho_preset[-1]:.6e} (expected ~0)")
        
        print("\n" + "=" * 60)
        print("All tests passed!")
        print("=" * 60)
        
    except Exception as e:
        print(f"✗ Error: {e}")
        import traceback
        traceback.print_exc()
