#!/usr/bin/env python3
"""
Lane-Emden equation solution module - Single Source of Truth (SSOT)

This module provides the CORRECT Lane-Emden solution loaded from pre-computed
data files. DO NOT re-implement ODE integration - use this module instead.

Lane-Emden equation:
    (1/ξ²) d/dξ (ξ² dθ/dξ) = -θⁿ

For n=1.5 (γ=5/3 polytrope):
    - First zero: ξ₁ = 3.6537540101
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
from typing import Optional, Tuple
from dataclasses import dataclass


# Constants from exact Lane-Emden n=1.5 solution
XI_1_N15 = 3.6537540101  # First zero of θ for n=1.5
DTHETA_1_N15 = -0.20331244769  # dθ/dξ at ξ₁ for n=1.5


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
        
        # Test density profile
        rho_c = 1.0
        R = 1.0
        r, rho = get_density_profile(solution, rho_c, R)
        print("\n✓ Generated density profile")
        print(f"  - ρ(0) = {rho[0]:.6f} (expected {rho_c:.6f})")
        print(f"  - ρ(R) = {rho[-1]:.6e} (expected ~0)")
        
        # Test mass calculation
        M = lane_emden_mass(rho_c, R)
        print("\n✓ Mass calculation")
        print(f"  - M = {M:.6f}")
        
    except Exception as e:
        print(f"✗ Error: {e}")
        import traceback
        traceback.print_exc()
