#!/usr/bin/env python3
"""
Single Source of Truth for SR Shock Tube Test Cases
====================================================

Kitajima et al. (2025) test case definitions for relativistic shock tubes.
All scripts should import from this module to ensure consistency.
"""

import os
import sys
import glob
import json
from pathlib import Path

# Add docs directory to path for relativistic_riemann_solver
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'docs'))
from relativistic_riemann_solver import RelativisiticRiemannSolver

# Kitajima et al. (2025) test case definitions
# Each test case has: (P, n, v) for left and right states
TEST_CASES = {
    'sod': {
        'name': 'Sod Problem (Section 3.1.1)',
        'left': (1.0, 1.0, 0.0),      # (P, n, v)
        'right': (0.1, 0.125, 0.0),
    },
    'blast_wave': {
        'name': 'Standard Blast Wave (Section 3.1.2)',
        'left': (40.0/3.0, 10.0, 0.0),
        'right': (1.0e-6, 1.0, 0.0),
    },
    'strong_blast': {
        'name': 'Strong Blast Wave (Section 3.1.3)',
        'left': (1000.0, 1.0, 0.0),
        'right': (0.01, 1.0, 0.0),
    },
    'ultra_relat': {
        'name': 'Ultra-Relativistic Shock (Section 3.2)',
        'left': (1.0, 1.0, 0.9),  # v=0.9c default, can be overridden
        'right': (1.0, 1.0, 0.0),
    },
}


def detect_test_type(results_dir):
    """
    Detect test type from directory name or config file.
    
    Args:
        results_dir: Path to results directory
        
    Returns:
        tuple: (test_type, v_left) where test_type is a key in TEST_CASES
    """
    results_dir = Path(results_dir)
    
    # Try to load config file if it exists
    config_patterns = [
        results_dir / '..' / '*.json',
        results_dir / '..' / 'config' / '*.json',
        results_dir / '..' / 'config' / 'presets' / '*.json',
    ]
    
    for pattern in config_patterns:
        configs = glob.glob(str(pattern))
        for config_file in configs:
            try:
                with open(config_file, 'r') as f:
                    config = json.load(f)
                    if 'testType' in config:
                        return config['testType'], config.get('v_left', 0.9)
            except:
                pass
    
    # Detect from directory name
    dir_name = results_dir.name.lower()
    parent_name = results_dir.parent.name.lower() if results_dir.parent else ''
    dir_lower = f"{parent_name}/{dir_name}"
    
    if 'strong' in dir_lower:
        return 'strong_blast', 0.0
    elif 'blast' in dir_lower:
        return 'blast_wave', 0.0
    elif 'ultra' in dir_lower:
        return 'ultra_relat', 0.9
    elif 'sod' in dir_lower:
        return 'sod', 0.0
    else:
        return 'sod', 0.0


def get_initial_conditions(test_type, v_left_override=None):
    """
    Get initial conditions for a test case.
    
    Args:
        test_type: Key in TEST_CASES
        v_left_override: Override left velocity (for ultra_relat tests)
        
    Returns:
        dict with keys: pL, rhoL, vL, pR, rhoR, vR, name
    """
    if test_type not in TEST_CASES:
        print(f"Warning: Unknown test type '{test_type}', using 'sod'")
        test_type = 'sod'
    
    test = TEST_CASES[test_type]
    pL, rhoL, vL = test['left']
    pR, rhoR, vR = test['right']
    
    # Override velocity for ultra-relativistic test
    if test_type == 'ultra_relat' and v_left_override is not None:
        vL = v_left_override
    
    return {
        'pL': pL, 'rhoL': rhoL, 'vL': vL,
        'pR': pR, 'rhoR': rhoR, 'vR': vR,
        'name': test['name']
    }


def get_exact_solution(t, test_type='sod', v_left_override=None, gamma=5.0/3.0, 
                       n_points=500, x_min=-0.5, x_max=0.5, x0=0.0):
    """
    Compute exact Riemann solution for any Kitajima test case at time t.
    
    Args:
        t: time
        test_type: 'sod', 'blast_wave', 'strong_blast', or 'ultra_relat'
        v_left_override: override left velocity (for ultra_relat tests)
        gamma: adiabatic index (default 5/3)
        n_points: number of points for solution
        x_min, x_max: domain bounds
        x0: initial discontinuity location
        
    Returns:
        dict with keys: x, pres, dens, vel, p_star, v_star
    """
    ic = get_initial_conditions(test_type, v_left_override)
    pL, rhoL, vL = ic['pL'], ic['rhoL'], ic['vL']
    pR, rhoR, vR = ic['pR'], ic['rhoR'], ic['vR']
    
    x = np.linspace(x_min, x_max, n_points)
    
    if t <= 0:
        # Initial condition
        pres = np.where(x < x0, pL, pR)
        dens = np.where(x < x0, rhoL, rhoR)
        vel = np.where(x < x0, vL, vR)
        return {
            'x': x, 'pres': pres, 'dens': dens, 'vel': vel,
            'p_star': 0.5*(pL+pR), 'v_star': 0.0
        }
    
    solver = RelativisiticRiemannSolver(gamma)
    solver.set_initial_states(pL, rhoL, vL, pR, rhoR, vR)
    
    # Solve - returns solution on [0, 1] with discontinuity at 0.5
    x_raw, pres_raw, dens_raw, vel_raw, _ = solver.solve(t, n=n_points)
    
    # Convert from [0, 1] domain to requested domain
    # The solver puts the initial discontinuity at x=0.5
    x_shifted = x_raw - 0.5 + x0
    
    # Scale to requested domain
    scale = (x_max - x_min)
    x_scaled = x_min + (x_shifted - x_min) * scale / (1.0)
    
    # Interpolate to uniform grid
    pres = np.interp(x, x_shifted, pres_raw)
    dens = np.interp(x, x_shifted, dens_raw)
    vel = np.interp(x, x_shifted, vel_raw)
    
    # Get star state values (from the plateau region)
    mid_idx = len(pres_raw) // 2
    p_star = pres_raw[mid_idx]
    v_star = vel_raw[mid_idx]
    
    return {
        'x': x, 'pres': pres, 'dens': dens, 'vel': vel,
        'p_star': p_star, 'v_star': v_star
    }


# Import numpy here to avoid issues if module is imported before numpy
import numpy as np
