#!/usr/bin/env python3
"""
Shock Tracker for BNS Cocoon Simulations

Tracks the shock position, velocity, and Lorentz factor as a function of time
by analyzing simulation snapshots.

Methods:
1. Density jump detection
2. Velocity gradient detection
3. Pressure gradient detection
"""

import argparse
import glob
from pathlib import Path
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass
class ShockState:
    """Shock state at a single time."""
    time: float              # Simulation time
    r_shock: float           # Shock radius
    v_shock: float           # Shock velocity [c]
    gamma_shock: float       # Shock Lorentz factor
    rho_upstream: float      # Upstream (unshocked) density
    rho_downstream: float    # Downstream (shocked) density
    compression_ratio: float # ρ_down / ρ_up
    

class ShockTracker:
    """
    Track shock position and properties in simulation snapshots.
    
    Usage:
        tracker = ShockTracker(results_dir)
        history = tracker.track_all()
        tracker.save_history(history, 'shock_history.csv')
    """
    
    def __init__(self, results_dir: str, 
                 dimension: int = 1,
                 density_threshold: float = 2.0):
        """
        Initialize shock tracker.
        
        Args:
            results_dir: Directory containing snapshot files
            dimension: Simulation dimension (1 or 2)
            density_threshold: Minimum compression ratio for shock detection
        """
        self.results_dir = Path(results_dir)
        self.dimension = dimension
        self.density_threshold = density_threshold
        
        self.snapshots = self._load_snapshot_list()
    
    def _load_snapshot_list(self) -> List[Path]:
        """Find all snapshot files."""
        pattern = str(self.results_dir / "snapshot_*.csv")
        files = sorted(glob.glob(pattern))
        return [Path(f) for f in files]
    
    def load_snapshot(self, filepath: Path) -> Tuple[pd.DataFrame, float]:
        """
        Load a snapshot and extract time.
        
        Returns:
            DataFrame with particle data, simulation time
        """
        # Read time from header
        time = 0.0
        with open(filepath) as f:
            for line in f:
                if line.startswith('# Time (code):'):
                    time = float(line.split(':')[1].strip())
                    break
                elif line.startswith('# time='):
                    time = float(line.split('=')[1].strip())
                    break
        
        # Read particle data
        df = pd.read_csv(filepath, comment='#')
        
        return df, time
    
    def compute_radial_profile(self, df: pd.DataFrame, 
                               n_bins: int = 100) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Compute radial density and velocity profiles.
        
        Returns:
            r: Radial bin centers
            rho: Mean density in each bin
            v_r: Mean radial velocity in each bin
        """
        # Compute radial position
        if 'pos_x' in df.columns:
            x, y = df['pos_x'].values, df['pos_y'].values
            z = df['pos_z'].values if 'pos_z' in df.columns else np.zeros_like(x)
        else:
            x, y = df['x'].values, df['y'].values
            z = df['z'].values if 'z' in df.columns else np.zeros_like(x)
        
        r = np.sqrt(x**2 + y**2 + z**2)
        
        # Compute radial velocity
        if 'vel_x' in df.columns:
            vx, vy = df['vel_x'].values, df['vel_y'].values
            vz = df['vel_z'].values if 'vel_z' in df.columns else np.zeros_like(vx)
        else:
            vx, vy = df['vx'].values, df['vy'].values
            vz = df['vz'].values if 'vz' in df.columns else np.zeros_like(vx)
        
        v_r = (x*vx + y*vy + z*vz) / np.maximum(r, 1e-10)
        
        # Density
        rho = df['dens'].values if 'dens' in df.columns else df['density'].values
        
        # Bin data
        r_edges = np.linspace(0, r.max() * 1.1, n_bins + 1)
        r_centers = 0.5 * (r_edges[:-1] + r_edges[1:])
        
        rho_binned = np.zeros(n_bins)
        vr_binned = np.zeros(n_bins)
        
        for i in range(n_bins):
            mask = (r >= r_edges[i]) & (r < r_edges[i+1])
            if np.sum(mask) > 0:
                rho_binned[i] = np.mean(rho[mask])
                vr_binned[i] = np.mean(v_r[mask])
        
        return r_centers, rho_binned, vr_binned
    
    def detect_shock(self, r: np.ndarray, rho: np.ndarray, 
                     v_r: np.ndarray) -> Optional[Tuple[float, float, float]]:
        """
        Detect shock position from density profile.
        
        Returns:
            (r_shock, rho_upstream, rho_downstream) or None if no shock found
        """
        # Look for maximum density gradient (shock front)
        drho_dr = np.gradient(rho, r)
        
        # Find peaks in density gradient
        i_max = np.argmax(np.abs(drho_dr))
        
        # Check if it's a real shock (compression ratio > threshold)
        if i_max > 0 and i_max < len(rho) - 1:
            rho_up = np.mean(rho[i_max+1:min(i_max+5, len(rho))])
            rho_down = np.mean(rho[max(0, i_max-4):i_max])
            
            if rho_up > 0 and rho_down / rho_up > self.density_threshold:
                return r[i_max], rho_up, rho_down
        
        return None
    
    def track_snapshot(self, filepath: Path) -> Optional[ShockState]:
        """
        Track shock in a single snapshot.
        
        Returns:
            ShockState or None if no shock detected
        """
        df, time = self.load_snapshot(filepath)
        r, rho, v_r = self.compute_radial_profile(df)
        
        result = self.detect_shock(r, rho, v_r)
        if result is None:
            return None
        
        r_shock, rho_up, rho_down = result
        
        # Estimate shock velocity from velocity profile
        idx = np.argmin(np.abs(r - r_shock))
        v_shock = v_r[idx]
        
        # Lorentz factor
        gamma_shock = 1.0 / np.sqrt(1.0 - min(v_shock**2, 0.9999))
        
        return ShockState(
            time=time,
            r_shock=r_shock,
            v_shock=v_shock,
            gamma_shock=gamma_shock,
            rho_upstream=rho_up,
            rho_downstream=rho_down,
            compression_ratio=rho_down / rho_up if rho_up > 0 else 0,
        )
    
    def track_all(self) -> List[ShockState]:
        """
        Track shock through all snapshots.
        
        Returns:
            List of ShockState objects
        """
        history = []
        
        for snap in self.snapshots:
            state = self.track_snapshot(snap)
            if state is not None:
                history.append(state)
                print(f"t={state.time:.3f}: r_sh={state.r_shock:.3f}, "
                      f"v_sh={state.v_shock:.3f}c, Γ={state.gamma_shock:.2f}")
        
        return history
    
    def save_history(self, history: List[ShockState], filepath: str):
        """Save shock history to CSV file."""
        data = {
            'time': [s.time for s in history],
            'r_shock': [s.r_shock for s in history],
            'v_shock': [s.v_shock for s in history],
            'gamma_shock': [s.gamma_shock for s in history],
            'rho_upstream': [s.rho_upstream for s in history],
            'rho_downstream': [s.rho_downstream for s in history],
            'compression_ratio': [s.compression_ratio for s in history],
        }
        
        df = pd.DataFrame(data)
        df.to_csv(filepath, index=False)
        print(f"Saved shock history to {filepath}")
    
    def estimate_shock_velocity(self, history: List[ShockState]) -> np.ndarray:
        """
        Estimate shock velocity from position vs time.
        
        Uses finite differences of r_shock(t).
        """
        times = np.array([s.time for s in history])
        radii = np.array([s.r_shock for s in history])
        
        # Central difference
        v_shock = np.gradient(radii, times)
        
        return v_shock


def main():
    parser = argparse.ArgumentParser(description='Track shock in simulation')
    parser.add_argument('results_dir', help='Directory containing snapshots')
    parser.add_argument('--output', '-o', default='shock_history.csv',
                       help='Output CSV file')
    parser.add_argument('--dimension', '-d', type=int, default=1,
                       help='Simulation dimension')
    parser.add_argument('--threshold', type=float, default=2.0,
                       help='Density compression threshold')
    
    args = parser.parse_args()
    
    tracker = ShockTracker(
        args.results_dir,
        dimension=args.dimension,
        density_threshold=args.threshold,
    )
    
    print(f"Found {len(tracker.snapshots)} snapshots")
    print("Tracking shock...")
    
    history = tracker.track_all()
    
    if history:
        output_path = Path(args.results_dir) / args.output
        tracker.save_history(history, str(output_path))
        
        # Summary
        print(f"\n=== Shock Tracking Summary ===")
        print(f"Time range: {history[0].time:.3f} - {history[-1].time:.3f}")
        print(f"Initial shock radius: {history[0].r_shock:.3f}")
        print(f"Final shock radius: {history[-1].r_shock:.3f}")
        print(f"Max Lorentz factor: {max(s.gamma_shock for s in history):.2f}")
    else:
        print("No shock detected in any snapshot")


if __name__ == '__main__':
    main()
