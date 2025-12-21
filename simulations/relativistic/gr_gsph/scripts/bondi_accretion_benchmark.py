#!/usr/bin/env python3
"""
Bondi Accretion Benchmark Suite

Reproduces the Bondi accretion tests from Liptai & Price (2019) MNRAS 485, 819:

1. Geodesic Bondi Flow (Fig 14)
   - Marginally bound, pressure-less radial infall
   - Analytic solution: v^r = sqrt(2M/r) * (1 - 2M/r)
                       rho = rho_0 / r^2 * sqrt(r / 2M)
                       u = u_0 * (r^2 * sqrt(2M/r))^(1-gamma)

2. Sonic Point Bondi Flow (Fig 18-19)
   - Includes pressure effects
   - Flow passes through sonic point at r_c = 8M
   - Analytic solution via implicit equation for T(r)

Reference: Hawley, Smarr & Wilson (1984)

Usage:
    python bondi_accretion_benchmark.py [--plot-only]
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import brentq
import argparse


# =============================================================================
# Analytic Solutions for Bondi Accretion
# =============================================================================

class BondiGeodesic:
    """
    Geodesic (pressure-less) Bondi flow in Schwarzschild spacetime.

    From Hawley, Smarr & Wilson (1984), valid when u << 1.
    """

    def __init__(self, M=1.0, rho_0=1.0, u_0=1e-9, gamma=5.0/3.0):
        self.M = M
        self.rho_0 = rho_0
        self.u_0 = u_0
        self.gamma = gamma

    def v_r(self, r):
        """Radial coordinate velocity (infall, negative)."""
        M = self.M
        return -np.sqrt(2*M/r) * (1 - 2*M/r)

    def rho(self, r):
        """Rest-mass density."""
        M = self.M
        return self.rho_0 / r**2 * np.sqrt(r / (2*M))

    def rho_star(self, r):
        """Conserved density rho* = rho * W * sqrt(gamma)."""
        v = self.v_r(r)
        W = 1.0 / np.sqrt(1 - v**2)  # Lorentz factor
        # In Schwarzschild, sqrt(gamma) = 1
        return self.rho(r) * W

    def u(self, r):
        """Specific internal energy."""
        M = self.M
        return self.u_0 * (r**2 * np.sqrt(2*M/r))**(1 - self.gamma)


class BondiSonicPoint:
    """
    Sonic point Bondi flow in Schwarzschild spacetime.

    From Michel (1972), Hawley, Smarr & Wilson (1984).
    Flow passes through critical point at r_c.
    """

    def __init__(self, M=1.0, r_c=8.0, K_0=1.0, gamma=5.0/3.0):
        self.M = M
        self.r_c = r_c  # Critical (sonic) point
        self.K_0 = K_0  # Entropy normalization
        self.gamma = gamma
        self.n = 1.0 / (gamma - 1.0)  # Polytropic index

        # Compute constants from critical point
        self._compute_constants()

    def _compute_constants(self):
        """Compute C_1 and C_2 from critical point conditions."""
        M, r_c, n = self.M, self.r_c, self.n

        # Critical point values
        U_c = np.sqrt(M / (2 * r_c))
        v_c = U_c / np.sqrt(1 - 3 * U_c**2)
        T_c = v_c**2 * n / ((1 + n) * (1 - n * v_c**2))

        self.U_c = U_c
        self.v_c = v_c
        self.T_c = T_c

        # Constants
        self.C_1 = U_c * r_c**2 * T_c**n
        self.C_2 = (1 + (1 + n) * T_c)**2 * (1 - 2*M/r_c + self.C_1**2 / (r_c**4 * T_c**(2*n)))

    def _solve_T(self, r):
        """Solve implicit equation for temperature T(r)."""
        M, n, C_1, C_2 = self.M, self.n, self.C_1, self.C_2

        def equation(T):
            if T <= 0:
                return 1e10
            term1 = (1 + (n + 1) * T)**2
            term2 = 1 - 2*M/r + C_1**2 / (r**4 * T**(2*n))
            return term1 * term2 - C_2

        # For inflow, we need the lower branch
        try:
            # Bracket for inflow solution
            T = brentq(equation, 1e-10, self.T_c * 2)
        except:
            T = self.T_c
        return T

    def T(self, r):
        """Temperature T = P/rho = (gamma-1)*u at radius r."""
        if np.isscalar(r):
            return self._solve_T(r)
        return np.array([self._solve_T(ri) for ri in r])

    def U_r(self, r):
        """Four-velocity radial component U^r."""
        T = self.T(r)
        return self.C_1 / (r**2 * T**self.n)

    def v_r(self, r):
        """Coordinate velocity v^r (negative for infall)."""
        M = self.M
        U_r = self.U_r(r)
        U_0 = np.sqrt(1 - 2*M/r + U_r**2) / (1 - 2*M/r)
        return -U_r / U_0  # Negative for infall

    def rho(self, r):
        """Rest-mass density."""
        T = self.T(r)
        return self.K_0 * T**self.n

    def rho_star(self, r):
        """Conserved density."""
        v = self.v_r(r)
        W = 1.0 / np.sqrt(np.maximum(1 - v**2, 1e-10))
        return self.rho(r) * W

    def u(self, r):
        """Specific internal energy."""
        T = self.T(r)
        return self.n * T


# =============================================================================
# Plotting Functions
# =============================================================================

def plot_bondi_geodesic(output_dir, sim_data=None):
    """
    Plot Bondi geodesic flow (Fig 14).

    If sim_data is provided, overlay C++ simulation results.
    """
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    M = 1.0
    bondi = BondiGeodesic(M=M, rho_0=1.0, u_0=1e-9)

    # Radial range
    r = np.linspace(2.5*M, 18*M, 200)

    # Analytic solutions (use rest-frame density, not conserved density)
    rho_ana = bondi.rho(r)  # Rest-frame density
    v_r_ana = bondi.v_r(r)
    u_ana = bondi.u(r)

    # Plot density (rest-frame)
    ax = axes[0]
    ax.plot(r, rho_ana, 'g-', lw=2.5, label='Analytic Solution')
    if sim_data is not None:
        # Use rest-frame density from simulation (dens column)
        ax.scatter(sim_data['r'], sim_data['rho'], s=5, c='blue',
                   alpha=0.7, label='C++ GR-GSPH')
    ax.set_xlabel('r/M', fontsize=12)
    ax.set_ylabel(r'$\rho$ (rest-frame)', fontsize=12)
    ax.set_title('Rest-Frame Density', fontweight='bold')
    ax.legend(fontsize=10, framealpha=0.9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(2, 20)

    # Plot velocity
    ax = axes[1]
    ax.plot(r, -v_r_ana, 'g-', lw=2.5, label='Analytic Solution')
    if sim_data is not None:
        ax.scatter(sim_data['r'], -sim_data['v_r'], s=5, c='blue',
                   alpha=0.7, label='C++ GR-GSPH')
    ax.set_xlabel('r/M', fontsize=12)
    ax.set_ylabel(r'$|v^r|$', fontsize=12)
    ax.set_title('Infall Velocity', fontweight='bold')
    ax.legend(fontsize=10, framealpha=0.9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(2, 20)

    # Plot thermal energy
    ax = axes[2]
    ax.semilogy(r, u_ana, 'g-', lw=2.5, label='Analytic Solution')
    if sim_data is not None:
        ax.scatter(sim_data['r'], sim_data['u'], s=5, c='blue',
                   alpha=0.7, label='C++ GR-GSPH')
    ax.set_xlabel('r/M', fontsize=12)
    ax.set_ylabel(r'$u$ (specific internal energy)', fontsize=12)
    ax.set_title('Thermal Energy', fontweight='bold')
    ax.legend(fontsize=10, framealpha=0.9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(2, 20)

    fig.suptitle('Bondi Geodesic Accretion: Liptai & Price (2019) Fig 14\n'
                 'Marginally bound, pressure-less radial infall (M=1)',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()

    output_file = output_dir / 'bondi_geodesic.png'
    fig.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()


def plot_bondi_sonic(output_dir, sim_data=None):
    """
    Plot Bondi sonic point flow (Fig 18).
    """
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    M = 1.0
    r_c = 8.0 * M
    bondi = BondiSonicPoint(M=M, r_c=r_c, K_0=1.0)

    # Radial range
    r = np.linspace(2.5*M, 18*M, 100)

    # Analytic solutions
    rho_star_ana = bondi.rho_star(r)
    v_r_ana = bondi.v_r(r)
    u_ana = bondi.u(r)

    # Plot density
    ax = axes[0]
    ax.plot(r, rho_star_ana, 'g-', lw=2.5, label='Analytic Solution')
    if sim_data is not None:
        ax.scatter(sim_data['r'], sim_data['rho_star'], s=5, c='blue',
                   alpha=0.7, label='C++ GR-GSPH')
    ax.axvline(r_c, color='r', ls='--', lw=1.5, alpha=0.7, label=f'Sonic point r={r_c:.0f}M')
    ax.set_xlabel('r/M', fontsize=12)
    ax.set_ylabel(r'$\rho^*$', fontsize=12)
    ax.set_title('Conserved Density', fontweight='bold')
    ax.legend(fontsize=10, framealpha=0.9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(2, 20)

    # Plot velocity
    ax = axes[1]
    ax.plot(r, -v_r_ana, 'g-', lw=2.5, label='Analytic Solution')
    if sim_data is not None:
        ax.scatter(sim_data['r'], -sim_data['v_r'], s=5, c='blue',
                   alpha=0.7, label='C++ GR-GSPH')
    ax.axvline(r_c, color='r', ls='--', lw=1.5, alpha=0.7)
    ax.set_xlabel('r/M', fontsize=12)
    ax.set_ylabel(r'$|v^r|$', fontsize=12)
    ax.set_title('Infall Velocity', fontweight='bold')
    ax.legend(fontsize=10, framealpha=0.9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(2, 20)

    # Plot thermal energy
    ax = axes[2]
    ax.plot(r, u_ana, 'g-', lw=2.5, label='Analytic Solution')
    if sim_data is not None:
        ax.scatter(sim_data['r'], sim_data['u'], s=5, c='blue',
                   alpha=0.7, label='C++ GR-GSPH')
    ax.axvline(r_c, color='r', ls='--', lw=1.5, alpha=0.7)
    ax.set_xlabel('r/M', fontsize=12)
    ax.set_ylabel(r'$u$ (specific internal energy)', fontsize=12)
    ax.set_title('Thermal Energy', fontweight='bold')
    ax.legend(fontsize=10, framealpha=0.9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(2, 20)

    fig.suptitle('Bondi Sonic Point Accretion: Liptai & Price (2019) Fig 18\n'
                 f'Transonic flow with sonic point at r_c = {r_c:.0f}M (M=1, K_0=1)',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()

    output_file = output_dir / 'bondi_sonic_point.png'
    fig.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()


def plot_combined_bondi(output_dir, geodesic_data=None, sonic_data=None):
    """Create combined figure showing both Bondi flows with C++ simulation overlay."""

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    M = 1.0
    r = np.linspace(2.5*M, 18*M, 200)

    # Top row: Geodesic flow
    bondi_geo = BondiGeodesic(M=M, rho_0=1.0, u_0=1e-9)

    ax = axes[0, 0]
    ax.plot(r, bondi_geo.rho(r), 'g-', lw=2.5, label='Analytic Solution')
    if geodesic_data is not None:
        ax.scatter(geodesic_data['r'], geodesic_data['rho'], s=5, c='blue',
                   alpha=0.7, label='C++ GR-GSPH')
    ax.set_xlabel('r/M')
    ax.set_ylabel(r'$\rho$ (rest-frame)')
    ax.set_title('Geodesic: Density', fontweight='bold')
    ax.legend(fontsize=9, framealpha=0.9)
    ax.grid(True, alpha=0.3)
    ax.text(0.02, 0.98, 'GEODESIC FLOW\n(pressure-less)', transform=ax.transAxes,
            fontsize=9, fontweight='bold', va='top', color='darkblue')

    ax = axes[0, 1]
    ax.plot(r, -bondi_geo.v_r(r), 'g-', lw=2.5, label='Analytic Solution')
    if geodesic_data is not None:
        ax.scatter(geodesic_data['r'], -geodesic_data['v_r'], s=5, c='blue',
                   alpha=0.7, label='C++ GR-GSPH')
    ax.set_xlabel('r/M')
    ax.set_ylabel(r'$|v^r|$')
    ax.set_title('Geodesic: Velocity', fontweight='bold')
    ax.legend(fontsize=9, framealpha=0.9)
    ax.grid(True, alpha=0.3)

    ax = axes[0, 2]
    ax.semilogy(r, bondi_geo.u(r), 'g-', lw=2.5, label='Analytic Solution')
    if geodesic_data is not None and 'u' in geodesic_data:
        ax.scatter(geodesic_data['r'], geodesic_data['u'], s=5, c='blue',
                   alpha=0.7, label='C++ GR-GSPH')
    ax.set_xlabel('r/M')
    ax.set_ylabel(r'$u$')
    ax.set_title('Geodesic: Thermal Energy', fontweight='bold')
    ax.legend(fontsize=9, framealpha=0.9)
    ax.grid(True, alpha=0.3)

    # Bottom row: Sonic point flow
    r_c = 8.0 * M
    r_sonic = np.linspace(2.5*M, 18*M, 100)
    bondi_sonic = BondiSonicPoint(M=M, r_c=r_c, K_0=1.0)

    ax = axes[1, 0]
    ax.plot(r_sonic, bondi_sonic.rho_star(r_sonic), 'g-', lw=2.5, label='Analytic Solution')
    if sonic_data is not None:
        ax.scatter(sonic_data['r'], sonic_data['rho_star'], s=5, c='blue',
                   alpha=0.7, label='C++ GR-GSPH')
    ax.axvline(r_c, color='r', ls='--', lw=1.5, alpha=0.7, label=f'Sonic r={r_c:.0f}M')
    ax.set_xlabel('r/M')
    ax.set_ylabel(r'$\rho^*$')
    ax.set_title('Sonic Point: Density', fontweight='bold')
    ax.legend(fontsize=9, framealpha=0.9)
    ax.grid(True, alpha=0.3)
    ax.text(0.02, 0.98, 'SONIC POINT FLOW\n(with pressure)', transform=ax.transAxes,
            fontsize=9, fontweight='bold', va='top', color='darkgreen')

    ax = axes[1, 1]
    ax.plot(r_sonic, -bondi_sonic.v_r(r_sonic), 'g-', lw=2.5, label='Analytic Solution')
    if sonic_data is not None:
        ax.scatter(sonic_data['r'], -sonic_data['v_r'], s=5, c='blue',
                   alpha=0.7, label='C++ GR-GSPH')
    ax.axvline(r_c, color='r', ls='--', lw=1.5, alpha=0.7)
    ax.set_xlabel('r/M')
    ax.set_ylabel(r'$|v^r|$')
    ax.set_title('Sonic Point: Velocity', fontweight='bold')
    ax.legend(fontsize=9, framealpha=0.9)
    ax.grid(True, alpha=0.3)

    ax = axes[1, 2]
    ax.plot(r_sonic, bondi_sonic.u(r_sonic), 'g-', lw=2.5, label='Analytic Solution')
    if sonic_data is not None and 'u' in sonic_data:
        ax.scatter(sonic_data['r'], sonic_data['u'], s=5, c='blue',
                   alpha=0.7, label='C++ GR-GSPH')
    ax.axvline(r_c, color='r', ls='--', lw=1.5, alpha=0.7)
    ax.set_xlabel('r/M')
    ax.set_ylabel(r'$u$')
    ax.set_title('Sonic Point: Thermal Energy', fontweight='bold')
    ax.legend(fontsize=9, framealpha=0.9)
    ax.grid(True, alpha=0.3)

    title = 'Bondi Accretion in Schwarzschild Spacetime\nLiptai & Price (2019) Figs 14 & 18'
    if geodesic_data is not None:
        title += ' - Blue: C++ GR-GSPH, Green: Analytic'
    fig.suptitle(title, fontsize=14, fontweight='bold')
    plt.tight_layout()

    output_file = output_dir / 'bondi_accretion_combined.png'
    fig.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()


def load_simulation_data(sim_dir, snapshot_num=0):
    """
    Load C++ simulation snapshot.

    Returns dict with 'r', 'rho_star', 'v_r', 'u' arrays, or None if not found.
    """
    import glob
    import pandas as pd

    if not sim_dir.exists():
        print(f"  Simulation directory not found: {sim_dir}")
        return None

    # Find snapshot file
    snapshot_file = sim_dir / f'snapshot_{snapshot_num:04d}.csv'
    if not snapshot_file.exists():
        # Try to find any snapshot
        files = sorted(glob.glob(str(sim_dir / 'snapshot_*.csv')))
        if not files:
            print(f"  No snapshot files found in {sim_dir}")
            return None
        snapshot_file = Path(files[0])

    print(f"  Loading: {snapshot_file.name}")

    # Load data
    data = pd.read_csv(snapshot_file, comment='#')

    # Filter out ghost particles
    if 'is_ghost' in data.columns:
        mask = data['is_ghost'] == 0
        data = data[mask]

    # Extract arrays
    result = {
        'r': data['pos_x'].values,  # radial position
    }

    # Rest-frame density
    # Note: For GR-GSPH, 'dens' is the SPH-estimated rest-frame density
    if 'dens' in data.columns:
        result['rho'] = data['dens'].values
        # Compute Lorentz factor from velocity for rho_star if needed
        if 'vel_x' in data.columns:
            v = data['vel_x'].values
            W = 1.0 / np.sqrt(1 - v**2)
            result['rho_star'] = result['rho'] * W
        else:
            result['rho_star'] = result['rho']

    if 'vel_x' in data.columns:
        result['v_r'] = data['vel_x'].values

    if 'ene' in data.columns:
        result['u'] = data['ene'].values

    return result


def main():
    parser = argparse.ArgumentParser(description='Bondi Accretion Benchmark')
    parser.add_argument('--plot-only', action='store_true',
                        help='Only generate analytic solution plots')
    args = parser.parse_args()

    script_dir = Path(__file__).parent
    output_dir = script_dir.parent / 'results'
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("  Bondi Accretion Benchmark Suite")
    print("  Based on Liptai & Price (2019) MNRAS 485, 819")
    print("=" * 70)

    # Try to load C++ simulation data
    project_root = script_dir.parent.parent.parent.parent
    geodesic_dir = project_root / 'simulations' / 'relativistic' / 'gr_gsph' / 'results' / 'bondi_geodesic'
    sonic_dir = project_root / 'simulations' / 'relativistic' / 'gr_gsph' / 'results' / 'bondi_sonic'

    print("\n[0] Loading C++ Simulation Data")
    geodesic_data = load_simulation_data(geodesic_dir, snapshot_num=0)
    if geodesic_data is not None:
        print(f"  Loaded {len(geodesic_data['r'])} particles for geodesic flow")
    else:
        print("  No geodesic simulation data found. Run:")
        print("  ./build/sph simulations/relativistic/gr_gsph/config/presets/gr_bondi_geodesic.json")

    sonic_data = load_simulation_data(sonic_dir, snapshot_num=0)
    if sonic_data is not None:
        print(f"  Loaded {len(sonic_data['r'])} particles for sonic point flow")
    else:
        print("  No sonic simulation data found. Run:")
        print("  ./build/sph simulations/relativistic/gr_gsph/config/presets/gr_bondi_sonic.json")

    print("\n[1] Bondi Geodesic Flow (Fig 14)")
    plot_bondi_geodesic(output_dir, sim_data=geodesic_data)

    print("\n[2] Bondi Sonic Point Flow (Fig 18)")
    plot_bondi_sonic(output_dir, sim_data=sonic_data)

    print("\n[3] Combined Bondi Figure")
    plot_combined_bondi(output_dir, geodesic_data=geodesic_data, sonic_data=sonic_data)

    print("\n" + "=" * 70)
    print("  BONDI BENCHMARKS COMPLETE")
    print("=" * 70)

    if geodesic_data is not None:
        print(f"\nBlue points: C++ GR-GSPH simulation ({len(geodesic_data['r'])} particles)")
        print("Green lines: Analytic solution")

    print("""
NOTES:
======
Bondi accretion in Schwarzschild spacetime.

For geodesic flow (pressure-less):
  v^r = -sqrt(2M/r) * (1 - 2M/r)
  rho = rho_0 / r^2 * sqrt(r / 2M)
  u = u_0 * (r^2 * sqrt(2M/r))^(1-gamma)

For sonic point flow:
  Critical point at r_c = 8M
  Flow becomes supersonic inside r_c
""")


if __name__ == '__main__':
    main()
