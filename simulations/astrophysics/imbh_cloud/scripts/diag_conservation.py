#!/usr/bin/env python3
"""
Conservation Diagnostics for IMBH-Cloud Simulations
====================================================

Verifies conservation laws in IMBH tidal disruption simulations, accounting
for particles accreted by the central black hole (sink accretion).

Quantities Tracked:
-------------------
1. Total Energy (E_tot = E_kin + E_therm + E_pot + E_BH)
   - E_kin: Kinetic energy = sum(0.5 * m * v^2)
   - E_therm: Thermal energy = sum(m * u)
   - E_pot: Self-gravitational potential = 0.5 * sum(m * phi_self)
   - E_BH: External BH potential = sum(m * phi_BH)  [NOT halved, external field]

2. Linear Momentum (P = sum(m * v) + P_accreted)
   - Tracks x, y, z components
   - Accounts for momentum removed by sink accretion

3. Angular Momentum (L = sum(m * r x v) + L_accreted)
   - About the BH position (typically origin)
   - Accounts for angular momentum removed by sink accretion

4. Mass Conservation (M_tot = M_gas + M_accreted)
   - Tracks mass removed by sink accretion

Output:
-------
- conservation_diagnostics.csv: Time series of all quantities
- conservation_plots.png: Multi-panel diagnostic plots

Usage:
------
    python3 diag_conservation.py <results_dir> [--output-dir OUTPUT] [--bh-mass M_BH]

    # Example:
    python3 diag_conservation.py results/CAT_OKA/A_61k/oka_corrected --bh-mass 100.0

Notes:
------
- BH mass is in code units (typically 1 code mass = 1000 M_sun)
- For fixed BH at origin, BH position = (0, 0, 0)
- Sink accretion is inferred from particle count changes between snapshots

Units (IMBH_ENCOUNTER):
-----------------------
- [L] = pc
- [M] = 1000 M_sun
- [V] = km/s
- [T] = 0.978 Myr
- [E] = 1000 M_sun * (km/s)^2

References:
-----------
- Oka et al. (2017) Nature Astronomy - CO-0.40-0.22 IMBH candidate
"""

import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, Tuple, Dict, List
import warnings

warnings.filterwarnings('ignore', category=RuntimeWarning)


# =============================================================================
# PHYSICAL CONSTANTS (CODE UNITS)
# =============================================================================

G_CODE = 1.0  # Gravitational constant in code units
GAMMA = 5.0 / 3.0  # Adiabatic index


# =============================================================================
# DATA I/O
# =============================================================================

def load_snapshot(filepath: Path) -> Optional[pd.DataFrame]:
    """
    Load a single CSV snapshot file.

    Returns DataFrame with columns: id, pos_x, pos_y, pos_z, vel_x, vel_y, vel_z,
    mass, dens, pres, ene, sml, phi, etc.
    """
    try:
        # Find header line (starts with "id,")
        with open(filepath, 'r') as f:
            skiprows = 0
            for line in f:
                if line.strip().startswith('id,'):
                    break
                skiprows += 1

        df = pd.read_csv(filepath, skiprows=skiprows)
        return df
    except Exception as e:
        print(f"  Warning: Failed to load {filepath}: {e}")
        return None


def get_snapshot_files(results_dir: Path) -> List[Path]:
    """Get sorted list of snapshot files."""
    snapshots = sorted(results_dir.glob('snapshot_*.csv'))
    return snapshots


def extract_time_from_snapshot(filepath: Path) -> Optional[float]:
    """Extract simulation time from snapshot header."""
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('# time,'):
                    parts = line.strip().split(',')
                    if len(parts) >= 2:
                        return float(parts[1])
                if line.startswith('id,'):
                    break
        return None
    except Exception:
        return None


def load_energy_file(filepath: Path) -> Optional[pd.DataFrame]:
    """Load energy.dat file if it exists."""
    if not filepath.exists():
        return None
    try:
        data = []
        with open(filepath, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    values = line.split()
                    if len(values) >= 5:
                        data.append([float(x) for x in values[:6]])

        if not data:
            return None

        df = pd.DataFrame(data, columns=['time', 'E_kin', 'E_therm', 'E_pot', 'E_tot', 'E_bh'][:len(data[0])])
        return df
    except Exception:
        return None


# =============================================================================
# CONSERVATION CALCULATIONS
# =============================================================================

def compute_kinetic_energy(df: pd.DataFrame) -> float:
    """Compute total kinetic energy: E_kin = sum(0.5 * m * v^2)"""
    v_squared = df['vel_x']**2 + df['vel_y']**2 + df['vel_z']**2
    return 0.5 * np.sum(df['mass'] * v_squared)


def compute_thermal_energy(df: pd.DataFrame) -> float:
    """Compute total thermal energy: E_therm = sum(m * u)"""
    return np.sum(df['mass'] * df['ene'])


def compute_self_gravity_potential(df: pd.DataFrame) -> float:
    """
    Compute self-gravitational potential energy: E_pot = 0.5 * sum(m * phi)

    The factor 0.5 accounts for double-counting in pairwise interactions.
    The phi field should contain the self-gravitational potential.
    """
    if 'phi' not in df.columns:
        return 0.0
    return 0.5 * np.sum(df['mass'] * df['phi'])


def compute_bh_potential_energy(
    df: pd.DataFrame,
    bh_mass: float,
    bh_pos: np.ndarray,
    softening: float = 0.001
) -> float:
    """
    Compute external BH potential energy: E_BH = sum(m * phi_BH)

    Using Plummer softening:
        phi_BH(r) = -G * M_BH / sqrt(r^2 + eps^2)

    No factor of 0.5 since this is an external field, not self-interaction.
    """
    dx = df['pos_x'].values - bh_pos[0]
    dy = df['pos_y'].values - bh_pos[1]
    dz = df['pos_z'].values - bh_pos[2]
    r_squared = dx**2 + dy**2 + dz**2

    phi_bh = -G_CODE * bh_mass / np.sqrt(r_squared + softening**2)
    return np.sum(df['mass'] * phi_bh)


def compute_linear_momentum(df: pd.DataFrame) -> np.ndarray:
    """Compute total linear momentum: P = sum(m * v)"""
    px = np.sum(df['mass'] * df['vel_x'])
    py = np.sum(df['mass'] * df['vel_y'])
    pz = np.sum(df['mass'] * df['vel_z'])
    return np.array([px, py, pz])


def compute_angular_momentum(
    df: pd.DataFrame,
    origin: np.ndarray = np.array([0.0, 0.0, 0.0])
) -> np.ndarray:
    """
    Compute total angular momentum about origin: L = sum(m * r x v)

    For IMBH simulations, origin is typically the BH position.
    """
    rx = df['pos_x'].values - origin[0]
    ry = df['pos_y'].values - origin[1]
    rz = df['pos_z'].values - origin[2]

    vx = df['vel_x'].values
    vy = df['vel_y'].values
    vz = df['vel_z'].values

    m = df['mass'].values

    # L = r x (m*v)
    Lx = np.sum(m * (ry * vz - rz * vy))
    Ly = np.sum(m * (rz * vx - rx * vz))
    Lz = np.sum(m * (rx * vy - ry * vx))

    return np.array([Lx, Ly, Lz])


def compute_total_mass(df: pd.DataFrame) -> float:
    """Compute total gas mass."""
    return np.sum(df['mass'])


def compute_center_of_mass(df: pd.DataFrame) -> np.ndarray:
    """Compute center of mass position."""
    M = compute_total_mass(df)
    if M == 0:
        return np.array([0.0, 0.0, 0.0])

    cx = np.sum(df['mass'] * df['pos_x']) / M
    cy = np.sum(df['mass'] * df['pos_y']) / M
    cz = np.sum(df['mass'] * df['pos_z']) / M
    return np.array([cx, cy, cz])


# =============================================================================
# MAIN DIAGNOSTICS ROUTINE
# =============================================================================

def run_conservation_diagnostics(
    results_dir: Path,
    output_dir: Path,
    bh_mass: float = 100.0,
    bh_pos: np.ndarray = np.array([0.0, 0.0, 0.0]),
    bh_softening: float = 0.001
) -> pd.DataFrame:
    """
    Run full conservation diagnostics on simulation results.

    Returns DataFrame with time series of all conservation quantities.
    """
    snapshots = get_snapshot_files(results_dir)

    if not snapshots:
        raise ValueError(f"No snapshots found in {results_dir}")

    print(f"Found {len(snapshots)} snapshots")
    print(f"BH parameters: M={bh_mass}, pos={bh_pos}, eps={bh_softening}")

    # Storage for results
    records = []

    # Track accreted quantities
    accreted_mass = 0.0
    accreted_momentum = np.array([0.0, 0.0, 0.0])
    accreted_angular_momentum = np.array([0.0, 0.0, 0.0])
    accreted_energy = 0.0

    prev_df = None
    prev_ids = set()

    for i, snap_path in enumerate(snapshots):
        df = load_snapshot(snap_path)
        if df is None:
            continue

        # Get time
        time = extract_time_from_snapshot(snap_path)
        if time is None:
            time = float(i)

        current_ids = set(df['id'].values)

        # Identify accreted particles (present before, missing now)
        if prev_df is not None and prev_ids:
            accreted_ids = prev_ids - current_ids
            if accreted_ids:
                # Get properties of accreted particles from previous snapshot
                accreted_mask = prev_df['id'].isin(accreted_ids)
                accreted_df = prev_df[accreted_mask]

                # Accumulate accreted quantities
                accreted_mass += compute_total_mass(accreted_df)
                accreted_momentum += compute_linear_momentum(accreted_df)
                accreted_angular_momentum += compute_angular_momentum(accreted_df, bh_pos)

                # Accreted energy (kinetic + thermal + BH potential)
                E_kin_acc = compute_kinetic_energy(accreted_df)
                E_therm_acc = compute_thermal_energy(accreted_df)
                E_bh_acc = compute_bh_potential_energy(accreted_df, bh_mass, bh_pos, bh_softening)
                accreted_energy += E_kin_acc + E_therm_acc + E_bh_acc

                if len(accreted_ids) > 0:
                    print(f"  t={time:.3f}: {len(accreted_ids)} particles accreted (total mass: {accreted_mass:.4f})")

        # Compute current quantities
        E_kin = compute_kinetic_energy(df)
        E_therm = compute_thermal_energy(df)
        E_pot = compute_self_gravity_potential(df)
        E_bh = compute_bh_potential_energy(df, bh_mass, bh_pos, bh_softening)
        E_tot = E_kin + E_therm + E_pot + E_bh

        # Energy including accreted contribution
        E_tot_corrected = E_tot + accreted_energy

        P = compute_linear_momentum(df)
        P_corrected = P + accreted_momentum

        L = compute_angular_momentum(df, bh_pos)
        L_corrected = L + accreted_angular_momentum

        M_gas = compute_total_mass(df)
        M_tot = M_gas + accreted_mass

        CoM = compute_center_of_mass(df)

        record = {
            'time': time,
            'n_particles': len(df),
            'n_accreted_total': len(prev_ids) - len(current_ids) if prev_ids else 0,

            # Energies
            'E_kin': E_kin,
            'E_therm': E_therm,
            'E_pot': E_pot,
            'E_bh': E_bh,
            'E_tot': E_tot,
            'E_accreted': accreted_energy,
            'E_tot_corrected': E_tot_corrected,

            # Momenta (uncorrected)
            'Px': P[0],
            'Py': P[1],
            'Pz': P[2],
            'P_mag': np.linalg.norm(P),

            # Momenta (corrected for accretion)
            'Px_corr': P_corrected[0],
            'Py_corr': P_corrected[1],
            'Pz_corr': P_corrected[2],
            'P_mag_corr': np.linalg.norm(P_corrected),

            # Angular momenta (uncorrected)
            'Lx': L[0],
            'Ly': L[1],
            'Lz': L[2],
            'L_mag': np.linalg.norm(L),

            # Angular momenta (corrected for accretion)
            'Lx_corr': L_corrected[0],
            'Ly_corr': L_corrected[1],
            'Lz_corr': L_corrected[2],
            'L_mag_corr': np.linalg.norm(L_corrected),

            # Mass
            'M_gas': M_gas,
            'M_accreted': accreted_mass,
            'M_tot': M_tot,

            # Center of mass
            'CoM_x': CoM[0],
            'CoM_y': CoM[1],
            'CoM_z': CoM[2],
        }
        records.append(record)

        prev_df = df.copy()
        prev_ids = current_ids

        if (i + 1) % 10 == 0:
            print(f"  Processed {i + 1}/{len(snapshots)} snapshots...")

    result_df = pd.DataFrame(records)

    # Compute conservation errors
    if len(result_df) > 0:
        E0 = result_df['E_tot_corrected'].iloc[0]
        result_df['E_drift_pct'] = (result_df['E_tot_corrected'] - E0) / np.abs(E0) * 100

        L0 = result_df['L_mag_corr'].iloc[0]
        if L0 > 0:
            result_df['L_drift_pct'] = (result_df['L_mag_corr'] - L0) / L0 * 100
        else:
            result_df['L_drift_pct'] = 0.0

        M0 = result_df['M_tot'].iloc[0]
        result_df['M_drift_pct'] = (result_df['M_tot'] - M0) / M0 * 100

    return result_df


# =============================================================================
# VISUALIZATION
# =============================================================================

def create_conservation_plots(
    df: pd.DataFrame,
    output_path: Path,
    title: str = "Conservation Diagnostics"
) -> None:
    """Create multi-panel conservation diagnostic plots."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("  Warning: matplotlib not available, skipping plots")
        return

    fig, axes = plt.subplots(3, 2, figsize=(14, 12))

    t = df['time']

    # --- Panel 1: Energy Components ---
    ax = axes[0, 0]
    ax.plot(t, df['E_kin'], 'g-', label='Kinetic', linewidth=1.5)
    ax.plot(t, df['E_therm'], 'r-', label='Thermal', linewidth=1.5)
    ax.plot(t, df['E_pot'], 'b-', label='Self-grav', linewidth=1.5)
    ax.plot(t, df['E_bh'], 'm-', label='BH potential', linewidth=1.5)
    ax.plot(t, df['E_tot'], 'k-', label='Total', linewidth=2)
    ax.set_xlabel('Time')
    ax.set_ylabel('Energy')
    ax.set_title('Energy Components')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # --- Panel 2: Energy Conservation Error ---
    ax = axes[0, 1]
    ax.plot(t, df['E_drift_pct'], 'k-', linewidth=1.5)
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.axhline(y=1, color='r', linestyle=':', alpha=0.5, label='1% threshold')
    ax.axhline(y=-1, color='r', linestyle=':', alpha=0.5)
    ax.set_xlabel('Time')
    ax.set_ylabel('Energy Drift (%)')
    ax.set_title('Energy Conservation (corrected for accretion)')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # --- Panel 3: Linear Momentum ---
    ax = axes[1, 0]
    ax.plot(t, df['Px_corr'], 'r-', label='Px', linewidth=1)
    ax.plot(t, df['Py_corr'], 'g-', label='Py', linewidth=1)
    ax.plot(t, df['Pz_corr'], 'b-', label='Pz', linewidth=1)
    ax.plot(t, df['P_mag_corr'], 'k-', label='|P|', linewidth=2)
    ax.set_xlabel('Time')
    ax.set_ylabel('Momentum')
    ax.set_title('Linear Momentum (corrected for accretion)')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # --- Panel 4: Angular Momentum ---
    ax = axes[1, 1]
    ax.plot(t, df['Lx_corr'], 'r-', label='Lx', linewidth=1)
    ax.plot(t, df['Ly_corr'], 'g-', label='Ly', linewidth=1)
    ax.plot(t, df['Lz_corr'], 'b-', label='Lz', linewidth=1)
    ax.plot(t, df['L_mag_corr'], 'k-', label='|L|', linewidth=2)
    ax.set_xlabel('Time')
    ax.set_ylabel('Angular Momentum')
    ax.set_title('Angular Momentum (corrected for accretion)')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # --- Panel 5: Mass Conservation ---
    ax = axes[2, 0]
    ax.plot(t, df['M_gas'], 'b-', label='Gas mass', linewidth=1.5)
    ax.plot(t, df['M_accreted'], 'r-', label='Accreted mass', linewidth=1.5)
    ax.plot(t, df['M_tot'], 'k--', label='Total mass', linewidth=2)
    ax.set_xlabel('Time')
    ax.set_ylabel('Mass')
    ax.set_title('Mass Conservation')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # --- Panel 6: Particle Count ---
    ax = axes[2, 1]
    ax.plot(t, df['n_particles'], 'b-', linewidth=1.5)
    ax.set_xlabel('Time')
    ax.set_ylabel('Number of Particles')
    ax.set_title('Particle Count (decreases due to sink accretion)')
    ax.grid(True, alpha=0.3)

    # Add particle loss info
    if len(df) > 0:
        n_initial = df['n_particles'].iloc[0]
        n_final = df['n_particles'].iloc[-1]
        n_lost = n_initial - n_final
        ax.text(0.95, 0.95, f'Initial: {n_initial}\nFinal: {n_final}\nAccreted: {n_lost}',
                transform=ax.transAxes, fontsize=9, verticalalignment='top',
                horizontalalignment='right', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.suptitle(title, fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(output_path, dpi=150)
    plt.close()

    print(f"  Saved: {output_path}")


def print_summary(df: pd.DataFrame) -> None:
    """Print conservation summary to console."""
    print("\n" + "=" * 70)
    print("  CONSERVATION DIAGNOSTICS SUMMARY")
    print("=" * 70)

    if len(df) == 0:
        print("  No data to summarize")
        return

    # Energy
    E0 = df['E_tot_corrected'].iloc[0]
    E_final = df['E_tot_corrected'].iloc[-1]
    E_drift = df['E_drift_pct'].iloc[-1]

    print(f"\n  Energy Conservation:")
    print(f"    Initial E:  {E0:.6e}")
    print(f"    Final E:    {E_final:.6e}")
    print(f"    Drift:      {E_drift:+.4f}%")
    status = "PASS" if abs(E_drift) < 1.0 else "MARGINAL" if abs(E_drift) < 5.0 else "FAIL"
    print(f"    Status:     {status}")

    # Angular momentum
    L0 = df['L_mag_corr'].iloc[0]
    L_final = df['L_mag_corr'].iloc[-1]
    L_drift = df['L_drift_pct'].iloc[-1] if 'L_drift_pct' in df else 0.0

    print(f"\n  Angular Momentum Conservation:")
    print(f"    Initial |L|: {L0:.6e}")
    print(f"    Final |L|:   {L_final:.6e}")
    print(f"    Drift:       {L_drift:+.4f}%")

    # Mass
    M0 = df['M_tot'].iloc[0]
    M_final = df['M_tot'].iloc[-1]
    M_accreted = df['M_accreted'].iloc[-1]

    print(f"\n  Mass Conservation:")
    print(f"    Initial M:   {M0:.6f}")
    print(f"    Final M:     {M_final:.6f}")
    print(f"    Accreted:    {M_accreted:.6f} ({M_accreted/M0*100:.2f}%)")
    print(f"    Status:      {'PASS' if abs(M_final - M0) < 1e-10 else 'FAIL'}")

    # Particles
    n_initial = df['n_particles'].iloc[0]
    n_final = df['n_particles'].iloc[-1]

    print(f"\n  Particle Tracking:")
    print(f"    Initial:     {n_initial}")
    print(f"    Final:       {n_final}")
    print(f"    Accreted:    {n_initial - n_final}")

    print("\n" + "=" * 70)


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Conservation diagnostics for IMBH-cloud simulations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('results_dir', type=Path,
                        help='Directory containing snapshot_*.csv files')
    parser.add_argument('--output-dir', '-o', type=Path, default=None,
                        help='Output directory (default: results_dir)')
    parser.add_argument('--bh-mass', type=float, default=100.0,
                        help='BH mass in code units (default: 100.0 = 10^5 Msun)')
    parser.add_argument('--bh-pos', type=float, nargs=3, default=[0.0, 0.0, 0.0],
                        help='BH position (default: 0 0 0)')
    parser.add_argument('--bh-softening', type=float, default=0.001,
                        help='BH softening length (default: 0.001 pc)')
    parser.add_argument('--no-plots', action='store_true',
                        help='Skip generating plots')

    args = parser.parse_args()

    if not args.results_dir.exists():
        print(f"Error: Results directory not found: {args.results_dir}")
        return 1

    output_dir = args.output_dir or args.results_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nConservation Diagnostics")
    print(f"========================")
    print(f"Results:    {args.results_dir}")
    print(f"Output:     {output_dir}")

    # Run diagnostics
    bh_pos = np.array(args.bh_pos)
    df = run_conservation_diagnostics(
        args.results_dir,
        output_dir,
        bh_mass=args.bh_mass,
        bh_pos=bh_pos,
        bh_softening=args.bh_softening
    )

    # Save CSV
    csv_path = output_dir / 'conservation_diagnostics.csv'
    df.to_csv(csv_path, index=False)
    print(f"\nSaved: {csv_path}")

    # Generate plots
    if not args.no_plots:
        plot_path = output_dir / 'conservation_plots.png'
        create_conservation_plots(df, plot_path, title=f"Conservation: {args.results_dir.name}")

    # Print summary
    print_summary(df)

    return 0


if __name__ == '__main__':
    exit(main())
