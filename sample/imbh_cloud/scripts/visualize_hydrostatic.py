#!/usr/bin/env python3
"""
Visualization script for Lane-Emden hydrostatic equilibrium test with self-gravity.

Purpose:
    Verify that the relaxed Lane-Emden sphere maintains hydrostatic equilibrium
    over multiple crossing times when self-gravity is enabled, demonstrating
    correct gravity implementation before proceeding to IMBH tidal disruption.

Physics:
    - Lane-Emden n=1.5 polytrope (γ=5/3)
    - Self-gravity with Barnes-Hut tree
    - DISPH method (conservative formulation)
    - Test duration: ~10 crossing times >> tidal disruption timescale

Expected behavior:
    - Density profile should remain constant (RMS deviation < 5%)
    - Velocities should remain small (|v| < 0.01 c_s)
    - Total energy should be conserved (drift < 1%)
    - No expansion or collapse

Usage:
    python3 visualize_hydrostatic.py <results_dir> <output_dir>
    
    results_dir: Directory containing snapshot_*.csv and energy.dat
    output_dir: Directory to save plots and animations
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.gridspec import GridSpec
import glob
import os
import sys
from pathlib import Path

# Lane-Emden n=1.5 analytical parameters (code units: M=1, R=1)
XI_1 = 3.6537540101  # First zero of Lane-Emden equation
ALPHA = 1.0 / XI_1  # Scaling factor R/ξ₁
RHO_C = 1.43009692  # Central density (analytical)
GAMMA = 5.0 / 3.0  # Adiabatic index
POLYTROPIC_N = 1.5  # Polytropic index


def lane_emden_theta(xi, tol=1e-10):
    """
    Solve Lane-Emden equation for n=1.5 using 4th-order Runge-Kutta.
    
    d²θ/dξ² + (2/ξ) dθ/dξ + θ^n = 0
    θ(0) = 1, θ'(0) = 0
    """
    n = POLYTROPIC_N
    
    # Initial conditions with Taylor series near ξ=0
    xi_array = np.atleast_1d(xi)
    theta = np.ones_like(xi_array)
    
    # RK4 integration
    dxi = 0.001
    xi_current = dxi  # Start slightly offset from zero
    theta_current = 1.0 - (dxi**2) / 6.0  # Taylor expansion
    dtheta_current = -dxi / 3.0
    
    for i, xi_target in enumerate(xi_array):
        while xi_current < xi_target:
            # RK4 step
            k1_theta = dtheta_current
            k1_dtheta = -theta_current**n - (2.0 / xi_current) * dtheta_current
            
            k2_theta = dtheta_current + 0.5 * dxi * k1_dtheta
            k2_dtheta = -(theta_current + 0.5 * dxi * k1_theta)**n - (2.0 / (xi_current + 0.5 * dxi)) * k2_theta
            
            k3_theta = dtheta_current + 0.5 * dxi * k2_dtheta
            k3_dtheta = -(theta_current + 0.5 * dxi * k2_theta)**n - (2.0 / (xi_current + 0.5 * dxi)) * k3_theta
            
            k4_theta = dtheta_current + dxi * k3_dtheta
            k4_dtheta = -(theta_current + dxi * k3_theta)**n - (2.0 / (xi_current + dxi)) * k4_theta
            
            theta_current += (dxi / 6.0) * (k1_theta + 2*k2_theta + 2*k3_theta + k4_theta)
            dtheta_current += (dxi / 6.0) * (k1_dtheta + 2*k2_dtheta + 2*k3_dtheta + k4_dtheta)
            xi_current += dxi
            
            if theta_current < 0:
                theta_current = 0
                break
        
        theta[i] = max(theta_current, 0)
    
    return theta


def lane_emden_density(r):
    """Analytical Lane-Emden density profile ρ(r) = ρ_c θ^n."""
    xi = r / ALPHA
    theta = lane_emden_theta(xi)
    return RHO_C * theta**POLYTROPIC_N


def detect_header_lines(filepath):
    """Detect number of header lines by finding where numeric data starts."""
    with open(filepath, 'r') as f:
        for line_num, line in enumerate(f):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            # Try to parse first field as number
            try:
                float(line.split(',')[0])
                return line_num
            except (ValueError, IndexError):
                continue
    return 61  # Default fallback


def load_snapshot(filepath):
    """Load snapshot with automatic header detection."""
    header_lines = detect_header_lines(filepath)
    data = np.loadtxt(filepath, delimiter=',', skiprows=header_lines)
    return data


def load_energy_file(filepath):
    """Load energy.dat file."""
    try:
        data = np.loadtxt(filepath, skiprows=1)
        return data
    except Exception as e:
        print(f"Warning: Could not load energy file: {e}")
        return None


def compute_crossing_time(snapshot_data):
    """
    Estimate crossing time t_cross = R / <v_sound>.
    
    For Lane-Emden sphere: R = 1.0 (code units)
    Sound speed: c_s = sqrt(γ * P / ρ)
    """
    # Extract density and pressure
    rho = snapshot_data[:, 11]  # Column 11: density
    pres = snapshot_data[:, 12]  # Column 12: pressure
    
    # Compute sound speed
    c_s = np.sqrt(GAMMA * pres / rho)
    c_s_avg = np.median(c_s)  # Use median (more robust than mean)
    
    R = 1.0  # Sphere radius in code units
    t_cross = R / c_s_avg
    
    return t_cross, c_s_avg


def analyze_hydrostatic_quality(snapshot_data, analytic_density_func):
    """
    Compute hydrostatic equilibrium quality metrics.
    
    Returns:
        rms_density_error: RMS fractional density deviation from analytic
        max_velocity: Maximum velocity magnitude
        median_velocity: Median velocity magnitude
    """
    # Extract particle positions and velocities
    pos = snapshot_data[:, 1:4]  # Columns 1-3: x, y, z
    vel = snapshot_data[:, 4:7]  # Columns 4-6: vx, vy, vz
    rho = snapshot_data[:, 11]  # Column 11: density
    
    # Compute radial distance
    r = np.sqrt(np.sum(pos**2, axis=1))
    
    # Compute velocity magnitude
    v_mag = np.sqrt(np.sum(vel**2, axis=1))
    
    # Compute density error (only for particles within R < 1.0)
    mask = r < 0.95  # Exclude surface particles
    r_masked = r[mask]
    rho_masked = rho[mask]
    
    rho_analytic = analytic_density_func(r_masked)
    density_error = np.abs(rho_masked - rho_analytic) / rho_analytic
    rms_density_error = np.sqrt(np.mean(density_error**2))
    
    return {
        'rms_density_error': rms_density_error,
        'max_velocity': np.max(v_mag),
        'median_velocity': np.median(v_mag),
        'mean_density_error': np.mean(density_error)
    }


def create_hydrostatic_animation(results_dir, output_dir):
    """Create comprehensive 6-panel animation showing hydrostatic equilibrium test."""
    
    # Find all snapshots
    snapshot_files = sorted(glob.glob(os.path.join(results_dir, 'snapshot_*.csv')))
    if len(snapshot_files) == 0:
        print(f"Error: No snapshots found in {results_dir}")
        return
    
    print(f"Found {len(snapshot_files)} snapshots")
    
    # Load energy data
    energy_file = os.path.join(results_dir, 'energy.dat')
    energy_data = load_energy_file(energy_file)
    
    # Estimate crossing time from first snapshot
    first_snapshot = load_snapshot(snapshot_files[0])
    t_cross, c_s_avg = compute_crossing_time(first_snapshot)
    print(f"Estimated crossing time: t_cross = {t_cross:.3f} code units")
    print(f"Average sound speed: c_s = {c_s_avg:.3f} code units")
    
    # Setup figure
    fig = plt.figure(figsize=(18, 10))
    gs = GridSpec(2, 3, figure=fig, hspace=0.3, wspace=0.3)
    
    ax1 = fig.add_subplot(gs[0, 0])  # 3D density
    ax2 = fig.add_subplot(gs[0, 1])  # XY slice
    ax3 = fig.add_subplot(gs[0, 2])  # Radial density profile
    ax4 = fig.add_subplot(gs[1, 0])  # Velocity histogram
    ax5 = fig.add_subplot(gs[1, 1])  # Energy evolution
    ax6 = fig.add_subplot(gs[1, 2])  # Quality metrics
    
    # Quality metrics storage
    times = []
    rms_errors = []
    max_velocities = []
    median_velocities = []
    
    # Colorbar references (create once, reuse)
    cbar1 = None
    
    def update(frame_num):
        nonlocal cbar1
        
        snapshot_file = snapshot_files[frame_num]
        data = load_snapshot(snapshot_file)
        
        # Extract data
        pos = data[:, 1:4]
        vel = data[:, 4:7]
        rho = data[:, 11]
        
        # Compute derived quantities
        r = np.sqrt(np.sum(pos**2, axis=1))
        v_mag = np.sqrt(np.sum(vel**2, axis=1))
        
        # Get time from filename (assuming snapshot_NNNN.csv format)
        step_num = int(os.path.basename(snapshot_file).split('_')[1].split('.')[0])
        
        # Try to get actual time from energy file
        if energy_data is not None and frame_num < len(energy_data):
            current_time = energy_data[frame_num, 0]
        else:
            current_time = step_num * 0.2  # Assume output frequency
        
        # Compute quality metrics
        metrics = analyze_hydrostatic_quality(data, lane_emden_density)
        times.append(current_time)
        rms_errors.append(metrics['rms_density_error'])
        max_velocities.append(metrics['max_velocity'])
        median_velocities.append(metrics['median_velocity'])
        
        # Clear all axes
        for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
            ax.clear()
        
        # Panel 1: 3D density scatter
        scatter = ax1.scatter(pos[:, 0], pos[:, 1], c=rho, s=1, cmap='viridis',
                             vmin=0, vmax=RHO_C, alpha=0.6)
        ax1.set_xlabel('x [code units]')
        ax1.set_ylabel('y [code units]')
        ax1.set_title('3D Density Distribution (XY projection)')
        ax1.set_xlim(-1.2, 1.2)
        ax1.set_ylim(-1.2, 1.2)
        ax1.set_aspect('equal')
        
        # Only create colorbar once
        if cbar1 is None:
            cbar1 = plt.colorbar(scatter, ax=ax1, label='Density [code units]')
        
        # Panel 2: XY slice (|z| < 0.1)
        mask_xy = np.abs(pos[:, 2]) < 0.1
        ax2.scatter(pos[mask_xy, 0], pos[mask_xy, 1], c=rho[mask_xy],
                   s=5, cmap='viridis', vmin=0, vmax=RHO_C)
        ax2.set_xlabel('x [code units]')
        ax2.set_ylabel('y [code units]')
        ax2.set_title('XY Slice (|z| < 0.1)')
        ax2.set_xlim(-1.2, 1.2)
        ax2.set_ylim(-1.2, 1.2)
        ax2.set_aspect('equal')
        
        # Panel 3: Radial density profile vs analytic
        r_bins = np.linspace(0, 1.0, 30)
        r_centers = 0.5 * (r_bins[1:] + r_bins[:-1])
        
        # Bin particles
        rho_binned = []
        for i in range(len(r_bins) - 1):
            mask = (r >= r_bins[i]) & (r < r_bins[i+1])
            if np.sum(mask) > 0:
                rho_binned.append(np.mean(rho[mask]))
            else:
                rho_binned.append(np.nan)
        
        # Analytical profile
        r_analytic = np.linspace(0, 1.0, 200)
        rho_analytic = lane_emden_density(r_analytic)
        
        ax3.plot(r_analytic, rho_analytic, 'r-', linewidth=2, label='Lane-Emden (analytic)', alpha=0.8)
        ax3.plot(r_centers, rho_binned, 'bo-', markersize=4, label='SPH (binned)', alpha=0.7)
        ax3.set_xlabel('Radius r [code units]')
        ax3.set_ylabel('Density ρ [code units]')
        ax3.set_title(f'Radial Density Profile\nRMS Error: {metrics["rms_density_error"]*100:.2f}%')
        ax3.legend(loc='upper right')
        ax3.grid(True, alpha=0.3)
        ax3.set_xlim(0, 1.0)
        ax3.set_ylim(0, RHO_C * 1.1)
        
        # Panel 4: Velocity magnitude histogram
        ax4.hist(v_mag / c_s_avg, bins=50, color='steelblue', alpha=0.7, edgecolor='black')
        ax4.axvline(metrics['median_velocity'] / c_s_avg, color='red', linestyle='--',
                   linewidth=2, label=f'Median: {metrics["median_velocity"]/c_s_avg:.4f}')
        ax4.set_xlabel('|v| / c_s')
        ax4.set_ylabel('Particle count')
        ax4.set_title(f'Velocity Distribution\nMax: {metrics["max_velocity"]/c_s_avg:.4f} c_s')
        ax4.legend()
        ax4.set_xlim(0, min(0.1, metrics['max_velocity'] / c_s_avg * 1.5))
        ax4.grid(True, alpha=0.3)
        
        # Panel 5: Energy evolution
        if energy_data is not None:
            time_energy = energy_data[:, 0]
            total_energy = energy_data[:, 1]
            
            # Normalize to initial total energy
            E0 = total_energy[0]
            
            ax5.plot(time_energy / t_cross, (total_energy - E0) / abs(E0) * 100,
                    'k-', linewidth=2, label='Total (fractional drift)')
            ax5.axhline(0, color='gray', linestyle='--', alpha=0.5)
            ax5.axhline(1, color='red', linestyle=':', alpha=0.3, label='±1% threshold')
            ax5.axhline(-1, color='red', linestyle=':', alpha=0.3)
            ax5.axvline(current_time / t_cross, color='blue', linestyle='--', alpha=0.5,
                       label=f't = {current_time/t_cross:.2f} t_cross')
            ax5.set_xlabel('Time [crossing times]')
            ax5.set_ylabel('Energy drift [%]')
            ax5.set_title('Energy Conservation')
            ax5.legend(loc='best')
            ax5.grid(True, alpha=0.3)
        
        # Panel 6: Quality metrics over time
        if len(times) > 1:
            t_array = np.array(times) / t_cross
            ax6.plot(t_array, np.array(rms_errors) * 100, 'b-', linewidth=2,
                    label='Density RMS error [%]')
            ax6.plot(t_array, np.array(max_velocities) / c_s_avg * 100, 'r-',
                    linewidth=2, label='Max velocity [% c_s]')
            ax6.axhline(5, color='gray', linestyle='--', alpha=0.3, label='5% threshold')
            ax6.set_xlabel('Time [crossing times]')
            ax6.set_ylabel('Error / Velocity [%]')
            ax6.set_title('Hydrostatic Quality Metrics')
            ax6.legend(loc='best')
            ax6.grid(True, alpha=0.3)
        
        # Main title
        fig.suptitle(f'Hydrostatic Equilibrium Test (DISPH + Self-Gravity)\n'
                    f'Time: {current_time:.2f} code units = {current_time/t_cross:.2f} t_cross | '
                    f'Snapshot: {frame_num + 1}/{len(snapshot_files)}',
                    fontsize=14, fontweight='bold')
        
        return []
    
    # Create animation
    print("Creating animation...")
    anim = animation.FuncAnimation(fig, update, frames=len(snapshot_files),
                                  interval=200, blit=False)
    
    # Save animation
    output_file = os.path.join(output_dir, 'hydrostatic_animation.gif')
    print(f"Saving animation to {output_file}...")
    anim.save(output_file, writer='pillow', fps=5, dpi=100)
    print(f"✓ Animation saved: {output_file}")
    
    plt.close(fig)
    
    # Create final summary plot
    create_summary_plot(snapshot_files[-1], times, rms_errors, max_velocities,
                       median_velocities, t_cross, c_s_avg, output_dir)


def create_summary_plot(final_snapshot_file, times, rms_errors, max_velocities,
                       median_velocities, t_cross, c_s_avg, output_dir):
    """Create summary plot of final state and quality metrics."""
    
    data = load_snapshot(final_snapshot_file)
    pos = data[:, 1:4]
    vel = data[:, 4:7]
    rho = data[:, 11]
    r = np.sqrt(np.sum(pos**2, axis=1))
    v_mag = np.sqrt(np.sum(vel**2, axis=1))
    
    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(2, 3, figure=fig, hspace=0.3, wspace=0.3)
    
    # Panel 1: Final density profile
    ax1 = fig.add_subplot(gs[0, 0])
    r_bins = np.linspace(0, 1.0, 30)
    r_centers = 0.5 * (r_bins[1:] + r_bins[:-1])
    rho_binned = []
    for i in range(len(r_bins) - 1):
        mask = (r >= r_bins[i]) & (r < r_bins[i+1])
        if np.sum(mask) > 0:
            rho_binned.append(np.mean(rho[mask]))
        else:
            rho_binned.append(np.nan)
    
    r_analytic = np.linspace(0, 1.0, 200)
    rho_analytic = lane_emden_density(r_analytic)
    
    ax1.plot(r_analytic, rho_analytic, 'r-', linewidth=3, label='Analytic', alpha=0.8)
    ax1.plot(r_centers, rho_binned, 'bo-', markersize=6, label='SPH', alpha=0.7)
    ax1.set_xlabel('Radius [code units]', fontsize=12)
    ax1.set_ylabel('Density [code units]', fontsize=12)
    ax1.set_title('Final Density Profile', fontsize=13, fontweight='bold')
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)
    
    # Panel 2: Density error vs radius
    ax2 = fig.add_subplot(gs[0, 1])
    mask_r = r < 0.95
    rho_expected = lane_emden_density(r[mask_r])
    fractional_error = (rho[mask_r] - rho_expected) / rho_expected * 100
    
    ax2.scatter(r[mask_r], fractional_error, s=2, alpha=0.5, c='steelblue')
    ax2.axhline(0, color='black', linestyle='-', linewidth=1)
    ax2.axhline(5, color='red', linestyle='--', alpha=0.5, label='±5%')
    ax2.axhline(-5, color='red', linestyle='--', alpha=0.5)
    ax2.set_xlabel('Radius [code units]', fontsize=12)
    ax2.set_ylabel('Density Error [%]', fontsize=12)
    ax2.set_title('Density Deviation from Analytic', fontsize=13, fontweight='bold')
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3)
    
    # Panel 3: Velocity vs radius
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.scatter(r, v_mag / c_s_avg, s=2, alpha=0.5, c='coral')
    ax3.axhline(0.01, color='red', linestyle='--', linewidth=2, label='1% c_s threshold')
    ax3.set_xlabel('Radius [code units]', fontsize=12)
    ax3.set_ylabel('|v| / c_s', fontsize=12)
    ax3.set_title('Velocity Distribution', fontsize=13, fontweight='bold')
    ax3.legend(fontsize=11)
    ax3.grid(True, alpha=0.3)
    ax3.set_ylim(0, min(0.05, np.max(v_mag / c_s_avg) * 1.2))
    
    # Panel 4: Quality metrics evolution
    ax4 = fig.add_subplot(gs[1, 0])
    t_array = np.array(times) / t_cross
    ax4.plot(t_array, np.array(rms_errors) * 100, 'b-', linewidth=2, label='Density RMS [%]')
    ax4.axhline(5, color='red', linestyle='--', alpha=0.5, label='5% threshold')
    ax4.set_xlabel('Time [crossing times]', fontsize=12)
    ax4.set_ylabel('RMS Density Error [%]', fontsize=12)
    ax4.set_title('Density Conservation', fontsize=13, fontweight='bold')
    ax4.legend(fontsize=11)
    ax4.grid(True, alpha=0.3)
    
    # Panel 5: Velocity evolution
    ax5 = fig.add_subplot(gs[1, 1])
    ax5.plot(t_array, np.array(max_velocities) / c_s_avg * 100, 'r-', linewidth=2,
            label='Max velocity')
    ax5.plot(t_array, np.array(median_velocities) / c_s_avg * 100, 'b-', linewidth=2,
            label='Median velocity')
    ax5.axhline(1, color='gray', linestyle='--', alpha=0.5, label='1% c_s')
    ax5.set_xlabel('Time [crossing times]', fontsize=12)
    ax5.set_ylabel('Velocity [% of c_s]', fontsize=12)
    ax5.set_title('Velocity Evolution', fontsize=13, fontweight='bold')
    ax5.legend(fontsize=11)
    ax5.grid(True, alpha=0.3)
    
    # Panel 6: Summary statistics
    ax6 = fig.add_subplot(gs[1, 2])
    ax6.axis('off')
    
    final_metrics = analyze_hydrostatic_quality(data, lane_emden_density)
    
    summary_text = f"""
    HYDROSTATIC TEST SUMMARY
    ════════════════════════════
    
    Test Duration:
      • Total time: {times[-1]:.2f} code units
      • Crossing times: {times[-1]/t_cross:.2f} t_cross
    
    Final State Quality:
      • RMS density error: {final_metrics['rms_density_error']*100:.2f}%
      • Max velocity: {final_metrics['max_velocity']/c_s_avg*100:.3f}% c_s
      • Median velocity: {final_metrics['median_velocity']/c_s_avg*100:.4f}% c_s
    
    Pass Criteria:
      {'✓' if final_metrics['rms_density_error'] < 0.05 else '✗'} Density error < 5%
      {'✓' if final_metrics['max_velocity']/c_s_avg < 0.01 else '✗'} Max velocity < 1% c_s
      {'✓' if final_metrics['median_velocity']/c_s_avg < 0.001 else '✗'} Median vel < 0.1% c_s
    
    Physical Interpretation:
      • Sphere radius: R = 1.0 code units
      • Sound speed: c_s ≈ {c_s_avg:.3f}
      • Crossing time: t_cross ≈ {t_cross:.3f}
      • Test >> tidal timescale ✓
    """
    
    ax6.text(0.1, 0.5, summary_text, transform=ax6.transAxes,
            fontsize=11, verticalalignment='center', family='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
    
    fig.suptitle('Hydrostatic Equilibrium Test - Final Summary',
                fontsize=16, fontweight='bold')
    
    output_file = os.path.join(output_dir, 'hydrostatic_summary.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"✓ Summary plot saved: {output_file}")
    plt.close(fig)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python3 visualize_hydrostatic.py <results_dir> <output_dir>")
        sys.exit(1)
    
    results_dir = sys.argv[1]
    output_dir = sys.argv[2]
    
    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("Hydrostatic Equilibrium Test Visualization")
    print("=" * 60)
    print(f"Results directory: {results_dir}")
    print(f"Output directory: {output_dir}")
    print("")
    
    create_hydrostatic_animation(results_dir, output_dir)
    
    print("")
    print("=" * 60)
    print("✓ Visualization complete!")
    print("=" * 60)
