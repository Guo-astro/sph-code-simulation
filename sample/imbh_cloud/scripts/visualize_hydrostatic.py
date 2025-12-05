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
    # Standard incremental mode (default - only processes new snapshots)
    python3 visualize_hydrostatic.py <results_dir> <output_dir>
    
    # Force full rebuild from scratch
    python3 visualize_hydrostatic.py <results_dir> <output_dir> --rebuild
    
    results_dir: Directory containing snapshot_*.csv and energy.dat
    output_dir: Directory to save plots and animations

Incremental Mode:
    By default, the script operates in incremental mode which:
    - Stores metadata in .hydrostatic_gif_metadata.json
    - Only renders new snapshot frames not previously processed
    - Appends new frames to the existing GIF
    - Much faster for long-running simulations with many snapshots
    
    Use --rebuild to force a complete regeneration of the GIF.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import glob
import io
import os
import sys
import json
import argparse
from pathlib import Path
from PIL import Image

# Lane-Emden n=1.5 analytical parameters (code units: M=1, R=1)
XI_1 = 3.6537540101  # First zero of Lane-Emden equation
ALPHA = 1.0 / XI_1  # Scaling factor R/ξ₁
RHO_C = 1.43009692  # Central density (analytical)
GAMMA = 5.0 / 3.0  # Adiabatic index
POLYTROPIC_N = 1.5  # Polytropic index
G = 1.0  # Gravitational constant

# Theoretical polytropic constant K for isentropic Lane-Emden
# K = 4*pi*G*alpha^2*rho_c^(1/3)/(n+1) for n=1.5
K_EXPECTED = 4.0 * np.pi * G * ALPHA**2 * RHO_C**(1.0/3.0) / (POLYTROPIC_N + 1)


def load_gif_metadata(output_dir):
    """Load metadata about previously processed frames."""
    metadata_file = os.path.join(output_dir, '.hydrostatic_gif_metadata.json')
    if os.path.exists(metadata_file):
        try:
            with open(metadata_file, 'r') as f:
                return json.load(f)
        except (json.JSONDecodeError, IOError):
            pass
    return {'processed_snapshots': [], 'metrics': {'times': [], 'rms_errors': [], 
            'max_velocities': [], 'median_velocities': [], 'K_variations': []}}


def save_gif_metadata(output_dir, metadata):
    """Save metadata about processed frames."""
    metadata_file = os.path.join(output_dir, '.hydrostatic_gif_metadata.json')
    with open(metadata_file, 'w') as f:
        json.dump(metadata, f, indent=2)


def append_frames_to_gif(existing_gif_path, new_frames, output_path, fps=5):
    """
    Append new frames to an existing GIF.
    
    Args:
        existing_gif_path: Path to existing GIF file
        new_frames: List of PIL Image objects to append
        output_path: Path for output GIF (can be same as existing)
        fps: Frames per second
    """
    frames = []
    duration = int(1000 / fps)
    
    # Load existing frames if GIF exists
    if os.path.exists(existing_gif_path):
        try:
            existing_gif = Image.open(existing_gif_path)
            # Extract all frames from existing GIF
            try:
                while True:
                    frames.append(existing_gif.copy())
                    existing_gif.seek(existing_gif.tell() + 1)
            except EOFError:
                pass
            print(f"  Loaded {len(frames)} existing frames from GIF")
        except Exception as e:
            print(f"  Warning: Could not load existing GIF: {e}")
            frames = []
    
    # Append new frames
    frames.extend(new_frames)
    
    if len(frames) > 0:
        # Ensure all frames have the same size (use first frame as reference)
        target_size = frames[0].size
        resized_frames = []
        for frame in frames:
            if frame.size != target_size:
                # Resize to match first frame, maintaining aspect ratio with padding
                frame = frame.resize(target_size, Image.Resampling.LANCZOS)
            resized_frames.append(frame)
        
        # Save combined GIF
        resized_frames[0].save(
            output_path,
            save_all=True,
            append_images=resized_frames[1:],
            duration=duration,
            loop=0
        )
        return len(resized_frames)
    return 0


def lane_emden_theta(xi, tol=1e-10):
    """
    Solve Lane-Emden equation for n=1.5 using 4th-order Runge-Kutta.
    
    d²θ/dξ² + (2/ξ) dθ/dξ + θ^n = 0
    θ(0) = 1, θ'(0) = 0
    
    Note: This function handles unsorted input arrays by sorting internally
    and restoring the original order before returning.
    """
    n = POLYTROPIC_N
    
    # Initial conditions with Taylor series near ξ=0
    xi_array = np.atleast_1d(xi)
    theta = np.ones_like(xi_array, dtype=float)
    
    # Sort the input and remember original order for restoration
    sort_idx = np.argsort(xi_array)
    xi_sorted = xi_array[sort_idx]
    
    # RK4 integration
    dxi = 0.001
    xi_current = dxi  # Start slightly offset from zero
    theta_current = 1.0 - (dxi**2) / 6.0  # Taylor expansion
    dtheta_current = -dxi / 3.0
    
    for sorted_i, xi_target in enumerate(xi_sorted):
        # Handle xi values near or below zero
        if xi_target <= dxi:
            theta[sort_idx[sorted_i]] = 1.0 - xi_target**2 / 6.0
            continue
            
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
        
        # Store result at the correct original position
        theta[sort_idx[sorted_i]] = max(theta_current, 0)
    
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
        data = np.loadtxt(filepath, skiprows=1, ndmin=2)
        # Handle case where file has only one row (returns shape (1, N))
        # or zero rows (returns shape (0, 0) or similar)
        if data.size == 0:
            print("Warning: Energy file is empty")
            return None
        return data
    except Exception as e:
        print(f"Warning: Could not load energy file: {e}")
        return None


def compute_crossing_time(snapshot_data):
    """
    Estimate crossing time t_cross = R / <v_sound>.
    
    For Lane-Emden sphere: R = 1.0 (code units)
    Sound speed is pre-computed in column 15 of snapshot data.
    Falls back to computing from P/ρ if needed.
    """
    # Try to use pre-computed sound speed (column 15)
    if snapshot_data.shape[1] > 15:
        c_s = snapshot_data[:, 15]  # Column 15: sound speed
        # Filter out invalid values
        valid_mask = np.isfinite(c_s) & (c_s > 0)
        if np.sum(valid_mask) > 0:
            c_s_avg = np.median(c_s[valid_mask])
            R = 1.0
            t_cross = R / c_s_avg
            return t_cross, c_s_avg
    
    # Fallback: compute from pressure and density
    rho = snapshot_data[:, 11]  # Column 11: density
    pres = snapshot_data[:, 12]  # Column 12: pressure
    
    # Filter out invalid values (zero or negative density/pressure)
    valid_mask = (rho > 0) & (pres > 0) & np.isfinite(rho) & np.isfinite(pres)
    
    if np.sum(valid_mask) == 0:
        # Last resort: use typical Lane-Emden values
        print("Warning: Could not compute sound speed from data, using default")
        c_s_avg = 0.7  # Typical for Lane-Emden n=1.5
        return 1.0 / c_s_avg, c_s_avg
    
    # Compute sound speed
    c_s = np.sqrt(GAMMA * pres[valid_mask] / rho[valid_mask])
    c_s_avg = np.median(c_s)
    
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


def create_hydrostatic_animation(results_dir, output_dir, append_mode=True):
    """
    Create comprehensive 6-panel animation showing hydrostatic equilibrium test.
    
    Args:
        results_dir: Directory containing snapshot CSV files
        output_dir: Directory to save output files
        append_mode: If True, only process new snapshots and append to existing GIF
    """
    
    # Find all snapshots
    snapshot_files = sorted(glob.glob(os.path.join(results_dir, 'snapshot_*.csv')))
    if len(snapshot_files) == 0:
        print(f"Error: No snapshots found in {results_dir}")
        return
    
    print(f"Found {len(snapshot_files)} total snapshots")
    
    # Load metadata for incremental processing
    metadata = load_gif_metadata(output_dir) if append_mode else {
        'processed_snapshots': [], 
        'metrics': {'times': [], 'rms_errors': [], 'max_velocities': [], 
                   'median_velocities': [], 'K_variations': []}
    }
    
    processed_set = set(metadata['processed_snapshots'])
    new_snapshot_files = [f for f in snapshot_files if os.path.basename(f) not in processed_set]
    
    if append_mode and len(new_snapshot_files) == 0:
        print("No new snapshots to process. GIF is up-to-date.")
        return
    
    if append_mode and len(processed_set) > 0:
        print(f"  Previously processed: {len(processed_set)} snapshots")
        print(f"  New snapshots to add: {len(new_snapshot_files)}")
    
    # Load energy data
    energy_file = os.path.join(results_dir, 'energy.dat')
    energy_data = load_energy_file(energy_file)
    
    # Estimate crossing time from first snapshot
    first_snapshot = load_snapshot(snapshot_files[0])
    t_cross, c_s_avg = compute_crossing_time(first_snapshot)
    print(f"Estimated crossing time: t_cross = {t_cross:.3f} code units")
    print(f"Average sound speed: c_s = {c_s_avg:.3f} code units")
    
    # Quality metrics storage - restore from metadata if appending
    times = metadata['metrics'].get('times', [])
    rms_errors = metadata['metrics'].get('rms_errors', [])
    max_velocities = metadata['metrics'].get('max_velocities', [])
    median_velocities = metadata['metrics'].get('median_velocities', [])
    K_variations = metadata['metrics'].get('K_variations', [])
    
    # Generate frames for new snapshots only
    new_frames = []
    
    def render_frame(snapshot_file, frame_idx, total_frames):
        """Render a single frame and return as PIL Image."""
        # Setup figure
        fig = plt.figure(figsize=(18, 10))
        gs = GridSpec(2, 3, figure=fig, hspace=0.3, wspace=0.3)
        
        ax1 = fig.add_subplot(gs[0, 0])  # 3D density
        ax2 = fig.add_subplot(gs[0, 1])  # XY slice
        ax3 = fig.add_subplot(gs[0, 2])  # Radial density profile
        ax4 = fig.add_subplot(gs[1, 0])  # Velocity histogram
        ax5 = fig.add_subplot(gs[1, 1])  # K(r) profile (polytropic constant)
        ax6 = fig.add_subplot(gs[1, 2])  # Quality metrics
        
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
        if energy_data is not None and frame_idx < len(energy_data):
            current_time = energy_data[frame_idx, 0]
        else:
            current_time = step_num * 0.2  # Assume output frequency
        
        # Compute quality metrics
        metrics = analyze_hydrostatic_quality(data, lane_emden_density)
        times.append(current_time)
        rms_errors.append(metrics['rms_density_error'])
        max_velocities.append(metrics['max_velocity'])
        median_velocities.append(metrics['median_velocity'])
        
        # Panel 1: 3D density scatter
        scatter = ax1.scatter(pos[:, 0], pos[:, 1], c=rho, s=1, cmap='viridis',
                             vmin=0, vmax=RHO_C, alpha=0.6)
        ax1.set_xlabel('x [code units]')
        ax1.set_ylabel('y [code units]')
        ax1.set_title('3D Density Distribution (XY projection)')
        ax1.set_xlim(-1.2, 1.2)
        ax1.set_ylim(-1.2, 1.2)
        ax1.set_aspect('equal')
        plt.colorbar(scatter, ax=ax1, label='Density [code units]')
        
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
        # Ensure xlim upper bound is > 0 to avoid singular transformation
        xlim_upper = max(0.001, min(0.1, metrics['max_velocity'] / c_s_avg * 1.5))
        ax4.set_xlim(0, xlim_upper)
        ax4.grid(True, alpha=0.3)
        
        # Panel 5: K(r) profile (polytropic constant)
        # K = P / rho^gamma should be CONSTANT for isentropic Lane-Emden
        pres = data[:, 12]  # Column 12: pressure
        K_values = pres / rho**GAMMA
        
        # Bin K by radius
        K_binned = []
        K_std_binned = []
        for i in range(len(r_bins) - 1):
            mask_bin = (r >= r_bins[i]) & (r < r_bins[i+1])
            if np.sum(mask_bin) > 0:
                K_binned.append(np.mean(K_values[mask_bin]))
                K_std_binned.append(np.std(K_values[mask_bin]))
            else:
                K_binned.append(np.nan)
                K_std_binned.append(np.nan)
        
        K_binned = np.array(K_binned)
        K_std_binned = np.array(K_std_binned)
        
        # Plot K(r) profile
        ax5.plot(r_centers, K_binned, 'b-o', markersize=4, linewidth=2, label='K(r) SPH')
        ax5.fill_between(r_centers, K_binned - K_std_binned, K_binned + K_std_binned,
                        alpha=0.3, color='blue')
        ax5.axhline(K_EXPECTED, color='red', linestyle='--', linewidth=2,
                   label=f'K expected = {K_EXPECTED:.4f}')
        ax5.set_xlabel('Radius r [code units]')
        ax5.set_ylabel('K = P/ρ^γ')
        
        # Compute K variation metric
        valid_K = K_binned[~np.isnan(K_binned) & (r_centers < 0.85)]
        K_variation = np.std(valid_K) / np.mean(valid_K) * 100 if len(valid_K) > 0 else 0
        K_mean = np.mean(valid_K) if len(valid_K) > 0 else 0
        K_variations.append(K_variation)
        
        ax5.set_title(f'Polytropic Constant K(r)\nK_mean={K_mean:.4f}, Variation={K_variation:.1f}%')
        ax5.legend(loc='best', fontsize=9)
        ax5.grid(True, alpha=0.3)
        ax5.set_xlim(0, 1.0)
        ax5.set_ylim(0, K_EXPECTED * 1.5)
        
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
        
        # Main title - use frame count (len(times)) which accumulates correctly across incremental runs
        # Don't show total since old frames would have stale totals when new snapshots are added
        fig.suptitle(f'Hydrostatic Equilibrium Test (DISPH + Self-Gravity)\n'
                    f'Time: {current_time:.2f} code units = {current_time/t_cross:.2f} t_cross | '
                    f'Frame #{len(times)}',
                    fontsize=14, fontweight='bold')
        
        # Convert figure to PIL Image using buffer (cross-platform compatible)
        # Use tight_layout but NOT bbox_inches='tight' to maintain consistent frame sizes
        fig.tight_layout(rect=[0, 0, 1, 0.95])  # Leave room for suptitle
        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=100)
        buf.seek(0)
        img = Image.open(buf).convert('RGB')
        plt.close(fig)
        return img
    
    # Process new snapshots
    output_file = os.path.join(output_dir, 'hydrostatic_animation.gif')
    
    if append_mode and len(new_snapshot_files) > 0:
        print(f"Processing {len(new_snapshot_files)} new frames (incremental mode)...")
        
        # Generate frames for new snapshots only
        for i, snapshot_file in enumerate(new_snapshot_files):
            # Find the global index of this snapshot
            global_idx = snapshot_files.index(snapshot_file)
            print(f"  Rendering frame {i+1}/{len(new_snapshot_files)}: {os.path.basename(snapshot_file)}")
            
            frame_img = render_frame(snapshot_file, global_idx, len(snapshot_files))
            new_frames.append(frame_img)
            
            # Track processed snapshots
            metadata['processed_snapshots'].append(os.path.basename(snapshot_file))
        
        # Save updated metrics
        metadata['metrics'] = {
            'times': times,
            'rms_errors': rms_errors,
            'max_velocities': max_velocities,
            'median_velocities': median_velocities,
            'K_variations': K_variations
        }
        
        # Append new frames to existing GIF
        total_frames = append_frames_to_gif(output_file, new_frames, output_file, fps=5)
        print(f"✓ Animation updated: {output_file} ({total_frames} total frames)")
        
        # Save metadata for next incremental update
        save_gif_metadata(output_dir, metadata)
        
        # Create/update interactive HTML viewer
        create_interactive_viewer(output_dir, metadata['processed_snapshots'])
        
    else:
        # Full rebuild mode (no existing data or forced rebuild)
        print(f"Creating full animation with {len(snapshot_files)} frames...")
        
        all_frames = []
        for i, snapshot_file in enumerate(snapshot_files):
            print(f"  Rendering frame {i+1}/{len(snapshot_files)}: {os.path.basename(snapshot_file)}")
            frame_img = render_frame(snapshot_file, i, len(snapshot_files))
            all_frames.append(frame_img)
            metadata['processed_snapshots'].append(os.path.basename(snapshot_file))
        
        # Save metrics
        metadata['metrics'] = {
            'times': times,
            'rms_errors': rms_errors,
            'max_velocities': max_velocities,
            'median_velocities': median_velocities,
            'K_variations': K_variations
        }
        
        # Save GIF
        if len(all_frames) > 0:
            all_frames[0].save(
                output_file,
                save_all=True,
                append_images=all_frames[1:],
                duration=200,
                loop=0
            )
            print(f"✓ Animation saved: {output_file}")
        
        # Save metadata
        save_gif_metadata(output_dir, metadata)
    
    # Create interactive HTML viewer with slider
    create_interactive_viewer(output_dir, metadata['processed_snapshots'])
    
    # Create final summary plot
    create_summary_plot(snapshot_files[-1], times, rms_errors, max_velocities,
                       median_velocities, t_cross, c_s_avg, output_dir)


def create_interactive_viewer(output_dir, snapshot_names):
    """
    Create an interactive HTML viewer with a slider to navigate frames.
    
    Args:
        output_dir: Directory containing the frame images
        snapshot_names: List of snapshot filenames that were processed
    """
    frames_dir = os.path.join(output_dir, 'frames')
    os.makedirs(frames_dir, exist_ok=True)
    
    # Check if we need to extract frames from the GIF
    gif_path = os.path.join(output_dir, 'hydrostatic_animation.gif')
    frame_files = sorted(glob.glob(os.path.join(frames_dir, 'frame_*.png')))
    
    # Extract frames from GIF if frames directory is empty or has fewer frames
    if os.path.exists(gif_path):
        try:
            gif = Image.open(gif_path)
            n_gif_frames = 0
            try:
                while True:
                    n_gif_frames += 1
                    gif.seek(gif.tell() + 1)
            except EOFError:
                pass
            
            # Re-extract if frame count doesn't match
            if len(frame_files) != n_gif_frames:
                print(f"Extracting {n_gif_frames} frames for interactive viewer...")
                gif = Image.open(gif_path)
                frame_idx = 0
                try:
                    while True:
                        frame_path = os.path.join(frames_dir, f'frame_{frame_idx:04d}.png')
                        gif.save(frame_path)
                        frame_idx += 1
                        gif.seek(gif.tell() + 1)
                except EOFError:
                    pass
                frame_files = sorted(glob.glob(os.path.join(frames_dir, 'frame_*.png')))
                print(f"  Extracted {len(frame_files)} frames")
        except Exception as e:
            print(f"Warning: Could not extract frames from GIF: {e}")
            return
    
    if len(frame_files) == 0:
        print("Warning: No frames available for interactive viewer")
        return
    
    # Generate HTML with embedded JavaScript slider
    html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Hydrostatic Equilibrium Test - Interactive Viewer</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, sans-serif;
            background: #1a1a2e;
            color: #eee;
            min-height: 100vh;
            padding: 20px;
        }}
        .container {{
            max-width: 1400px;
            margin: 0 auto;
        }}
        h1 {{
            text-align: center;
            margin-bottom: 20px;
            color: #4fc3f7;
        }}
        .viewer {{
            background: #16213e;
            border-radius: 12px;
            padding: 20px;
            box-shadow: 0 4px 20px rgba(0,0,0,0.3);
        }}
        .frame-container {{
            text-align: center;
            margin-bottom: 20px;
        }}
        .frame-container img {{
            max-width: 100%;
            height: auto;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.2);
        }}
        .controls {{
            background: #0f3460;
            padding: 20px;
            border-radius: 8px;
        }}
        .slider-container {{
            display: flex;
            align-items: center;
            gap: 15px;
            margin-bottom: 15px;
        }}
        .slider-container label {{
            min-width: 80px;
            font-weight: 500;
        }}
        .slider-container input[type="range"] {{
            flex: 1;
            height: 8px;
            -webkit-appearance: none;
            background: #1a1a2e;
            border-radius: 4px;
            outline: none;
        }}
        .slider-container input[type="range"]::-webkit-slider-thumb {{
            -webkit-appearance: none;
            width: 20px;
            height: 20px;
            background: #4fc3f7;
            border-radius: 50%;
            cursor: pointer;
            transition: background 0.2s;
        }}
        .slider-container input[type="range"]::-webkit-slider-thumb:hover {{
            background: #81d4fa;
        }}
        .frame-info {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            flex-wrap: wrap;
            gap: 10px;
        }}
        .frame-number {{
            font-size: 1.2em;
            font-weight: bold;
            color: #4fc3f7;
        }}
        .playback-controls {{
            display: flex;
            gap: 10px;
        }}
        .playback-controls button {{
            padding: 10px 20px;
            font-size: 1em;
            border: none;
            border-radius: 6px;
            cursor: pointer;
            transition: all 0.2s;
            font-weight: 500;
        }}
        .btn-play {{
            background: #4caf50;
            color: white;
        }}
        .btn-play:hover {{
            background: #66bb6a;
        }}
        .btn-pause {{
            background: #ff9800;
            color: white;
        }}
        .btn-pause:hover {{
            background: #ffb74d;
        }}
        .btn-step {{
            background: #2196f3;
            color: white;
        }}
        .btn-step:hover {{
            background: #64b5f6;
        }}
        .speed-control {{
            display: flex;
            align-items: center;
            gap: 10px;
        }}
        .speed-control select {{
            padding: 8px 12px;
            font-size: 1em;
            border: none;
            border-radius: 6px;
            background: #1a1a2e;
            color: #eee;
            cursor: pointer;
        }}
        .keyboard-hints {{
            margin-top: 15px;
            padding: 10px;
            background: #1a1a2e;
            border-radius: 6px;
            font-size: 0.9em;
            color: #888;
        }}
        .keyboard-hints kbd {{
            background: #0f3460;
            padding: 2px 8px;
            border-radius: 4px;
            font-family: monospace;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🔬 Hydrostatic Equilibrium Test - Interactive Viewer</h1>
        <div class="viewer">
            <div class="frame-container">
                <img id="frameImage" src="frames/frame_0000.png" alt="Frame">
            </div>
            <div class="controls">
                <div class="slider-container">
                    <label for="frameSlider">Frame:</label>
                    <input type="range" id="frameSlider" min="0" max="{len(frame_files) - 1}" value="0">
                </div>
                <div class="frame-info">
                    <span class="frame-number">Frame: <span id="frameNum">1</span> / {len(frame_files)}</span>
                    <div class="playback-controls">
                        <button class="btn-step" onclick="stepBackward()">⏮ Prev</button>
                        <button class="btn-play" id="playBtn" onclick="togglePlay()">▶ Play</button>
                        <button class="btn-step" onclick="stepForward()">Next ⏭</button>
                    </div>
                    <div class="speed-control">
                        <label for="speedSelect">Speed:</label>
                        <select id="speedSelect" onchange="updateSpeed()">
                            <option value="500">0.5x</option>
                            <option value="200" selected>1x</option>
                            <option value="100">2x</option>
                            <option value="50">4x</option>
                            <option value="25">8x</option>
                        </select>
                    </div>
                </div>
                <div class="keyboard-hints">
                    Keyboard: <kbd>Space</kbd> Play/Pause | <kbd>←</kbd> Previous | <kbd>→</kbd> Next | <kbd>Home</kbd> First | <kbd>End</kbd> Last
                </div>
            </div>
        </div>
    </div>

    <script>
        const totalFrames = {len(frame_files)};
        let currentFrame = 0;
        let isPlaying = false;
        let playInterval = null;
        let playSpeed = 200;

        const frameImage = document.getElementById('frameImage');
        const frameSlider = document.getElementById('frameSlider');
        const frameNum = document.getElementById('frameNum');
        const playBtn = document.getElementById('playBtn');

        // Preload images for smoother playback
        const images = [];
        for (let i = 0; i < totalFrames; i++) {{
            const img = new Image();
            img.src = `frames/frame_${{String(i).padStart(4, '0')}}.png`;
            images.push(img);
        }}

        function updateFrame(frame) {{
            currentFrame = Math.max(0, Math.min(totalFrames - 1, frame));
            frameImage.src = `frames/frame_${{String(currentFrame).padStart(4, '0')}}.png`;
            frameSlider.value = currentFrame;
            frameNum.textContent = currentFrame + 1;
        }}

        function stepForward() {{
            updateFrame(currentFrame + 1);
        }}

        function stepBackward() {{
            updateFrame(currentFrame - 1);
        }}

        function togglePlay() {{
            isPlaying = !isPlaying;
            if (isPlaying) {{
                playBtn.textContent = '⏸ Pause';
                playBtn.className = 'btn-pause';
                playInterval = setInterval(() => {{
                    if (currentFrame >= totalFrames - 1) {{
                        currentFrame = -1;  // Loop back
                    }}
                    stepForward();
                }}, playSpeed);
            }} else {{
                playBtn.textContent = '▶ Play';
                playBtn.className = 'btn-play';
                clearInterval(playInterval);
            }}
        }}

        function updateSpeed() {{
            playSpeed = parseInt(document.getElementById('speedSelect').value);
            if (isPlaying) {{
                clearInterval(playInterval);
                playInterval = setInterval(() => {{
                    if (currentFrame >= totalFrames - 1) {{
                        currentFrame = -1;
                    }}
                    stepForward();
                }}, playSpeed);
            }}
        }}

        frameSlider.addEventListener('input', (e) => {{
            updateFrame(parseInt(e.target.value));
        }});

        // Keyboard controls
        document.addEventListener('keydown', (e) => {{
            switch(e.key) {{
                case ' ':
                    e.preventDefault();
                    togglePlay();
                    break;
                case 'ArrowRight':
                    stepForward();
                    break;
                case 'ArrowLeft':
                    stepBackward();
                    break;
                case 'Home':
                    updateFrame(0);
                    break;
                case 'End':
                    updateFrame(totalFrames - 1);
                    break;
            }}
        }});
    </script>
</body>
</html>
'''
    
    html_path = os.path.join(output_dir, 'interactive_viewer.html')
    with open(html_path, 'w') as f:
        f.write(html_content)
    print(f"✓ Interactive viewer saved: {html_path}")
    print(f"  Open in browser: file://{os.path.abspath(html_path)}")


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


def find_data_directory(results_dir):
    """
    Auto-detect the correct data directory containing snapshots.
    
    Handles cases where:
    1. Snapshots are directly in results_dir
    2. Snapshots are in a 'continue' subdirectory
    3. Snapshots are in the parent GSPH folder when using GSPH_CONTINUE
    4. Need to search sibling directories
    
    Returns:
        Tuple of (data_dir, energy_file_path) or (None, None) if not found
    """
    results_path = Path(results_dir)
    
    # Strategy 1: Check if snapshots exist in the provided directory
    snapshots = list(results_path.glob('snapshot_*.csv'))
    if snapshots:
        energy_file = results_path / 'energy.dat'
        return str(results_path), str(energy_file) if energy_file.exists() else None
    
    # Strategy 2: Check 'continue' subdirectory
    continue_dir = results_path / 'continue'
    if continue_dir.exists():
        snapshots = list(continue_dir.glob('snapshot_*.csv'))
        if snapshots:
            energy_file = continue_dir / 'energy.dat'
            return str(continue_dir), str(energy_file) if energy_file.exists() else None
    
    # Strategy 3: Handle *_CONTINUE pattern - look for base directory without _CONTINUE
    dir_name = results_path.name
    if dir_name.endswith('_CONTINUE') or dir_name.endswith('_continue'):
        base_name = dir_name.rsplit('_', 1)[0]  # Remove _CONTINUE suffix
        base_dir = results_path.parent / base_name
        
        if base_dir.exists():
            # Check if data is in base_dir directly
            snapshots = list(base_dir.glob('snapshot_*.csv'))
            if snapshots:
                energy_file = base_dir / 'energy.dat'
                return str(base_dir), str(energy_file) if energy_file.exists() else None
            
            # Check if data is in base_dir/continue
            continue_subdir = base_dir / 'continue'
            if continue_subdir.exists():
                snapshots = list(continue_subdir.glob('snapshot_*.csv'))
                if snapshots:
                    energy_file = continue_subdir / 'energy.dat'
                    return str(continue_subdir), str(energy_file) if energy_file.exists() else None
    
    # Strategy 4: Search all sibling directories
    parent_dir = results_path.parent
    if parent_dir.exists():
        for sibling in parent_dir.iterdir():
            if sibling.is_dir():
                snapshots = list(sibling.glob('snapshot_*.csv'))
                if snapshots:
                    print(f"  Found data in sibling directory: {sibling}")
                    energy_file = sibling / 'energy.dat'
                    return str(sibling), str(energy_file) if energy_file.exists() else None
                
                # Also check continue subdirectory
                continue_subdir = sibling / 'continue'
                if continue_subdir.exists():
                    snapshots = list(continue_subdir.glob('snapshot_*.csv'))
                    if snapshots:
                        print(f"  Found data in: {continue_subdir}")
                        energy_file = continue_subdir / 'energy.dat'
                        return str(continue_subdir), str(energy_file) if energy_file.exists() else None
    
    return None, None


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Visualization script for Lane-Emden hydrostatic equilibrium test.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Standard run (incremental by default - only processes new snapshots)
  python3 visualize_hydrostatic.py results/ output/
  
  # Force full rebuild from scratch
  python3 visualize_hydrostatic.py results/ output/ --rebuild
  
  # Explicitly enable incremental mode (default behavior)
  python3 visualize_hydrostatic.py results/ output/ --append
  
The incremental mode stores metadata in .hydrostatic_gif_metadata.json
and only processes new snapshot files, appending them to the existing GIF.
This is much faster for long-running simulations.
        """
    )
    parser.add_argument('results_dir', help='Directory containing snapshot_*.csv and energy.dat')
    parser.add_argument('output_dir', help='Directory to save plots and animations')
    parser.add_argument('--rebuild', action='store_true', 
                       help='Force full rebuild of GIF from scratch (ignore cached metadata)')
    parser.add_argument('--append', action='store_true', default=True,
                       help='Incremental mode: only process new snapshots (default)')
    
    args = parser.parse_args()
    
    results_dir = args.results_dir
    output_dir = args.output_dir
    append_mode = not args.rebuild  # Default to append mode unless --rebuild is specified
    
    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("Hydrostatic Equilibrium Test Visualization")
    print("=" * 60)
    print(f"Results directory: {results_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Mode: {'Incremental (append new frames)' if append_mode else 'Full rebuild'}")
    print("")
    
    # Auto-detect data directory
    data_dir, energy_file = find_data_directory(results_dir)
    
    if data_dir is None:
        print(f"Error: No snapshots found in {results_dir}")
        print("")
        print("Searched locations:")
        print(f"  - {results_dir}")
        print(f"  - {results_dir}/continue")
        dir_name = Path(results_dir).name
        if '_CONTINUE' in dir_name.upper():
            base_name = dir_name.rsplit('_', 1)[0]
            print(f"  - {Path(results_dir).parent / base_name}")
            print(f"  - {Path(results_dir).parent / base_name / 'continue'}")
        print("")
        print("Please check that simulations have been run and snapshots exist.")
        sys.exit(1)
    
    if data_dir != results_dir:
        print(f"Auto-detected data directory: {data_dir}")
        print("")
    
    create_hydrostatic_animation(data_dir, output_dir, append_mode=append_mode)
    
    print("")
    print("=" * 60)
    print("✓ Visualization complete!")
    print("=" * 60)
