#!/usr/bin/env python3
"""
Visualize Lane-Emden relaxation results
Creates animations and diagnostic plots for IMBH-cloud research setup validation
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import sys
import os
from pathlib import Path

# Add shared module path for SSOT Lane-Emden solution
_script_dir = Path(__file__).parent
_shared_path = _script_dir.parent.parent.parent / "scripts" / "shared"
if _shared_path.exists():
    sys.path.insert(0, str(_shared_path))
    from lane_emden import load_lane_emden_solution as _load_le
    
    def load_lane_emden_solution(filepath='data/lane_emden/n1.5_3d.dat'):
        """Load the exact Lane-Emden n=1.5 solution using shared module (SSOT)"""
        try:
            solution = _load_le(n=1.5, dim=3)
            return {
                'xi': solution.xi,
                'theta': solution.theta,
                'dtheta': solution.dtheta
            }
        except Exception:
            return None
else:
    def load_lane_emden_solution(filepath='data/lane_emden/n1.5_3d.dat'):
        """Load the exact Lane-Emden n=1.5 solution (fallback)"""
        try:
            data = np.loadtxt(filepath, skiprows=4)
            return {
                'xi': data[:, 0],
                'theta': data[:, 1],
                'dtheta': data[:, 2]
            }
        except Exception:
            return None

def load_snapshot(filepath):
    """Load a single CSV snapshot"""
    try:
        # Find the line with "id," header dynamically
        with open(filepath, 'r') as f:
            skiprows = 0
            for line in f:
                if line.startswith('id,'):
                    break
                skiprows += 1
        
        # Skip header lines + 1 to get to data
        data = np.loadtxt(filepath, delimiter=',', skiprows=skiprows+1)
        # Format: id, pos_x, pos_y, pos_z, vel_x, vel_y, vel_z, acc_x, acc_y, acc_z, 
        #         mass, dens, pres, ene, sml, sound, alpha, balsara, gradh, phi, neighbor, is_ghost
        return {
            'id': data[:, 0].astype(int),
            'pos': data[:, 1:4],
            'vel': data[:, 4:7],
            'mass': data[:, 10],
            'rho': data[:, 11],
            'P': data[:, 12],
            'u': data[:, 13],
            'h': data[:, 14]
        }
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return None

def load_energy_file(filepath):
    """Load energy.dat file"""
    try:
        data = []
        with open(filepath, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    values = line.split()
                    if len(values) >= 5:
                        data.append([float(x) for x in values[:5]])
        
        if not data:
            return None
            
        data = np.array(data)
        return {
            'time': data[:, 0],
            'kinetic': data[:, 1],
            'thermal': data[:, 2],
            'potential': data[:, 3],
            'total': data[:, 4]
        }
    except Exception as e:
        print(f"Error loading energy file: {e}")
        return None

def create_relaxation_animation(results_dir, output_dir):
    """Create comprehensive relaxation animation"""
    
    results_path = Path(results_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Find all snapshots
    snapshots = sorted(results_path.glob('snapshot_*.csv'))
    if not snapshots:
        print(f"No snapshots found in {results_dir}")
        return
    
    print(f"Found {len(snapshots)} snapshots")
    
    # Load energy data
    energy_file = results_path / 'energy.dat'
    energy = load_energy_file(energy_file) if energy_file.exists() else None
    
    # Load all snapshots
    print("Loading snapshots...")
    data_list = []
    for snap in snapshots:
        data = load_snapshot(snap)
        if data is not None:
            data_list.append(data)
    
    if not data_list:
        print("No valid snapshot data loaded!")
        return
    
    n_frames = len(data_list)
    print(f"Loaded {n_frames} valid snapshots")
    
    # Calculate global bounds for consistent plotting
    all_pos = np.vstack([d['pos'] for d in data_list])
    pos_min, pos_max = all_pos.min(axis=0), all_pos.max(axis=0)
    pos_range = max(pos_max - pos_min)
    
    all_rho = np.concatenate([d['rho'] for d in data_list])
    rho_min, rho_max = all_rho.min(), all_rho.max()
    
    all_vel_mag = np.concatenate([np.linalg.norm(d['vel'], axis=1) for d in data_list])
    vel_max = all_vel_mag.max()
    
    print(f"Position range: [{pos_min}, {pos_max}]")
    print(f"Density range: [{rho_min:.3f}, {rho_max:.3f}]")
    print(f"Velocity range: [0, {vel_max:.3e}]")
    
    # Create figure with subplots
    fig = plt.figure(figsize=(16, 10))
    
    # 3D density view (top left)
    ax1 = fig.add_subplot(2, 3, 1, projection='3d')
    
    # XY slice (top middle)
    ax2 = fig.add_subplot(2, 3, 2)
    
    # XZ slice (top right)
    ax3 = fig.add_subplot(2, 3, 3)
    
    # Radial density profile (bottom left)
    ax4 = fig.add_subplot(2, 3, 4)
    
    # Velocity evolution (bottom middle)
    ax5 = fig.add_subplot(2, 3, 5)
    
    # Energy conservation (bottom right)
    ax6 = fig.add_subplot(2, 3, 6)
    
    # Create colorbars once and reuse them
    colorbar_ax2 = None
    colorbar_ax3 = None
    
    def init():
        """Initialize animation"""
        return []
    
    def update(frame):
        """Update animation frame"""
        nonlocal colorbar_ax2, colorbar_ax3
        
        data = data_list[frame]
        
        # Clear all axes
        ax1.clear()
        ax2.clear()
        ax3.clear()
        ax4.clear()
        ax5.clear()
        ax6.clear()
        
        pos = data['pos']
        rho = data['rho']
        vel = data['vel']
        vel_mag = np.linalg.norm(vel, axis=1)
        
        # Calculate radial distance from center
        r = np.linalg.norm(pos, axis=1)
        
        # 1. 3D scatter plot (density-colored) - NO colorbar to reduce clutter
        scatter = ax1.scatter(pos[:, 0], pos[:, 1], pos[:, 2], 
                             c=rho, s=2, cmap='viridis', 
                             vmin=rho_min, vmax=rho_max, alpha=0.6)
        ax1.set_xlabel('X [code units]')
        ax1.set_ylabel('Y [code units]')
        ax1.set_zlabel('Z [code units]')
        ax1.set_title(f'3D Density (Step {frame}/{n_frames-1})')
        ax1.set_xlim(pos_min[0], pos_max[0])
        ax1.set_ylim(pos_min[1], pos_max[1])
        ax1.set_zlim(pos_min[2], pos_max[2])
        
        # 2. XY slice (z ≈ 0)
        z_slice = np.abs(pos[:, 2]) < 0.1
        if z_slice.any():
            scatter2 = ax2.scatter(pos[z_slice, 0], pos[z_slice, 1], 
                                  c=rho[z_slice], s=10, cmap='viridis',
                                  vmin=rho_min, vmax=rho_max)
            # Only create colorbar on first frame
            if colorbar_ax2 is None:
                colorbar_ax2 = plt.colorbar(scatter2, ax=ax2, label='ρ [code units]')
        ax2.set_xlabel('X [code units]')
        ax2.set_ylabel('Y [code units]')
        ax2.set_title('XY Slice (|z| < 0.1)')
        ax2.set_aspect('equal')
        ax2.grid(True, alpha=0.3)
        
        # 3. XZ slice (y ≈ 0)
        y_slice = np.abs(pos[:, 1]) < 0.1
        if y_slice.any():
            scatter3 = ax3.scatter(pos[y_slice, 0], pos[y_slice, 2], 
                                  c=rho[y_slice], s=10, cmap='viridis',
                                  vmin=rho_min, vmax=rho_max)
            # Only create colorbar on first frame
            if colorbar_ax3 is None:
                colorbar_ax3 = plt.colorbar(scatter3, ax=ax3, label='ρ [code units]')
        ax3.set_xlabel('X [code units]')
        ax3.set_ylabel('Z [code units]')
        ax3.set_title('XZ Slice (|y| < 0.1)')
        ax3.set_aspect('equal')
        ax3.grid(True, alpha=0.3)
        
        # 4. Radial density profile
        r_bins = np.linspace(0, r.max(), 50)
        r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])
        rho_profile = []
        for i in range(len(r_bins) - 1):
            mask = (r >= r_bins[i]) & (r < r_bins[i+1])
            if mask.any():
                rho_profile.append(np.mean(rho[mask]))
            else:
                rho_profile.append(np.nan)
        
        ax4.plot(r_centers, rho_profile, 'b-', linewidth=2, label='SPH')
        
        # Exact Lane-Emden n=1.5 solution
        lane_emden = load_lane_emden_solution()
        if lane_emden is not None:
            # Lane-Emden parameters (from lane_emden.cpp)
            xi_1 = 3.6537540101
            alpha = 1.0 / xi_1  # R=1, so α = R/ξ₁
            rho_c = 1.43009692  # Central density in code units
            
            # Compute analytic density: ρ(r) = ρ_c × θ(r/α)^1.5
            r_theory = np.linspace(0, r.max(), 200)
            xi_theory = r_theory / alpha
            
            # Interpolate θ(ξ) from exact solution
            theta_theory = np.interp(xi_theory, lane_emden['xi'], lane_emden['theta'], 
                                    left=1.0, right=0.0)
            rho_theory = rho_c * theta_theory**1.5
            
            ax4.plot(r_theory, rho_theory, 'r--', linewidth=1.5, alpha=0.7, 
                    label=f'Lane-Emden n=1.5 (ρ_c={rho_c:.2f})')
        else:
            # Fallback to approximate formula if data file not found
            r_theory = np.linspace(0, r.max(), 100)
            rho_c = 1.43
            xi = r_theory / 0.273691
            rho_theory = rho_c / (1 + xi**2 / 3)**1.5
            ax4.plot(r_theory, rho_theory, 'r--', linewidth=1.5, alpha=0.7, 
                    label='Lane-Emden n=1.5 (approx)')
        
        ax4.set_xlabel('Radius [code units]')
        ax4.set_ylabel('Density [code units]')
        ax4.set_title('Radial Density Profile')
        ax4.grid(True, alpha=0.3)
        ax4.legend()
        ax4.set_xlim(0, r.max())
        
        # 5. Velocity magnitude histogram
        ax5.hist(vel_mag, bins=50, color='orange', alpha=0.7, edgecolor='black')
        ax5.axvline(np.mean(vel_mag), color='red', linestyle='--', linewidth=2, 
                   label=f'Mean: {np.mean(vel_mag):.2e}')
        ax5.axvline(np.median(vel_mag), color='blue', linestyle='--', linewidth=2,
                   label=f'Median: {np.median(vel_mag):.2e}')
        ax5.set_xlabel('Velocity Magnitude [pc/Myr]')
        ax5.set_ylabel('Particle Count')
        ax5.set_title('Velocity Distribution (Should → 0)')
        ax5.legend()
        ax5.grid(True, alpha=0.3)
        if vel_mag.max() > 0:
            ax5.set_yscale('log')
        
        # 6. Energy conservation (if available)
        if energy is not None and len(energy['time']) > 0:
            steps_so_far = min(frame + 1, len(energy['time']))
            ax6.plot(range(steps_so_far), energy['total'][:steps_so_far], 
                    'b-', linewidth=2, label='Total Energy')
            ax6.plot(range(steps_so_far), energy['thermal'][:steps_so_far], 
                    'r--', linewidth=1.5, label='Thermal')
            if energy['kinetic'].max() > 0:
                ax6.plot(range(steps_so_far), energy['kinetic'][:steps_so_far], 
                        'g:', linewidth=1.5, label='Kinetic')
            
            # Calculate energy drift
            if steps_so_far > 1:
                E0 = energy['total'][0]
                E_current = energy['total'][steps_so_far-1]
                drift = abs(E_current - E0) / abs(E0) * 100 if E0 != 0 else 0
                ax6.text(0.05, 0.95, f'ΔE/E₀ = {drift:.3f}%', 
                        transform=ax6.transAxes, fontsize=10,
                        verticalalignment='top',
                        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        ax6.set_xlabel('Step')
        ax6.set_ylabel('Energy [code units]')
        ax6.set_title('Energy Conservation')
        ax6.legend()
        ax6.grid(True, alpha=0.3)
        
        # Overall title
        fig.suptitle(f'Lane-Emden Relaxation: Step {frame}/{n_frames-1} | N={len(pos)} particles', 
                    fontsize=14, fontweight='bold')
        
        plt.tight_layout()
        return []
    
    # Create animation
    print("Creating animation...")
    anim = FuncAnimation(fig, update, frames=n_frames, init_func=init, 
                        blit=False, interval=200, repeat=True)
    
    # Save as GIF
    output_file = output_path / 'relaxation_animation.gif'
    print(f"Saving animation to {output_file}...")
    writer = PillowWriter(fps=5)
    anim.save(output_file, writer=writer, dpi=100)
    print(f"✓ Animation saved: {output_file}")
    
    # Also save final frame as static image
    update(n_frames - 1)
    summary_file = output_path / 'relaxation_summary.png'
    plt.savefig(summary_file, dpi=150, bbox_inches='tight')
    print(f"✓ Summary plot saved: {summary_file}")
    
    plt.close()
    
    # Create velocity decay plot
    print("Creating velocity decay analysis...")
    fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    mean_vel = []
    max_vel = []
    for data in data_list:
        vel_mag = np.linalg.norm(data['vel'], axis=1)
        mean_vel.append(np.mean(vel_mag))
        max_vel.append(np.max(vel_mag))
    
    steps = np.arange(len(mean_vel))
    
    ax1.semilogy(steps, mean_vel, 'b-o', linewidth=2, markersize=4, label='Mean velocity')
    ax1.semilogy(steps, max_vel, 'r-^', linewidth=2, markersize=4, label='Max velocity')
    ax1.set_xlabel('Relaxation Step')
    ax1.set_ylabel('Velocity Magnitude [pc/Myr]')
    ax1.set_title('Velocity Decay (Should → 0)')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Density profile comparison (first vs last)
    r_initial = np.linalg.norm(data_list[0]['pos'], axis=1)
    r_final = np.linalg.norm(data_list[-1]['pos'], axis=1)
    
    ax2.hist(r_initial, bins=50, alpha=0.5, color='blue', label='Initial', density=True)
    ax2.hist(r_final, bins=50, alpha=0.5, color='red', label='Final', density=True)
    ax2.set_xlabel('Radius [pc]')
    ax2.set_ylabel('Normalized Particle Count')
    ax2.set_title('Radial Distribution (Should remain constant)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    fig2.suptitle(f'Relaxation Quality Metrics (N={len(data_list[0]["pos"])} particles)', 
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    velocity_file = output_path / 'velocity_decay.png'
    plt.savefig(velocity_file, dpi=150, bbox_inches='tight')
    print(f"✓ Velocity analysis saved: {velocity_file}")
    plt.close()
    
    # Print summary statistics
    print("\n" + "="*60)
    print("RELAXATION SUMMARY")
    print("="*60)
    print(f"Number of particles: {len(data_list[0]['pos'])}")
    print(f"Number of steps: {len(data_list)}")
    print(f"Initial mean velocity: {mean_vel[0]:.3e} pc/Myr")
    print(f"Final mean velocity: {mean_vel[-1]:.3e} pc/Myr")
    print(f"Velocity reduction: {(1 - mean_vel[-1]/mean_vel[0])*100:.2f}%")
    
    if energy is not None and len(energy['time']) > 0:
        E0 = energy['total'][0]
        Ef = energy['total'][-1]
        drift = abs(Ef - E0) / abs(E0) * 100 if E0 != 0 else 0
        print(f"\nEnergy conservation:")
        print(f"  Initial E_total: {E0:.6f}")
        print(f"  Final E_total: {Ef:.6f}")
        print(f"  Drift: {drift:.4f}%")
    
    print("\n✓ Relaxation appears " + ("GOOD" if mean_vel[-1]/mean_vel[0] < 0.1 else "NEEDS MORE STEPS"))
    print("="*60)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python3 visualize_relaxation.py <results_dir> [output_dir]")
        print("Example: python3 visualize_relaxation.py simulations/astrophysics/imbh_cloud/results/lane_emden_50k_relax")
        sys.exit(1)
    
    results_dir = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else results_dir
    
    create_relaxation_animation(results_dir, output_dir)
