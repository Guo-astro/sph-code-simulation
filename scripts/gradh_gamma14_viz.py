#!/usr/bin/env python3
"""
Visualization script for γ=1.4 grad-h comparison test suite.
Generates comparison plots and animated GIF.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from pathlib import Path
import glob

# Directory paths
BASE_DIR = Path("results/gradh_gamma14")
OUTPUT_DIR = BASE_DIR

# Simulation configs
CONFIGS = [
    ("gsph_gradh", "GSPH + grad-h", "blue"),
    ("ssph_gradh", "SSPH + grad-h", "green"),
    ("ssph_nogradh", "SSPH (no grad-h)", "orange"),
    ("gsph_nogradh", "GSPH (no grad-h)", "red"),
]

def lane_emden_planar(n: float):
    """
    Solve 1D planar Lane-Emden equation: d²θ/dξ² = -θ^n
    
    Args:
        n: Polytropic index
        
    Returns:
        xi: Dimensionless coordinate array
        theta: Dimensionless density (ρ/ρ_c)^(1/n) array
    """
    xi_max = 5.0
    dxi = 0.001
    
    xi_vals = [0.0]
    theta_vals = [1.0]
    
    xi = 0.0
    theta = 1.0
    dtheta = 0.0
    
    while xi < xi_max and theta > 1e-10:
        # RK4 integration
        k1_t = dtheta
        k1_dt = -theta**n if theta > 0 else 0
        
        t2 = theta + 0.5*dxi*k1_t
        dt2 = dtheta + 0.5*dxi*k1_dt
        k2_t = dt2
        k2_dt = -(max(t2, 0)**n)
        
        t3 = theta + 0.5*dxi*k2_t
        dt3 = dtheta + 0.5*dxi*k2_dt
        k3_t = dt3
        k3_dt = -(max(t3, 0)**n)
        
        t4 = theta + dxi*k3_t
        dt4 = dtheta + dxi*k3_dt
        k4_t = dt4
        k4_dt = -(max(t4, 0)**n)
        
        theta += dxi * (k1_t + 2*k2_t + 2*k3_t + k4_t) / 6
        dtheta += dxi * (k1_dt + 2*k2_dt + 2*k3_dt + k4_dt) / 6
        xi += dxi
        
        if theta > 0:
            xi_vals.append(xi)
            theta_vals.append(theta)
    
    return np.array(xi_vals), np.array(theta_vals)


def load_csv_data(filepath: str):
    """
    Load position and density data from CSV snapshot.
    
    Args:
        filepath: Path to CSV file
        
    Returns:
        x: Position array
        rho: Density array
        time: Simulation time
    """
    with open(filepath) as f:
        lines = f.readlines()
    
    header_idx = None
    for i, line in enumerate(lines):
        if line.startswith("id,"):
            header_idx = i
            break
    
    if header_idx is None:
        return None, None, None
    
    x_vals, rho_vals = [], []
    time = 0.0
    
    for line in lines[:header_idx]:
        if "Time (code):" in line:
            time = float(line.split(":")[1].strip())
    
    for line in lines[header_idx+1:]:
        if not line.strip():
            continue
        parts = line.split(",")
        x_vals.append(float(parts[1]))
        rho_vals.append(float(parts[6]))
    
    return np.array(x_vals), np.array(rho_vals), time


def get_snapshots(sim_dir: str):
    """Get sorted list of snapshot files."""
    pattern = str(BASE_DIR / sim_dir / "snapshot_*.csv")
    return sorted(glob.glob(pattern))


def main():
    """Generate all visualizations."""
    print("=" * 60)
    print("γ=1.4 Grad-h Comparison Visualization")
    print("=" * 60)
    
    # Compute Lane-Emden solution (n=2.5 for γ=1.4)
    n = 2.5  # n = 1/(γ-1) = 1/0.4 = 2.5
    xi_le, theta_le = lane_emden_planar(n)
    xi_1 = xi_le[-1]
    print(f"Lane-Emden n={n}: ξ₁ = {xi_1:.4f}")
    
    # Get relaxed IC properties
    relax_file = BASE_DIR / "relaxation" / "snapshot_0006.csv"
    if not relax_file.exists():
        print(f"Error: Relaxation file not found: {relax_file}")
        return
    
    x_relax, rho_relax, _ = load_csv_data(str(relax_file))
    rho_c = rho_relax.max()
    x_max = np.max(np.abs(x_relax))
    print(f"Relaxed IC: ρ_c = {rho_c:.4f}, x_max = {x_max:.4f}")
    
    # Scale Lane-Emden solution
    alpha = x_max / xi_1
    x_analytic = alpha * xi_le
    rho_analytic = rho_c * theta_le**n
    
    # Get frame count
    n_frames = min(len(get_snapshots(c[0])) for c, *_ in [CONFIGS[0]])
    for sim_dir, _, _ in CONFIGS:
        n_frames = min(n_frames, len(get_snapshots(sim_dir)))
    print(f"Using {n_frames} frames")
    
    # Create animation
    print("\nGenerating animation...")
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()
    
    lines = []
    for idx, (sim_dir, title, color) in enumerate(CONFIGS):
        ax = axes[idx]
        ax.set_xlim(-0.8, 0.8)
        ax.set_ylim(0, 2.0)
        ax.set_xlabel('Position x', fontsize=11)
        ax.set_ylabel('Density ρ', fontsize=11)
        ax.set_title(title, fontsize=13, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        # Analytic Lane-Emden
        ax.fill_between(x_analytic, 0, rho_analytic, alpha=0.15, color='green', 
                       label='Lane-Emden n=2.5')
        ax.fill_between(-x_analytic, 0, rho_analytic, alpha=0.15, color='green')
        ax.plot(x_analytic, rho_analytic, 'g--', linewidth=1, alpha=0.7)
        ax.plot(-x_analytic, rho_analytic, 'g--', linewidth=1, alpha=0.7)
        
        # Simulation line
        line, = ax.plot([], [], color=color, linewidth=1.5, alpha=0.8, label='Simulation')
        lines.append(line)
        ax.legend(loc='upper right', fontsize=9)
    
    time_text = fig.text(0.5, 0.02, '', ha='center', fontsize=12)
    fig.suptitle('Grad-h Effect: γ=1.4 Polytrope (Relaxed IC)', fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0.04, 1, 0.96])
    
    def init():
        for line in lines:
            line.set_data([], [])
        time_text.set_text('')
        return lines + [time_text]
    
    def animate(frame):
        for idx, (sim_dir, _, _) in enumerate(CONFIGS):
            files = get_snapshots(sim_dir)
            if frame < len(files):
                x, rho, time = load_csv_data(files[frame])
                if x is not None:
                    sort_idx = np.argsort(x)
                    lines[idx].set_data(x[sort_idx], rho[sort_idx])
        
        files = get_snapshots(CONFIGS[0][0])
        if frame < len(files):
            _, _, time = load_csv_data(files[frame])
            time_text.set_text(f't = {time:.2f}')
        
        return lines + [time_text]
    
    anim = animation.FuncAnimation(fig, animate, init_func=init, 
                                   frames=n_frames, interval=100, blit=True)
    anim.save(str(OUTPUT_DIR / 'gradh_gamma14_comparison.gif'), 
              writer='pillow', fps=10, dpi=100)
    print(f"✓ Animation saved: {OUTPUT_DIR / 'gradh_gamma14_comparison.gif'}")
    plt.close()
    
    # Create static comparison plot
    print("\nGenerating static comparison...")
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()
    
    for idx, (sim_dir, title, color) in enumerate(CONFIGS):
        ax = axes[idx]
        
        # Initial state
        x0, rho0, t0 = load_csv_data(str(BASE_DIR / sim_dir / "snapshot_0000.csv"))
        sort_idx = np.argsort(x0)
        ax.plot(x0[sort_idx], rho0[sort_idx], 'k--', linewidth=1.5, alpha=0.5, 
               label=f't={t0:.1f} (initial)')
        
        # Final state
        xf, rhof, tf = load_csv_data(str(BASE_DIR / sim_dir / "snapshot_0100.csv"))
        sort_idx = np.argsort(xf)
        ax.plot(xf[sort_idx], rhof[sort_idx], color=color, linewidth=2, 
               label=f't={tf:.1f} (final)')
        
        # Calculate change
        change = (rhof.max() - rho0.max()) / rho0.max() * 100
        
        ax.set_xlim(-0.8, 0.8)
        ax.set_ylim(0, 2.0)
        ax.set_xlabel('Position x', fontsize=11)
        ax.set_ylabel('Density ρ', fontsize=11)
        ax.set_title(f'{title}\nΔρ_c = {change:+.1f}%', fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper right', fontsize=9)
    
    plt.suptitle('Initial vs Final: γ=1.4 Polytrope (t=20)', fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(str(OUTPUT_DIR / 'gradh_gamma14_initial_vs_final.png'), dpi=150)
    print(f"✓ Saved: {OUTPUT_DIR / 'gradh_gamma14_initial_vs_final.png'}")
    plt.close()
    
    # Summary table
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"{'Method':15} | {'ρ_c(t=0)':>10} | {'ρ_c(t=20)':>10} | {'Change':>10}")
    print("-" * 52)
    
    for sim_dir, title, _ in CONFIGS:
        x0, rho0, _ = load_csv_data(str(BASE_DIR / sim_dir / "snapshot_0000.csv"))
        xf, rhof, _ = load_csv_data(str(BASE_DIR / sim_dir / "snapshot_0100.csv"))
        change = (rhof.max() - rho0.max()) / rho0.max() * 100
        print(f"{sim_dir:15} | {rho0.max():10.4f} | {rhof.max():10.4f} | {change:+9.1f}%")


if __name__ == "__main__":
    main()
