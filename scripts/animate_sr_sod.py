#!/usr/bin/env python3
"""
Create animation of SR Sod shock tube evolution
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.gridspec import GridSpec
import sys
import os
import glob

def load_snapshot(filename):
    """Load CSV snapshot"""
    data = []
    time = 0.0

    with open(filename, 'r') as f:
        for line in f:
            if '# Time (code):' in line:
                time = float(line.split(':')[1].strip())
            if line.startswith('#') or line.strip().startswith('id'):
                continue
            parts = line.strip().split(',')
            if len(parts) >= 15:
                try:
                    data.append([float(x) for x in parts])
                except ValueError:
                    continue

    if len(data) == 0:
        return None

    data = np.array(data)
    result = {
        'time': time,
        'pos_x': data[:, 1],
        'vel_x': data[:, 4],
        'dens': data[:, 11],
        'pres': data[:, 12],
    }

    idx = np.argsort(result['pos_x'])
    for key in ['pos_x', 'vel_x', 'dens', 'pres']:
        result[key] = result[key][idx]

    # Compute lab-frame density
    c = 1.0
    v = result['vel_x']
    gamma = 1.0 / np.sqrt(1.0 - v**2 / c**2)
    result['N'] = gamma * result['dens']

    return result

def create_animation(snapshot_dir, output_file, skip=5):
    """Create animation from snapshots"""

    # Find all snapshots
    pattern = os.path.join(snapshot_dir, "snapshot_*.csv")
    files = sorted(glob.glob(pattern))

    if len(files) == 0:
        print(f"Error: No snapshots found in {snapshot_dir}")
        return

    print(f"Found {len(files)} snapshots")

    # Load every N-th snapshot
    snapshots = []
    for i, f in enumerate(files):
        if i % skip == 0:
            data = load_snapshot(f)
            if data is not None:
                snapshots.append(data)

    print(f"Using {len(snapshots)} snapshots (skip={skip})")

    if len(snapshots) == 0:
        print("Error: Could not load any snapshots")
        return

    # Set up figure
    fig = plt.figure(figsize=(12, 10))
    gs = GridSpec(3, 1, figure=fig, hspace=0.3)

    ax1 = fig.add_subplot(gs[0])  # Lab-frame density
    ax2 = fig.add_subplot(gs[1])  # Pressure
    ax3 = fig.add_subplot(gs[2])  # Velocity

    # Find global limits
    all_N = np.concatenate([s['N'] for s in snapshots])
    all_P = np.concatenate([s['pres'] for s in snapshots])
    all_v = np.concatenate([s['vel_x'] for s in snapshots])
    all_x = np.concatenate([s['pos_x'] for s in snapshots])

    x_min, x_max = all_x.min(), all_x.max()
    N_min, N_max = all_N.min() * 0.9, all_N.max() * 1.1
    P_min, P_max = all_P.min() * 0.9, all_P.max() * 1.1
    v_min, v_max = all_v.min() * 1.1, all_v.max() * 1.1

    # Initialize lines
    line1, = ax1.plot([], [], 'r-o', markersize=3, linewidth=1.5)
    line2, = ax2.plot([], [], 'b-o', markersize=3, linewidth=1.5)
    line3, = ax3.plot([], [], 'g-o', markersize=3, linewidth=1.5)

    # Set limits and labels
    ax1.set_xlim(x_min, x_max)
    ax1.set_ylim(N_min, N_max)
    ax1.set_xlabel('Position x')
    ax1.set_ylabel('Lab-frame Density N')
    ax1.grid(True, alpha=0.3)

    ax2.set_xlim(x_min, x_max)
    ax2.set_ylim(P_min, P_max)
    ax2.set_xlabel('Position x')
    ax2.set_ylabel('Pressure P')
    ax2.grid(True, alpha=0.3)

    ax3.set_xlim(x_min, x_max)
    ax3.set_ylim(v_min, v_max)
    ax3.axhline(0, color='k', linestyle='--', alpha=0.3)
    ax3.set_xlabel('Position x')
    ax3.set_ylabel('Velocity v')
    ax3.grid(True, alpha=0.3)

    time_text = fig.text(0.5, 0.95, '', ha='center', fontsize=12, weight='bold')

    def init():
        line1.set_data([], [])
        line2.set_data([], [])
        line3.set_data([], [])
        time_text.set_text('')
        return line1, line2, line3, time_text

    def animate(i):
        data = snapshots[i]

        line1.set_data(data['pos_x'], data['N'])
        line2.set_data(data['pos_x'], data['pres'])
        line3.set_data(data['pos_x'], data['vel_x'])

        time_text.set_text(f'SR Sod Shock Tube: t = {data["time"]:.4f}')

        return line1, line2, line3, time_text

    # Create animation
    print(f"Creating animation with {len(snapshots)} frames...")
    anim = animation.FuncAnimation(
        fig, animate, init_func=init,
        frames=len(snapshots), interval=100, blit=True
    )

    # Save animation
    print(f"Saving to {output_file}...")
    if output_file.endswith('.mp4'):
        anim.save(output_file, writer='ffmpeg', fps=10, dpi=100)
    elif output_file.endswith('.gif'):
        anim.save(output_file, writer='pillow', fps=10)
    else:
        print("Warning: Unknown format, saving as MP4")
        anim.save(output_file + '.mp4', writer='ffmpeg', fps=10, dpi=100)

    print(f"Animation saved: {output_file}")
    plt.close()

def main():
    if len(sys.argv) < 2:
        print("Usage: python animate_sr_sod.py <snapshot_directory> [output.mp4] [skip]")
        print("\nExample:")
        print("  python animate_sr_sod.py sample/sr_sod/results/sharp sr_sod.mp4 5")
        print("\nFormats: .mp4 (requires ffmpeg) or .gif")
        sys.exit(1)

    snapshot_dir = sys.argv[1]

    if len(sys.argv) > 2:
        output_file = sys.argv[2]
    else:
        output_file = "sr_sod_animation.mp4"

    skip = 5
    if len(sys.argv) > 3:
        skip = int(sys.argv[3])

    create_animation(snapshot_dir, output_file, skip)

if __name__ == "__main__":
    main()
