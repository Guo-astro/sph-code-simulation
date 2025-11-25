#!/usr/bin/env python3
"""
Create side-by-side animation comparing two SR Sod simulations
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

    c = 1.0
    v = result['vel_x']
    gamma = 1.0 / np.sqrt(1.0 - v**2 / c**2)
    result['N'] = gamma * result['dens']

    return result

def load_snapshots(directory, skip=5):
    """Load all snapshots from directory"""
    pattern = os.path.join(directory, "snapshot_*.csv")
    files = sorted(glob.glob(pattern))

    snapshots = []
    for i, f in enumerate(files):
        if i % skip == 0:
            data = load_snapshot(f)
            if data is not None:
                snapshots.append(data)

    return snapshots

def create_comparison_animation(dir1, dir2, label1, label2, output_file, skip=5):
    """Create side-by-side comparison animation"""

    print(f"Loading snapshots from {dir1}...")
    snapshots1 = load_snapshots(dir1, skip)

    print(f"Loading snapshots from {dir2}...")
    snapshots2 = load_snapshots(dir2, skip)

    print(f"Loaded {len(snapshots1)} and {len(snapshots2)} snapshots")

    # Use minimum number of frames
    n_frames = min(len(snapshots1), len(snapshots2))
    if n_frames == 0:
        print("Error: No snapshots loaded")
        return

    # Set up figure with side-by-side layout
    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(3, 2, figure=fig, hspace=0.3, wspace=0.25)

    # Left column (before)
    ax1 = fig.add_subplot(gs[0, 0])  # Density
    ax2 = fig.add_subplot(gs[1, 0])  # Pressure
    ax3 = fig.add_subplot(gs[2, 0])  # Velocity

    # Right column (after)
    ax4 = fig.add_subplot(gs[0, 1])
    ax5 = fig.add_subplot(gs[1, 1])
    ax6 = fig.add_subplot(gs[2, 1])

    # Find global limits
    all_x1 = np.concatenate([s['pos_x'] for s in snapshots1])
    all_x2 = np.concatenate([s['pos_x'] for s in snapshots2])
    x_min = min(all_x1.min(), all_x2.min())
    x_max = max(all_x1.max(), all_x2.max())

    all_N = np.concatenate([s['N'] for s in snapshots1 + snapshots2])
    all_P = np.concatenate([s['pres'] for s in snapshots1 + snapshots2])
    all_v = np.concatenate([s['vel_x'] for s in snapshots1 + snapshots2])

    N_min, N_max = all_N.min() * 0.9, all_N.max() * 1.1
    P_min, P_max = all_P.min() * 0.9, all_P.max() * 1.1
    v_min, v_max = all_v.min() * 1.1, all_v.max() * 1.1

    # Initialize lines
    line1, = ax1.plot([], [], 'r-', linewidth=2, label=label1)
    line2, = ax2.plot([], [], 'b-', linewidth=2)
    line3, = ax3.plot([], [], 'g-', linewidth=2)

    line4, = ax4.plot([], [], 'r-', linewidth=2, label=label2)
    line5, = ax5.plot([], [], 'b-', linewidth=2)
    line6, = ax6.plot([], [], 'g-', linewidth=2)

    # Set up axes - LEFT
    ax1.set_xlim(x_min, x_max)
    ax1.set_ylim(N_min, N_max)
    ax1.set_ylabel('Lab-frame Density N')
    ax1.set_title(label1)
    ax1.legend(loc='upper right')
    ax1.grid(True, alpha=0.3)

    ax2.set_xlim(x_min, x_max)
    ax2.set_ylim(P_min, P_max)
    ax2.set_ylabel('Pressure P')
    ax2.grid(True, alpha=0.3)

    ax3.set_xlim(x_min, x_max)
    ax3.set_ylim(v_min, v_max)
    ax3.axhline(0, color='k', linestyle='--', alpha=0.3)
    ax3.set_xlabel('Position x')
    ax3.set_ylabel('Velocity v')
    ax3.grid(True, alpha=0.3)

    # Set up axes - RIGHT
    ax4.set_xlim(x_min, x_max)
    ax4.set_ylim(N_min, N_max)
    ax4.set_ylabel('Lab-frame Density N')
    ax4.set_title(label2)
    ax4.legend(loc='upper right')
    ax4.grid(True, alpha=0.3)

    ax5.set_xlim(x_min, x_max)
    ax5.set_ylim(P_min, P_max)
    ax5.set_ylabel('Pressure P')
    ax5.grid(True, alpha=0.3)

    ax6.set_xlim(x_min, x_max)
    ax6.set_ylim(v_min, v_max)
    ax6.axhline(0, color='k', linestyle='--', alpha=0.3)
    ax6.set_xlabel('Position x')
    ax6.set_ylabel('Velocity v')
    ax6.grid(True, alpha=0.3)

    time_text = fig.text(0.5, 0.96, '', ha='center', fontsize=14, weight='bold')

    def init():
        line1.set_data([], [])
        line2.set_data([], [])
        line3.set_data([], [])
        line4.set_data([], [])
        line5.set_data([], [])
        line6.set_data([], [])
        time_text.set_text('')
        return line1, line2, line3, line4, line5, line6, time_text

    def animate(i):
        # Interpolate index for second dataset if different lengths
        i2 = int(i * len(snapshots2) / len(snapshots1))
        i2 = min(i2, len(snapshots2) - 1)

        data1 = snapshots1[i]
        data2 = snapshots2[i2]

        line1.set_data(data1['pos_x'], data1['N'])
        line2.set_data(data1['pos_x'], data1['pres'])
        line3.set_data(data1['pos_x'], data1['vel_x'])

        line4.set_data(data2['pos_x'], data2['N'])
        line5.set_data(data2['pos_x'], data2['pres'])
        line6.set_data(data2['pos_x'], data2['vel_x'])

        time_text.set_text(f'SR Sod Comparison: t = {data1["time"]:.4f}')

        return line1, line2, line3, line4, line5, line6, time_text

    print(f"Creating animation with {n_frames} frames...")
    anim = animation.FuncAnimation(
        fig, animate, init_func=init,
        frames=n_frames, interval=100, blit=True
    )

    print(f"Saving to {output_file}...")
    if output_file.endswith('.mp4'):
        anim.save(output_file, writer='ffmpeg', fps=10, dpi=100)
    elif output_file.endswith('.gif'):
        anim.save(output_file, writer='pillow', fps=10)
    else:
        anim.save(output_file + '.mp4', writer='ffmpeg', fps=10, dpi=100)

    print(f"Animation saved: {output_file}")
    plt.close()

def main():
    if len(sys.argv) < 3:
        print("Usage: python compare_animations.py <dir1> <dir2> [output.mp4] [skip]")
        print("\nExample:")
        print("  python compare_animations.py \\")
        print("    sample/sr_sod/results/sharp \\")
        print("    sample/sr_sod/results/sharp_fixed \\")
        print("    comparison.mp4 5")
        sys.exit(1)

    dir1 = sys.argv[1]
    dir2 = sys.argv[2]

    output_file = sys.argv[3] if len(sys.argv) > 3 else "comparison_animation.mp4"
    skip = int(sys.argv[4]) if len(sys.argv) > 4 else 5

    label1 = "Before: 100p, C_cd=0.2, 1st order"
    label2 = "After: 400p, C_cd=1.0, 2nd order"

    create_comparison_animation(dir1, dir2, label1, label2, output_file, skip)

if __name__ == "__main__":
    main()
