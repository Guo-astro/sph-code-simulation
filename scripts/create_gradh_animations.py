#!/usr/bin/env python3
"""
Create GIF animations for all 4 grad-h test scenarios.

Scenarios:
1. GSPH without grad-h (collapses)
2. GSPH with grad-h (stable)
3. SSPH without grad-h (stable)
4. SSPH with grad-h (stable)
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.gridspec import GridSpec
import os
import glob
import re
import pandas as pd

def load_all_snapshots(method_dir):
    """Load all snapshots for a method."""
    if not os.path.exists(method_dir):
        return None
    
    files = sorted(glob.glob(os.path.join(method_dir, "snapshot_*.csv")))
    if not files:
        return None
    
    snapshots = []
    for f in files:
        try:
            # Read time from header
            time_val = None
            with open(f) as fp:
                for line in fp:
                    if 'Time (physical):' in line:
                        match = re.search(r'Time \(physical\):\s*([\d.e+-]+)', line)
                        if match:
                            time_val = float(match.group(1))
                        break
            
            if time_val is None:
                idx = int(re.search(r'snapshot_(\d+)', f).group(1))
                time_val = idx * 0.2
            
            # Read data
            df = pd.read_csv(f, comment='#')
            
            if 'pos_x' in df.columns and 'dens' in df.columns:
                snapshots.append({
                    'time': time_val,
                    'x': df['pos_x'].values,
                    'rho': df['dens'].values,
                    'rho_max': df['dens'].max(),
                    'file': f
                })
        except Exception as e:
            print(f"Error loading {f}: {e}")
            continue
    
    # Sort by time
    snapshots.sort(key=lambda s: s['time'])
    return snapshots


def create_single_method_animation(snapshots, method_name, output_file, fps=10):
    """Create animation for a single method."""
    if not snapshots:
        print(f"No data for {method_name}")
        return
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    
    # Find global limits
    all_x = np.concatenate([s['x'] for s in snapshots])
    x_min, x_max = all_x.min(), all_x.max()
    
    # For density, use log scale for collapsing case
    rho_max_all = max(s['rho_max'] for s in snapshots)
    use_log = rho_max_all > 100
    
    if use_log:
        rho_min, rho_max = 0.1, rho_max_all * 1.5
    else:
        rho_min, rho_max = 0, min(5, rho_max_all * 1.2)
    
    # Setup top panel (density profile)
    line1, = ax1.plot([], [], 'b-', linewidth=1.5)
    scatter1 = ax1.scatter([], [], c='blue', s=10, alpha=0.5)
    
    ax1.set_xlim(x_min - 0.1, x_max + 0.1)
    if use_log:
        ax1.set_yscale('log')
        ax1.set_ylim(rho_min, rho_max)
    else:
        ax1.set_ylim(rho_min, rho_max)
    ax1.set_xlabel('Position x', fontsize=11)
    ax1.set_ylabel('Density ρ', fontsize=11)
    ax1.grid(True, alpha=0.3)
    
    title1 = ax1.set_title('', fontsize=12)
    
    # Setup bottom panel (max density vs time)
    times = [s['time'] for s in snapshots]
    rho_maxes = [s['rho_max'] for s in snapshots]
    
    ax2.plot(times, rho_maxes, 'k-', linewidth=1, alpha=0.3)
    point2, = ax2.plot([], [], 'ro', markersize=10)
    
    ax2.set_xlim(0, max(times) * 1.05)
    if use_log:
        ax2.set_yscale('log')
        ax2.set_ylim(0.5, rho_max_all * 2)
    else:
        ax2.set_ylim(0, min(5, rho_max_all * 1.2))
    ax2.set_xlabel('Time', fontsize=11)
    ax2.set_ylabel('Maximum density ρ_max', fontsize=11)
    ax2.grid(True, alpha=0.3)
    ax2.set_title(f'{method_name}', fontsize=13, fontweight='bold')
    
    def init():
        line1.set_data([], [])
        scatter1.set_offsets(np.empty((0, 2)))
        point2.set_data([], [])
        title1.set_text('')
        return line1, scatter1, point2, title1
    
    def animate(i):
        s = snapshots[i]
        
        # Sort by x for line plot
        idx = np.argsort(s['x'])
        x_sorted = s['x'][idx]
        rho_sorted = s['rho'][idx]
        
        line1.set_data(x_sorted, rho_sorted)
        scatter1.set_offsets(np.column_stack([s['x'], s['rho']]))
        
        point2.set_data([s['time']], [s['rho_max']])
        
        title1.set_text(f't = {s["time"]:.2f}, ρ_max = {s["rho_max"]:.2f}')
        
        return line1, scatter1, point2, title1
    
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=len(snapshots), interval=1000//fps,
                                   blit=True)
    
    # Save
    print(f"Saving {output_file}...")
    anim.save(output_file, writer='pillow', fps=fps)
    plt.close(fig)
    print(f"  Done: {output_file}")


def create_4panel_animation(all_data, output_file, fps=8):
    """Create 4-panel animation comparing all methods."""
    
    methods = ['gsph_nogradh', 'gsph_gradh', 'ssph_nogradh', 'ssph_gradh']
    titles = ['GSPH without grad-h', 'GSPH with grad-h', 
              'SSPH without grad-h', 'SSPH with grad-h']
    colors = ['red', 'blue', 'green', 'purple']
    
    # Find common time points
    all_times = set()
    for method in methods:
        if method in all_data and all_data[method]:
            for s in all_data[method]:
                all_times.add(round(s['time'], 2))
    
    common_times = sorted(all_times)
    
    # Limit frames for reasonable GIF size
    if len(common_times) > 100:
        step = len(common_times) // 100
        common_times = common_times[::step]
    
    print(f"Creating 4-panel animation with {len(common_times)} frames...")
    
    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(2, 2, figure=fig, hspace=0.3, wspace=0.25)
    
    axes = [fig.add_subplot(gs[i//2, i%2]) for i in range(4)]
    
    # Global limits
    x_min, x_max = -1.5, 1.5
    rho_max_global = 1
    for method in methods:
        if method in all_data and all_data[method]:
            rho_max_global = max(rho_max_global, 
                                 max(s['rho_max'] for s in all_data[method]))
    
    # Setup axes
    lines = []
    scatters = []
    title_texts = []
    
    for i, (ax, method, title, color) in enumerate(zip(axes, methods, titles, colors)):
        ax.set_xlim(x_min, x_max)
        
        # Use log scale only for collapsing case
        if method == 'gsph_nogradh':
            ax.set_yscale('log')
            ax.set_ylim(0.1, 10000)
        else:
            ax.set_ylim(0, 4)
        
        ax.set_xlabel('Position x', fontsize=10)
        ax.set_ylabel('Density ρ', fontsize=10)
        ax.grid(True, alpha=0.3)
        
        line, = ax.plot([], [], color=color, linewidth=1.5)
        scatter = ax.scatter([], [], c=color, s=8, alpha=0.4)
        
        lines.append(line)
        scatters.append(scatter)
        
        txt = ax.set_title(title, fontsize=11, fontweight='bold')
        title_texts.append(txt)
    
    # Time display
    time_text = fig.suptitle('', fontsize=14, fontweight='bold')
    
    def get_snapshot_at_time(snapshots, target_time):
        """Find snapshot closest to target time."""
        if not snapshots:
            return None
        
        best = min(snapshots, key=lambda s: abs(s['time'] - target_time))
        if abs(best['time'] - target_time) < 0.5:
            return best
        return None
    
    def init():
        for line in lines:
            line.set_data([], [])
        for scatter in scatters:
            scatter.set_offsets(np.empty((0, 2)))
        time_text.set_text('')
        return lines + scatters + [time_text]
    
    def animate(frame):
        t = common_times[frame]
        
        for i, method in enumerate(methods):
            if method in all_data:
                s = get_snapshot_at_time(all_data[method], t)
                if s:
                    idx = np.argsort(s['x'])
                    x_sorted = s['x'][idx]
                    rho_sorted = s['rho'][idx]
                    
                    lines[i].set_data(x_sorted, rho_sorted)
                    scatters[i].set_offsets(np.column_stack([s['x'], s['rho']]))
                    
                    title_texts[i].set_text(f'{titles[i]}\nρ_max = {s["rho_max"]:.2f}')
        
        time_text.set_text(f't = {t:.2f}')
        
        return lines + scatters + title_texts + [time_text]
    
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=len(common_times), interval=1000//fps,
                                   blit=False)
    
    print(f"Saving {output_file}...")
    anim.save(output_file, writer='pillow', fps=fps)
    plt.close(fig)
    print(f"  Done: {output_file}")


def create_comparison_timeline(all_data, output_file, fps=8):
    """Create animation showing density evolution timeline for all methods."""
    
    methods = ['gsph_nogradh', 'gsph_gradh', 'ssph_nogradh', 'ssph_gradh']
    labels = ['GSPH no-gradh', 'GSPH gradh', 'SSPH no-gradh', 'SSPH gradh']
    colors = ['red', 'blue', 'green', 'purple']
    markers = ['o', '^', 's', 'D']
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Collect all time series
    time_series = {}
    for method in methods:
        if method in all_data and all_data[method]:
            times = [s['time'] for s in all_data[method]]
            rho_maxes = [s['rho_max'] for s in all_data[method]]
            time_series[method] = {'t': times, 'rho': rho_maxes}
    
    # Find time range
    all_times = []
    for method in methods:
        if method in time_series:
            all_times.extend(time_series[method]['t'])
    t_max = max(all_times) if all_times else 10
    
    # Plot background lines
    for method, label, color in zip(methods, labels, colors):
        if method in time_series:
            ax1.plot(time_series[method]['t'], time_series[method]['rho'],
                    color=color, linewidth=1, alpha=0.3)
            ax2.semilogy(time_series[method]['t'], time_series[method]['rho'],
                        color=color, linewidth=1, alpha=0.3)
    
    ax1.set_xlim(0, t_max * 1.05)
    ax1.set_ylim(0, 50)
    ax1.set_xlabel('Time', fontsize=12)
    ax1.set_ylabel('Maximum density ρ_max', fontsize=12)
    ax1.set_title('Linear Scale', fontsize=13)
    ax1.grid(True, alpha=0.3)
    
    ax2.set_xlim(0, t_max * 1.05)
    ax2.set_ylim(0.5, 10000)
    ax2.set_xlabel('Time', fontsize=12)
    ax2.set_ylabel('Maximum density ρ_max (log)', fontsize=12)
    ax2.set_title('Log Scale', fontsize=13)
    ax2.grid(True, alpha=0.3, which='both')
    
    # Moving points
    points1 = []
    points2 = []
    for method, color, marker in zip(methods, colors, markers):
        p1, = ax1.plot([], [], color=color, marker=marker, markersize=10, linestyle='none')
        p2, = ax2.plot([], [], color=color, marker=marker, markersize=10, linestyle='none')
        points1.append(p1)
        points2.append(p2)
    
    # Legend
    ax1.legend(points1, labels, loc='upper left', fontsize=9)
    ax2.legend(points2, labels, loc='upper left', fontsize=9)
    
    time_text = fig.suptitle('', fontsize=14, fontweight='bold')
    
    # Common time points
    common_times = sorted(set(all_times))
    if len(common_times) > 100:
        step = len(common_times) // 100
        common_times = common_times[::step]
    
    def init():
        for p in points1 + points2:
            p.set_data([], [])
        time_text.set_text('')
        return points1 + points2 + [time_text]
    
    def animate(frame):
        t = common_times[frame]
        
        for i, method in enumerate(methods):
            if method in time_series:
                ts = time_series[method]
                # Find closest time
                idx = np.argmin(np.abs(np.array(ts['t']) - t))
                if abs(ts['t'][idx] - t) < 0.5:
                    points1[i].set_data([ts['t'][idx]], [ts['rho'][idx]])
                    points2[i].set_data([ts['t'][idx]], [ts['rho'][idx]])
        
        time_text.set_text(f't = {t:.2f}')
        
        return points1 + points2 + [time_text]
    
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=len(common_times), interval=1000//fps,
                                   blit=False)
    
    print(f"Saving {output_file}...")
    anim.save(output_file, writer='pillow', fps=fps)
    plt.close(fig)
    print(f"  Done: {output_file}")


def main():
    print("=" * 70)
    print("CREATING GIF ANIMATIONS FOR GRAD-H COMPARISON")
    print("=" * 70)
    
    results_dir = "results/gradh_comparison"
    output_dir = results_dir
    os.makedirs(output_dir, exist_ok=True)
    
    # Load all data
    print("\nLoading simulation data...")
    all_data = {}
    
    methods = {
        'gsph_nogradh': 'GSPH without grad-h (COLLAPSES)',
        'gsph_gradh': 'GSPH with grad-h (stable)',
        'ssph_nogradh': 'SSPH without grad-h (stable)',
        'ssph_gradh': 'SSPH with grad-h (stable)'
    }
    
    for method, desc in methods.items():
        method_dir = os.path.join(results_dir, method)
        snapshots = load_all_snapshots(method_dir)
        if snapshots:
            all_data[method] = snapshots
            print(f"  {method}: {len(snapshots)} snapshots loaded")
        else:
            print(f"  {method}: NO DATA FOUND")
    
    if not all_data:
        print("\nERROR: No data found!")
        return
    
    # Create individual animations
    print("\n" + "=" * 70)
    print("CREATING INDIVIDUAL ANIMATIONS")
    print("=" * 70)
    
    for method, desc in methods.items():
        if method in all_data:
            output_file = os.path.join(output_dir, f'{method}_animation.gif')
            create_single_method_animation(all_data[method], desc, output_file, fps=8)
    
    # Create 4-panel comparison
    print("\n" + "=" * 70)
    print("CREATING 4-PANEL COMPARISON")
    print("=" * 70)
    
    output_file = os.path.join(output_dir, 'gradh_4panel_comparison.gif')
    create_4panel_animation(all_data, output_file, fps=6)
    
    # Create timeline comparison
    print("\n" + "=" * 70)
    print("CREATING TIMELINE COMPARISON")
    print("=" * 70)
    
    output_file = os.path.join(output_dir, 'gradh_timeline_comparison.gif')
    create_comparison_timeline(all_data, output_file, fps=6)
    
    print("\n" + "=" * 70)
    print("ALL ANIMATIONS COMPLETE!")
    print("=" * 70)
    print(f"\nOutput files in {output_dir}/:")
    for f in sorted(glob.glob(os.path.join(output_dir, '*.gif'))):
        size_kb = os.path.getsize(f) / 1024
        print(f"  {os.path.basename(f)} ({size_kb:.0f} KB)")


if __name__ == "__main__":
    main()
