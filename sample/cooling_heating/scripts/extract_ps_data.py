#!/usr/bin/env python3
"""
Extract data from PostScript figure files.
Parse the actual curve coordinates from f1a.ps
"""

import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def parse_ps_file(filename):
    """Extract moveto/lineto coordinates from PostScript file."""
    
    with open(filename, 'r') as f:
        content = f.read()
    
    # Find all moveto (M) and lineto (L) commands
    # Format: "x y M" or "x y L"
    pattern = r'(\d+)\s+(\d+)\s+([ML])'
    matches = re.findall(pattern, content)
    
    print(f"Found {len(matches)} drawing commands in {filename}")
    
    # Separate into different paths (each M starts a new path)
    paths = []
    current_path = []
    
    for x, y, cmd in matches:
        if cmd == 'M':
            if current_path:
                paths.append(current_path)
            current_path = [(int(x), int(y))]
        else:  # L
            current_path.append((int(x), int(y)))
    
    if current_path:
        paths.append(current_path)
    
    print(f"Extracted {len(paths)} paths")
    
    # Find the longest paths (likely the data curves)
    path_lengths = [(i, len(p)) for i, p in enumerate(paths)]
    path_lengths.sort(key=lambda x: x[1], reverse=True)
    
    print("\nLongest paths:")
    for i, length in path_lengths[:10]:
        print(f"  Path {i}: {length} points")
        if length > 5:
            path = paths[i]
            print(f"    Start: {path[0]}, End: {path[-1]}")
            print(f"    X range: {min(p[0] for p in path)} - {max(p[0] for p in path)}")
            print(f"    Y range: {min(p[1] for p in path)} - {max(p[1] for p in path)}")
    
    return paths

def plot_ps_coordinates(paths, output_file):
    """Plot the extracted coordinates to visualize."""
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Plot each path with different color
    colors = plt.cm.tab10(np.linspace(0, 1, min(len(paths), 10)))
    
    for i, path in enumerate(paths[:20]):  # Plot first 20 paths
        if len(path) > 5:  # Only plot substantial paths
            xs = [p[0] for p in path]
            ys = [p[1] for p in path]
            ax.plot(xs, ys, 'o-', color=colors[i % 10], 
                   linewidth=2, markersize=3, alpha=0.7,
                   label=f'Path {i} ({len(path)} pts)')
    
    ax.set_xlabel('PostScript X coordinate')
    ax.set_ylabel('PostScript Y coordinate')
    ax.set_title('Extracted Paths from f1a.ps')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=200, bbox_inches='tight')
    print(f"\n✓ Plot saved: {output_file}")

def main():
    ps_file = '/Users/guo/Downloads/sphcode/docs/papers/cooling-heating/f1a.ps'
    output_plot = '../results/ps_extraction.png'
    
    print("=" * 70)
    print("Extracting Data from PostScript Figure")
    print("=" * 70)
    print(f"Input: {ps_file}")
    
    paths = parse_ps_file(ps_file)
    plot_ps_coordinates(paths, output_plot)
    
    # Save the longest path (likely the main T(n) curve)
    if paths:
        longest_idx = max(range(len(paths)), key=lambda i: len(paths[i]))
        longest_path = paths[longest_idx]
        
        print(f"\nLongest path (index {longest_idx}):")
        print(f"  {len(longest_path)} points")
        print(f"  First 10 points:")
        for i, (x, y) in enumerate(longest_path[:10]):
            print(f"    {i}: ({x}, {y})")
    
    print("\n" + "=" * 70)

if __name__ == '__main__':
    main()
