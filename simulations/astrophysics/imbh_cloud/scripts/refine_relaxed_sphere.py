#!/usr/bin/env python3
"""
Refine a relaxed sphere by adding particles at midpoints between neighbors.

Key insight: The relaxed particles have OPTIMIZED local arrangement.
Instead of random sampling (which creates clustering), we PRESERVE
the good arrangement by adding particles at midpoints.

Algorithm:
1. Load relaxed particles
2. Build KD-tree for neighbor finding
3. For each particle, find k nearest neighbors
4. Add new particles at midpoints (avoiding duplicates)
5. Adjust masses to preserve total mass and density profile

This creates ~3-5x more particles with inherited good arrangement,
requiring only ~10-20 relaxation steps instead of 100+.

Usage:
    python refine_relaxed_sphere.py <template.h5> <output.h5> [--neighbors 6]
"""

import numpy as np
import h5py
import sys
from scipy.spatial import cKDTree


def load_relaxed(filename):
    """Load relaxed particle configuration."""
    with h5py.File(filename, 'r') as f:
        x = np.array(f['particles/pos_x'])
        y = np.array(f['particles/pos_y'])
        z = np.array(f['particles/pos_z'])
        rho = np.array(f['particles/dens'])
        pres = np.array(f['particles/pres'])
        ene = np.array(f['particles/ene'])
        sml = np.array(f['particles/sml'])
        mass = np.array(f['particles/mass'])
        is_ghost = np.array(f['particles/is_ghost'])

    # Filter real particles only
    real = is_ghost == 0
    pos = np.column_stack([x[real], y[real], z[real]])

    return {
        'pos': pos,
        'rho': rho[real],
        'pres': pres[real],
        'ene': ene[real],
        'sml': sml[real],
        'mass': mass[real][0],
        'n_particles': len(pos)
    }


def refine_particles(data, n_neighbors=6):
    """
    Refine particle distribution by adding midpoint particles.

    Args:
        data: Dictionary with particle data
        n_neighbors: Number of nearest neighbors to use for midpoints

    Returns:
        New positions array
    """
    pos = data['pos']
    n_orig = len(pos)

    print(f"Building KD-tree for {n_orig:,} particles...")
    tree = cKDTree(pos)

    # Find k+1 neighbors (includes self)
    print(f"Finding {n_neighbors} neighbors per particle...")
    distances, indices = tree.query(pos, k=n_neighbors+1)

    # Collect midpoint positions (use set to avoid duplicates)
    print("Generating midpoint particles...")
    midpoints_set = set()

    for i in range(n_orig):
        for j_idx in range(1, n_neighbors+1):  # Skip self (index 0)
            j = indices[i, j_idx]

            # Create canonical edge key to avoid duplicates
            edge = (min(i, j), max(i, j))
            if edge not in midpoints_set:
                midpoints_set.add(edge)

    print(f"  Found {len(midpoints_set):,} unique edges")

    # Compute midpoint positions
    midpoints = []
    for i, j in midpoints_set:
        mid = 0.5 * (pos[i] + pos[j])
        midpoints.append(mid)

    midpoints = np.array(midpoints)

    # Combine original + midpoints
    new_pos = np.vstack([pos, midpoints])

    print(f"  Original: {n_orig:,}")
    print(f"  Midpoints: {len(midpoints):,}")
    print(f"  Total: {len(new_pos):,}")
    print(f"  Refinement factor: {len(new_pos)/n_orig:.2f}x")

    return new_pos, data


def interpolate_quantities(new_pos, data):
    """Interpolate SPH quantities to new particle positions."""
    old_pos = data['pos']

    # Build tree for interpolation
    tree = cKDTree(old_pos)

    # For each new particle, find nearest old particle and copy quantities
    distances, indices = tree.query(new_pos, k=1)

    # Compute radii
    r_new = np.sqrt(new_pos[:, 0]**2 + new_pos[:, 1]**2 + new_pos[:, 2]**2)
    r_old = np.sqrt(old_pos[:, 0]**2 + old_pos[:, 1]**2 + old_pos[:, 2]**2)

    # Interpolate based on radius (for spherically symmetric profile)
    sort_idx = np.argsort(r_old)
    r_sorted = r_old[sort_idx]

    rho_new = np.interp(r_new, r_sorted, data['rho'][sort_idx])
    pres_new = np.interp(r_new, r_sorted, data['pres'][sort_idx])
    ene_new = np.interp(r_new, r_sorted, data['ene'][sort_idx])
    sml_new = np.interp(r_new, r_sorted, data['sml'][sort_idx])

    # Scale smoothing length for higher resolution
    scale = (len(new_pos) / len(old_pos)) ** (1/3)
    sml_new = sml_new / scale

    return rho_new, pres_new, ene_new, sml_new


def write_snapshot(filename, pos, rho, pres, ene, sml, total_mass, template_file):
    """Write refined configuration to HDF5."""
    n = len(pos)
    mass_per_particle = total_mass / n

    with h5py.File(template_file, 'r') as f_in:
        with h5py.File(filename, 'w') as f:
            # Copy metadata
            if 'metadata' in f_in:
                f_in.copy('metadata', f)

            p = f.create_group('particles')

            # Positions
            p.create_dataset('pos_x', data=pos[:, 0])
            p.create_dataset('pos_y', data=pos[:, 1])
            p.create_dataset('pos_z', data=pos[:, 2])

            # Velocities (zero)
            p.create_dataset('vel_x', data=np.zeros(n))
            p.create_dataset('vel_y', data=np.zeros(n))
            p.create_dataset('vel_z', data=np.zeros(n))
            p.create_dataset('vel_t', data=np.zeros(n))

            # Accelerations (zero)
            p.create_dataset('acc_x', data=np.zeros(n))
            p.create_dataset('acc_y', data=np.zeros(n))
            p.create_dataset('acc_z', data=np.zeros(n))

            # SPH quantities
            p.create_dataset('mass', data=np.full(n, mass_per_particle))
            p.create_dataset('dens', data=rho)
            p.create_dataset('pres', data=pres)
            p.create_dataset('ene', data=ene)
            p.create_dataset('sml', data=sml)
            p.create_dataset('sound', data=np.sqrt(1.6667 * pres / rho))

            # Other fields
            p.create_dataset('id', data=np.arange(n, dtype=np.int32))
            p.create_dataset('is_ghost', data=np.zeros(n, dtype=np.int32))
            p.create_dataset('alpha', data=np.ones(n))
            p.create_dataset('balsara', data=np.ones(n))
            p.create_dataset('gradh', data=np.ones(n))
            p.create_dataset('neighbor', data=np.full(n, 50, dtype=np.int32))
            p.create_dataset('phi', data=np.zeros(n))

            # Energy
            e = f.create_group('energy')
            e.create_dataset('kinetic', data=[0.0])
            e.create_dataset('potential', data=[0.0])
            e.create_dataset('thermal', data=[ene.sum() * mass_per_particle])
            e.create_dataset('total', data=[ene.sum() * mass_per_particle])

    print(f"Wrote: {filename}")
    print(f"  Particles: {n:,}")
    print(f"  Mass/particle: {mass_per_particle:.6f}")
    print(f"  Total mass: {total_mass:.4f}")


def main():
    if len(sys.argv) < 3:
        print("Usage: python refine_relaxed_sphere.py <template.h5> <output.h5> [--neighbors N]")
        print("Example: python refine_relaxed_sphere.py snapshot_0128.h5 refined_1M.h5 --neighbors 8")
        return

    template_file = sys.argv[1]
    output_file = sys.argv[2]

    # Parse optional neighbors argument
    n_neighbors = 6  # Default: 6 neighbors → ~3.5x particles
    for i, arg in enumerate(sys.argv):
        if arg == '--neighbors' and i+1 < len(sys.argv):
            n_neighbors = int(sys.argv[i+1])

    print(f"Loading template: {template_file}")
    data = load_relaxed(template_file)

    total_mass = data['mass'] * data['n_particles']
    print(f"Original: {data['n_particles']:,} particles, total mass = {total_mass:.4f}")

    # Refine
    print(f"\nRefining with {n_neighbors} neighbors per particle...")
    new_pos, data = refine_particles(data, n_neighbors)

    # Interpolate quantities
    print("\nInterpolating SPH quantities...")
    rho, pres, ene, sml = interpolate_quantities(new_pos, data)

    # Write
    print(f"\nWriting output...")
    write_snapshot(output_file, new_pos, rho, pres, ene, sml, total_mass, template_file)

    print("\n=== NEXT STEPS ===")
    print("1. Run SHORT relaxation (~20-50 steps) to optimize local arrangement")
    print("2. Then run hydrostatic simulation")
    print("Expected: Much faster than random sampling (10-20 steps vs 100+)")


if __name__ == "__main__":
    main()
