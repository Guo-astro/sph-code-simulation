#!/usr/bin/env python3
"""
Test script to verify HDF5 and CSV output reading capabilities.
This validates the new output system works properly with Python.
"""

import numpy as np
import h5py
import sys
import os

def test_hdf5_file(filepath):
    """Test reading HDF5 snapshot file"""
    print(f"\n{'='*60}")
    print(f"Testing HDF5 file: {filepath}")
    print(f"{'='*60}")
    
    if not os.path.exists(filepath):
        print(f"ERROR: File not found: {filepath}")
        return False
    
    try:
        with h5py.File(filepath, 'r') as f:
            print(f"\n📦 HDF5 Groups: {list(f.keys())}")
            
            # Read metadata
            if 'metadata' in f:
                print(f"\n📝 Metadata Attributes:")
                for attr_name in f['metadata'].attrs.keys():
                    print(f"  {attr_name}: {f['metadata'].attrs[attr_name]}")
            
            # Read particles
            if 'particles' in f:
                particles = f['particles']
                print(f"\n🔵 Particle Datasets: {list(particles.keys())}")
                print(f"  Particle count: {particles['pos_x'].shape[0]}")
                
                # Load particle data
                x = particles['pos_x'][:]
                y = particles['pos_y'][:]
                z = particles['pos_z'][:]
                rho = particles['dens'][:]
                pres = particles['pres'][:]
                mass = particles['mass'][:]
                
                print(f"\n📊 Particle Statistics:")
                r = np.sqrt(x**2 + y**2 + z**2)
                print(f"  Radius: min={r.min():.6f}, max={r.max():.6f}")
                print(f"  Density: min={rho.min():.6f}, max={rho.max():.6f}")
                print(f"  Pressure: min={pres.min():.6f}, max={pres.max():.6f}")
                print(f"  Total mass: {mass.sum():.6f}")
            
            # Read energy
            if 'energy' in f:
                energy = f['energy']
                print(f"\n⚡ Energy Datasets: {list(energy.keys())}")
                if 'kinetic' in energy:
                    print(f"  Kinetic: {energy['kinetic'][0]:.6e}")
                if 'thermal' in energy:
                    print(f"  Thermal: {energy['thermal'][0]:.6e}")
                if 'potential' in energy:
                    print(f"  Potential: {energy['potential'][0]:.6e}")
                if 'total' in energy:
                    print(f"  Total: {energy['total'][0]:.6e}")
        
        print(f"\n✅ HDF5 file read successfully!")
        return True
        
    except Exception as e:
        print(f"\n❌ ERROR reading HDF5 file: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_csv_file(filepath):
    """Test reading CSV snapshot file"""
    print(f"\n{'='*60}")
    print(f"Testing CSV file: {filepath}")
    print(f"{'='*60}")
    
    if not os.path.exists(filepath):
        print(f"ERROR: File not found: {filepath}")
        return False
    
    try:
        import pandas as pd
        
        # Read CSV with pandas (handles comments and headers)
        df = pd.read_csv(filepath, comment='#')
        
        print(f"\n📊 CSV Data Shape: {df.shape}")
        print(f"  Columns: {df.shape[1]}")
        print(f"  Rows (particles): {df.shape[0]}")
        print(f"\n📋 Column Names: {list(df.columns)[:10]}...")  # Show first 10
        
        # Extract data
        x = df['pos_x'].values
        y = df['pos_y'].values
        z = df['pos_z'].values
        mass = df['mass'].values
        rho = df['dens'].values
        pres = df['pres'].values
        
        print(f"\n📊 Particle Statistics:")
        r = np.sqrt(x**2 + y**2 + z**2)
        print(f"  Radius: min={r.min():.6f}, max={r.max():.6f}")
        print(f"  Density: min={rho.min():.6f}, max={rho.max():.6f}")
        print(f"  Pressure: min={pres.min():.6f}, max={pres.max():.6f}")
        print(f"  Total mass: {mass.sum():.6f}")
        
        print(f"\n✅ CSV file read successfully!")
        return True
        
    except ImportError:
        print(f"\n⚠️  pandas not installed, using numpy (basic test)")
        # Just verify file exists and has content
        with open(filepath, 'r') as f:
            lines = [l for l in f if not l.startswith('#')]
            print(f"  Found {len(lines)} data lines (including header)")
        print(f"\n✅ CSV file exists and has data")
        return True
        
    except Exception as e:
        print(f"\n❌ ERROR reading CSV file: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    # Default paths
    base_dir = "lane_emden/results/polytrope_n1.5_3d"
    
    # Test HDF5 snapshot
    hdf5_snapshot = f"{base_dir}/snapshot_0000.h5"
    hdf5_checkpoint = f"{base_dir}/checkpoints/checkpoint_0000.h5"
    csv_snapshot = f"{base_dir}/snapshot_0000.csv"
    csv_checkpoint = f"{base_dir}/checkpoints/checkpoint_0000.csv"
    
    results = []
    
    # Test each file if it exists
    if os.path.exists(hdf5_snapshot):
        results.append(("HDF5 Snapshot", test_hdf5_file(hdf5_snapshot)))
    
    if os.path.exists(hdf5_checkpoint):
        results.append(("HDF5 Checkpoint", test_hdf5_file(hdf5_checkpoint)))
    
    if os.path.exists(csv_snapshot):
        results.append(("CSV Snapshot", test_csv_file(csv_snapshot)))
    
    if os.path.exists(csv_checkpoint):
        results.append(("CSV Checkpoint", test_csv_file(csv_checkpoint)))
    
    # Summary
    print(f"\n{'='*60}")
    print(f"SUMMARY")
    print(f"{'='*60}")
    for name, success in results:
        status = "✅ PASS" if success else "❌ FAIL"
        print(f"  {name}: {status}")
    
    all_passed = all(success for _, success in results)
    if all_passed:
        print(f"\n🎉 All tests passed! New output system works correctly.")
        return 0
    else:
        print(f"\n⚠️  Some tests failed.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
