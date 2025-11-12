#!/usr/bin/env python3
"""
Wrapper script to run SPH simulations with a detailed progress bar.
Monitors simulation output and displays real-time progress.
"""

import sys
import subprocess
import re
import os
from datetime import datetime

def format_time(seconds):
    """Format seconds into HH:MM:SS."""
    hours = int(seconds // 3600)
    minutes = int((seconds % 3600) // 60)
    secs = int(seconds % 60)
    return f"{hours:02d}:{minutes:02d}:{secs:02d}"

def format_scientific(value):
    """Format value in scientific notation."""
    if value < 0.01 or value > 1000:
        return f"{value:.3e}"
    else:
        return f"{value:.4f}"

def run_with_progress(command, label="Simulation"):
    """Run command and display progress bar based on output parsing."""
    
    print(f"{'='*70}")
    print(f"{label}")
    print(f"{'='*70}")
    
    # Start process
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1
    )
    
    params = {
        'particle_num': None,
        'end_time': None,
        'current_time': 0.0,
        'timestep': 0,
        'dt': 0.0,
        'snapshot_count': 0,
        'start_real_time': datetime.now(),
        'last_snapshot_time': 0.0
    }
    
    simulation_started = False
    last_progress_line = ""
    show_setup = True
    
    # Read output line by line
    for line in iter(process.stdout.readline, ''):
        line = line.rstrip()
        
        # Parse "loop:" lines from SPH output
        loop_match = re.match(r'loop:\s*(\d+),\s*time:\s*([\d.eE+-]+),\s*dt:\s*([\d.eE+-]+),\s*num:\s*(\d+)', line)
        
        if loop_match:
            # Extract loop information
            params['timestep'] = int(loop_match.group(1))
            params['current_time'] = float(loop_match.group(2))
            params['dt'] = float(loop_match.group(3))
            if params['particle_num'] is None:
                params['particle_num'] = int(loop_match.group(4))
            
            # Start progress display after first loop
            if not simulation_started:
                simulation_started = True
                print(f"\n{'─'*70}")
                print("Simulation Progress:")
                print(f"{'─'*70}")
            
            # Calculate progress
            if params['end_time'] and params['end_time'] > 0:
                progress_pct = min(100.0, (params['current_time'] / params['end_time']) * 100)
                elapsed = (datetime.now() - params['start_real_time']).total_seconds()
                
                # Build progress bar
                bar_width = 30
                filled = int(bar_width * progress_pct / 100)
                bar = '█' * filled + '░' * (bar_width - filled)
                
                # Estimate time remaining
                if progress_pct > 1.0:
                    total_time = elapsed / (progress_pct / 100)
                    remaining = total_time - elapsed
                    eta_str = format_time(remaining)
                    elapsed_str = format_time(elapsed)
                else:
                    eta_str = "--:--:--"
                    elapsed_str = format_time(elapsed)
                
                # Build detailed progress line
                progress_line = (
                    f"[{bar}] {progress_pct:5.1f}% │ "
                    f"Step: {params['timestep']:>6} │ "
                    f"t: {format_scientific(params['current_time'])}/{params['end_time']:.1f} │ "
                    f"dt: {format_scientific(params['dt'])} │ "
                    f"Snaps: {params['snapshot_count']:>2} │ "
                    f"Time: {elapsed_str} │ "
                    f"ETA: {eta_str}"
                )
                
                # Update progress bar (overwrite previous line)
                # Only update if enough change (avoid spam, update every 10 steps or 1% progress)
                if (params['timestep'] % 10 == 0 or 
                    int(progress_pct * 10) != int(getattr(run_with_progress, '_last_pct', 0) * 10)):
                    
                    # Clear line and print
                    print(f"\r{progress_line}", end='', flush=True)
                    run_with_progress._last_pct = progress_pct
                    last_progress_line = progress_line
        
        elif not simulation_started:
            # Show initial setup output
            if show_setup:
                # Parse important initial parameters
                if 'particle_num =' in line or 'Initial condition made' in line:
                    match = re.search(r'particle_num\s*=\s*(\d+)', line)
                    if match:
                        params['particle_num'] = int(match.group(1))
                    print(line)
                elif '* end time' in line.lower() or 'end time' in line.lower():
                    match = re.search(r'[*\s]*end time\s*[=:]\s*([\d.]+)', line, re.IGNORECASE)
                    if match:
                        params['end_time'] = float(match.group(1))
                    print(line)
                elif 'output time' in line.lower() and params['end_time'] is None:
                    # Sometimes end time is in output time section
                    print(line)
                elif any(keyword in line for keyword in ['SPH type:', 'Sample:', 'Lane-Emden:', 'parameters', 'time', 'Wrote snapshot']):
                    print(line)
                    if 'Wrote snapshot' in line:
                        params['snapshot_count'] += 1
        
        else:
            # During simulation, check for snapshots
            if 'Wrote snapshot' in line or 'snapshot_' in line:
                params['snapshot_count'] += 1
            
            # Show important messages on new lines
            if any(keyword in line.lower() for keyword in ['error', 'warning', 'failed', 'simulation complete']):
                print(f"\n{line}")
    
    # Wait for process to complete
    return_code = process.wait()
    
    # Final status
    print()  # New line after progress bar
    print(f"{'─'*70}")
    
    if return_code == 0:
        elapsed = (datetime.now() - params['start_real_time']).total_seconds()
        print(f"✓ Simulation completed successfully")
        print(f"  Total time:    {format_time(elapsed)}")
        print(f"  Particles:     {params['particle_num']:,}")
        print(f"  Timesteps:     {params['timestep']:,}")
        print(f"  Snapshots:     {params['snapshot_count']}")
        print(f"  Final time:    {params['current_time']:.6f} / {params['end_time']:.1f}")
        if params['timestep'] > 0 and elapsed > 0:
            steps_per_sec = params['timestep'] / elapsed
            print(f"  Performance:   {steps_per_sec:.1f} steps/sec")
    else:
        print(f"✗ Simulation failed with exit code {return_code}")
    
    print(f"{'='*70}\n")
    
    return return_code

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: run_with_progress.py <label> <command> [args...]")
        print("Example: run_with_progress.py 'SSPH Test' ./build/sph lane_emden")
        sys.exit(1)
    
    label = sys.argv[1]
    command = sys.argv[2:]
    
    exit_code = run_with_progress(command, label)
    sys.exit(exit_code)
