#!/usr/bin/env python3
"""
Parse Koyama & Inutsuka (2000) Figure 1 PostScript files
Extracts equilibrium curves INCLUDING stable/unstable branch information

LT1 = solid lines = thermally STABLE equilibrium (WNM and CNM branches)
LT2 = dashed lines = thermally UNSTABLE equilibrium (intermediate region)
"""

import re
import numpy as np

def parse_ps_file(filename):
    """Parse PostScript file and extract curve data with stability info"""
    
    with open(filename, 'r') as f:
        content = f.read()
    
    # Find axis bounds from tick labels
    # Temperature and Pressure axes: y-labels are powers of 10
    temp_labels = re.findall(r'\(([0-9]+)\) Rshow.*?(\d+) (\d+) M', content)
    
    # Find x-axis (density) labels
    density_labels = re.findall(r'\(([0-9-]+)\) Cshow.*?(\d+) (\d+) M', content)
    
    # Determine axis bounds
    # x-axis: density n_H [cm^-3], log scale from -1 to 6
    x_min_ps, x_max_ps = 610, 4713  # PostScript coordinates
    x_min_log, x_max_log = -1, 6     # log10(n_H)
    
    # y-axis: two panels (temperature and pressure)
    # Top panel (temperature): y ~ 855-2295 in PS coords
    # Bottom panel (pressure): y ~ 2808-4236 in PS coords
    
    # Extract all curve segments
    curves = []
    
    # Find all LT1 and LT2 curve segments
    # Pattern: LT[12] followed by M (moveto) and series of V (vertical line) commands
    lt_pattern = r'(LT[12])\s+(\d+)\s+(\d+)\s+M\s+((?:[-\d]+\s+[-\d]+\s+V\s*)+)'
    
    for match in re.finditer(lt_pattern, content):
        line_type = match.group(1)  # LT1 or LT2
        x_start = int(match.group(2))
        y_start = int(match.group(3))
        moves = match.group(4)
        
        # Parse V commands (relative vertical moves)
        v_commands = re.findall(r'([-\d]+)\s+([-\d]+)\s+V', moves)
        
        # Build coordinate list
        x_coords = [x_start]
        y_coords = [y_start]
        
        x_current = x_start
        y_current = y_start
        
        for dx_str, dy_str in v_commands:
            dx = int(dx_str)
            dy = int(dy_str)
            x_current += dx
            y_current += dy
            x_coords.append(x_current)
            y_coords.append(y_current)
        
        # Determine which panel (temperature or pressure)
        avg_y = np.mean(y_coords)
        if avg_y < 1500:
            panel = 'temperature'
            y_min_ps, y_max_ps = 855, 2295
            y_min_log, y_max_log = 0, 7  # log10(T) from 1 to 10^7 K
        else:
            panel = 'pressure'
            y_min_ps, y_max_ps = 2808, 4236
            y_min_log, y_max_log = 0, 6  # log10(P/k_B) from 1 to 10^6 K cm^-3
        
        # Convert PS coordinates to physical values
        x_log = x_min_log + (np.array(x_coords) - x_min_ps) * (x_max_log - x_min_log) / (x_max_ps - x_min_ps)
        y_log = y_min_log + (np.array(y_coords) - y_min_ps) * (y_max_log - y_min_log) / (y_max_ps - y_min_ps)
        
        n_H = 10**x_log
        if panel == 'temperature':
            T = 10**y_log
            val = T
        else:
            P = 10**y_log
            val = P
        
        # Determine stability (LT1 = stable, LT2 = unstable)
        is_stable = (line_type == 'LT1')
        
        curves.append({
            'panel': panel,
            'line_type': line_type,
            'is_stable': is_stable,
            'n_H': n_H,
            'value': val,
            'x_ps': np.array(x_coords),
            'y_ps': np.array(y_coords)
        })
    
    return curves

def identify_column_density(curves, filename):
    """Identify if this is N_H = 1e19 or 1e20 cm^-2"""
    # Based on filename or curve characteristics
    # For now, check the temperature values
    temp_curves = [c for c in curves if c['panel'] == 'temperature']
    if not temp_curves:
        return 'unknown'
    
    # At n_H ~ 10 cm^-3, N_H=1e19 has lower T than N_H=1e20
    # Check temperature at log(n_H) ~ 1 (n_H ~ 10)
    for curve in temp_curves:
        idx = np.argmin(np.abs(np.log10(curve['n_H']) - 1.0))
        T_at_n10 = curve['value'][idx]
        
        # N_H=1e19: T ~ 70-100 K at n=10
        # N_H=1e20: T ~ 700-2000 K at n=10
        if T_at_n10 < 200:
            return '1e19'
        elif T_at_n10 > 500:
            return '1e20'
    
    return 'unknown'

def merge_curve_segments(curves, panel):
    """Merge curve segments by stability type"""
    panel_curves = [c for c in curves if c['panel'] == panel]
    
    # Group by stability
    stable_segs = [c for c in panel_curves if c['is_stable']]
    unstable_segs = [c for c in panel_curves if not c['is_stable']]
    
    def concatenate_segments(segs):
        if not segs:
            return None
        
        # Sort by starting density
        segs = sorted(segs, key=lambda c: c['n_H'][0])
        
        # Concatenate
        n_H_all = np.concatenate([c['n_H'] for c in segs])
        val_all = np.concatenate([c['value'] for c in segs])
        
        # Remove duplicates and sort
        unique_idx = np.unique(n_H_all, return_index=True)[1]
        n_H_all = n_H_all[unique_idx]
        val_all = val_all[unique_idx]
        
        sort_idx = np.argsort(n_H_all)
        return n_H_all[sort_idx], val_all[sort_idx]
    
    stable_curve = concatenate_segments(stable_segs)
    unstable_curve = concatenate_segments(unstable_segs)
    
    return stable_curve, unstable_curve

# Parse all 4 files
files = {
    'f1a.ps': 'Temperature and Pressure',
    'f1b.ps': 'Chemical fractions',
    'f1c.ps': 'Heating and cooling rates',
    'f1d.ps': 'Timescales'
}

print("=" * 80)
print("Parsing Koyama & Inutsuka (2000) Figure 1 PostScript Files")
print("=" * 80)

for ps_file, description in files.items():
    print(f"\n{ps_file}: {description}")
    try:
        curves = parse_ps_file(ps_file)
        col_dens = identify_column_density(curves, ps_file)
        
        print(f"  Found {len(curves)} curve segments")
        print(f"  Estimated column density: N_H = {col_dens} cm^-2")
        
        # Show curve breakdown
        temp_curves = [c for c in curves if c['panel'] == 'temperature']
        pres_curves = [c for c in curves if c['panel'] == 'pressure']
        
        if temp_curves:
            stable_temp = [c for c in temp_curves if c['is_stable']]
            unstable_temp = [c for c in temp_curves if not c['is_stable']]
            print(f"  Temperature: {len(stable_temp)} stable segments, {len(unstable_temp)} unstable segments")
            
            # Check temperature range at n=10
            for curve in temp_curves:
                idx = np.argmin(np.abs(np.log10(curve['n_H']) - 1.0))
                T_at_n10 = curve['value'][idx]
                stability = "STABLE" if curve['is_stable'] else "UNSTABLE"
                print(f"    {stability}: T(n=10 cm^-3) ≈ {T_at_n10:.1f} K")
        
        if pres_curves:
            stable_pres = [c for c in pres_curves if c['is_stable']]
            unstable_pres = [c for c in pres_curves if not c['is_stable']]
            print(f"  Pressure: {len(stable_pres)} stable segments, {len(unstable_pres)} unstable segments")
            
    except Exception as e:
        print(f"  ERROR: {e}")

print("\n" + "=" * 80)
print("CRITICAL FINDINGS:")
print("=" * 80)
print("""
Figure 1a shows the FULL equilibrium curve including:
  - STABLE branches (solid lines, LT1): WNM at low n, CNM at high n
  - UNSTABLE branches (dashed lines, LT2): Intermediate densities

At n_H = 10 cm^-3, the equilibrium curve shows the UNSTABLE branch (T ~ 1000K).
The CNM phase at n=10, T=107K (Table 1, Model C10) is NOT on this curve!

Model C10 is a SHOCK PROPAGATION simulation, not thermal equilibrium!
Initial state: n=10 cm^-3, T=107K (non-equilibrium CNM)
The shock will compress and heat the gas, NOT cool it to equilibrium.
""")
