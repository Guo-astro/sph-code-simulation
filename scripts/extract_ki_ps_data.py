#!/usr/bin/env python3
"""
Extract curve data from K&I 2000 PostScript figure files.

This script parses gnuplot-generated PostScript files to extract
the numerical data points for each curve with pixel-perfect accuracy.

Usage:
    python extract_ki_ps_data.py

PostScript coordinate system:
- M (moveto): absolute coordinates
- L (lineto): absolute coordinates  
- V (rlineto): relative coordinates (dx, dy from current position)
- R (rmoveto): relative move

Figure Contents (from K&I 2000 ApJ 532, 980):
- f1a.ps: Temperature T(n) and Pressure P/k(n) vs density
- f1b.ps: Electron fraction x_e(n) and H2 fraction x_H2(n) vs density
- f1c.ps: Heating and cooling rates Gamma(n), Lambda(n) vs density
- f1d.ps: Thermal equilibrium timescale t_eq(n) vs density
"""

import re
import numpy as np
from pathlib import Path


def parse_ps_curve(ps_content: str) -> list[dict]:
    """
    Parse gnuplot PostScript content and extract curve data.
    
    Returns list of curves, each containing:
    - style: LT1, LT2, etc.
    - points: list of (x, y) in PS coordinates
    """
    curves = []
    current_curve = None
    current_x, current_y = 0, 0
    
    lines = ps_content.split('\n')
    
    for line in lines:
        line = line.strip()
        
        # Detect line style change (start of new curve)
        lt_match = re.match(r'^LT(\d+)$', line)
        if lt_match:
            if current_curve and len(current_curve['points']) > 0:
                curves.append(current_curve)
            current_curve = {'style': line, 'points': []}
            continue
        
        # Parse moveto command (absolute)
        match = re.match(r'^(\d+)\s+(\d+)\s+M$', line)
        if match:
            current_x = int(match.group(1))
            current_y = int(match.group(2))
            if current_curve is not None:
                current_curve['points'].append((current_x, current_y))
            continue
        
        # Parse relative lineto (V = rlineto)
        match = re.match(r'^(-?\d+)\s+(-?\d+)\s+V$', line)
        if match:
            dx = int(match.group(1))
            dy = int(match.group(2))
            current_x += dx
            current_y += dy
            if current_curve is not None:
                current_curve['points'].append((current_x, current_y))
            continue
        
        # Parse absolute lineto (L)
        match = re.match(r'^(\d+)\s+(\d+)\s+L$', line)
        if match:
            current_x = int(match.group(1))
            current_y = int(match.group(2))
            if current_curve is not None:
                current_curve['points'].append((current_x, current_y))
            continue
    
    # Don't forget the last curve
    if current_curve and len(current_curve['points']) > 0:
        curves.append(current_curve)
    
    return curves


def convert_coordinates(ps_coords: list[tuple], 
                        x_ps_range: tuple, y_ps_range: tuple,
                        x_log_range: tuple, y_log_range: tuple,
                        log_x: bool = True, log_y: bool = True) -> list[tuple]:
    """
    Convert PS coordinates to physical values.
    
    Args:
        ps_coords: List of (x, y) PS coordinates
        x_ps_range: (x_min_ps, x_max_ps)
        y_ps_range: (y_min_ps, y_max_ps)
        x_log_range: (log_x_min, log_x_max) physical range
        y_log_range: (log_y_min, log_y_max) physical range
        log_x: If True, x is log scale
        log_y: If True, y is log scale
    """
    x_min_ps, x_max_ps = x_ps_range
    y_min_ps, y_max_ps = y_ps_range
    log_x_min, log_x_max = x_log_range
    log_y_min, log_y_max = y_log_range
    
    physical_points = []
    for px, py in ps_coords:
        # Linear interpolation in PS space
        frac_x = (px - x_min_ps) / (x_max_ps - x_min_ps)
        frac_y = (py - y_min_ps) / (y_max_ps - y_min_ps)
        
        # Convert to physical coordinates
        if log_x:
            log_x_val = log_x_min + frac_x * (log_x_max - log_x_min)
            x_val = 10**log_x_val
        else:
            x_val = log_x_min + frac_x * (log_x_max - log_x_min)
            
        if log_y:
            log_y_val = log_y_min + frac_y * (log_y_max - log_y_min)
            y_val = 10**log_y_val
        else:
            y_val = log_y_min + frac_y * (log_y_max - log_y_min)
        
        physical_points.append((x_val, y_val))
    
    return physical_points


def extract_f1a_data(ps_file: Path) -> dict:
    """
    Extract temperature and pressure curves from f1a.ps.
    
    f1a.ps axis system (from tick marks in PS file):
    - X axis: log(n [cm^-3]) from -2 to 4
    - Y axis: log(T [K]) or log(P/k [cm^-3 K]) from 0 to 7
    
    PS coordinates (from axis tick positions):
    - X: 1320 = log(n) = -2, 4713 = log(n) = 4
    - Y: 551 = 0, 2913 = 7
    
    Curves (identified by starting y-coordinate):
    - Temperature curves start at y ~ 800-900 (high n, low T)
    - Pressure curves start at y ~ 2700-2800 (high n, high P)
    
    LT1 = solid line = N_H = 10^19 cm^-2
    LT2 = dashed line = N_H = 10^20 cm^-2
    """
    content = ps_file.read_text()
    curves = parse_ps_curve(content)
    
    print(f"\nFound {len(curves)} curves in {ps_file.name}:")
    
    # PS coordinate ranges
    x_ps_range = (1320, 4713)
    y_ps_range = (551, 2913)
    x_log_range = (-2.0, 4.0)  # log(n)
    y_log_range = (0.0, 7.0)   # log(T) or log(P/k)
    
    result = {}
    
    for i, curve in enumerate(curves):
        if not curve['points']:
            continue
            
        y_start = curve['points'][0][1]
        style = curve['style']
        
        # Convert to physical coordinates
        physical = convert_coordinates(
            curve['points'], x_ps_range, y_ps_range,
            x_log_range, y_log_range
        )
        
        # Identify curve type by starting y position
        if y_start < 1500:  # Temperature curves (start at low y = low log T)
            curve_type = 'temperature'
            if style == 'LT1':
                result['temperature_N19'] = physical
            elif style == 'LT2':
                result['temperature_N20'] = physical
        else:  # Pressure curves (start at high y = high log P/k)
            curve_type = 'pressure'
            if style == 'LT1':
                result['pressure_N19'] = physical
            elif style == 'LT2':
                result['pressure_N20'] = physical
        
        # Report
        n_start = physical[0][0]
        n_end = physical[-1][0]
        y_start_phys = physical[0][1]
        y_end_phys = physical[-1][1]
        print(f"  Curve {i}: {style} -> {curve_type}")
        print(f"    n: {n_start:.2e} -> {n_end:.2e}")
        print(f"    {curve_type}: {y_start_phys:.1f} -> {y_end_phys:.1f}")
    
    return result


def extract_f1b_data(ps_file: Path) -> dict:
    """
    Extract electron fraction, H2 fraction, and CO fraction curves from f1b.ps.

    f1b.ps axis system (from tick marks in PS file):
    - X axis: log(n [cm^-3]) from -2 to 6
    - Y axis: log(x_i) from -8 to 0

    The 6 curves in order are:
    - Curve 0 (LT1): x_e for N_H = 10^19 cm^-2
    - Curve 1 (LT2): x_e for N_H = 10^20 cm^-2
    - Curve 2 (LT1): x_H2 for N_H = 10^19 cm^-2
    - Curve 3 (LT2): x_H2 for N_H = 10^20 cm^-2
    - Curve 4 (LT1): x_CO for N_H = 10^19 cm^-2
    - Curve 5 (LT2): x_CO for N_H = 10^20 cm^-2

    Note: Curves are drawn from high n to low n (right to left).

    LT1 = solid = N_H = 10^19 cm^-2
    LT2 = dashed = N_H = 10^20 cm^-2
    """
    content = ps_file.read_text()
    curves = parse_ps_curve(content)

    print(f"\nFound {len(curves)} curves in {ps_file.name}:")

    # PS coordinate ranges (from tick positions)
    x_ps_range = (1320, 4713)
    y_ps_range = (551, 2913)
    x_log_range = (-2.0, 6.0)  # log(n)
    y_log_range = (-8.0, 0.0)  # log(x_i)

    result = {}

    # Map curve index to species based on K&I Fig 1b structure
    # The curves appear in order: x_e (2 curves), x_H2 (2 curves), x_CO (2 curves)
    curve_map = {
        0: ('x_e', 'N19'),
        1: ('x_e', 'N20'),
        2: ('x_H2', 'N19'),
        3: ('x_H2', 'N20'),
        4: ('x_CO', 'N19'),
        5: ('x_CO', 'N20'),
    }

    for i, curve in enumerate(curves):
        if not curve['points'] or i not in curve_map:
            continue

        species, nh_label = curve_map[i]
        style = curve['style']

        # Convert to physical coordinates
        physical = convert_coordinates(
            curve['points'], x_ps_range, y_ps_range,
            x_log_range, y_log_range
        )

        # Store with proper key
        key = f'{species}_{nh_label}'
        result[key] = physical

        # Report
        n_start = physical[0][0]
        n_end = physical[-1][0]
        y_start_phys = physical[0][1]
        y_end_phys = physical[-1][1]
        print(f"  Curve {i}: {style} -> {species} ({nh_label})")
        print(f"    n: {n_start:.2e} -> {n_end:.2e}")
        print(f"    {species}: {y_start_phys:.2e} -> {y_end_phys:.2e}")

    return result


def extract_f1c_data(ps_file: Path) -> dict:
    """
    Extract heating and cooling rate curves from f1c.ps.
    
    NOTE: For accurate f1c extraction, use extract_f1c.py instead.
    This function only provides basic curve parsing.
    
    f1c.ps shows heating and cooling rates vs density:
    - X axis: log(n [cm^-3]) from -2 to 6
    - Y axis: log(rate [erg s^-1 H^-1]) from -28 to -25
    
    PS coordinate system (from tick marks):
    - X: 1320 = log(n)=-2, 4713 = log(n)=6 (424 PS units per log unit)
    - Y: 551 = log(rate)=-28, tick spacing = 675 PS units per log unit
    """
    content = ps_file.read_text()
    curves = parse_ps_curve(content)
    
    print(f"\nFound {len(curves)} curves in {ps_file.name}:")
    for i, curve in enumerate(curves):
        print(f"  Curve {i}: {curve['style']}, {len(curve['points'])} points")
        if curve['points']:
            p0, pn = curve['points'][0], curve['points'][-1]
            print(f"    PS range: ({p0[0]},{p0[1]}) -> ({pn[0]},{pn[1]})")
    
    return {'curves': curves}


def extract_f1d_data(ps_file: Path) -> dict:
    """
    Extract timescale curves from f1d.ps.
    
    f1d.ps shows:
    - Thermal equilibrium timescale t_eq(n)
    - Other characteristic timescales
    
    X axis: log(n [cm^-3]) from -2 to 4
    Y axis: log(t [yr])
    """
    content = ps_file.read_text()
    curves = parse_ps_curve(content)
    
    print(f"\nFound {len(curves)} curves in {ps_file.name}:")
    
    # Try to identify axis ranges from the PS file
    # Look for typical tick values
    
    # For timescales, y-axis is typically log(t) in years
    # t_eq ~ 10^4 to 10^8 years depending on density
    
    for i, curve in enumerate(curves):
        print(f"  Curve {i}: {curve['style']}, {len(curve['points'])} points")
        if curve['points']:
            p0, pn = curve['points'][0], curve['points'][-1]
            print(f"    PS range: ({p0[0]},{p0[1]}) -> ({pn[0]},{pn[1]})")
    
    return {'curves': curves}


def save_data(result: dict, output_dir: Path, prefix: str):
    """Save extracted data to text files."""
    for name, points in result.items():
        if isinstance(points, list) and len(points) > 0:
            arr = np.array(points)
            filename = f'{prefix}_{name}.txt'
            np.savetxt(output_dir / filename, arr, 
                       header=f'{name}: col1=n [cm^-3], col2=value', fmt='%.6e')
            print(f"  Saved {filename}: {len(arr)} points")


def print_sample_data(result: dict, name: str, ylabel: str):
    """Print sample data points for verification."""
    if name not in result:
        return
        
    points = result[name]
    print(f"\n{name}: {len(points)} points")
    print(f"{'n [cm^-3]':>12s}  {ylabel:>12s}")
    print("-" * 30)
    
    # Print at specific densities
    densities = [1e4, 1e3, 1e2, 1e1, 1e0, 1e-1, 1e-2]
    arr = np.array(points)
    for target_n in densities:
        idx = np.argmin(np.abs(arr[:, 0] - target_n))
        n, val = arr[idx]
        if abs(np.log10(n) - np.log10(target_n)) < 0.5:
            print(f"{n:12.2e}  {val:12.4e}")


def main():
    """Main extraction routine."""
    ps_dir = Path("/Users/guo/Downloads/sphcode/docs/papers/cooling-heating")
    output_dir = Path("/Users/guo/Downloads/sphcode/data/ki2000_extracted")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("="*70)
    print("K&I 2000 PostScript Figure Data Extraction")
    print("="*70)
    
    # Extract f1a.ps (Temperature and Pressure)
    f1a = ps_dir / "f1a.ps"
    if f1a.exists():
        print("\n" + "="*70)
        print("f1a.ps: Temperature T(n) and Pressure P/k(n)")
        print("="*70)
        data_f1a = extract_f1a_data(f1a)
        
        # Print sample data
        print_sample_data(data_f1a, 'temperature_N19', 'T [K]')
        print_sample_data(data_f1a, 'temperature_N20', 'T [K]')
        print_sample_data(data_f1a, 'pressure_N19', 'P/k')
        print_sample_data(data_f1a, 'pressure_N20', 'P/k')
        
        # Save
        print("\nSaving f1a data:")
        save_data(data_f1a, output_dir, 'f1a')
    
    # Extract f1b.ps (Electron and H2 fractions)
    f1b = ps_dir / "f1b.ps"
    if f1b.exists():
        print("\n" + "="*70)
        print("f1b.ps: Electron fraction x_e(n) and H2 fraction x_H2(n)")
        print("="*70)
        data_f1b = extract_f1b_data(f1b)
        
        print_sample_data(data_f1b, 'x_e_N19', 'x_e')
        print_sample_data(data_f1b, 'x_e_N20', 'x_e')
        print_sample_data(data_f1b, 'x_H2_N19', 'x_H2')
        print_sample_data(data_f1b, 'x_H2_N20', 'x_H2')
        
        print("\nSaving f1b data:")
        save_data(data_f1b, output_dir, 'f1b')
    
    # Extract f1c.ps (Heating/Cooling rates)
    f1c = ps_dir / "f1c.ps"
    if f1c.exists():
        print("\n" + "="*70)
        print("f1c.ps: Heating and Cooling rates")
        print("="*70)
        data_f1c = extract_f1c_data(f1c)
    
    # Extract f1d.ps (Timescales)
    f1d = ps_dir / "f1d.ps"
    if f1d.exists():
        print("\n" + "="*70)
        print("f1d.ps: Thermal timescales")
        print("="*70)
        data_f1d = extract_f1d_data(f1d)
    
    print("\n" + "="*70)
    print("Extraction complete!")
    print(f"Output directory: {output_dir}")
    print("="*70)


if __name__ == '__main__':
    main()
