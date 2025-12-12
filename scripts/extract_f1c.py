#!/usr/bin/env python3
"""
Extract heating/cooling curves from K&I 2000 f1c.ps.

Panel (c) shows heating and cooling rates vs density n.
- X axis: log(n [cm^-3]) from 0 to 6
- Y axis: log(Gamma, Lambda [erg s^-1 H^-1]) from -28 to -25

LT1 = solid lines (N_H = 10^19 cm^-2)
LT2 = dashed lines (N_H = 10^20 cm^-2)

From the PS file, after skipping grid lines (first 26 M commands),
there are multiple curves for different processes.
"""

import numpy as np
from pathlib import Path
import re


def parse_f1c_ps(ps_file: Path) -> list[dict]:
    """Parse f1c.ps and extract all data curves."""
    content = ps_file.read_text()
    lines = content.split('\n')

    curves = []
    current_curve = None
    current_style = None
    x, y = 0, 0
    m_count = 0  # Count moveto commands

    for line in lines:
        line = line.strip()

        # Track line style
        if re.match(r'^LT[0-9]$', line):
            current_style = line
            continue

        # Moveto - start of curve segment
        match = re.match(r'^(\d+)\s+(\d+)\s+M$', line)
        if match:
            m_count += 1
            # Skip first 26 M commands (grid lines, labels, etc.)
            if m_count <= 26:
                continue

            # Save previous curve if exists
            if current_curve and len(current_curve['points']) > 1:
                curves.append(current_curve)

            x, y = int(match.group(1)), int(match.group(2))
            current_curve = {
                'style': current_style,
                'points': [(x, y)]
            }
            continue

        # Relative lineto
        match = re.match(r'^(-?\d+)\s+(-?\d+)\s+V$', line)
        if match and current_curve:
            dx, dy = int(match.group(1)), int(match.group(2))
            x += dx
            y += dy
            current_curve['points'].append((x, y))
            continue

        # Absolute lineto
        match = re.match(r'^(\d+)\s+(\d+)\s+L$', line)
        if match and current_curve:
            x, y = int(match.group(1)), int(match.group(2))
            current_curve['points'].append((x, y))
            continue

        # Relative moveto (jump within curve)
        match = re.match(r'^(-?\d+)\s+(-?\d+)\s+R$', line)
        if match and current_curve:
            dx, dy = int(match.group(1)), int(match.group(2))
            x += dx
            y += dy
            # This is a discontinuity - save current and start new
            if len(current_curve['points']) > 1:
                curves.append(current_curve)
            current_curve = {
                'style': current_style,
                'points': [(x, y)]
            }
            continue

    # Don't forget last curve
    if current_curve and len(current_curve['points']) > 1:
        curves.append(current_curve)

    return curves


def convert_to_physical(curves: list[dict]) -> list[dict]:
    """Convert PS coordinates to physical units.

    PS coordinate system (from tick marks in f1c.ps):
    - X-axis ticks: 1320=-2, 1744=-1, 2168=0, 2592=1, 3017=2, 3441=3, 3865=4, 4289=5, 4713=6
    - Y-axis ticks: 551=-28, 1226=-27, 1901=-26, 2576=-25
    - Tick spacing: X = 424 PS units per log unit, Y = 675 PS units per log unit
    
    Note: Some data points extend above Y=2576 (above the -25 tick), so we use
    the tick spacing to extrapolate rather than clamping to the tick range.
    """
    # X-axis: linear mapping from tick marks
    x_ps_min, x_ps_max = 1320, 4713
    log_n_min, log_n_max = -2.0, 6.0
    
    # Y-axis: use tick spacing for accurate conversion
    # Y=551 corresponds to log(rate)=-28
    # Tick spacing: 675 PS units per log unit
    y_ps_ref = 551  # Reference point
    log_rate_ref = -28.0  # log(rate) at reference
    y_ps_per_log = 675.0  # PS units per log unit (from tick spacing)

    physical_curves = []

    for curve in curves:
        phys_points = []
        for px, py in curve['points']:
            # X: Linear interpolation
            frac_x = (px - x_ps_min) / (x_ps_max - x_ps_min)
            log_n = log_n_min + frac_x * (log_n_max - log_n_min)
            
            # Y: Use tick spacing for accurate conversion
            log_rate = log_rate_ref + (py - y_ps_ref) / y_ps_per_log

            n = 10**log_n
            rate = 10**log_rate

            phys_points.append((n, rate))

        physical_curves.append({
            'style': curve['style'],
            'points': phys_points
        })

    return physical_curves


def identify_curves(curves: list[dict]) -> dict:
    """Identify which curve corresponds to which process.

    Based on K&I 2000 Figure 1c labels and curve shapes:
    LT1 (solid, N19):
      - Curve starting high right, going down-left: Total heating/cooling
      - Curve starting at n~1000: Lya cooling (high T region)
      - Curve going up then down: CII cooling
      - Curve at bottom: grain/CR heating

    LT2 (dashed, N20):
      - Similar structure but for higher column density
    """
    result = {
        'N19': {'total_heat': None, 'total_cool': None, 'curves': []},
        'N20': {'total_heat': None, 'total_cool': None, 'curves': []},
    }

    for i, curve in enumerate(curves):
        case = 'N19' if curve['style'] == 'LT1' else 'N20'
        result[case]['curves'].append({
            'index': i,
            'n_points': len(curve['points']),
            'n_range': (curve['points'][0][0], curve['points'][-1][0]),
            'rate_range': (curve['points'][0][1], curve['points'][-1][1]),
            'points': curve['points']
        })

    return result


def main():
    ps_file = Path('/Users/guo/Downloads/sphcode/docs/papers/cooling-heating/f1c.ps')
    output_dir = Path('/Users/guo/Downloads/sphcode/data/ki2000_extracted')

    print("Extracting curves from f1c.ps...")
    curves = parse_f1c_ps(ps_file)
    print(f"Found {len(curves)} data curves")

    # Convert to physical units
    physical_curves = convert_to_physical(curves)

    # Identify curves
    identified = identify_curves(physical_curves)

    # Print summary
    print("\nLT1 (N_H = 10^19 cm^-2) curves:")
    for c in identified['N19']['curves']:
        print(f"  Curve {c['index']}: {c['n_points']} pts, "
              f"n=[{c['n_range'][0]:.1e}, {c['n_range'][1]:.1e}], "
              f"rate=[{c['rate_range'][0]:.1e}, {c['rate_range'][1]:.1e}]")

    print("\nLT2 (N_H = 10^20 cm^-2) curves:")
    for c in identified['N20']['curves']:
        print(f"  Curve {c['index']}: {c['n_points']} pts, "
              f"n=[{c['n_range'][0]:.1e}, {c['n_range'][1]:.1e}], "
              f"rate=[{c['rate_range'][0]:.1e}, {c['rate_range'][1]:.1e}]")

    # Save all curves
    for case in ['N19', 'N20']:
        for i, c in enumerate(identified[case]['curves']):
            arr = np.array(c['points'])
            filename = f'f1c_{case}_curve{i}.txt'
            header = f'n [cm^-3], rate [erg/s/H] - {case} curve {i}'
            np.savetxt(output_dir / filename, arr, header=header, fmt='%.6e')
            print(f"Saved: {filename}")

    return identified


if __name__ == '__main__':
    main()
