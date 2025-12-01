#!/usr/bin/env python3
"""
Digitize Koyama & Inutsuka (2000) Figure 1 from PostScript files.

Extracts pixel-perfect curve data from the original PostScript figures
and converts to physical units (density, temperature, pressure, etc.).

Based on the paper: Koyama & Inutsuka (2000) ApJ 532:980-993
Figure 1 shows equilibrium curves for two column densities (1e19, 1e20 cm^-2).
"""

import re
import numpy as np
from pathlib import Path
from typing import List, Tuple, Dict

class PostScriptParser:
    """Parse PostScript files to extract curve data."""
    
    def __init__(self, ps_file: str):
        self.ps_file = Path(ps_file)
        with open(self.ps_file, 'r') as f:
            self.content = f.read()
        
        # Extract scale and translation from PostScript
        self.parse_transform()
    
    def parse_transform(self):
        """Extract coordinate transformation from PostScript."""
        # Look for "translate" and "scale" commands
        translate_match = re.search(r'(\d+)\s+(\d+)\s+translate', self.content)
        scale_match = re.search(r'([\d.]+)\s+([\d.]+)\s+scale', self.content)
        
        if translate_match:
            self.translate_x = float(translate_match.group(1))
            self.translate_y = float(translate_match.group(2))
        else:
            self.translate_x = 0
            self.translate_y = 0
        
        if scale_match:
            self.scale_x = float(scale_match.group(1))
            self.scale_y = float(scale_match.group(2))
        else:
            self.scale_x = 1.0
            self.scale_y = 1.0
        
        print(f"Transform: translate({self.translate_x}, {self.translate_y}), scale({self.scale_x}, {self.scale_y})")
    
    def extract_curves(self) -> Dict[str, List[Tuple[float, float]]]:
        """Extract all curves from PostScript file.
        
        Returns:
            Dictionary mapping curve names to list of (x, y) points
        """
        curves = {}
        curve_counter = {}  # Track multiple curves with same line type
        
        # Find all curve sections (marked by line type commands like LT0, LT1, etc.)
        # Pattern: LT<number> followed by M (moveto) and V (vertical line) commands
        curve_pattern = r'(LT\d+)\s*\n(\d+)\s+(\d+)\s+M\s*\n((?:[-\d]+\s+[-\d]+\s+V\s*\n?)*)'
        
        matches = re.finditer(curve_pattern, self.content, re.MULTILINE)
        
        for match in matches:
            line_type = match.group(1)
            start_x = int(match.group(2))
            start_y = int(match.group(3))
            v_commands = match.group(4)
            
            # Build curve from start point and relative moves
            points = [(start_x, start_y)]
            x, y = start_x, start_y
            
            # Parse all "dx dy V" commands
            v_matches = re.findall(r'([-\d]+)\s+([-\d]+)\s+V', v_commands)
            for dx, dy in v_matches:
                x += int(dx)
                y += int(dy)
                points.append((x, y))
            
            # Create unique curve name by adding counter
            if line_type not in curve_counter:
                curve_counter[line_type] = 0
            curve_counter[line_type] += 1
            
            curve_name = f"{line_type}_{curve_counter[line_type]}"
            curves[curve_name] = points
            
            print(f"  {curve_name}: {len(points)} points, start=({start_x}, {start_y})")
        
        return curves
    
    def get_axis_bounds(self) -> Tuple[float, float, float, float]:
        """Extract axis bounds from PostScript comments or defaults.
        
        Returns:
            (x_min, x_max, y_min, y_max) in PostScript coordinates
        """
        # Try to find axis labels or use defaults
        # For K&I Figure 1, x-axis is log(n_H), y-axis is log(T) or log(P)
        
        # Default bounds (will be refined based on actual data)
        x_min, x_max = 1320, 4713  # Typical PS coordinates for log(n) = -1 to 6
        y_min, y_max = 551, 2913   # Typical PS coordinates for log(T) = 1 to 5
        
        return x_min, x_max, y_min, y_max


def ps_to_physical(x_ps: float, y_ps: float, 
                   x_bounds_ps: Tuple[float, float],
                   y_bounds_ps: Tuple[float, float],
                   x_bounds_phys: Tuple[float, float],
                   y_bounds_phys: Tuple[float, float]) -> Tuple[float, float]:
    """Convert PostScript coordinates to physical units.
    
    Args:
        x_ps, y_ps: PostScript coordinates
        x_bounds_ps: (x_min, x_max) in PostScript
        y_bounds_ps: (y_min, y_max) in PostScript
        x_bounds_phys: (x_min, x_max) in physical units (log scale)
        y_bounds_phys: (y_min, y_max) in physical units (log scale)
    
    Returns:
        (x_phys, y_phys) in physical units (linear scale if log bounds given)
    """
    # Linear interpolation
    x_frac = (x_ps - x_bounds_ps[0]) / (x_bounds_ps[1] - x_bounds_ps[0])
    y_frac = (y_ps - y_bounds_ps[0]) / (y_bounds_ps[1] - y_bounds_ps[0])
    
    # Map to physical range (assumes logarithmic)
    log_x = x_bounds_phys[0] + x_frac * (x_bounds_phys[1] - x_bounds_phys[0])
    log_y = y_bounds_phys[0] + y_frac * (y_bounds_phys[1] - y_bounds_phys[0])
    
    # Convert from log to linear
    x_phys = 10**log_x
    y_phys = 10**log_y
    
    return x_phys, y_phys


def digitize_figure1a(ps_file: str) -> Dict[str, np.ndarray]:
    """Digitize Figure 1a: Temperature and Pressure vs density.
    
    Panel (a) shows:
      - Temperature T vs n_H for N_H = 1e19, 1e20 cm^-2
      - Pressure P/k vs n_H for N_H = 1e19, 1e20 cm^-2
    
    Args:
        ps_file: Path to f1a.ps
    
    Returns:
        Dictionary with keys: 'T_1e19', 'T_1e20', 'P_1e19', 'P_1e20'
        Each value is (n_H, value) array
    """
    parser = PostScriptParser(ps_file)
    curves = parser.extract_curves()
    
    # Figure 1a axis bounds (from examining the PostScript):
    # X-axis: n_H from 10^-1 to 10^6 cm^-3
    # Y-axis (top): Temperature from 10^1 to 10^5 K
    # Y-axis (bottom): Pressure from 10^1 to 10^9 K cm^-3
    
    # PostScript coordinate bounds (measured from file):
    # X: 1320 to 4713 corresponds to log(n_H) = -1 to 6
    # Y: Two different scales for T and P
    
    x_bounds_ps = (1320, 4713)
    x_bounds_phys = (-1, 6)  # log(n_H / cm^-3)
    
    result = {}
    
    # Identify curves by their starting y-position:
    # Around y=850: Temperature curves (top panel, log(T) = 1 to 5)
    # Around y=2800: Pressure curves (bottom panel, log(P) = 1 to 9)
    
    temp_y_bounds_ps = (855, 2238)  # Approximate from tick positions
    temp_y_bounds_phys = (5, 1)  # log(T/K), NOTE: Y-axis is INVERTED in PS!
    
    pres_y_bounds_ps = (2808, 551)  # Pressure axis, also inverted
    pres_y_bounds_phys = (9, 1)  # log(P/k / K cm^-3), INVERTED
    
    # Process each curve
    curve_assignments = []
    for curve_name, points in curves.items():
        if not points:
            continue
        
        # Determine curve type by first y-coordinate
        first_y = points[0][1]
        
        # Convert all points
        n_H_vals = []
        phys_vals = []
        
        # Determine if this is temperature or pressure based on y-position
        if first_y < 1500:  # Temperature curves (upper panel)
            y_bounds_ps = temp_y_bounds_ps
            y_bounds_phys = temp_y_bounds_phys
            curve_type = 'T'
        else:  # Pressure curves (lower panel)
            y_bounds_ps = pres_y_bounds_ps
            y_bounds_phys = pres_y_bounds_phys
            curve_type = 'P'
        
        for x_ps, y_ps in points:
            n_H, phys_val = ps_to_physical(x_ps, y_ps,
                                            x_bounds_ps, y_bounds_ps,
                                            x_bounds_phys, y_bounds_phys)
            n_H_vals.append(n_H)
            phys_vals.append(phys_val)
        
        curve_assignments.append({
            'name': curve_name,
            'type': curve_type,
            'first_y': first_y,
            'n_H': n_H_vals,
            'vals': phys_vals
        })
    
    # Sort curves within each type by first_y to distinguish 1e19 vs 1e20
    temp_curves = sorted([c for c in curve_assignments if c['type'] == 'T'], 
                         key=lambda x: x['first_y'])
    pres_curves = sorted([c for c in curve_assignments if c['type'] == 'P'],
                         key=lambda x: x['first_y'])
    
    # Assign column densities (lower curve first = 1e20, upper = 1e19 typically)
    if len(temp_curves) >= 2:
        result['T_1e19'] = np.column_stack([temp_curves[1]['n_H'], temp_curves[1]['vals']])
        result['T_1e20'] = np.column_stack([temp_curves[0]['n_H'], temp_curves[0]['vals']])
        print(f"  T_1e19: {len(temp_curves[1]['n_H'])} points")
        print(f"    n_H: [{min(temp_curves[1]['n_H']):.2e}, {max(temp_curves[1]['n_H']):.2e}]")
        print(f"    T: [{min(temp_curves[1]['vals']):.2e}, {max(temp_curves[1]['vals']):.2e}]")
        print(f"  T_1e20: {len(temp_curves[0]['n_H'])} points")
        print(f"    n_H: [{min(temp_curves[0]['n_H']):.2e}, {max(temp_curves[0]['n_H']):.2e}]")
        print(f"    T: [{min(temp_curves[0]['vals']):.2e}, {max(temp_curves[0]['vals']):.2e}]")
    
    if len(pres_curves) >= 2:
        result['P_1e19'] = np.column_stack([pres_curves[1]['n_H'], pres_curves[1]['vals']])
        result['P_1e20'] = np.column_stack([pres_curves[0]['n_H'], pres_curves[0]['vals']])
        print(f"  P_1e19: {len(pres_curves[1]['n_H'])} points")
        print(f"  P_1e20: {len(pres_curves[0]['n_H'])} points")
    
    return result


def main():
    """Main digitization routine."""
    ps_dir = Path('/Users/guo/Downloads/sphcode/docs/papers/cooling-heating')
    
    print("="*70)
    print("Digitizing Koyama & Inutsuka (2000) Figure 1")
    print("="*70)
    
    # Figure 1a: Temperature and Pressure
    print("\nProcessing Figure 1a (Temperature & Pressure)...")
    f1a_data = digitize_figure1a(str(ps_dir / 'f1a.ps'))
    
    # Save results
    output_dir = Path('/Users/guo/Downloads/sphcode/data/lane_emden')
    output_dir.mkdir(parents=True, exist_ok=True)
    
    for name, data in f1a_data.items():
        output_file = output_dir / f'koyama_inutsuka_{name}.dat'
        np.savetxt(output_file, data, 
                   header=f'{name}: Column 1 = n_H [cm^-3], Column 2 = {name.split("_")[0]}',
                   fmt='%.6e')
        print(f"  Saved: {output_file}")
    
    print("\n" + "="*70)
    print("Digitization complete!")
    print("="*70)


if __name__ == '__main__':
    main()
