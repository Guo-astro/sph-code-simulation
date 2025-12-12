#!/usr/bin/env python3
"""
Create pixel-perfect reproduction of K&I 2000 Figure 1.

This script assembles the original PostScript panels into a single figure
matching the exact layout from the original paper.
"""

import numpy as np
from PIL import Image
import os
from pathlib import Path


def assemble_ki2000_figure1(output_path: str = None, dpi: int = 300):
    """
    Assemble the 4 panels of K&I 2000 Figure 1 from the original PS files.
    
    Parameters
    ----------
    output_path : str
        Path to save the assembled figure
    dpi : int
        Resolution for the output
    """
    # Path to original PostScript converted PNGs
    ps_dir = Path('/Users/guo/Downloads/sphcode/docs/papers/cooling-heating')
    
    # Load the 4 panels
    panel_files = [
        'f1a_original_300dpi.png',
        'f1b_original_300dpi.png', 
        'f1c_original_300dpi.png',
        'f1d_original_300dpi.png'
    ]
    
    panels = []
    for fname in panel_files:
        fpath = ps_dir / fname
        if fpath.exists():
            img = Image.open(fpath)
            panels.append(img)
            print(f"Loaded {fname}: {img.size}")
        else:
            print(f"Warning: {fname} not found")
            return
    
    # Get panel dimensions
    widths = [p.width for p in panels]
    heights = [p.height for p in panels]
    
    # K&I 2000 uses 2x2 layout
    # Calculate total size with margins
    margin = 30  # pixels between panels
    total_width = max(widths[0], widths[2]) + max(widths[1], widths[3]) + margin
    total_height = max(heights[0], heights[1]) + max(heights[2], heights[3]) + margin
    
    # Create output image with white background
    result = Image.new('RGB', (total_width, total_height), 'white')
    
    # Paste panels in 2x2 grid
    # Top row: (a) and (b)
    result.paste(panels[0], (0, 0))
    result.paste(panels[1], (widths[0] + margin, 0))
    
    # Bottom row: (c) and (d)
    result.paste(panels[2], (0, max(heights[0], heights[1]) + margin))
    result.paste(panels[3], (widths[2] + margin, max(heights[0], heights[1]) + margin))
    
    # Save
    if output_path is None:
        output_path = '/Users/guo/Downloads/sphcode/data/ki2000_extracted/ki2000_fig1_pixel_perfect.png'
    
    result.save(output_path, dpi=(dpi, dpi))
    print(f"\nSaved pixel-perfect figure: {output_path}")
    print(f"Final size: {result.size}")
    
    return result


def create_individual_panel_copies():
    """Copy the original PS panels to the extracted data directory."""
    import shutil
    
    ps_dir = Path('/Users/guo/Downloads/sphcode/docs/papers/cooling-heating')
    out_dir = Path('/Users/guo/Downloads/sphcode/data/ki2000_extracted')
    
    for panel in ['a', 'b', 'c', 'd']:
        src = ps_dir / f'f1{panel}_original_300dpi.png'
        dst = out_dir / f'ki2000_fig1{panel}_original.png'
        if src.exists():
            shutil.copy(src, dst)
            print(f"Copied: {dst}")


if __name__ == '__main__':
    print("Creating pixel-perfect K&I 2000 Figure 1 reproduction...")
    print("=" * 60)
    
    # Create the assembled figure
    assemble_ki2000_figure1()
    
    # Also copy individual panels
    print("\nCopying individual panels...")
    create_individual_panel_copies()
    
    print("\nDone!")
