"""
K&I 2000 Analysis Package

Tools for reproducing and analyzing heating/cooling curves from 
Koyama & Inutsuka (2000, ApJ 532, 980).

Subpackages:
    extraction - Extract data from PostScript figures
    physics - First-principles heating/cooling calculations
    visualization - Plotting and comparison tools

Usage:
    # Calculate thermal equilibrium from first principles
    from ki2000.physics import ThermalEquilibrium
    eq = ThermalEquilibrium(N_H=1e19)
    result = eq.equilibrium_curve()
    
    # Plot panel C reproduction
    from ki2000.physics import reproduce_ki2000_panel_c
    reproduce_ki2000_panel_c(output_path='panel_c.png', N_H=1e19)
"""

__version__ = '1.0.0'
