#!/usr/bin/env python3
"""
Generate C++ source file with hardcoded Koyama & Inutsuka (2000) data.
All 19 curves embedded as static arrays - no external file dependencies.
"""

import numpy as np

def generate_cpp_array(varname, description, n_data, val_data):
    """Generate C++ static array initialization."""
    n_points = len(n_data)
    
    # Start array
    code = f"    // {description} ({n_points} points)\n"
    code += f"    static constexpr size_t N_{varname} = {n_points};\n"
    code += f"    static constexpr real n_{varname}[N_{varname}] = {{\n"
    
    # Density values (8 per line)
    for i in range(0, n_points, 8):
        end = min(i + 8, n_points)
        values = ', '.join([f'{n:.6e}' for n in n_data[i:end]])
        code += f"        {values}"
        if end < n_points:
            code += ","
        code += "\n"
    code += "    };\n"
    
    # Value array
    code += f"    static constexpr real val_{varname}[N_{varname}] = {{\n"
    for i in range(0, n_points, 8):
        end = min(i + 8, n_points)
        values = ', '.join([f'{v:.6e}' for v in val_data[i:end]])
        code += f"        {values}"
        if end < n_points:
            code += ","
        code += "\n"
    code += "    };\n\n"
    
    return code

def main():
    print("="*70)
    print("GENERATING C++ HARDCODED KOYAMA-INUTSUKA DATA")
    print("="*70)
    
    # Define all curves to process
    curves = [
        ('f1a_0', 'T_1e19', 'Temperature N_H=1e19 cm^-2 [K]'),
        ('f1a_1', 'T_1e20', 'Temperature N_H=1e20 cm^-2 [K]'),
        ('f1a_2', 'P_1e19', 'Pressure N_H=1e19 cm^-2 [K cm^-3]'),
        ('f1a_3', 'P_1e20', 'Pressure N_H=1e20 cm^-2 [K cm^-3]'),
        ('f1b_0', 'xe_1e19', 'Electron fraction N_H=1e19'),
        ('f1b_1', 'xe_1e20', 'Electron fraction N_H=1e20'),
        ('f1b_2', 'xH2', 'H2 molecular fraction'),
        ('f1b_3', 'xCO', 'CO molecular fraction'),
        ('f1c_0', 'Gamma_PE', 'Photoelectric heating [erg/s/H]'),
        ('f1c_1', 'Gamma_CR', 'Cosmic ray heating [erg/s/H]'),
        ('f1c_2', 'Gamma_H2', 'H2 heating [erg/s/H]'),
        ('f1c_3', 'Lambda_CII', 'CII cooling [erg/s/H]'),
        ('f1c_4', 'Lambda_OI', 'OI cooling [erg/s/H]'),
        ('f1c_5', 'Lambda_Lya', 'Lya cooling [erg/s/H]'),
        ('f1c_6', 'Lambda_CO', 'CO cooling [erg/s/H]'),
        ('f1d_0', 't_cool', 'Cooling timescale [years]'),
        ('f1d_1', 't_rec', 'Recombination timescale [years]'),
        ('f1d_2', 't_ff', 'Free-fall timescale [years]'),
        ('f1d_3', 't_H2', 'H2 formation timescale [years]'),
    ]
    
    # Generate all arrays
    arrays_code = ""
    for file_key, varname, desc in curves:
        panel, idx = file_key.split('_')
        filename = f'../results/{panel}_curve_{idx}.txt'
        
        try:
            data = np.loadtxt(filename)
            # Reverse if needed (data might be high-to-low)
            if data[0,0] > data[-1,0]:
                data = data[::-1]
            
            n_data = data[:, 0]
            val_data = data[:, 1]
            
            arrays_code += generate_cpp_array(varname, desc, n_data, val_data)
            print(f"✓ {varname}: {len(n_data)} points")
            
        except Exception as e:
            print(f"✗ Failed to load {filename}: {e}")
    
    # Write complete C++ file
    cpp_code = f'''/**
 * @file koyama_inutsuka_data.hpp
 * @brief Hardcoded interpolation data from Koyama & Inutsuka (2000) Figure 1
 * 
 * All 19 curves from the paper are embedded as static constexpr arrays.
 * No external file dependencies - completely self-contained.
 * 
 * Data extracted pixel-perfectly from original PostScript files:
 * - Panel (a): Temperature and Pressure vs density (4 curves)
 * - Panel (b): Chemical fractions vs density (4 curves) 
 * - Panel (c): Heating/cooling rates vs density (7 curves)
 * - Panel (d): Timescales vs density (4 curves)
 * 
 * Total: {sum([len(np.loadtxt(f'../results/{c[0].split("_")[0]}_curve_{c[0].split("_")[1]}.txt')) for c in curves])} data points across 19 curves
 */

#pragma once

#include "../../defines.hpp"
#include <cstddef>

namespace sph {{
namespace thermal {{
namespace koyama_inutsuka_data {{

{arrays_code}
}}  // namespace koyama_inutsuka_data
}}  // namespace thermal
}}  // namespace sph
'''
    
    output_file = '../../../include/thermal/koyama_inutsuka_data.hpp'
    with open(output_file, 'w') as f:
        f.write(cpp_code)
    
    print("\n" + "="*70)
    print(f"✓ Generated: {output_file}")
    print("="*70)
    print("\nUsage in C++:")
    print("  #include \"thermal/koyama_inutsuka_data.hpp\"")
    print("  using namespace sph::thermal::koyama_inutsuka_data;")
    print("  // Access arrays: n_T_1e19[i], val_T_1e19[i], etc.")
    print("="*70)

if __name__ == '__main__':
    main()
