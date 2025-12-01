#!/usr/bin/env python3
"""
Generate C++ header file with digitized Koyama & Inutsuka (2000) data.

WARNING: The digitized data contains the unstable/WNM branch at n~10 cm⁻³,
NOT the CNM branch! See KOYAMA_INUTSUKA_DIGITIZATION_REPORT.md for details.
"""

import numpy as np
from pathlib import Path

def format_array(name, values, items_per_line=8, comment=""):
    """Format a C++ array with proper line breaks."""
    lines = [f"    static constexpr real {name}[{len(values)}] = {{"]
    
    for i, val in enumerate(values):
        if i % items_per_line == 0:
            lines.append("        ")
        
        lines[-1] += f"{val:.6e}"
        
        if i < len(values) - 1:
            lines[-1] += ", "
        
        if (i + 1) % items_per_line == 0 and i < len(values) - 1:
            pass  # Line already ended
    
    lines.append("    };")
    
    if comment:
        return ["    // " + comment] + lines
    return lines

def main():
    data_dir = Path('/Users/guo/Downloads/sphcode/data/lane_emden')
    
    # Load digitized data
    T_1e20 = np.loadtxt(data_dir / 'koyama_inutsuka_T_1e20.dat')
    T_1e19 = np.loadtxt(data_dir / 'koyama_inutsuka_T_1e19.dat')
    P_1e20 = np.loadtxt(data_dir / 'koyama_inutsuka_P_1e20.dat')
    P_1e19 = np.loadtxt(data_dir / 'koyama_inutsuka_P_1e19.dat')
    
    # Generate header file
    header = []
    header.append("/**")
    header.append(" * @file koyama_inutsuka_data_digitized.hpp")
    header.append(" * @brief Pixel-perfect digitization of Koyama & Inutsuka (2000) Figure 1")
    header.append(" * ")
    header.append(" * CRITICAL WARNING:")
    header.append(" * ================")
    header.append(" * At intermediate densities (1 < n_H < 100 cm⁻³), these curves show the")
    header.append(" * UNSTABLE BRANCH, not the Cold Neutral Medium (CNM) branch!")
    header.append(" * ")
    header.append(" * Example: At n_H = 10 cm⁻³:")
    header.append(" *   - This data:  T_eq ≈ 1200 K  (unstable/WNM branch)")
    header.append(" *   - CNM branch: T_eq ≈ 20-50 K (NOT in this data!)")
    header.append(" * ")
    header.append(" * The K&I equilibrium curves are multi-valued (S-shaped) with three")
    header.append(" * equilibrium states at intermediate densities:")
    header.append(" *   1. Warm Neutral Medium (WNM): T ≈ 6000-8000 K")
    header.append(" *   2. Unstable branch: T ≈ 500-2000 K  ← THIS IS WHAT'S IN THE DATA")
    header.append(" *   3. Cold Neutral Medium (CNM): T ≈ 20-100 K  ← NOT in this data!")
    header.append(" * ")
    header.append(" * For CNM simulations, either:")
    header.append(" *   - Use n_H > 100 cm⁻³ where only CNM exists, OR")
    header.append(" *   - Implement multi-phase cooling with explicit branch selection")
    header.append(" * ")
    header.append(" * Data extracted from original PostScript files using digitize_koyama_inutsuka.py")
    header.append(" * See docs/KOYAMA_INUTSUKA_DIGITIZATION_REPORT.md for full analysis.")
    header.append(" */")
    header.append("")
    header.append("#pragma once")
    header.append("")
    header.append("#include \"../defines.hpp\"")
    header.append("#include <cstddef>")
    header.append("")
    header.append("namespace sph {")
    header.append("namespace thermal {")
    header.append("namespace koyama_inutsuka_data_digitized {")
    header.append("")
    
    # Temperature N_H=1e20
    n_H = T_1e20[:, 0]
    T = T_1e20[:, 1]
    idx_n10 = np.argmin(np.abs(n_H - 10.0))
    
    header.append("    // ====== Temperature for N_H = 1e20 cm⁻² ======")
    header.append(f"    // WARNING: At n_H={n_H[idx_n10]:.1e} cm⁻³, T={T[idx_n10]:.0f}K (UNSTABLE BRANCH, not CNM!)")
    header.append(f"    static constexpr size_t N_T_1e20 = {len(T_1e20)};")
    header.extend(format_array("n_T_1e20", n_H, comment="Density [cm⁻³]"))
    header.extend(format_array("val_T_1e20", T, comment="Temperature [K] - CAUTION: unstable branch at n~10!"))
    header.append("")
    
    # Temperature N_H=1e19
    n_H_19 = T_1e19[:, 0]
    T_19 = T_1e19[:, 1]
    
    header.append("    // ====== Temperature for N_H = 1e19 cm⁻² ======")
    header.append(f"    static constexpr size_t N_T_1e19 = {len(T_1e19)};")
    header.extend(format_array("n_T_1e19", n_H_19, comment="Density [cm⁻³]"))
    header.extend(format_array("val_T_1e19", T_19, comment="Temperature [K] - CAUTION: unstable branch at n~10!"))
    header.append("")
    
    # Pressure N_H=1e20
    n_H_P = P_1e20[:, 0]
    P = P_1e20[:, 1]
    
    header.append("    // ====== Pressure for N_H = 1e20 cm⁻² ======")
    header.append(f"    static constexpr size_t N_P_1e20 = {len(P_1e20)};")
    header.extend(format_array("n_P_1e20", n_H_P, comment="Density [cm⁻³]"))
    header.extend(format_array("val_P_1e20", P, comment="Pressure/k [K cm⁻³]"))
    header.append("")
    
    # Pressure N_H=1e19
    n_H_P19 = P_1e19[:, 0]
    P19 = P_1e19[:, 1]
    
    header.append("    // ====== Pressure for N_H = 1e19 cm⁻² ======")
    header.append(f"    static constexpr size_t N_P_1e19 = {len(P_1e19)};")
    header.extend(format_array("n_P_1e19", n_H_P19, comment="Density [cm⁻³]"))
    header.extend(format_array("val_P_1e19", P19, comment="Pressure/k [K cm⁻³]"))
    header.append("")
    
    header.append("} // namespace koyama_inutsuka_data_digitized")
    header.append("} // namespace thermal")
    header.append("} // namespace sph")
    
    # Write to file
    output_file = Path('/Users/guo/Downloads/sphcode/include/thermal/koyama_inutsuka_data_digitized.hpp')
    with open(output_file, 'w') as f:
        f.write('\n'.join(header))
        f.write('\n')
    
    print(f"✓ Generated: {output_file}")
    print(f"  - 4 curves digitized (T_1e19, T_1e20, P_1e19, P_1e20)")
    print(f"  - {len(T_1e20)} points per curve")
    print(f"\n⚠️  WARNING: Data contains unstable branch at n~10 cm⁻³, not CNM!")
    print(f"   At n={n_H[idx_n10]:.1e} cm⁻³: T={T[idx_n10]:.0f}K (should be ~25K for CNM)")
    print(f"\n📖 See docs/KOYAMA_INUTSUKA_DIGITIZATION_REPORT.md for details")

if __name__ == '__main__':
    main()
