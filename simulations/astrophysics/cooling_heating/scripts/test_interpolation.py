#!/usr/bin/env python3
"""
Test interpolation-based thermal equilibrium for SPH simulations.
Compares interpolated T(n) and P(n) with exact digitized data.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

class KoyamaInutsukaCoolingPython:
    """Python version of the C++ cooling class for testing."""
    
    def __init__(self, data_dir, N_H_column=1e19):
        self.N_H_column = N_H_column
        
        # Load appropriate curves
        if N_H_column >= 5e19:
            T_file = f'{data_dir}/f1a_curve_1.txt'
            P_file = f'{data_dir}/f1a_curve_3.txt'
            print(f'Loading N_H = 1e20 cm^-2 curves')
        else:
            T_file = f'{data_dir}/f1a_curve_0.txt'
            P_file = f'{data_dir}/f1a_curve_2.txt'
            print(f'Loading N_H = 1e19 cm^-2 curves')
        
        # Load data
        self.T_data = np.loadtxt(T_file)
        self.P_data = np.loadtxt(P_file)
        
        # Reverse if needed (data is high-to-low)
        if self.T_data[0,0] > self.T_data[-1,0]:
            self.T_data = self.T_data[::-1]
        if self.P_data[0,0] > self.P_data[-1,0]:
            self.P_data = self.P_data[::-1]
        
        print(f'  Temperature: {len(self.T_data)} points, n=[{self.T_data[0,0]:.2e}, {self.T_data[-1,0]:.2e}] cm^-3')
        print(f'  T range: [{self.T_data[0,1]:.1f}, {self.T_data[-1,1]:.1f}] K')
        print(f'  Pressure: {len(self.P_data)} points')
        print(f'  P range: [{self.P_data[0,1]:.1f}, {self.P_data[-1,1]:.1f}] K cm^-3')
    
    def temperature(self, n_H):
        """Interpolate temperature in log-log space."""
        if np.isscalar(n_H):
            return np.interp(np.log10(n_H), 
                           np.log10(self.T_data[:,0]), 
                           np.log10(self.T_data[:,1]))
        else:
            log_T = np.interp(np.log10(n_H), 
                            np.log10(self.T_data[:,0]), 
                            np.log10(self.T_data[:,1]))
            return 10**log_T
    
    def pressure(self, n_H):
        """Interpolate pressure in log-log space."""
        if np.isscalar(n_H):
            return np.interp(np.log10(n_H), 
                           np.log10(self.P_data[:,0]), 
                           np.log10(self.P_data[:,1]))
        else:
            log_P = np.interp(np.log10(n_H), 
                            np.log10(self.P_data[:,0]), 
                            np.log10(self.P_data[:,1]))
            return 10**log_P

def test_interpolation():
    """Test the interpolation and create comparison plots."""
    
    print("="*70)
    print("TESTING KOYAMA-INUTSUKA COOLING INTERPOLATION")
    print("="*70)
    
    # Load cooling curves
    data_dir = '../results'
    cooling_1e19 = KoyamaInutsukaCoolingPython(data_dir, N_H_column=1e19)
    cooling_1e20 = KoyamaInutsukaCoolingPython(data_dir, N_H_column=1e20)
    
    # Test points
    n_test = np.logspace(-1, 6, 1000)
    
    # Interpolate
    T_interp_1e19 = cooling_1e19.temperature(n_test)
    P_interp_1e19 = cooling_1e19.pressure(n_test)
    T_interp_1e20 = cooling_1e20.temperature(n_test)
    P_interp_1e20 = cooling_1e20.pressure(n_test)
    
    # Create plots
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    
    # Temperature - N_H = 1e19
    ax1 = axes[0, 0]
    ax1.loglog(cooling_1e19.T_data[:,0], cooling_1e19.T_data[:,1], 'bo', 
              markersize=4, alpha=0.6, label='Exact data points')
    ax1.loglog(n_test, T_interp_1e19, 'r-', linewidth=2, 
              label='Interpolation')
    ax1.set_xlabel(r'$n$ [cm$^{-3}$]', fontsize=12, fontweight='bold')
    ax1.set_ylabel(r'$T$ [K]', fontsize=12, fontweight='bold')
    ax1.set_xlim(0.1, 1e6)
    ax1.set_ylim(10, 1e5)
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3, which='both')
    ax1.set_title(r'Temperature: $N_H=10^{19}$ cm$^{-2}$', fontsize=13, fontweight='bold')
    
    # Temperature - N_H = 1e20
    ax2 = axes[0, 1]
    ax2.loglog(cooling_1e20.T_data[:,0], cooling_1e20.T_data[:,1], 'bo', 
              markersize=4, alpha=0.6, label='Exact data points')
    ax2.loglog(n_test, T_interp_1e20, 'r-', linewidth=2, 
              label='Interpolation')
    ax2.set_xlabel(r'$n$ [cm$^{-3}$]', fontsize=12, fontweight='bold')
    ax2.set_ylabel(r'$T$ [K]', fontsize=12, fontweight='bold')
    ax2.set_xlim(0.1, 1e6)
    ax2.set_ylim(10, 1e5)
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3, which='both')
    ax2.set_title(r'Temperature: $N_H=10^{20}$ cm$^{-2}$', fontsize=13, fontweight='bold')
    
    # Pressure - N_H = 1e19
    ax3 = axes[1, 0]
    ax3.loglog(cooling_1e19.P_data[:,0], cooling_1e19.P_data[:,1], 'bo', 
              markersize=4, alpha=0.6, label='Exact data points')
    ax3.loglog(n_test, P_interp_1e19, 'r-', linewidth=2, 
              label='Interpolation')
    ax3.set_xlabel(r'$n$ [cm$^{-3}$]', fontsize=12, fontweight='bold')
    ax3.set_ylabel(r'$P/k_B$ [K cm$^{-3}$]', fontsize=12, fontweight='bold')
    ax3.set_xlim(0.1, 1e6)
    ax3.set_ylim(100, 1e8)
    ax3.legend(fontsize=11)
    ax3.grid(True, alpha=0.3, which='both')
    ax3.set_title(r'Pressure: $N_H=10^{19}$ cm$^{-2}$', fontsize=13, fontweight='bold')
    
    # Pressure - N_H = 1e20
    ax4 = axes[1, 1]
    ax4.loglog(cooling_1e20.P_data[:,0], cooling_1e20.P_data[:,1], 'bo', 
              markersize=4, alpha=0.6, label='Exact data points')
    ax4.loglog(n_test, P_interp_1e20, 'r-', linewidth=2, 
              label='Interpolation')
    ax4.set_xlabel(r'$n$ [cm$^{-3}$]', fontsize=12, fontweight='bold')
    ax4.set_ylabel(r'$P/k_B$ [K cm$^{-3}$]', fontsize=12, fontweight='bold')
    ax4.set_xlim(0.1, 1e6)
    ax4.set_ylim(100, 1e8)
    ax4.legend(fontsize=11)
    ax4.grid(True, alpha=0.3, which='both')
    ax4.set_title(r'Pressure: $N_H=10^{20}$ cm$^{-2}$', fontsize=13, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig('../results/interpolation_test.png', dpi=300, bbox_inches='tight')
    print("\n✓ Saved: ../results/interpolation_test.png")
    
    # Print test values
    print("\n" + "="*70)
    print("INTERPOLATION TEST VALUES (N_H = 1e19 cm^-2)")
    print("="*70)
    print(f"{'n [cm^-3]':>12} {'T [K]':>12} {'P/k_B [K cm^-3]':>20}")
    print("-"*70)
    for n in [0.1, 1.0, 10.0, 100.0, 1000.0]:
        T = cooling_1e19.temperature(n)
        P = cooling_1e19.pressure(n)
        print(f"{n:12.1f} {T:12.1f} {P:20.1f}")
    
    print("\n" + "="*70)
    print("✓ INTERPOLATION TEST COMPLETE")
    print("="*70)

if __name__ == '__main__':
    test_interpolation()
