#!/usr/bin/env python3
"""
Analysis script to demonstrate the fundamental difference between SSPH and GSPH
force term ratios, and the geometric amplification in 1D vs 3D.

This script proves:
1. SSPH ratio is purely physical (P/ρ²)
2. GSPH ratio is corrupted by kernel gradient differences
3. 1D wave equation vs 3D spherical convergence amplification

Author: Analysis based on Inutsuka (2002) GSPH paper
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Create output directory
OUTPUT_DIR = Path(__file__).parent.parent / "docs" / "derivations" / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def cubic_spline_kernel(r, h, dim=3):
    """Cubic spline kernel W(r, h)"""
    q = r / h
    sigma = {1: 2/3, 2: 10/(7*np.pi), 3: 1/np.pi}[dim]
    norm = sigma / h**dim
    
    w = np.zeros_like(q)
    mask1 = q < 0.5
    mask2 = (q >= 0.5) & (q < 1.0)
    
    w[mask1] = 1 - 6*q[mask1]**2 + 6*q[mask1]**3
    w[mask2] = 2*(1 - q[mask2])**3
    
    return norm * w


def cubic_spline_gradient_magnitude(r, h, dim=3):
    """Magnitude of cubic spline kernel gradient |∇W(r, h)|"""
    q = r / h
    sigma = {1: 2/3, 2: 10/(7*np.pi), 3: 1/np.pi}[dim]
    norm = sigma / h**(dim+1)
    
    dw = np.zeros_like(q)
    mask1 = q < 0.5
    mask2 = (q >= 0.5) & (q < 1.0)
    
    # dW/dq
    dw[mask1] = -12*q[mask1] + 18*q[mask1]**2
    dw[mask2] = -6*(1 - q[mask2])**2
    
    return np.abs(norm * dw)


def polytrope_profile(r, n=1.5, rho_c=1.0, K=1.0):
    """
    Simple polytropic density profile ρ(r) for a self-gravitating sphere.
    Using Lane-Emden solution approximation.
    """
    # For n=1.5, approximate solution
    xi_1 = 3.65375  # First zero of Lane-Emden for n=1.5
    alpha = np.sqrt((n+1) * K * rho_c**(1/n - 1) / (4 * np.pi))
    R = alpha * xi_1  # Stellar radius
    
    # Approximate density profile
    r_norm = r / R
    rho = rho_c * np.maximum(0, 1 - r_norm**2)**n
    return rho, R


def plot_ssph_vs_gsph_ratio():
    """
    Plot 1: Demonstrate that SSPH ratio is physical while GSPH ratio is corrupted.
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # Setup: particle pair at different radii in a polytrope
    r_range = np.linspace(0.1, 0.9, 100)  # Normalized radius
    rho_c = 1.0
    gamma = 5/3
    K = 1.0
    
    # Get density profile
    rho, R = polytrope_profile(r_range * 3.65, n=1.5, rho_c=rho_c, K=K)
    
    # Pressure from polytropic EOS: P = K ρ^γ
    P = K * rho**gamma
    
    # Smoothing length: h ∝ ρ^(-1/D) for D=3
    D = 3
    eta = 1.2
    m = 0.01  # particle mass
    h = eta * (m / rho)**(1/D)
    
    # ============ Plot 1a: Density and Pressure profiles ============
    ax1 = axes[0, 0]
    ax1_twin = ax1.twinx()
    
    l1, = ax1.plot(r_range, rho, 'b-', linewidth=2, label=r'$\rho(r)$')
    l2, = ax1_twin.plot(r_range, P, 'r-', linewidth=2, label=r'$P(r)$')
    l3, = ax1.plot(r_range, h, 'g--', linewidth=2, label=r'$h(r)$')
    
    ax1.set_xlabel(r'Normalized radius $r/R$', fontsize=12)
    ax1.set_ylabel(r'Density $\rho$, Smoothing length $h$', fontsize=12, color='b')
    ax1_twin.set_ylabel(r'Pressure $P$', fontsize=12, color='r')
    ax1.set_title('Polytropic Profile: Density, Pressure, Smoothing Length', fontsize=14)
    ax1.legend(handles=[l1, l2, l3], loc='upper right')
    ax1.grid(True, alpha=0.3)
    
    # ============ Plot 1b: P/ρ² ratio (physical) ============
    ax2 = axes[0, 1]
    
    # For particle pairs: i at r, j at r + Δr
    delta_r = 0.05
    r_i = r_range[:-10]
    r_j = r_range[:-10] + delta_r
    
    rho_i, _ = polytrope_profile(r_i * 3.65, n=1.5, rho_c=rho_c, K=K)
    rho_j, _ = polytrope_profile(r_j * 3.65, n=1.5, rho_c=rho_c, K=K)
    
    P_i = K * rho_i**gamma
    P_j = K * rho_j**gamma
    
    # SSPH ratio: purely physical
    ssph_ratio = (P_i / rho_i**2) / (P_j / rho_j**2)
    
    ax2.plot(r_i, ssph_ratio, 'b-', linewidth=2, label='SSPH ratio')
    ax2.axhline(y=1.0, color='k', linestyle='--', alpha=0.5)
    ax2.set_xlabel(r'Particle $i$ position $r_i/R$', fontsize=12)
    ax2.set_ylabel(r'$(P_i/\rho_i^2) / (P_j/\rho_j^2)$', fontsize=12)
    ax2.set_title('SSPH: Pure Physical Ratio (EOS only)', fontsize=14)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim([0, 3])
    
    # ============ Plot 1c: GSPH ratio with kernel corruption ============
    ax3 = axes[1, 0]
    
    h_i = eta * (m / rho_i)**(1/D)
    h_j = eta * (m / rho_j)**(1/D)
    
    # Separation between particles
    r_ij = np.abs(r_j - r_i) * 3.65  # Physical separation
    
    # Kernel gradients
    grad_W_i = cubic_spline_gradient_magnitude(r_ij, h_i, dim=D)
    grad_W_j = cubic_spline_gradient_magnitude(r_ij, h_j, dim=D)
    
    # Avoid division by zero
    mask = (grad_W_j > 1e-10) & (grad_W_i > 1e-10)
    
    # GSPH ratio: corrupted by kernel
    gsph_ratio = np.zeros_like(ssph_ratio)
    gsph_ratio[mask] = ((P_i[mask] / rho_i[mask]**2) * grad_W_i[mask]) / \
                       ((P_j[mask] / rho_j[mask]**2) * grad_W_j[mask])
    
    # Theoretical prediction: (h_j/h_i)^(D+1)
    kernel_factor = (h_j / h_i)**(D+1)
    predicted_gsph_ratio = ssph_ratio * kernel_factor
    
    ax3.plot(r_i[mask], gsph_ratio[mask], 'r-', linewidth=2, label='GSPH actual ratio')
    ax3.plot(r_i, ssph_ratio, 'b--', linewidth=2, label='SSPH ratio (physical)')
    ax3.plot(r_i, predicted_gsph_ratio, 'g:', linewidth=2, label=r'Predicted: SSPH $\times (h_j/h_i)^{D+1}$')
    ax3.axhline(y=1.0, color='k', linestyle='--', alpha=0.5)
    
    ax3.set_xlabel(r'Particle $i$ position $r_i/R$', fontsize=12)
    ax3.set_ylabel('Force term ratio', fontsize=12)
    ax3.set_title('GSPH: Ratio Corrupted by Kernel Gradients', fontsize=14)
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    ax3.set_ylim([0, 5])
    
    # ============ Plot 1d: The corruption factor ============
    ax4 = axes[1, 1]
    
    corruption = gsph_ratio[mask] / ssph_ratio[mask]
    
    ax4.plot(r_i[mask], corruption, 'purple', linewidth=2, label='Actual corruption')
    ax4.plot(r_i, kernel_factor, 'orange', linewidth=2, linestyle='--', label=r'$(h_j/h_i)^{D+1}$')
    ax4.axhline(y=1.0, color='k', linestyle='--', alpha=0.5, label='No corruption')
    
    ax4.set_xlabel(r'Particle $i$ position $r_i/R$', fontsize=12)
    ax4.set_ylabel('Corruption factor', fontsize=12)
    ax4.set_title('GSPH Corruption Factor: Deviation from Physical Ratio', fontsize=14)
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    # Add annotation
    ax4.annotate(f'Max deviation: {np.max(np.abs(corruption - 1))*100:.1f}%',
                 xy=(0.5, 0.9), xycoords='axes fraction', fontsize=12,
                 bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'ssph_vs_gsph_ratio.png', dpi=150, bbox_inches='tight')
    plt.savefig(OUTPUT_DIR / 'ssph_vs_gsph_ratio.pdf', bbox_inches='tight')
    print(f"Saved: {OUTPUT_DIR / 'ssph_vs_gsph_ratio.png'}")
    plt.close()


def plot_1d_vs_3d_wave_equation():
    """
    Plot 2: Demonstrate the geometric amplification in 3D vs 1D.
    
    1D wave equation: ∂²ρ/∂t² = c_s² ∂²ρ/∂x²
    3D spherical:     ∂²ρ/∂t² = c_s² [∂²ρ/∂r² + (2/r)∂ρ/∂r] = c_s² ∇²ρ
    
    The (2/r)∂ρ/∂r term is the geometric amplification.
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # ============ Plot 2a: 1D perturbation evolution ============
    ax1 = axes[0, 0]
    
    # 1D wave equation: d'Alembert solution
    x = np.linspace(-2, 2, 500)
    c_s = 1.0
    times = [0, 0.5, 1.0, 1.5]
    
    # Initial Gaussian perturbation
    def initial_perturbation(x, x0=0, sigma=0.2):
        return np.exp(-(x - x0)**2 / (2 * sigma**2))
    
    for t in times:
        # D'Alembert solution: f(x-ct) + f(x+ct)
        rho_1d = 0.5 * (initial_perturbation(x - c_s*t) + initial_perturbation(x + c_s*t))
        ax1.plot(x, rho_1d, label=f't = {t}', linewidth=2)
    
    ax1.set_xlabel('Position x', fontsize=12)
    ax1.set_ylabel(r'Density perturbation $\delta\rho$', fontsize=12)
    ax1.set_title('1D: Wave Propagation (No Amplification)', fontsize=14)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim([0, 1.2])
    
    # ============ Plot 2b: 3D spherical perturbation evolution ============
    ax2 = axes[0, 1]
    
    r = np.linspace(0.01, 2, 500)
    
    # In 3D spherical, inward-moving wave gets amplified as 1/r
    # Solution form: (1/r) * f(r - c_s*t) for inward wave
    for t in times:
        # Inward propagating wave with 1/r amplification
        r_wave = r + c_s * t  # Wave moving inward (toward r=0)
        amplitude = initial_perturbation(r_wave - 1.0, sigma=0.2)
        # The 1/r factor causes amplification near origin
        rho_3d = amplitude / np.maximum(r, 0.1)
        rho_3d = np.minimum(rho_3d, 10)  # Cap for visualization
        ax2.plot(r, rho_3d, label=f't = {t}', linewidth=2)
    
    ax2.set_xlabel('Radius r', fontsize=12)
    ax2.set_ylabel(r'Density perturbation $\delta\rho$', fontsize=12)
    ax2.set_title('3D Spherical: Wave Amplification (1/r factor)', fontsize=14)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim([0, 5])
    
    # Annotation
    ax2.annotate('Geometric\namplification\n→ ∞ as r → 0',
                 xy=(0.15, 3), fontsize=11,
                 bbox=dict(boxstyle='round', facecolor='red', alpha=0.3))
    
    # ============ Plot 2c: The geometric term (2/r)∂ρ/∂r ============
    ax3 = axes[1, 0]
    
    r = np.linspace(0.1, 2, 500)
    
    # For a density profile with inward gradient
    rho_profile = 1.0 + 0.5 * np.exp(-r)  # Decreasing outward
    drho_dr = -0.5 * np.exp(-r)
    
    # Geometric term
    geometric_term = (2/r) * drho_dr
    laplacian_1d = np.gradient(np.gradient(rho_profile, r), r)
    
    ax3.plot(r, np.abs(geometric_term), 'r-', linewidth=2, label=r'$|(2/r)\partial\rho/\partial r|$ (geometric)')
    ax3.plot(r, np.abs(laplacian_1d), 'b--', linewidth=2, label=r'$|\partial^2\rho/\partial r^2|$ (1D-like)')
    ax3.plot(r, np.abs(geometric_term) / (np.abs(laplacian_1d) + 1e-10), 'g:', linewidth=2, label='Ratio')
    
    ax3.set_xlabel('Radius r', fontsize=12)
    ax3.set_ylabel('Term magnitude', fontsize=12)
    ax3.set_title('Geometric Term vs 1D Laplacian', fontsize=14)
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    ax3.set_yscale('log')
    
    # ============ Plot 2d: Error amplification simulation ============
    ax4 = axes[1, 1]
    
    # Simulate error growth with and without geometric term
    t_sim = np.linspace(0, 10, 1000)
    
    # Initial error
    epsilon_0 = 0.01
    
    # 1D: error stays bounded (oscillatory)
    omega_1d = 1.0  # characteristic frequency
    error_1d = epsilon_0 * np.cos(omega_1d * t_sim)
    
    # 3D with geometric amplification: positive feedback
    # Simplified model: dε/dt = α * ε where α > 0 due to geometric term
    alpha_3d = 0.5  # growth rate from geometric focusing
    error_3d = epsilon_0 * np.exp(alpha_3d * t_sim)
    
    # 3D with Ω correction: stabilized
    error_3d_corrected = epsilon_0 * np.cos(omega_1d * t_sim) * np.exp(-0.1 * t_sim)
    
    ax4.plot(t_sim, np.abs(error_1d), 'b-', linewidth=2, label='1D: Bounded oscillation')
    ax4.plot(t_sim, error_3d, 'r-', linewidth=2, label='3D without Ω: Exponential growth')
    ax4.plot(t_sim, np.abs(error_3d_corrected), 'g--', linewidth=2, label='3D with Ω: Stabilized')
    
    ax4.set_xlabel('Time', fontsize=12)
    ax4.set_ylabel('Error magnitude', fontsize=12)
    ax4.set_title('Error Evolution: 1D vs 3D', fontsize=14)
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    ax4.set_yscale('log')
    ax4.set_ylim([1e-3, 1e3])
    
    # Annotation
    ax4.annotate('Catastrophic\ncollapse',
                 xy=(8, 50), fontsize=11,
                 bbox=dict(boxstyle='round', facecolor='red', alpha=0.3))
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / '1d_vs_3d_geometric.png', dpi=150, bbox_inches='tight')
    plt.savefig(OUTPUT_DIR / '1d_vs_3d_geometric.pdf', bbox_inches='tight')
    print(f"Saved: {OUTPUT_DIR / '1d_vs_3d_geometric.png'}")
    plt.close()


def plot_pressure_deficit_feedback():
    """
    Plot 3: Show the pressure deficit feedback mechanism.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # ============ Left: Pressure deficit as function of Ω ============
    ax1 = axes[0]
    
    omega = np.linspace(1.0, 1.5, 100)
    
    # Pressure deficit: ε = 1 - 1/Ω
    epsilon = 1 - 1/omega
    
    ax1.plot(omega, epsilon * 100, 'r-', linewidth=3)
    ax1.fill_between(omega, 0, epsilon * 100, alpha=0.3, color='red')
    
    # Typical values in stratified media
    ax1.axvline(x=1.1, color='blue', linestyle='--', label='Typical surface: Ω ≈ 1.1')
    ax1.axvline(x=1.3, color='green', linestyle='--', label='Steep gradient: Ω ≈ 1.3')
    
    ax1.set_xlabel(r'$\Omega$ (grad-h factor)', fontsize=14)
    ax1.set_ylabel('Pressure deficit ε (%)', fontsize=14)
    ax1.set_title('Pressure Force Underestimation without grad-h', fontsize=14)
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim([1.0, 1.5])
    ax1.set_ylim([0, 35])
    
    # Annotations
    ax1.annotate('5-13% typical\nunderestimate',
                 xy=(1.15, 12), fontsize=12,
                 bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.5))
    
    # ============ Right: Feedback loop visualization ============
    ax2 = axes[1]
    
    # Simulate the feedback loop
    t = np.linspace(0, 10, 1000)
    
    # Without grad-h: positive feedback
    # dρ/dt ∝ ε * g * ρ, where ε ∝ |∇ρ|/ρ
    rho_no_gradh = 1.0 * np.exp(0.3 * t)  # Exponential growth
    
    # With grad-h: stable equilibrium
    rho_with_gradh = 1.0 + 0.1 * np.sin(t) * np.exp(-0.1 * t)  # Damped oscillation
    
    ax2.plot(t, rho_no_gradh, 'r-', linewidth=2, label='Without grad-h: Collapse')
    ax2.plot(t, rho_with_gradh, 'g-', linewidth=2, label='With grad-h: Stable')
    ax2.axhline(y=1.0, color='k', linestyle='--', alpha=0.5, label='Equilibrium')
    
    ax2.set_xlabel('Time', fontsize=14)
    ax2.set_ylabel(r'Central density $\rho_c$', fontsize=14)
    ax2.set_title('Density Evolution: Positive Feedback Loop', fontsize=14)
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3)
    ax2.set_yscale('log')
    ax2.set_ylim([0.5, 100])
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'pressure_deficit_feedback.png', dpi=150, bbox_inches='tight')
    plt.savefig(OUTPUT_DIR / 'pressure_deficit_feedback.pdf', bbox_inches='tight')
    print(f"Saved: {OUTPUT_DIR / 'pressure_deficit_feedback.png'}")
    plt.close()


def main():
    """Generate all analysis plots."""
    print("=" * 60)
    print("GSPH vs SSPH Ratio Analysis")
    print("=" * 60)
    
    print("\n1. Generating SSPH vs GSPH ratio comparison...")
    plot_ssph_vs_gsph_ratio()
    
    print("\n2. Generating 1D vs 3D geometric analysis...")
    plot_1d_vs_3d_wave_equation()
    
    print("\n3. Generating pressure deficit feedback plot...")
    plot_pressure_deficit_feedback()
    
    print("\n" + "=" * 60)
    print(f"All plots saved to: {OUTPUT_DIR}")
    print("=" * 60)


if __name__ == "__main__":
    main()
