"""
SRGSPH Particle Class

Represents a single SPH particle with Special Relativistic Godunov formulation.
Uses the Kitajima formulation with baryon number density.
"""

import numpy as np


class Particle:
    """
    Special Relativistic SPH Particle.
    
    Attributes:
        # Particle properties
        nu: Baryon number in this particle
        x: Position in 1D
        
        # Conserved quantities (time-evolved)
        S: Canonical momentum per baryon = γHv
        e: Canonical energy per baryon = γH - P/(Nc²)
        
        # Primitive variables (recovered from S, e)
        v: Velocity
        n: Baryon number density in rest frame
        P: Pressure
        u: Thermal energy per baryon
        
        # Derived quantities
        gamma: Lorentz factor = 1/sqrt(1 - v²/c²)
        H: Enthalpy per baryon = 1 + u/c² + P/(nc²)
        N: Baryon number density in lab frame = γn
        cs: Sound speed
        
        # SPH quantities
        h: Smoothing length
        Vp: Particle volume
        neighbors: List of neighbor particle indices
        
        # Gradient quantities (for MUSCL reconstruction)
        grad_n: Gradient of rest frame density
        grad_P: Gradient of pressure
        grad_v: Gradient of velocity
    """
    
    def __init__(self, nu, x, v, n, P, gamma_c, c=1.0):
        """
        Initialize particle with primitive variables.
        
        Args:
            nu: Baryon number in particle
            x: Position
            v: Velocity
            n: Rest frame baryon density
            P: Pressure
            gamma_c: Adiabatic index
            c: Speed of light
        """
        self.nu = nu
        self.x = x
        self.gamma_c = gamma_c
        self.c = c
        
        # Set primitive variables
        self.v = v
        self.n = n
        self.P = P
        
        # Compute thermal energy: u = P/((γ_c - 1)n)
        self.u = P / ((gamma_c - 1.0) * n)
        
        # Compute derived quantities
        self.gamma = 1.0 / np.sqrt(1.0 - (v/c)**2)
        self.H = 1.0 + self.u/(c**2) + P/(n * c**2)
        self.N = self.gamma * n
        self.cs = np.sqrt(gamma_c * P / (n * self.H))
        
        # Compute conserved quantities
        # S = γHv
        self.S = self.gamma * self.H * v
        # e = γH - P/(Nc²)
        self.e = self.gamma * self.H - P / (self.N * c**2)
        
        # SPH quantities (to be set later)
        self.h = 0.0
        self.Vp = 0.0
        self.neighbors = []
        
        # Gradients (for MUSCL)
        self.grad_n = 0.0
        self.grad_P = 0.0
        self.grad_v = 0.0
        
    def recover_primitives(self):
        """
        Recover primitive variables (v, n, P, u) from conserved (S, e).
        
        Uses the quartic equation for γ from Eq. 67 in the paper:
        (γ² - 1)(Xe γ - 1)² - S²(Xγ² - 1)² = 0
        where X = γ_c/(γ_c - 1)
        
        Then recovers velocity from Eq. 69.
        """
        X = self.gamma_c / (self.gamma_c - 1.0)
        c = self.c
        
        # Solve quartic equation for γ
        # Using numpy polynomial solver
        # Expand: (γ² - 1)(Xe γ - 1)² - S²(Xγ² - 1)² = 0
        
        # Coefficients for polynomial in γ
        # After expansion, we get a quartic in γ
        a4 = X**2 * self.e**2 - self.S**2 * X**2
        a3 = -2.0 * X * self.e
        a2 = 1.0 - self.S**2 * (2.0 * X - 1.0)
        a1 = 0.0
        a0 = self.S**2
        
        coeffs = [a4, a3, a2, a1, a0]
        roots = np.roots(coeffs)
        
        # Select physical root (real and γ >= 1)
        gamma_new = None
        for root in roots:
            if np.isreal(root) and np.real(root) >= 1.0:
                gamma_new = np.real(root)
                break
        
        if gamma_new is None:
            # Fallback: use Newton iteration
            gamma_new = self._solve_gamma_newton(X)
        
        self.gamma = gamma_new
        
        # Recover velocity from Eq. 69
        # v = (Xγ² - 1) / (γ(Xe γ - 1)) * S
        denominator = self.gamma * (X * self.e * self.gamma - 1.0)
        if abs(denominator) < 1e-15:
            self.v = 0.0
        else:
            self.v = (X * self.gamma**2 - 1.0) / denominator * self.S
        
        # Enforce causality: |v| < c
        v_max = 0.99 * c  # Slightly below c for numerical stability
        if abs(self.v) > v_max:
            self.v = np.sign(self.v) * v_max
            # Recompute gamma with limited velocity
            self.gamma = 1.0 / np.sqrt(1.0 - (self.v/c)**2)
        
        # Recover enthalpy
        # H = (Xe γ - 1) / (Xγ² - 1)
        self.H = (X * self.e * self.gamma - 1.0) / (X * self.gamma**2 - 1.0)
        
        # Recover rest frame density from N = γn
        self.n = self.N / self.gamma
        
        # Recover thermal energy from H = 1 + u/c² + P/(nc²)
        # We need to use equation of state: P = (γ_c - 1)nu
        # H = 1 + u/c² + (γ_c - 1)u/(c²) = 1 + γ_c u/c²
        self.u = (self.H - 1.0) * c**2 / self.gamma_c
        
        # Recover pressure
        self.P = (self.gamma_c - 1.0) * self.n * self.u
        
        # Ensure positive pressure and density
        if self.P < 1e-10:
            self.P = 1e-10
        if self.n < 1e-10:
            self.n = 1e-10
        
        # Update sound speed
        # cs² = γ_c P / (n H)
        cs_squared = self.gamma_c * self.P / (self.n * self.H)
        if cs_squared < 0:
            self.cs = 0.0
        else:
            self.cs = np.sqrt(cs_squared)
        
    def _solve_gamma_newton(self, X, max_iter=100, tol=1e-12):
        """
        Solve for γ using Newton-Raphson iteration.
        
        Args:
            X: γ_c/(γ_c - 1)
            max_iter: Maximum iterations
            tol: Tolerance
            
        Returns:
            γ (Lorentz factor)
        """
        # Initial guess
        gamma = 1.0 + abs(self.S) / self.c
        
        for _ in range(max_iter):
            # f(γ) = (γ² - 1)(Xe γ - 1)² - S²(Xγ² - 1)²
            f = ((gamma**2 - 1.0) * (X * self.e * gamma - 1.0)**2 -
                 self.S**2 * (X * gamma**2 - 1.0)**2)
            
            # f'(γ)
            df = (2.0 * gamma * (X * self.e * gamma - 1.0)**2 +
                  2.0 * (gamma**2 - 1.0) * (X * self.e * gamma - 1.0) * X * self.e -
                  2.0 * self.S**2 * (X * gamma**2 - 1.0) * 2.0 * X * gamma)
            
            if abs(df) < 1e-20:
                break
            
            gamma_new = gamma - f / df
            
            # Ensure γ >= 1
            gamma_new = max(gamma_new, 1.0)
            
            if abs(gamma_new - gamma) < tol:
                return gamma_new
            
            gamma = gamma_new
        
        return gamma
    
    def __repr__(self):
        return (f"Particle(x={self.x:.4f}, v={self.v:.4f}, n={self.n:.4f}, "
                f"P={self.P:.4f}, gamma={self.gamma:.4f})")
