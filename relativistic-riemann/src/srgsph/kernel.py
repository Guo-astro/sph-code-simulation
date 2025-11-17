"""
SPH Kernel Functions for SRGSPH

Implements Gaussian kernel as specified in Eq. 11 of the paper.
"""

import numpy as np


class GaussianKernel:
    """
    Gaussian kernel for SPH.
    
    For 1D: W(x, h) = (1/(h√π)) exp(-x²/h²)
    
    From Eq. 11 in the paper.
    """
    
    def __init__(self, dim=1):
        """
        Initialize kernel.
        
        Args:
            dim: Number of dimensions (1, 2, or 3)
        """
        self.dim = dim
        
        # Normalization constant
        self.norm = 1.0 / np.sqrt(np.pi)**dim
    
    def W(self, r, h):
        """
        Kernel function W(r, h).
        
        Args:
            r: Distance (scalar or array)
            h: Smoothing length (scalar or array, same shape as r)
            
        Returns:
            Kernel value
        """
        q = r / h
        return self.norm / h**self.dim * np.exp(-q**2)
    
    def gradW(self, dx, h):
        """
        Gradient of kernel ∇W(x, h).
        
        For 1D: dW/dx = -2x/h³ W(x, h) = -2x/(h²) * W(x,h) / h
        
        Args:
            dx: Displacement vector (for 1D, this is just Δx)
            h: Smoothing length
            
        Returns:
            Gradient (same shape as dx)
        """
        r = abs(dx)
        q = r / h
        
        # dW/dr
        dW_dr = -2.0 * q / h * self.W(r, h)
        
        # For 1D: ∇W = (dW/dr) * (dx/r)
        # Handle r=0 case
        if np.isscalar(dx):
            if abs(dx) < 1e-15:
                return 0.0
            else:
                return dW_dr * np.sign(dx)
        else:
            grad = np.zeros_like(dx)
            mask = abs(dx) > 1e-15
            grad[mask] = dW_dr[mask] * np.sign(dx[mask])
            return grad
    
    def kernel_at_midpoint(self, xi, xj, hi, hj):
        """
        Kernel at midpoint with averaged smoothing length.
        
        This is W(xi - xj, 2h̄) where h̄ = (hi + hj)/2
        Used in SRGSPH force calculations (Eq. 64-65).
        
        Args:
            xi, xj: Particle positions
            hi, hj: Smoothing lengths
            
        Returns:
            W(xi - xj, 2h̄)
        """
        h_avg = 0.5 * (hi + hj)
        dx = xi - xj
        r = abs(dx)
        return self.W(r, 2.0 * h_avg)
    
    def grad_kernel_pair(self, xi, xj, hi, hj):
        """
        Gradient terms for particle pair interaction.
        
        Returns [∇i W(xi - xj, 2hi) - ∇j W(xi - xj, 2hj)]
        as used in Eq. 64-65 for variable smoothing length.
        
        Args:
            xi, xj: Particle positions
            hi, hj: Smoothing lengths
            
        Returns:
            (grad_i, grad_j): Gradient contributions
        """
        dx = xi - xj
        
        # ∇i W(xi - xj, 2hi)
        grad_i = self.gradW(dx, 2.0 * hi)
        
        # ∇j W(xi - xj, 2hj) = -∇i W(xj - xi, 2hj)
        grad_j = self.gradW(-dx, 2.0 * hj)
        
        return grad_i, grad_j


class WendlandC2Kernel:
    """
    Wendland C2 kernel (compact support).
    
    Alternative to Gaussian kernel, mentioned in the paper.
    For 1D: W(q) = (21/(16π)) (1-q)³(1+3q) for q < 1
    """
    
    def __init__(self, dim=1):
        """
        Initialize kernel.
        
        Args:
            dim: Number of dimensions
        """
        self.dim = dim
        
        # Normalization constants
        if dim == 1:
            self.sigma = 5.0 / 4.0
        elif dim == 2:
            self.sigma = 7.0 / (4.0 * np.pi)
        elif dim == 3:
            self.sigma = 21.0 / (16.0 * np.pi)
        else:
            raise ValueError(f"Unsupported dimension: {dim}")
    
    def W(self, r, h):
        """
        Kernel function.
        
        Args:
            r: Distance
            h: Smoothing length
            
        Returns:
            Kernel value
        """
        q = r / h
        
        if np.isscalar(q):
            if q >= 1.0:
                return 0.0
            else:
                return self.sigma / h**self.dim * (1.0 - q)**3 * (1.0 + 3.0 * q)
        else:
            W = np.zeros_like(q)
            mask = q < 1.0
            W[mask] = self.sigma / h**self.dim * ((1.0 - q[mask])**3 * 
                                                   (1.0 + 3.0 * q[mask]))
            return W
    
    def gradW(self, dx, h):
        """
        Gradient of kernel.
        
        Args:
            dx: Displacement
            h: Smoothing length
            
        Returns:
            Gradient
        """
        r = abs(dx)
        q = r / h
        
        # dW/dq = sigma * [-4(1-q)²(1+3q) + 3(1-q)³]
        #       = sigma * (1-q)² [-4(1+3q) + 3(1-q)]
        #       = sigma * (1-q)² [-4 - 12q + 3 - 3q]
        #       = sigma * (1-q)² [-1 - 15q]
        
        if np.isscalar(q):
            if q >= 1.0:
                return 0.0
            else:
                dW_dq = self.sigma / h**self.dim * (1.0 - q)**2 * (-1.0 - 15.0 * q)
                dW_dr = dW_dq / h
                if abs(dx) < 1e-15:
                    return 0.0
                else:
                    return dW_dr * np.sign(dx)
        else:
            grad = np.zeros_like(dx)
            mask = (abs(dx) > 1e-15) & (q < 1.0)
            dW_dq = self.sigma / h**self.dim * (1.0 - q[mask])**2 * (-1.0 - 15.0 * q[mask])
            dW_dr = dW_dq / h
            grad[mask] = dW_dr * np.sign(dx[mask])
            return grad
    
    def kernel_at_midpoint(self, xi, xj, hi, hj):
        """Kernel at midpoint with averaged smoothing length."""
        h_avg = 0.5 * (hi + hj)
        dx = xi - xj
        r = abs(dx)
        return self.W(r, 2.0 * h_avg)
    
    def grad_kernel_pair(self, xi, xj, hi, hj):
        """Gradient terms for particle pair interaction."""
        dx = xi - xj
        grad_i = self.gradW(dx, 2.0 * hi)
        grad_j = self.gradW(-dx, 2.0 * hj)
        return grad_i, grad_j
