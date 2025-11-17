"""
Volume-based Density Calculation for SRGSPH

Implements the volume-based approach to particle volume and density calculation
as described in Section 2.3-2.4 of the paper.

Key equations:
- Vp(x) = [Σj W(x - xj, h(x))]⁻¹  (Eq. 33)
- h(x) = η Vp*(x)^(1/d)  (Eq. 35)
- Vp*(x) = [Σj W(x - xj, C_smooth h(x))]⁻¹  (Eq. 36)
- N(x) = ν(x) / Vp(x)  (Eq. 37)
"""

import numpy as np


class VolumeCalculator:
    """
    Calculate particle volumes and densities using volume-based approach.
    
    This avoids discontinuities in smoothing length at contact discontinuities
    when particles have different baryon numbers (see Section 2.4 and Fig. 2-3).
    """
    
    def __init__(self, kernel, eta=1.0, C_smooth=2.0, dim=1):
        """
        Initialize volume calculator.
        
        Args:
            kernel: Kernel object (GaussianKernel or WendlandC2Kernel)
            eta: Smoothing length coefficient (default: 1.0)
            C_smooth: Smoothing coefficient > 1 to ensure smooth h(x) (default: 2.0)
            dim: Number of dimensions
        """
        self.kernel = kernel
        self.eta = eta
        self.C_smooth = C_smooth
        self.dim = dim
    
    def compute_volumes_and_smoothing_lengths(self, particles, max_iter=10, tol=1e-6):
        """
        Compute particle volumes and smoothing lengths iteratively.
        
        For first timestep: iterate Eqs. 35-36 until convergence
        For subsequent timesteps: use previous h as initial guess
        
        Args:
            particles: List of Particle objects
            max_iter: Maximum iterations for convergence
            tol: Convergence tolerance
            
        Returns:
            None (modifies particles in place)
        """
        n_particles = len(particles)
        positions = np.array([p.x for p in particles])
        
        # Initialize or get previous smoothing lengths
        h = np.array([p.h if p.h > 0 else self._estimate_initial_h(particles, i) 
                      for i, p in enumerate(particles)])
        
        # Iterate to find consistent h and Vp*
        for iteration in range(max_iter):
            h_old = h.copy()
            
            # Compute Vp* for each particle using C_smooth * h
            # Optimized: vectorize distance calculations
            Vp_star = np.zeros(n_particles)
            for i in range(n_particles):
                dx = positions[i] - positions
                r = np.abs(dx)
                # Only consider particles within support
                mask = r < 3.0 * self.C_smooth * h[i]
                sum_W = np.sum(self.kernel.W(r[mask], self.C_smooth * h[i]))
                
                if sum_W > 1e-15:
                    Vp_star[i] = 1.0 / sum_W
                else:
                    Vp_star[i] = h[i]**self.dim  # Fallback
            
            # Update h from Eq. 35: h = η Vp*^(1/d)
            h = self.eta * Vp_star**(1.0 / self.dim)
            
            # Check convergence
            max_change = np.max(np.abs(h - h_old) / np.maximum(h_old, 1e-10))
            if max_change < tol:
                break
        
        # Compute actual Vp using converged h (Eq. 33)
        # Optimized: vectorize distance calculations
        Vp = np.zeros(n_particles)
        for i in range(n_particles):
            dx = positions[i] - positions
            r = np.abs(dx)
            # Only consider particles within support
            mask = r < 3.0 * h[i]
            sum_W = np.sum(self.kernel.W(r[mask], h[i]))
            
            if sum_W > 1e-15:
                Vp[i] = 1.0 / sum_W
            else:
                Vp[i] = h[i]**self.dim
        
        # Update particles
        for i, p in enumerate(particles):
            p.h = h[i]
            p.Vp = Vp[i]
            # Update lab frame density: N = ν / Vp  (Eq. 37)
            p.N = p.nu / Vp[i]
            # Update rest frame density: n = N / γ
            p.n = p.N / p.gamma
    
    def _estimate_initial_h(self, particles, i):
        """
        Estimate initial smoothing length from particle spacing.
        
        Args:
            particles: List of particles
            i: Index of particle
            
        Returns:
            Estimated smoothing length
        """
        # Find average distance to nearest neighbors
        x_i = particles[i].x
        distances = []
        
        for j, p in enumerate(particles):
            if j != i:
                distances.append(abs(p.x - x_i))
        
        if len(distances) == 0:
            return 1.0
        
        distances.sort()
        # Use average of 3 nearest neighbors (or fewer if not enough particles)
        n_neighbors = min(3, len(distances))
        avg_dist = np.mean(distances[:n_neighbors])
        
        # h should be roughly the particle spacing
        return self.eta * avg_dist
    
    def find_neighbors(self, particles, kernel_radius_factor=3.0):
        """
        Find neighbors for each particle.
        
        Neighbors are particles within kernel_radius_factor * h.
        
        Args:
            particles: List of Particle objects
            kernel_radius_factor: Multiplier for smoothing length
            
        Returns:
            None (modifies particles.neighbors in place)
        """
        n_particles = len(particles)
        
        for i, p_i in enumerate(particles):
            neighbors = []
            search_radius = kernel_radius_factor * p_i.h
            
            for j, p_j in enumerate(particles):
                if i == j:
                    continue
                
                dx = abs(p_i.x - p_j.x)
                # Use max of both smoothing lengths for neighbor search
                max_h = max(p_i.h, p_j.h)
                
                if dx <= kernel_radius_factor * max_h:
                    neighbors.append(j)
            
            p_i.neighbors = neighbors
    
    def compute_Vij_squared(self, particles, i, j):
        """
        Compute V²ij for particle pair (Eq. 63).
        
        V²ij = (V²ij(hi) + V²ij(hj)) / 2
        
        This is used in the force calculation as V²ij,interp.
        
        Args:
            particles: List of particles
            i, j: Particle indices
            
        Returns:
            V²ij value
        """
        p_i = particles[i]
        p_j = particles[j]
        
        # For simplicity, we approximate V²ij ≈ Vp,i * Vp,j
        # This is consistent with the volume-based approach
        # More sophisticated interpolation could be implemented
        
        return np.sqrt(p_i.Vp * p_j.Vp)


class GradientCalculator:
    """
    Calculate gradients for MUSCL reconstruction.
    
    Implements monotonicity-preserving gradient calculation with
    shock detection from Eq. 66.
    """
    
    def __init__(self, C_shock=3.0, C_cd=1.0):
        """
        Initialize gradient calculator.
        
        Args:
            C_shock: Shock detection constant (default: 3.0)
            C_cd: Contact discontinuity detection constant (default: 1.0)
        """
        self.C_shock = C_shock
        self.C_cd = C_cd
    
    def compute_gradients(self, particles, kernel):
        """
        Compute gradients of primitive variables for each particle.
        
        Uses simple SPH gradient formula:
        ∇f_i ≈ Σj (f_j - f_i) ∇W_ij
        
        Args:
            particles: List of Particle objects
            kernel: Kernel object
            
        Returns:
            None (modifies particle gradients in place)
        """
        for i, p_i in enumerate(particles):
            grad_n = 0.0
            grad_P = 0.0
            grad_v = 0.0
            sum_W = 0.0
            
            for j in p_i.neighbors:
                p_j = particles[j]
                
                dx = p_i.x - p_j.x
                h_avg = 0.5 * (p_i.h + p_j.h)
                grad_W = kernel.gradW(dx, h_avg)
                
                grad_n += (p_j.n - p_i.n) * grad_W
                grad_P += (p_j.P - p_i.P) * grad_W
                grad_v += (p_j.v - p_i.v) * grad_W
                sum_W += abs(grad_W)
            
            # Normalize
            if sum_W > 1e-15:
                p_i.grad_n = grad_n / sum_W
                p_i.grad_P = grad_P / sum_W
                p_i.grad_v = grad_v / sum_W
            else:
                p_i.grad_n = 0.0
                p_i.grad_P = 0.0
                p_i.grad_v = 0.0
    
    def apply_limiters(self, particles, i, j):
        """
        Apply monotonicity constraints for MUSCL reconstruction.
        
        From Eq. 66: Set gradients to zero if:
        1. C_shock * e_ij · (v_i - v_j) > min(cs,i, cs,j)  (shock)
        2. |log10(P_i / P_j)| > C_cd  (contact discontinuity)
        
        Args:
            particles: List of particles
            i, j: Particle indices
            
        Returns:
            (use_gradient_i, use_gradient_j): Booleans indicating if gradients should be used
        """
        p_i = particles[i]
        p_j = particles[j]
        
        # Unit vector from i to j
        dx = p_j.x - p_i.x
        if abs(dx) < 1e-15:
            return True, True
        
        e_ij = np.sign(dx)
        
        # Shock detection
        dv = p_i.v - p_j.v
        cs_min = min(p_i.cs, p_j.cs)
        
        is_shock = self.C_shock * e_ij * dv > cs_min
        
        # Contact discontinuity detection
        if abs(p_j.P) > 1e-15:
            pressure_ratio = abs(np.log10(p_i.P / p_j.P))
            is_cd = pressure_ratio > self.C_cd
        else:
            is_cd = True
        
        # If shock or CD detected, disable gradients
        use_gradient = not (is_shock or is_cd)
        
        return use_gradient, use_gradient
