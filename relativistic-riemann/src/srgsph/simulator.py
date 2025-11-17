"""
SRGSPH 1D Simulator

Main simulation class implementing Special Relativistic Godunov SPH
based on the Kitajima formulation (arXiv:2510.18251v1).

Key features:
- Volume-based density calculation (Eqs. 33-37)
- Riemann solver for particle interactions
- MUSCL reconstruction for second-order accuracy
- Euler time integration (Eqs. 70-72)
- CFL timestep criterion (Eq. 73)
"""

import numpy as np
from .particle import Particle
from .kernel import GaussianKernel
from .density import VolumeCalculator, GradientCalculator
import sys
import os

# Add parent directory to path to import kitajima_solver
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from kitajima_solver import KitajimaRiemannSolver


class SRGSPH1D:
    """
    1D Special Relativistic Godunov SPH Simulator.
    
    Implements the SRGSPH method with:
    - Variable smoothing length
    - Volume-based particle definition
    - Riemann solver for shock capturing
    - MUSCL reconstruction
    """
    
    def __init__(self, gamma_c=5.0/3.0, c=1.0, C_CFL=0.3, 
                 eta=1.0, C_smooth=2.0, C_shock=3.0, C_cd=1.0,
                 use_variable_h=True, use_muscl=True):
        """
        Initialize simulator.
        
        Args:
            gamma_c: Adiabatic index (default: 5/3)
            c: Speed of light (default: 1.0)
            C_CFL: CFL constant (default: 0.3)
            eta: Smoothing length coefficient (default: 1.0)
            C_smooth: Smoothing coefficient for h calculation (default: 2.0)
            C_shock: Shock detection constant (default: 3.0)
            C_cd: Contact discontinuity detection constant (default: 1.0)
            use_variable_h: Use variable smoothing length (default: True)
            use_muscl: Use MUSCL reconstruction (default: True)
        """
        self.gamma_c = gamma_c
        self.c = c
        self.C_CFL = C_CFL
        self.use_variable_h = use_variable_h
        self.use_muscl = use_muscl
        
        # Initialize components
        self.kernel = GaussianKernel(dim=1)
        self.volume_calc = VolumeCalculator(self.kernel, eta=eta, C_smooth=C_smooth, dim=1)
        self.gradient_calc = GradientCalculator(C_shock=C_shock, C_cd=C_cd)
        
        # Riemann solver
        self.riemann_solver = KitajimaRiemannSolver(gamma_c, c)
        
        # Simulation state
        self.particles = []
        self.time = 0.0
        self.step = 0
        
    def setup_particles(self, positions, velocities, densities, pressures, baryon_numbers=None):
        """
        Initialize particles with primitive variables.
        
        Args:
            positions: Array of particle positions
            velocities: Array of particle velocities
            densities: Array of rest frame baryon densities
            pressures: Array of pressures
            baryon_numbers: Array of baryon numbers (default: uniform)
            
        Returns:
            None
        """
        n_particles = len(positions)
        
        if baryon_numbers is None:
            # Uniform baryon number
            baryon_numbers = np.ones(n_particles)
        
        self.particles = []
        for i in range(n_particles):
            p = Particle(
                nu=baryon_numbers[i],
                x=positions[i],
                v=velocities[i],
                n=densities[i],
                P=pressures[i],
                gamma_c=self.gamma_c,
                c=self.c
            )
            self.particles.append(p)
        
        # Initialize volumes and smoothing lengths
        if self.use_variable_h:
            self.volume_calc.compute_volumes_and_smoothing_lengths(self.particles)
        else:
            # Fixed smoothing length
            avg_spacing = (positions[-1] - positions[0]) / (n_particles - 1)
            for p in self.particles:
                p.h = avg_spacing
                p.Vp = p.h  # For 1D
        
        # Find neighbors
        self.volume_calc.find_neighbors(self.particles)
        
    def compute_timestep(self):
        """
        Compute timestep from CFL condition (Eq. 73).
        
        Δt = C_CFL * min_i[h_i / cs,i]
        
        Returns:
            Timestep
        """
        dt_min = np.inf
        
        for p in self.particles:
            if p.cs > 1e-15:
                dt_i = p.h / p.cs
                dt_min = min(dt_min, dt_i)
        
        return self.C_CFL * dt_min
    
    def compute_forces(self):
        """
        Compute SPH forces using Riemann solver.
        
        Implements Eqs. 64-65:
        ⟨ν_i Ṡ_i⟩ = -Σj P*_ij V²ij,interp [∇i W(xi-xj, 2hi) - ∇j W(xi-xj, 2hj)]
        ⟨ν_i ė_i⟩ = -Σj P*_ij v*_ij · V²ij,interp [∇i W(xi-xj, 2hi) - ∇j W(xi-xj, 2hj)]
        
        Returns:
            (dS_dt, de_dt): Arrays of time derivatives
        """
        n_particles = len(self.particles)
        dS_dt = np.zeros(n_particles)
        de_dt = np.zeros(n_particles)
        
        # Compute gradients for MUSCL
        if self.use_muscl:
            self.gradient_calc.compute_gradients(self.particles, self.kernel)
        
        # Loop over particle pairs
        for i, p_i in enumerate(self.particles):
            for j in p_i.neighbors:
                if j <= i:
                    continue  # Avoid double counting
                
                p_j = self.particles[j]
                
                # Check if gradients should be used (Eq. 66)
                if self.use_muscl:
                    use_grad_i, use_grad_j = self.gradient_calc.apply_limiters(
                        self.particles, i, j)
                else:
                    use_grad_i, use_grad_j = False, False
                
                # Get states at interface (with MUSCL reconstruction)
                dx = p_j.x - p_i.x
                s = 0.5 * dx  # Distance to midpoint
                
                if use_grad_i:
                    n_L = p_i.n + p_i.grad_n * s
                    P_L = p_i.P + p_i.grad_P * s
                    v_L = p_i.v + p_i.grad_v * s
                else:
                    n_L = p_i.n
                    P_L = p_i.P
                    v_L = p_i.v
                
                if use_grad_j:
                    n_R = p_j.n - p_j.grad_n * s
                    P_R = p_j.P - p_j.grad_P * s
                    v_R = p_j.v - p_j.grad_v * s
                else:
                    n_R = p_j.n
                    P_R = p_j.P
                    v_R = p_j.v
                
                # Ensure positive pressure and density
                n_L = max(n_L, 1e-10)
                n_R = max(n_R, 1e-10)
                P_L = max(P_L, 1e-10)
                P_R = max(P_R, 1e-10)
                
                # Enforce causality: |v| < c
                v_max = 0.99 * self.c
                v_L = np.clip(v_L, -v_max, v_max)
                v_R = np.clip(v_R, -v_max, v_max)
                
                # If negative values, fall back to first-order
                if n_L <= 0 or n_R <= 0 or P_L <= 0 or P_R <= 0:
                    n_L, P_L, v_L = p_i.n, p_i.P, p_i.v
                    n_R, P_R, v_R = p_j.n, p_j.P, p_j.v
                
                # Solve Riemann problem
                self.riemann_solver.set_initial_states(P_L, n_L, v_L, P_R, n_R, v_R)
                
                try:
                    # Get star state
                    pmin = 0.5 * min(P_L, P_R)
                    pmax = 2.0 * max(P_L, P_R)
                    P_star = self.riemann_solver.get_pressure(pmin, pmax)
                    v_star = 0.5 * (self.riemann_solver.vls + self.riemann_solver.vrs)
                except:
                    # Fallback if Riemann solver fails
                    P_star = 0.5 * (P_L + P_R)
                    v_star = 0.5 * (v_L + v_R)
                
                # Compute V²ij,interp (Eq. 63)
                V2_ij = self.volume_calc.compute_Vij_squared(self.particles, i, j)
                
                # Compute kernel gradients
                grad_i, grad_j = self.kernel.grad_kernel_pair(p_i.x, p_j.x, p_i.h, p_j.h)
                
                # Force term: [∇i W - ∇j W]
                grad_diff = grad_i - grad_j
                
                # Update forces (anti-symmetric)
                # dS/dt = -Σ P* V²ij [∇i W - ∇j W]
                force_momentum = -P_star * V2_ij * grad_diff
                dS_dt[i] += force_momentum
                dS_dt[j] -= force_momentum  # Newton's third law
                
                # de/dt = -Σ P* v* · V²ij [∇i W - ∇j W]
                force_energy = -P_star * v_star * V2_ij * grad_diff
                de_dt[i] += force_energy
                de_dt[j] -= force_energy
        
        # Divide by baryon number (Eqs. 64-65 have ⟨ν_i Ṡ_i⟩ on left side)
        for i, p in enumerate(self.particles):
            if p.nu > 1e-15:
                dS_dt[i] /= p.nu
                de_dt[i] /= p.nu
        
        return dS_dt, de_dt
    
    def step_forward(self, dt):
        """
        Advance simulation by one timestep using Euler method (Eqs. 70-72).
        
        Args:
            dt: Timestep
            
        Returns:
            None
        """
        # Compute forces
        dS_dt, de_dt = self.compute_forces()
        
        # Update conserved variables (Eqs. 70-71)
        for i, p in enumerate(self.particles):
            p.S += dS_dt[i] * dt
            p.e += de_dt[i] * dt
        
        # Recover primitive variables
        for p in self.particles:
            p.recover_primitives()
        
        # Update positions (Eq. 72)
        for p in self.particles:
            p.x += p.v * dt
        
        # Update volumes and smoothing lengths
        if self.use_variable_h:
            self.volume_calc.compute_volumes_and_smoothing_lengths(
                self.particles, max_iter=5)
        
        # Update neighbors
        self.volume_calc.find_neighbors(self.particles)
        
        # Update time
        self.time += dt
        self.step += 1
    
    def run(self, t_end, dt=None, output_freq=10):
        """
        Run simulation until t_end.
        
        Args:
            t_end: End time
            dt: Fixed timestep (if None, use CFL condition)
            output_freq: Output frequency (print every N steps)
            
        Returns:
            None
        """
        print(f"Starting SRGSPH simulation")
        print(f"  Particles: {len(self.particles)}")
        print(f"  Variable h: {self.use_variable_h}")
        print(f"  MUSCL: {self.use_muscl}")
        print(f"  Target time: {t_end}")
        print()
        
        while self.time < t_end:
            # Compute timestep
            if dt is None:
                dt_cfl = self.compute_timestep()
                # Don't overshoot end time
                dt_step = min(dt_cfl, t_end - self.time)
            else:
                dt_step = min(dt, t_end - self.time)
            
            # Take step
            self.step_forward(dt_step)
            
            # Output
            if self.step % output_freq == 0:
                print(f"Step {self.step:6d}: t = {self.time:.6f}, dt = {dt_step:.6e}")
        
        print(f"\nSimulation complete at t = {self.time:.6f}")
    
    def get_solution(self):
        """
        Extract current solution.
        
        Returns:
            Dictionary with arrays of particle properties
        """
        return {
            'x': np.array([p.x for p in self.particles]),
            'v': np.array([p.v for p in self.particles]),
            'n': np.array([p.n for p in self.particles]),
            'N': np.array([p.N for p in self.particles]),
            'P': np.array([p.P for p in self.particles]),
            'u': np.array([p.u for p in self.particles]),
            'gamma': np.array([p.gamma for p in self.particles]),
            'S': np.array([p.S for p in self.particles]),
            'e': np.array([p.e for p in self.particles]),
            'h': np.array([p.h for p in self.particles]),
            'cs': np.array([p.cs for p in self.particles]),
            'time': self.time,
            'step': self.step
        }
    
    def save_solution(self, filename):
        """
        Save current solution to file.
        
        Args:
            filename: Output filename
        """
        sol = self.get_solution()
        
        header = (f"SRGSPH Solution at t={sol['time']:.6f}, step={sol['step']}\n"
                  f"Columns: x, v, n, N, P, u, gamma, S, e, h, cs")
        
        data = np.column_stack([
            sol['x'], sol['v'], sol['n'], sol['N'], sol['P'],
            sol['u'], sol['gamma'], sol['S'], sol['e'], sol['h'], sol['cs']
        ])
        
        np.savetxt(filename, data, header=header, comments='# ')
        print(f"Solution saved to {filename}")
