"""
Relativistic Riemann Problem Solver

This program computes the solution of a 1D relativistic Riemann problem 
with initial data UL if X<0.5 and UR if X>0.5 in the spatial domain [0, 1].

Converted from Fortran 77 code (41114_2016_3_MOESM6_ESM.f)
Based on Marti and Mueller, J. Fluid Mech., (1994)
"""

import numpy as np
import sys


class RelativisiticRiemannSolver:
    """Solver for relativistic Riemann problems in hydrodynamics"""
    
    def __init__(self, gamma):
        """
        Initialize the solver with adiabatic index.
        
        Args:
            gamma: Adiabatic index of the gas
        """
        self.gamma = gamma
        
        # Left state
        self.rhol = None
        self.pl = None
        self.ul = None
        self.hl = None
        self.csl = None
        self.vell = None
        self.wl = None
        
        # Right state
        self.rhor = None
        self.pr = None
        self.ur = None
        self.hr = None
        self.csr = None
        self.velr = None
        self.wr = None
        
        # Left star state
        self.rhols = None
        self.uls = None
        self.hls = None
        self.csls = None
        self.vells = None
        self.vshockl = None
        
        # Right star state
        self.rhors = None
        self.urs = None
        self.hrs = None
        self.csrs = None
        self.velrs = None
        self.vshockr = None
        
    def set_initial_states(self, pl, rhol, vell, pr, rhor, velr):
        """
        Set the initial left and right states.
        
        Args:
            pl: Left pressure
            rhol: Left density
            vell: Left flow velocity
            pr: Right pressure
            rhor: Right density
            velr: Right flow velocity
        """
        # Left state
        self.pl = pl
        self.rhol = rhol
        self.vell = vell
        
        # Right state
        self.pr = pr
        self.rhor = rhor
        self.velr = velr
        
        # Compute derived quantities for left state
        self.ul = self.pl / (self.gamma - 1.0) / self.rhol
        self.hl = 1.0 + self.ul + self.pl / self.rhol
        self.csl = np.sqrt(self.gamma * self.pl / self.rhol / self.hl)
        self.wl = 1.0 / np.sqrt(1.0 - self.vell**2)
        
        # Compute derived quantities for right state
        self.ur = self.pr / (self.gamma - 1.0) / self.rhor
        self.hr = 1.0 + self.ur + self.pr / self.rhor
        self.csr = np.sqrt(self.gamma * self.pr / self.rhor / self.hr)
        self.wr = 1.0 / np.sqrt(1.0 - self.velr**2)
        
    def get_velocity(self, p, rhoa, pa, ua, ha, csa, vela, wa, direction):
        """
        Compute flow velocity behind a rarefaction or shock.
        
        Args:
            p: Post-wave pressure
            rhoa: Density ahead of wave
            pa: Pressure ahead of wave
            ua: Specific internal energy ahead of wave
            ha: Specific enthalpy ahead of wave
            csa: Local sound speed ahead of wave
            vela: Flow velocity ahead of wave
            wa: Flow Lorentz factor ahead of wave
            direction: 'L' for left or 'R' for right propagating wave
            
        Returns:
            Tuple of (rho, u, h, cs, vel, vshock) in post-wave state
        """
        sign = -1.0 if direction == 'L' else 1.0
        
        if p > pa:
            # Shock wave
            a = 1.0 + (self.gamma - 1.0) * (pa - p) / self.gamma / p
            b = 1.0 - a
            c = ha * (pa - p) / rhoa - ha**2
            
            # Check for unphysical enthalpies
            if c > (b**2 / 4.0 / a):
                raise ValueError("Unphysical specific enthalpy in intermediate state")
            
            # Specific enthalpy in post-wave state
            h = (-b + np.sqrt(b**2 - 4.0 * a * c)) / 2.0 / a
            
            # Density in post-wave state
            rho = self.gamma * p / (self.gamma - 1.0) / (h - 1.0)
            
            # Specific internal energy in post-wave state
            u = p / (self.gamma - 1.0) / rho
            
            # Mass flux across the wave
            j = sign * np.sqrt((p - pa) / (ha / rhoa - h / rho))
            
            # Shock velocity
            a = j**2 + (rhoa * wa)**2
            b = -vela * rhoa**2 * wa**2
            vshock = (-b + sign * j**2 * np.sqrt(1.0 + rhoa**2 / j**2)) / a
            wshock = 1.0 / np.sqrt(1.0 - vshock**2)
            
            # Flow velocity in post-shock state
            a = wshock * (p - pa) / j + ha * wa * vela
            b = ha * wa + (p - pa) * (wshock * vela / j + 1.0 / rhoa / wa)
            vel = a / b
            
            # Local sound speed in post-shock state
            cs = np.sqrt(self.gamma * p / rho / h)
            
        else:
            # Rarefaction wave
            # Polytropic constant across rarefaction
            k = pa / rhoa**self.gamma
            
            # Density behind rarefaction
            rho = (p / k)**(1.0 / self.gamma)
            
            # Specific internal energy behind rarefaction
            u = p / (self.gamma - 1.0) / rho
            
            # Specific enthalpy behind rarefaction
            h = 1.0 + u + p / rho
            
            # Local sound speed behind rarefaction
            cs = np.sqrt(self.gamma * p / (rho + self.gamma * p / (self.gamma - 1.0)))
            
            # Flow velocity behind rarefaction
            sqgl1 = np.sqrt(self.gamma - 1.0)
            a = ((1.0 + vela) / (1.0 - vela) *
                 ((sqgl1 + csa) / (sqgl1 - csa) *
                  (sqgl1 - cs) / (sqgl1 + cs))**(-sign * 2.0 / sqgl1))
            
            vel = (a - 1.0) / (a + 1.0)
            vshock = 0.0
            
        return rho, u, h, cs, vel, vshock
        
    def get_dvel(self, p):
        """
        Compute difference in flow speed between left and right intermediate states.
        
        Args:
            p: Pressure in intermediate state
            
        Returns:
            Difference in flow velocity (vells - velrs)
        """
        # Left wave
        (self.rhols, self.uls, self.hls, self.csls, 
         self.vells, self.vshockl) = self.get_velocity(
            p, self.rhol, self.pl, self.ul, self.hl, 
            self.csl, self.vell, self.wl, 'L')
        
        # Right wave
        (self.rhors, self.urs, self.hrs, self.csrs, 
         self.velrs, self.vshockr) = self.get_velocity(
            p, self.rhor, self.pr, self.ur, self.hr, 
            self.csr, self.velr, self.wr, 'R')
        
        return self.vells - self.velrs
        
    def get_pressure(self, pmin, pmax, tol=0.0):
        """
        Find pressure in intermediate state using bisection method.
        
        Args:
            pmin: Left endpoint of interval
            pmax: Right endpoint of interval
            tol: Tolerance for convergence (>= 0)
            
        Returns:
            Pressure in intermediate state
        """
        # Use simple bisection for robustness
        a = pmin
        b = pmax
        fa = self.get_dvel(a)
        fb = self.get_dvel(b)
        
        # Ensure fa and fb have opposite signs
        if fa * fb > 0:
            # Try to find a valid bracket
            for _ in range(50):
                mid = 0.5 * (a + b)
                fmid = self.get_dvel(mid)
                if fa * fmid < 0:
                    b = mid
                    fb = fmid
                    break
                elif fmid * fb < 0:
                    a = mid
                    fa = fmid
                    break
                a *= 0.5
                b *= 2.0
                fa = self.get_dvel(a)
                fb = self.get_dvel(b)
        
        # Bisection iteration
        max_iter = 100
        eps = 1e-12
        
        for _ in range(max_iter):
            c = 0.5 * (a + b)
            fc = self.get_dvel(c)
            
            if abs(fc) < eps or abs(b - a) < eps * abs(c):
                return c
            
            if fa * fc < 0:
                b = c
                fb = fc
            else:
                a = c
                fa = fc
        
        return 0.5 * (a + b)
                
    def rarefaction(self, xi, rhoa, pa, ua, csa, vela, direction):
        """
        Compute flow state within a rarefaction.
        
        Args:
            xi: Similarity variable
            rhoa: Density ahead of rarefaction
            pa: Pressure ahead of rarefaction
            ua: Specific internal energy ahead of rarefaction
            csa: Local sound speed ahead of rarefaction
            vela: Flow velocity ahead of rarefaction
            direction: 'L' for left or 'R' for right propagating wave
            
        Returns:
            Tuple of (rho, p, u, vel) at point within rarefaction
        """
        sign = 1.0 if direction == 'L' else -1.0
        
        b = np.sqrt(self.gamma - 1.0)
        c = (b + csa) / (b - csa)
        d = -sign * b / 2.0
        k = (1.0 + xi) / (1.0 - xi)
        l = c * k**d
        v = ((1.0 - vela) / (1.0 + vela))**d
        
        # Newton-Raphson iteration for sound speed
        ocs2 = csa
        max_iter = 100
        for _ in range(max_iter):
            fcs2 = (l * v * (1.0 + sign * ocs2)**d * (ocs2 - b) +
                    (1.0 - sign * ocs2)**d * (ocs2 + b))
            
            dfdcs2 = (l * v * (1.0 + sign * ocs2)**d *
                      (1.0 + sign * d * (ocs2 - b) / (1.0 + sign * ocs2)) +
                      (1.0 - sign * ocs2)**d *
                      (1.0 - sign * d * (ocs2 + b) / (1.0 - sign * ocs2)))
            
            cs2 = ocs2 - fcs2 / dfdcs2
            
            if abs(cs2 - ocs2) / ocs2 <= 5e-7:
                break
            ocs2 = cs2
        
        vel = (xi + sign * cs2) / (1.0 + sign * xi * cs2)
        
        rho = rhoa * ((cs2**2 * (self.gamma - 1.0 - csa**2)) /
                      (csa**2 * (self.gamma - 1.0 - cs2**2)))**(1.0 / (self.gamma - 1.0))
        
        p = cs2**2 * (self.gamma - 1.0) * rho / (self.gamma - 1.0 - cs2**2) / self.gamma
        
        u = p / (self.gamma - 1.0) / rho
        
        return rho, p, u, vel
        
    def solve(self, t, n=400, tol=0.0):
        """
        Solve the Riemann problem at time t.
        
        Args:
            t: Time for the solution
            n: Number of spatial points
            tol: Tolerance for pressure solver
            
        Returns:
            Tuple of arrays (x, pressure, density, velocity, internal_energy)
        """
        # Find pressure bounds
        pmin = (self.pl + self.pr) / 2.0
        pmax = pmin
        
        iloop = 0
        max_loops = 100
        
        for _ in range(max_loops):
            iloop += 1
            pmin = 0.5 * max(pmin, 0.0)
            pmax = 2.0 * pmax
            
            dvel1 = self.get_dvel(pmin)
            dvel2 = self.get_dvel(pmax)
            
            if dvel1 * dvel2 <= 0.0:
                break
        
        # Find pressure in intermediate state
        ps = self.get_pressure(pmin, pmax, tol)
        vels = 0.5 * (self.vells + self.velrs)
        
        # Compute wave positions
        if self.pl >= ps:
            # Left rarefaction
            x1 = 0.5 + (self.vell - self.csl) / (1.0 - self.vell * self.csl) * t
            x2 = 0.5 + (vels - self.csls) / (1.0 - vels * self.csls) * t
        else:
            # Left shock
            x1 = 0.5 + self.vshockl * t
            x2 = x1
        
        x3 = 0.5 + vels * t
        
        if self.pr >= ps:
            # Right rarefaction
            x4 = 0.5 + (vels + self.csrs) / (1.0 + vels * self.csrs) * t
            x5 = 0.5 + (self.velr + self.csr) / (1.0 + self.velr * self.csr) * t
        else:
            # Right shock
            x4 = 0.5 + self.vshockr * t
            x5 = x4
        
        # Compute solution on mesh
        x = np.linspace(0, 1, n)
        pressure = np.zeros(n)
        density = np.zeros(n)
        velocity = np.zeros(n)
        internal_energy = np.zeros(n)
        
        for i in range(n):
            if x[i] <= x1:
                # Left state
                pressure[i] = self.pl
                density[i] = self.rhol
                velocity[i] = self.vell
                internal_energy[i] = self.ul
                
            elif x[i] <= x2:
                # Left rarefaction fan
                xi = (x[i] - 0.5) / t
                density[i], pressure[i], internal_energy[i], velocity[i] = \
                    self.rarefaction(xi, self.rhol, self.pl, self.ul, 
                                   self.csl, self.vell, 'L')
                
            elif x[i] <= x3:
                # Left star state
                pressure[i] = ps
                density[i] = self.rhols
                velocity[i] = vels
                internal_energy[i] = self.uls
                
            elif x[i] <= x4:
                # Right star state
                pressure[i] = ps
                density[i] = self.rhors
                velocity[i] = vels
                internal_energy[i] = self.urs
                
            elif x[i] <= x5:
                # Right rarefaction fan
                xi = (x[i] - 0.5) / t
                density[i], pressure[i], internal_energy[i], velocity[i] = \
                    self.rarefaction(xi, self.rhor, self.pr, self.ur, 
                                   self.csr, self.velr, 'R')
                
            else:
                # Right state
                pressure[i] = self.pr
                density[i] = self.rhor
                velocity[i] = self.velr
                internal_energy[i] = self.ur
        
        return x, pressure, density, velocity, internal_energy
        
    def save_solution(self, filename, x, pressure, density, velocity, internal_energy, t):
        """
        Save solution to file.
        
        Args:
            filename: Output filename
            x: Spatial positions
            pressure: Pressure array
            density: Density array
            velocity: Velocity array
            internal_energy: Internal energy array
            t: Time of solution
        """
        n = len(x)
        with open(filename, 'w') as f:
            f.write(f"{n:5d} {t:10.5f}\n")
            for i in range(n):
                f.write(f"{x[i]:15.8e} {pressure[i]:15.8e} {density[i]:15.8e} "
                       f"{velocity[i]:15.8e} {internal_energy[i]:15.8e}\n")


def main():
    """Main program - interactive mode"""
    print("Relativistic Riemann Problem Solver")
    print("=" * 50)
    
    # Get input parameters
    gamma = float(input("Adiabatic index of the gas: "))
    t = float(input("Time for the solution: "))
    
    print("\n-- LEFT STATE --")
    pl = float(input("  Pressure: "))
    rhol = float(input("  Density: "))
    vell = float(input("  Flow velocity: "))
    
    print("\n-- RIGHT STATE --")
    pr = float(input("  Pressure: "))
    rhor = float(input("  Density: "))
    velr = float(input("  Flow velocity: "))
    
    # Create solver and solve
    solver = RelativisiticRiemannSolver(gamma)
    solver.set_initial_states(pl, rhol, vell, pr, rhor, velr)
    
    print("\nSolving...")
    x, pressure, density, velocity, internal_energy = solver.solve(t, n=400)
    
    # Save solution
    solver.save_solution('solution.dat', x, pressure, density, velocity, internal_energy, t)
    
    print("Solution saved to solution.dat")
    

if __name__ == "__main__":
    main()
