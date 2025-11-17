"""
Core solver module for relativistic Riemann problems.

This module implements the exact Riemann solver for special relativistic 
hydrodynamics following Marti and Mueller (1994).
"""

import numpy as np


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
        Find pressure in intermediate state using Brent's method.
        
        Args:
            pmin: Left endpoint of interval
            pmax: Right endpoint of interval
            tol: Tolerance for convergence (>= 0)
            
        Returns:
            Pressure in intermediate state
        """
        # Compute machine precision
        eps = 1.0
        while (1.0 + eps / 2.0) > 1.0:
            eps = eps / 2.0
        
        # Initialization
        a = pmin
        b = pmax
        fa = self.get_dvel(a)
        fb = self.get_dvel(b)
        
        # Begin iteration
        c = a
        fc = fa
        d = b - a
        e = d
        
        while True:
            if abs(fc) >= abs(fb):
                a, b, c = b, c, a
                fa, fb, fc = fb, fc, fa
            
            # Convergence test
            tol1 = 2.0 * eps * abs(b) + 0.5 * tol
            xm = 0.5 * (c - b)
            
            if abs(xm) <= tol1 or fb == 0.0:
                return b
            
            # Is bisection necessary?
            if abs(e) < tol1 or abs(fa) <= abs(fb):
                # Bisection
                d = xm
                e = d
            else:
                # Is quadratic interpolation possible?
                if a != c:
                    # Inverse quadratic interpolation
                    q = fa / fc
                    r = fb / fc
                    s = fb / fa
                    p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0))
                    q = (q - 1.0) * (r - 1.0) * (s - 1.0)
                else:
                    # Linear interpolation
                    s = fb / fa
                    p = 2.0 * xm * s
                    q = 1.0 - s
                
                # Adjust signs
                if p > 0.0:
                    q = -q
                p = abs(p)
                
                # Is interpolation acceptable?
                if (2.0 * p) >= (3.0 * xm * q - abs(tol1 * q)) or p >= abs(0.5 * e * q):
                    d = xm
                    e = d
                else:
                    e = d
                    d = p / q
            
            # Complete step
            a = b
            fa = fb
            if abs(d) > tol1:
                b = b + d
            else:
                b = b + np.sign(xm) * tol1
            
            fb = self.get_dvel(b)
            
            if (fb * (fc / abs(fc))) > 0.0:
                c = a
                fc = fa
                d = b - a
                e = d
                
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
        
        max_loops = 100
        
        for _ in range(max_loops):
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
