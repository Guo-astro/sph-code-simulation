"""
Kitajima Formulation: Special Relativistic Riemann Solver

This module implements the relativistic Riemann solver using the Kitajima formulation
(Kitajima et al., 2025, arXiv:2510.18251v1) which uses:
- Baryon number density N (not mass density)
- Explicit light speed c (not c=1)
- Canonical momentum per baryon S
- Canonical energy per baryon e

Based on equations from Section 2.1 of the paper.
"""

import numpy as np


class KitajimaRiemannSolver:
    """
    Relativistic Riemann solver using Kitajima formulation with baryon number density.
    
    Physical variables:
    - N: Baryon number density in lab frame
    - n: Baryon number density in rest frame  
    - S: Relativistic canonical momentum per baryon = γHv
    - e: Canonical energy per baryon = γH - P/(Nc²)
    - H: Enthalpy per baryon = 1 + u/c² + P/(nc²)
    - γ: Lorentz factor = 1/sqrt(1 - v²/c²)
    - c: Speed of light (explicit, not set to 1)
    """
    
    def __init__(self, gamma_c, c=1.0):
        """
        Initialize solver.
        
        Args:
            gamma_c: Ratio of specific heats (adiabatic index)
            c: Speed of light (default=1.0 for natural units, 299792.458 km/s for SI)
        """
        self.gamma_c = gamma_c
        self.c = c
        
        # Left state
        self.nl = None  # Rest frame baryon number density
        self.Pl = None  # Pressure
        self.vl = None  # Velocity
        self.Nl = None  # Lab frame baryon number density
        self.Hl = None  # Enthalpy per baryon
        self.csl = None # Sound speed
        self.gammal = None  # Lorentz factor
        
        # Right state  
        self.nr = None
        self.Pr = None
        self.vr = None
        self.Nr = None
        self.Hr = None
        self.csr = None
        self.gammar = None
        
        # Left star state
        self.nls = None
        self.Pls = None
        self.vls = None
        self.Nls = None
        self.Hls = None
        self.csls = None
        self.vshockl = None
        
        # Right star state
        self.nrs = None
        self.Prs = None
        self.vrs = None
        self.Nrs = None
        self.Hrs = None
        self.csrs = None
        self.vshockr = None
        
    def set_initial_states(self, Pl, nl, vl, Pr, nr, vr):
        """
        Set initial left and right states.
        
        Args:
            Pl, Pr: Pressure (left, right)
            nl, nr: Baryon number density in rest frame (left, right)
            vl, vr: Velocity (left, right)
        """
        self.Pl = Pl
        self.nl = nl
        self.vl = vl
        
        self.Pr = Pr
        self.nr = nr
        self.vr = vr
        
        # Compute derived quantities for left state
        # Thermal energy per baryon: u = P/((gamma_c-1)*n)
        ul = Pl / ((self.gamma_c - 1.0) * nl)
        
        # Enthalpy per baryon: H = 1 + u/c² + P/(nc²)
        self.Hl = 1.0 + ul/(self.c**2) + Pl/(nl * self.c**2)
        
        # Lorentz factor: γ = 1/sqrt(1 - v²/c²)
        self.gammal = 1.0 / np.sqrt(1.0 - (vl/self.c)**2)
        
        # Lab frame baryon number density: N = γn
        self.Nl = self.gammal * nl
        
        # Sound speed: cs = sqrt(gamma_c * P / (n*H))
        self.csl = np.sqrt(self.gamma_c * Pl / (nl * self.Hl))
        
        # Compute derived quantities for right state
        ur = Pr / ((self.gamma_c - 1.0) * nr)
        self.Hr = 1.0 + ur/(self.c**2) + Pr/(nr * self.c**2)
        self.gammar = 1.0 / np.sqrt(1.0 - (vr/self.c)**2)
        self.Nr = self.gammar * nr
        self.csr = np.sqrt(self.gamma_c * Pr / (nr * self.Hr))
        
    def get_velocity(self, P, na, Pa, Ha, csa, va, gammaa, direction):
        """
        Compute post-wave state from pre-wave state and pressure.
        
        Args:
            P: Post-wave pressure
            na: Pre-wave rest frame baryon density
            Pa: Pre-wave pressure
            Ha: Pre-wave enthalpy per baryon
            csa: Pre-wave sound speed
            va: Pre-wave velocity
            gammaa: Pre-wave Lorentz factor
            direction: 'L' or 'R'
            
        Returns:
            (n, u, H, cs, v, vshock) in post-wave state
        """
        sign = -1.0 if direction == 'L' else 1.0
        c = self.c
        
        if P > Pa:
            # Shock wave
            a = 1.0 + (self.gamma_c - 1.0) * (Pa - P) / (self.gamma_c * P)
            b = 1.0 - a
            c_term = Ha * (Pa - P) / na - Ha**2
            
            if c_term > (b**2 / (4.0 * a)):
                raise ValueError("Unphysical enthalpy in intermediate state")
            
            # Post-shock enthalpy per baryon
            H = (-b + np.sqrt(b**2 - 4.0 * a * c_term)) / (2.0 * a)
            
            # Post-shock rest frame density
            n = self.gamma_c * P / ((self.gamma_c - 1.0) * (H - 1.0))
            
            # Thermal energy per baryon
            u = P / ((self.gamma_c - 1.0) * n)
            
            # Mass flux across shock
            j = sign * np.sqrt((P - Pa) / (Ha/na - H/n))
            
            # Shock velocity
            Na = gammaa * na
            a = j**2 + Na**2
            b = -va * Na**2
            vshock = (-b + sign * j**2 * np.sqrt(1.0 + na**2 / j**2)) / a
            gamma_shock = 1.0 / np.sqrt(1.0 - (vshock/c)**2)
            
            # Post-shock velocity
            a = gamma_shock * (P - Pa) / j + Ha * gammaa * va
            b = Ha * gammaa + (P - Pa) * (gamma_shock * va / j + 1.0 / (na * gammaa))
            v = a / b
            
            # Sound speed
            cs = np.sqrt(self.gamma_c * P / (n * H))
            
        else:
            # Rarefaction wave
            # Polytropic constant
            K = Pa / na**self.gamma_c
            
            # Post-rarefaction rest frame density
            n = (P / K)**(1.0 / self.gamma_c)
            
            # Thermal energy per baryon
            u = P / ((self.gamma_c - 1.0) * n)
            
            # Enthalpy per baryon
            H = 1.0 + u/(c**2) + P/(n * c**2)
            
            # Sound speed
            cs = np.sqrt(self.gamma_c * P / (n * H))
            
            # Velocity from rarefaction relation
            sqgl1 = np.sqrt(self.gamma_c - 1.0)
            A = ((1.0 + va/c) / (1.0 - va/c) *
                 ((sqgl1 + csa/c) / (sqgl1 - csa/c) *
                  (sqgl1 - cs/c) / (sqgl1 + cs/c))**(-sign * 2.0 / sqgl1))
            
            v = c * (A - 1.0) / (A + 1.0)
            vshock = 0.0
            
        return n, u, H, cs, v, vshock
        
    def get_dvel(self, P):
        """
        Compute velocity difference between left and right star states.
        
        Args:
            P: Star pressure
            
        Returns:
            Velocity difference vls - vrs
        """
        # Left wave
        (self.nls, uls, self.Hls, self.csls,
         self.vls, self.vshockl) = self.get_velocity(
            P, self.nl, self.Pl, self.Hl, self.csl, self.vl, self.gammal, 'L')
        
        # Right wave
        (self.nrs, urs, self.Hrs, self.csrs,
         self.vrs, self.vshockr) = self.get_velocity(
            P, self.nr, self.Pr, self.Hr, self.csr, self.vr, self.gammar, 'R')
        
        return self.vls - self.vrs
        
    def get_pressure(self, pmin, pmax, tol=0.0):
        """
        Find star pressure using Brent's method.
        
        Args:
            pmin, pmax: Pressure bracket
            tol: Tolerance
            
        Returns:
            Star pressure
        """
        # Machine precision
        eps = np.finfo(float).eps
        
        # Initialize
        a, b = pmin, pmax
        fa = self.get_dvel(a)
        fb = self.get_dvel(b)
        
        c, fc = a, fa
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
            
            # Interpolation
            if abs(e) < tol1 or abs(fa) <= abs(fb):
                d = xm
                e = d
            else:
                if a != c:
                    # Inverse quadratic
                    q, r, s = fa/fc, fb/fc, fb/fa
                    p = s * (2.0*xm*q*(q-r) - (b-a)*(r-1.0))
                    q = (q-1.0) * (r-1.0) * (s-1.0)
                else:
                    # Linear
                    s = fb / fa
                    p = 2.0 * xm * s
                    q = 1.0 - s
                
                if p > 0.0:
                    q = -q
                p = abs(p)
                
                if 2.0*p >= 3.0*xm*q - abs(tol1*q) or p >= abs(0.5*e*q):
                    d = xm
                    e = d
                else:
                    e = d
                    d = p / q
            
            a, fa = b, fb
            if abs(d) > tol1:
                b += d
            else:
                b += np.sign(xm) * tol1
            
            fb = self.get_dvel(b)
            
            if fb * (fc / abs(fc)) > 0.0:
                c, fc = a, fa
                d = b - a
                e = d
                
    def rarefaction(self, xi, na, Pa, ua, csa, va, direction):
        """
        Compute state within rarefaction fan.
        
        Args:
            xi: Similarity variable (x-x0)/t
            na, Pa, ua, csa, va: Pre-rarefaction state
            direction: 'L' or 'R'
            
        Returns:
            (n, P, u, v) at xi
        """
        sign = 1.0 if direction == 'L' else -1.0
        c = self.c
        
        b = np.sqrt(self.gamma_c - 1.0)
        C = (b + csa/c) / (b - csa/c)
        d = -sign * b / 2.0
        k = (1.0 + xi/c) / (1.0 - xi/c)
        l = C * k**d
        v_factor = ((1.0 - va/c) / (1.0 + va/c))**d
        
        # Newton iteration for sound speed
        ocs2 = csa / c
        for _ in range(100):
            fcs2 = (l * v_factor * (1.0 + sign*ocs2)**d * (ocs2 - b) +
                    (1.0 - sign*ocs2)**d * (ocs2 + b))
            
            dfdcs2 = (l * v_factor * (1.0 + sign*ocs2)**d *
                      (1.0 + sign*d*(ocs2-b)/(1.0+sign*ocs2)) +
                      (1.0 - sign*ocs2)**d *
                      (1.0 - sign*d*(ocs2+b)/(1.0-sign*ocs2)))
            
            cs2 = ocs2 - fcs2 / dfdcs2
            
            if abs(cs2 - ocs2) / max(abs(ocs2), 1e-10) <= 5e-7:
                break
            ocs2 = cs2
        
        # Velocity
        v = c * (xi/c + sign*cs2) / (1.0 + sign*xi*cs2/c)
        
        # Density
        n = na * ((cs2**2 * (self.gamma_c - 1.0 - (csa/c)**2)) /
                  ((csa/c)**2 * (self.gamma_c - 1.0 - cs2**2)))**(1.0/(self.gamma_c-1.0))
        
        # Pressure
        P = cs2**2 * (self.gamma_c - 1.0) * n / (self.gamma_c - 1.0 - cs2**2) / self.gamma_c * c**2
        
        # Thermal energy
        u = P / ((self.gamma_c - 1.0) * n)
        
        return n, P, u, v
        
    def solve(self, t, x0=0.5, n_points=400, tol=0.0):
        """
        Solve Riemann problem at time t.
        
        Args:
            t: Time
            x0: Initial discontinuity position
            n_points: Number of grid points
            tol: Tolerance for pressure solver
            
        Returns:
            (x, P, n, N, v, u, gamma, S, e) arrays
        """
        # Find pressure bounds
        pmin = (self.Pl + self.Pr) / 2.0
        pmax = pmin
        
        for _ in range(100):
            pmin = 0.5 * max(pmin, 0.0)
            pmax = 2.0 * pmax
            
            dvel1 = self.get_dvel(pmin)
            dvel2 = self.get_dvel(pmax)
            
            if dvel1 * dvel2 <= 0.0:
                break
        
        # Find star pressure
        Ps = self.get_pressure(pmin, pmax, tol)
        vs = 0.5 * (self.vls + self.vrs)
        
        # Wave positions
        if self.Pl >= Ps:
            # Left rarefaction
            x1 = x0 + (self.vl - self.csl) / (1.0 - self.vl*self.csl/self.c**2) * t
            gamma_ls = 1.0 / np.sqrt(1.0 - (vs/self.c)**2)
            x2 = x0 + (vs - self.csls) / (1.0 - vs*self.csls/self.c**2) * t
        else:
            # Left shock
            x1 = x0 + self.vshockl * t
            x2 = x1
        
        x3 = x0 + vs * t
        
        if self.Pr >= Ps:
            # Right rarefaction
            gamma_rs = 1.0 / np.sqrt(1.0 - (vs/self.c)**2)
            x4 = x0 + (vs + self.csrs) / (1.0 + vs*self.csrs/self.c**2) * t
            x5 = x0 + (self.vr + self.csr) / (1.0 + self.vr*self.csr/self.c**2) * t
        else:
            # Right shock
            x4 = x0 + self.vshockr * t
            x5 = x4
        
        # Solution on grid
        x = np.linspace(0, 1, n_points)
        P = np.zeros(n_points)
        n = np.zeros(n_points)
        N = np.zeros(n_points)
        v = np.zeros(n_points)
        u = np.zeros(n_points)
        gamma_array = np.zeros(n_points)
        S = np.zeros(n_points)
        e = np.zeros(n_points)
        
        for i in range(n_points):
            if x[i] <= x1:
                # Left state
                P[i] = self.Pl
                n[i] = self.nl
                v[i] = self.vl
                u[i] = self.Pl / ((self.gamma_c - 1.0) * self.nl)
                gamma_array[i] = self.gammal
                N[i] = self.Nl
                
            elif x[i] <= x2:
                # Left rarefaction fan
                xi = (x[i] - x0) / t
                n[i], P[i], u[i], v[i] = self.rarefaction(
                    xi, self.nl, self.Pl, self.Pl/((self.gamma_c-1)*self.nl),
                    self.csl, self.vl, 'L')
                gamma_array[i] = 1.0 / np.sqrt(1.0 - (v[i]/self.c)**2)
                N[i] = gamma_array[i] * n[i]
                
            elif x[i] <= x3:
                # Left star state
                P[i] = Ps
                n[i] = self.nls
                v[i] = vs
                u[i] = Ps / ((self.gamma_c - 1.0) * self.nls)
                gamma_array[i] = 1.0 / np.sqrt(1.0 - (vs/self.c)**2)
                N[i] = gamma_array[i] * n[i]
                
            elif x[i] <= x4:
                # Right star state
                P[i] = Ps
                n[i] = self.nrs
                v[i] = vs
                u[i] = Ps / ((self.gamma_c - 1.0) * self.nrs)
                gamma_array[i] = 1.0 / np.sqrt(1.0 - (vs/self.c)**2)
                N[i] = gamma_array[i] * n[i]
                
            elif x[i] <= x5:
                # Right rarefaction fan
                xi = (x[i] - x0) / t
                n[i], P[i], u[i], v[i] = self.rarefaction(
                    xi, self.nr, self.Pr, self.Pr/((self.gamma_c-1)*self.nr),
                    self.csr, self.vr, 'R')
                gamma_array[i] = 1.0 / np.sqrt(1.0 - (v[i]/self.c)**2)
                N[i] = gamma_array[i] * n[i]
                
            else:
                # Right state
                P[i] = self.Pr
                n[i] = self.nr
                v[i] = self.vr
                u[i] = self.Pr / ((self.gamma_c - 1.0) * self.nr)
                gamma_array[i] = self.gammar
                N[i] = self.Nr
            
            # Compute canonical momentum and energy
            # S = γHv
            H = 1.0 + u[i]/(self.c**2) + P[i]/(n[i]*self.c**2)
            S[i] = gamma_array[i] * H * v[i]
            # e = γH - P/(Nc²)
            e[i] = gamma_array[i] * H - P[i]/(N[i]*self.c**2)
        
        return x, P, n, N, v, u, gamma_array, S, e
        
    def save_solution(self, filename, x, P, n, N, v, u, gamma, S, e, t):
        """Save solution to file."""
        np.savetxt(filename, 
                   np.column_stack([x, P, n, N, v, u, gamma, S, e]),
                   header=f'Time: {t}\nColumns: x, P, n, N, v, u, gamma, S, e')
