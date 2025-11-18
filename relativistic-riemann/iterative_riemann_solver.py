"""
Iterative Relativistic Riemann Solver

This module implements an iterative relativistic Riemann solver based on
the van Leer formulation, adapted from the Kitajima solver framework.
Uses Newton-Raphson iteration to solve for the star pressure.

Physical variables:
- n: Baryon number density in rest frame  
- N: Baryon number density in lab frame = γn
- P: Pressure
- v: Velocity
- γ: Lorentz factor = 1/sqrt(1 - v²/c²)
- H: Enthalpy per baryon = 1 + u/c² + P/(nc²)
- c: Speed of light (explicit)
"""

import numpy as np


class IterativeRiemannSolver:
    """
    Iterative relativistic Riemann solver using Newton-Raphson method.
    
    This solver iteratively finds the star pressure that satisfies
    the velocity continuity condition across the contact discontinuity.
    """
    
    def __init__(self, gamma_c, c=1.0, max_iter=100, tol=1e-10, verbose=False):
        """
        Initialize solver.
        
        Args:
            gamma_c: Ratio of specific heats (adiabatic index)
            c: Speed of light (default=1.0 for natural units)
            max_iter: Maximum number of Newton-Raphson iterations
            tol: Convergence tolerance for pressure
            verbose: Print detailed iteration information
        """
        self.gamma_c = gamma_c
        self.c = c
        self.max_iter = max_iter
        self.tol = tol
        self.verbose = verbose
        
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
        
        # Star states
        self.nls = None
        self.Pls = None
        self.vls = None
        self.Nls = None
        self.Hls = None
        self.csls = None
        self.vshockl = None
        
        self.nrs = None
        self.Prs = None
        self.vrs = None
        self.Nrs = None
        self.Hrs = None
        self.csrs = None
        self.vshockr = None
        
        # Iteration tracking
        self.iterations = 0
        self.pressure_history = []
        self.residual_history = []
        self.convergence_info = {}  # Store detailed convergence information
        
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
        ul = Pl / ((self.gamma_c - 1.0) * nl)
        self.Hl = 1.0 + ul/(self.c**2) + Pl/(nl * self.c**2)
        self.gammal = 1.0 / np.sqrt(1.0 - (vl/self.c)**2)
        self.Nl = self.gammal * nl
        self.csl = np.sqrt(self.gamma_c * Pl / (nl * self.Hl))
        
        # Compute derived quantities for right state
        ur = Pr / ((self.gamma_c - 1.0) * nr)
        self.Hr = 1.0 + ur/(self.c**2) + Pr/(nr * self.c**2)
        self.gammar = 1.0 / np.sqrt(1.0 - (vr/self.c)**2)
        self.Nr = self.gammar * nr
        self.csr = np.sqrt(self.gamma_c * Pr / (nr * self.Hr))
        
    def get_velocity_and_derivative(self, P, na, Pa, Ha, csa, va, gammaa, direction):
        """
        Compute post-wave velocity and its derivative with respect to pressure.
        Uses adaptive finite difference for robust derivative calculation.
        
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
            (n, v, dv_dP, H, cs, vshock) in post-wave state
        """
        sign = -1.0 if direction == 'L' else 1.0
        c = self.c
        
        if P > Pa:
            # Shock wave - use Rankine-Hugoniot relations
            a = 1.0 + (self.gamma_c - 1.0) * (Pa - P) / (self.gamma_c * P)
            b = 1.0 - a
            c_term = Ha * (Pa - P) / na - Ha**2
            
            if c_term > (b**2 / (4.0 * a)):
                raise ValueError("Unphysical enthalpy in intermediate state")
            
            # Post-shock enthalpy
            H = (-b + np.sqrt(b**2 - 4.0 * a * c_term)) / (2.0 * a)
            
            # Post-shock density
            n = self.gamma_c * P / ((self.gamma_c - 1.0) * (H - 1.0))
            
            # Mass flux
            j = sign * np.sqrt((P - Pa) / (Ha/na - H/n))
            
            # Shock velocity
            Na = gammaa * na
            a_shock = j**2 + Na**2
            b_shock = -va * Na**2
            vshock = (-b_shock + sign * j**2 * np.sqrt(1.0 + na**2 / j**2)) / a_shock
            gamma_shock = 1.0 / np.sqrt(1.0 - (vshock/c)**2)
            
            # Post-shock velocity
            a_vel = gamma_shock * (P - Pa) / j + Ha * gammaa * va
            b_vel = Ha * gammaa + (P - Pa) * (gamma_shock * va / j + 1.0 / (na * gammaa))
            v = a_vel / b_vel
            
            # Sound speed
            cs = np.sqrt(self.gamma_c * P / (n * H))
            
            # Adaptive finite difference for derivative
            # Scale perturbation with pressure magnitude and convergence state
            rel_p = max(abs(P - Pa) / Pa, 0.01)
            dP = max(1e-8 * P * rel_p, 1e-12)
            P_pert = P + dP
            
            # Perturbed enthalpy
            a_pert = 1.0 + (self.gamma_c - 1.0) * (Pa - P_pert) / (self.gamma_c * P_pert)
            b_pert = 1.0 - a_pert
            c_term_pert = Ha * (Pa - P_pert) / na - Ha**2
            
            if c_term_pert > (b_pert**2 / (4.0 * a_pert)):
                # Use smaller perturbation
                dP = dP * 0.1
                P_pert = P + dP
                a_pert = 1.0 + (self.gamma_c - 1.0) * (Pa - P_pert) / (self.gamma_c * P_pert)
                b_pert = 1.0 - a_pert
                c_term_pert = Ha * (Pa - P_pert) / na - Ha**2
            
            H_pert = (-b_pert + np.sqrt(b_pert**2 - 4.0 * a_pert * c_term_pert)) / (2.0 * a_pert)
            n_pert = self.gamma_c * P_pert / ((self.gamma_c - 1.0) * (H_pert - 1.0))
            
            j_pert = sign * np.sqrt((P_pert - Pa) / (Ha/na - H_pert/n_pert))
            vshock_pert = (-b_shock + sign * j_pert**2 * np.sqrt(1.0 + na**2 / j_pert**2)) / (j_pert**2 + Na**2)
            gamma_shock_pert = 1.0 / np.sqrt(1.0 - (vshock_pert/c)**2)
            
            a_vel_pert = gamma_shock_pert * (P_pert - Pa) / j_pert + Ha * gammaa * va
            b_vel_pert = Ha * gammaa + (P_pert - Pa) * (gamma_shock_pert * va / j_pert + 1.0 / (na * gammaa))
            v_pert = a_vel_pert / b_vel_pert
            
            dv_dP = (v_pert - v) / dP
            
        else:
            # Rarefaction wave - use isentropic relations
            K = Pa / na**self.gamma_c
            n = (P / K)**(1.0 / self.gamma_c)
            
            u = P / ((self.gamma_c - 1.0) * n)
            H = 1.0 + u/(c**2) + P/(n * c**2)
            cs = np.sqrt(self.gamma_c * P / (n * H))
            
            # Velocity from rarefaction invariant
            sqgl1 = np.sqrt(self.gamma_c - 1.0)
            A = ((1.0 + va/c) / (1.0 - va/c) *
                 ((sqgl1 + csa/c) / (sqgl1 - csa/c) *
                  (sqgl1 - cs/c) / (sqgl1 + cs/c))**(-sign * 2.0 / sqgl1))
            
            v = c * (A - 1.0) / (A + 1.0)
            vshock = 0.0
            
            # Adaptive derivative for rarefaction
            rel_p = max(abs(P - Pa) / Pa, 0.01)
            dP = max(1e-8 * P * rel_p, 1e-12)
            P_pert = P + dP
            n_pert = (P_pert / K)**(1.0 / self.gamma_c)
            u_pert = P_pert / ((self.gamma_c - 1.0) * n_pert)
            H_pert = 1.0 + u_pert/(c**2) + P_pert/(n_pert * c**2)
            cs_pert = np.sqrt(self.gamma_c * P_pert / (n_pert * H_pert))
            
            A_pert = ((1.0 + va/c) / (1.0 - va/c) *
                      ((sqgl1 + csa/c) / (sqgl1 - csa/c) *
                       (sqgl1 - cs_pert/c) / (sqgl1 + cs_pert/c))**(-sign * 2.0 / sqgl1))
            v_pert = c * (A_pert - 1.0) / (A_pert + 1.0)
            
            dv_dP = (v_pert - v) / dP
            
        return n, v, dv_dP, H, cs, vshock
        
    def iterate_pressure(self, P_guess):
        """
        Single Newton-Raphson iteration to improve pressure estimate.
        Follows robust formulation inspired by classic Newtonian solver.
        
        Args:
            P_guess: Current pressure guess
            
        Returns:
            (P_new, residual, converged)
        """
        # Get velocities and derivatives from both sides
        (self.nls, self.vls, dvl_dP, self.Hls, 
         self.csls, self.vshockl) = self.get_velocity_and_derivative(
            P_guess, self.nl, self.Pl, self.Hl, self.csl, 
            self.vl, self.gammal, 'L')
        
        (self.nrs, self.vrs, dvr_dP, self.Hrs, 
         self.csrs, self.vshockr) = self.get_velocity_and_derivative(
            P_guess, self.nr, self.Pr, self.Hr, self.csr, 
            self.vr, self.gammar, 'R')
        
        # Residual: velocity difference (f = vls - vrs should be zero)
        f = self.vls - self.vrs
        
        # Derivative of residual with respect to pressure
        df_dP = dvl_dP - dvr_dP
        
        # Robust Newton-Raphson update with safeguards
        if abs(df_dP) < 1e-20:
            # Derivative too small - likely at solution or singular point
            P_new = P_guess
            converged = True
        else:
            # Standard Newton-Raphson step
            delta_P = -f / df_dP
            
            # Damping factor for stability (inspired by classic solver's implicit damping)
            # Limit pressure change to prevent overshooting
            max_factor = 2.0
            if abs(delta_P) > max_factor * P_guess:
                delta_P = max_factor * P_guess * np.sign(delta_P)
            
            P_new = P_guess + delta_P
            
            # Enforce physical bounds (similar to classic solver's smallp)
            smallp = 1e-25
            P_new = max(P_new, smallp)
            P_new = max(P_new, 0.01 * min(self.Pl, self.Pr))
            P_new = min(P_new, 100.0 * max(self.Pl, self.Pr))
            
            # Check convergence using both relative pressure change and residual
            # (dual criterion like classic solver checking both iteration change and function value)
            rel_change = abs(P_new - P_guess) / max(P_guess, smallp)
            rel_residual = abs(f) / max(abs(self.vls), abs(self.vrs), 1e-10)
            
            converged = (rel_change < self.tol) or (rel_residual < self.tol)
        
        return P_new, f, converged
        
    def solve_star_pressure(self):
        """
        Iteratively solve for star pressure using Newton-Raphson method.
        Uses robust initial guess based on acoustic approximation.
        
        Returns:
            Star pressure
        """
        # Robust initial guess using acoustic approximation (inspired by classic solver)
        # This is the relativistic equivalent of the Newtonian acoustic guess
        c = self.c
        
        # Compute acoustic impedances (relativistic sound speeds times densities)
        Zl = self.nl * self.csl  # Left acoustic impedance
        Zr = self.nr * self.csr  # Right acoustic impedance
        
        # Acoustic approximation: P* = Pl + (Pr - Pl - Zr*(vr - vl)) * Zl/(Zl + Zr)
        # This assumes small velocity difference and linear wave propagation
        P_guess = self.Pr - self.Pl - Zr * (self.vr - self.vl)
        P_guess = self.Pl + P_guess * Zl / (Zl + Zr)
        
        # Ensure initial guess is positive and reasonable
        smallp = 1e-25
        P_guess = max(P_guess, smallp)
        P_guess = max(P_guess, 0.1 * min(self.Pl, self.Pr))
        P_guess = min(P_guess, 10.0 * max(self.Pl, self.Pr))
        
        # Reset iteration tracking
        self.iterations = 0
        self.pressure_history = [P_guess]
        self.residual_history = []
        
        if self.verbose:
            print(f"\n=== Robust Iterative Solver ===")
            print(f"Initial states:")
            print(f"  Left:  P={self.Pl:.6e}, n={self.nl:.6e}, v={self.vl:.6e}, cs={self.csl:.6e}")
            print(f"  Right: P={self.Pr:.6e}, n={self.nr:.6e}, v={self.vr:.6e}, cs={self.csr:.6e}")
            print(f"Acoustic impedances: Zl={Zl:.6e}, Zr={Zr:.6e}")
            print(f"Initial guess: P*={P_guess:.6e}")
            print(f"\nIteration history:")
        
        # Newton-Raphson iteration
        for i in range(self.max_iter):
            P_new, residual, converged = self.iterate_pressure(P_guess)
            
            self.iterations = i + 1
            self.pressure_history.append(P_new)
            self.residual_history.append(abs(residual))
            
            if self.verbose:
                print(f"  {i+1:3d}: P*={P_new:.6e}, |Δv|={abs(residual):.3e}, "
                      f"ΔP={abs(P_new-P_guess):.3e}")
            
            if converged:
                if self.verbose:
                    print(f"✓ Converged in {self.iterations} iterations")
                return P_new
            
            P_guess = P_new
        
        # Maximum iterations reached
        print(f"Warning: Maximum iterations ({self.max_iter}) reached")
        print(f"  Final residual: {abs(residual):.3e}")
        return P_new
    
    def get_convergence_diagnostics(self):
        """
        Get detailed convergence diagnostics.
        
        Returns:
            Dictionary with convergence information
        """
        if len(self.pressure_history) < 2:
            return {"error": "No iteration history available"}
        
        diagnostics = {
            "iterations": self.iterations,
            "initial_pressure": self.pressure_history[0],
            "final_pressure": self.pressure_history[-1],
            "pressure_change": abs(self.pressure_history[-1] - self.pressure_history[0]),
            "final_residual": self.residual_history[-1] if self.residual_history else 0.0,
            "convergence_rate": [],
        }
        
        # Compute convergence rate (quadratic for Newton-Raphson in ideal case)
        if len(self.residual_history) > 2:
            rates = []
            for i in range(1, len(self.residual_history) - 1):
                if self.residual_history[i] > 1e-15:
                    rate = np.log(self.residual_history[i+1]) - np.log(self.residual_history[i])
                    rates.append(rate)
            if rates:
                diagnostics["convergence_rate"] = rates
                diagnostics["avg_convergence_rate"] = np.mean(rates)
        
        return diagnostics
        
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
        
    def solve(self, t, x0=0.5, n_points=400):
        """
        Solve Riemann problem at time t using iterative method.
        
        Args:
            t: Time
            x0: Initial discontinuity position
            n_points: Number of grid points
            
        Returns:
            (x, P, n, N, v, u, gamma, S, e) arrays
        """
        # Solve for star pressure iteratively
        Ps = self.solve_star_pressure()
        vs = 0.5 * (self.vls + self.vrs)
        
        print(f"Iterative solver converged in {self.iterations} iterations")
        print(f"  Star pressure: Ps = {Ps:.6f}")
        print(f"  Star velocity: vs = {vs:.6f}")
        print(f"  Final residual: {self.residual_history[-1]:.3e}")
        
        # Compute wave positions (same as Kitajima solver)
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
            H = 1.0 + u[i]/(self.c**2) + P[i]/(n[i]*self.c**2)
            S[i] = gamma_array[i] * H * v[i]
            e[i] = gamma_array[i] * H - P[i]/(N[i]*self.c**2)
        
        return x, P, n, N, v, u, gamma_array, S, e
        
    def save_solution(self, filename, x, P, n, N, v, u, gamma, S, e, t):
        """Save solution to file."""
        np.savetxt(filename, 
                   np.column_stack([x, P, n, N, v, u, gamma, S, e]),
                   header=f'Time: {t}\nIterations: {self.iterations}\nColumns: x, P, n, N, v, u, gamma, S, e')
