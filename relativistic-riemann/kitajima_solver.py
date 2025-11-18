"""
Kitajima Formulation: Special Relativistic Riemann Solver

This module implements the relativistic Riemann solver using the Kitajima formulation
(Kitajima et al., 2025, arXiv:2510.18251v1) which uses:
- Baryon number density N (not mass density)
- Explicit light speed c (not c=1)
- Canonical momentum per baryon S
- Canonical energy per baryon e

Enhanced with tangential velocity support following Pons, Martí & Müller (1999)
"The exact solution of the Riemann problem with non-zero tangential velocities
in relativistic hydrodynamics" (J. Fluid Mech. 1999)

TANGENTIAL VELOCITY IMPLEMENTATION DETAILS
==========================================

Physical Background (from Pons et al. 1999):
--------------------------------------------
In relativistic hydrodynamics, tangential velocities (perpendicular to shock normal)
are coupled to the flow through:
1. The Lorentz factor: γ = 1/sqrt(1 - v²/c²) where v² = (v^x)² + (v^y)² + (v^z)²
2. The specific enthalpy: h (or H in Kitajima notation)

Unlike Newtonian hydrodynamics where tangential velocities are constant across
waves, relativistic flows exhibit tangential velocity changes governed by
energy-momentum conservation.

Conservation Laws Across Waves:
-------------------------------

1. RAREFACTION WAVES (Pons et al. eqs. 27-28):
   - Conserved quantity: h*W*v^y = constant
   - Conserved quantity: h*W*v^z = constant
   - Direction preserved: v^y/v^z = constant
   - Magnitude changes: v^t = sqrt((v^y)² + (v^z)²) varies
   
   These arise from the simple wave ODE (eq. 30) and the algebraic
   constraints (eqs. 27-28) for self-similar rarefaction solutions.

2. SHOCK WAVES (Pons et al. eqs. 54-55):
   - Conserved quantity: S^y/D = (ρhW²v^y)/(ρW) = hWv^y = constant
   - Conserved quantity: S^z/D = (ρhW²v^z)/(ρW) = hWv^z = constant
   - Direction preserved: v^y/v^z = constant
   - Magnitude changes: v^t varies
   
   These follow from Rankine-Hugoniot jump conditions across shock fronts.

3. CONTACT DISCONTINUITY:
   - Normal velocity continuous: v^x_L* = v^x_R*
   - Pressure continuous: P_L* = P_R*
   - Tangential velocities can jump: v^y_L* ≠ v^y_R*, v^z_L* ≠ v^z_R*

Implementation Strategy:
-----------------------
For both rarefactions and shocks, we solve iteratively for tangential
velocities because:
1. Conservation laws give v^t as function of enthalpy and Lorentz factor
2. Lorentz factor depends on total velocity including v^t
3. This creates implicit equation requiring iterative solution

Iteration procedure:
1. Start with approximate v^t (based on enthalpy ratio)
2. Compute Lorentz factor γ = 1/sqrt(1 - (v^x² + v^t²)/c²)
3. Update v^t from conservation law
4. Repeat until convergence (typically <50 iterations)
5. Decompose v^t into (v^y, v^z) preserving direction

Notes on formulation equivalence:
- Kitajima: N = γn (lab frame baryon density), H (enthalpy per baryon)
- Pons et al: ρ = n (proper rest-mass density), h (specific enthalpy)  
- Relation: ρ in Pons et al. ≡ n in Kitajima; h ≡ H (for ideal gas)
- D = ρW = nγ (lab frame mass density)
- S^i = ρhW²v^i (momentum density)

Physical Interpretation:
-----------------------
The tangential velocity changes across waves because:
1. Energy-momentum tensor couples all velocity components
2. In ultrarelativistic regime (γ >> 1), even small tangential velocities
   significantly affect the energy density
3. Thermodynamically relativistic flows (h > 1) also show coupling even
   at moderate Lorentz factors

Testing:
--------
Test cases from Pons et al. (1999) Table 1 show that tangential velocities
can significantly modify:
- Star region pressure P* (increases with v^t_R)
- Normal velocity v^x_* (decreases with v^t_L)
- Wave speeds (all tend to zero as v^t → sqrt(1-v^x²))
"""

import numpy as np


class KitajimaRiemannSolver:
    """
    Relativistic Riemann solver using Kitajima formulation with baryon number density.
    
    Physical variables (Kitajima formulation):
    - N: Baryon number density in lab frame = γn
    - n: Baryon number density in rest frame  
    - S: Relativistic canonical momentum per baryon = γHv
    - e: Canonical energy per baryon = γH - P/(Nc²)
    - H: Enthalpy per baryon = 1 + u/c² + P/(nc²)
    - γ: Lorentz factor = 1/sqrt(1 - v²/c²)
    - c: Speed of light (explicit, not set to 1)
    
    Velocity components:
    - v^x: Normal velocity (along shock/discontinuity normal direction)
    - v^y: Tangential velocity (first perpendicular component)
    - v^z: Tangential velocity (second perpendicular component)
    - v²: Total velocity squared = (v^x)² + (v^y)² + (v^z)²
    
    Tangential velocity conservation (Pons, Martí & Müller 1999):
    -------------------------------------------------------------
    
    Rarefaction waves (eqs. 27-28, rar2, rar3):
        h*W*v^y = const
        h*W*v^z = const
        → Direction v^y/v^z preserved, magnitude varies
    
    Shock waves (eqs. 54-55, rhvy, rhvz):
        S^y/D = const  ⟺  (ρhW²v^y)/(ρW) = hWv^y = const
        S^z/D = const  ⟺  (ρhW²v^z)/(ρW) = hWv^z = const
        → Direction v^y/v^z preserved, magnitude varies
    
    Contact discontinuity:
        v^x continuous, P continuous
        v^y, v^z can jump
    
    Key physical insight:
    - Tangential velocities couple through Lorentz factor and enthalpy
    - Even thermodynamically relativistic flows (h>1) with γ≈1 show coupling
    - Ultrarelativistic tangential flow (v^t→c) strongly affects solution
    - All wave speeds → 0 as v^t → sqrt(1-(v^x)²) (maximum tangential velocity)
    
    Implementation notes:
    - Iterative solution for tangential velocities (γ depends on total v)
    - Under-relaxation for numerical stability (factor 0.3)
    - Convergence tolerance: 1e-10
    - Maximum iterations: 50
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
        self.vl = None  # Normal velocity (v^x)
        self.vyl = None # Tangential velocity y-component (v^y)
        self.vzl = None # Tangential velocity z-component (v^z)
        self.Nl = None  # Lab frame baryon number density
        self.Hl = None  # Enthalpy per baryon
        self.csl = None # Sound speed
        self.gammal = None  # Lorentz factor
        
        # Right state  
        self.nr = None
        self.Pr = None
        self.vr = None
        self.vyr = None
        self.vzr = None
        self.Nr = None
        self.Hr = None
        self.csr = None
        self.gammar = None
        
        # Left star state
        self.nls = None
        self.Pls = None
        self.vls = None
        self.vyls = None
        self.vzls = None
        self.Nls = None
        self.Hls = None
        self.csls = None
        self.vshockl = None
        
        # Right star state
        self.nrs = None
        self.Prs = None
        self.vrs = None
        self.vyrs = None
        self.vzrs = None
        self.Nrs = None
        self.Hrs = None
        self.csrs = None
        self.vshockr = None
        
    def set_initial_states(self, Pl, nl, vl, Pr, nr, vr, vyl=0.0, vzl=0.0, vyr=0.0, vzr=0.0):
        """
        Set initial left and right states.
        
        Args:
            Pl, Pr: Pressure (left, right)
            nl, nr: Baryon number density in rest frame (left, right)
            vl, vr: Normal velocity v^x (left, right)
            vyl, vyr: Tangential velocity v^y (left, right) [default=0]
            vzl, vzr: Tangential velocity v^z (left, right) [default=0]
        """
        self.Pl = Pl
        self.nl = nl
        self.vl = vl
        self.vyl = vyl
        self.vzl = vzl
        
        self.Pr = Pr
        self.nr = nr
        self.vr = vr
        self.vyr = vyr
        self.vzr = vzr
        
        # Compute derived quantities for left state
        # Total velocity squared: v² = (v^x)² + (v^y)² + (v^z)²
        v2_l = vl**2 + vyl**2 + vzl**2
        
        # Thermal energy per baryon: u = P/((gamma_c-1)*n)
        ul = Pl / ((self.gamma_c - 1.0) * nl)
        
        # Enthalpy per baryon: H = 1 + u/c² + P/(nc²)
        self.Hl = 1.0 + ul/(self.c**2) + Pl/(nl * self.c**2)
        
        # Lorentz factor: γ = 1/sqrt(1 - v²/c²)
        self.gammal = 1.0 / np.sqrt(1.0 - v2_l/(self.c**2))
        
        # Lab frame baryon number density: N = γn
        self.Nl = self.gammal * nl
        
        # Sound speed: cs = sqrt(gamma_c * P / (n*H))
        self.csl = np.sqrt(self.gamma_c * Pl / (nl * self.Hl))
        
        # Compute derived quantities for right state
        v2_r = vr**2 + vyr**2 + vzr**2
        ur = Pr / ((self.gamma_c - 1.0) * nr)
        self.Hr = 1.0 + ur/(self.c**2) + Pr/(nr * self.c**2)
        self.gammar = 1.0 / np.sqrt(1.0 - v2_r/(self.c**2))
        self.Nr = self.gammar * nr
        self.csr = np.sqrt(self.gamma_c * Pr / (nr * self.Hr))
        
    def get_velocity(self, P, na, Pa, Ha, csa, va, vya, vza, gammaa, direction):
        """
        Compute post-wave state from pre-wave state and pressure.
        
        Following Pons et al. (1999, J. Fluid Mech.) with tangential velocity handling:
        - For rarefactions: h*W*v^y = const, h*W*v^z = const (eqs. 27-28, rar2, rar3)
        - For shocks: S^y/D = const, S^z/D = const (eqs. 54-55, rhvy, rhvz)
          which translates to: (h*W*v^y)/(γn) = const, (h*W*v^z)/(γn) = const
        
        Key insight from paper (section 3):
        - Tangential velocity direction is preserved: v^y/v^z = const
        - Only magnitude changes according to energy-momentum conservation
        
        Args:
            P: Post-wave pressure
            na: Pre-wave rest frame baryon density
            Pa: Pre-wave pressure
            Ha: Pre-wave enthalpy per baryon
            csa: Pre-wave sound speed
            va: Pre-wave normal velocity (v^x)
            vya: Pre-wave tangential velocity (v^y)
            vza: Pre-wave tangential velocity (v^z)
            gammaa: Pre-wave Lorentz factor (includes tangential velocities)
            direction: 'L' or 'R'
            
        Returns:
            (n, u, H, cs, v, vy, vz, vshock) in post-wave state
        """
        sign = -1.0 if direction == 'L' else 1.0
        c = self.c
        
        # Tangential velocity magnitude and direction (direction preserved across waves)
        vta = np.sqrt(vya**2 + vza**2)
        
        if P > Pa:
            # Shock wave - use Taub adiabat (Pons et al. eq. 63)
            a = 1.0 + (self.gamma_c - 1.0) * (Pa - P) / (self.gamma_c * P)
            b_term = -(self.gamma_c - 1.0) * (Pa - P) / (self.gamma_c * P)
            c_term = Ha * (Pa - P) / na - Ha**2
            
            # Solve quadratic equation for post-shock enthalpy H_b
            discriminant = b_term**2 - 4.0 * a * c_term
            if discriminant < 0:
                raise ValueError("Unphysical enthalpy in shock state (negative discriminant)")
            
            H = (-b_term + np.sqrt(discriminant)) / (2.0 * a)
            
            if H <= 0:
                raise ValueError("Unphysical enthalpy in shock state (H <= 0)")
            
            # Post-shock rest frame density from ideal gas EOS
            n = self.gamma_c * P / ((self.gamma_c - 1.0) * (H - 1.0))
            
            # Thermal energy per baryon
            u = P / ((self.gamma_c - 1.0) * n)
            
            # Mass flux from eq. (64): j^2 = -[p]/[h/ρ]
            j2 = -(P - Pa) / (H/n - Ha/na)
            if j2 < 0:
                raise ValueError("Negative mass flux squared in shock")
            j = sign * np.sqrt(j2)
            
            # Shock velocity from eq. (61)
            Na = gammaa * na
            a_shock = j**2 + Na**2
            b_shock = -va * Na**2
            vshock = (-b_shock + sign * j**2 * np.sqrt(1.0 + na**2 / j**2)) / a_shock
            gamma_shock = 1.0 / np.sqrt(1.0 - (vshock/c)**2)
            
            # Post-shock normal velocity (from Rankine-Hugoniot conditions)
            a_v = gamma_shock * (P - Pa) / j + Ha * gammaa * va
            b_v = Ha * gammaa + (P - Pa) * (gamma_shock * va / j + 1.0 / (na * gammaa))
            v = a_v / b_v
            
            # Sound speed
            cs = np.sqrt(self.gamma_c * P / (n * H))
            
            # Post-shock tangential velocities (Pons et al. eqs. 54-55, rhvy, rhvz)
            # Rankine-Hugoniot conditions: S^y/D = const, S^z/D = const
            # where S^y = ρhW²v^y = nhW²v^y and D = ρW = nW
            # Therefore: (hWv^y)/W = hv^y = const across shock
            # But need to account for proper coupling through energy-momentum tensor
            # 
            # From eq. (54-55): [S^y/D] = 0 and [S^z/D] = 0
            # S^y = ρhW²v^y, D = ρW, so S^y/D = hWv^y
            # Therefore: (hWv^y)_a = (hWv^y)_b and (hWv^z)_a = (hWv^z)_b
            #
            # Iterative solution needed since W_b depends on total velocity including v^y, v^z
            
            if vta > 1e-15:
                # Initial guess for tangential velocity magnitude
                vt_mag = min(vta * (Ha * gammaa) / H, 0.9 * c)  # Start with safe value
                
                # Iterative solution for self-consistent tangential velocity
                for iteration in range(50):
                    v2_tot = v**2 + vt_mag**2
                    if v2_tot >= c**2 * 0.9999:
                        # Velocity approaching speed of light - reduce tangential component
                        vt_mag *= 0.90
                        continue
                    
                    # Compute Lorentz factor with current tangential velocity
                    gamma_b = 1.0 / np.sqrt(1.0 - v2_tot/(c**2))
                    
                    # Apply conservation law: (h*W*v^t)_a = (h*W*v^t)_b (Pons eq. 54-55)
                    vt_new = (Ha * gammaa * vta) / (H * gamma_b)
                    
                    # Enforce causality: total velocity must be < c
                    v2_new = v**2 + vt_new**2
                    if v2_new >= c**2 * 0.9999:
                        vt_new = np.sqrt(max(0, c**2 * 0.9999 - v**2))
                    
                    # Check convergence
                    if abs(vt_new - vt_mag) < 1e-10:
                        break
                    
                    # Under-relaxation for stability
                    vt_mag = 0.3 * vt_new + 0.7 * vt_mag
                
                # Preserve direction: v^y/v^z = const (Pons et al. section 3)
                if vta > 1e-15:
                    vy = (vya / vta) * vt_mag
                    vz = (vza / vta) * vt_mag
                else:
                    vy = 0.0
                    vz = 0.0
            else:
                vy = 0.0
                vz = 0.0
            
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
            
            # Normal velocity from rarefaction integral (Pons et al. eq. 30, ODE integration)
            # This equation relates v^x across rarefaction through isentropic flow
            # For ideal gas: integrating eq. (30) gives the velocity-pressure relation
            sqgl1 = np.sqrt(self.gamma_c - 1.0)
            A = ((1.0 + va/c) / (1.0 - va/c) *
                 ((sqgl1 + csa/c) / (sqgl1 - csa/c) *
                  (sqgl1 - cs/c) / (sqgl1 + cs/c))**(-sign * 2.0 / sqgl1))
            
            v = c * (A - 1.0) / (A + 1.0)
            vshock = 0.0
            
            # Post-rarefaction tangential velocities (Pons et al. eqs. 27-28, rar2, rar3)
            # Conservation law: h*W*v^y = const and h*W*v^z = const
            # Equivalent to: (h*W*v^t) = const where v^t = sqrt((v^y)² + (v^z)²)
            # Direction v^y/v^z is preserved across rarefaction
            #
            # From ODE integration (Pons eq. 30), we get v^x as function of P
            # Tangential velocity then follows from h*W*v^t = const
            # Iterative solution needed since W depends on total velocity
            
            if vta > 1e-15:
                # Initial guess for tangential velocity magnitude
                vt_mag = min(vta * (Ha * gammaa) / H, 0.9 * c)  # Start with safe value
                
                # Iterative solution for self-consistent tangential velocity
                for iteration in range(50):
                    v2_tot = v**2 + vt_mag**2
                    if v2_tot >= c**2 * 0.9999:
                        # Velocity approaching speed of light - reduce tangential component
                        vt_mag *= 0.90
                        continue
                    
                    # Compute Lorentz factor with current tangential velocity
                    gamma_post = 1.0 / np.sqrt(1.0 - v2_tot/(c**2))
                    
                    # Apply conservation law: (h*W*v^t)_a = (h*W*v^t)_b (Pons eq. 27-28)
                    vt_new = (Ha * gammaa * vta) / (H * gamma_post)
                    
                    # Enforce causality: total velocity must be < c
                    v2_new = v**2 + vt_new**2
                    if v2_new >= c**2 * 0.9999:
                        vt_new = np.sqrt(max(0, c**2 * 0.9999 - v**2))
                    
                    # Check convergence
                    if abs(vt_new - vt_mag) < 1e-10:
                        break
                    
                    # Under-relaxation for stability
                    vt_mag = 0.3 * vt_new + 0.7 * vt_mag
                
                # Preserve direction: v^y/v^z = const (Pons et al. section 3)
                if vta > 1e-15:
                    vy = (vya / vta) * vt_mag
                    vz = (vza / vta) * vt_mag
                else:
                    vy = 0.0
                    vz = 0.0
            else:
                vy = 0.0
                vz = 0.0
            
        return n, u, H, cs, v, vy, vz, vshock
        
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
         self.vls, self.vyls, self.vzls, self.vshockl) = self.get_velocity(
            P, self.nl, self.Pl, self.Hl, self.csl, self.vl, self.vyl, self.vzl, self.gammal, 'L')
        
        # Right wave
        (self.nrs, urs, self.Hrs, self.csrs,
         self.vrs, self.vyrs, self.vzrs, self.vshockr) = self.get_velocity(
            P, self.nr, self.Pr, self.Hr, self.csr, self.vr, self.vyr, self.vzr, self.gammar, 'R')
        
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
                
    def rarefaction(self, xi, na, Pa, ua, csa, va, vya, vza, direction):
        """
        Compute state within rarefaction fan.
        
        Following Pons et al. (1999, J. Fluid Mech.) with tangential velocity coupling.
        
        Key equations:
        - Self-similarity: ξ = (x-x0)/t relates position to wave structure
        - Characteristic speed (eq. 6): ξ = [v^x(1-c_s²) ± c_s*sqrt(...)]/[1-v²c_s²]
        - Conservation: h*W*v^y = const, h*W*v^z = const (eqs. 27-28)
        - Direction: v^y/v^z = const (tangential velocity direction preserved)
        
        The rarefaction solution involves:
        1. Finding sound speed c_s(ξ) from characteristic equation
        2. Computing normal velocity v^x from isentropic relation
        3. Computing tangential velocities from h*W*v^t = const
        
        Args:
            xi: Similarity variable (x-x0)/t
            na, Pa, ua, csa, va: Pre-rarefaction state (normal velocity)
            vya, vza: Pre-rarefaction tangential velocities
            direction: 'L' or 'R'
            
        Returns:
            (n, P, u, v, vy, vz) at xi
        """
        sign = 1.0 if direction == 'L' else -1.0
        c = self.c
        
        # Tangential velocity magnitude and direction (direction preserved)
        vta = np.sqrt(vya**2 + vza**2)
        
        # Pre-rarefaction total velocity and Lorentz factor
        v2a = va**2 + vta**2
        gammaa = 1.0 / np.sqrt(1.0 - v2a/(c**2))
        
        # Pre-rarefaction enthalpy
        Ha = 1.0 + ua/(c**2) + Pa/(na * c**2)
        
        # Constants for rarefaction integral (Pons et al., based on eq. 30)
        b = np.sqrt(self.gamma_c - 1.0)
        C = (b + csa/c) / (b - csa/c)
        d = -sign * b / 2.0
        k = (1.0 + xi/c) / (1.0 - xi/c)
        l = C * k**d
        v_factor = ((1.0 - va/c) / (1.0 + va/c))**d
        
        # Newton iteration to find sound speed c_s from characteristic equation (eq. 6)
        # The equation couples ξ, v^x, and c_s through the relativistic dispersion relation
        ocs2 = csa / c  # Initial guess: sound speed from pre-rarefaction state
        for iteration in range(100):
            # Characteristic equation residual
            fcs2 = (l * v_factor * (1.0 + sign*ocs2)**d * (ocs2 - b) +
                    (1.0 - sign*ocs2)**d * (ocs2 + b))
            
            # Derivative for Newton step
            dfdcs2 = (l * v_factor * (1.0 + sign*ocs2)**d *
                      (1.0 + sign*d*(ocs2-b)/(1.0+sign*ocs2)) +
                      (1.0 - sign*ocs2)**d *
                      (1.0 - sign*d*(ocs2+b)/(1.0-sign*ocs2)))
            
            cs2 = ocs2 - fcs2 / dfdcs2
            
            # Convergence check
            if abs(cs2 - ocs2) / max(abs(ocs2), 1e-10) <= 5e-7:
                break
            ocs2 = cs2
        
        # Normal velocity from similarity solution (Pons et al., derived from eq. 30)
        # Relates velocity to position through self-similar structure
        v = c * (xi/c + sign*cs2) / (1.0 + sign*xi*cs2/c)
        
        # Density from isentropic relation (polytropic along rarefaction)
        n = na * ((cs2**2 * (self.gamma_c - 1.0 - (csa/c)**2)) /
                  ((csa/c)**2 * (self.gamma_c - 1.0 - cs2**2)))**(1.0/(self.gamma_c-1.0))
        
        # Pressure from sound speed definition
        P = cs2**2 * (self.gamma_c - 1.0) * n / (self.gamma_c - 1.0 - cs2**2) / self.gamma_c * c**2
        
        # Thermal energy
        u = P / ((self.gamma_c - 1.0) * n)
        
        # Enthalpy in rarefaction state
        H = 1.0 + u/(c**2) + P/(n * c**2)
        
        # Tangential velocities from conservation laws (Pons et al. eqs. 27-28)
        # h*W*v^y = const and h*W*v^z = const across rarefaction
        # Equivalently: h*W*v^t = const where v^t = sqrt((v^y)² + (v^z)²)
        # Direction v^y/v^z is preserved
        
        if vta > 1e-15:
            # Initial guess for tangential velocity magnitude
            vt_mag = min(vta * (Ha * gammaa) / H, 0.9 * c)  # Start with safe value
            
            # Iterative solution for self-consistent tangential velocity
            # Needed because Lorentz factor W depends on total velocity including v^t
            for iteration in range(50):
                v2_post = v**2 + vt_mag**2
                if v2_post >= c**2 * 0.9999:
                    # Velocity approaching speed of light - reduce tangential component
                    vt_mag *= 0.90
                    continue
                
                # Lorentz factor with current total velocity
                gamma_post = 1.0 / np.sqrt(1.0 - v2_post/(c**2))
                
                # Apply conservation: (h*W*v^t)_a = (h*W*v^t)_post
                vt_new = (Ha * gammaa * vta) / (H * gamma_post)
                
                # Enforce causality: total velocity must be < c
                v2_new = v**2 + vt_new**2
                if v2_new >= c**2 * 0.9999:
                    vt_new = np.sqrt(max(0, c**2 * 0.9999 - v**2))
                
                # Check convergence
                if abs(vt_new - vt_mag) < 1e-10:
                    break
                
                # Under-relaxation for numerical stability
                vt_mag = 0.3 * vt_new + 0.7 * vt_mag
            
            # Decompose into y and z components (direction preserved)
            if vta > 1e-15:
                vy = (vya / vta) * vt_mag
                vz = (vza / vta) * vt_mag
            else:
                vy = 0.0
                vz = 0.0
        else:
            vy = 0.0
            vz = 0.0
        
        return n, P, u, v, vy, vz
        
    def solve(self, t, x0=0.5, n_points=400, tol=0.0):
        """
        Solve Riemann problem at time t.
        
        Args:
            t: Time
            x0: Initial discontinuity position
            n_points: Number of grid points
            tol: Tolerance for pressure solver
            
        Returns:
            (x, P, n, N, v, vy, vz, u, gamma, S, e) arrays
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
            vta_l = np.sqrt(self.vyl**2 + self.vzl**2)
            v2a_l = self.vl**2 + vta_l**2
            gamma_a_l = 1.0 / np.sqrt(1.0 - v2a_l/(self.c**2))
            
            x1 = x0 + (self.vl - self.csl) / (1.0 - self.vl*self.csl/self.c**2) * t
            
            vts_l = np.sqrt(self.vyls**2 + self.vzls**2)
            v2s_l = vs**2 + vts_l**2
            gamma_ls = 1.0 / np.sqrt(1.0 - v2s_l/(self.c**2))
            x2 = x0 + (vs - self.csls) / (1.0 - vs*self.csls/self.c**2) * t
        else:
            # Left shock
            x1 = x0 + self.vshockl * t
            x2 = x1
        
        x3 = x0 + vs * t
        
        if self.Pr >= Ps:
            # Right rarefaction
            vts_r = np.sqrt(self.vyrs**2 + self.vzrs**2)
            v2s_r = vs**2 + vts_r**2
            gamma_rs = 1.0 / np.sqrt(1.0 - v2s_r/(self.c**2))
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
        vy = np.zeros(n_points)
        vz = np.zeros(n_points)
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
                vy[i] = self.vyl
                vz[i] = self.vzl
                u[i] = self.Pl / ((self.gamma_c - 1.0) * self.nl)
                gamma_array[i] = self.gammal
                N[i] = self.Nl
                
            elif x[i] <= x2:
                # Left rarefaction fan
                xi = (x[i] - x0) / t
                n[i], P[i], u[i], v[i], vy[i], vz[i] = self.rarefaction(
                    xi, self.nl, self.Pl, self.Pl/((self.gamma_c-1)*self.nl),
                    self.csl, self.vl, self.vyl, self.vzl, 'L')
                v2_tot = v[i]**2 + vy[i]**2 + vz[i]**2
                gamma_array[i] = 1.0 / np.sqrt(1.0 - v2_tot/(self.c**2))
                N[i] = gamma_array[i] * n[i]
                
            elif x[i] <= x3:
                # Left star state
                P[i] = Ps
                n[i] = self.nls
                v[i] = vs
                vy[i] = self.vyls
                vz[i] = self.vzls
                u[i] = Ps / ((self.gamma_c - 1.0) * self.nls)
                v2_tot = vs**2 + self.vyls**2 + self.vzls**2
                gamma_array[i] = 1.0 / np.sqrt(1.0 - v2_tot/(self.c**2))
                N[i] = gamma_array[i] * n[i]
                
            elif x[i] <= x4:
                # Right star state
                P[i] = Ps
                n[i] = self.nrs
                v[i] = vs
                vy[i] = self.vyrs
                vz[i] = self.vzrs
                u[i] = Ps / ((self.gamma_c - 1.0) * self.nrs)
                v2_tot = vs**2 + self.vyrs**2 + self.vzrs**2
                gamma_array[i] = 1.0 / np.sqrt(1.0 - v2_tot/(self.c**2))
                N[i] = gamma_array[i] * n[i]
                
            elif x[i] <= x5:
                # Right rarefaction fan
                xi = (x[i] - x0) / t
                n[i], P[i], u[i], v[i], vy[i], vz[i] = self.rarefaction(
                    xi, self.nr, self.Pr, self.Pr/((self.gamma_c-1)*self.nr),
                    self.csr, self.vr, self.vyr, self.vzr, 'R')
                v2_tot = v[i]**2 + vy[i]**2 + vz[i]**2
                gamma_array[i] = 1.0 / np.sqrt(1.0 - v2_tot/(self.c**2))
                N[i] = gamma_array[i] * n[i]
                
            else:
                # Right state
                P[i] = self.Pr
                n[i] = self.nr
                v[i] = self.vr
                vy[i] = self.vyr
                vz[i] = self.vzr
                u[i] = self.Pr / ((self.gamma_c - 1.0) * self.nr)
                gamma_array[i] = self.gammar
                N[i] = self.Nr
            
            # Compute canonical momentum and energy
            # S = γHv (total momentum magnitude)
            H = 1.0 + u[i]/(self.c**2) + P[i]/(n[i]*self.c**2)
            vtot2 = v[i]**2 + vy[i]**2 + vz[i]**2
            vtot = np.sqrt(vtot2)
            S[i] = gamma_array[i] * H * vtot
            # e = γH - P/(Nc²)
            e[i] = gamma_array[i] * H - P[i]/(N[i]*self.c**2)
        
        return x, P, n, N, v, vy, vz, u, gamma_array, S, e
        
    def save_solution(self, filename, x, P, n, N, v, vy, vz, u, gamma, S, e, t):
        """Save solution to file."""
        np.savetxt(filename, 
                   np.column_stack([x, P, n, N, v, vy, vz, u, gamma, S, e]),
                   header=f'Time: {t}\nColumns: x, P, n, N, v, vy, vz, u, gamma, S, e')
