"""
Relativistic Riemann Problem Solver with Tangential Velocity

Based on Pons, Marti, Muller (2000) "The exact solution of the Riemann problem 
with non-zero tangential velocities in relativistic hydrodynamics"
J. Fluid Mech. (2000), vol. 422, pp. 125-139.

Key physics:
- Tangential velocity v^t affects the solution through the Lorentz factor
- The quantity K = h * W * v^t is conserved across waves (rarefactions and shocks)
- For rarefactions: integrate Eq. (12) dv^x/dp = ± 1/(rho*h*W^2*cs) * 1/sqrt(1+g)
  where g = (v^t)^2 * (xi^2 - 1) / (1 - xi * v^x)^2
- For shocks: use Rankine-Hugoniot with Eq. (25) for tangent velocity
"""

import numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp


class RelativisticRiemannSolverTangent:
    """
    Relativistic Riemann solver with tangential velocities.
    Based on Pons et al. (2000).
    
    The key insight is that the tangent velocity modifies the ODE for
    velocity within rarefactions through the factor 1/sqrt(1+g) where
    g depends on (v^t)^2.
    """
    
    def __init__(self, gamma):
        """Initialize with adiabatic index."""
        self.gamma = gamma
        self.gm1 = gamma - 1.0
        
    def compute_enthalpy(self, P, rho):
        """Compute specific enthalpy h = 1 + eps + P/rho."""
        eps = P / (self.gm1 * rho)
        return 1.0 + eps + P / rho
    
    def compute_sound_speed(self, P, rho, h):
        """Compute relativistic sound speed."""
        cs2 = self.gamma * P / (rho * h)
        return np.sqrt(max(cs2, 1e-20))
    
    def characteristic_speed_pons(self, vx, vt, cs, sign):
        """
        Compute characteristic velocity from Eq (6) of Pons et al. (2000).
        
        xi = [v^x(1-cs^2) ± cs*sqrt((1-v^2)[1 - v^2*cs^2 - (v^x)^2*(1-cs^2)])] / (1 - v^2*cs^2)
        
        Args:
            vx: normal velocity
            vt: tangent velocity (scalar, not vector)
            cs: sound speed
            sign: -1 for left-going, +1 for right-going
            
        Returns:
            xi: characteristic velocity
        """
        v2 = vx * vx + vt * vt
        cs2 = cs * cs
        
        numerator_1 = vx * (1.0 - cs2)
        
        sqrt_arg = (1.0 - v2) * (1.0 - v2 * cs2 - vx * vx * (1.0 - cs2))
        if sqrt_arg < 0:
            # Fallback to simple formula if sqrt argument is negative
            if sign < 0:
                return (vx - cs) / (1.0 - vx * cs)
            else:
                return (vx + cs) / (1.0 + vx * cs)
        
        sqrt_term = cs * np.sqrt(sqrt_arg)
        
        denom = 1.0 - v2 * cs2
        if abs(denom) < 1e-15:
            if sign < 0:
                return (vx - cs) / (1.0 - vx * cs)
            else:
                return (vx + cs) / (1.0 + vx * cs)
        
        # sign = -1 for left-going (minus in Eq 6), +1 for right-going (plus)
        return (numerator_1 + sign * sqrt_term) / denom
    
    def compute_lorentz_factor(self, vx, vt):
        """Compute Lorentz factor W = 1/sqrt(1 - vx^2 - vt^2)."""
        v2 = vx * vx + vt * vt
        if v2 >= 1.0:
            return 1000.0
        return 1.0 / np.sqrt(1.0 - v2)
    
    def compute_vt_from_K(self, K, h, vx):
        """
        Compute tangential velocity from conserved K = h * W * v^t.
        
        From Eq. (25) in Pons:
        v^t = K * sqrt((1 - v^x^2) / (h^2 + K^2))
        """
        if abs(K) < 1e-15:
            return 0.0
        
        vx2 = min(vx * vx, 0.9999)
        denom = h * h + K * K
        
        if denom <= 0:
            return 0.0
            
        vt = K * np.sqrt((1.0 - vx2) / denom)
        
        # Ensure subluminal
        if vx2 + vt * vt >= 0.9999:
            vt = np.sqrt(max(0, 0.9998 - vx2))
            
        return vt
    
    def get_velocity_shock(self, p_b, rho_a, p_a, h_a, vx_a, vt_a, W_a, K, sign):
        """
        Compute post-shock velocity using Rankine-Hugoniot conditions.
        
        Based on Eqs (17)-(21) from Pons et al.
        
        Returns:
            (vx_b, h_b, rho_b) post-shock state
        """
        gamma = self.gamma
        
        # Taub adiabat: solve quadratic for h_b (Eq. 17-19 from Pons)
        a = 1.0 + self.gm1 * (p_a - p_b) / (gamma * p_b)
        b = -self.gm1 * (p_a - p_b) / (gamma * p_b)
        c = h_a * (p_a - p_b) / rho_a - h_a * h_a
        
        disc = b*b - 4.0*a*c
        if disc < 0:
            h_b = h_a
        else:
            h_b = (-b + np.sqrt(disc)) / (2.0 * a)
        
        h_b = max(h_b, 1.0 + 1e-10)
        rho_b = gamma * p_b / (self.gm1 * (h_b - 1.0))
        
        # Mass flux j^2 (Eq. 21)
        denom = h_a / rho_a - h_b / rho_b
        if abs(denom) < 1e-15:
            return vx_a, h_a, rho_a
            
        j2 = (p_b - p_a) / denom
        if j2 < 0:
            j2 = 0.0
        j = sign * np.sqrt(j2)
        
        # Shock velocity (Eq. 20 from Pons / standard relativistic shock formula)
        rhoW_a = rho_a * W_a
        a_shock = j*j + rhoW_a*rhoW_a
        if a_shock < 1e-15:
            V_shock = vx_a
        else:
            # Standard formula: V_shock = (-b + sign * |j| * sqrt(1 + rho^2/j^2)) / a
            # where b = -vx * rho^2 * W^2
            # Simplifies to: (vx * rho^2 * W^2 + sign * |j| * sqrt(j^2 + rho^2*W^2)) / (j^2 + rho^2*W^2)
            sqrt_term = np.sqrt(j*j + rhoW_a*rhoW_a)
            V_shock = (rhoW_a*rhoW_a * vx_a + sign * abs(j) * sqrt_term) / a_shock
        
        V_shock = max(min(V_shock, 0.9999), -0.9999)
        W_shock = 1.0 / np.sqrt(max(1.0 - V_shock*V_shock, 1e-10))
        
        # Post-shock normal velocity (from conservation laws)
        # Using the standard relativistic shock formulas:
        # a = W_shock * (p_b - p_a) / j + h_a * W_a * vx_a
        # b = h_a * W_a + (p_b - p_a) * (W_shock * vx_a / j + 1.0 / (rho_a * W_a))
        # vx_b = a / b
        if abs(j) < 1e-15:
            vx_b = vx_a
        else:
            a_vel = W_shock * (p_b - p_a) / j + h_a * W_a * vx_a
            b_vel = h_a * W_a + (p_b - p_a) * (W_shock * vx_a / j + 1.0 / (rho_a * W_a))
            vx_b = a_vel / b_vel if abs(b_vel) > 1e-15 else vx_a
        
        # Clamp velocity
        vt_b = self.compute_vt_from_K(K, h_b, vx_b)
        vx_max = np.sqrt(max(0, 0.9999 - vt_b*vt_b))
        vx_b = max(min(vx_b, vx_max), -vx_max)
        
        return vx_b, h_b, rho_b
    
    def get_velocity_rarefaction_ode(self, p_b, rho_a, p_a, h_a, cs_a, vx_a, vt_a, W_a, K, sign):
        """
        Compute post-rarefaction velocity by integrating the ODE (Eq. 11/12).
        
        The ODE from Eq (11) is:
        rho*h*W^2*(v^x - xi)*dv^x + (1 - xi*v^x)*dp = 0
        
        So: dv^x/dp = -(1 - xi*v^x) / (rho*h*W^2*(v^x - xi))
        
        where xi is the characteristic velocity from Eq (6).
        
        The key correction from the Pons paper is that xi must be computed
        using the FULL Eq (6) formula which includes tangent velocity effects,
        not the simple (vx ± cs)/(1 ± vx*cs) formula.
        
        Returns:
            (vx_b, h_b, rho_b) post-rarefaction state
        """
        gamma = self.gamma
        
        # Isentropic constant
        k_poly = p_a / (rho_a ** gamma)
        
        def rhs(p, vx):
            """Right-hand side of the ODE dv^x/dp from Eq (11)."""
            if p <= 0:
                return 0.0
                
            # Density from isentropic relation
            rho = (p / k_poly) ** (1.0 / gamma)
            
            # Enthalpy
            eps = p / (self.gm1 * rho)
            h = 1.0 + eps + p / rho
            
            # Sound speed
            cs = np.sqrt(gamma * p / (rho * h))
            
            # Tangent velocity from K = h * W * v^t
            vt = self.compute_vt_from_K(K, h, vx)
            
            # Lorentz factor
            v2 = vx * vx + vt * vt
            if v2 >= 0.9999:
                return 0.0
            W2 = 1.0 / (1.0 - v2)
            
            # Use the CORRECT characteristic velocity from Eq (6)!
            # This is the key fix - the characteristic speed depends on vt
            xi = self.characteristic_speed_pons(vx, vt, cs, sign)
            
            # From Eq (11): dv^x/dp = -(1 - xi*v^x) / (rho*h*W^2*(v^x - xi))
            numerator = -(1.0 - xi * vx)
            denominator = rho * h * W2 * (vx - xi)
            
            if abs(denominator) < 1e-30:
                return 0.0
            
            return numerator / denominator
        
        # Integrate from p_a to p_b
        if abs(p_b - p_a) < 1e-15:
            rho_b = (p_b / k_poly) ** (1.0 / gamma)
            h_b = self.compute_enthalpy(p_b, rho_b)
            return vx_a, h_b, rho_b
        
        # Use simple Euler integration with small steps for robustness
        n_steps = 2000
        dp = (p_b - p_a) / n_steps
        vx = vx_a
        p = p_a
        
        for _ in range(n_steps):
            dvdp = rhs(p, vx)
            vx += dvdp * dp
            p += dp
            
            # Clamp velocity to ensure subluminal
            rho = (p / k_poly) ** (1.0 / gamma) if p > 0 else rho_a
            h = self.compute_enthalpy(max(p, 1e-10), max(rho, 1e-10))
            vt = self.compute_vt_from_K(K, h, vx)
            vx_max = np.sqrt(max(0, 0.9999 - vt * vt))
            vx = max(min(vx, vx_max), -vx_max)
        
        # Final state
        rho_b = (p_b / k_poly) ** (1.0 / gamma)
        h_b = self.compute_enthalpy(p_b, rho_b)
        
        # Clamp velocity
        vt_b = self.compute_vt_from_K(K, h_b, vx)
        vx_max = np.sqrt(max(0, 0.9999 - vt_b * vt_b))
        vx = max(min(vx, vx_max), -vx_max)
        
        return vx, h_b, rho_b
    
    def velocity_curve(self, p, state, sign):
        """
        Compute velocity behind a wave as function of pressure.
        
        Returns:
            vx_b: Post-wave normal velocity
        """
        rho_a, p_a, h_a, cs_a, vx_a, vt_a, W_a, K = state
        
        if p > p_a:
            # Shock
            vx_b, _, _ = self.get_velocity_shock(p, rho_a, p_a, h_a, vx_a, vt_a, W_a, K, sign)
        else:
            # Rarefaction
            vx_b, _, _ = self.get_velocity_rarefaction_ode(p, rho_a, p_a, h_a, cs_a, vx_a, vt_a, W_a, K, sign)
        
        return vx_b
    
    def velocity_difference(self, p, state_L, state_R):
        """Compute difference in post-wave velocities."""
        vx_L = self.velocity_curve(p, state_L, -1)
        vx_R = self.velocity_curve(p, state_R, +1)
        return vx_L - vx_R
    
    def solve(self, P_L, rho_L, vx_L, vt_L, P_R, rho_R, vx_R, vt_R, t, n_points=500):
        """
        Solve the Riemann problem with tangential velocities.
        
        Returns:
            Dictionary with x, pressure, density, vel_x, vel_t arrays
        """
        gamma = self.gamma
        
        # Compute derived quantities
        h_L = self.compute_enthalpy(P_L, rho_L)
        h_R = self.compute_enthalpy(P_R, rho_R)
        W_L = self.compute_lorentz_factor(vx_L, vt_L)
        W_R = self.compute_lorentz_factor(vx_R, vt_R)
        cs_L = self.compute_sound_speed(P_L, rho_L, h_L)
        cs_R = self.compute_sound_speed(P_R, rho_R, h_R)
        
        # Conserved tangent momentum K = h * W * v^t
        K_L = h_L * W_L * vt_L
        K_R = h_R * W_R * vt_R
        
        state_L = (rho_L, P_L, h_L, cs_L, vx_L, vt_L, W_L, K_L)
        state_R = (rho_R, P_R, h_R, cs_R, vx_R, vt_R, W_R, K_R)
        
        # Find P_s using bisection
        # For relativistic problems, velocity curves can saturate at high pressure,
        # so we need sensible initial bounds
        P_min = min(P_L, P_R) * 0.01
        P_max = max(P_L, P_R) * 2.0
        
        # Ensure we have a valid bracket
        dvel_min = self.velocity_difference(P_min, state_L, state_R)
        dvel_max = self.velocity_difference(P_max, state_L, state_R)
        
        # Adjust bounds to ensure bracketing
        for _ in range(30):
            if dvel_min * dvel_max < 0 and abs(dvel_min) > 1e-10 and abs(dvel_max) > 1e-10:
                break
            # Expand in appropriate direction
            if dvel_min > 0 and dvel_max >= 0:
                # Need to find where it goes negative - increase P_max
                P_max *= 3.0
            elif dvel_min <= 0 and dvel_max < 0:
                # Need to find where it goes positive - decrease P_min
                P_min *= 0.3
            elif abs(dvel_max) < 1e-10:
                # dvel_max is ~0, reduce P_max
                P_max *= 0.5
            elif abs(dvel_min) < 1e-10:
                # dvel_min is ~0, increase P_min
                P_min *= 2.0
            else:
                # Try shrinking from both sides
                P_min *= 1.5
                P_max *= 0.7
            
            dvel_min = self.velocity_difference(P_min, state_L, state_R)
            dvel_max = self.velocity_difference(P_max, state_L, state_R)
        
        # Find pressure using bisection
        try:
            P_s = brentq(lambda P: self.velocity_difference(P, state_L, state_R), 
                        P_min, P_max, xtol=1e-10)
        except ValueError:
            # Fallback: try geometric mean as starting point
            P_s = np.sqrt(P_L * P_R)
        
        # Compute star state velocities and densities
        vx_Ls = self.velocity_curve(P_s, state_L, -1)
        vx_Rs = self.velocity_curve(P_s, state_R, +1)
        vx_star = 0.5 * (vx_Ls + vx_Rs)
        
        if P_s > P_L:
            _, h_Ls, rho_Ls = self.get_velocity_shock(P_s, rho_L, P_L, h_L, vx_L, vt_L, W_L, K_L, -1)
            left_is_shock = True
        else:
            _, h_Ls, rho_Ls = self.get_velocity_rarefaction_ode(P_s, rho_L, P_L, h_L, cs_L, vx_L, vt_L, W_L, K_L, -1)
            left_is_shock = False
        
        if P_s > P_R:
            _, h_Rs, rho_Rs = self.get_velocity_shock(P_s, rho_R, P_R, h_R, vx_R, vt_R, W_R, K_R, +1)
            right_is_shock = True
        else:
            _, h_Rs, rho_Rs = self.get_velocity_rarefaction_ode(P_s, rho_R, P_R, h_R, cs_R, vx_R, vt_R, W_R, K_R, +1)
            right_is_shock = False
        
        # Compute tangent velocities in star states
        vt_Ls = self.compute_vt_from_K(K_L, h_Ls, vx_star)
        vt_Rs = self.compute_vt_from_K(K_R, h_Rs, vx_star)
        
        # Sound speeds in star states
        cs_Ls = self.compute_sound_speed(P_s, rho_Ls, h_Ls)
        cs_Rs = self.compute_sound_speed(P_s, rho_Rs, h_Rs)
        
        # Wave positions
        x_min, x_max = -0.5, 0.5
        x = np.linspace(x_min, x_max, n_points)
        
        if t <= 0:
            pres = np.where(x < 0, P_L, P_R)
            dens = np.where(x < 0, rho_L, rho_R)
            vel_x = np.where(x < 0, vx_L, vx_R)
            vel_t = np.where(x < 0, vt_L, vt_R)
            return {'x': x, 'pres': pres, 'dens': dens, 'vel_x': vel_x, 'vel_t': vel_t,
                    'P_star': P_s, 'vx_star': vx_star}
        
        # Left wave positions
        if left_is_shock:
            j2 = (P_s - P_L) / max(abs(h_L / rho_L - h_Ls / rho_Ls), 1e-15)
            j = -np.sqrt(max(j2, 0))
            rhoW_L = rho_L * W_L
            a_shock = j*j + rhoW_L*rhoW_L
            V_shock_L = (rhoW_L*rhoW_L * vx_L + j*j - abs(j) * np.sqrt(j*j + rhoW_L*rhoW_L)) / max(a_shock, 1e-15)
            x1 = V_shock_L * t
            x2 = x1
        else:
            x1 = (vx_L - cs_L) / (1.0 - vx_L * cs_L) * t
            x2 = (vx_star - cs_Ls) / (1.0 - vx_star * cs_Ls) * t
        
        x3 = vx_star * t
        
        if right_is_shock:
            j2 = (P_s - P_R) / max(abs(h_R / rho_R - h_Rs / rho_Rs), 1e-15)
            j = np.sqrt(max(j2, 0))
            rhoW_R = rho_R * W_R
            a_shock = j*j + rhoW_R*rhoW_R
            V_shock_R = (rhoW_R*rhoW_R * vx_R + j*j + abs(j) * np.sqrt(j*j + rhoW_R*rhoW_R)) / max(a_shock, 1e-15)
            x4 = V_shock_R * t
            x5 = x4
        else:
            x4 = (vx_star + cs_Rs) / (1.0 + vx_star * cs_Rs) * t
            x5 = (vx_R + cs_R) / (1.0 + vx_R * cs_R) * t
        
        # Build solution
        pres = np.zeros(n_points)
        dens = np.zeros(n_points)
        vel_x = np.zeros(n_points)
        vel_t = np.zeros(n_points)
        
        for i in range(n_points):
            xi = x[i]
            
            if xi <= x1:
                pres[i], dens[i], vel_x[i], vel_t[i] = P_L, rho_L, vx_L, vt_L
            elif xi <= x2:
                if not left_is_shock:
                    frac = (xi - x1) / max(abs(x2 - x1), 1e-10)
                    pres[i] = P_L + frac * (P_s - P_L)
                    dens[i] = rho_L + frac * (rho_Ls - rho_L)
                    vel_x[i] = vx_L + frac * (vx_star - vx_L)
                    vel_t[i] = vt_L + frac * (vt_Ls - vt_L)
                else:
                    pres[i], dens[i], vel_x[i], vel_t[i] = P_s, rho_Ls, vx_star, vt_Ls
            elif xi <= x3:
                pres[i], dens[i], vel_x[i], vel_t[i] = P_s, rho_Ls, vx_star, vt_Ls
            elif xi <= x4:
                pres[i], dens[i], vel_x[i], vel_t[i] = P_s, rho_Rs, vx_star, vt_Rs
            elif xi <= x5:
                if not right_is_shock:
                    frac = (xi - x4) / max(abs(x5 - x4), 1e-10)
                    pres[i] = P_s + frac * (P_R - P_s)
                    dens[i] = rho_Rs + frac * (rho_R - rho_Rs)
                    vel_x[i] = vx_star + frac * (vx_R - vx_star)
                    vel_t[i] = vt_Rs + frac * (vt_R - vt_Rs)
                else:
                    pres[i], dens[i], vel_x[i], vel_t[i] = P_R, rho_R, vx_R, vt_R
            else:
                pres[i], dens[i], vel_x[i], vel_t[i] = P_R, rho_R, vx_R, vt_R
        
        return {
            'x': x, 'pres': pres, 'dens': dens, 'vel_x': vel_x, 'vel_t': vel_t,
            'P_star': P_s, 'vx_star': vx_star,
            'rho_Ls': rho_Ls, 'rho_Rs': rho_Rs,
            'vt_Ls': vt_Ls, 'vt_Rs': vt_Rs,
            'x1': x1, 'x2': x2, 'x3': x3, 'x4': x4, 'x5': x5,
            'left_is_shock': left_is_shock, 'right_is_shock': right_is_shock
        }


def test_table1():
    """Test against Table 1 from Pons et al. (2000)."""
    print("Testing against Pons et al. (2000) Table 1")
    print("=" * 70)
    print("Initial data: P_L=1000, P_R=0.01, rho=1, vx=0, gamma=5/3")
    print()
    
    gamma = 5.0/3.0
    solver = RelativisticRiemannSolverTangent(gamma)
    
    # Test cases from Table 1
    # Format: (vt_L, vt_R, expected_p_star, expected_v_star, expected_rho_L_star, expected_rho_R_star)
    test_cases = [
        (0.0, 0.0, 18.6, 0.960, 0.0916, 10.4),
        (0.0, 0.9, 42.8, 0.913, 0.151, 14.6),
        (0.9, 0.0, 0.189, 0.328, 0.00583, 3.44),
        (0.9, 0.9, 0.904, 0.319, 0.0149, 4.46),
        (0.99, 0.99, 0.706, 0.095, 0.0129, 4.29),
    ]
    
    P_L, rho_L, vx_L = 1000.0, 1.0, 0.0
    P_R, rho_R, vx_R = 0.01, 1.0, 0.0
    
    print(f"{'vt_L':>5} {'vt_R':>5} | {'P*':>8} {'v*':>8} {'rho_L*':>10} {'rho_R*':>8} | Expected P* v*")
    print("-" * 70)
    
    for vt_L, vt_R, exp_p, exp_v, exp_rhoL, exp_rhoR in test_cases:
        result = solver.solve(P_L, rho_L, vx_L, vt_L, P_R, rho_R, vx_R, vt_R, t=0.4)
        
        wave_L = "R" if not result['left_is_shock'] else "S"
        wave_R = "S" if result['right_is_shock'] else "R"
        
        p_err = abs(result['P_star'] - exp_p) / exp_p * 100
        v_err = abs(result['vx_star'] - exp_v) / max(exp_v, 0.001) * 100
        
        status = "OK" if p_err < 5 and v_err < 5 else "CHECK"
        
        print(f"{vt_L:5.2f} {vt_R:5.2f} | {result['P_star']:8.4f} {result['vx_star']:8.4f} "
              f"{result['rho_Ls']:10.5f} {result['rho_Rs']:8.4f} | "
              f"{exp_p:6.3f} {exp_v:5.3f} [{status}]")


if __name__ == "__main__":
    test_table1()
