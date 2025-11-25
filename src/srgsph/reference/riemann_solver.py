import numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp

class State:
    def __init__(self, p, rho, vx, vy, vz, gamma_ad=5.0/3.0):
        self.p = p
        self.rho = rho
        self.vx = vx
        self.vy = vy
        self.vz = vz
        self.gamma_ad = gamma_ad
        
        self.v2 = vx**2 + vy**2 + vz**2
        if self.v2 >= 0.99999999:
            # Cap velocity to avoid NaN
            factor = np.sqrt(0.99999999 / self.v2)
            self.vx *= factor
            self.vy *= factor
            self.vz *= factor
            self.v2 = 0.99999999
            
        self.W = 1.0 / np.sqrt(1.0 - self.v2)
        # h = 1 + epsilon + p/rho
        # p = (gamma - 1) * rho * epsilon => epsilon = p / ((gamma-1)*rho)
        # h = 1 + p/((gamma-1)*rho) + p/rho = 1 + p/rho * (1/(gamma-1) + 1) = 1 + p/rho * (gamma/(gamma-1))
        self.h = 1.0 + (gamma_ad / (gamma_ad - 1.0)) * (p / rho)
        self.cs2 = (gamma_ad * p) / (self.rho * self.h) # Sound speed squared
        self.cs = np.sqrt(self.cs2)

def solve_shock(p_target, state, direction):
    # direction: -1 for Left wave (moving left), +1 for Right wave (moving right)
    # But in Riemann problem context:
    # Left state connects to L* via Left Wave (could be shock or rarefaction)
    # Right state connects to R* via Right Wave
    # For Left Wave (from L to L*), wave moves Left relative to fluid? 
    # Actually, standard convention:
    # Left Wave propagates to the left (relative to contact).
    # Right Wave propagates to the right.
    
    # Shock relations from Pons et al. 2000
    # Eq 27: Quadratic for h_b
    # h_b^2 (1 + A) - A h_b + B - h_a^2 = 0
    # A = (gamma-1)(p_a - p_b) / (gamma * p_b)
    # B = h_a (p_a - p_b) / rho_a
    
    # Note: p_a is state.p, p_b is p_target
    
    if p_target <= 0:
        return None, None # Invalid pressure

    gamma = state.gamma_ad
    pa = max(state.p, 1e-10)
    pb = p_target
    ha = state.h
    rhoa = state.rho
    
    A = (gamma - 1.0) * (pa - pb) / (gamma * pb)
    B = ha * (pa - pb) / rhoa
    
    # Quadratic: a x^2 + b x + c = 0
    # a = 1 + A
    # b = -A
    # c = B - ha**2
    
    qa = 1.0 + A
    qb = -A
    qc = B - ha**2
    
    delta = qb**2 - 4*qa*qc
    if delta < 0:
        return None, None # Should not happen for physical shocks
        
    hb = (-qb + np.sqrt(delta)) / (2*qa)
    
    # Recover rho_b
    # h = 1 + gamma/(gamma-1) * p/rho
    # h - 1 = ...
    # rho = gamma/(gamma-1) * p / (h-1)
    rhob = (gamma / (gamma - 1.0)) * pb / (hb - 1.0)
    
    # Mass flux j^2
    # Eq 28: j^2 = - (pb - pa) / (hb/rhob - ha/rhoa)
    denom = (hb/rhob - ha/rhoa)
    if abs(denom) < 1e-15:
        # Weak shock limit or error
        # For weak shock, j^2 approx (rho u)^2?
        # Use acoustic approximation?
        # j = rho W cs?
        # Let's just return None to force fallback or handle gracefully
        # But for p > pa, denom should be non-zero.
        # If pb -> pa, denom -> 0.
        # Let's use a small epsilon.
        j2 = 0.0
    else:
        j2 = -(pb - pa) / denom
        
    if j2 <= 0:
        # Unphysical or weak shock limit
        # If pb approx pa, j2 might be noise.
        # Assume vxb = vxa
        return state.vx, (rhob, hb, state.vy, state.vz)
        
    j = np.sqrt(j2)
    
    Da = rhoa * state.W
    
    # Term under sqrt in Eq 26
    term = j**2 + Da**2 * (1.0 - state.vx**2)
    sqrt_term = np.sqrt(term)
    
    denom_vs = Da**2 + j**2
    
    if direction == 'left':
        # Left wave (S_left)
        Vs = (Da**2 * state.vx - j * sqrt_term) / denom_vs
    else:
        # Right wave (S_right)
        Vs = (Da**2 * state.vx + j * sqrt_term) / denom_vs
        
    Ws = 1.0 / np.sqrt(1.0 - Vs**2)
    
    j_signed = j if direction == 'right' else -j
    
    if abs(j_signed) < 1e-15:
         return state.vx, (rhob, hb, state.vy, state.vz)

    num = ha * state.W * state.vx + Ws * (pb - pa) / j_signed
    den = ha * state.W + (pb - pa) * (Ws * state.vx / j_signed + 1.0 / Da)
    
    if abs(den) < 1e-15:
        return None, None

    vxb = num / den
    
    # Tangential velocity
    # Eq 25
    # vt_b = ha Wa vt_a * sqrt( (1 - vxb^2) / (hb^2 + (ha Wa vt_a)^2) )
    vta = np.sqrt(state.vy**2 + state.vz**2)
    if vta > 0:
        factor = np.sqrt((1.0 - vxb**2) / (hb**2 + (ha * state.W * vta)**2))
        vyb = state.vy * ha * state.W * factor / vta # Scaling original components?
        # Eq 25 says v^{y,z}_b = ... so yes, components scale.
        # Wait, Eq 25 is for components.
        vyb = state.vy * ha * state.W * np.sqrt((1.0 - vxb**2) / (hb**2 + (ha * state.W * state.vy)**2 + (ha * state.W * state.vz)**2 )) # Wait, denominator has vt_a
        # The denominator is h_b^2 + (h_a W_a v_a^t)^2.
        # So the factor is common for both y and z.
        vyb = state.vy * ha * state.W * np.sqrt((1.0 - vxb**2) / (hb**2 + (ha * state.W * vta)**2))
        vzb = state.vz * ha * state.W * np.sqrt((1.0 - vxb**2) / (hb**2 + (ha * state.W * vta)**2))
    else:
        vyb = 0.0
        vzb = 0.0
        
    return vxb, (rhob, hb, vyb, vzb)

def solve_rarefaction_analytical(p_target, state, direction):
    # Analytical solution for vt=0 (Marti & Muller 1994)
    gamma = state.gamma_ad
    pa = max(state.p, 1e-10)
    rhoa = state.rho
    
    # Isentropic relation
    const_entropy = pa / (rhoa**gamma)
    rhob = (p_target / const_entropy)**(1.0/gamma)
    
    # Sound speeds
    # csa
    ua = pa / ((gamma - 1.0) * rhoa)
    ha = 1.0 + ua + pa/rhoa
    csa = np.sqrt(gamma * pa / (rhoa * ha))
    
    # csb
    ub = p_target / ((gamma - 1.0) * rhob)
    hb = 1.0 + ub + p_target/rhob
    csb = np.sqrt(gamma * p_target / (rhob * hb))
    
    sign = 1.0 if direction == 'left' else -1.0
    
    sqrt_g_minus_1 = np.sqrt(gamma - 1.0)
    
    term_v = (1.0 + state.vx) / (1.0 - state.vx)
    term_ca = (sqrt_g_minus_1 + csa) / (sqrt_g_minus_1 - csa)
    term_cb = (sqrt_g_minus_1 + csb) / (sqrt_g_minus_1 - csb)
    
    # Exponent sign:
    # Left wave (direction='left'): v increases as p decreases.
    # term_ca / term_cb > 1 (since csa > csb).
    # We need A > term_v (v > va).
    # So exponent must be positive.
    # sign = 1.0 for left.
    # So exponent = sign * 2.0 / ...
    
    exponent = sign * 2.0 / sqrt_g_minus_1
    
    base = term_ca * (1.0 / term_cb)
    
    A = term_v * (base ** exponent)
    
    vxb = (A - 1.0) / (A + 1.0)
    
    return vxb, (rhob, hb, 0.0, 0.0)

def solve_rarefaction(p_target, state, direction):
    # Check for tangential velocity
    vta = np.sqrt(state.vy**2 + state.vz**2)
    if vta < 1e-6:
        return solve_rarefaction_analytical(p_target, state, direction)

    # direction: 'left' or 'right'
    # ODE: dvx/dp = +/- 1/(rho h W^2 cs) * 1/sqrt(1+g)
    # Sign: 
    # Eq 17: dvx/dp = +/- ...
    # For Left wave (L -> L*), pressure decreases? Or increases?
    # Rarefaction: p_b <= p_a.
    # If we go from L to L* (p* < pL), we integrate from pL to p*.
    # Which sign?
    # Standard hydro: Left rarefaction, v decreases as p increases? 
    # Left rarefaction: head moves at v - c, tail at v* - c*.
    # v increases across left rarefaction (v* > vL usually).
    # Since dp < 0, we need dvx > 0 => dvx/dp < 0.
    # So for Left wave, sign is MINUS.
    # Right wave: v decreases (v* < vR). dp < 0 => dvx < 0 => dvx/dp > 0.
    # So for Right wave, sign is PLUS.
    
    # Wait, Pons Eq 17 says +/-.
    # Let's verify with Eq 30 in Paper I (vt=0): W^2 dvx = +/- cs/rho drho = +/- cs/(gamma p) dp
    # Left wave: J- invariant? v + 2c/(g-1)?
    # In relativistic:
    # Left wave: dx/dt = (v - cs)/(1 - v cs).
    # Across Left rarefaction, Riemann invariant is constant.
    # Usually J+ is const across C-, J- across C+.
    # Let's stick to the sign deduction:
    # Left wave: v increases as p decreases. dvx/dp < 0. -> Minus.
    # Right wave: v decreases as p decreases. dvx/dp > 0. -> Plus.
    
    sign = -1.0 if direction == 'left' else 1.0
    
    # We need to integrate.
    # State variables change along the curve.
    # We need rho(p) along the adiabat.
    # Adiabatic relation: p / rho^gamma = const
    # rho = rho_a * (p / p_a)^(1/gamma)
    
    pa = max(state.p, 1e-10)
    rhoa = state.rho
    gamma = state.gamma_ad
    const_entropy = pa / (rhoa**gamma)
    
    # Invariants for tangential velocity:
    # h W vy = const, h W vz = const
    # Let K_y = h_a W_a vy_a
    # Let K_z = h_a W_a vz_a
    # Let K_t = sqrt(K_y^2 + K_z^2) = h W vt
    
    vta = np.sqrt(state.vy**2 + state.vz**2)
    Kt = state.h * state.W * vta
    
    def deriv(p, y):
        # y = [vx]
        vx = y[0]
        
        # Calculate local state
        rho = (p / const_entropy)**(1.0/gamma)
        h = 1.0 + (gamma / (gamma - 1.0)) * (p / rho)
        
        # Calculate vt from invariant
        # h W vt = Kt => vt = Kt / (h W)
        # W = 1/sqrt(1 - vx^2 - vt^2)
        # W^2 (1 - vx^2 - vt^2) = 1
        # W^2 (1 - vx^2) - (W vt)^2 = 1
        # W^2 (1 - vx^2) - (Kt/h)^2 = 1
        # W^2 = (1 + (Kt/h)^2) / (1 - vx^2)
        # W = sqrt(...)
        
        W2 = (1.0 + (Kt/h)**2) / (1.0 - vx**2)
        W = np.sqrt(W2)
        
        vt = Kt / (h * W)
        
        cs2 = (gamma * p) / (rho * h)
        cs = np.sqrt(cs2)
        
        # Calculate xi (similarity variable)
        # Eq 14 in Pons
        # xi = (vx(1-cs^2) +/- cs sqrt((1-v2)(1-v2 cs^2 - vx^2(1-cs^2)))) / (1 - v2 cs^2)
        # Sign in xi:
        # Left wave: minus
        # Right wave: plus
        
        v2 = vx**2 + vt**2
        
        term1 = vx * (1.0 - cs2)
        term2_inner = (1.0 - v2) * (1.0 - v2*cs2 - vx**2 * (1.0 - cs2))
        if term2_inner < 0: term2_inner = 0 # Numerical safety
        term2 = cs * np.sqrt(term2_inner)
        denom = 1.0 - v2 * cs2
        
        xi_sign = -1.0 if direction == 'left' else 1.0
        xi = (term1 + xi_sign * term2) / denom
        
        # g function
        # g = vt^2 (xi^2 - 1) / (1 - xi vx)^2
        g = (vt**2 * (xi**2 - 1.0)) / ((1.0 - xi * vx)**2)
        
        # ODE
        # dvx/dp = sign * 1/(rho h W^2 cs) * 1/sqrt(1+g)
        dvdP = sign * (1.0 / (rho * h * W2 * cs)) * (1.0 / np.sqrt(1.0 + g))
        
        return [dvdP]

    # Integrate from pa to p_target
    res = solve_ivp(deriv, [pa, p_target], [state.vx], rtol=1e-6, atol=1e-8)
    
    vxb = res.y[0][-1]
    
    # Recover other variables
    pb = p_target
    rhob = (pb / const_entropy)**(1.0/gamma)
    hb = 1.0 + (gamma / (gamma - 1.0)) * (pb / rhob)
    
    W2 = (1.0 + (Kt/hb)**2) / (1.0 - vxb**2)
    Wb = np.sqrt(W2)
    vtb = Kt / (hb * Wb)
    
    if vta > 0:
        vyb = state.vy * (vtb / vta) # Direction preserved
        vzb = state.vz * (vtb / vta)
    else:
        vyb = 0.0
        vzb = 0.0
        
    return vxb, (rhob, hb, vyb, vzb)

def wave_curve(p, state, direction):
    if p > state.p:
        res = solve_shock(p, state, direction)
        if res[0] is None:
            return np.nan
        return res[0]
    else:
        res = solve_rarefaction(p, state, direction)
        if res[0] is None:
            return np.nan
        return res[0]

def solve_riemann(left_state, right_state):
    # Optimization for identical states
    if abs(left_state.p - right_state.p) < 1e-6 and \
       abs(left_state.rho - right_state.rho) < 1e-6 and \
       abs(left_state.vx - right_state.vx) < 1e-6 and \
       abs(left_state.vy - right_state.vy) < 1e-6 and \
       abs(left_state.vz - right_state.vz) < 1e-6:
        return {
            'p': left_state.p,
            'rho': left_state.rho,
            'vx': left_state.vx,
            'vy': left_state.vy,
            'vz': left_state.vz,
            'h': left_state.h
        }

    # Find p_star such that v_left(p_star) = v_right(p_star)
    
    def func(p):
        vl = wave_curve(p, left_state, 'left')
        vr = wave_curve(p, right_state, 'right')
        if np.isnan(vl) or np.isnan(vr):
            return np.nan
        return vl - vr
        
    # Bracket the root
    p_min = min(left_state.p, right_state.p) * 0.001
    p_max = max(left_state.p, right_state.p) * 1000.0
    
    # Check signs
    try:
        f_min = func(p_min)
        f_max = func(p_max)
    except Exception as e:
        print(f"Error in Riemann solver bracketing: {e}")
        return None
        
    if np.isnan(f_min) or np.isnan(f_max):
        # print(f"NaN in bracket endpoints: f_min={f_min}, f_max={f_max}, p_min={p_min}, p_max={p_max}")
        return None

    if f_min * f_max > 0:
        # Try to expand bracket
        p_min *= 0.01
        p_max *= 100.0
        f_min = func(p_min)
        f_max = func(p_max)
        if np.isnan(f_min) or np.isnan(f_max) or f_min * f_max > 0:
            # Vacuum or error
            # print(f"Cannot bracket root: f_min={f_min}, f_max={f_max}, p_min={p_min}, p_max={p_max}")
            # print(f"Left: p={left_state.p}, rho={left_state.rho}, vx={left_state.vx}")
            # print(f"Right: p={right_state.p}, rho={right_state.rho}, vx={right_state.vx}")
            return None

    try:
        p_star = brentq(func, p_min, p_max)
    except ValueError:
        # print("Brentq failed")
        return None
    
    # Compute star states
    vl_star, (rhol, hl, vyl, vzl) = (None, (None, None, None, None))
    vr_star, (rhor, hr, vyr, vzr) = (None, (None, None, None, None))
    
    if p_star > left_state.p:
        vl_star, (rhol, hl, vyl, vzl) = solve_shock(p_star, left_state, 'left')
    else:
        vl_star, (rhol, hl, vyl, vzl) = solve_rarefaction(p_star, left_state, 'left')
        
    if p_star > right_state.p:
        vr_star, (rhor, hr, vyr, vzr) = solve_shock(p_star, right_state, 'right')
    else:
        vr_star, (rhor, hr, vyr, vzr) = solve_rarefaction(p_star, right_state, 'right')
        
    # Return the interface state
    # If v_star > 0, use Left Star state
    # If v_star < 0, use Right Star state
    
    v_star = 0.5 * (vl_star + vr_star)
    
    if v_star > 0:
        return {
            'p': p_star,
            'rho': rhol,
            'vx': vl_star,
            'vy': vyl,
            'vz': vzl,
            'h': hl
        }
    else:
        return {
            'p': p_star,
            'rho': rhor,
            'vx': vr_star,
            'vy': vyr,
            'vz': vzr,
            'h': hr
        }

