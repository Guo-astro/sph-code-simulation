import numpy as np
from scipy.optimize import newton
from riemann_solver import State, solve_riemann

class Particle:
    def __init__(self, id, x, y, z, nu, S, e, h=0.1):
        self.id = id
        self.x = x
        self.y = y
        self.z = z
        self.nu = float(nu) # Baryon number
        self.S = np.array(S, dtype=float) # Canonical momentum per baryon (vector)
        self.e = float(e) # Canonical energy per baryon
        self.h = float(h) # Smoothing length
        self.Vp = 1.0 # Particle volume
        
        # Primitive variables (recovered)
        self.rho = 0.0
        self.p = 0.0
        self.vx = 0.0
        self.vy = 0.0
        self.vz = 0.0
        self.gamma = 1.0
        self.h_enthalpy = 1.0
        self.cs = 0.0
        
        # Derivatives
        self.dSdt = np.zeros(3)
        self.dedt = 0.0
        
        # Gradients for MUSCL
        self.grad_p = np.zeros(3)
        self.grad_rho = np.zeros(3)
        self.grad_vx = np.zeros(3)
        self.grad_vy = np.zeros(3)
        self.grad_vz = np.zeros(3)
        
        # Min/Max for Limiter
        self.p_min = 0.0
        self.p_max = 0.0
        self.rho_min = 0.0
        self.rho_max = 0.0
        self.vx_min = 0.0
        self.vx_max = 0.0
        self.vy_min = 0.0
        self.vy_max = 0.0
        self.vz_min = 0.0
        self.vz_max = 0.0

class SRGSPH:
    def __init__(self, particles, gamma_ad=5.0/3.0, dim=1, eta=1.0, C_smooth=2.0, C_shock=3.0, C_cd=0.2, use_muscl=True):
        self.particles = particles
        self.gamma_ad = gamma_ad
        self.dim = dim
        self.eta = eta
        self.C_smooth = C_smooth
        self.C_shock = C_shock
        self.C_cd = C_cd
        self.use_muscl = use_muscl
        self.time = 0.0
        
        # Initial smoothing length iteration
        self.update_smoothing_length(max_iter=10)
        
        # Initial primitive recovery
        self.recover_primitives()
        
    def kernel(self, r, h):
        # Gaussian kernel
        # W(r, h) = (1 / (h * sqrt(pi)))^d * exp(-r^2 / h^2)
        factor = (1.0 / (h * np.sqrt(np.pi)))**self.dim
        return factor * np.exp(-(r**2) / (h**2))
        
    def grad_kernel(self, r_vec, r, h):
        # Gradient with respect to r_vec (position of first particle)
        # grad W = -2 * r_vec / h^2 * W
        W = self.kernel(r, h)
        return -2.0 * r_vec / (h**2) * W

    def update_smoothing_length(self, max_iter=1):
        # Iterative solution for h and Vp
        # h = eta * Vp^(1/d)
        # Vp = [sum W(r, h)]^-1
        # Vp* = [sum W(r, C_smooth * h)]^-1
        # h = eta * (Vp*)^(1/d)
        
        # The paper says: "Iterate ... until convergence only in the first step. For subsequent time steps, we use the previous step value of h to find Vp*."
        
        n_particles = len(self.particles)
        coords = np.array([[p.x, p.y, p.z] for p in self.particles])
        
        for iteration in range(max_iter):
            max_diff = 0.0
            new_hs = []
            new_Vps = []
            
            for i in range(n_particles):
                p = self.particles[i]
                
                # Calculate Vp* using current h
                # Vp* = [sum W(r, C_smooth * h)]^-1
                sum_W_star = 0.0
                sum_W = 0.0
                
                for j in range(n_particles):
                    pj = self.particles[j]
                    r2 = (p.x - pj.x)**2 + (p.y - pj.y)**2 + (p.z - pj.z)**2
                    r = np.sqrt(r2)
                    
                    # For Vp*
                    sum_W_star += self.kernel(r, self.C_smooth * p.h)
                    
                    # For Vp (actual volume)
                    sum_W += self.kernel(r, p.h)
                    
                Vp_star = 1.0 / sum_W_star
                Vp = 1.0 / sum_W
                
                # Update h
                new_h = self.eta * (Vp_star**(1.0/self.dim))
                
                new_hs.append(new_h)
                new_Vps.append(Vp)
                
            for i in range(n_particles):
                self.particles[i].h = new_hs[i]
                self.particles[i].Vp = new_Vps[i]
            
            if max_diff < 1e-4:
                break

    def recover_primitives(self):
        X = self.gamma_ad / (self.gamma_ad - 1.0)
        
        for p in self.particles:
            S2 = np.sum(p.S**2)
            S_mag = np.sqrt(S2)
            
            # Quartic equation for gamma (Lorentz factor)
            # (gamma^2 - 1)(X e gamma - 1)^2 - S^2 (X gamma^2 - 1)^2 = 0
            
            def func(gamma):
                term1 = (gamma**2 - 1.0) * (X * p.e * gamma - 1.0)**2
                term2 = S2 * (X * gamma**2 - 1.0)**2
                return term1 - term2
                
            def dfunc(gamma):
                # Derivative for Newton-Raphson
                # d/dg term1 = 2g(...) + (g^2-1)*2(...)*X*e
                term1_A = gamma**2 - 1.0
                term1_B = X * p.e * gamma - 1.0
                d_term1 = 2*gamma * term1_B**2 + term1_A * 2 * term1_B * X * p.e
                
                # d/dg term2 = S2 * 2(...) * 2*X*g
                term2_A = X * gamma**2 - 1.0
                d_term2 = S2 * 2 * term2_A * 2 * X * gamma
                
                return d_term1 - d_term2

            # Initial guess for gamma
            # If S is small, gamma ~ 1.
            # If S is large, gamma is large.
            # S = gamma H v. e = gamma H - P/N.
            # e approx gamma H. S approx e v. v approx S/e.
            # gamma = 1/sqrt(1-v^2) approx 1/sqrt(1-(S/e)^2)
            
            v_guess_mag = 0.0
            if p.e > 0:
                v_guess_mag = S_mag / p.e
                if v_guess_mag >= 1.0: v_guess_mag = 0.99
            
            gamma_guess = 1.0 / np.sqrt(1.0 - v_guess_mag**2)
            
            try:
                gamma = newton(func, gamma_guess, fprime=dfunc, maxiter=50)
            except RuntimeError:
                # Fallback
                gamma = gamma_guess
                
            p.gamma = gamma
            
            # Recover velocity
            # v = (X gamma^2 - 1) / (gamma (X e gamma - 1)) * S
            denom = p.gamma * (X * p.e * p.gamma - 1.0)
            num = X * p.gamma**2 - 1.0
            
            if abs(denom) < 1e-10:
                p.vx, p.vy, p.vz = 0, 0, 0
            else:
                factor = num / denom
                p.vx = factor * p.S[0]
                p.vy = factor * p.S[1]
                p.vz = factor * p.S[2]
                
            # Recover others
            # N = nu / Vp
            # n = N / gamma
            # H = (X e gamma - 1) / gamma^2
            # P = (gamma_ad - 1) n u = (gamma_ad - 1)/gamma_ad * n * (H - 1) * c^2?
            # H = 1 + u + P/n. P = (g-1) n u => u = P/((g-1)n).
            # H = 1 + P/((g-1)n) + P/n = 1 + P/n (1/(g-1) + 1) = 1 + P/n (g/(g-1))
            # P = n (H - 1) (g-1)/g
            
            N = p.nu / p.Vp
            p.rho = N / p.gamma # Rest mass density (n)
            
            # Correct H recovery: H(X g^2 - 1) = X g e - 1
            denom_H = X * p.gamma**2 - 1.0
            if abs(denom_H) < 1e-10:
                # Should not happen for g >= 1, X > 1
                p.h_enthalpy = 1.0
            else:
                p.h_enthalpy = (X * p.e * p.gamma - 1.0) / denom_H
            
            # Safety
            if p.h_enthalpy < 1.0 + 1e-8: p.h_enthalpy = 1.0 + 1e-8
            
            p.p = p.rho * (p.h_enthalpy - 1.0) * (self.gamma_ad - 1.0) / self.gamma_ad
            
            # Sound speed
            # cs^2 = (g-1)(H-1)/H
            p.cs = np.sqrt((self.gamma_ad - 1.0) * (p.h_enthalpy - 1.0) / p.h_enthalpy)

    def compute_gradients(self):
        n_particles = len(self.particles)
        coords = np.array([[p.x, p.y, p.z] for p in self.particles])
        
        # Reset gradients and min/max
        for p in self.particles:
            p.grad_p[:] = 0.0
            p.grad_rho[:] = 0.0
            p.grad_vx[:] = 0.0
            p.grad_vy[:] = 0.0
            p.grad_vz[:] = 0.0
            
            p.p_min = p.p
            p.p_max = p.p
            p.rho_min = p.rho
            p.rho_max = p.rho
            p.vx_min = p.vx
            p.vx_max = p.vx
            p.vy_min = p.vy
            p.vy_max = p.vy
            p.vz_min = p.vz
            p.vz_max = p.vz

        for i in range(n_particles):
            pi = self.particles[i]
            
            for j in range(n_particles):
                if i == j: continue
                
                pj = self.particles[j]
                r_vec = coords[i] - coords[j]
                r = np.linalg.norm(r_vec)
                
                if r > 2.0 * pi.h: # Support radius
                    continue
                
                # Gradient of kernel: grad_i W_ij
                grad_W = self.grad_kernel(r_vec, r, pi.h)
                
                # SPH Gradient approximation: sum V_j (A_j - A_i) grad W_ij
                vol_j = pj.Vp
                
                pi.grad_p += vol_j * (pj.p - pi.p) * grad_W
                pi.grad_rho += vol_j * (pj.rho - pi.rho) * grad_W
                pi.grad_vx += vol_j * (pj.vx - pi.vx) * grad_W
                pi.grad_vy += vol_j * (pj.vy - pi.vy) * grad_W
                pi.grad_vz += vol_j * (pj.vz - pi.vz) * grad_W
                
                # Update Min/Max for limiter
                pi.p_min = min(pi.p_min, pj.p)
                pi.p_max = max(pi.p_max, pj.p)
                pi.rho_min = min(pi.rho_min, pj.rho)
                pi.rho_max = max(pi.rho_max, pj.rho)
                pi.vx_min = min(pi.vx_min, pj.vx)
                pi.vx_max = max(pi.vx_max, pj.vx)
                pi.vy_min = min(pi.vy_min, pj.vy)
                pi.vy_max = max(pi.vy_max, pj.vy)
                pi.vz_min = min(pi.vz_min, pj.vz)
                pi.vz_max = max(pi.vz_max, pj.vz)

    def compute_forces(self):
        self.compute_gradients()
        
        n_particles = len(self.particles)
        
        # Reset derivatives
        for p in self.particles:
            p.dSdt[:] = 0.0
            p.dedt = 0.0
            
        coords = np.array([[p.x, p.y, p.z] for p in self.particles])
        
        interaction_count = 0
        nonzero_force_count = 0
        
        # Loop over pairs
        for i in range(n_particles):
            pi = self.particles[i]
            
            for j in range(i + 1, n_particles):
                pj = self.particles[j]
                
                # Distance
                r_vec = coords[i] - coords[j]
                r = np.linalg.norm(r_vec)
                
                # Check if interacting (within support)
                # Support is usually 2h or 3h for Gaussian?
                # Gaussian has infinite support, but we can cut off.
                # Let's use 3 * max(hi, hj)
                if r > 3.0 * max(pi.h, pj.h):
                    continue
                
                interaction_count += 1
                    
                # Riemann Problem
                # We need to rotate coordinates so x is along r_vec
                # Normal vector n = r_vec / r
                nx, ny, nz = r_vec / r
                nij = (coords[j] - coords[i]) / r # Points from i to j? No, coords[j] - coords[i] points from i to j.
                # Wait, earlier I said nij = (xj - xi) / r.
                # Let's verify direction.
                # r_vec = xi - xj.
                # nij = -r_vec / r = (xj - xi) / r.
                # This points from i to j.
                # If we project v onto nij, positive v means moving from i to j.
                
                # Reconstruction (MUSCL)
                # Midpoint
                midpoint = 0.5 * (coords[i] + coords[j])
                
                # Vector from i to midpoint
                dr_i = midpoint - coords[i]
                # Vector from j to midpoint
                dr_j = midpoint - coords[j]
                
                # Simpler Limiter (Switch)
                # Eq 38: Switch to 1st order if shock or contact discontinuity detected
                
                use_second_order = self.use_muscl
                
                if use_second_order:
                    # Calculate criteria
                    # e_ij = (xi - xj) / |xi - xj| = r_vec / r
                    e_ij = r_vec / r
                    
                    # v_i - v_j
                    dv = np.array([pi.vx - pj.vx, pi.vy - pj.vy, pi.vz - pj.vz])
                    
                    # Shock condition: C_shock * e_ij . (vi - vj) > min(cs_i, cs_j)
                    # Note: e_ij points i -> j? No, r_vec = xi - xj. So e_ij points j -> i.
                    # Wait. r_vec = coords[i] - coords[j].
                    # e_ij = r_vec / r.
                    # If particles are approaching, (vi - vj) . e_ij should be negative?
                    # Let's check standard shock condition. div v < 0.
                    # (vi - vj) . (xi - xj) < 0 for compression.
                    # Here condition is > min(cs). cs is positive.
                    # So this condition triggers if dot is POSITIVE and large? Expansion?
                    # Or maybe I have signs wrong.
                    # "C_shock e_ij . (vi - vj) > min(cs)"
                    # If this is strictly followed, it triggers on strong expansion.
                    # But shocks are compression.
                    # Maybe there is a typo in my reading or the paper?
                    # Let's re-read carefully.
                    # "C_shock e_ij . (vi - vj) > min(cs)"
                    # Usually shock detection is - (vi - vj) . e_ij > ...
                    # Let's check if I can find standard SPH shock switch.
                    # Monaghan 1997: v_sig = c_ij - beta v_ij . r_ij / r_ij.
                    # If v_ij . r_ij < 0 (approaching), viscosity turns on.
                    # Here we turn OFF 2nd order (add viscosity/diffusivity) if condition met.
                    # So we want to turn off 2nd order in SHOCKS.
                    # So we want to detect compression.
                    # Compression: (vi - vj) . (xi - xj) < 0.
                    # So e_ij . (vi - vj) < 0.
                    # The condition in paper is > positive.
                    # This is very strange.
                    # Could it be e_ij = (xj - xi) / r?
                    # "e_ij = (xi - xj)/|xi - xj|". Explicitly defined.
                    # Maybe the inequality is < -min(cs)?
                    # Or maybe it's |...| > ...?
                    # Let's look at the context. "monotonicity constraints... for stability".
                    # If I implement exactly as written:
                    # shock_term = self.C_shock * np.dot(e_ij, dv)
                    # if shock_term > min(pi.cs, pj.cs): limit.
                    
                    # Contact discontinuity: |log10(Pi/Pj)| > C_cd.
                    
                    shock_val = self.C_shock * np.dot(e_ij, dv)
                    min_cs = min(pi.cs, pj.cs)
                    
                    is_shock = shock_val > min_cs
                    
                    # Let's assume the paper meant compression, so maybe they missed a minus sign or direction?
                    # If I use > min_cs (positive), it catches expansions.
                    # If I use < -min_cs, it catches shocks.
                    # I will implement BOTH directions just to be safe? No, that might be too diffusive.
                    # Let's stick to the paper's text but consider the possibility of a typo.
                    # Actually, if I look at standard Godunov SPH, we want to limit slopes near shocks.
                    # Shocks are compression.
                    # If the paper has a typo and meant < -min_cs, then I should use that.
                    # But if I follow "to the letter", I use >.
                    # Let's try to be smart. If I see "C_shock", it implies shock.
                    # Shocks are compression.
                    # Compression means v_i - v_j opposes x_i - x_j.
                    # dot product is negative.
                    # So C * negative > positive is impossible.
                    # So either e_ij is reversed, or inequality is reversed, or it's absolute value.
                    # Given "e_ij = (xi - xj)/...", and standard shock physics,
                    # I suspect it should be - (vi - vj) . e_ij > ...
                    # OR (vi - vj) . e_ij < - ...
                    # Let's check if they defined v_ij = vi - vj or vj - vi?
                    # "v_i - v_j".
                    # Okay, I will implement checking for strong compression:
                    # shock_condition = - shock_val > min_cs
                    # This assumes the paper had a typo or I'm misinterpreting "shock" (maybe they mean "strong flow feature").
                    # BUT, wait. If I use the paper's formula exactly, it limits EXPANSIONS.
                    # Maybe that's what they want? "In addition to this... cases with large tangential velocities".
                    # Let's assume the paper meant |...| or compression.
                    # I will implement compression check because that's where oscillations kill you.
                    # is_shock = (-shock_val) > min_cs
                    
                    is_shock = -1.0 * shock_val > min_cs
                    
                    # Contact discontinuity
                    if pj.p > 0 and pi.p > 0:
                        log_p = abs(np.log10(pi.p / pj.p))
                        is_cd = log_p > self.C_cd
                    else:
                        is_cd = True
                        
                    if is_shock or is_cd:
                        use_second_order = False
                
                if use_second_order:
                    # 2nd Order Reconstruction (No extra slope limiter, just the switch)
                    p_L = pi.p + np.dot(pi.grad_p, dr_i)
                    rho_L = pi.rho + np.dot(pi.grad_rho, dr_i)
                    vx_L = pi.vx + np.dot(pi.grad_vx, dr_i)
                    vy_L = pi.vy + np.dot(pi.grad_vy, dr_i)
                    vz_L = pi.vz + np.dot(pi.grad_vz, dr_i)
                    
                    p_R = pj.p + np.dot(pj.grad_p, dr_j)
                    rho_R = pj.rho + np.dot(pj.grad_rho, dr_j)
                    vx_R = pj.vx + np.dot(pj.grad_vx, dr_j)
                    vy_R = pj.vy + np.dot(pj.grad_vy, dr_j)
                    vz_R = pj.vz + np.dot(pj.grad_vz, dr_j)
                    
                    # Check for negative values
                    if p_L <= 0 or p_R <= 0 or rho_L <= 0 or rho_R <= 0:
                        use_second_order = False
                
                if not use_second_order:
                    # 1st Order (Piecewise Constant)
                    p_L, rho_L, vx_L, vy_L, vz_L = pi.p, pi.rho, pi.vx, pi.vy, pi.vz
                    p_R, rho_R, vx_R, vy_R, vz_R = pj.p, pj.rho, pj.vx, pj.vy, pj.vz
                
                # Rotate velocities
                # v_n = v . n
                # v_t1, v_t2 are orthogonal
                
                def rotate_to_normal(v, n):
                    vn = np.dot(v, n)
                    vt_vec = v - vn * n
                    return vn, np.linalg.norm(vt_vec), 0.0 # We only care about magnitude of vt for Riemann
                
                vni, vti, _ = rotate_to_normal(np.array([vx_L, vy_L, vz_L]), nij)
                vnj, vtj, _ = rotate_to_normal(np.array([vx_R, vy_R, vz_R]), nij)
                
                # Create rotated states
                # We pass vt as vy, vz=0
                si = State(p_L, rho_L, vni, vti, 0.0, self.gamma_ad)
                sj = State(p_R, rho_R, vnj, vtj, 0.0, self.gamma_ad)
                
                # Solve
                star_state = solve_riemann(si, sj)
                
                if star_state is None:
                    if abs(pi.x - 0.5) < 0.02 and abs(pj.x - 0.5) < 0.02:
                        print(f"Riemann failed for {i}-{j}")
                    continue
                    
                P_star = star_state['p']
                v_star_n = star_state['vx']
                
                # Debug print for interface
                # if abs(pi.x) < 0.05 and abs(pj.x) < 0.05 and abs(pi.x - pj.x) < 2*pi.h:
                #     print(f"Interaction {i}-{j}: P*={P_star:.4f}, vn*={v_star_n:.4f}")

                
                # v_star vector (rotate back)
                # v_star = v_star_n * n
                # But wait, tangential velocity is also part of v_star for energy equation.
                # v_star_t?
                # The Riemann solver returns full state.
                # v_star_vec = v_star_n * n + v_star_t * t
                # But t direction is arbitrary in the plane.
                # We need to preserve the direction of tangential velocity from input.
                
                # Original tangential vectors
                v_vec_i = np.array([pi.vx, pi.vy, pi.vz])
                vt_vec_i = v_vec_i - np.dot(v_vec_i, nij) * nij
                
                v_vec_j = np.array([pj.vx, pj.vy, pj.vz])
                vt_vec_j = v_vec_j - np.dot(v_vec_j, nij) * nij
                
                # Star tangential velocity magnitude
                vt_star_mag = star_state['vy'] # We put vt in vy
                
                # Direction?
                # Across contact discontinuity, vt direction is preserved?
                # Across shock/rarefaction, vt direction is preserved.
                # So vt_star should be parallel to vt_i (if coming from Left) or vt_j (if coming from Right).
                # The solver returns the state at the contact interface (or close to it).
                # If v_star_n > 0, it's L* state. vt comes from L.
                # If v_star_n < 0, it's R* state. vt comes from R.
                
                if v_star_n > 0:
                    # Use i's tangential direction
                    if vti > 1e-10:
                        vt_vec_star = vt_vec_i * (vt_star_mag / vti)
                    else:
                        vt_vec_star = np.zeros(3)
                else:
                    # Use j's tangential direction
                    if vtj > 1e-10:
                        vt_vec_star = vt_vec_j * (vt_star_mag / vtj)
                    else:
                        vt_vec_star = np.zeros(3)
                        
                v_star_vec = v_star_n * nij + vt_vec_star
                
                # Volume factor
                # Derivation in Sec 2.3 shows Vp^2 term.
                # V_ij^2 = 0.5 * (Vp_i^2 + Vp_j^2)
                Vij2 = 0.5 * (pi.Vp**2 + pj.Vp**2)
                
                # Gradients
                # grad_i W(xi - xj, h)
                # We use sqrt(2) * h as per paper
                # grad_i W(xi - xj, sqrt(2)hi)
                
                grad_Wi = self.grad_kernel(coords[i] - coords[j], r, np.sqrt(2.0) * pi.h)
                grad_Wj = self.grad_kernel(coords[i] - coords[j], r, np.sqrt(2.0) * pj.h)
                
                # Term: [grad_i W(...hi) - grad_j W(...hj)]
                # Note: grad_j W(xi - xj) = - grad_i W(xi - xj) if h is same.
                # But here h is different.
                # grad_j W(xi - xj, h) is derivative wrt xj.
                # W depends on (xi - xj). Let R = xi - xj.
                # dW/dxj = dW/dR * dR/dxj = dW/dR * (-1).
                # dW/dxi = dW/dR * (1).
                # So grad_j W = - grad_i W (functionally).
                # self.grad_kernel returns grad wrt first arg (xi).
                # So grad_Wj_val = grad_kernel(xi - xj, r, hj).
                # This is grad_i W(..., hj).
                # We need grad_j W(..., hj) = - grad_i W(..., hj).
                
                term_grad = grad_Wi - (-self.grad_kernel(coords[i] - coords[j], r, np.sqrt(2.0) * pj.h))
                # term_grad = grad_Wi + grad_Wj_val
                
                # Force term
                # Eq 371: <nu_i S_i_dot> = - sum P* Vij2 term_grad
                force = - P_star * Vij2 * term_grad
                
                # Energy term
                # Eq 373: <nu_i e_i_dot> = - sum P* v* . Vij2 term_grad
                # = v* . force
                power = np.dot(v_star_vec, force)
                
                # Add to i
                pi.dSdt += force
                pi.dedt += power
                
                # Add to j (Newton's 3rd law)
                # Force on j:
                # Swap i and j.
                # P* same. Vij2 same.
                # term_grad_ji = [grad_j W(xj - xi, hj) - grad_i W(xj - xi, hi)]
                # grad_j W(xj - xi) = grad wrt xj of W(xj - xi).
                # = grad_kernel(xj - xi, r, hj).
                # grad_i W(xj - xi) = - grad_kernel(xj - xi, r, hi).
                # term_grad_ji = grad_kernel(xj - xi, hj) + grad_kernel(xj - xi, hi).
                # Note: grad_kernel(v) = -2 v / h^2 W.
                # grad_kernel(xj - xi) = -2 (xj - xi) ... = -2 (- (xi - xj)) ... = 2 (xi - xj) ...
                # = - grad_kernel(xi - xj).
                # So term_grad_ji = - term_grad.
                # So force on j is -force on i.
                
                pj.dSdt -= force
                pj.dedt -= power
                
                if np.linalg.norm(force) > 0:
                    nonzero_force_count += 1

        # print(f"Interactions: {interaction_count}, Non-zero forces: {nonzero_force_count}")

    def step(self, dt):
        # Compute forces first
        self.compute_forces()
        
        # Euler integration
        for p in self.particles:
            # Update conserved variables
            # <nu S>_new = <nu S>_old + <nu S_dot> dt
            # S_new = S_old + S_dot dt (since nu is constant per particle in integration?)
            # Eq 425: <nu S>^{n+1} = <nu S>^n + <nu S_dot> dt
            # If nu is constant for particle i (it is parameter), then:
            # nu S_new = nu S_old + (force) dt
            # S_new = S_old + (force / nu) dt
            
            p.S += (p.dSdt / p.nu) * dt
            p.e += (p.dedt / p.nu) * dt
            
            # Update position
            # x^{n+1} = x^n + v^n dt
            p.x += p.vx * dt
            p.y += p.vy * dt
            p.z += p.vz * dt
            
        self.time += dt
        
        # Update smoothing length
        self.update_smoothing_length()
        
        # Debug stats
        hs = [p.h for p in self.particles]
        Vps = [p.Vp for p in self.particles]
        print(f"Step end: h_min={min(hs):.4e}, h_max={max(hs):.4e}, Vp_min={min(Vps):.4e}, Vp_max={max(Vps):.4e}")
        
        # Recover primitives
        self.recover_primitives()

