import numpy as np
import matplotlib.pyplot as plt
from srg_sph import Particle, SRGSPH
import time

def run_sod_test():
    # Parameters
    gamma_ad = 5.0/3.0
    n_particles_side = 100 # Reduced for speed in this demo
    domain_size = 1.0
    x_min = -0.5
    x_max = 0.5
    t_end = 0.35
    
    # Left state
    PL = 1.0
    nL = 1.0
    vxL = 0.0
    
    # Right state
    PR = 0.1
    nR = 0.125
    vxR = 0.0
    
    particles = []
    
    # Spacing
    dx = domain_size / (2 * n_particles_side)
    
    # Initialize Left
    X_factor = gamma_ad / (gamma_ad - 1.0)
    
    for i in range(n_particles_side):
        x = x_min + (i + 0.5) * dx
        
        # Primitive to Conservative
        u = PL / ((gamma_ad - 1.0) * nL)
        H = 1.0 + u + PL/nL
        gamma = 1.0 / np.sqrt(1.0 - vxL**2)
        N = gamma * nL
        
        S = gamma * H * vxL
        # e consistent with recover_primitives
        # H(X g^2 - 1) = X g e - 1 => e = (H(X g^2 - 1) + 1) / (X g)
        e = (H * (X_factor * gamma**2 - 1.0) + 1.0) / (X_factor * gamma)
        
        # Baryon number nu
        # N = nu / Vp => nu = N * Vp
        # Initially Vp approx dx (in 1D)
        nu = N * dx
        
        p = Particle(i, x, 0, 0, nu, [S, 0, 0], e, h=dx)
        p.rho = nL
        p.p = PL
        p.vx = vxL
        p.gamma = gamma
        p.h_enthalpy = H
        p.cs = np.sqrt((gamma_ad - 1.0) * (H - 1.0) / H)
        
        particles.append(p)
        
    # Initialize Right
    for i in range(n_particles_side):
        x = 0.0 + (i + 0.5) * dx
        
        u = PR / ((gamma_ad - 1.0) * nR)
        H = 1.0 + u + PR/nR
        gamma = 1.0 / np.sqrt(1.0 - vxR**2)
        N = gamma * nR
        
        S = gamma * H * vxR
        # e consistent with recover_primitives
        e = (H * (X_factor * gamma**2 - 1.0) + 1.0) / (X_factor * gamma)
        
        nu = N * dx
        
        p = Particle(i + n_particles_side, x, 0, 0, nu, [S, 0, 0], e, h=dx)
        p.rho = nR
        p.p = PR
        p.vx = vxR
        p.gamma = gamma
        p.h_enthalpy = H
        p.cs = np.sqrt((gamma_ad - 1.0) * (H - 1.0) / H)
        
        particles.append(p)
        
    # Simulation
    sim = SRGSPH(particles, gamma_ad=gamma_ad, dim=1, eta=1.2, C_smooth=2.0)
    
    # Initial smoothing length update
    print("Initializing smoothing lengths...")
    # Iterate a few times to converge h
    for k in range(5):
        sim.update_smoothing_length()
        hs = [p.h for p in sim.particles]
        print(f"Iter {k}: h_min={min(hs):.4f}, h_max={max(hs):.4f}")
        
    # Check Vp and nu
    Vps = [p.Vp for p in sim.particles]
    nus = [p.nu for p in sim.particles]
    print(f"Vp_min={min(Vps):.4e}, Vp_max={max(Vps):.4e}")
    print(f"nu_min={min(nus):.4e}, nu_max={max(nus):.4e}")
    Vij2_est = 0.5 * ((min(Vps)/min(nus))**2 + (min(Vps)/min(nus))**2)
    print(f"Est Vij2 (min) = {Vij2_est:.4e}")
    
    print("Starting simulation...")
    step = 0
    while sim.time < t_end:
        # CFL condition
        min_dt = 1e9
        for p in sim.particles:
            if p.cs > 0:
                dt_p = p.h / p.cs
                if dt_p < min_dt: min_dt = dt_p
        
        dt = 0.3 * min_dt # CFL = 0.3
        if sim.time + dt > t_end:
            dt = t_end - sim.time
            
        # Compute forces
        sim.compute_forces()
        
        # Integrate
        sim.step(dt)
        
        step += 1
        if step % 10 == 0:
            print(f"Step {step}, Time: {sim.time:.4f}, dt: {dt:.4e}")
            
    # Output results
    print("Simulation finished.")
    
    # Collect data
    xs = [p.x for p in sim.particles]
    ps = [p.p for p in sim.particles]
    rhos = [p.rho for p in sim.particles]
    vxs = [p.vx for p in sim.particles]
    
    # Save to file
    with open('sod_results.txt', 'w') as f:
        f.write("x,p,rho,vx\n")
        for i in range(len(xs)):
            f.write(f"{xs[i]},{ps[i]},{rhos[i]},{vxs[i]}\n")
            
    print("Results saved to sod_results.txt")

if __name__ == "__main__":
    run_sod_test()
