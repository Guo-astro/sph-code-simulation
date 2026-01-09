# Email to Professor

---

**Subject:** IMBH-Cloud Simulation: Switching from Isothermal to γ=5/3 with Cooling

---

Dear Professor [Name],

I am writing to explain a necessary change in my simulation approach for the IMBH-cloud tidal interaction project.

## The Problem

The isothermal equation of state (γ=1) cannot produce stable initial conditions for our cloud parameters:

- For M = 40 M_sun at T = 15 K, the Bonnor-Ebert critical mass is only ~1-2 M_sun
- Our cloud is **20-40× super-critical** → no equilibrium exists
- Numerical test confirmed this: the cloud collapsed immediately with timestep dropping to dt ~ 10⁻¹⁰ Myr

## The Solution: Hybrid Approach

I propose using **γ = 5/3 polytrope with cooling switched ON during the flyby**:

| Phase | γ | Cooling | Purpose |
|-------|---|---------|---------|
| 1-2. Relaxation & Test | 5/3 | OFF | Stable equilibrium (Lane-Emden n=1.5) |
| 3. IMBH Flyby | 5/3 | **ON** | Isothermal shocks form dense clump |

## Why This Works

**Isothermal shocks are required to explain observations:**

The observed dense clump (n ~ 10^6.5 cm⁻³) requires compression of ~100-1000×.

- Adiabatic shocks: max compression ρ₂/ρ₁ ≤ 4 (insufficient)
- Isothermal shocks: ρ₂/ρ₁ = M² → can reach 100-1000× (matches observations)

**The physics during flyby:**
1. Tidal compression heats gas adiabatically (T → 10³ K)
2. Cooling rapidly removes thermal energy (t_cool << t_compress)
3. Loss of pressure support triggers further collapse
4. Dense, cold clump forms via isothermal shock (n → 10^6.5 cm⁻³, T → 60 K)

## Summary

| Approach | Stability | Dense Clump | Observed T=60K |
|----------|-----------|-------------|----------------|
| Isothermal (γ=1) | Unstable | - | - |
| Adiabatic only (γ=5/3) | Stable | Limited (ρ×4 max) | Too hot |
| **γ=5/3 + Cooling** | **Stable** | **Yes (ρ×1000)** | **Yes** |

The hybrid approach resolves the stability problem while naturally producing the observed dense clump through cooling-triggered isothermal shocks.

I would appreciate your feedback on this approach.

Best regards,
[Your Name]
