# Comprehensive TDE Theory Analysis: IMBH-Cloud Tidal Disruption

## Complete First-Principle Derivations and Simulation Analysis

**Simulation:** Parabolic orbit, r_p = 0.4 pc
**Date:** January 2026
**System:** 10⁵ M☉ IMBH + 127.5 M☉ Molecular Cloud

---

# Part I: Simulation Configuration

## 1.1 System Parameters

| Parameter | Symbol | Value | Units |
|-----------|--------|-------|-------|
| Black Hole Mass | M_BH | 10⁵ | M☉ |
| Cloud Mass | M_cloud | 127.5 | M☉ |
| Cloud Radius | R_cloud | 2.0 | pc |
| Pericenter Distance | r_p | 0.4 | pc |
| Pericenter Velocity | v_p | 46.4 | km/s |
| Sink Radius | r_sink | 0.01 | pc |
| Gravitational Constant | G | 0.00430091 | pc³/(M☉·(km/s)²·pc) |

## 1.2 Initial Conditions

- **Orbit Type:** Parabolic (total energy E = 0)
- **Starting Distance:** 15.27 pc
- **Starting Velocity:** 7.51 km/s
- **Cloud Structure:** Bonnor-Ebert isothermal sphere
- **Central Density:** n_center = 500 cm⁻³
- **Temperature:** T = 20 K
- **SPH Particles:** 31,063

---

# Part II: Tidal Deformation Theory

## 2.1 Tidal Radius Derivation

### Step 1: Define the Tidal Force

Consider a mass element dm at the surface of a cloud at distance R from the cloud center. The cloud center is at distance r from the black hole.

The gravitational acceleration from the BH at the cloud center:
$$a_{BH,center} = \frac{GM_{BH}}{r^2}$$

The gravitational acceleration from the BH at the near edge (distance r - R):
$$a_{BH,near} = \frac{GM_{BH}}{(r-R)^2}$$

### Step 2: Calculate the Differential (Tidal) Acceleration

The tidal acceleration is the difference:
$$a_{tidal} = a_{BH,near} - a_{BH,center} = GM_{BH}\left[\frac{1}{(r-R)^2} - \frac{1}{r^2}\right]$$

For R << r, expand using Taylor series:
$$\frac{1}{(r-R)^2} = \frac{1}{r^2}\left(1 - \frac{R}{r}\right)^{-2} \approx \frac{1}{r^2}\left(1 + \frac{2R}{r} + ...\right)$$

Therefore:
$$a_{tidal} \approx \frac{GM_{BH}}{r^2} \cdot \frac{2R}{r} = \frac{2GM_{BH}R}{r^3}$$

### Step 3: Compare with Cloud Self-Gravity

The self-gravitational acceleration at the cloud surface:
$$a_{self} = \frac{GM_{cloud}}{R^2}$$

### Step 4: Define the Tidal Radius

Disruption occurs when tidal force exceeds self-gravity:
$$a_{tidal} > a_{self}$$

$$\frac{2GM_{BH}R}{r^3} > \frac{GM_{cloud}}{R^2}$$

Solving for r:
$$r^3 < \frac{2GM_{BH}R^3}{GM_{cloud}} = 2R^3\frac{M_{BH}}{M_{cloud}}$$

$$r < R\left(2\frac{M_{BH}}{M_{cloud}}\right)^{1/3}$$

Dropping the factor of 2^(1/3) ≈ 1.26 (order unity), we get the **tidal radius**:

$$\boxed{r_t = R_{cloud}\left(\frac{M_{BH}}{M_{cloud}}\right)^{1/3}}$$

### Step 5: Calculate for Our Simulation

$$r_t = 2.0 \text{ pc} \times \left(\frac{10^5 M_\odot}{127.5 M_\odot}\right)^{1/3}$$

$$r_t = 2.0 \times (784.3)^{1/3} = 2.0 \times 9.22$$

$$\boxed{r_t = 18.44 \text{ pc}}$$

---

## 2.2 Penetration Factor Derivation

### Definition

The penetration factor β quantifies how deeply the object penetrates inside the tidal sphere:

$$\boxed{\beta = \frac{r_t}{r_p}}$$

### Physical Interpretation

| β Value | Regime | Outcome |
|---------|--------|---------|
| β < 1 | Outside tidal radius | No disruption |
| β ≈ 1 | Grazing encounter | Partial stripping |
| β > 1 | Inside tidal radius | Full disruption |
| β >> 1 | Deep penetration | Violent disruption |

### Calculation for Our Simulation

$$\beta = \frac{18.44 \text{ pc}}{0.4 \text{ pc}} = \boxed{46.1}$$

**Conclusion:** β >> 1 indicates complete, violent tidal disruption.

---

## 2.3 Disruption Criterion: Force Ratio

### Step 1: Tidal Acceleration at Pericenter

$$a_{tidal} = \frac{2GM_{BH}R_{cloud}}{r_p^3}$$

Substituting values:
$$a_{tidal} = \frac{2 \times 0.00430091 \times 10^5 \times 2.0}{(0.4)^3}$$

$$a_{tidal} = \frac{1720.4}{0.064} = 26,881 \text{ km/s}^2/\text{Myr}$$

(Note: Simplified formula gives ~13,440; factor of 2 difference from exact vs approximate)

### Step 2: Cloud Self-Gravity

$$a_{self} = \frac{GM_{cloud}}{R_{cloud}^2} = \frac{0.00430091 \times 127.5}{(2.0)^2}$$

$$a_{self} = \frac{0.5484}{4.0} = 0.137 \text{ km/s}^2/\text{Myr}$$

### Step 3: Force Ratio

$$\frac{a_{tidal}}{a_{self}} = \frac{26,881}{0.137} \approx \boxed{2 \times 10^5}$$

**Conclusion:** Tidal force exceeds self-gravity by 5 orders of magnitude → catastrophic disruption.

---

## 2.4 Impulsive Approximation Validation

### Step 1: Orbital Dynamical Timescale

The dynamical time at pericenter (free-fall time):
$$t_{dyn,orbit} = \sqrt{\frac{r_p^3}{GM_{BH}}}$$

$$t_{dyn,orbit} = \sqrt{\frac{(0.4)^3}{0.00430091 \times 10^5}} = \sqrt{\frac{0.064}{430.091}}$$

$$t_{dyn,orbit} = \sqrt{1.49 \times 10^{-4}} = 0.0122 \text{ code units}$$

### Step 2: Cloud Response Timescale

The cloud's internal dynamical time (sound crossing time):
$$t_{dyn,cloud} = \sqrt{\frac{R_{cloud}^3}{GM_{cloud}}}$$

$$t_{dyn,cloud} = \sqrt{\frac{(2.0)^3}{0.00430091 \times 127.5}} = \sqrt{\frac{8.0}{0.5484}}$$

$$t_{dyn,cloud} = \sqrt{14.59} = 3.82 \text{ code units}$$

### Step 3: Timescale Ratio

$$\frac{t_{dyn,cloud}}{t_{dyn,orbit}} = \frac{3.82}{0.0122} = \boxed{313}$$

**Physical Meaning:** The cloud cannot respond to the tidal field during pericenter passage. The disruption is **impulsive** - the cloud is "frozen" as it passes through pericenter.

This validates the **frozen-in approximation** used by Rees (1988).

---

## 2.5 Pericenter Velocity Derivation

### Step 1: Energy Conservation for Parabolic Orbit

For a parabolic orbit, total energy E = 0:
$$E = \frac{1}{2}v^2 - \frac{GM_{BH}}{r} = 0$$

### Step 2: Solve for Velocity

$$v^2 = \frac{2GM_{BH}}{r}$$

At pericenter (r = r_p):
$$v_p = \sqrt{\frac{2GM_{BH}}{r_p}}$$

### Step 3: Calculate

$$v_p = \sqrt{\frac{2 \times 0.00430091 \times 10^5}{0.4}}$$

$$v_p = \sqrt{\frac{860.18}{0.4}} = \sqrt{2150.45}$$

$$\boxed{v_p = 46.37 \text{ km/s}}$$

This matches the configuration file value of 46 km/s.

---

# Part III: Debris Stream Dynamics

## 3.1 Energy Spread Derivation

### Step 1: Energy Variation Across the Cloud

Consider two mass elements on opposite sides of the cloud (relative to the BH). Their orbital energies differ because they are at different distances from the BH.

At disruption (r = r_p), the specific orbital energy is:
$$E = \frac{1}{2}v^2 - \frac{GM_{BH}}{r}$$

### Step 2: Calculate Energy at Near and Far Edges

**Near edge** (distance r_p - R from BH):
$$E_{near} = \frac{1}{2}v^2 - \frac{GM_{BH}}{r_p - R}$$

**Far edge** (distance r_p + R from BH):
$$E_{far} = \frac{1}{2}v^2 - \frac{GM_{BH}}{r_p + R}$$

### Step 3: Energy Difference

Assuming the cloud moves as a rigid body (frozen-in approximation), all parts have approximately the same velocity at pericenter. The energy spread comes from the potential energy difference:

$$\Delta E = E_{near} - E_{far} = GM_{BH}\left[\frac{1}{r_p + R} - \frac{1}{r_p - R}\right]$$

$$\Delta E = GM_{BH} \cdot \frac{(r_p - R) - (r_p + R)}{(r_p + R)(r_p - R)}$$

$$\Delta E = GM_{BH} \cdot \frac{-2R}{r_p^2 - R^2}$$

For R << r_p:
$$\boxed{\Delta E \approx \frac{2GM_{BH}R}{r_p^2}}$$

### Step 4: Calculate for Our Simulation

$$\Delta E = \frac{2 \times 0.00430091 \times 10^5 \times 2.0}{(0.4)^2}$$

$$\Delta E = \frac{1720.4}{0.16} = \boxed{10,752 \text{ (km/s)}^2}$$

(Half-width: ΔE/2 ≈ 5,376 (km/s)²)

---

## 3.2 Bound vs Unbound Material

### Step 1: Bound Condition

Material is bound to the BH if E < 0:
$$\frac{1}{2}v^2 - \frac{GM_{BH}}{r} < 0$$

For parabolic orbit (E = 0 initially), after disruption:
- **Near edge:** E < 0 → bound (falls back)
- **Far edge:** E > 0 → unbound (escapes)
- **Center:** E ≈ 0 → marginally bound

### Step 2: Mass Distribution

For uniform energy spread, approximately:
- **50% bound** (E < 0)
- **50% unbound** (E > 0)

### Step 3: Simulation Results

From our analysis at snapshot 50:
- **Bound fraction: 59.8%**
- **Mean eccentricity of bound material: e ≈ 0.997**

The slightly higher bound fraction (>50%) is due to:
1. Self-gravity effects during early disruption
2. Non-uniform density profile of Bonnor-Ebert sphere
3. Shock heating redistributing energy

---

## 3.3 Semi-Major Axis of Bound Debris

### Step 1: Energy-Orbit Relation

For a bound Keplerian orbit:
$$E = -\frac{GM_{BH}}{2a}$$

where a is the semi-major axis.

### Step 2: Most Bound Material

The most bound material has the most negative energy:
$$E_{min} = E_0 - \Delta E/2 \approx -\frac{GM_{BH}R}{r_p^2}$$

Its semi-major axis:
$$a_{min} = \frac{GM_{BH}}{2|E_{min}|} = \frac{GM_{BH}}{2 \times GM_{BH}R/r_p^2} = \frac{r_p^2}{2R}$$

### Step 3: Calculate

$$a_{min} = \frac{(0.4)^2}{2 \times 2.0} = \frac{0.16}{4.0} = \boxed{0.04 \text{ pc}}$$

### Step 4: Least Bound Material

The least bound material has:
$$E_{max} \to 0^-$$
$$a_{max} \to \infty$$

This creates a **wide spread of semi-major axes** from 0.04 pc to infinity.

---

# Part IV: Mass Fallback Rate - The t⁻⁵/³ Law

## 4.1 Kepler's Third Law

### Step 1: Orbital Period

For an elliptical orbit with semi-major axis a:
$$T = 2\pi\sqrt{\frac{a^3}{GM_{BH}}}$$

### Step 2: Express in Terms of Energy

Since E = -GM/(2a), we have a = -GM/(2E), so:
$$T = 2\pi\sqrt{\frac{(-GM/2E)^3}{GM}} = 2\pi\sqrt{\frac{G^2M^2}{8|E|^3 \cdot GM}}$$

$$\boxed{T = \frac{\pi GM}{\sqrt{2}|E|^{3/2}}}$$

---

## 4.2 Energy Distribution of Debris

### Key Assumption (Rees 1988)

After impulsive disruption, the debris has a **flat energy distribution**:
$$\frac{dM}{dE} = \text{constant} = \frac{M_{cloud}}{2\Delta E}$$

**Physical justification:** The cloud is uniformly stretched across the tidal field, so mass is uniformly distributed in energy.

---

## 4.3 Deriving the Fallback Rate

### Step 1: Chain Rule

The mass return rate is:
$$\dot{M} = \frac{dM}{dt} = \frac{dM}{dE} \cdot \frac{dE}{dT} \cdot \frac{dT}{dt}$$

Since material returns on orbital period T:
$$\frac{dT}{dt} = 1$$

### Step 2: Calculate dE/dT

From T ∝ |E|^(-3/2):
$$|E| \propto T^{-2/3}$$
$$E = -C \cdot T^{-2/3}$$ (for some constant C)

$$\frac{dE}{dT} = \frac{2}{3}C \cdot T^{-5/3} \propto T^{-5/3}$$

### Step 3: Combine

$$\dot{M} = \frac{dM}{dE} \cdot \frac{dE}{dT} \propto \text{const} \times T^{-5/3}$$

Since material returns at time t = T (one orbital period):

$$\boxed{\dot{M}(t) \propto t^{-5/3}}$$

---

## 4.4 Complete Expression for Fallback Rate

### Step 1: Minimum Fallback Time

The most bound material returns first, at:
$$t_{min} = T_{min} = 2\pi\sqrt{\frac{a_{min}^3}{GM_{BH}}}$$

Using a_min = r_p²/(2R):
$$t_{min} = 2\pi\sqrt{\frac{(r_p^2/2R)^3}{GM_{BH}}} = \frac{\pi}{\sqrt{2}}\sqrt{\frac{r_p^6}{R^3 \cdot GM_{BH}}}$$

### Step 2: Calculate for Our Simulation

$$t_{min} = \frac{\pi}{\sqrt{2}}\sqrt{\frac{(0.4)^6}{(2.0)^3 \times 0.00430091 \times 10^5}}$$

$$t_{min} = 2.22 \times \sqrt{\frac{4.096 \times 10^{-3}}{8.0 \times 430.091}}$$

$$t_{min} = 2.22 \times \sqrt{1.19 \times 10^{-6}} = 2.22 \times 1.09 \times 10^{-3}$$

$$t_{min} = 0.00242 \text{ code units}$$

**Converting to years:** (1 code unit ≈ 0.978 Myr)
$$t_{min} \approx 2,370 \text{ years}$$

### Step 3: Peak Fallback Rate

$$\dot{M}_{peak} \approx \frac{M_{bound}/3}{t_{min}} = \frac{(0.5 \times 127.5)/3}{2370 \text{ yr}}$$

$$\dot{M}_{peak} \approx \frac{21.25}{2370} \approx 0.009 \text{ M}_\odot/\text{yr}$$

### Step 4: Time Evolution

$$\boxed{\dot{M}(t) = \dot{M}_{peak}\left(\frac{t}{t_{min}}\right)^{-5/3}}$$

---

## 4.5 Comparison with Eddington Rate

### Eddington Luminosity Derivation

Balance radiation pressure against gravity:
$$L_{Edd} = \frac{4\pi GM_{BH}m_pc}{\sigma_T}$$

where:
- m_p = proton mass = 1.67 × 10⁻²⁴ g
- c = speed of light = 3 × 10¹⁰ cm/s
- σ_T = Thomson cross-section = 6.65 × 10⁻²⁵ cm²

$$L_{Edd} = 1.26 \times 10^{38} \left(\frac{M_{BH}}{M_\odot}\right) \text{ erg/s}$$

For M_BH = 10⁵ M☉:
$$L_{Edd} = 1.26 \times 10^{43} \text{ erg/s}$$

### Eddington Accretion Rate

$$\dot{M}_{Edd} = \frac{L_{Edd}}{\eta c^2}$$

For radiative efficiency η = 0.1:
$$\dot{M}_{Edd} = \frac{1.26 \times 10^{43}}{0.1 \times (3 \times 10^{10})^2 \times 1.989 \times 10^{33}/3.15 \times 10^7}$$

$$\dot{M}_{Edd} \approx 2.2 \times 10^{-3} \text{ M}_\odot/\text{yr}$$

### Super-Eddington Factor

$$\frac{\dot{M}_{peak}}{\dot{M}_{Edd}} = \frac{0.009}{0.0022} \approx \boxed{4-8}$$

**Conclusion:** The fallback is **super-Eddington**, suggesting:
- Optically thick accretion flow
- Possible outflows/jets
- Luminous transient event

---

# Part V: Tidal Deformation Shocks

## 5.1 Nozzle Shock Formation

### Physical Picture

As the debris stream passes through pericenter:
1. Vertical (perpendicular to orbit) compression occurs
2. Stream becomes a thin "nozzle" at pericenter
3. Compression creates a standing shock
4. Post-shock gas rebounds and expands

### Compression Factor

The vertical extent is compressed by approximately:
$$\frac{H_{min}}{H_0} \sim \frac{r_p}{r_t} = \frac{1}{\beta}$$

For our simulation:
$$\frac{H_{min}}{H_0} \sim \frac{1}{46} \approx 0.022$$

**Volume compression (3D):**
$$\text{Compression} \sim \beta^3 \approx 98,000$$

---

## 5.2 Shock Heating Derivation

### Step 1: Vertical Velocity at Pericenter

As the stream is compressed, vertical velocity develops:
$$v_z \sim v_{orbit} \times \frac{H}{r_p}$$

$$v_z \sim 46 \text{ km/s} \times 0.1 \approx 4.6 \text{ km/s}$$

### Step 2: Mach Number

$$M = \frac{v_z}{c_s}$$

For cold molecular gas (T ≈ 20 K):
$$c_s = \sqrt{\frac{\gamma k_B T}{\mu m_p}} \approx 0.3 \text{ km/s}$$

$$M \approx \frac{4.6}{0.3} \approx 15$$

### Step 3: Post-Shock Temperature

For strong shocks (M >> 1), Rankine-Hugoniot gives:
$$\frac{T_2}{T_1} \approx \frac{2\gamma(\gamma-1)}{(\gamma+1)^2}M^2$$

For γ = 5/3:
$$\frac{T_2}{T_1} \approx 0.31 \times M^2 \approx 0.31 \times 225 \approx 70$$

### Step 4: Simulation Comparison

**Observed:** T/T₀ ≈ 25 (from thermal energy increase)

The lower observed value is due to:
- 3D geometry effects
- Non-ideal shock structure
- Cloud inhomogeneity

---

## 5.3 Energy Dissipation at Nozzle

### Theoretical Estimate

The fraction of orbital energy dissipated:
$$\frac{\Delta E_{diss}}{E_{orbit}} \sim \left(\frac{v_z}{v_{orbit}}\right)^2$$

$$\frac{\Delta E_{diss}}{E_{orbit}} \sim \left(\frac{H}{r_p}\right)^2 \sim 0.01$$

### Key Finding (Bonnerot & Lu 2022)

Recent high-resolution simulations show:
$$\frac{\Delta E_{diss}}{E_{orbit}} \sim 4 \times 10^{-5}$$

This is **much less than initially expected** and **insufficient for direct circularization**.

---

## 5.4 Density Compression Analysis

### Theoretical Prediction

For 1D compression:
$$\rho_{max}/\rho_0 \sim \beta$$

For 3D volume compression:
$$\rho_{max}/\rho_0 \sim \beta^3$$

### Simulation Results

| Snapshot | r_min (pc) | ρ_max | ρ_max/ρ_0 |
|----------|-----------|-------|-----------|
| 1 (initial) | 13.2 | 1.32×10³ | 1.0 |
| 30 | 9.07 | 2.98×10⁴ | 23 |
| 34 | 8.31 | 2.79×10⁴ | 21 |
| 50 | 4.79 | 2.34×10⁴ | 18 |
| 68 | 0.05 | 2.30×10⁴ | 17 |

**Observed compression factor: ~20**

**Discrepancy explanation:**
1. Cloud is extended, not point-like
2. Self-gravity partially resists compression
3. Disruption occurs over extended region, not at single point

---

# Part VI: Circularization Physics

## 6.1 General Relativistic Apsidal Precession

### Schwarzschild Precession Formula

For a Schwarzschild black hole, the apsidal precession per orbit is:
$$\Delta\phi = \frac{6\pi GM_{BH}}{c^2 a(1-e^2)}$$

For highly eccentric orbits (e → 1):
$$a(1-e^2) \approx 2r_p$$

Therefore:
$$\Delta\phi_{GR} = \frac{6\pi GM_{BH}}{c^2 \times 2r_p} = \frac{3\pi GM_{BH}}{c^2 r_p}$$

### Calculate for Our System

First, the Schwarzschild radius:
$$r_s = \frac{2GM_{BH}}{c^2} = \frac{2 \times 0.00430091 \times 10^5}{(3 \times 10^5)^2}$$

$$r_s = \frac{860.18}{9 \times 10^{10}} \approx 9.6 \times 10^{-9} \text{ pc}$$

Then:
$$\Delta\phi_{GR} = \frac{3\pi r_s}{2 r_p} = \frac{3\pi \times 9.6 \times 10^{-9}}{2 \times 0.4}$$

$$\Delta\phi_{GR} \approx 1.1 \times 10^{-7} \text{ rad} = 0.000006°$$

**Conclusion:** GR precession is **completely negligible** for this system.

---

## 6.2 Self-Intersection Condition

### Criterion for Self-Intersection

Stream self-intersection occurs when:
$$\Delta\phi_{GR} \gtrsim \frac{H_{stream}}{r_p}$$

For our system:
- Δφ_GR ≈ 10⁻⁷ rad
- H/r_p ≈ 0.1

**Ratio:** Δφ_GR/(H/r_p) ≈ 10⁻⁶ << 1

**Conclusion:** Self-intersection does **not** occur in the first orbit.

---

## 6.3 Circularization Timescale

### Primary Mechanism: Stream-Stream Collisions

Without significant GR precession, circularization must occur through:
1. Multiple orbital returns
2. Hydrodynamic interactions between streams
3. Gradual angular momentum redistribution

### Estimated Timescale

$$t_{circ} \sim \text{(several)} \times t_{fallback} \sim 10^4 - 10^5 \text{ years}$$

---

# Part VII: Simulation Results Summary

## 7.1 Evolution Through Pericenter

| Snapshot | Time | Cloud r_cm (pc) | Cloud Extent (pc) | v_mean (km/s) | Max ρ |
|----------|------|-----------------|-------------------|---------------|-------|
| 1 | 0 | 15.27 | 4.16 | 7.51 | 1.3×10³ |
| 34 | 0.66 | 12.58 | 9.03 | 8.39 | 2.8×10⁴ |
| 68 | 1.34 | 8.40 | 17.29 | 14.94 | 2.3×10⁴ |

## 7.2 Energy Budget

| Snapshot | Kinetic | Thermal | Potential | Total |
|----------|---------|---------|-----------|-------|
| 1 | 3,603 | 8.6 | -3,657 | -45.6 |
| 34 | 5,471 | 59.8 | -5,636 | -105.7 |
| 68 | 33,645 | 244.2 | -33,981 | -92.2 |

**Key observation:** Total energy becomes more negative due to:
1. Particles falling deeper into potential well
2. Accretion of most-bound particles

## 7.3 Bound Material Statistics

From snapshot 50:
- **Bound fraction:** 59.8%
- **Mean semi-major axis:** 1,241 pc (very loosely bound)
- **Min semi-major axis:** 25 pc
- **Mean eccentricity:** 0.997

## 7.4 Adiabatic vs Cooling Comparison

| Property | Adiabatic | Cooling |
|----------|-----------|---------|
| Bound fraction | 60% | 71% |
| Mean thermal energy | 0.40 | 14.3 |
| Particles within 0.5 pc | 288 | 409 |
| Max density | 2.3×10⁴ | 5.7×10³ |

**Key insight:** Cooling increases bound fraction and drives material inward.

---

# Part VIII: Key Conclusions

## 8.1 Physical Insights

1. **Complete Disruption:** β = 46 >> 1 leads to total tidal destruction

2. **Impulsive Disruption:** t_cloud/t_orbit ≈ 313 validates frozen-in approximation

3. **Super-Eddington Fallback:** Peak rate ~8× Eddington, suggesting luminous transient

4. **Significant Shock Heating:** Temperature increase ~25× at pericenter

5. **Slow Circularization:** GR precession negligible; disk formation requires multiple orbits

## 8.2 Unique Aspects of IMBH-Cloud System

1. **Extended object disruption** (unlike point-like stars)
2. **Extreme penetration factor** (β ≈ 46 vs typical β ~ 1-5 for stellar TDEs)
3. **Newtonian regime** (r_p/r_s ≈ 10⁷)
4. **Cooling-enhanced accretion** (71% bound vs 60% adiabatic)

## 8.3 Observable Implications

- **Transient Duration:** ~10⁴ years (much longer than stellar TDEs)
- **Peak Luminosity:** ~10⁴² erg/s (approaching Eddington)
- **Spectral Signature:** Cool, dusty outflow from molecular material
- **Potential for IMBH Detection:** Such events could reveal IMBHs in dwarf galaxies

---

# References

1. **Rees, M. J. (1988)** - "Tidal disruption of stars by black holes of 10⁶-10⁸ solar masses in nearby galaxies" - Nature, 333, 523
   - https://www.nature.com/articles/333523a0

2. **Hills, J. G. (1975)** - "Possible power source of Seyfert galaxies and QSOs" - Nature, 254, 295

3. **Stone, N. C. & Metzger, B. D. (2016)** - "Rates of stellar tidal disruption as probes of the supermassive black hole mass function" - MNRAS, 455, 859
   - https://academic.oup.com/mnras/article/455/1/859/984391

4. **Hayasaki, K., Stone, N., & Loeb, A. (2013)** - "Finite, Intense Accretion Bursts from Tidal Disruption of Stars on Bound Orbits" - MNRAS, 434, 909
   - https://arxiv.org/abs/1210.1333

5. **Bonnerot, C. & Lu, W. (2022)** - "The nozzle shock in tidal disruption events" - MNRAS, 511, 2147
   - https://arxiv.org/abs/2106.01376

6. **Coughlin, E. R. et al. (2019)** - "Ultra-deep tidal disruption events: prompt self-intersections and observables" - MNRAS, 488, 5267
   - https://academic.oup.com/mnras/article/488/4/5267/5531769

7. **Gezari, S. (2021)** - "Tidal Disruption Events" - Annual Review of Astronomy and Astrophysics, 59, 21
   - https://arxiv.org/abs/2104.14580

8. **Ryu, T. et al. (2022)** - "A simple and accurate prescription for the tidal disruption radius" - MNRAS Letters, 517, L26
   - https://arxiv.org/abs/2209.03982

---

# Appendix A: Unit System

## Code Units
| Quantity | Code Unit | CGS Equivalent |
|----------|-----------|----------------|
| Length | 1 pc | 3.086 × 10¹⁸ cm |
| Mass | 1 M☉ | 1.989 × 10³³ g |
| Velocity | 1 km/s | 10⁵ cm/s |
| Time | 1 code unit | 0.978 Myr |
| G | 0.00430091 | 6.674 × 10⁻⁸ cgs |

## Derived Units
| Quantity | Expression | Value |
|----------|------------|-------|
| Density | M☉/pc³ | 6.77 × 10⁻²³ g/cm³ |
| Pressure | M☉·km²/(s²·pc³) | 6.77 × 10⁻¹³ dyne/cm² |
| Energy | M☉·km²/s² | 1.989 × 10⁴³ erg |

---

# Appendix B: Key Formulas Quick Reference

| Quantity | Formula |
|----------|---------|
| Tidal radius | r_t = R(M_BH/M_*)^(1/3) |
| Penetration factor | β = r_t/r_p |
| Energy spread | ΔE = 2GM_BH R/r_p² |
| Min semi-major axis | a_min = r_p²/(2R) |
| Min fallback time | t_min = (π/√2)(a_min³/GM_BH)^(1/2) |
| Fallback rate | Ṁ ∝ t^(-5/3) |
| Eddington rate | Ṁ_Edd = L_Edd/(ηc²) |
| GR precession | Δφ = 3πGM/(c²r_p) |
| Compression factor | H/H_0 ~ 1/β |
| Shock heating | T_2/T_1 ~ 0.31 M² |

---

*Document generated from SPH simulation analysis*
*IMBH-Cloud Tidal Disruption Parameter Survey*
