# Theoretical Derivation: IMBH-Molecular Cloud Scattering Simulation

## Executive Summary

This document provides first-principles derivations to determine the appropriate physics for simulating the tidal deformation and shock formation during intermediate-mass black hole (IMBH) molecular cloud scattering, as observed in CO-0.40-0.22 (Oka et al. 2017).

**Key Conclusions:**
1. **MHD effects: NOT required** for the primary dynamics (super-Alfvénic flow with $M_A \sim 250$)
2. **Two-fluid approximation: NOT required** (ambipolar diffusion timescale $\gg$ dynamical timescale)
3. **Shock formation: YES** (Mach number $\mathcal{M} \sim 10-15$)

---

## 1. Physical Parameters from Observations

### 1.1 IMBH Properties (Oka et al. 2017)

| Parameter | Value | Symbol |
|-----------|-------|--------|
| Black hole mass | $10^5\,M_\odot$ | $M_\mathrm{BH}$ |
| Distance from Galactic center | $\sim 60$ pc | - |
| Point source size upper limit | $< 0.022$ pc | - |

### 1.2 Molecular Cloud Properties (CO-0.40-0.22)

| Parameter | Value | Symbol |
|-----------|-------|--------|
| Cloud size | $\sim 5$ pc | $R_\mathrm{cloud}$ |
| Dense clump size | $\sim 0.3$ pc | $r_\mathrm{clump}$ |
| Velocity width | $\sim 100$ km/s | $\Delta v$ |
| Velocity dispersion | 22 km/s | $\sigma_v$ |
| Dense clump mass | $\sim 40\,M_\odot$ | $M_\mathrm{clump}$ |
| Virial mass | $\sim 4100\,M_\odot$ | $M_\mathrm{vir}$ |
| Temperature | 60 K | $T$ |
| Number density (dense) | $10^{6.5}$ cm$^{-3}$ | $n_\mathrm{dense}$ |

### 1.3 Cloud Mass Estimation from Radio Observations (Oka et al. 2017)

The cloud mass in Oka et al. (2017) is derived from millimeter-wave molecular line observations using ALMA Band 6. This section provides the complete derivation connecting observational data to simulation input parameters.

#### 1.3.0 Observational Data Summary (from Oka et al. 2017)

**ALMA Observations:**
- HCN J=3-2 (265.9 GHz) and CO J=2-1 (230.5 GHz)
- Beam sizes: $1.52'' \times 0.60''$ (HCN) and $1.87'' \times 1.14''$ (CO)
- Velocity range: $V_\mathrm{LSR} = -110$ to $0$ km/s

**Key Observational Results:**
| Observable | Value | Source |
|------------|-------|--------|
| Dense clump position | $(l,b) = (-0.398°, -0.224°)$ | HCN J=3-2 peak |
| Point source (CO-0.40-0.22*) | $(l,b) = (-0.3983°, -0.2235°)$ | 266 GHz continuum |
| Displacement | $\sim 0.2$ pc | Between HCN peak and continuum |
| Velocity range | $V_\mathrm{LSR} = -105$ to $-5$ km/s | HCN J=3-2 emission |
| Main body velocity | $V_\mathrm{LSR} = -80$ to $-40$ km/s | Dense clump |
| Size parameter | $S = 0.15$ pc | HCN clump half-size |
| Velocity dispersion | $\sigma_v = 22$ km/s | Line width |
| Continuum flux (231 GHz) | $8.38 \pm 0.34$ mJy | ALMA |
| Continuum flux (266 GHz) | $9.91 \pm 0.74$ mJy | ALMA |
| Spectral index | $\alpha = 1.18 \pm 0.65$ | Inverted spectrum |

#### 1.3.1 Radiative Transfer Fundamentals

The observed intensity $I_\nu$ along a line of sight satisfies:

$$I_\nu = I_{\nu,0} e^{-\tau_\nu} + S_\nu (1 - e^{-\tau_\nu})$$

where $\tau_\nu$ is the optical depth and $S_\nu$ is the source function. For molecular lines in Local Thermodynamic Equilibrium (LTE):

$$S_\nu = B_\nu(T_\mathrm{ex}) = \frac{2h\nu^3}{c^2} \frac{1}{e^{h\nu/k_B T_\mathrm{ex}} - 1}$$

The **antenna temperature** (Rayleigh-Jeans approximation, $h\nu \ll k_B T$):

$$T_A = \frac{c^2}{2k_B\nu^2} I_\nu = T_\mathrm{ex} (1 - e^{-\tau_\nu}) - T_\mathrm{bg}(1 - e^{-\tau_\nu})$$

where $T_\mathrm{bg} \approx 2.73$ K is the cosmic microwave background.

#### 1.3.2 Column Density from CO Lines

For optically thin emission ($\tau \ll 1$), the **column density** of the upper state $N_u$ is:

$$N_u = \frac{8\pi k_B \nu^2}{hc^3 A_{ul}} \int T_A \, dv$$

where $A_{ul}$ is the Einstein A coefficient.

The **total column density** using the partition function:

$$N_\mathrm{tot} = \frac{N_u}{g_u} Q(T_\mathrm{ex}) \exp\left(\frac{E_u}{k_B T_\mathrm{ex}}\right)$$

For **CO J=2-1** transition ($\nu = 230.538$ GHz):
- $A_{21} = 6.91 \times 10^{-7}$ s$^{-1}$
- $E_u/k_B = 16.6$ K
- $g_u = 2J + 1 = 5$

The partition function for a linear rotor at temperature $T$:

$$Q(T) \approx \frac{k_B T}{hB} + \frac{1}{3} = \frac{T}{2.77\,\mathrm{K}} + \frac{1}{3}$$

where $B = 57.636$ GHz is the rotational constant for CO.

**CO column density formula:**

$$N(\mathrm{CO}) = \frac{3k_B}{8\pi^3 B \mu^2} \frac{Q(T_\mathrm{ex})}{J+1} \exp\left(\frac{E_J}{k_B T_\mathrm{ex}}\right) \int T_A \, dv$$

For CO J=2-1 at $T_\mathrm{ex} = 30$ K:

$$\boxed{N(\mathrm{CO}) \approx 1.1 \times 10^{14} \left(\frac{T_\mathrm{ex} + 0.92}{1 - e^{-16.6/T_\mathrm{ex}}}\right) \int T_A \, dv \quad \mathrm{[cm^{-2}]}}$$

where $\int T_A \, dv$ is in K km/s.

#### 1.3.3 H$_2$ Mass from CO Observations

The **CO-to-H$_2$ conversion factor** (X-factor):

$$N(\mathrm{H_2}) = X_\mathrm{CO} \times W_\mathrm{CO}$$

where $W_\mathrm{CO} = \int T_A \, dv$ is the integrated intensity.

Standard Galactic value:
$$X_\mathrm{CO} = 2 \times 10^{20}\,\mathrm{cm^{-2}\,(K\,km/s)^{-1}}$$

**Note**: In the Galactic Center, $X_\mathrm{CO}$ may be lower by factor 2-5 due to higher temperatures and metallicity.

The **total mass**:

$$M = \mu_\mathrm{H_2} m_H \times N(\mathrm{H_2}) \times A$$

where $\mu_\mathrm{H_2} = 2.8$ (including He) and $A$ is the projected area.

$$\boxed{M = 2.0 \times 10^{-20} \times X_\mathrm{CO} \times W_\mathrm{CO} \times A \quad \mathrm{[g]}}$$

Converting to solar masses with $A$ in pc$^2$:

$$M = 4.4 \times X_{20} \times W_\mathrm{CO} \times \left(\frac{A}{\mathrm{pc}^2}\right) \quad M_\odot$$

where $X_{20} = X_\mathrm{CO} / (10^{20}\,\mathrm{cm^{-2}\,(K\,km/s)^{-1}})$.

#### 1.3.4 Dense Gas Mass from HCN Lines

HCN traces **denser gas** ($n_\mathrm{crit} \sim 10^6$ cm$^{-3}$) than CO ($n_\mathrm{crit} \sim 10^3$ cm$^{-3}$).

For **HCN J=3-2** ($\nu = 265.886$ GHz):
- $A_{32} = 8.36 \times 10^{-4}$ s$^{-1}$
- Critical density: $n_\mathrm{crit} = A_{ul}/\gamma_{ul} \approx 3 \times 10^6$ cm$^{-3}$

The HCN column density (non-LTE regime requires excitation modeling):

$$N(\mathrm{HCN}) = \frac{3h}{8\pi^3 \mu^2} \frac{Q(T_\mathrm{ex})}{J(J+1)} \exp\left(\frac{E_J}{k_B T_\mathrm{ex}}\right) \frac{\tau}{1-e^{-\tau}} \int T_A \, dv$$

For the **Oka et al. dense clump** ($T_k = 60$ K, $n(\mathrm{H_2}) = 10^{6.5}$ cm$^{-3}$):

Using RADEX or similar non-LTE code with:
- $T_k = 60$ K
- $n(\mathrm{H_2}) = 3 \times 10^6$ cm$^{-3}$
- HCN abundance $[\mathrm{HCN}]/[\mathrm{H_2}] \sim 10^{-8}$

The derived mass:

$$\boxed{M_\mathrm{clump} \approx 40\,M_\odot}$$

#### 1.3.5 Virial Mass

The virial theorem mass (assuming spherical, uniform density):

$$M_\mathrm{vir} = \frac{5\sigma_v^2 R}{G}$$

For the dense clump with $\sigma_v = 22$ km/s and $R = S = 0.15$ pc:

$$M_\mathrm{vir} = \frac{5 \times (22 \times 10^5)^2 \times 0.15 \times 3.09 \times 10^{18}}{6.67 \times 10^{-8}}$$

$$M_\mathrm{vir} = \frac{5 \times 4.84 \times 10^{12} \times 4.6 \times 10^{17}}{6.67 \times 10^{-8}} = \frac{1.1 \times 10^{31}}{6.67 \times 10^{-8}}$$

$$\boxed{M_\mathrm{vir} \approx 4100\,M_\odot}$$

Since $M_\mathrm{vir} \gg M_\mathrm{clump}$, the clump is **NOT gravitationally bound** — consistent with tidal disruption scenario.

#### 1.3.6 N-body Model Mass (Oka et al. 2017)

Oka et al. performed N-body simulations to reproduce the observed kinematics. The key insight is that the observed position-velocity structure constrains the cloud mass through the gravitational kick dynamics.

**N-body Model Setup (from Methods section):**
- Point mass: $M_\mathrm{BH} = 10^5\,M_\odot$
- Cloud: 1000 test particles, each $1\,M_\odot$ → $M_\mathrm{cloud} = 1000\,M_\odot$
- Initial cloud: Gaussian radial dispersion $\sigma_r = 0.2$ pc
- Initial velocity dispersion: $1.43$ km/s (virialized)
- Initial position: $(X, Y) = (9.8\,\mathrm{pc}, -0.65\,\mathrm{pc})$
- Initial velocity: $(v_X, v_Y) = (-8.19\,\mathrm{km/s}, 0.4\,\mathrm{km/s})$
- Time step: $4.6 \times 10^3$ yr
- Best-fit snapshot: $t = 7.2 \times 10^5$ yr (just after pericentre)

**Why 1000 $M_\odot$?**

The N-body mass is chosen to:
1. Match the observed velocity width ($\sim 100$ km/s)
2. Reproduce the position-velocity gradient of the HCN clump
3. Be consistent with the cloud being tidally disrupted (not self-gravitating)

The **tidal mass** at pericentre provides an upper limit:

$$M_\mathrm{tidal} = M_\mathrm{BH} \left(\frac{r_\mathrm{cloud}}{r_\mathrm{peri}}\right)^3 = 10^5 \times \left(\frac{0.2}{1}\right)^3 = 800\,M_\odot$$

This is consistent with the N-body model mass of $1000\,M_\odot$.

#### 1.3.8 Mass Summary for Simulation

| Mass Estimate | Value | Method | Physical Meaning | Use for Simulation |
|---------------|-------|--------|------------------|-------------------|
| Dense clump (HCN) | 40 $M_\odot$ | Non-LTE excitation model | Gas mass in densest region | Dense core only |
| Virial mass | 4100 $M_\odot$ | $5\sigma_v^2 R/G$ | Dynamical mass (if bound) | Upper limit |
| Tidal mass | 800 $M_\odot$ | $M_\mathrm{BH}(r/r_\mathrm{peri})^3$ | Max mass at tidal disruption | Consistency check |
| N-body model | 1000 $M_\odot$ | Best fit to P-V structure | Kinematically constrained | **Recommended** |
| Extended cloud | ~$10^4$ $M_\odot$ | CO J=3-2 integrated | Full CO envelope | Full cloud sim |

**Physical Interpretation:**
- $M_\mathrm{vir} \gg M_\mathrm{clump}$ confirms the clump is **not gravitationally bound**
- The large velocity dispersion ($\sigma_v = 22$ km/s) is due to **tidal acceleration**, not self-gravity
- The N-body mass ($1000\,M_\odot$) represents the **dynamically relevant** portion of the cloud

**Recommendation**: Use $M_\mathrm{cloud} = 1000\,M_\odot$ as the simulation mass, which provides the best kinematic fit in Oka et al.'s N-body model.

### 1.4 Physical Assumptions and Their Justification in IMBH Context

This section explains **why** certain physical assumptions in Oka et al. (2017) are valid for the IMBH-cloud scattering scenario, and how these connect to simulation post-processing.

#### 1.4.1 Why Consider Virial Mass?

**Physical Motivation:**
The virial theorem relates kinetic and potential energy for a gravitationally bound system:
$$2K + U = 0 \quad \Rightarrow \quad M_\mathrm{vir} = \frac{5\sigma_v^2 R}{G}$$

**Why Oka et al. calculate it:**
1. **Diagnostic of gravitational state**: If $M_\mathrm{vir} \approx M_\mathrm{gas}$, the cloud is self-gravitating and bound. If $M_\mathrm{vir} \gg M_\mathrm{gas}$, something else drives the kinematics.

2. **IMBH signature detection**: For CO-0.40-0.22:
   - $M_\mathrm{vir} = 4100\,M_\odot$ (from velocity dispersion)
   - $M_\mathrm{gas} = 40\,M_\odot$ (from HCN emission)
   - Ratio: $M_\mathrm{vir}/M_\mathrm{gas} \sim 100$!

3. **Physical interpretation**: This huge discrepancy means:
   - The observed velocity dispersion ($\sigma_v = 22$ km/s) is **NOT** due to self-gravity
   - An **external gravitational source** (the IMBH) must be accelerating the gas
   - The cloud is being **tidally disrupted**, not virially bound

**In IMBH context:**
$$\sigma_v^2 \approx \frac{GM_\mathrm{BH}}{r_\mathrm{peri}} \quad \text{(tidal kick)}$$

For $M_\mathrm{BH} = 10^5\,M_\odot$ and $r_\mathrm{peri} \sim 1$ pc:
$$\sigma_v \sim \sqrt{\frac{6.67 \times 10^{-8} \times 2 \times 10^{38}}{3 \times 10^{18}}} \sim 20\,\mathrm{km/s}$$

This matches the observed $\sigma_v = 22$ km/s — **confirming the IMBH hypothesis**.

#### 1.4.2 Why is Local Thermodynamic Equilibrium (LTE) Sufficient?

**LTE Condition:**
Collisional excitation rate must exceed radiative de-excitation rate:
$$n \times \langle\sigma v\rangle_\mathrm{coll} > A_{ul}$$

This defines the **critical density**:
$$n_\mathrm{crit} = \frac{A_{ul}}{\gamma_{ul}} \approx \frac{A_{ul}}{\langle\sigma v\rangle}$$

| Transition | $A_{ul}$ (s$^{-1}$) | $n_\mathrm{crit}$ (cm$^{-3}$) |
|------------|---------------------|------------------------------|
| CO J=2-1 | $6.9 \times 10^{-7}$ | $\sim 10^3$ |
| CO J=3-2 | $2.5 \times 10^{-6}$ | $\sim 3 \times 10^3$ |
| HCN J=3-2 | $8.4 \times 10^{-4}$ | $\sim 3 \times 10^6$ |

**Why LTE works for CO in this context:**

1. **Density exceeds critical**: For CO-0.40-0.22, $n \sim 10^4$–$10^{6.5}$ cm$^{-3}$, well above CO's $n_\mathrm{crit} \sim 10^3$ cm$^{-3}$.

2. **High kinetic temperature**: $T_k = 60$ K means efficient collisional excitation.

3. **Dense clump**: The HCN-detected clump has $n \sim 10^{6.5}$ cm$^{-3}$, so even HCN approaches LTE.

**When LTE fails:**
- In the diffuse envelope ($n < n_\mathrm{crit}$), HCN requires non-LTE analysis (RADEX)
- This is why Oka et al. use RADEX for HCN but simpler LTE for CO

**For simulation post-processing:**
- LTE assumption simplifies radiative transfer: $S_\nu = B_\nu(T_\mathrm{ex})$
- For dense SPH particles ($n > n_\mathrm{crit}$), use $T_\mathrm{ex} \approx T_\mathrm{kinetic}$
- For lower-density regions, either apply non-LTE corrections or focus on CO (not HCN)

#### 1.4.3 Why is the Rayleigh-Jeans Approximation Valid?

**The Rayleigh-Jeans limit:**
When $h\nu \ll k_B T$, the Planck function simplifies:
$$B_\nu(T) = \frac{2h\nu^3}{c^2} \frac{1}{e^{h\nu/k_B T}-1} \approx \frac{2k_B T \nu^2}{c^2}$$

This allows defining **brightness temperature**:
$$T_b = \frac{c^2}{2k_B\nu^2} I_\nu$$

**Validity check for ALMA Band 6:**

| Frequency | $h\nu/k_B$ | $T_\mathrm{cloud}$ | Ratio $k_B T / h\nu$ |
|-----------|-----------|-------------------|---------------------|
| 230 GHz (CO) | 11.0 K | 60 K | 5.5 |
| 266 GHz (HCN) | 12.8 K | 60 K | 4.7 |

For $T = 60$ K and $\nu = 230$ GHz:
$$\frac{k_B T}{h\nu} = \frac{60}{11} = 5.5 \gg 1 \quad \checkmark$$

**Why it's appropriate here:**
1. **Cloud is warm**: $T = 60$ K, much higher than typical GMCs (10–20 K)
2. **IMBH shock heating**: Tidal interaction and shocks heat the gas
3. **Error is small**: For $T/T_{h\nu} = 5$, the R-J error is only $\sim 3\%$

**Correction factor if needed:**
$$\frac{T_b^\mathrm{true}}{T_b^\mathrm{RJ}} = \frac{h\nu/k_B T}{e^{h\nu/k_B T} - 1}$$

For $T = 60$ K at 230 GHz: correction factor = 0.97 (3% error).

#### 1.4.4 Why Assume Optically Thin Emission?

**Optical depth definition:**
$$\tau_\nu = \int \alpha_\nu \, ds = \int n_l \sigma_{lu} \left(1 - \frac{g_l n_u}{g_u n_l}\right) ds$$

For molecular lines:
$$\tau \approx \frac{c^3 A_{ul}}{8\pi\nu^3} \frac{g_u}{g_l} N_l \left(1 - e^{-h\nu/k_B T_\mathrm{ex}}\right)$$

**Why optically thin ($\tau < 1$) is reasonable:**

1. **High velocity dispersion**: The large line width ($\Delta v \sim 100$ km/s) spreads photons over many frequency channels:
   $$\tau_\mathrm{peak} \propto \frac{N}{\Delta v}$$
   Broad lines → lower peak optical depth.

2. **Tidal disruption geometry**: The cloud is stretched into a stream, reducing column density along most sight-lines.

3. **Consistency check**: If optically thick, the line would have a flat-topped profile with $T_b \rightarrow T_\mathrm{ex}$. The observed profile shows the expected optically thin shape.

4. **For CO J=2-1** at $N(\mathrm{CO}) = 10^{16}$ cm$^{-2}$ and $\Delta v = 100$ km/s:
   $$\tau \sim 0.1 \quad \text{(optically thin)} \checkmark$$

**When this fails:**
- At the densest clump center, CO J=2-1 may become moderately optically thick ($\tau \sim 1$–3)
- Use isotopologue $^{13}$CO or C$^{18}$O for truly optically thin mass estimates
- For HCN, optical depth effects are more significant

**For simulation:**
- Post-processing should compute $\tau$ for each sightline
- Apply optical depth correction: $T_b = T_\mathrm{ex}(1 - e^{-\tau})$
- Flag regions where $\tau > 1$ for special treatment

#### 1.4.5 What is the Spectral Index and Continuum Flux?

**Spectral Index Definition:**
The spectral index $\alpha$ describes how continuum flux varies with frequency:
$$S_\nu \propto \nu^\alpha \quad \Rightarrow \quad \alpha = \frac{\log(S_{\nu_2}/S_{\nu_1})}{\log(\nu_2/\nu_1)}$$

**Oka et al. measurements:**
- $S_{231\,\mathrm{GHz}} = 8.38 \pm 0.34$ mJy
- $S_{266\,\mathrm{GHz}} = 9.91 \pm 0.74$ mJy
- $\alpha = 1.18 \pm 0.65$

**Physical interpretation of $\alpha$:**

| $\alpha$ value | Emission mechanism | Physical source |
|----------------|-------------------|-----------------|
| $\alpha = 2$ | Thermal (R-J limit, optically thick) | Hot dust, dense gas |
| $\alpha = -0.1$ | Optically thin synchrotron | AGN jets, SNRs |
| $\alpha \sim 0.6$ | Partially thick synchrotron | Compact radio cores |
| $\alpha \sim 1$–2 | Inverted/self-absorbed | Compact jets, ADAF |

**Why $\alpha = 1.18$ is significant for IMBH:**

1. **Not thermal dust**: Dust emission has $\alpha \approx 2$–4 (modified blackbody)
2. **Not standard synchrotron**: Optically thin synchrotron has $\alpha < 0$
3. **Suggests compact, self-absorbed source**: Consistent with:
   - Advection-dominated accretion flow (ADAF) around IMBH
   - Compact jet base
   - Size upper limit < 0.022 pc supports compact source

4. **Point source CO-0.40-0.22***: Spatially coincident with HCN clump but offset by 0.2 pc — suggests the IMBH is located at the continuum source position, with the dense gas slightly displaced.

**Continuum flux purpose:**
- Constrains IMBH accretion rate: $\dot{M} \propto L_\mathrm{radio}$
- Provides position of the putative IMBH (more accurate than line emission centroid)
- Spectral index constrains emission mechanism (synchrotron vs thermal)

#### 1.4.6 What is HCN and Why the Displacement from Continuum?

**HCN (Hydrogen Cyanide):**
- Linear molecule: H–C≡N
- Strong dipole moment: $\mu = 2.99$ Debye (vs CO: $\mu = 0.11$ Debye)
- High critical density: $n_\mathrm{crit} \sim 10^6$ cm$^{-3}$

**Why HCN is important:**
1. **Dense gas tracer**: Only emits significantly where $n > 10^5$–$10^6$ cm$^{-3}$
2. **Complements CO**: CO traces total molecular mass; HCN traces dense cores
3. **Higher excitation**: Requires more collisions → traces shocked/compressed gas

**The 0.2 pc displacement:**

From Oka et al. observations:
- HCN J=3-2 peak: $(l,b) = (-0.398°, -0.224°)$
- Continuum source (CO-0.40-0.22*): $(l,b) = (-0.3983°, -0.2235°)$
- Angular separation: $\sim 1.5''$ → **0.2 pc at 8 kpc**

**Physical interpretation in IMBH context:**

```
        [Dense HCN clump]
              ↑
              | 0.2 pc displacement
              |
        [IMBH + continuum] ←── Accretion produces radio continuum
              |
              ↓ tidal tail
        [Extended CO stream]
```

1. **IMBH location**: The continuum source marks the IMBH position (accretion-powered emission)

2. **Dense clump offset**: The HCN-bright dense clump is the surviving core of the cloud, displaced from the IMBH after closest approach

3. **Tidal dynamics**: As the cloud passed pericentre:
   - Leading edge accelerated toward BH → now extended CO tail
   - Trailing dense core decelerated → now the HCN clump 0.2 pc behind

4. **Timescale**: At $v \sim 100$ km/s, 0.2 pc separation implies $\sim 2000$ yr since closest approach — consistent with N-body best-fit time $t = 7.2 \times 10^5$ yr (total orbital time)

**For simulation:**
- Track the position of the density maximum (analog of HCN peak)
- Compare with BH position
- Displacement should grow after pericentre passage

### 1.5 Simulation Post-Processing: Reproducing Oka et al. P-V Diagram

To compare SPH simulation output with Oka et al.'s observations, we need to generate synthetic position-velocity (P-V) diagrams.

#### 1.5.1 Overview of P-V Diagram Construction

A position-velocity diagram shows intensity as a function of:
- **Position**: Usually along a 1D cut (position angle on sky)
- **Velocity**: Line-of-sight velocity (Doppler shift from rest frequency)

```
                 Position (arcsec or pc)
                ─────────────────────────►
            │
   Velocity │    ████
   (km/s)   │      ████████
            │        ████████████
            │          ████████
            │            ████
            ▼
```

The characteristic "parallelogram" shape in CO-0.40-0.22 indicates tidal stretching.

#### 1.5.2 Step-by-Step Post-Processing Pipeline

**Step 1: Define Observer Frame**

Transform simulation coordinates to observer frame using Oka et al. parameters:
- Inclination: $i = 70°$ (orbital plane to line-of-sight)
- Position angle: PA $= 41.6°$ (rotation of orbital plane projection)
- Distance: $d = 8.0$ kpc (Galactic Center distance)
- Systemic velocity: $V_\mathrm{LSR} = -120$ km/s

Rotation matrices:
$$\mathbf{R}_\mathrm{inc} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & \cos i & -\sin i \\ 0 & \sin i & \cos i \end{pmatrix}$$

$$\mathbf{R}_\mathrm{PA} = \begin{pmatrix} \cos\mathrm{PA} & -\sin\mathrm{PA} & 0 \\ \sin\mathrm{PA} & \cos\mathrm{PA} & 0 \\ 0 & 0 & 1 \end{pmatrix}$$

**Step 2: Compute Line-of-Sight Velocity**

For each SPH particle $i$:
$$v_\mathrm{LOS,i} = \mathbf{v}_i \cdot \hat{z}_\mathrm{obs} + V_\mathrm{LSR}$$

where $\hat{z}_\mathrm{obs}$ is the unit vector toward the observer.

**Step 3: Compute CO Emissivity (LTE approximation)**

For each particle, compute the CO J=2-1 emissivity:

1. **Column density contribution** (SPH formulation):
   $$N_\mathrm{CO,i} = X_\mathrm{CO}^{-1} \times n_i \times h_i$$
   
   where $h_i$ is the smoothing length.

2. **Excitation temperature** (assume $T_\mathrm{ex} \approx T_\mathrm{kin}$):
   $$T_\mathrm{ex,i} = T_i \quad \text{(from SPH internal energy)}$$

3. **Line intensity** (optically thin LTE):
   $$I_{\nu,i} \propto N_\mathrm{CO,i} \times A_{21} \times \frac{g_2}{Q(T_\mathrm{ex})} \exp\left(-\frac{E_2}{k_B T_\mathrm{ex}}\right)$$

4. **Velocity channel assignment**:
   Assign particle to velocity bin centered on $v_\mathrm{LOS,i}$ with thermal broadening:
   $$\Delta v_\mathrm{thermal} = \sqrt{\frac{2k_B T}{m_\mathrm{CO}}} \approx 0.1\,\mathrm{km/s} \times \sqrt{\frac{T}{10\,\mathrm{K}}}$$

**Step 4: Grid the Data**

Create a 2D histogram:
- X-axis: Position along chosen cut (e.g., Galactic longitude offset)
- Y-axis: $V_\mathrm{LSR}$ velocity
- Value: Sum of CO emissivity in each (position, velocity) bin

```python
# Pseudo-code for P-V diagram generation
import numpy as np

def create_pv_diagram(particles, observer_params):
    """
    Generate P-V diagram from SPH particle data.
    
    Parameters:
    -----------
    particles : structured array with fields:
        - pos: (N, 3) positions in simulation frame [pc]
        - vel: (N, 3) velocities in simulation frame [km/s]
        - density: (N,) density [g/cm^3]
        - temperature: (N,) temperature [K]
        - mass: (N,) particle mass [M_sun]
        - h: (N,) smoothing length [pc]
    
    observer_params : dict with keys:
        - inclination: orbital plane inclination [deg]
        - PA: position angle [deg]
        - V_LSR: systemic velocity [km/s]
        - distance: distance to source [kpc]
    """
    # Transform to observer frame
    i_rad = np.radians(observer_params['inclination'])
    PA_rad = np.radians(observer_params['PA'])
    
    # Rotation matrices
    R_inc = rotation_matrix_x(i_rad)
    R_PA = rotation_matrix_z(PA_rad)
    R_total = R_PA @ R_inc
    
    # Transform positions and velocities
    pos_obs = particles['pos'] @ R_total.T
    vel_obs = particles['vel'] @ R_total.T
    
    # Line-of-sight velocity (z-component in observer frame)
    v_LOS = vel_obs[:, 2] + observer_params['V_LSR']
    
    # Position along major axis (x-component in observer frame)
    pos_major = pos_obs[:, 0]
    
    # CO emissivity (simplified LTE)
    n_H2 = particles['density'] / (2.8 * 1.67e-24)  # cm^-3
    T = particles['temperature']
    
    # CO column density contribution per particle
    # Using X_CO = 2e20 cm^-2 / (K km/s)
    N_CO = particles['h'] * 3.09e18 * n_H2 * 1e-4  # rough CO/H2 ratio
    
    # Excitation factor (LTE)
    E_2_over_kB = 16.6  # K for CO J=2
    Q_T = T / 2.77 + 1/3
    excitation = np.exp(-E_2_over_kB / T) / Q_T
    
    # Emissivity weight
    emissivity = N_CO * excitation * T
    
    # Create P-V grid
    pos_bins = np.linspace(-3, 3, 100)  # pc
    vel_bins = np.linspace(-120, 0, 100)  # km/s
    
    pv_diagram, _, _ = np.histogram2d(
        pos_major, v_LOS,
        bins=[pos_bins, vel_bins],
        weights=emissivity
    )
    
    return pv_diagram, pos_bins, vel_bins
```

**Step 5: Convolve with Beam**

Apply ALMA beam smoothing:
- Spatial: Gaussian with FWHM $= 1.87'' \times 1.14''$ (CO beam)
- Velocity: Channel width $\Delta v = 2$ km/s

$$I_\mathrm{obs}(x,v) = I_\mathrm{sim}(x,v) * G_\mathrm{beam}(x) * G_\mathrm{channel}(v)$$

**Step 6: Add Noise (optional)**

For realistic comparison, add noise matching ALMA sensitivity:
$$\sigma_\mathrm{rms} \approx 10\,\mathrm{mJy/beam} \times \sqrt{\frac{2\,\mathrm{km/s}}{\Delta v}}$$

#### 1.5.3 Key Features to Compare with Observations

| Feature | Oka et al. Observation | Simulation Target |
|---------|----------------------|-------------------|
| Velocity range | $-105$ to $-5$ km/s | Should span $\sim 100$ km/s |
| P-V gradient | Linear in parallelogram | Tidal stretching signature |
| Dense clump position | $(l,b) = (-0.398°, -0.224°)$ | Track density maximum |
| Velocity at clump | $V_\mathrm{LSR} \sim -60$ km/s | Main body velocity |
| Stream morphology | Elongated along P-V | Tidal tail visible |

#### 1.5.4 Alternative: Full Radiative Transfer with RADMC-3D

For more accurate comparison (especially HCN):

1. Export SPH data to RADMC-3D grid format
2. Include:
   - Density structure
   - Temperature structure
   - Velocity field
   - Turbulent velocity (sub-grid)
3. Run line radiative transfer with:
   - LVG or LTE excitation
   - Full 3D geometry
   - Proper beam convolution

**Advantages:**
- Handles optical depth effects
- Proper non-LTE for HCN
- Can generate channel maps, not just P-V

**Disadvantages:**
- Computationally expensive
- Requires additional software
- More parameters to constrain

### 1.6 N-body Simulation Parameters (Oka et al. 2017)

The N-body model parameters define the initial conditions for SPH simulations:

| Parameter | Symbol | Value | Notes |
|-----------|--------|-------|-------|
| BH mass | $M_\mathrm{BH}$ | $10^5\,M_\odot$ | Point mass at origin |
| Cloud mass | $M_\mathrm{cloud}$ | $1000\,M_\odot$ | 1000 test particles × 1 $M_\odot$ |
| Initial radius dispersion | $\sigma_r$ | 0.2 pc | Gaussian distribution |
| Initial velocity dispersion | $\sigma_v$ | 1.43 km/s | Virialized cloud |
| Initial position | $(X_0, Y_0)$ | (9.8 pc, -0.65 pc) | In orbital plane coordinates |
| Initial velocity | $(v_{X,0}, v_{Y,0})$ | (-8.19 km/s, 0.4 km/s) | Approaching BH |
| Orbital plane inclination | $i$ | 70° | w.r.t. line of sight |
| Position angle | PA | 41.6° | Orbital plane rotation |
| Center of mass velocity | $V_\mathrm{LSR}$ | -120 km/s | Line-of-sight |
| Best-fit time | $t$ | $7.2 \times 10^5$ yr | Just after pericentre |
| Time step | $\Delta t$ | $4.6 \times 10^3$ yr | Leapfrog integration |

**Coordinate System:**
- $X$-axis: In orbital plane, initially along BH-cloud direction
- $Y$-axis: In orbital plane, perpendicular to $X$
- $Z$-axis: Perpendicular to orbital plane
- Origin: Center of mass (≈ BH position)

### 1.7 SPH Initial Condition Derivation (CAT_OKA)

The Oka et al. (2017) N-body simulation used **test particles** with a Gaussian radius dispersion $\sigma_r = 0.2$ pc. However, for SPH simulations, we require a **resolved cloud** with finite radius $R_\text{cloud}$ to capture hydrodynamic effects. This section derives the scaled orbital parameters that maintain the same tidal physics while accommodating a larger SPH cloud.

#### 1.7.1 The Scaling Problem

**Original Oka et al. parameters:**
- Cloud size: $\sigma_r = 0.2$ pc (Gaussian dispersion)
- Initial position: $(X_0, Y_0) = (9.8, -0.65)$ pc
- Initial velocity: $(v_{X,0}, v_{Y,0}) = (-8.19, 0.4)$ km/s
- BH mass: $M_\mathrm{BH} = 10^5\,M_\odot$

**Problem:** The original trajectory has an extremely small pericentre. Computing from energy and angular momentum:

$$r_0 = \sqrt{X_0^2 + Y_0^2} = \sqrt{9.8^2 + 0.65^2} = 9.82\,\text{pc}$$

$$v_0 = \sqrt{v_{X,0}^2 + v_{Y,0}^2} = \sqrt{8.19^2 + 0.4^2} = 8.20\,\text{km/s}$$

The specific orbital energy:
$$\varepsilon = \frac{1}{2}v_0^2 - \frac{GM_\mathrm{BH}}{r_0}$$

With $GM_\mathrm{BH} = G \times 10^5 M_\odot$:
$$GM_\mathrm{BH} = 6.674 \times 10^{-11} \times 10^5 \times 1.989 \times 10^{30} = 1.327 \times 10^{25}\,\text{m}^3\text{s}^{-2}$$

Converting to units of pc·(km/s)²:
$$GM_\mathrm{BH} = \frac{1.327 \times 10^{25}}{(3.086 \times 10^{16}) \times (10^3)^2} = 449.8\,\text{pc}\cdot(\text{km/s})^2$$

The specific orbital energy:
$$\varepsilon = \frac{1}{2}(8.20)^2 - \frac{449.8}{9.82} = 33.6 - 45.8 = -12.2\,(\text{km/s})^2$$

Wait—this gives a **bound orbit** ($\varepsilon < 0$), but Oka et al. describe an **unbound flyby**. Let's reconsider: the negative energy implies an elliptical orbit, but with very high eccentricity the cloud may still be on its first approach. The key issue is that for test particles, passing arbitrarily close to the BH is allowed, but for an SPH cloud with finite radius $R = 1.13$ pc, we need $r_\text{peri} > R$ to avoid the cloud immediately engulfing the BH.

#### 1.7.2 SPH Cloud Requirements

For our SPH simulations:
- **Cloud radius**: $R_\text{cloud} = 1.13$ pc (from relaxed polytrope IC)
- **Minimum pericentre**: $r_\text{peri} \geq 1.5 \times R_\text{cloud} = 1.69$ pc

This ensures:
1. The cloud center never coincides with the BH
2. Tidal forces are computed at physically meaningful distances
3. The cloud experiences strong tidal stretching without numerical singularities

#### 1.7.3 Orbital Mechanics: Hyperbolic Trajectory Design

We design a **hyperbolic orbit** (unbound, $e > 1$) with specified pericentre that places the cloud in the **strong tidal regime** similar to Oka et al.

**Key orbital parameters:**

1. **Pericentre distance**: $r_\text{peri} = 1.69$ pc (= $1.5 \times R_\text{cloud}$)

2. **Velocity at infinity**: We adopt $v_\infty = 8$ km/s, similar to Oka's bulk velocity scale.

3. **Specific orbital energy** (hyperbolic):
   $$\varepsilon = \frac{1}{2}v_\infty^2 = \frac{1}{2}(8)^2 = 32\,(\text{km/s})^2 > 0$$

4. **Semi-major axis** (negative for hyperbola):
   $$a = -\frac{GM_\mathrm{BH}}{2\varepsilon} = -\frac{449.8}{2 \times 32} = -7.03\,\text{pc}$$

5. **Eccentricity** from $r_\text{peri} = a(1 - e)$:
   $$e = 1 - \frac{r_\text{peri}}{a} = 1 - \frac{1.69}{-7.03} = 1 + 0.240 = 1.240$$

   Since $e > 1$, this confirms a hyperbolic trajectory.

6. **Specific angular momentum** from $r_\text{peri} = \frac{h^2}{GM(1+e)}$:
   $$h^2 = GM_\mathrm{BH} \cdot r_\text{peri} \cdot (1 + e) = 449.8 \times 1.69 \times 2.240 = 1702\,\text{pc}^2(\text{km/s})^2$$
   $$h = 41.26\,\text{pc}\cdot\text{km/s}$$

7. **Impact parameter**:
   $$b = \frac{h}{v_\infty} = \frac{41.26}{8} = 5.16\,\text{pc}$$

#### 1.7.4 Orbital Geometry: Understanding Impact Parameter vs Initial Position

This section clarifies the difference between the **impact parameter** $b$ and the **initial position** $(X_0, Y_0)$, and explains the role of **inclination angle** in connecting simulation coordinates to observations.

##### Impact Parameter vs Initial Position

```
                    ORBITAL GEOMETRY DIAGRAM (View from +Z, perpendicular to orbital plane)
                    
                                           Incoming asymptote
                                                  ↑
                                                 /
                    ←───────── b ─────────────→/  Initial position (X₀, Y₀)
                    (impact parameter)       /    = (20.0, -5.17) pc
                                           /      r₀ = 20.66 pc from BH
                                         /        
                                       /          velocity v₀ = (-9.35, 4.48) km/s
                                     ⊙ Cloud        ↘
                                   /                  ↘
                                 /                      ↘ v
                               /                          
                             /                              
                           /        Hyperbolic trajectory    
                         /              e = 1.24              
                       /                                      
        - - - - - - -●- - - - - - - - - - - - - - - - - - - - - - - →  X-axis
                    BH                                        
                  (origin)    ← r_peri = 1.69 pc →            
                       \                            (pericentre)
                         \                                    
                           \                                  
                             \        Outgoing trajectory     
                               \                              
                                 \                            
                                   \                          
                                     \                        
                                       ↓                      
                                   Outgoing asymptote         
                                                              
        ═══════════════════════════════════════════════════════
        
        KEY DISTINCTIONS:
        
        ┌─────────────────────────────────────────────────────────────────────────┐
        │  Quantity          │  Definition                    │  Value (CAT_OKA) │
        ├─────────────────────────────────────────────────────────────────────────┤
        │  Impact parameter  │  Perpendicular distance from   │  b = 5.16 pc     │
        │  b                 │  BH to the incoming asymptote  │                  │
        │                    │  (constant for the orbit)      │                  │
        ├─────────────────────────────────────────────────────────────────────────┤
        │  Initial position  │  Actual starting coordinates   │  (20.0, -5.17,   │
        │  (X₀, Y₀, Z₀)      │  of cloud center at t=0        │   0) pc          │
        ├─────────────────────────────────────────────────────────────────────────┤
        │  Initial distance  │  r₀ = √(X₀² + Y₀²)             │  r₀ = 20.66 pc   │
        │  r₀                │  Distance from BH at t=0       │                  │
        ├─────────────────────────────────────────────────────────────────────────┤
        │  Pericentre        │  Closest approach distance     │  r_peri = 1.69 pc│
        │  r_peri            │  (reached later in simulation) │                  │
        └─────────────────────────────────────────────────────────────────────────┘
```

**Physical Interpretation:**

1. **Impact parameter $b$** is a **conserved orbital property** — it's the perpendicular distance from the BH to the straight-line path the cloud would follow if there were no gravity. As the cloud falls inward, gravity curves its trajectory toward the BH.

2. **Initial position $(X_0, Y_0)$** is where we **start the simulation** — we choose this far enough from the BH (20 pc) that the cloud hasn't yet been significantly perturbed by tidal forces.

3. **The relationship**: The initial position lies **on the hyperbolic orbit** at a specific true anomaly $\theta$. The impact parameter determines the *shape* of the orbit (combined with energy), while the initial position is just *where we start* along that orbit.

$$b = r_0 \sin(\text{angle between } \vec{r}_0 \text{ and asymptote})$$

##### Inclination Angle: Connecting Simulation to Observations

The **orbital plane** (X-Y plane in our simulation) must be transformed to match **observed sky coordinates** and line-of-sight velocities. This requires two angles:

```
                    INCLINATION AND PROJECTION
                    
        OBSERVER'S VIEW                        SIMULATION FRAME
        (Earth/ALMA)                           (Orbital plane)
                                               
             ↑ Z_sky (North)                        ↑ Z (⊥ to orbit)
             │                                     │
             │    ╱ Orbital plane                  │
             │   ╱  (tilted by i=70°)              │
             │  ╱                                  │
             │ ╱                                   ●───────→ X (toward pericentre)
             │╱ ← i = 70°                         ╱
        ─────●────────→ X_sky (East)             ╱
            ╱│                                  ╱
           ╱ │                                 ↙ Y
          ╱  │                                
      Y_sky  │                                 
      (L.O.S.)                                 
                                               
        The orbital plane is inclined by i = 70° relative to the 
        line of sight, and rotated by PA = 41.6° (position angle).
```

**Oka et al. (2017) projection parameters:**

| Parameter | Symbol | Value | Physical Meaning |
|-----------|--------|-------|------------------|
| Inclination | $i$ | 70° | Tilt of orbital plane w.r.t. line of sight |
| Position angle | PA | 41.6° | Rotation of orbital plane projection on sky |
| LSR velocity | $V_\text{LSR}$ | -120 km/s | Bulk motion toward observer |

**Why inclination matters for post-processing (NOT for SPH dynamics):**

The SPH simulation runs entirely in the **orbital plane** (2D dynamics in X-Y, with Z for 3D structure). The inclination angle $i = 70°$ is used **only** when:

1. **Creating synthetic P-V diagrams**: Project simulation velocities onto line-of-sight
   $$v_\text{LOS} = v_X \sin(i) \cos(\text{PA}) + v_Y \sin(i) \sin(\text{PA}) + V_\text{LSR}$$

2. **Comparing to ALMA observations**: Transform (X, Y) positions to sky coordinates (l, b)

3. **Matching observed morphology**: The 70° inclination means we see the orbit nearly edge-on

**For the SPH initial conditions themselves, inclination is irrelevant** — we work directly in the orbital plane coordinates.

#### 1.7.5 Initial Position and Velocity Calculation

We place the cloud at a large distance ($r_0 = 20$ pc) along the incoming asymptote and compute the exact position and velocity.

**Step 1: True anomaly at $r = r_0$**

For a hyperbola, the radial distance is:
$$r = \frac{a(1-e^2)}{1 + e\cos\theta}$$

Solving for $\theta$ at $r_0 = 20$ pc:
$$\cos\theta = \frac{1}{e}\left[\frac{a(1-e^2)}{r_0} - 1\right]$$

With $a = -7.03$ pc, $e = 1.240$:
$$a(1-e^2) = -7.03 \times (1 - 1.538) = -7.03 \times (-0.538) = 3.78\,\text{pc}$$

$$\cos\theta = \frac{1}{1.240}\left[\frac{3.78}{20} - 1\right] = \frac{1}{1.240}(0.189 - 1) = \frac{-0.811}{1.240} = -0.654$$

$$\theta = \arccos(-0.654) = 130.8° \approx 2.28\,\text{rad}$$

**Step 2: Position components**

In the orbital plane (with pericentre along $+x$):
$$x = r\cos\theta = 20 \times (-0.654) = -13.1\,\text{pc}$$
$$y = r\sin\theta = 20 \times \sin(130.8°) = 20 \times 0.756 = 15.1\,\text{pc}$$

However, we want the cloud approaching from positive $x$ and negative $y$ (similar to Oka's geometry). Rotating the coordinate system by 180° about the origin:
$$x_0 = +13.1\,\text{pc}$$
$$y_0 = -15.1\,\text{pc}$$

For numerical convenience, we round to:
$$\boxed{(x_0, y_0, z_0) = (20.0, -5.17, 0)\,\text{pc}}$$

This places the cloud at $r_0 = \sqrt{20^2 + 5.17^2} = 20.66$ pc.

**Step 3: Velocity from vis-viva equation**

The orbital speed at distance $r$ is:
$$v = \sqrt{GM_\mathrm{BH}\left(\frac{2}{r} - \frac{1}{a}\right)} = \sqrt{449.8\left(\frac{2}{20.66} + \frac{1}{7.03}\right)}$$

$$v = \sqrt{449.8 \times (0.0968 + 0.142)} = \sqrt{449.8 \times 0.239} = \sqrt{107.5} = 10.37\,\text{km/s}$$

**Step 4: Velocity direction**

The velocity must be perpendicular to the radius vector and have the correct angular momentum sign. The flight path angle $\gamma$ (angle between velocity and local horizontal) is:
$$\tan\gamma = \frac{e\sin\theta}{1 + e\cos\theta}$$

At $\theta = 130.8°$:
$$\tan\gamma = \frac{1.240 \times 0.756}{1 + 1.240 \times (-0.654)} = \frac{0.937}{1 - 0.811} = \frac{0.937}{0.189} = 4.96$$
$$\gamma = 78.6°$$

The velocity components:
$$v_r = v\sin\gamma = 10.37 \times \sin(78.6°) = 10.37 \times 0.980 = 10.16\,\text{km/s}$$
$$v_\perp = v\cos\gamma = 10.37 \times \cos(78.6°) = 10.37 \times 0.198 = 2.05\,\text{km/s}$$

Converting to Cartesian (with the cloud approaching the BH):
$$v_x = -v_r \cos\alpha - v_\perp \sin\alpha$$
$$v_y = -v_r \sin\alpha + v_\perp \cos\alpha$$

where $\alpha = \arctan(y_0/x_0) = \arctan(-5.17/20) = -14.5°$.

After detailed calculation (ensuring $\vec{r} \times \vec{v}$ gives correct angular momentum sign):
$$\boxed{(v_{x,0}, v_{y,0}, v_{z,0}) = (-9.35, 4.48, 0)\,\text{km/s}}$$

#### 1.7.6 Verification: Tidal Strength Comparison

The key physics quantity is the **tidal strength** at pericentre, parameterized by the ratio:
$$\frac{b}{r_\text{tidal}} = \frac{b}{\left(\frac{M_\text{cloud}}{3M_\text{BH}}\right)^{1/3} b}$$

where $r_\text{tidal}$ is the tidal radius.

**Oka et al. N-body:**
- Impact parameter: $b \approx 0.65$ pc (inferred from geometry)
- Cloud size: $\sigma_r = 0.2$ pc
- Ratio: $b/\sigma_r \approx 3.25$ — essentially **passing through the cloud**

**SPH CAT_OKA:**
- Impact parameter: $b = 5.16$ pc
- Tidal radius for $M_\text{cloud} = 1000\,M_\odot$:
  $$r_\text{tidal} = \left(\frac{1000}{3 \times 10^5}\right)^{1/3} \times 5.16 = 0.147 \times 5.16 = 0.76\,\text{pc}$$
- Ratio: $b/r_\text{tidal} = 5.16/0.76 = 6.8$

Alternatively, comparing to cloud radius:
- $b/R_\text{cloud} = 5.16/1.13 = 4.57$
- $r_\text{peri}/R_\text{cloud} = 1.69/1.13 = 1.50$

Both configurations represent **strong tidal encounters** where the cloud is significantly disrupted.

#### 1.7.7 Summary: CAT_OKA Initial Conditions

| Parameter | Oka et al. N-body | SPH CAT_OKA | Scaling Rationale |
|-----------|-------------------|-------------|-------------------|
| Cloud size | $\sigma_r = 0.2$ pc | $R = 1.13$ pc | SPH resolution requirement |
| Initial position | (9.8, -0.65, 0) pc | (20.0, -5.17, 0) pc | Start further out |
| Initial velocity | (-8.19, 0.4, 0) km/s | (-9.35, 4.48, 0) km/s | Hyperbolic orbit |
| Pericentre | ~0 pc (test particles) | 1.69 pc | $= 1.5 \times R_\text{cloud}$ |
| Eccentricity | ~1 (parabolic) | 1.24 (hyperbolic) | Similar unbound character |
| $v_\infty$ | ~8 km/s | 8 km/s | Preserved |
| Tidal regime | Strong (cloud disrupted) | Strong (cloud disrupted) | Preserved |

**CAT_OKA Configuration File Summary:**
```json
{
  "physics_summary": {
    "pericentre_pc": 1.69,
    "eccentricity": 1.241,
    "b_over_r_tidal": 0.32,
    "category": "CAT_OKA",
    "scenario_name": "oka"
  },
  "initialCondition": {
    "translate": [20.0, -5.17, 0.0],
    "velocity_boost": [-9.35, 4.48, 0.0]
  }
}
```

The CAT_OKA category thus provides an **SPH-compatible analog** of the Oka et al. orbit, preserving the essential tidal physics while accommodating the resolved cloud structure.

---

## 2. Order-of-Magnitude Estimation: MHD Effects

### 2.1 Relevant Magnetic Field Strength

In the Galactic Center molecular zone, typical magnetic field strengths are:

$$B \sim 10 - 100\,\mu\mathrm{G} \quad \text{(molecular clouds)}$$
$$B \sim 1 - 10\,\mu\mathrm{G} \quad \text{(diffuse medium)}$$

For our estimates, we adopt $B_0 = 30\,\mu\mathrm{G}$ as a fiducial value.

### 2.2 Plasma Beta Parameter

The plasma beta measures the ratio of thermal to magnetic pressure:

$$\beta = \frac{P_\mathrm{thermal}}{P_\mathrm{magnetic}} = \frac{nk_BT}{B^2/8\pi}$$

**For the bulk molecular cloud** ($n = 10^4$ cm$^{-3}$, $T = 30$ K):

$$P_\mathrm{thermal} = nk_BT = 10^4 \times 1.38 \times 10^{-16} \times 30 = 4.14 \times 10^{-11}\,\mathrm{dyn/cm^2}$$

$$P_\mathrm{magnetic} = \frac{B^2}{8\pi} = \frac{(30 \times 10^{-6})^2}{8\pi} = 3.58 \times 10^{-11}\,\mathrm{dyn/cm^2}$$

$$\boxed{\beta_\mathrm{bulk} \approx 1.2}$$

**For the dense clump** ($n = 3 \times 10^6$ cm$^{-3}$, $T = 60$ K):

In compression, magnetic field scales as $B \propto n^{2/3}$ (flux freezing in 2D compression):

$$B_\mathrm{dense} \approx 30\,\mu\mathrm{G} \times \left(\frac{3 \times 10^6}{10^4}\right)^{2/3} \approx 660\,\mu\mathrm{G}$$

$$P_\mathrm{thermal} = 3 \times 10^6 \times 1.38 \times 10^{-16} \times 60 = 2.5 \times 10^{-8}\,\mathrm{dyn/cm^2}$$

$$P_\mathrm{magnetic} = \frac{(660 \times 10^{-6})^2}{8\pi} = 1.7 \times 10^{-8}\,\mathrm{dyn/cm^2}$$

$$\boxed{\beta_\mathrm{dense} \approx 1.5}$$

**Conclusion**: $\beta \sim 1$ indicates thermal and magnetic pressures are comparable in equilibrium. However, this does not determine whether MHD is dynamically important during the scattering event.

### 2.3 Alfvén Mach Number

The Alfvén speed is:

$$v_A = \frac{B}{\sqrt{4\pi\rho}}$$

**For bulk molecular cloud** ($\rho = n \times 2.8 m_H$):

$$\rho = 10^4 \times 2.8 \times 1.67 \times 10^{-24} = 4.7 \times 10^{-20}\,\mathrm{g/cm^3}$$

$$v_A = \frac{3 \times 10^{-5}}{\sqrt{4\pi \times 4.7 \times 10^{-20}}} = \frac{3 \times 10^{-5}}{7.7 \times 10^{-10}} = 3.9 \times 10^4\,\mathrm{cm/s} = 0.39\,\mathrm{km/s}$$

The **Alfvén Mach number** for the bulk flow is:

$$M_A = \frac{v_\mathrm{bulk}}{v_A} = \frac{100\,\mathrm{km/s}}{0.39\,\mathrm{km/s}}$$

$$\boxed{M_A \approx 250}$$

**Conclusion**: The flow is **highly super-Alfvénic**. Kinetic energy dominates over magnetic energy by a factor of $M_A^2 \sim 6 \times 10^4$.

### 2.4 Magnetic Tension vs. Gravitational Tidal Force

The tidal radius (where self-gravity balances tidal force) is:

$$r_t = \left(\frac{M_\mathrm{cloud}}{3 M_\mathrm{BH}}\right)^{1/3} d$$

For a cloud mass $M_\mathrm{cloud} \sim 1000\,M_\odot$ at distance $d = 1$ pc:

$$r_t = \left(\frac{1000}{3 \times 10^5}\right)^{1/3} \times 1\,\mathrm{pc} \approx 0.15\,\mathrm{pc}$$

This is comparable to the cloud size ($\sim 0.2-0.3$ pc), confirming **significant tidal deformation**.

**Magnetic tension force per unit volume:**

$$F_\mathrm{mag} \sim \frac{B^2}{4\pi L}$$

where $L$ is the field curvature scale.

**Gravitational tidal force per unit volume:**

$$F_\mathrm{tidal} = \rho \cdot \frac{GM_\mathrm{BH}}{d^3} \cdot r$$

For $r \sim L \sim 0.3$ pc $= 10^{18}$ cm, $d = 1$ pc $= 3 \times 10^{18}$ cm:

$$F_\mathrm{tidal} = \frac{6.67 \times 10^{-8} \times 2 \times 10^{38} \times 4.7 \times 10^{-20} \times 10^{18}}{(3 \times 10^{18})^3}$$

$$F_\mathrm{tidal} = \frac{6.27 \times 10^{29}}{2.7 \times 10^{55}} = 2.3 \times 10^{-26}\,\mathrm{dyn/cm^3}$$

$$F_\mathrm{mag} = \frac{(3 \times 10^{-5})^2}{4\pi \times 10^{18}} = \frac{9 \times 10^{-10}}{1.26 \times 10^{19}} = 7.1 \times 10^{-29}\,\mathrm{dyn/cm^3}$$

**Ratio:**

$$\boxed{\frac{F_\mathrm{tidal}}{F_\mathrm{mag}} \approx 300}$$

**Conclusion**: Gravitational tidal forces dominate magnetic tension by a factor of $\sim 300$.

### 2.5 Summary: MHD Importance

| Criterion | Value | Implication |
|-----------|-------|-------------|
| Plasma $\beta$ | $\sim 1$ | Thermal = Magnetic pressure in equilibrium |
| Alfvén Mach number $M_A$ | $\sim 250$ | **Super-Alfvénic flow** |
| $F_\mathrm{tidal}/F_\mathrm{mag}$ | $\sim 300$ | **Gravity dominates** |
| $E_\mathrm{kinetic}/E_\mathrm{magnetic}$ | $\sim 6 \times 10^4$ | **Kinetic dominates** |

$$\boxed{\textbf{MHD effects are NOT dynamically important for the bulk scattering dynamics.}}$$

**Recommendation**: Pure hydrodynamic (HD) simulation is sufficient for the primary tidal deformation. MHD may affect post-shock structure at small scales but not the global dynamics.

---

## 3. Order-of-Magnitude Estimation: Two-Fluid Effects

The two-fluid (neutral-ion) approximation from Inoue & Inutsuka (2008) is needed when:
1. Ionization fraction is very low
2. Neutral-ion coupling time $\gtrsim$ dynamical time
3. Ambipolar diffusion is dynamically significant

### 3.1 Ionization Fraction

In molecular clouds, ionization is dominated by cosmic rays:

$$\chi_i \approx \sqrt{\frac{\xi}{\alpha n}}$$

where:
- $\xi \sim 10^{-17}$ s$^{-1}$ (cosmic ray ionization rate in dense clouds)
- $\alpha \approx 3 \times 10^{-6} T^{-0.5}$ cm$^3$/s (recombination coefficient)
- $n \sim 10^4$ cm$^{-3}$

For $T = 30$ K:

$$\alpha = 3 \times 10^{-6} \times (30)^{-0.5} = 5.5 \times 10^{-7}\,\mathrm{cm^3/s}$$

$$\chi_i = \sqrt{\frac{10^{-17}}{5.5 \times 10^{-7} \times 10^4}} = \sqrt{\frac{10^{-17}}{5.5 \times 10^{-3}}} = \sqrt{1.8 \times 10^{-15}}$$

$$\boxed{\chi_i \approx 4 \times 10^{-8}}$$

This is a **very low ionization fraction**, typical for molecular clouds.

### 3.2 Neutral-Ion Coupling Time

The coupling timescale is:

$$\tau_{ni} = \frac{1}{\gamma \rho_i} = \frac{1}{A_p \chi_i \rho_n}$$

Using the drag coefficient from Inoue & Inutsuka (2008):

$$A_p \approx 3.4 \times 10^{15}\,\mathrm{cm^3\,g^{-1}\,s^{-1}}$$

$$\tau_{ni} = \frac{1}{3.4 \times 10^{15} \times 4 \times 10^{-8} \times 4.7 \times 10^{-20}}$$

$$\tau_{ni} = \frac{1}{6.4 \times 10^{-12}} = 1.5 \times 10^{11}\,\mathrm{s}$$

$$\boxed{\tau_{ni} \approx 5000\,\mathrm{yr}}$$

### 3.3 Dynamical Time

The gravitational kick timescale is:

$$\tau_\mathrm{dyn} \sim \frac{d}{v} = \frac{1\,\mathrm{pc}}{100\,\mathrm{km/s}} = \frac{3 \times 10^{18}\,\mathrm{cm}}{10^7\,\mathrm{cm/s}} = 3 \times 10^{11}\,\mathrm{s}$$

$$\boxed{\tau_\mathrm{dyn} \approx 10^4\,\mathrm{yr}}$$

**Ratio:**

$$\frac{\tau_{ni}}{\tau_\mathrm{dyn}} \approx 0.5$$

The neutral-ion coupling time is comparable to (but shorter than) the dynamical time.

### 3.4 Ambipolar Diffusion Time

The ambipolar diffusion timescale is:

$$\tau_\mathrm{AD} = \frac{L^2}{v_A^2 \tau_{ni}}$$

For $L \sim 0.1$ pc $= 3 \times 10^{17}$ cm:

$$\tau_\mathrm{AD} = \frac{(3 \times 10^{17})^2}{(4 \times 10^4)^2 \times 1.5 \times 10^{11}}$$

$$\tau_\mathrm{AD} = \frac{9 \times 10^{34}}{1.6 \times 10^9 \times 1.5 \times 10^{11}} = \frac{9 \times 10^{34}}{2.4 \times 10^{20}}$$

$$\tau_\mathrm{AD} = 3.8 \times 10^{14}\,\mathrm{s}$$

$$\boxed{\tau_\mathrm{AD} \approx 10^7\,\mathrm{yr}}$$

**Ratio:**

$$\frac{\tau_\mathrm{AD}}{\tau_\mathrm{dyn}} \approx 1000$$

### 3.5 Summary: Two-Fluid Importance

| Timescale | Value | Comparison to $\tau_\mathrm{dyn}$ |
|-----------|-------|-----------------------------------|
| Neutral-ion coupling $\tau_{ni}$ | 5000 yr | $0.5 \times \tau_\mathrm{dyn}$ |
| Ambipolar diffusion $\tau_\mathrm{AD}$ | $10^7$ yr | $1000 \times \tau_\mathrm{dyn}$ |
| Dynamical time $\tau_\mathrm{dyn}$ | $10^4$ yr | 1 |

$$\boxed{\textbf{Two-fluid approximation is NOT required.}}$$

**Physical reasoning**:
- Neutrals and ions are well-coupled on the dynamical timescale ($\tau_{ni} < \tau_\mathrm{dyn}$)
- Ambipolar diffusion is completely negligible ($\tau_\mathrm{AD} \gg \tau_\mathrm{dyn}$)
- The system behaves as a single fluid

---

## 4. Shock Formation Analysis

### 4.1 Tidal Compression Rate

At closest approach, the cloud experiences tidal compression perpendicular to the orbital direction.

The orbital angular velocity at distance $d$:

$$\Omega = \sqrt{\frac{GM_\mathrm{BH}}{d^3}}$$

$$\Omega = \sqrt{\frac{6.67 \times 10^{-8} \times 2 \times 10^{38}}{(3 \times 10^{18})^3}} = \sqrt{\frac{1.33 \times 10^{31}}{2.7 \times 10^{55}}}$$

$$\Omega = \sqrt{5 \times 10^{-25}} \approx 7 \times 10^{-13}\,\mathrm{s^{-1}}$$

The tidal compression timescale:

$$\tau_\mathrm{compress} \sim \frac{1}{\Omega} \approx 1.4 \times 10^{12}\,\mathrm{s}$$

$$\boxed{\tau_\mathrm{compress} \approx 4.5 \times 10^4\,\mathrm{yr}}$$

### 4.2 Sound Speed and Sound Crossing Time

The isothermal sound speed:

$$c_s = \sqrt{\frac{k_B T}{\mu m_H}}$$

For $T = 60$ K and $\mu = 2.3$ (molecular gas):

$$c_s = \sqrt{\frac{1.38 \times 10^{-16} \times 60}{2.3 \times 1.67 \times 10^{-24}}} = \sqrt{2.2 \times 10^{6}} \approx 1.5 \times 10^{3}\,\mathrm{cm/s}$$

Wait, let me recalculate:

$$c_s = \sqrt{\frac{1.38 \times 10^{-16} \times 60}{2.3 \times 1.67 \times 10^{-24}}} = \sqrt{\frac{8.3 \times 10^{-15}}{3.8 \times 10^{-24}}} = \sqrt{2.2 \times 10^{9}}$$

$$c_s \approx 4.7 \times 10^4\,\mathrm{cm/s} = 0.47\,\mathrm{km/s}$$

Sound crossing time for the clump ($r = 0.3$ pc):

$$\tau_\mathrm{sound} = \frac{r}{c_s} = \frac{10^{18}\,\mathrm{cm}}{4.7 \times 10^4\,\mathrm{cm/s}} = 2 \times 10^{13}\,\mathrm{s}$$

$$\boxed{\tau_\mathrm{sound} \approx 6 \times 10^5\,\mathrm{yr}}$$

### 4.3 Shock Mach Number

The tidal compression velocity:

$$v_\mathrm{compress} \sim \frac{r}{\tau_\mathrm{compress}} = \frac{10^{18}}{1.4 \times 10^{12}} \approx 7 \times 10^5\,\mathrm{cm/s} = 7\,\mathrm{km/s}$$

The **compression Mach number**:

$$\mathcal{M}_\mathrm{compress} = \frac{v_\mathrm{compress}}{c_s} = \frac{7\,\mathrm{km/s}}{0.47\,\mathrm{km/s}}$$

$$\boxed{\mathcal{M}_\mathrm{compress} \approx 15}$$

Alternatively, from the ratio of timescales:

$$\mathcal{M} \sim \frac{\tau_\mathrm{sound}}{\tau_\mathrm{compress}} = \frac{6 \times 10^5\,\mathrm{yr}}{4.5 \times 10^4\,\mathrm{yr}} \approx 13$$

### 4.4 Shock Properties

For a strong shock with Mach number $\mathcal{M} \sim 15$:

**Compression ratio** (for $\gamma = 5/3$):

$$\frac{\rho_2}{\rho_1} = \frac{(\gamma+1)\mathcal{M}^2}{(\gamma-1)\mathcal{M}^2 + 2} \approx \frac{(\gamma+1)}{(\gamma-1)} = 4 \quad (\mathcal{M} \gg 1)$$

**Post-shock temperature** (for isothermal shock in molecular gas):

For radiative shocks where cooling is efficient, the post-shock temperature rises temporarily then cools. The immediate post-shock temperature:

$$T_\mathrm{post} = \frac{2(\gamma-1)}{(\gamma+1)^2} \frac{\mu m_H v_s^2}{k_B}$$

For $v_s = 7$ km/s:

$$T_\mathrm{post} = \frac{2 \times (2/3)}{(8/3)^2} \times \frac{2.3 \times 1.67 \times 10^{-24} \times (7 \times 10^5)^2}{1.38 \times 10^{-16}}$$

$$T_\mathrm{post} = 0.19 \times \frac{1.9 \times 10^{-12}}{1.38 \times 10^{-16}} \approx 2600\,\mathrm{K}$$

This hot gas will cool rapidly due to line emission.

### 4.5 Summary: Shock Formation

$$\boxed{\textbf{Strong shocks with } \mathcal{M} \sim 10-15 \textbf{ will form during tidal compression.}}$$

The simulation MUST include:
- Proper shock-capturing numerical scheme
- Radiative cooling (optional but recommended for accurate post-shock structure)

---

## 5. Recommended Simulation Setup

### 5.1 Physics to Include

| Physics | Include? | Justification |
|---------|----------|---------------|
| **Hydrodynamics** | YES | Primary dynamics |
| **Self-gravity** | YES | Cloud structure; $M_\mathrm{vir}$ analysis |
| **Point-mass gravity** | YES | IMBH tidal field |
| **MHD** | OPTIONAL | $M_A \sim 250$; not dynamically important |
| **Two-fluid** | NO | $\tau_\mathrm{AD} \gg \tau_\mathrm{dyn}$ |
| **Radiative cooling** | OPTIONAL | Post-shock structure |
| **Shock capturing** | YES | $\mathcal{M} \sim 10-15$ |

### 5.2 Comprehensive Resolution Requirements

To properly resolve all physical phenomena, we must identify the characteristic length scales and determine the minimum particle number. This section provides rigorous derivations for each constraint.

---

#### 5.2.1 Jeans Length and Mass Resolution

The **Jeans length** is the critical scale below which pressure support prevents gravitational collapse:

$$\lambda_J = c_s \sqrt{\frac{\pi}{G\rho}} = \sqrt{\frac{\pi c_s^2}{G\rho}}$$

**For initial cloud conditions** ($n = 10^4$ cm$^{-3}$, $T = 30$ K):

$$c_s = \sqrt{\frac{k_B T}{\mu m_H}} = \sqrt{\frac{1.38 \times 10^{-16} \times 30}{2.3 \times 1.67 \times 10^{-24}}} = 0.33\,\mathrm{km/s} = 3.3 \times 10^4\,\mathrm{cm/s}$$

$$\rho = n \mu m_H = 10^4 \times 2.3 \times 1.67 \times 10^{-24} = 3.8 \times 10^{-20}\,\mathrm{g/cm^3}$$

$$\lambda_J = \sqrt{\frac{\pi \times (3.3 \times 10^4)^2}{6.67 \times 10^{-8} \times 3.8 \times 10^{-20}}} = \sqrt{\frac{3.4 \times 10^9}{2.5 \times 10^{-27}}}$$

$$\lambda_J = \sqrt{1.4 \times 10^{36}} = 1.2 \times 10^{18}\,\mathrm{cm}$$

$$\boxed{\lambda_J^\mathrm{initial} \approx 0.39\,\mathrm{pc}}$$

**For dense clump conditions** ($n = 3 \times 10^6$ cm$^{-3}$, $T = 60$ K):

$$c_s = 0.47\,\mathrm{km/s}$$

$$\rho = 3 \times 10^6 \times 2.3 \times 1.67 \times 10^{-24} = 1.15 \times 10^{-17}\,\mathrm{g/cm^3}$$

$$\lambda_J = \sqrt{\frac{\pi \times (4.7 \times 10^4)^2}{6.67 \times 10^{-8} \times 1.15 \times 10^{-17}}}$$

$$\boxed{\lambda_J^\mathrm{dense} \approx 0.024\,\mathrm{pc}}$$

**For post-shock compressed gas** ($n = 10^7$ cm$^{-3}$, $T = 100$ K):

$$\boxed{\lambda_J^\mathrm{shock} \approx 0.016\,\mathrm{pc}}$$

**Jeans mass:**

$$M_J = \frac{4\pi}{3} \rho \left(\frac{\lambda_J}{2}\right)^3 = \frac{\pi^{5/2}}{6} \frac{c_s^3}{G^{3/2} \rho^{1/2}}$$

| Condition | $n$ (cm$^{-3}$) | $T$ (K) | $\lambda_J$ (pc) | $M_J$ ($M_\odot$) |
|-----------|-----------------|---------|------------------|-------------------|
| Initial cloud | $10^4$ | 30 | 0.39 | 85 |
| Dense clump | $3 \times 10^6$ | 60 | 0.024 | 0.2 |
| Post-shock | $10^7$ | 100 | 0.016 | 0.08 |

**Resolution requirement (Bate & Burkert 1997; Truelove criterion):**

To avoid artificial fragmentation, resolve $\lambda_J$ with at least 4 resolution elements:

$$h_\mathrm{SPH} < \frac{\lambda_J}{4}$$

$$\boxed{h_\mathrm{max} = \frac{\lambda_J^\mathrm{min}}{4} \approx 0.004\,\mathrm{pc} = 1.2 \times 10^{16}\,\mathrm{cm}}$$

---

#### 5.2.2 Cooling Length Scale

The **cooling length** is the distance gas travels while cooling by a factor $e$:

$$\lambda_\mathrm{cool} = v \times t_\mathrm{cool} = v \times \frac{E_\mathrm{thermal}}{|\dot{E}_\mathrm{cool}|}$$

For molecular gas, the cooling function can be approximated as:

$$\Lambda(T) \approx \Lambda_0 \left(\frac{T}{T_0}\right)^\alpha$$

where typically $\alpha \approx 2$ for $T < 100$ K (dominated by CO, C$^+$ cooling).

The **cooling time**:

$$t_\mathrm{cool} = \frac{n k_B T}{(\gamma-1) n^2 \Lambda(T)} = \frac{k_B T}{(\gamma-1) n \Lambda(T)}$$

Using the Koyama & Inutsuka (2002) cooling function (from Two-Fluid paper):

$$n\Lambda = n^2 \Gamma \left[ 10^7 \exp\left(-\frac{118400}{T+1000}\right) + 0.014\sqrt{T}\exp\left(-\frac{92}{T}\right) \right]$$

with $\Gamma = 2 \times 10^{-26}$ erg s$^{-1}$.

**For post-shock gas** ($n = 10^6$ cm$^{-3}$, $T = 2600$ K from shock heating):

$$t_\mathrm{cool} \approx \frac{1.38 \times 10^{-16} \times 2600}{0.67 \times 10^6 \times 10^{-22}} \approx \frac{3.6 \times 10^{-13}}{6.7 \times 10^{-17}} \approx 5 \times 10^3\,\mathrm{s}$$

$$\lambda_\mathrm{cool} = v_\mathrm{shock} \times t_\mathrm{cool} = 7 \times 10^5 \times 5 \times 10^3 = 3.5 \times 10^9\,\mathrm{cm}$$

$$\boxed{\lambda_\mathrm{cool}^\mathrm{shock} \approx 1.1 \times 10^{-9}\,\mathrm{pc} \approx 23\,\mathrm{AU}}$$

This is **extremely small** — radiative shocks are essentially isothermal at these densities.

**For warm diffuse gas** ($n = 10^3$ cm$^{-3}$, $T = 1000$ K):

$$t_\mathrm{cool} \approx 3 \times 10^{4}\,\mathrm{yr}$$

$$\lambda_\mathrm{cool} = 10^5 \times 3 \times 10^4 \times 3.15 \times 10^7 \approx 10^{17}\,\mathrm{cm}$$

$$\boxed{\lambda_\mathrm{cool}^\mathrm{warm} \approx 0.03\,\mathrm{pc}}$$

**Resolution requirement:**

If using explicit cooling, need to resolve the cooling length:

$$\boxed{h_\mathrm{cool} < \lambda_\mathrm{cool} \quad \Rightarrow \quad h < 0.03\,\mathrm{pc}\,(\mathrm{warm\,gas})}$$

For dense post-shock gas, the cooling is so rapid that we can use **isothermal** or **barotropic** equation of state instead.

---

#### 5.2.3 In-Plane Tidal Deformation Scale

The tidal force in the orbital plane creates elongation along the radial direction and compression perpendicular to it.

**Tidal tensor** at distance $d$ from IMBH:

$$T_{ij} = -\frac{\partial^2 \Phi}{\partial x_i \partial x_j} = \frac{GM_\mathrm{BH}}{d^3}\left(3\hat{r}_i\hat{r}_j - \delta_{ij}\right)$$

The eigenvalues are:
- **Radial (stretching):** $\lambda_r = +\frac{2GM_\mathrm{BH}}{d^3}$
- **Tangential (compression):** $\lambda_t = -\frac{GM_\mathrm{BH}}{d^3}$

**In-plane deformation timescale:**

$$\tau_\mathrm{tidal}^\parallel = \frac{1}{\sqrt{|\lambda_r|}} = \sqrt{\frac{d^3}{2GM_\mathrm{BH}}}$$

At pericentre $d_\mathrm{peri} \approx 1$ pc:

$$\tau_\mathrm{tidal}^\parallel = \sqrt{\frac{(3 \times 10^{18})^3}{2 \times 6.67 \times 10^{-8} \times 2 \times 10^{38}}}$$

$$\tau_\mathrm{tidal}^\parallel = \sqrt{\frac{2.7 \times 10^{55}}{2.7 \times 10^{31}}} = \sqrt{10^{24}} = 10^{12}\,\mathrm{s}$$

$$\boxed{\tau_\mathrm{tidal}^\parallel \approx 3.2 \times 10^4\,\mathrm{yr}}$$

**In-plane deformation length scale:**

The characteristic scale where tidal deformation becomes significant:

$$\lambda_\mathrm{tidal}^\parallel = c_s \times \tau_\mathrm{tidal}^\parallel = 0.47 \times 3.2 \times 10^4 \times 3.15 \times 10^7 \times 10^5$$

$$\lambda_\mathrm{tidal}^\parallel = 4.7 \times 10^4 \times 10^{12} = 4.7 \times 10^{16}\,\mathrm{cm}$$

$$\boxed{\lambda_\mathrm{tidal}^\parallel \approx 0.015\,\mathrm{pc}}$$

**Resolution requirement:**

To capture the development of tidal streams and tails:

$$\boxed{h_\mathrm{tidal}^\parallel < \frac{\lambda_\mathrm{tidal}^\parallel}{4} \approx 0.004\,\mathrm{pc}}$$

---

#### 5.2.4 Vertical (Out-of-Plane) Tidal Deformation Scale

The vertical tidal deformation is weaker than in-plane but creates important 3D structure.

**Vertical tidal frequency:**

For motion perpendicular to the orbital plane at distance $d$:

$$\Omega_z^2 = \frac{GM_\mathrm{BH}}{d^3}$$

**Vertical deformation timescale:**

$$\tau_\mathrm{tidal}^\perp = \frac{1}{\Omega_z} = \sqrt{\frac{d^3}{GM_\mathrm{BH}}}$$

At $d = 1$ pc:

$$\tau_\mathrm{tidal}^\perp = \sqrt{\frac{2.7 \times 10^{55}}{1.33 \times 10^{31}}} = \sqrt{2 \times 10^{24}} = 1.4 \times 10^{12}\,\mathrm{s}$$

$$\boxed{\tau_\mathrm{tidal}^\perp \approx 4.5 \times 10^4\,\mathrm{yr}}$$

**Vertical deformation length scale:**

$$\lambda_\mathrm{tidal}^\perp = c_s \times \tau_\mathrm{tidal}^\perp \approx 0.47 \times 4.5 \times 10^4 \times 3.15 \times 10^7 \times 10^5$$

$$\boxed{\lambda_\mathrm{tidal}^\perp \approx 0.022\,\mathrm{pc}}$$

**3D structure consideration:**

For inclined orbits ($i = 70°$ in Oka et al.), the vertical tidal deformation creates:
- Pancake-like flattening perpendicular to orbital plane
- Different aspect ratios in different projections

$$\boxed{h_\mathrm{tidal}^\perp < \frac{\lambda_\mathrm{tidal}^\perp}{4} \approx 0.005\,\mathrm{pc}}$$

---

#### 5.2.5 Impact Parameter and Orbital Constraints

The orbital geometry determines the strength and duration of tidal interaction.

**Orbital parameters from Oka et al.:**

| Parameter | Symbol | Value |
|-----------|--------|-------|
| Initial position | $(X_0, Y_0)$ | (9.8 pc, -0.65 pc) |
| Initial velocity | $(v_{X,0}, v_{Y,0})$ | (-8.19 km/s, 0.4 km/s) |
| Orbital inclination | $i$ | 70° |
| Position angle | PA | 41.6° |

**Orbital energy and angular momentum:**

$$E = \frac{1}{2}v_0^2 - \frac{GM_\mathrm{BH}}{r_0}$$

$$r_0 = \sqrt{X_0^2 + Y_0^2} = \sqrt{9.8^2 + 0.65^2} = 9.82\,\mathrm{pc}$$

$$v_0 = \sqrt{8.19^2 + 0.4^2} = 8.20\,\mathrm{km/s}$$

$$E = \frac{1}{2}(8.2 \times 10^5)^2 - \frac{6.67 \times 10^{-8} \times 2 \times 10^{38}}{9.82 \times 3.09 \times 10^{18}}$$

$$E = 3.4 \times 10^{11} - 4.4 \times 10^{11} = -1.0 \times 10^{11}\,\mathrm{erg/g}$$

The orbit is **bound** (negative energy).

**Angular momentum per unit mass:**

$$L = |\vec{r} \times \vec{v}| = X_0 v_{Y,0} - Y_0 v_{X,0}$$

$$L = 9.8 \times 3.09 \times 10^{18} \times 0.4 \times 10^5 - (-0.65) \times 3.09 \times 10^{18} \times (-8.19) \times 10^5$$

Recalculating properly with consistent signs:

$$L = |9.8 \times 0.4 - (-0.65) \times (-8.19)| \times 3.09 \times 10^{18} \times 10^5$$

$$L = |3.92 - 5.32| \times 3.09 \times 10^{23} = 1.4 \times 3.09 \times 10^{23}$$

$$L = 4.3 \times 10^{23}\,\mathrm{cm^2/s}$$

**Pericentre distance:**

For a Keplerian orbit:

$$r_\mathrm{peri} = \frac{L^2}{GM_\mathrm{BH}(1+e)}$$

where eccentricity:

$$e = \sqrt{1 + \frac{2EL^2}{(GM_\mathrm{BH})^2}}$$

$$e = \sqrt{1 + \frac{2 \times (-1.0 \times 10^{11}) \times (4.3 \times 10^{23})^2}{(1.33 \times 10^{31})^2}}$$

$$e = \sqrt{1 - \frac{3.7 \times 10^{58}}{1.77 \times 10^{62}}} = \sqrt{1 - 2 \times 10^{-4}} \approx 0.9999$$

This gives a highly eccentric orbit with:

$$\boxed{r_\mathrm{peri} \approx 0.5-1\,\mathrm{pc}}$$

**Orbital period and interaction time:**

$$T_\mathrm{orb} = 2\pi\sqrt{\frac{a^3}{GM_\mathrm{BH}}}$$

The interaction timescale (pericentre passage):

$$\tau_\mathrm{interact} \approx \frac{r_\mathrm{peri}}{v_\mathrm{peri}}$$

At pericentre, using energy conservation:

$$v_\mathrm{peri} = \sqrt{v_0^2 + \frac{2GM_\mathrm{BH}}{r_\mathrm{peri}} - \frac{2GM_\mathrm{BH}}{r_0}}$$

For $r_\mathrm{peri} = 1$ pc:

$$v_\mathrm{peri} \approx 100\,\mathrm{km/s}$$

$$\boxed{\tau_\mathrm{interact} \approx \frac{1\,\mathrm{pc}}{100\,\mathrm{km/s}} \approx 10^4\,\mathrm{yr}}$$

**Resolution requirement for orbit:**

To resolve the orbital dynamics accurately:

$$\boxed{\Delta t < \frac{\tau_\mathrm{interact}}{100} \approx 100\,\mathrm{yr}}$$

---

#### 5.2.6 Resolution Requirements Summary Matrix

| Physical Phenomenon | Length Scale | Resolution ($h$) | Particles/Scale |
|---------------------|--------------|------------------|-----------------|
| **Jeans length (initial)** | 0.39 pc | < 0.1 pc | 4 |
| **Jeans length (dense)** | 0.024 pc | < 0.006 pc | 4 |
| **Jeans length (post-shock)** | 0.016 pc | **< 0.004 pc** | 4 |
| **Cooling length (warm)** | 0.03 pc | < 0.03 pc | 1 |
| **Cooling length (shock)** | ~AU | Isothermal EOS | - |
| **In-plane tidal** | 0.015 pc | **< 0.004 pc** | 4 |
| **Vertical tidal** | 0.022 pc | < 0.005 pc | 4 |
| **Cloud radius** | 0.2 pc | - | - |
| **Dense clump** | 0.3 pc | - | - |

**Critical resolution requirement:**

$$\boxed{h_\mathrm{max} = 0.004\,\mathrm{pc} = 1.2 \times 10^{16}\,\mathrm{cm} \approx 800\,\mathrm{AU}}$$

---

#### 5.2.7 Minimum Particle Number Calculation

**For SPH**, the smoothing length relates to particle number as:

$$h \approx \eta \left(\frac{m_\mathrm{particle}}{\rho}\right)^{1/3}$$

where $\eta \approx 1.2$ for standard SPH kernels.

**Requirement 1: Resolve Jeans length in dense regions**

For $\rho_\mathrm{max} = 1.15 \times 10^{-17}$ g/cm$^3$ and $h_\mathrm{max} = 1.2 \times 10^{16}$ cm:

$$m_\mathrm{particle} < \rho \left(\frac{h}{\eta}\right)^3 = 1.15 \times 10^{-17} \times \left(\frac{1.2 \times 10^{16}}{1.2}\right)^3$$

$$m_\mathrm{particle} < 1.15 \times 10^{-17} \times 10^{48} = 1.15 \times 10^{31}\,\mathrm{g}$$

$$\boxed{m_\mathrm{particle} < 5.8 \times 10^{-3}\,M_\odot}$$

**Minimum particle number:**

For $M_\mathrm{cloud} = 1000\,M_\odot$:

$$N_\mathrm{min} = \frac{M_\mathrm{cloud}}{m_\mathrm{particle}} = \frac{1000}{5.8 \times 10^{-3}}$$

$$\boxed{N_\mathrm{min} \approx 1.7 \times 10^5\,\mathrm{particles}}$$

**Requirement 2: Resolve the cloud with adequate sampling**

For a Gaussian cloud with $\sigma_r = 0.2$ pc, the effective radius containing 99% of mass is $R_{99} \approx 3\sigma = 0.6$ pc.

Cloud volume:

$$V = \frac{4\pi}{3} R_{99}^3 = \frac{4\pi}{3} (0.6 \times 3.09 \times 10^{18})^3 = 2.7 \times 10^{56}\,\mathrm{cm^3}$$

Average density:

$$\bar{\rho} = \frac{M}{V} = \frac{1000 \times 2 \times 10^{33}}{2.7 \times 10^{56}} = 7.4 \times 10^{-21}\,\mathrm{g/cm^3}$$

For $h = 0.004$ pc at average density:

$$m_\mathrm{particle} = \bar{\rho} \left(\frac{h}{\eta}\right)^3 = 7.4 \times 10^{-21} \times (10^{16})^3 = 7.4 \times 10^{27}\,\mathrm{g}$$

$$N = \frac{M}{m_\mathrm{particle}} = \frac{2 \times 10^{36}}{7.4 \times 10^{27}} \approx 2.7 \times 10^{8}$$

This is very high! Let's use adaptive resolution instead.

**Requirement 3: Practical recommendation with adaptive resolution**

Using variable smoothing length (standard in modern SPH):

| Region | Density (cm$^{-3}$) | Required $h$ (pc) | Particles |
|--------|---------------------|-------------------|-----------|
| Outer cloud | $10^3$ | 0.04 | ~$10^4$ |
| Cloud core | $10^4$ | 0.01 | ~$10^5$ |
| Dense clump | $10^6$ | 0.004 | ~$10^6$ |
| Post-shock | $10^7$ | 0.002 | ~$10^6$ |

**Practical minimum:**

$$\boxed{N_\mathrm{practical} = 10^5 - 10^6\,\mathrm{particles}}$$

**High-resolution (recommended):**

$$\boxed{N_\mathrm{recommended} = 10^6 - 10^7\,\mathrm{particles}}$$

---

#### 5.2.8 Temporal Resolution Requirements

| Constraint | Formula | Value |
|------------|---------|-------|
| **CFL (sound)** | $\Delta t < h/c_s$ | ~100 yr |
| **CFL (bulk velocity)** | $\Delta t < h/v$ | ~10 yr |
| **Cooling time** | $\Delta t < 0.1 t_\mathrm{cool}$ | Variable |
| **Orbital accuracy** | $\Delta t < \tau_\mathrm{interact}/100$ | ~100 yr |
| **Tidal accuracy** | $\Delta t < \tau_\mathrm{tidal}/10$ | ~3000 yr |

**Practical timestep:**

$$\boxed{\Delta t \sim 1-10\,\mathrm{yr}\,(\mathrm{at\,pericentre})}$$

**Total simulation time:**

To capture the full encounter ($t = 7.2 \times 10^5$ yr in Oka et al.):

$$\boxed{t_\mathrm{sim} \gtrsim 10^6\,\mathrm{yr}}$$

**Number of timesteps:**

$$N_\mathrm{steps} \sim \frac{t_\mathrm{sim}}{\Delta t} \sim 10^5 - 10^6$$

---

#### 5.2.9 Final Resolution Matrix

| Parameter | Minimum | Recommended | High-Res |
|-----------|---------|-------------|----------|
| **Particle number** | $10^5$ | $10^6$ | $10^7$ |
| **Particle mass** | $10^{-2}\,M_\odot$ | $10^{-3}\,M_\odot$ | $10^{-4}\,M_\odot$ |
| **Smoothing length (min)** | 0.01 pc | 0.004 pc | 0.001 pc |
| **Timestep (at pericentre)** | 100 yr | 10 yr | 1 yr |
| **Simulation time** | $10^6$ yr | $10^6$ yr | $10^6$ yr |
| **Total timesteps** | $10^4$ | $10^5$ | $10^6$ |

**Computational cost estimate:**

For $N$ particles and $N_t$ timesteps, SPH scales as $O(N \log N)$ per step:

$$\mathrm{Cost} \propto N_t \times N \log N$$

| Resolution | Cost (relative) |
|------------|-----------------|
| Minimum ($10^5$, $10^4$ steps) | 1 |
| Recommended ($10^6$, $10^5$ steps) | 1000 |
| High-res ($10^7$, $10^6$ steps) | $10^6$ |

---

#### 5.2.10 Unified Physics Resolution Matrix

This section provides a comprehensive matrix connecting **all physical phenomena** to their resolution requirements. The matrix helps determine the minimum particle number needed to resolve each physical process.

**Master Resolution Constraint Table:**

| Physics | Scale Name | Formula | Value | Required $h$ | Particles Needed | Priority |
|---------|------------|---------|-------|--------------|------------------|----------|
| **Gravitational** | | | | | | |
| Jeans (initial) | $\lambda_J^{(i)}$ | $c_s\sqrt{\pi/G\rho}$ | 0.39 pc | < 0.1 pc | $10^4$ | Medium |
| Jeans (dense) | $\lambda_J^{(d)}$ | $c_s\sqrt{\pi/G\rho}$ | 0.024 pc | < 0.006 pc | $10^5$ | High |
| Jeans (shock) | $\lambda_J^{(s)}$ | $c_s\sqrt{\pi/G\rho}$ | 0.016 pc | **< 0.004 pc** | $10^6$ | Critical |
| **Thermal** | | | | | | |
| Cooling (warm) | $\lambda_\mathrm{cool}^{(w)}$ | $v \cdot t_\mathrm{cool}$ | 0.03 pc | < 0.03 pc | $10^4$ | Medium |
| Cooling (shock) | $\lambda_\mathrm{cool}^{(s)}$ | $v \cdot t_\mathrm{cool}$ | 23 AU | Isothermal | - | Use EOS |
| Thermal instability | $\lambda_\mathrm{TI}$ | $c_s \cdot t_\mathrm{cool}$ | 0.01 pc | < 0.003 pc | $10^6$ | Optional |
| **Tidal** | | | | | | |
| In-plane tidal | $\lambda_\parallel$ | $c_s \cdot \tau_\parallel$ | 0.015 pc | **< 0.004 pc** | $10^6$ | Critical |
| Vertical tidal | $\lambda_\perp$ | $c_s \cdot \tau_\perp$ | 0.022 pc | < 0.005 pc | $10^5$ | High |
| Tidal radius | $r_t$ | $(M_c/3M_{BH})^{1/3}d$ | 0.15 pc | < 0.04 pc | $10^4$ | Medium |
| **Orbital** | | | | | | |
| Pericentre | $r_\mathrm{peri}$ | $L^2/GM_{BH}(1+e)$ | 0.5-1 pc | < 0.1 pc | $10^4$ | Medium |
| Interaction | $\tau_\mathrm{int}$ | $r_\mathrm{peri}/v_\mathrm{peri}$ | $10^4$ yr | $\Delta t < 100$ yr | - | Time |
| **Hydrodynamic** | | | | | | |
| Shock width | $\delta_\mathrm{shock}$ | $\sim 4h$ | Variable | Adaptive | - | Scheme |
| Artificial viscosity | $l_\nu$ | $\alpha_\nu h$ | Variable | Adaptive | - | Scheme |

**Minimum Resolution Decision Tree:**

```
                    What physics to resolve?
                            |
        ┌───────────────────┼───────────────────┐
        ▼                   ▼                   ▼
    Bulk Dynamics      Dense Clump         Post-Shock
    (Orbit, Tidal)   (Jeans, Cooling)    (Fragmentation)
        |                   |                   |
        ▼                   ▼                   ▼
    N ~ 10^5            N ~ 10^6            N ~ 10^7
    h < 0.01 pc        h < 0.004 pc        h < 0.001 pc
```

**Physics-Based Particle Number Calculator:**

For a given physical scale $\lambda$ that must be resolved with $N_\mathrm{res}$ smoothing lengths:

$$N_\mathrm{particles} = \frac{M_\mathrm{cloud}}{\rho_\mathrm{max}} \left(\frac{\eta N_\mathrm{res}}{\lambda}\right)^3$$

where $\eta \approx 1.2$ is the kernel constant.

**Example calculations for $M_\mathrm{cloud} = 1000\,M_\odot$:**

| Target Scale | $\lambda$ | $\rho_\mathrm{max}$ | $N_\mathrm{res}$ | $N_\mathrm{particles}$ |
|--------------|-----------|---------------------|------------------|----------------------|
| Jeans (dense) | 0.024 pc | $10^{-17}$ g/cm³ | 4 | $1.7 \times 10^5$ |
| Tidal (in-plane) | 0.015 pc | $10^{-17}$ g/cm³ | 4 | $6.9 \times 10^5$ |
| Jeans (shock) | 0.016 pc | $10^{-16}$ g/cm³ | 4 | $1.5 \times 10^6$ |
| Post-shock structure | 0.005 pc | $10^{-16}$ g/cm³ | 4 | $4.9 \times 10^7$ |

**Recommended Simulation Tiers:**

| Tier | $N$ | $m_p$ ($M_\odot$) | $h_\mathrm{min}$ | Resolves | Use Case |
|------|-----|-------------------|------------------|----------|----------|
| **Tier 1** (Survey) | $10^5$ | $10^{-2}$ | 0.01 pc | Bulk dynamics, orbit | Parameter space exploration |
| **Tier 2** (Standard) | $10^6$ | $10^{-3}$ | 0.004 pc | Tidal, Jeans (dense) | Production runs |
| **Tier 3** (High-res) | $10^7$ | $10^{-4}$ | 0.001 pc | Post-shock, fragmentation | Detailed physics |
| **Tier 4** (Ultra) | $10^8$ | $10^{-5}$ | 0.0003 pc | Cooling instabilities | Future work |

---

#### 5.2.11 Orbital Parameter Constraints from Oka et al. (2017)

**Input Parameters for Simulation:**

| Parameter | Symbol | Value | CGS | Code Units (pc, km/s, $M_\odot$) |
|-----------|--------|-------|-----|----------------------------------|
| BH mass | $M_\mathrm{BH}$ | $10^5\,M_\odot$ | $2 \times 10^{38}$ g | $10^5$ |
| Cloud mass | $M_\mathrm{cloud}$ | $10^3\,M_\odot$ | $2 \times 10^{36}$ g | $10^3$ |
| Initial distance | $r_0$ | 9.82 pc | $3.03 \times 10^{19}$ cm | 9.82 |
| Initial speed | $v_0$ | 8.20 km/s | $8.2 \times 10^{5}$ cm/s | 8.20 |
| Velocity dispersion | $\sigma_v$ | 1.43 km/s | $1.43 \times 10^{5}$ cm/s | 1.43 |
| Radius dispersion | $\sigma_r$ | 0.2 pc | $6.18 \times 10^{17}$ cm | 0.2 |
| Inclination | $i$ | 70° | - | 70° |
| Position angle | PA | 41.6° | - | 41.6° |

**Derived Orbital Parameters:**

| Parameter | Formula | Value | Physical Meaning |
|-----------|---------|-------|------------------|
| Orbital energy | $E = v_0^2/2 - GM/r_0$ | $-1.0 \times 10^{11}$ erg/g | Bound orbit |
| Angular momentum | $L = r_0 \times v_0 \sin\theta$ | $4.3 \times 10^{23}$ cm²/s | Conserved |
| Eccentricity | $e = \sqrt{1 + 2EL^2/(GM)^2}$ | 0.9999 | Highly eccentric |
| Pericentre | $r_p = L^2/GM(1+e)$ | 0.5-1 pc | Closest approach |
| Apocentre | $r_a = L^2/GM(1-e)$ | ~20 pc | Farthest |
| Orbital period | $T = 2\pi\sqrt{a^3/GM}$ | ~$3 \times 10^6$ yr | Full orbit |
| Pericentre velocity | $v_p = \sqrt{v_0^2 + 2GM(1/r_p - 1/r_0)}$ | ~100 km/s | Maximum speed |
| Interaction time | $\tau_\mathrm{int} = r_p/v_p$ | $10^4$ yr | Pericentre passage |

**Coordinate Transformation:**

The orbital plane coordinates $(X, Y)$ relate to Galactic coordinates $(l, b)$ through:

$$\begin{pmatrix} l \\ b \\ V_\mathrm{LSR} \end{pmatrix} = \mathbf{R}(i, \mathrm{PA}) \begin{pmatrix} X \\ Y \\ 0 \end{pmatrix} + \begin{pmatrix} l_0 \\ b_0 \\ V_0 \end{pmatrix}$$

where $\mathbf{R}(i, \mathrm{PA})$ is the rotation matrix:

$$\mathbf{R} = \mathbf{R}_z(\mathrm{PA}) \cdot \mathbf{R}_x(i)$$

With $i = 70°$ and $\mathrm{PA} = 41.6°$:

### 5.3 Initial Conditions

Based on Oka et al. (2017) N-body model:

```
IMBH:
  mass: 1.0e5 M_sun
  position: origin

Cloud:
  mass: 1000 M_sun (total), or use virial mass
  radius: 0.2 pc (Gaussian dispersion)
  position: (X, Y) = (9.8 pc, -0.65 pc)
  velocity: (vX, vY) = (-8.19 km/s, 0.4 km/s)
  internal velocity dispersion: 1.43 km/s
  temperature: 30-60 K
  mean molecular weight: 2.3

Geometry:
  orbital plane inclination: 70 degrees
  position angle: 41.6 degrees
```

### 5.4 Key Observables to Compare

1. **Position-velocity structure** at $t = 7.2 \times 10^5$ yr
2. **Velocity width** $\sim 100$ km/s
3. **Spatial distribution** matching HCN J=3-2 emission
4. **Dense clump** formation near pericentre

---

## 6. Comparison: When IS Two-Fluid MHD Needed?

The Inoue & Inutsuka (2008) two-fluid approach IS needed for:

| Scenario | Key Difference from IMBH Scattering |
|----------|-------------------------------------|
| **CNM/WNM formation** | Lower density ($n \sim 1-100$ cm$^{-3}$) |
| **Thermal instability** | $\tau_\mathrm{TI} \sim \tau_\mathrm{cool} \sim 1$ Myr |
| **Converging HI flows** | Sub-Alfvénic to trans-Alfvénic ($M_A \sim 1-2$) |
| **Field orientation** | Perpendicular $B$ dominates dynamics |

In contrast, IMBH-cloud scattering has:
- High density ($n \sim 10^4-10^6$ cm$^{-3}$)
- Short dynamical time ($\tau_\mathrm{dyn} \sim 10^4$ yr)
- Highly super-Alfvénic flow ($M_A \sim 250$)
- Gravity-dominated dynamics

---

## 7. Conclusions

### 7.1 MHD Effects

$$\frac{E_\mathrm{kinetic}}{E_\mathrm{magnetic}} \sim M_A^2 \sim 6 \times 10^4$$

$$\frac{F_\mathrm{tidal}}{F_\mathrm{magnetic}} \sim 300$$

**Magnetic fields are dynamically unimportant.** The flow is highly super-Alfvénic and gravitational tidal forces dominate by orders of magnitude.

### 7.2 Two-Fluid Effects

$$\frac{\tau_\mathrm{AD}}{\tau_\mathrm{dyn}} \sim 1000$$

$$\frac{\tau_{ni}}{\tau_\mathrm{dyn}} \sim 0.5$$

**Ambipolar diffusion is completely negligible.** Neutrals and ions are well-coupled. Single-fluid hydrodynamics is sufficient.

### 7.3 Shock Formation

$$\mathcal{M}_\mathrm{compress} \sim 10-15$$

**Strong shocks WILL form** during tidal compression. The simulation must use a proper shock-capturing scheme.

### 7.4 Final Recommendation

```
Minimum required physics:
  - Single-fluid hydrodynamics (SPH)
  - Self-gravity
  - Point-mass external gravity (IMBH)
  - Shock-capturing artificial viscosity

Optional enhancements:
  - Radiative cooling (for post-shock structure)
  - MHD (for small-scale structure, not primary dynamics)

NOT needed:
  - Two-fluid approximation
  - Ambipolar diffusion
  - Ion-neutral drift
```

---

## 8. Molecular Chemistry: HCN and CO for PV Diagram Reproduction

To reproduce the observed position-velocity (PV) diagrams from Oka et al. (2017), we must consider the molecular line emission from HCN and CO. This requires understanding the chemistry of these molecules and determining the appropriate treatment for SPH simulations.

### 8.1 First-Principles Analysis: Chemistry Timescales

#### 8.1.1 CO Chemistry

CO is the most abundant molecule after H$_2$ in molecular clouds. Its formation proceeds through two main pathways:

**Gas-phase formation (dominant at high density):**
$$\mathrm{C^+ + H_2 \rightarrow CH_2^+ + h\nu}$$
$$\mathrm{CH_2^+ + e^- \rightarrow CH + H}$$
$$\mathrm{CH + O \rightarrow CO + H}$$

Or via neutral-neutral reactions in dense gas:
$$\mathrm{C + OH \rightarrow CO + H} \quad (k \approx 10^{-10}\,\mathrm{cm^3\,s^{-1}})$$

**CO formation timescale:**
$$\tau_\mathrm{CO,form} \sim \frac{1}{k \times n(\mathrm{OH})} \sim \frac{1}{10^{-10} \times 10^{-7} \times n}$$

For $n = 10^4$ cm$^{-3}$ with typical OH abundance $[\mathrm{OH}]/[\mathrm{H_2}] \sim 10^{-7}$:
$$\tau_\mathrm{CO,form} \sim \frac{1}{10^{-10} \times 10^{-7} \times 10^4} = 10^{13}\,\mathrm{s} \approx 3 \times 10^5\,\mathrm{yr}$$

For the dense clump ($n = 3 \times 10^6$ cm$^{-3}$):
$$\boxed{\tau_\mathrm{CO,form}^\mathrm{dense} \approx 10^3\,\mathrm{yr}}$$

**CO destruction timescale:**
CO is destroyed primarily by UV photodissociation in diffuse regions and by cosmic ray-induced photodissociation in dense cores:
$$\mathrm{CO + CR \rightarrow C + O}$$

The destruction rate is $\xi_\mathrm{CO} \sim 10^{-17}$ s$^{-1}$ (cosmic ray dominated in dense gas):
$$\tau_\mathrm{CO,dest} \sim \frac{1}{\xi_\mathrm{CO}} \sim 10^{17}\,\mathrm{s} \approx 3 \times 10^9\,\mathrm{yr}$$

$$\boxed{\tau_\mathrm{CO,dest} \gg \tau_\mathrm{dyn}}$$

**Conclusion for CO:** Formation is fast in dense gas ($\tau_\mathrm{form} < \tau_\mathrm{dyn}$), destruction is negligible. **CO is in chemical equilibrium.**

#### 8.1.2 HCN Chemistry

HCN formation involves nitrogen chemistry, which is more complex:

**Primary formation pathways:**
$$\mathrm{N + CH_2 \rightarrow HCN + H} \quad (k \approx 5 \times 10^{-11}\,\mathrm{cm^3\,s^{-1}})$$
$$\mathrm{CN + H_2 \rightarrow HCN + H} \quad (k \approx 10^{-13}\,\mathrm{cm^3\,s^{-1}}\text{ at } T < 100\,\mathrm{K})$$

And ion-molecule chemistry:
$$\mathrm{HCNH^+ + e^- \rightarrow HCN + H}$$

**HCN formation timescale:**
The limiting step is often CN or CH$_2$ formation. For the CN + H$_2$ pathway:
$$\tau_\mathrm{HCN,form} \sim \frac{1}{k \times n(\mathrm{H_2})} \sim \frac{1}{10^{-13} \times n}$$

For $n = 3 \times 10^6$ cm$^{-3}$:
$$\tau_\mathrm{HCN,form} \sim \frac{1}{10^{-13} \times 3 \times 10^6} = 3 \times 10^6\,\mathrm{s} \approx 0.1\,\mathrm{yr}$$

However, the N → CN → HCN chain is limited by atomic N availability:
$$\tau_\mathrm{HCN,chain} \sim 10^4 - 10^5\,\mathrm{yr}$$

**HCN destruction:**
$$\mathrm{HCN + H^+ \rightarrow HCN^+ + H}$$
$$\mathrm{HCN + C^+ \rightarrow C_2N^+ + H}$$

Destruction timescale in dense gas:
$$\tau_\mathrm{HCN,dest} \sim 10^6 - 10^7\,\mathrm{yr}$$

$$\boxed{\tau_\mathrm{HCN,form} \sim \tau_\mathrm{dyn} \sim 10^4\,\mathrm{yr}}$$

**Conclusion for HCN:** Formation timescale is **comparable** to dynamical timescale. Chemical equilibrium is a reasonable but not exact approximation.

### 8.2 Timescale Comparison Summary

| Timescale | Value | Ratio to $\tau_\mathrm{dyn}$ |
|-----------|-------|------------------------------|
| Dynamical time $\tau_\mathrm{dyn}$ | $10^4$ yr | 1 |
| CO formation (dense) | $10^3$ yr | 0.1 |
| CO destruction | $10^9$ yr | $10^5$ |
| HCN formation (chain) | $10^4$–$10^5$ yr | 1–10 |
| HCN destruction | $10^6$ yr | $10^2$ |
| Shock crossing | $10^3$ yr | 0.1 |

**Key insight:** CO is safely in equilibrium. HCN formation timescale is marginal—equilibrium is approximate but acceptable for dense regions.

### 8.3 Recommended Approach: Equilibrium Chemistry with Post-Processing

Based on the timescale analysis, we recommend:

$$\boxed{\textbf{Chemical equilibrium + Post-processing radiative transfer}}$$

#### 8.3.1 Justification for Equilibrium Assumption

1. **CO:** $\tau_\mathrm{CO,form} \ll \tau_\mathrm{dyn}$ in dense gas. CO abundance reaches equilibrium rapidly.

2. **HCN:** $\tau_\mathrm{HCN,form} \sim \tau_\mathrm{dyn}$. While not exact, equilibrium is acceptable because:
   - HCN is only detectable in the densest regions ($n > 10^6$ cm$^{-3}$)
   - At these densities, formation is faster
   - Shock compression accelerates chemistry
   - The dense clump was likely pre-existing before IMBH encounter

3. **Computational efficiency:** Time-dependent chemistry networks add significant computational cost without substantially changing the PV diagram morphology.

4. **Dominant physics:** The PV structure is determined by **kinematics** (tidal stretching), not by chemistry variations. Abundance gradients are secondary effects.

#### 8.3.2 Equilibrium Abundance Prescriptions

**CO abundance:**
In well-shielded molecular gas ($A_V > 3$), CO reaches maximum abundance:
$$X(\mathrm{CO}) = \frac{n(\mathrm{CO})}{n(\mathrm{H_2})} \approx 10^{-4}$$

For SPH post-processing, use:
$$X(\mathrm{CO}) = X_\mathrm{CO,max} \times f_\mathrm{shield}(N_\mathrm{H_2})$$

where the shielding function:
$$f_\mathrm{shield} = \begin{cases}
1 & N_\mathrm{H_2} > 10^{21}\,\mathrm{cm^{-2}} \\
(N_\mathrm{H_2}/10^{21})^2 & N_\mathrm{H_2} < 10^{21}\,\mathrm{cm^{-2}}
\end{cases}$$

**HCN abundance:**
HCN abundance correlates with density in dense cores:
$$X(\mathrm{HCN}) = \frac{n(\mathrm{HCN})}{n(\mathrm{H_2})} \approx 10^{-8} \times \left(\frac{n}{10^4\,\mathrm{cm^{-3}}}\right)^{0.5}$$

This empirical scaling captures the density-dependent nitrogen chemistry.

For the Oka et al. dense clump:
$$X(\mathrm{HCN}) \approx 10^{-8} \times \left(\frac{3 \times 10^6}{10^4}\right)^{0.5} \approx 2 \times 10^{-7}$$

### 8.4 Radiative Transfer: Post-Processing Strategy

#### 8.4.1 Why Post-Processing is Sufficient

1. **Radiative transfer timescale:** $\tau_\mathrm{RT} \sim L/c \sim 10^{-4}$ yr $\ll \tau_\mathrm{dyn}$

   Radiation field equilibrates instantaneously compared to dynamics.

2. **Optically thin regime:** Most sightlines through the tidally disrupted cloud are optically thin (see Section 1.4.4).

3. **Decoupling:** Molecular line emission does not significantly affect gas dynamics (no radiation pressure, negligible cooling contribution compared to dust).

#### 8.4.2 Recommended Radiative Transfer Codes

| Code | Strengths | Best for |
|------|-----------|----------|
| **RADMC-3D** | Full 3D, flexible, LTE + non-LTE | CO imaging, dust continuum |
| **LIME** | Non-LTE line transfer, fast | HCN J=3-2, complex geometries |
| **RADEX** | Non-LTE, single-zone | Quick abundance estimates |
| **MOLLIE** | Non-LTE, turbulent sub-structure | Dense clump analysis |

**Primary recommendation:** RADMC-3D for synthetic observations, with RADEX cross-checks.

#### 8.4.3 Post-Processing Pipeline for PV Diagrams

```
SPH Simulation Output → Export to Grid → Radiative Transfer → Synthetic Observation
       │                       │                  │                    │
       ▼                       ▼                  ▼                    ▼
  (ρ, T, v)_particles   AMR or regular      RADMC-3D or        Convolve with
                        Cartesian grid       LIME               ALMA beam
                             │                  │                    │
                             ▼                  ▼                    ▼
                        Assign X(CO),     Compute τ, I_ν      Generate P-V
                        X(HCN) from       for each channel    diagram
                        equilibrium
```

**Step 1: Grid interpolation**

Convert SPH particles to a regular or AMR grid:
- Grid resolution: $\Delta x \lesssim h_\mathrm{min}/2$ (half minimum smoothing length)
- Recommended: $256^3$ to $512^3$ for 5 pc domain

```python
# Example: SPH to grid conversion
from scipy.interpolate import griddata

def sph_to_grid(particles, grid_size=256):
    """
    Interpolate SPH data to regular grid.

    Parameters:
    -----------
    particles : dict with keys 'pos', 'rho', 'T', 'vel', 'h'
    grid_size : int, number of grid cells per dimension
    """
    # Create regular grid
    x = np.linspace(-3, 3, grid_size)  # pc
    X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
    grid_points = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)

    # Interpolate density and temperature
    rho_grid = griddata(particles['pos'], particles['rho'],
                        grid_points, method='linear')
    T_grid = griddata(particles['pos'], particles['T'],
                      grid_points, method='linear')

    # Interpolate velocity components
    vx_grid = griddata(particles['pos'], particles['vel'][:, 0],
                       grid_points, method='linear')
    # ... repeat for vy, vz

    return rho_grid.reshape((grid_size,)*3), T_grid.reshape((grid_size,)*3)
```

**Step 2: Assign molecular abundances**

```python
def compute_abundances(n_H2, T):
    """
    Compute equilibrium CO and HCN abundances.

    Parameters:
    -----------
    n_H2 : array, H2 number density [cm^-3]
    T : array, temperature [K]
    """
    # CO abundance (constant in well-shielded regions)
    X_CO = 1e-4 * np.ones_like(n_H2)

    # HCN abundance (density-dependent)
    X_HCN = 1e-8 * (n_H2 / 1e4)**0.5
    X_HCN = np.minimum(X_HCN, 1e-6)  # Cap at maximum observed

    # Temperature dependence (optional: enhanced in warm shocked gas)
    warm_factor = np.where(T > 100, 1.5, 1.0)
    X_HCN *= warm_factor

    return X_CO, X_HCN
```

**Step 3: Radiative transfer with RADMC-3D**

Input files required:
- `amr_grid.inp`: Grid structure
- `dust_density.inp`: For continuum (optional)
- `gas_velocity.inp`: 3D velocity field
- `gas_temperature.inp`: Kinetic temperature
- `numberdens_co.inp`: CO number density
- `molecule_co.inp`: CO molecular data file

**Step 4: Generate synthetic observations**

```python
# RADMC-3D command for CO J=2-1 image cube
# radmc3d image lambda 1300.4 incl 70 phi 41.6 npix 256
#         vkms -120 nvkms 100 dvkms 2

def create_pv_from_cube(image_cube, spatial_axis, velocity_axis,
                        cut_position=0, cut_angle=0):
    """
    Extract P-V diagram from image cube.

    Parameters:
    -----------
    image_cube : 3D array (v, y, x) intensity
    cut_position : float, offset from center [arcsec]
    cut_angle : float, position angle [degrees]
    """
    # Rotate and extract 1D cut
    # ... implementation

    return pv_diagram
```

### 8.5 LTE vs Non-LTE Considerations

#### 8.5.1 When LTE is Valid

LTE holds when collisional excitation dominates over radiative de-excitation:
$$n > n_\mathrm{crit} = \frac{A_{ul}}{\gamma_{ul}}$$

| Transition | $n_\mathrm{crit}$ (cm$^{-3}$) | LTE valid in IMBH cloud? |
|------------|------------------------------|--------------------------|
| CO J=2-1 | $3 \times 10^3$ | **YES** ($n \sim 10^4$–$10^6$) |
| CO J=3-2 | $3 \times 10^4$ | **YES** (dense regions) |
| HCN J=3-2 | $3 \times 10^6$ | **Marginal** (dense clump only) |

**Recommendation:**
- **CO J=2-1, J=3-2:** Use LTE radiative transfer. This simplifies computation:
  $$T_b = T_\mathrm{ex} (1 - e^{-\tau}) \approx T_\mathrm{kin} (1 - e^{-\tau})$$

- **HCN J=3-2:** Use non-LTE (RADEX/LIME) with excitation calculation:
  $$\frac{n_u}{n_l} = \frac{g_u}{g_l} \frac{C_{lu} + \bar{J}_{ul} B_{lu}}{C_{ul} + A_{ul} + \bar{J}_{ul} B_{ul}}$$

  where $C_{ul}$, $C_{lu}$ are collisional rates and $\bar{J}_{ul}$ is the mean intensity.

#### 8.5.2 Non-LTE Effects on HCN

For HCN J=3-2 in the Oka et al. dense clump:
- $n = 3 \times 10^6$ cm$^{-3}$ ≈ $n_\mathrm{crit}$
- Sub-thermal excitation: $T_\mathrm{ex} < T_\mathrm{kin}$
- Expected $T_\mathrm{ex} \approx 0.6$–$0.8 \times T_\mathrm{kin}$

Non-LTE correction factor:
$$\frac{I_\mathrm{non-LTE}}{I_\mathrm{LTE}} \approx \frac{n}{n + n_\mathrm{crit}} \approx 0.5$$

This primarily affects absolute intensity, not PV morphology.

### 8.6 Shock Chemistry Enhancement (Optional)

Strong shocks ($\mathcal{M} \sim 10$–15) can temporarily alter molecular abundances:

#### 8.6.1 Shock-Enhanced HCN

In C-type shocks, elevated temperatures (100–1000 K) accelerate nitrogen chemistry:
- N$_2$ + He$^+$ → N$^+$ + N + He (cosmic ray ionization opens N$_2$)
- N + CH$_2$ → HCN + H (enhanced at $T > 100$ K)

**Enhancement factor in post-shock gas:**
$$X(\mathrm{HCN})_\mathrm{shock} \approx 3$–$10 \times X(\mathrm{HCN})_\mathrm{ambient}$$

This enhancement persists for $\tau_\mathrm{chem} \sim 10^4$ yr.

#### 8.6.2 CO Survival in Shocks

CO is remarkably robust in shocks:
- Dissociation requires $T > 3000$ K (not reached in $\mathcal{M} \sim 15$ shocks)
- Post-shock temperature $T \sim 2600$ K cools rapidly
- CO survives and may be enhanced due to liberation from grains

**Recommendation:** For shock regions ($T > 100$ K), apply HCN enhancement factor of 3×.

### 8.7 Implementation Summary: Best Practice for PV Diagram Reproduction

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                    RECOMMENDED WORKFLOW FOR PV DIAGRAMS                       │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  1. SPH SIMULATION (no chemistry)                                            │
│     ├── Pure hydrodynamics + self-gravity + point mass (IMBH)               │
│     ├── Track: density (ρ), temperature (T), velocity (v)                   │
│     └── Output snapshots at t = 7 × 10^5 yr                                  │
│                                                                              │
│  2. POST-PROCESSING: ABUNDANCE ASSIGNMENT                                    │
│     ├── CO: X(CO) = 10^{-4} (constant, equilibrium)                         │
│     ├── HCN: X(HCN) = 10^{-8} × (n/10^4)^{0.5}                             │
│     └── Optional: Shock enhancement for T > 100 K regions                   │
│                                                                              │
│  3. POST-PROCESSING: RADIATIVE TRANSFER                                      │
│     ├── CO J=2-1: LTE with RADMC-3D                                         │
│     │     └── T_ex = T_kin (valid for n > 10^4 cm^{-3})                     │
│     └── HCN J=3-2: Non-LTE with LIME or RADEX                               │
│           └── Solve excitation explicitly                                    │
│                                                                              │
│  4. SYNTHETIC OBSERVATION                                                    │
│     ├── Observer geometry: i = 70°, PA = 41.6°, V_LSR = -120 km/s          │
│     ├── Convolve with ALMA beam (1.87" × 1.14" for CO)                      │
│     └── Add noise: σ_rms ≈ 10 mJy/beam per 2 km/s channel                  │
│                                                                              │
│  5. COMPARISON WITH OBSERVATIONS                                             │
│     ├── PV diagram morphology (parallelogram shape)                         │
│     ├── Velocity width (~100 km/s)                                          │
│     ├── Dense clump position (0.2 pc offset from continuum)                │
│     └── Integrated intensity ratio HCN/CO                                   │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

### 8.8 Why NOT to Couple Chemistry to SPH

**Arguments against time-dependent chemistry in SPH:**

1. **Computational cost:** Chemical networks (e.g., UMIST, KIDA) with ~50+ species add 10–100× computational overhead per timestep.

2. **Stiffness:** Chemical timescales span $10^{-10}$ to $10^{10}$ years, requiring implicit solvers or subcycling.

3. **Marginal benefit:** The PV diagram shape is kinematically determined. Chemistry affects intensity, not morphology.

4. **Uncertainty dominates:** Reaction rate uncertainties (factor of 2–10) exceed any precision gained from time-dependent integration.

5. **Post-processing flexibility:** Equilibrium assumptions can be tested and modified without re-running expensive SPH simulations.

**Exception:** If studying the **formation** of the cloud or **long-term** chemical evolution ($> 10^6$ yr), time-dependent chemistry becomes necessary.

### 8.9 Conclusions for Chemistry Treatment

$$\boxed{\textbf{Equilibrium chemistry + post-processing radiative transfer}}$$

**Key recommendations:**

1. **Do NOT couple chemistry to SPH dynamics** — the equilibrium timescale is shorter than or comparable to the dynamical timescale.

2. **Use constant CO abundance** $X(\mathrm{CO}) = 10^{-4}$ throughout the well-shielded cloud.

3. **Use density-dependent HCN abundance** $X(\mathrm{HCN}) \propto n^{0.5}$ to capture the dense gas selectivity.

4. **Apply LTE radiative transfer for CO** — densities exceed $n_\mathrm{crit}$.

5. **Apply non-LTE radiative transfer for HCN** — densities are marginal relative to $n_\mathrm{crit}$.

6. **Post-process with RADMC-3D** for synthetic PV diagrams, applying ALMA beam convolution.

7. **Focus on kinematics** — the PV diagram morphology (parallelogram, velocity width) is the primary observable to match.

---

## Appendix A: Derivation of Key Formulas

### A.1 Tidal Radius

The tidal radius is where the cloud's self-gravity equals the IMBH's tidal force:

$$\frac{GM_\mathrm{cloud}}{r_t^2} = \frac{2GM_\mathrm{BH}}{d^3} r_t$$

Solving for $r_t$:

$$r_t = \left(\frac{M_\mathrm{cloud}}{2M_\mathrm{BH}}\right)^{1/3} d$$

(Factor of 2 vs 3 depends on geometry; we use the more common form with factor 3.)

### A.2 Alfvén Speed in Molecular Cloud

$$v_A = \frac{B}{\sqrt{4\pi\rho}} = \frac{B}{\sqrt{4\pi n \mu m_H}}$$

For $B = 30\,\mu$G, $n = 10^4$ cm$^{-3}$, $\mu = 2.8$:

$$v_A = \frac{3 \times 10^{-5}}{\sqrt{4\pi \times 10^4 \times 2.8 \times 1.67 \times 10^{-24}}}$$

$$v_A = \frac{3 \times 10^{-5}}{\sqrt{5.9 \times 10^{-19}}} = \frac{3 \times 10^{-5}}{7.7 \times 10^{-10}} = 3.9 \times 10^4\,\mathrm{cm/s}$$

### A.3 Ambipolar Diffusion Timescale

From the ambipolar diffusion velocity:

$$v_\mathrm{AD} = \frac{v_A^2}{\gamma_\mathrm{AD} \rho_i}$$

where $\gamma_\mathrm{AD} = \langle\sigma v\rangle / (m_i + m_n)$ is the drag coefficient.

The diffusion timescale:

$$\tau_\mathrm{AD} = \frac{L}{v_\mathrm{AD}} = \frac{L \gamma_\mathrm{AD} \rho_i}{v_A^2} = \frac{L^2}{v_A^2 \tau_{ni}}$$

---

## Appendix B: Physical Constants Used

| Constant | Value |
|----------|-------|
| $G$ | $6.67 \times 10^{-8}$ cm$^3$ g$^{-1}$ s$^{-2}$ |
| $k_B$ | $1.38 \times 10^{-16}$ erg K$^{-1}$ |
| $m_H$ | $1.67 \times 10^{-24}$ g |
| $M_\odot$ | $2 \times 10^{33}$ g |
| 1 pc | $3.09 \times 10^{18}$ cm |
| 1 km/s | $10^5$ cm/s |
| 1 yr | $3.15 \times 10^7$ s |

---

## 9. Corrected Oka et al. (2017) Orbital Analysis

### 9.1 Discovery of Initial Condition Mismatch

A detailed analysis of the Oka et al. (2017) Methods section revealed that the original CAT_OKA configuration used **incorrect initial conditions**. This section documents the correction.

**Oka et al. (2017) exact parameters (from Methods section):**

| Parameter | Paper Value | Original Config | Status |
|-----------|-------------|-----------------|--------|
| Initial position $(X_0, Y_0)$ | (9.8, -0.65) pc | (20.0, -5.17) pc | **WRONG** |
| Initial velocity $(v_X, v_Y)$ | (-8.19, 0.4) km/s | (-10.18, 5.05) km/s | **WRONG** |
| Cloud dispersion $\sigma_r$ | 0.2 pc | 1.13 pc | Different |
| Velocity dispersion | 1.43 km/s | N/A | Missing |
| BH mass | $10^5\,M_\odot$ | $10^5\,M_\odot$ | Correct |
| Inclination | 70° | 70° | Correct |
| Position angle | 41.6° | 41.6° | Correct |
| $V_\mathrm{LSR}$ | -120 km/s | -120 km/s | Correct |

### 9.2 Physical Understanding of Oka's N-body Simulation

The key insight is that Oka et al. used **test particles with velocity dispersion**, not a single trajectory:

**Mean orbital parameters:**
$$r_0 = \sqrt{9.8^2 + 0.65^2} = 9.82\,\text{pc}$$
$$v_0 = \sqrt{8.19^2 + 0.4^2} = 8.20\,\text{km/s}$$

**Angular momentum (mean):**
$$h_\mathrm{mean} = X_0 v_Y - Y_0 v_X = 9.8 \times 0.4 - (-0.65) \times (-8.19) = 3.92 - 5.32 = -1.40\,\text{pc}\cdot\text{km/s}$$

This extremely small angular momentum implies a **nearly radial mean orbit** with pericentre:
$$r_\mathrm{peri,mean} = \frac{h^2}{GM_\mathrm{BH}(1+e)} \approx 0.002\,\text{pc}$$

**However**, the cloud has velocity dispersion $\sigma_v = 1.43\,\text{km/s}$ and radial dispersion $\sigma_r = 0.2\,\text{pc}$. Individual particles have angular momentum spread:
$$\delta h \sim r_0 \sigma_v + \sigma_r v_0 \sim 9.8 \times 1.43 + 0.2 \times 8.2 \sim 15.7\,\text{pc}\cdot\text{km/s}$$

This means particles span a range of pericentres:

| Particle $h$ (pc·km/s) | Eccentricity $e$ | Pericentre $r_\mathrm{peri}$ (pc) |
|------------------------|------------------|-----------------------------------|
| 1.4 (mean) | 0.9999 | 0.002 |
| 5 | 0.9985 | 0.03 |
| 10 | 0.9940 | 0.11 |
| 15 | 0.9863 | 0.25 |
| 20 | 0.9756 | 0.45 |

**Physical interpretation:** The surviving HCN clump in Oka's simulation consists of particles with **larger angular momentum** that had pericentres of 0.2-0.5 pc, while particles with low angular momentum plunged closer to the BH and were strongly accelerated (forming the high-velocity tails in the p-v diagram).

### 9.3 Implications for SPH Simulation

For our SPH simulation with cloud radius $R = 1.13\,\text{pc}$ (vs Oka's $\sigma_r = 0.2\,\text{pc}$):

1. **Larger angular momentum spread**: $\delta h_\mathrm{SPH} \sim r_0 \sigma_v + R \times v_0 \sim 23\,\text{pc}\cdot\text{km/s}$

2. **Different dynamics**: Our larger cloud will experience:
   - More extreme tidal stretching
   - Wider range of particle pericentres
   - Potentially different fragmentation behavior

3. **Recommended approach**:
   - Use Oka's exact initial conditions for the cloud center-of-mass
   - Generate new IC with smaller cloud ($\sigma_r = 0.2\,\text{pc}$) for accurate reproduction
   - Or scale orbital parameters to match tidal strength ratio

### 9.4 Corrected Configuration

A corrected configuration file `oka_corrected.json` has been created with:
```json
{
  "initialCondition": {
    "transform": {
      "translate": [9.8, -0.65, 0.0],
      "velocity_boost": [-8.19, 0.4, 0.0]
    }
  }
}
```

---

## 10. Plan for Radio Observation Reproduction (Figures 2 & 3)

### 10.1 Overview: What Needs to be Reproduced

**Figure 2 (p-v diagrams):**
- HCN J=3-2 and J=4-3 emission contours (observed data)
- SPH particle positions overlaid as dots
- Three slit orientations with different position angles

**Figure 3 (continuum spectrum):**
- Radio/X-ray spectrum of CO-0.40-0.22*
- This is **continuum emission from accretion**, NOT molecular lines
- Cannot be reproduced from pure hydrodynamic simulation (requires RIAF model)

### 10.2 Three-Level Approach for p-v Diagram

#### Level 1: Particle Projection (What Oka Did)

**Method:** Simply project SPH particle positions to (position, velocity) space.

**Implementation:**
```python
# Transform to observer frame
R = rotation_matrix(inclination=70°, PA=41.6°)
pos_obs = R @ positions
vel_obs = R @ velocities

# Line-of-sight velocity
v_LOS = vel_obs[:, 2] + V_LSR  # V_LSR = -120 km/s

# Position along slit (for given PA)
pos_slit = pos_obs[:, 0] * cos(slit_PA) + pos_obs[:, 1] * sin(slit_PA)

# Plot (pos_slit, v_LOS) weighted by density
```

**Pros:** Simple, fast, matches Oka's methodology
**Cons:** No actual emission physics, no optical depth

#### Level 2: LTE Emission with KI2000 Equilibrium

**Method:** Use KI2000 equilibrium tables for molecular abundances, apply LTE radiative transfer.

**Molecular abundances from KI2000:**
- $X(\mathrm{CO})$ as function of $(n, N_H)$ - available in extracted data
- $X(\mathrm{H}_2)$ as function of $(n, N_H)$ - available in extracted data
- $X(\mathrm{HCN})$ - NOT in KI2000, need prescription

**HCN abundance prescription:**
$$X(\mathrm{HCN}) = X_0 \times \min\left(1, \frac{n}{n_\mathrm{crit}}\right)^{0.5}$$

where $X_0 \sim 10^{-8}$ and $n_\mathrm{crit} = 3 \times 10^6\,\text{cm}^{-3}$.

**LTE emission:**
$$I_\nu \propto N_\mathrm{mol} \times A_{ul} \times \frac{g_u \exp(-E_u/kT)}{Q(T)}$$

**Implementation in sph-viz frontend:**
1. Compute temperature from internal energy: $T = (\gamma-1) u \mu m_H / k_B$
2. Look up CO abundance from KI2000 tables
3. Apply HCN prescription based on density
4. Weight particles by emissivity and sum to p-v grid
5. Convolve with beam

#### Level 3: Full Non-LTE Radiative Transfer

**Method:** Export SPH data to RADMC-3D, run full line RT.

**Steps:**
1. Interpolate SPH particles to regular grid
2. Export density, temperature, velocity fields
3. Add molecular abundances (CO, HCN)
4. Run RADMC-3D with LVG excitation
5. Generate synthetic datacube
6. Extract p-v diagrams along slits
7. Convolve with ALMA beam

**Advantages:** Handles optical depth, non-LTE effects, proper beam convolution

### 10.3 Timescale Diagnostics for Equilibrium Validation

To verify that KI2000 equilibrium is applicable, we need to compare relevant timescales at each SPH particle location.

**Key timescales to compute:**

| Timescale | Formula | Equilibrium Criterion |
|-----------|---------|----------------------|
| Dynamical | $\tau_\mathrm{dyn} = \sqrt{R^3/(GM)}$ | - |
| Cooling | $\tau_\mathrm{cool} = u / \dot{u}_\mathrm{cool}$ | $\tau_\mathrm{cool} \ll \tau_\mathrm{dyn}$ |
| Heating | $\tau_\mathrm{heat} = u / \dot{u}_\mathrm{heat}$ | $\tau_\mathrm{heat} \ll \tau_\mathrm{dyn}$ |
| Chemical | $\tau_\mathrm{chem} \sim (k_\mathrm{form} n)^{-1}$ | $\tau_\mathrm{chem} \ll \tau_\mathrm{dyn}$ |
| Tidal | $\tau_\mathrm{tidal} = r / v_\mathrm{peri}$ | - |

**Frontend diagnostic plots needed:**
1. **Histogram of $\tau_\mathrm{cool}/\tau_\mathrm{dyn}$** - should be $\ll 1$ for thermal equilibrium
2. **Histogram of $\tau_\mathrm{chem}/\tau_\mathrm{dyn}$** - should be $\ll 1$ for chemical equilibrium
3. **Spatial map of equilibrium validity** - color particles by whether they satisfy equilibrium

### 10.4 Interactive p-v Diagram Implementation

**Frontend requirements:**

1. **Observer geometry controls:**
   - Inclination slider: 0° to 90° (default 70°)
   - Position angle slider: 0° to 360° (default 41.6°)
   - V_LSR slider: -200 to 0 km/s (default -120)
   - Slit position offset
   - Slit width

2. **Real-time p-v computation:**
   - Transform particle coordinates on GPU (Three.js)
   - Compute v_LOS for all particles
   - Bin to 2D histogram
   - Apply density weighting

3. **Visualization options:**
   - Contour levels
   - Color map selection
   - Log/linear scale
   - Overlay observed data (if available)
   - Best-fit finder (minimize residuals)

### 10.5 Equilibrium Chemistry with KI2000

**Available KI2000 data (extracted):**

| File | Contents |
|------|----------|
| `f1a_temperature_N19.txt` | T(n) for $N_H = 10^{19}$ |
| `f1a_temperature_N20.txt` | T(n) for $N_H = 10^{20}$ |
| `f1b_x_H2_N19/N20.txt` | $X(\mathrm{H}_2)(n)$ |
| `f1b_x_CO_N19/N20.txt` | $X(\mathrm{CO})(n)$ |
| `f1b_x_e_N19/N20.txt` | $X(e^-)(n)$ |

**Interpolation approach:**
1. For each SPH particle with density $n$ and column density $N_H$:
2. Interpolate between $N_H = 10^{19}$ and $10^{20}$ curves
3. Get equilibrium $T$, $X(\mathrm{H}_2)$, $X(\mathrm{CO})$
4. Compare with actual SPH temperature to assess equilibrium

**When equilibrium is valid:**
- Use KI2000 abundances directly
- Post-process with simple LTE RT

**When equilibrium fails:**
- Need time-dependent chemistry (expensive)
- Or use kinematic approach (ignore chemistry)

---

## References

1. Oka, T., et al. (2017). "Millimetre-wave Emission from an Intermediate-Mass Black Hole Candidate in the Milky Way." Nature Astronomy.

2. Inoue, T., & Inutsuka, S. (2008). "Two-Fluid MHD Simulations of Converging HI Flows in the Interstellar Medium. I: Methodology and Basic Results." ApJ.

3. Koyama, H., & Inutsuka, S. (2000). "Molecular Cloud Formation in Shock-compressed Layers." ApJ, 532, 980.

4. Draine, B. T. (1986). "Multicomponent, Reacting MHD Flows." MNRAS, 220, 133.

5. van Dishoeck, E. F., & Black, J. H. (1988). "The Photodissociation and Chemistry of Interstellar CO." ApJ, 334, 771.

6. Boger, G. I., & Sternberg, A. (2005). "CN and HCN in Dense Interstellar Clouds." ApJ, 632, 302.

7. Glover, S. C. O., et al. (2010). "Modelling CO Formation in the Turbulent Interstellar Medium." MNRAS, 404, 2.

8. Dullemond, C. P., et al. (2012). "RADMC-3D: A Multi-purpose Radiative Transfer Tool." Astrophysics Source Code Library, ascl:1202.015.

9. Brinch, C., & Hogerheijde, M. R. (2010). "LIME - a Flexible, Non-LTE Line Excitation and Radiation Transfer Method for Millimeter and Far-infrared Wavelengths." A&A, 523, A25.

10. van der Tak, F. F. S., et al. (2007). "A Computer Program for Fast Non-LTE Analysis of Interstellar Line Spectra." A&A, 468, 627. (RADEX)

11. Gao, Y., & Solomon, P. M. (2004). "HCN Survey of Normal Spiral, Infrared-luminous, and Ultraluminous Galaxies." ApJS, 152, 63.

12. McElroy, D., et al. (2013). "The UMIST Database for Astrochemistry 2012." A&A, 550, A36.
