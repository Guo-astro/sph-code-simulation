# 1D Diffusive Instability Study

## Overview

This test case studies the **secular diffusive instability** in GSPH without grad-h correction.

The theory predicts that GSPH without the Ω (grad-h) correction exhibits effective mass diffusion:

$$D_{\text{eff}} = \epsilon \cdot c_s \cdot h$$

with growth rate:

$$\Gamma = \frac{\epsilon \cdot c_s \cdot h}{L^2}$$

## Physical Setup

### 1D Isothermal Slab (Self-Gravitating)

A self-gravitating isothermal slab in hydrostatic equilibrium:

- **Density profile**: $\rho(x) = \rho_0 \cdot \text{sech}^2(x/H)$
- **Half-thickness**: $H$ (scale height)
- **Sound speed**: $c_s = \sqrt{P/\rho}$ (constant for isothermal)
- **1D Slab Gravity**: Uses Gauss's Law: $g(x) = -2\pi G \cdot \text{sign}(x) \cdot \Sigma(|x|)$

where $\Sigma(|x|) = \int_0^{|x|} \rho(x') dx'$ is the enclosed surface density.

**Note:** The code implements proper 1D slab gravity (not 3D point-mass gravity). This is computed by:
1. Sorting particles by $|x|$
2. Computing cumulative enclosed mass
3. Applying $g_i = -2\pi G \cdot \text{sign}(x_i) \cdot M_{\text{enclosed}}$

### Hydrostatic Equilibrium

For an isothermal slab:
$$\frac{dP}{dx} = -\rho g = -2\pi G \rho \Sigma(<x)$$

With $P = c_s^2 \rho$:
$$c_s^2 \frac{d\rho}{dx} = -2\pi G \rho \int_0^x \rho dx'$$

Solution: $\rho(x) = \rho_0 \, \text{sech}^2(x/H)$ where $H = c_s^2/(2\pi G \rho_0 H)$

## Test Cases

### 1. Stability Comparison

Compare GSPH with and without grad-h correction:

| Configuration | Expected Behavior |
|--------------|-------------------|
| `use_gradh: true` | Stable equilibrium |
| `use_gradh: false` | Exponential density growth |

### 2. Growth Rate Measurement

Measure $\Gamma$ from central density evolution:
$$\rho_c(t) = \rho_c(0) \cdot e^{\Gamma t}$$

### 3. Resolution Scaling

Verify $\Gamma \propto h \propto N^{-1}$ (in 1D):

| N particles | Expected $\Gamma / \Gamma_0$ |
|-------------|------------------------------|
| 100 | 1.0 |
| 200 | 0.5 |
| 400 | 0.25 |

### 4. Sound Speed Scaling

Verify $\Gamma \propto c_s$:

| $c_s / c_{s,0}$ | Expected $\Gamma / \Gamma_0$ |
|-----------------|------------------------------|
| 0.5 | 0.5 |
| 1.0 | 1.0 |
| 2.0 | 2.0 |

## Usage

### Quick Start

```bash
# Build with DIM=1
cd /Users/guo/Downloads/sphcode
mkdir -p build && cd build
cmake -DDIM=1 ..
make -j4

# Run stability comparison
make -f sample/diffusive_instability/Makefile.diffusive_instability compare_all
```

### Individual Tests

```bash
# With grad-h (stable)
./build/sph sample/diffusive_instability/config/presets/slab_gsph_gradh.json

# Without grad-h (unstable)
./build/sph sample/diffusive_instability/config/presets/slab_gsph_no_gradh.json

# Analyze results
python sample/diffusive_instability/scripts/analyze_growth_rate.py
```

## Expected Results

### Central Density Evolution

```
                   │
ln(ρ_c/ρ_c0)       │     Without grad-h: exponential growth
                   │        /
                   │       /  slope = Γ
                   │      /
                   │     /
                   │────/────────── With grad-h: stable
                   │   /
                   └───────────────────────────── time
```

### Theoretical vs Measured

For typical parameters ($\epsilon \approx 0.4$, $c_s = 1$, $h/L \approx 0.1$):

$$\Gamma_{\text{theory}} = \frac{0.4 \times 1.0 \times 0.1 \times L}{L^2} = \frac{0.04}{L}$$

## File Structure

```
diffusive_instability/
├── README.md                    # This file
├── Makefile.diffusive_instability
├── config/
│   └── presets/
│       ├── slab_gsph_gradh.json      # Stable (with grad-h)
│       ├── slab_gsph_no_gradh.json   # Unstable (without grad-h)
│       ├── slab_resolution_N100.json # Resolution study
│       ├── slab_resolution_N200.json
│       └── slab_resolution_N400.json
├── scripts/
│   ├── generate_slab.py              # Generate initial conditions
│   ├── analyze_growth_rate.py        # Measure Γ from simulation
│   ├── plot_comparison.py            # Compare stable vs unstable
│   └── validate_diffusion_theory.py  # Full theory validation
├── results/
│   └── (simulation outputs)
└── docs/
    └── THEORY.md                     # Detailed derivation
```

## References

1. Springel & Hernquist (2002) - Grad-h correction for SPH
2. Hopkins (2013) - A general class of SPH methods
3. This derivation: `sample/gradh_study/docs/SECULAR_INSTABILITY_DERIVATION.md`
