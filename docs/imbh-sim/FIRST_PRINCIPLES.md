# First Principles: KI2000 Bonnor-Ebert Cloud

## 1. Koyama & Inutsuka (2000) Thermal Equilibrium

The ISM thermal equilibrium is determined by the balance of heating and cooling.

### Embedded Data (N_H = 10^19 cm^-2):

| n [cm^-3] | T_eq [K] | P/k_B [K cm^-3] | Phase |
|-----------|----------|-----------------|-------|
| 0.1 | ~630 | ~63 | WNM |
| 1 | ~60 | ~60 | Unstable |
| 10 | ~20 | ~200 | Unstable |
| 50 | ~9.4 | ~470 | CNM |
| 100 | ~8.5 | ~850 | CNM |
| 1000 | ~7.3 | ~7300 | CNM/Molecular |

**Key insight**: T ~ 10 K requires n > 50 cm^-3 (cold neutral medium).

## 2. Cloud Size Scaling

For a pressure-truncated Bonnor-Ebert sphere:
```
R_cloud ~ xi_s * r_0
r_0 = c_s / sqrt(4 * pi * G * rho_c)
```

where xi_s ~ 3-6 is the dimensionless truncation radius.

### Implication for T ~ 10 K clouds:

At n_c = 1000 cm^-3, T ~ 7 K:
```
c_s ~ 0.17 km/s
rho_c ~ 31 M_sun/pc^3 (code units)
r_0 ~ 0.17 / sqrt(4 * pi * 1 * 31) ~ 0.009 pc
R_cloud ~ 3 * 0.009 ~ 0.03 pc
```

**Result**: A cold (T ~ 10 K) KI2000 cloud has R ~ 0.03 pc, not 1 pc!

### To get R ~ 1 pc with KI2000:

Need r_0 ~ 0.3 pc, which requires:
```
rho_c ~ (0.009/0.3)^2 * 31 ~ 0.03 M_sun/pc^3
n_c ~ 1 cm^-3
T_eq ~ 60 K (warm neutral medium, NOT molecular!)
```

## 3. Physical Reality

There are TWO physically distinct scenarios:

### Scenario A: Cold Molecular Cloud (T ~ 10 K)
- n_center ~ 1000 cm^-3
- R ~ 0.03 pc (very compact core)
- M ~ 0.1-1 M_sun
- Embedded in larger warm envelope

### Scenario B: Warm Diffuse Cloud (T ~ 60 K)
- n_center ~ 1 cm^-3
- R ~ 1 pc
- M ~ 1000 M_sun
- WNM conditions

The original Oka (2017) HVCC has:
- M ~ 40 M_sun
- R ~ 0.17 pc
- n_center ~ 10^5 cm^-3 (molecular core)

This is Scenario A - a compact molecular core.

## 4. Practical Implementation

### For molecular cloud simulation (our target):

Use the `isothermal_bonnor_ebert` sample with FIXED T = 10 K:
- Specify M_cloud, T_cloud, xi_s
- R_cloud is computed from Jeans analysis
- Does NOT use KI2000 variable T_eq(n)

### For full KI2000 equilibrium:

Use `bonnor_ebert_ki2000`:
- Specify n_center, P_ext
- R and M are determined by profile integration
- Accept the resulting (small) cloud size for T ~ 10 K

## 5. Recommended Configuration

For a T ~ 10 K molecular cloud simulation:

**Option 1: Compact KI2000 cloud**
```json
{
    "sample": "bonnor_ebert_ki2000",
    "rho_center_nH": 1000.0,
    "P_ext_K_cm3": 2000.0
}
```
Result: R ~ 0.03 pc, M ~ 0.001 M_sun, T_center ~ 7 K

**Option 2: Isothermal sphere (fixed T)**
```json
{
    "sample": "isothermal_bonnor_ebert",
    "M_cloud": 40.0,
    "T_cloud": 10.0,
    "xi_s": 5.0
}
```
Result: R determined by Jeans length, T = 10 K fixed

## 6. Conclusion

The KI2000 barotropic EOS naturally produces compact (R << 1 pc) clouds at cold temperatures (T ~ 10 K). This is physically correct - cold molecular cores are embedded in larger warm envelopes.

For large (R ~ 1 pc) clouds, either:
1. Accept warm (T ~ 60 K) temperatures from KI2000, or
2. Use isothermal (fixed T) approximation
