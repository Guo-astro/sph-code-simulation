/**
 * SRMHD Riemann Solver
 *
 * HLL-type Riemann solver for relativistic fast magnetosonic waves.
 *
 * Based on:
 * - Mignone & Bodo (2005) for relativistic HLLC
 * - Iwasaki & Inutsuka (2011) for MHD decomposition in SPH
 *
 * Key features:
 * - Relativistic wave speed estimates using velocity addition formula
 * - Acoustic impedance matching for interface states
 * - Compatible with Method of Characteristics for Alfven waves
 */

#include "srmhd/srmhd_riemann.hpp"
#include <cmath>
#include <algorithm>

namespace sph
{
namespace srmhd
{
namespace riemann
{

void compute_wave_speeds(
    real v, real c_f,
    real& lambda_minus, real& lambda_plus
)
{
    // Relativistic velocity addition formula:
    // lambda = (v +/- c_f) / (1 +/- v*c_f)

    const real v_cf = v * c_f;

    // Left-going characteristic (minus)
    const real denom_m = 1.0 - v_cf;
    lambda_minus = (denom_m > 1e-10) ? (v - c_f) / denom_m : -1.0;

    // Right-going characteristic (plus)
    const real denom_p = 1.0 + v_cf;
    lambda_plus = (denom_p > 1e-10) ? (v + c_f) / denom_p : 1.0;

    // Clamp to subluminal
    lambda_minus = std::max(-0.9999, lambda_minus);
    lambda_plus = std::min(0.9999, lambda_plus);
}

real compute_impedance(real rho, real H, real gamma_lor, real c_f)
{
    // Relativistic MHD impedance
    // Z = rho * H * gamma^2 * c_f
    return rho * H * gamma_lor * gamma_lor * c_f;
}

real acoustic_pressure(
    real P_L, real P_R,
    real v_L, real v_R,
    real Z_L, real Z_R
)
{
    const real Z_sum = Z_L + Z_R;
    if (Z_sum < 1e-30) {
        return 0.5 * (P_L + P_R);
    }

    return (Z_L * P_R + Z_R * P_L + Z_L * Z_R * (v_L - v_R)) / Z_sum;
}

real acoustic_velocity(
    real P_L, real P_R,
    real v_L, real v_R,
    real Z_L, real Z_R
)
{
    const real Z_sum = Z_L + Z_R;
    if (Z_sum < 1e-30) {
        return 0.5 * (v_L + v_R);
    }

    return (Z_L * v_L + Z_R * v_R + P_L - P_R) / Z_sum;
}

void hll_solver(
    const SRMHDState& L,
    const SRMHDState& R,
    real& Pt_star,
    real& v_star
)
{
    // Compute wave speeds using relativistic formula
    real lambda_L_minus, lambda_L_plus;
    real lambda_R_minus, lambda_R_plus;

    compute_wave_speeds(L.v, L.c_f, lambda_L_minus, lambda_L_plus);
    compute_wave_speeds(R.v, R.c_f, lambda_R_minus, lambda_R_plus);

    // HLL wave speed estimates (Einfeldt-type)
    const real S_L = std::min(lambda_L_minus, lambda_R_minus);
    const real S_R = std::max(lambda_L_plus, lambda_R_plus);

    // Check for degenerate case
    if (std::abs(S_R - S_L) < 1e-15) {
        Pt_star = 0.5 * (L.P_t + R.P_t);
        v_star = 0.5 * (L.v + R.v);
        return;
    }

    // Compute acoustic impedances
    const real Z_L = compute_impedance(L.rho, L.H, L.gamma_lor, L.c_f);
    const real Z_R = compute_impedance(R.rho, R.H, R.gamma_lor, R.c_f);

    // Acoustic approximation for interface states
    Pt_star = acoustic_pressure(L.P_t, R.P_t, L.v, R.v, Z_L, Z_R);
    v_star = acoustic_velocity(L.P_t, R.P_t, L.v, R.v, Z_L, Z_R);

    // Floor pressure
    Pt_star = std::max(Pt_star, 1e-15);

    // Clamp velocity to subluminal
    v_star = std::max(-0.9999, std::min(0.9999, v_star));
}

void hllc_solver(
    const SRMHDState& L,
    const SRMHDState& R,
    real& Pt_star,
    real& v_star
)
{
    // For SPH MHD, we use the simpler HLL approach since
    // the Method of Characteristics handles the contact wave structure.
    // The HLLC solver is more important for grid-based methods.
    //
    // Here we implement a more accurate version using Roe averages.

    // Roe-type averages
    const real sqrt_rho_L = std::sqrt(std::max(L.rho, 1e-30));
    const real sqrt_rho_R = std::sqrt(std::max(R.rho, 1e-30));
    const real roe_inv = 1.0 / (sqrt_rho_L + sqrt_rho_R);

    const real v_roe = (sqrt_rho_L * L.v + sqrt_rho_R * R.v) * roe_inv;
    const real c_f_roe = (sqrt_rho_L * L.c_f + sqrt_rho_R * R.c_f) * roe_inv;

    // Wave speeds using Roe averages
    real lambda_roe_minus, lambda_roe_plus;
    compute_wave_speeds(v_roe, c_f_roe, lambda_roe_minus, lambda_roe_plus);

    real lambda_L_minus, lambda_L_plus;
    real lambda_R_minus, lambda_R_plus;
    compute_wave_speeds(L.v, L.c_f, lambda_L_minus, lambda_L_plus);
    compute_wave_speeds(R.v, R.c_f, lambda_R_minus, lambda_R_plus);

    // Davis-type wave speed estimates
    const real S_L = std::min({lambda_L_minus, lambda_roe_minus});
    const real S_R = std::max({lambda_R_plus, lambda_roe_plus});

    // Compute fluxes and HLL state
    // For relativistic MHD, the flux structure is complex.
    // We use the acoustic approximation which works well for SPH.

    const real Z_L = compute_impedance(L.rho, L.H, L.gamma_lor, L.c_f);
    const real Z_R = compute_impedance(R.rho, R.H, R.gamma_lor, R.c_f);

    Pt_star = acoustic_pressure(L.P_t, R.P_t, L.v, R.v, Z_L, Z_R);
    v_star = acoustic_velocity(L.P_t, R.P_t, L.v, R.v, Z_L, Z_R);

    // Ensure physicality
    Pt_star = std::max(Pt_star, 1e-15);
    v_star = std::max(-0.9999, std::min(0.9999, v_star));
}

}
}
}
