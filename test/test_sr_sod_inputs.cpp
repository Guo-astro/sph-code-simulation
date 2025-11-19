// Test C++ Riemann solver with SR Sod initial conditions
// Compare directly with Python results

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>

using real = double;

// Identical to test_iterative_riemann_solver.cpp implementation
class SRRiemannSolver {
private:
    const real m_gamma;
    const real m_c_speed;
    
    void get_velocity_behind_wave(
        real P, real na, real Pa, real Ha, real csa, real va, real vta, real wa, real sign,
        real& n, real& H, real& cs, real& vel)
    {
        const real c = m_c_speed;
        const real gamma_c = m_gamma;
        
        if (P > Pa) {
            // SHOCK
            const real a = 1.0 + (gamma_c - 1.0) * (Pa - P) / (gamma_c * P);
            const real b_term = -(gamma_c - 1.0) * (Pa - P) / (gamma_c * P);
            const real c_term = Ha * (Pa - P) / na - Ha * Ha;
            const real discriminant = b_term * b_term - 4.0 * a * c_term;
            
            if (discriminant < 0.0) {
                n = na;
                H = Ha;
                cs = csa;
                vel = va;
                return;
            }
            
            H = (-b_term + std::sqrt(discriminant)) / (2.0 * a);
            if (H <= 0.0) {
                n = na;
                H = Ha;
                cs = csa;
                vel = va;
                return;
            }
            
            n = gamma_c * P / ((gamma_c - 1.0) * (H - 1.0));
            const real j2 = -(P - Pa) / (H / n - Ha / na);
            if (j2 < 0.0) {
                n = na;
                H = Ha;
                cs = csa;
                vel = va;
                return;
            }
            const real j = sign * std::sqrt(j2);
            
            const real Na = wa * na;
            const real a_shock = j2 + Na * Na;
            const real b_shock = -va * Na * Na;
            const real sqrt_term = std::sqrt(1.0 + na * na / j2);
            const real vshock = (-b_shock + sign * j2 * sqrt_term) / a_shock;
            const real gamma_shock = 1.0 / std::sqrt(1.0 - (vshock / c) * (vshock / c));
            
            const real a_vel = gamma_shock * (P - Pa) / j + Ha * wa * va;
            const real b_vel = Ha * wa + (P - Pa) * (gamma_shock * va / j + 1.0 / (na * wa));
            vel = a_vel / b_vel;
            
            const real v_limit = 0.99 * c;
            if (std::abs(vel) > v_limit) {
                vel = std::copysign(v_limit, vel);
            }
            
            cs = std::sqrt(gamma_c * P / (n * H));
        }
        else {
            // RAREFACTION
            const real K = Pa / std::pow(na, gamma_c);
            n = std::pow(P / K, 1.0 / gamma_c);
            const real u = P / ((gamma_c - 1.0) * n);
            H = 1.0 + u / (c * c) + P / (n * c * c);
            cs = std::sqrt(gamma_c * P / (n * H));
            
            const real sqgl1 = std::sqrt(gamma_c - 1.0);
            const real term1 = (1.0 + va / c) / (1.0 - va / c);
            const real term2 = (sqgl1 + csa / c) / (sqgl1 - csa / c);
            const real term3 = (sqgl1 - cs / c) / (sqgl1 + cs / c);
            const real exponent = -sign * 2.0 / sqgl1;
            const real A = term1 * std::pow(term2 * term3, exponent);
            
            vel = c * (A - 1.0) / (A + 1.0);
            
            const real v_limit = 0.99 * c;
            if (std::abs(vel) > v_limit || !std::isfinite(vel)) {
                vel = std::copysign(std::min(std::abs(va), v_limit), va);
            }
        }
    }
    
public:
    SRRiemannSolver(real gamma, real c) : m_gamma(gamma), m_c_speed(c) {}
    
    void solve(real Pl, real nl, real vl, real Pr, real nr, real vr,
               real& pstar, real& vstar, int& iterations)
    {
        const real c = m_c_speed;
        const real gamma_c = m_gamma;
        
        // Derived quantities
        const real c2 = c * c;
        const real ul = Pl / ((gamma_c - 1.0) * nl);
        const real ur = Pr / ((gamma_c - 1.0) * nr);
        const real Hl = 1.0 + ul / c2 + Pl / (nl * c2);
        const real Hr = 1.0 + ur / c2 + Pr / (nr * c2);
        const real csl = std::sqrt(gamma_c * Pl / (nl * Hl));
        const real csr = std::sqrt(gamma_c * Pr / (nr * Hr));
        const real wl = 1.0 / std::sqrt(1.0 - (vl / c) * (vl / c));
        const real wr = 1.0 / std::sqrt(1.0 - (vr / c) * (vr / c));
        
        // dvel function
        auto dvel = [&](real p_test) -> real {
            real nl_star, Hl_star, csl_star, vl_star;
            real nr_star, Hr_star, csr_star, vr_star;
            
            get_velocity_behind_wave(p_test, nl, Pl, Hl, csl, vl, 0.0, wl, -1.0,
                                     nl_star, Hl_star, csl_star, vl_star);
            get_velocity_behind_wave(p_test, nr, Pr, Hr, csr, vr, 0.0, wr, +1.0,
                                     nr_star, Hr_star, csr_star, vr_star);
            return vl_star - vr_star;
        };
        
        // Adaptive bracketing
        real p_min = 0.5 * (Pl + Pr);
        real p_max = p_min;
        
        for (int i = 0; i < 100; ++i) {
            p_min = 0.5 * std::max(p_min, 0.0);
            p_max = 2.0 * p_max;
            
            real f_min = dvel(p_min);
            real f_max = dvel(p_max);
            
            if (f_min * f_max <= 0.0) {
                break;
            }
        }
        
        // Brent's method
        constexpr int max_iter = 100;
        constexpr real tol = 1.0e-10;
        
        real a = p_min, b = p_max;
        real fa = dvel(a), fb = dvel(b);
        
        if (std::abs(fa) < std::abs(fb)) {
            std::swap(a, b);
            std::swap(fa, fb);
        }
        
        real brent_c = a, fc = fa;
        bool mflag = true;
        real d = 0.0;
        
        for (int iter = 0; iter < max_iter; ++iter) {
            if (std::abs(fb) < tol || std::abs(b - a) < tol) {
                pstar = b;
                
                real nl_star, Hl_star, csl_star, vl_star;
                real nr_star, Hr_star, csr_star, vr_star;
                get_velocity_behind_wave(b, nl, Pl, Hl, csl, vl, 0.0, wl, -1.0,
                                         nl_star, Hl_star, csl_star, vl_star);
                get_velocity_behind_wave(b, nr, Pr, Hr, csr, vr, 0.0, wr, +1.0,
                                         nr_star, Hr_star, csr_star, vr_star);
                vstar = 0.5 * (vl_star + vr_star);
                iterations = iter + 1;
                return;
            }
            
            real s;
            if (fa != fc && fb != fc) {
                s = a * fb * fc / ((fa - fb) * (fa - fc))
                  + b * fa * fc / ((fb - fa) * (fb - fc))
                  + brent_c * fa * fb / ((fc - fa) * (fc - fb));
            } else {
                s = b - fb * (b - a) / (fb - fa);
            }
            
            const real bisect = 0.25 * (3.0 * a + b);
            bool use_bisect = false;
            
            if ((s < bisect && s < b) || (s > bisect && s > b)) use_bisect = true;
            if (mflag && std::abs(s - b) >= 0.5 * std::abs(b - brent_c)) use_bisect = true;
            if (!mflag && std::abs(s - b) >= 0.5 * std::abs(brent_c - d)) use_bisect = true;
            if (mflag && std::abs(b - brent_c) < tol) use_bisect = true;
            if (!mflag && std::abs(brent_c - d) < tol) use_bisect = true;
            
            if (use_bisect) {
                s = 0.5 * (a + b);
                mflag = true;
            } else {
                mflag = false;
            }
            
            const real fs = dvel(s);
            d = brent_c;
            brent_c = b;
            fc = fb;
            
            if (fa * fs < 0.0) {
                b = s;
                fb = fs;
            } else {
                a = s;
                fa = fs;
            }
            
            if (std::abs(fa) < std::abs(fb)) {
                std::swap(a, b);
                std::swap(fa, fb);
            }
        }
        
        pstar = b;
        real nl_star, Hl_star, csl_star, vl_star;
        real nr_star, Hr_star, csr_star, vr_star;
        get_velocity_behind_wave(b, nl, Pl, Hl, csl, vl, 0.0, wl, -1.0,
                                 nl_star, Hl_star, csl_star, vl_star);
        get_velocity_behind_wave(b, nr, Pr, Hr, csr, vr, 0.0, wr, +1.0,
                                 nr_star, Hr_star, csr_star, vr_star);
        vstar = 0.5 * (vl_star + vr_star);
        iterations = max_iter;
    }
};

int main() {
    std::cout << std::setprecision(16);
    
    std::cout << "================================================================================" << std::endl;
    std::cout << "C++ RIEMANN SOLVER: SR Sod Initial Conditions" << std::endl;
    std::cout << "================================================================================" << std::endl;
    
    SRRiemannSolver solver(5.0/3.0, 1.0);
    
    // Test 1: Sharp discontinuity
    std::cout << "\n1. SHARP DISCONTINUITY (Theoretical SR Sod):" << std::endl;
    std::cout << "--------------------------------------------------------------------------------" << std::endl;
    {
        real Pl = 1.0, nl = 1.0, vl = 0.0;
        real Pr = 0.1, nr = 0.125, vr = 0.0;
        real pstar, vstar;
        int iterations;
        
        solver.solve(Pl, nl, vl, Pr, nr, vr, pstar, vstar, iterations);
        
        std::cout << "  Input: Pl=" << Pl << ", nl=" << nl << ", vl=" << vl << std::endl;
        std::cout << "         Pr=" << Pr << ", nr=" << nr << ", vr=" << vr << std::endl;
        std::cout << "  P* = " << pstar << std::endl;
        std::cout << "  v* = " << vstar << std::endl;
        std::cout << "  Iterations: " << iterations << std::endl;
    }
    
    // Test 2: Smoothed states
    std::cout << "\n2. SMOOTHED STATES (SPH Initial Conditions):" << std::endl;
    std::cout << "--------------------------------------------------------------------------------" << std::endl;
    {
        real Pl = 0.73, nl = 0.7375, vl = 0.0;
        real Pr = 0.37, nr = 0.3875, vr = 0.0;
        real pstar, vstar;
        int iterations;
        
        solver.solve(Pl, nl, vl, Pr, nr, vr, pstar, vstar, iterations);
        
        std::cout << "  Input: Pl=" << Pl << ", nl=" << nl << ", vl=" << vl << std::endl;
        std::cout << "         Pr=" << Pr << ", nr=" << nr << ", vr=" << vr << std::endl;
        std::cout << "  P* = " << pstar << std::endl;
        std::cout << "  v* = " << vstar << std::endl;
        std::cout << "  Iterations: " << iterations << std::endl;
    }
    
    // Test 3: Moderate ratio
    std::cout << "\n3. MODERATE RATIO (Pressure ratio ~2x):" << std::endl;
    std::cout << "--------------------------------------------------------------------------------" << std::endl;
    {
        real Pl = 0.8, nl = 0.6, vl = 0.0;
        real Pr = 0.4, nr = 0.4, vr = 0.0;
        real pstar, vstar;
        int iterations;
        
        solver.solve(Pl, nl, vl, Pr, nr, vr, pstar, vstar, iterations);
        
        std::cout << "  Input: Pl=" << Pl << ", nl=" << nl << ", vl=" << vl << std::endl;
        std::cout << "         Pr=" << Pr << ", nr=" << nr << ", vr=" << vr << std::endl;
        std::cout << "  P* = " << pstar << std::endl;
        std::cout << "  v* = " << vstar << std::endl;
        std::cout << "  Iterations: " << iterations << std::endl;
    }
    
    return 0;
}
