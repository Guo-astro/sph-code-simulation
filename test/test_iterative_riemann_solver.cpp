/**
 * Test program for iterative relativistic Riemann solver
 * 
 * This standalone program tests the iterative solver from sr_fluid_force.cpp
 * with specific initial conditions and outputs detailed debug information
 * for comparison with the Python kitajima_solver.py implementation.
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <functional>
#include <vector>
#include <string>

// Type definitions matching SPH code
using real = double;

/**
 * Iterative relativistic Riemann solver with Newton-Raphson method
 * 
 * This is extracted from sr_fluid_force.cpp for standalone testing.
 * Extended for tangential velocities following Pons, Marti & Muller (1999).
 */
class IterativeRiemannSolver {
private:
    real m_gamma;  // Adiabatic index
    real m_c_speed; // Speed of light
    
    /**
     * Compute post-wave state using relativistic jump conditions
     * (Extracted from FluidForce::get_velocity_behind_wave)
     */
    void get_velocity_behind_wave(
        const real P,
        const real na, const real Pa, const real Ha, const real csa,
        const real va, const real vta, const real wa,
        const real sign,
        real &n, real &H, real &cs, real &vel) const
    {
        const real c = m_c_speed;
        const real gamma_c = m_gamma;

        if (P > Pa)
        {
            // SHOCK WAVE - Use Taub adiabat and Rankine-Hugoniot relations
            const real a = 1.0 + (gamma_c - 1.0) * (Pa - P) / (gamma_c * P);
            const real b_term = -(gamma_c - 1.0) * (Pa - P) / (gamma_c * P);
            const real c_term = Ha * (Pa - P) / na - Ha * Ha;

            // Discriminant check for physical solution
            const real discriminant = b_term * b_term - 4.0 * a * c_term;
            if (discriminant < 0.0)
            {
                n = na;
                H = Ha;
                cs = csa;
                vel = va;
                std::cout << "    WARNING: Unphysical shock state (negative discriminant)" << std::endl;
                return;
            }

            // Post-shock enthalpy
            H = (-b_term + std::sqrt(discriminant)) / (2.0 * a);

            if (H <= 0.0)
            {
                n = na;
                H = Ha;
                cs = csa;
                vel = va;
                std::cout << "    WARNING: Unphysical enthalpy (H <= 0)" << std::endl;
                return;
            }

            // Post-shock rest-frame density
            n = gamma_c * P / ((gamma_c - 1.0) * (H - 1.0));

            // Mass flux squared
            const real j2 = -(P - Pa) / (H / n - Ha / na);
            if (j2 < 0.0)
            {
                n = na;
                H = Ha;
                cs = csa;
                vel = va;
                std::cout << "    WARNING: Negative mass flux squared" << std::endl;
                return;
            }
            const real j = sign * std::sqrt(j2);

            // Shock velocity
            const real Na = wa * na;
            const real a_shock = j2 + Na * Na;
            const real b_shock = -va * Na * Na;
            const real sqrt_term = std::sqrt(1.0 + na * na / j2);
            const real vshock = (-b_shock + sign * j2 * sqrt_term) / a_shock;

            // Shock Lorentz factor
            const real gamma_shock = 1.0 / std::sqrt(1.0 - (vshock / c) * (vshock / c));

            // Post-shock normal velocity
            const real a_vel = gamma_shock * (P - Pa) / j + Ha * wa * va;
            const real b_vel = Ha * wa + (P - Pa) * (gamma_shock * va / j + 1.0 / (na * wa));
            vel = a_vel / b_vel;

            // SAFETY: Ensure post-shock velocity is subluminal
            const real v_limit = 0.99 * c;
            if (std::abs(vel) > v_limit)
            {
                vel = std::copysign(v_limit, vel);
            }

            // Post-shock sound speed
            cs = std::sqrt(gamma_c * P / (n * H));
        }
        else
        {
            // RAREFACTION WAVE - Use isentropic relations
            const real K = Pa / std::pow(na, gamma_c);
            n = std::pow(P / K, 1.0 / gamma_c);
            const real u = P / ((gamma_c - 1.0) * n);
            H = 1.0 + u / (c * c) + P / (n * c * c);
            cs = std::sqrt(gamma_c * P / (n * H));

            // Velocity from rarefaction invariant
            const real sqgl1 = std::sqrt(gamma_c - 1.0);
            const real term1 = (1.0 + va / c) / (1.0 - va / c);
            const real term2 = (sqgl1 + csa / c) / (sqgl1 - csa / c);
            const real term3 = (sqgl1 - cs / c) / (sqgl1 + cs / c);
            const real exponent = -sign * 2.0 / sqgl1;
            const real A = term1 * std::pow(term2 * term3, exponent);

            vel = c * (A - 1.0) / (A + 1.0);

            // SAFETY: Ensure post-rarefaction velocity is subluminal
            const real v_limit = 0.99 * c;
            if (std::abs(vel) > v_limit || !std::isfinite(vel))
            {
                vel = std::copysign(std::min(std::abs(va), v_limit), va);
            }
        }
    }

public:
    IterativeRiemannSolver(real gamma, real c_speed) 
        : m_gamma(gamma), m_c_speed(c_speed) {}
    
    /**
     * Solve Riemann problem for interface pressure and velocity
     * (Extracted from FluidForce::iterative_solver lambda)
     * 
     * Returns convergence info for debugging
     */
    struct SolverResult {
        real pstar;
        real vstar;
        int iterations;
        bool converged;
        real final_residual;
        std::vector<real> p_history;
        std::vector<real> f_history;
    };
    
    SolverResult solve(const real left[5], const real right[5])
    {
        SolverResult result;
        
        // Extract states
        const real vl = left[0];
        const real nl = left[1];
        const real Pl = left[2];
        const real csl = left[3];
        const real vtl = left[4];

        const real vr = right[0];
        const real nr = right[1];
        const real Pr = right[2];
        const real csr = right[3];
        const real vtr = right[4];

        const real c = m_c_speed;
        const real gamma_c = m_gamma;

        std::cout << "\nInput states:" << std::endl;
        std::cout << "  Left:  v=" << vl << ", n=" << nl << ", P=" << Pl 
                  << ", cs=" << csl << ", vt=" << vtl << std::endl;
        std::cout << "  Right: v=" << vr << ", n=" << nr << ", P=" << Pr 
                  << ", cs=" << csr << ", vt=" << vtr << std::endl;

        // Compute enthalpy H = 1 + u/c² + P/(n·c²)
        const real ul = Pl / ((gamma_c - 1.0) * nl);
        const real ur = Pr / ((gamma_c - 1.0) * nr);
        const real Hl = 1.0 + ul / (c * c) + Pl / (nl * c * c);
        const real Hr = 1.0 + ur / (c * c) + Pr / (nr * c * c);

        std::cout << "  Computed: Hl=" << Hl << ", Hr=" << Hr << std::endl;

        // Compute Lorentz factors including tangential velocities
        const real v2l = vl * vl + vtl * vtl;
        const real v2r = vr * vr + vtr * vtr;
        const real wl = 1.0 / std::sqrt(1.0 - v2l / (c * c));
        const real wr = 1.0 / std::sqrt(1.0 - v2r / (c * c));

        std::cout << "  Lorentz factors: wl=" << wl << ", wr=" << wr << std::endl;

        // Check for trivial solution
        constexpr real trivial_tol = 1.0e-12;
        if (std::abs(Pl - Pr) < trivial_tol * std::max(Pl, Pr) &&
            std::abs(vl - vr) < trivial_tol * c)
        {
            std::cout << "  Trivial solution detected" << std::endl;
            result.pstar = Pl;
            result.vstar = vl;
            result.iterations = 0;
            result.converged = true;
            result.final_residual = 0.0;
            return result;
        }

        // ============================================================
        // BRENT'S METHOD with ADAPTIVE BRACKETING (matches Python)
        // ============================================================
        
        // Lambda to compute velocity difference dvel = vl_star - vr_star
        auto dvel = [&](real p_test) -> real {
            real nl_star, Hl_star, csl_star, vl_star;
            real nr_star, Hr_star, csr_star, vr_star;
            
            get_velocity_behind_wave(p_test, nl, Pl, Hl, csl, vl, vtl, wl, -1.0,
                                     nl_star, Hl_star, csl_star, vl_star);
            get_velocity_behind_wave(p_test, nr, Pr, Hr, csr, vr, vtr, wr, +1.0,
                                     nr_star, Hr_star, csr_star, vr_star);
            return vl_star - vr_star;
        };
        
        // Adaptive bracketing: find [p_min, p_max] where dvel changes sign
        real p_min = 0.5 * (Pl + Pr);
        real p_max = p_min;
        
        constexpr int max_bracket = 100;
        bool bracket_found = false;
        real f_min = 0.0, f_max = 0.0;
        
        for (int i = 0; i < max_bracket; ++i) {
            p_min = 0.5 * std::max(p_min, 0.0);
            p_max = 2.0 * p_max;
            
            f_min = dvel(p_min);
            f_max = dvel(p_max);
            
            if (f_min * f_max <= 0.0) {
                bracket_found = true;
                break;
            }
        }
        
        if (!bracket_found) {
            std::cout << "  ERROR: Bracketing failed after " << max_bracket << " iterations" << std::endl;
            result.converged = false;
            result.pstar = 0.5 * (Pl + Pr);
            result.vstar = 0.5 * (vl + vr);
            result.iterations = max_bracket;
            result.final_residual = std::min(std::abs(f_min), std::abs(f_max));
            return result;
        }
        
        // Brent's method (Brent 1973)
        constexpr int max_iter = 100;
        constexpr real tol = 1.0e-10;
        
        real a = p_min, b = p_max;
        real fa = f_min, fb = f_max;
        
        if (std::abs(fa) < std::abs(fb)) {
            std::swap(a, b);
            std::swap(fa, fb);
        }
        
        real brent_c = a, fc = fa;
        bool mflag = true;
        real d = 0.0;
        
        for (int iter = 0; iter < max_iter; ++iter) {
            result.p_history.push_back(b);
            result.f_history.push_back(fb);
            
            // Convergence check
            if (std::abs(fb) < tol || std::abs(b - a) < tol) {
                result.pstar = b;
                
                // Compute final velocity
                real nl_star, Hl_star, csl_star, vl_star;
                real nr_star, Hr_star, csr_star, vr_star;
                get_velocity_behind_wave(b, nl, Pl, Hl, csl, vl, vtl, wl, -1.0,
                                         nl_star, Hl_star, csl_star, vl_star);
                get_velocity_behind_wave(b, nr, Pr, Hr, csr, vr, vtr, wr, +1.0,
                                         nr_star, Hr_star, csr_star, vr_star);
                result.vstar = 0.5 * (vl_star + vr_star);
                result.iterations = iter + 1;
                result.converged = true;
                result.final_residual = fb;
                return result;
            }
            
            // Inverse quadratic interpolation or secant
            real s;
            if (fa != fc && fb != fc) {
                // Inverse quadratic interpolation
                s = a * fb * fc / ((fa - fb) * (fa - fc))
                  + b * fa * fc / ((fb - fa) * (fb - fc))
                  + brent_c * fa * fb / ((fc - fa) * (fc - fb));
            } else {
                // Secant method
                s = b - fb * (b - a) / (fb - fa);
            }
            
            // Bisection fallback conditions
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
        
        // Convergence failure
        std::cout << "  WARNING: Brent's method did not converge" << std::endl;
        result.pstar = b;
        
        real nl_star, Hl_star, csl_star, vl_star;
        real nr_star, Hr_star, csr_star, vr_star;
        get_velocity_behind_wave(b, nl, Pl, Hl, csl, vl, vtl, wl, -1.0,
                                 nl_star, Hl_star, csl_star, vl_star);
        get_velocity_behind_wave(b, nr, Pr, Hr, csr, vr, vtr, wr, +1.0,
                                 nr_star, Hr_star, csr_star, vr_star);
        result.vstar = 0.5 * (vl_star + vr_star);
        result.iterations = max_iter;
        result.converged = false;
        result.final_residual = fb;
        
        return result;
    }
};

/**
 * Test case structure
 */
struct TestCase {
    std::string name;
    real gamma;
    real c;
    real left[5];   // [v^n, n, P, cs, |v^t|]
    real right[5];
};

int main(int argc, char* argv[])
{
    std::cout << std::setprecision(16);
    
    // Define test cases
    std::vector<TestCase> test_cases;
    
    // Test 1: SR Sod shock tube (from sr_sod.json)
    test_cases.push_back({
        "SR_Sod_Shock_Tube",
        5.0/3.0,  // gamma
        1.0,      // c (natural units)
        {0.0, 10.0, 13.33, 0.0, 0.0},  // left: v, n, P, cs (computed), vt
        {0.0, 1.0, 1.0e-8, 0.0, 0.0}   // right
    });
    
    // Test 2: Two rarefactions
    test_cases.push_back({
        "Two_Rarefactions",
        5.0/3.0,
        1.0,
        {-0.5, 1.0, 1.0, 0.0, 0.0},
        {0.5, 1.0, 1.0, 0.0, 0.0}
    });
    
    // Test 3: Two shocks
    test_cases.push_back({
        "Two_Shocks",
        5.0/3.0,
        1.0,
        {0.5, 1.0, 1.0, 0.0, 0.0},
        {-0.5, 1.0, 1.0, 0.0, 0.0}
    });
    
    // Test 4: Left shock, right rarefaction
    test_cases.push_back({
        "Left_Shock_Right_Rarefaction",
        5.0/3.0,
        1.0,
        {0.3, 1.0, 10.0, 0.0, 0.0},
        {0.0, 1.0, 1.0, 0.0, 0.0}
    });
    
    // Test 5: With tangential velocity (Pons et al. test)
    test_cases.push_back({
        "With_Tangential_Velocity",
        5.0/3.0,
        1.0,
        {0.0, 1.0, 1.0, 0.0, 0.3},  // vt = 0.3
        {0.0, 1.0, 0.1, 0.0, 0.0}
    });
    
    // Compute sound speeds for test cases
    for (auto& tc : test_cases) {
        // Left sound speed
        real ul = tc.left[2] / ((tc.gamma - 1.0) * tc.left[1]);
        real Hl = 1.0 + ul/(tc.c*tc.c) + tc.left[2]/(tc.left[1]*tc.c*tc.c);
        tc.left[3] = std::sqrt(tc.gamma * tc.left[2] / (tc.left[1] * Hl));
        
        // Right sound speed
        real ur = tc.right[2] / ((tc.gamma - 1.0) * tc.right[1]);
        real Hr = 1.0 + ur/(tc.c*tc.c) + tc.right[2]/(tc.right[1]*tc.c*tc.c);
        tc.right[3] = std::sqrt(tc.gamma * tc.right[2] / (tc.right[1] * Hr));
    }
    
    // Run all test cases
    std::ofstream outfile("test_results_cpp.txt");
    outfile << std::setprecision(16);
    
    for (const auto& tc : test_cases) {
        std::cout << "\n" << std::string(80, '=') << std::endl;
        std::cout << "TEST CASE: " << tc.name << std::endl;
        std::cout << std::string(80, '=') << std::endl;
        
        outfile << "\n" << std::string(80, '=') << std::endl;
        outfile << "TEST: " << tc.name << std::endl;
        outfile << "gamma=" << tc.gamma << ", c=" << tc.c << std::endl;
        outfile << "LEFT:  v=" << tc.left[0] << ", n=" << tc.left[1] 
                << ", P=" << tc.left[2] << ", cs=" << tc.left[3] 
                << ", vt=" << tc.left[4] << std::endl;
        outfile << "RIGHT: v=" << tc.right[0] << ", n=" << tc.right[1] 
                << ", P=" << tc.right[2] << ", cs=" << tc.right[3] 
                << ", vt=" << tc.right[4] << std::endl;
        
        IterativeRiemannSolver solver(tc.gamma, tc.c);
        auto result = solver.solve(tc.left, tc.right);
        
        outfile << "\nRESULT:" << std::endl;
        outfile << "pstar=" << result.pstar << std::endl;
        outfile << "vstar=" << result.vstar << std::endl;
        outfile << "iterations=" << result.iterations << std::endl;
        outfile << "converged=" << (result.converged ? "true" : "false") << std::endl;
        outfile << "final_residual=" << result.final_residual << std::endl;
        
        outfile << "\nITERATION_HISTORY:" << std::endl;
        for (size_t i = 0; i < result.p_history.size(); ++i) {
            outfile << i << " " << result.p_history[i] << " " << result.f_history[i] << std::endl;
        }
        outfile << std::string(80, '=') << std::endl;
        
        std::cout << "\nFINAL RESULT:" << std::endl;
        std::cout << "  pstar = " << result.pstar << std::endl;
        std::cout << "  vstar = " << result.vstar << std::endl;
        std::cout << "  Converged: " << (result.converged ? "YES" : "NO") 
                  << " in " << result.iterations << " iterations" << std::endl;
    }
    
    outfile.close();
    std::cout << "\n\nResults written to test_results_cpp.txt" << std::endl;
    
    return 0;
}
