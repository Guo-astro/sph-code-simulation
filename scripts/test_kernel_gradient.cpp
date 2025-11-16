// Test cubic spline kernel gradient
#include <iostream>
#include <cmath>

constexpr int DIM = 1;
using real = double;

inline real powh(const real h) {
#if DIM == 1
    return h;
#elif DIM == 2
    return h * h;
#elif DIM == 3
    return h * h * h;
#endif
}

inline real sqr(real x) { return x*x; }
inline real pow3(real x) { return x*x*x; }

constexpr real sigma_cubic = 2.0 / 3.0;

real dw_magnitude(const real r, const real h) {
    if(r == 0.0) return 0.0;
    
    const real h_ = h * 0.5;
    const real q = r / h_;
    
    // Cubic spline gradient formula
    real term1 = 0.5 * (2.0 - q + std::abs(2.0 - q));  // This is (2-q) if q<2, else 0
    real term2 = 0.5 * (1.0 - q + std::abs(1.0 - q));  // This is (1-q) if q<1, else 0
    
    real dW_dq = 0.75 * term1 * term1 - 3.0 * term2 * term2;
    
    // Gradient magnitude: |dW/dr| = |dW/dq| / h_
    const real dW_dr = std::abs(dW_dq) / h_;
    
    // But the code multiplies by extra factor!
    const real c = -sigma_cubic / (powh(h_) * h_ * r) * dW_dq;
    
    printf("  Debug: q=%.4f, term1=%.4f, term2=%.4f, dW_dq=%.4f, dW_dr=%.2e, code_c=%.2e\n",
           q, term1, term2, dW_dq, dW_dr, std::abs(c));
    
    return std::abs(c);
}

int main() {
    // Test with actual h values from simulation
    // From simulation: h_i=0.00397758, passing 2*h_i to kernel
    const real dx = 0.00015625;
    const real h_i = 0.00397758;
    const real h_arg = 2.0 * h_i;  // Argument passed to kernel (2h)
    
    std::cout << "Testing cubic spline gradient\n";
    std::cout << "h_i = " << h_i << " (smoothing length)\n";
    std::cout << "h_arg = 2*h_i = " << h_arg << " (passed to kernel)\n";
    std::cout << "dx = " << dx << "\n\n";
    std::cout << "Gradient magnitude vs distance:\n";
    std::cout << "j\tr\t\t|dW/dr|\t\tq=r/h_eff\n";
    std::cout << "--\t----\t\t-------\t\t---------\n";
    
    for(int i = 1; i <= 10; i++) {
        real r = i * dx;
        real dw_mag = dw_magnitude(r, h_arg);
        real h_eff = h_arg * 0.5;  // After h_=h*0.5
        real q = r / h_eff;
        
        printf("%d\t%.6e\t%.6e\t%.4f\n", i, r, dw_mag, q);
    }
    
    std::cout << "\nExpected from simulation:" << std::endl;
    std::cout << "j=1 (r=0.00015625): |dw|=1522.53" << std::endl;
    std::cout << "j=2 (r=0.0003125): |dw|=6472.62" << std::endl;
    std::cout << "j=3 (r=0.00046875): |dw|=15384.8" << std::endl;
    
    return 0;
}
