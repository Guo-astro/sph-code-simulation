#pragma once

#include <cmath>
#include "defines.hpp"
#include "kernel_function.hpp"

// Gaussian kernel (as in paper Eq. 11)
// W(x, h) = [1/(h√π)]^d × exp(-x²/h²)
namespace sph
{
namespace Gaussian
{
#if DIM == 1
    constexpr real sigma_gauss = 0.564189583547756;  // 1/√π
#elif DIM == 2
    constexpr real sigma_gauss = 1.0 / M_PI;       // 1/π
#else
    constexpr real sigma_gauss = 0.179587122125167;  // 1/(π^(3/2))
#endif

class Gauss : public KernelFunction {
public:
    Gauss() {}

    real w(const real r, const real h) const
    {
        const real q = r / h;
        const real q2 = q * q;
        return sigma_gauss / powh(h) * std::exp(-q2);
    }

    vec_t dw(const vec_t &rij, const real r, const real h) const
    {
        if(r == 0.0) {
            return vec_t(0);
        }
        const real q = r / h;
        const real q2 = q * q;
        // dW/dr = W(r) × (-2r/h²) = W(r) × (-2q/h)
        // ∇W = (r_ij/r) × (dW/dr) = (r_ij/r) × W(r) × (-2q/h)
        const real w_val = sigma_gauss / powh(h) * std::exp(-q2);
        const real c = -2.0 * w_val * q / (h * r);  // Divide by r to get unit vector when multiplied by r_ij
        return rij * c;
    }

    real dhw(const real r, const real h) const
    {
        const real q = r / h;
        const real q2 = q * q;
        // ∂W/∂h = W × [(-d/h) + (2r²/h³)] = W × [(-d + 2q²)/h]
        const real w_val = sigma_gauss / powh(h) * std::exp(-q2);
        return w_val * (2.0 * q2 - static_cast<real>(DIM)) / h;
    }
};

}
}
