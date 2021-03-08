#ifndef Foptim_Gauss_BFGS_hpp
#define Foptim_Gauss_BFGS_hpp

#include <cstdint>

namespace { extern "C" {
    void gauss_bfgs_(
        // Required arguments
        void (*residue)(double *, const double *, const int32_t &, const int32_t &),
        void (*Jacobian)(double *, const double *, const int32_t &, const int32_t &),
        double * x, const int32_t & M, const int32_t & N,
        // Optional arguments
        const int32_t & max_iteration,
        const double & precision, const double & min_StepLength
    );
} }

namespace Foptim {

inline void Gauss_BFGS(
    void (*residue)(double *, const double *, const int32_t &, const int32_t &),
    void (*Jacobian)(double *, const double *, const int32_t &, const int32_t &),
    double * x, const int32_t & M, const int32_t & N,
    const int32_t & max_iteration = 100,
    const double & precision = 1e-15, const double & min_StepLength = 1e-15
) {
    gauss_bfgs_(
        residue, Jacobian, x, M, N,
        max_iteration, precision, min_StepLength
    );
}

} // namespace optim

#endif