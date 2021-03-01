#ifndef Foptim_trust_region_hpp
#define Foptim_trust_region_hpp

#include <cstdint>

namespace { extern "C" {
    void trust_region_(
        // Required arguments
        void (*residue)(double *, const double *, const int32_t &, const int32_t &),
        void (*Jacobian)(double *, const double *, const int32_t &, const int32_t &),
        double * x, const int32_t & M, const int32_t & N,
        // Optional arguments
        const int32_t & max_iteration, const int32_t & max_StepIteration,
        const double & precision, const double & min_StepLength
    );
} }

namespace Foptim {

inline void trust_region(
    void (*residue)(double *, const double *, const int32_t &, const int32_t &),
    void (*Jacobian)(double *, const double *, const int32_t &, const int32_t &),
    double * x, const int32_t & M, const int32_t & N,
    const int32_t & max_iteration = 100, const int32_t & max_StepIteration = 100,
    const double & precision = 1e-15, const double & min_StepLength = 1e-15
) {
    trust_region_(
        residue, Jacobian, x, M, N,
        max_iteration, max_StepIteration, precision, min_StepLength
    );
}

} // namespace optim

#endif