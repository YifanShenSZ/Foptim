#ifndef Foptim_augmented_Lagrangian_hpp
#define Foptim_augmented_Lagrangian_hpp

#include <cstdint>

namespace { extern "C" {
    void augmented_lagrangian_(
        // Required arguments
        void (*f)(double &, const double *, const int32_t &),
        void (*f_fd)(double &, double *, const double *, const int32_t &),
        void (*fdd)(double *, const double *, const int32_t &),
        void (*c)(double *, const double *, const int32_t &, const int32_t &),
        void (*c_cd)(double *, double *, const double *, const int32_t &, const int32_t &),
        void (*c_cd_cdd)(double *, double *, double *, const double *, const int32_t &, const int32_t &),
        double * x, const int32_t & N, const int32_t & M,
        // Optional arguments
        const double * lambda0, const double & miu0,
        const int32_t & max_iteration,
        const int32_t & Hessian_step, const int32_t & max_StepIteration,
        const double & precision, const double & min_StepLength
    );
} }

namespace Foptim {

inline void augmented_Lagrangian(
    void (*f)(double &, const double *, const int32_t &),
    void (*f_fd)(double &, double *, const double *, const int32_t &),
    void (*fdd)(double *, const double *, const int32_t &),
    void (*c)(double *, const double *, const int32_t &, const int32_t &),
    void (*c_cd)(double *, double *, const double *, const int32_t &, const int32_t &),
    void (*c_cd_cdd)(double *, double *, double *, const double *, const int32_t &, const int32_t &),
    double * x, const int32_t & N, const int32_t & M,
    const double * lambda0 = nullptr, const double & miu0 = 1.0,
    const int32_t & max_iteration = 100,
    const int32_t & Hessian_step = 10, const int32_t & max_StepIteration = 100,
    const double & precision = 1e-15, const double & min_StepLength = 1e-15
) {
    bool empty_lambda0 = nullptr == lambda0;
    if (empty_lambda0) lambda0 = new double[M]();
    augmented_lagrangian_(
        f, f_fd, fdd, c, c_cd, c_cd_cdd, x, N, M,
        lambda0, miu0,
        max_iteration,
        Hessian_step, max_StepIteration, precision, min_StepLength
    );
    if (empty_lambda0) delete [] lambda0;
}

} // namespace optim

#endif