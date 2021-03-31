#ifndef Foptim_ALagrangian_Newton_Raphson_hpp
#define Foptim_ALagrangian_Newton_Raphson_hpp

#include <cstdint>

namespace { extern "C" {
    void alagrangian_newton_raphson_(
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
        const int32_t & max_StepIteration,
        const double & precision, const double & min_StepLength
    );
} }

namespace Foptim {

inline void ALagrangian_Newton_Raphson(
    void (*f)(double &, const double *, const int32_t &),
    void (*f_fd)(double &, double *, const double *, const int32_t &),
    void (*fdd)(double *, const double *, const int32_t &),
    void (*c)(double *, const double *, const int32_t &, const int32_t &),
    void (*c_cd)(double *, double *, const double *, const int32_t &, const int32_t &),
    void (*c_cd_cdd)(double *, double *, double *, const double *, const int32_t &, const int32_t &),
    double * x, const int32_t & N, const int32_t & M,
    const double * lambda0 = nullptr, const double & miu0 = 1.0,
    const int32_t & max_iteration = 100,
    const int32_t & max_StepIteration = 100,
    const double & precision = 1e-12, const double & min_StepLength = 1e-12
) {
    bool empty_lambda0 = nullptr == lambda0;
    if (empty_lambda0) lambda0 = new double[M]();
    alagrangian_newton_raphson_(
        f, f_fd, fdd, c, c_cd, c_cd_cdd, x, N, M,
        lambda0, miu0,
        max_iteration,
        max_StepIteration, precision, min_StepLength
    );
    if (empty_lambda0) delete [] lambda0;
}

} // namespace optim

#endif