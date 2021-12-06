#ifndef Foptim_line_search_2nd_NewtonRaphson_hpp
#define Foptim_line_search_2nd_NewtonRaphson_hpp

#include <cstdint>

namespace { extern "C" {
    void newtonraphson_(
        // Required arguments
        void (*f)(double &, const double *, const int32_t &),
        void (*f_fd)(double &, double *, const double *, const int32_t &),
        void (*fdd)(double *, const double *, const int32_t &),
        double * x, const int32_t & dim,
        // Optional arguments
        const int32_t & max_iteration,
        const double & precision, const double & min_StepLength
    );
} }

namespace Foptim {

inline void NewtonRaphson(
    void (*f)(double &, const double *, const int32_t &),
    void (*f_fd)(double &, double *, const double *, const int32_t &),
    void (*fdd)(double *, const double *, const int32_t &),
    double * x, const int32_t & dim,
    const int32_t & max_iteration = 100,
    const double & precision = 1e-12, const double & min_StepLength = 1e-12
) {
    newtonraphson_(
        f, f_fd, fdd, x, dim,
        max_iteration, precision, min_StepLength
    );
}

} // namespace optim

#endif