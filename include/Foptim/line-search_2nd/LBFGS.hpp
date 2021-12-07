#ifndef Foptim_line_search_2nd_LBFGS_hpp
#define Foptim_line_search_2nd_LBFGS_hpp

#include <cstdint>

namespace { extern "C" {
    void lbfgs_(
        // Required arguments
        void (*f)(double &, const double *, const int32_t &),
        void (*f_fd)(double &, double *, const double *, const int32_t &),
        double * x, const int32_t & dim,
        // Optional arguments
        const double & initial_learning_rate,
        const int32_t & memory, const int32_t & max_iteration,
        const double & precision, const double & min_StepLength
    );
} }

namespace Foptim {

inline void LBFGS(
    void (*f)(double &, const double *, const int32_t &),
    void (*f_fd)(double &, double *, const double *, const int32_t &),
    double * x, const int32_t & dim,
    const double & initial_learning_rate = 1e-3,
    const int32_t & memory = 10, const int32_t & max_iteration = 100,
    const double & precision = 1e-12, const double & min_StepLength = 1e-12
) {
    lbfgs_(
        f, f_fd, x, dim,
        initial_learning_rate, memory, max_iteration, precision, min_StepLength
    );
}

} // namespace optim

#endif