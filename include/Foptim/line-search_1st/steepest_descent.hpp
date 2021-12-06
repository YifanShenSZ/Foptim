#ifndef Foptim_line_search_1st_steepest_descent_hpp
#define Foptim_line_search_1st_steepest_descent_hpp

#include <cstdint>

namespace { extern "C" {
    void steepest_descent_(
        // Required arguments
        void (*f)(double &, const double *, const int32_t &),
        void (*f_fd)(double &, double *, const double *, const int32_t &),
        double * x, const int32_t & dim,
        // Optional arguments
        const double & initial_learning_rate,
        const int32_t & max_iteration,
        const double & precision, const double & min_StepLength
    );
} }

namespace Foptim {

inline void steepest_descent(
    void (*f)(double &, const double *, const int32_t &),
    void (*f_fd)(double &, double *, const double *, const int32_t &),
    double * x, const int32_t & dim,
    const double & initial_learning_rate = 1e-3,
    const int32_t & max_iteration = 100,
    const double & precision = 1e-12, const double & min_StepLength = 1e-12
) {
    steepest_descent_(
        f, f_fd, x, dim,
        initial_learning_rate, max_iteration, precision, min_StepLength
    );
}

} // namespace optim

#endif