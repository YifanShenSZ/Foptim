#include <iostream>
#include <cmath>
#include <time.h>

#include <Foptim/Foptim.hpp>

void f(double & fx, const double * x, const int32_t & dim) {
    fx = 0.0;
    for (int32_t i = 0; i < dim; i++) fx += x[i] * x[i] * x[i] * x[i];
}

void f_fd(double & fx, double * fdx, const double * x, const int32_t & dim) {
    fx = 0.0;
    for (int32_t i = 0; i < dim; i++) {
        fx += x[i] * x[i] * x[i] * x[i];
        fdx[i] = 4.0 * x[i] * x[i] * x[i];
    }
}

void fdd(double * fddx, const double * x, const int32_t & dim) {
    for (int32_t i = 0; i < dim; i++)
    for (int32_t j = 0; j < dim; j++) {
        int32_t location = i * dim + j;
        if (i == j) fddx[location] = 12.0 * x[i] * x[i];
        else        fddx[location] = 0.0;
    }
}

void fd_residue(double * fdx, const double * x, const int32_t & M, const int32_t & N) {
    for (int32_t i = 0; i < M; i++) fdx[i] = 4.0 * x[i] * x[i] * x[i];
}

void fdd_Jacobian(double * fddx, const double * x, const int32_t & M, const int32_t & N) {
    for (int32_t i = 0; i < M; i++)
    for (int32_t j = 0; j < N; j++) {
        int32_t location = i * M + j;
        if (i == j) fddx[location] = 12.0 * x[i] * x[i];
        else        fddx[location] = 0.0;
    }
}

void c(double * cx, const double * x, const int32_t & M, const int32_t & N) {
    cx[0] = -1.0;
    for (int32_t i = 0; i < N; i++) cx[0] += x[i] * x[i];
}

void c_cd(double * cx, double * cdx, const double * x, const int32_t & M, const int32_t & N) {
    cx[0] = -1.0;
    for (int32_t i = 0; i < N; i++) {
        cx[0] += x[i] * x[i];
        cdx[i] = 2.0 * x[i];
    }
}

void c_cd_cdd(double * cx, double * cdx, double * cddx, const double * x, const int32_t & M, const int32_t & N) {
    cx[0] = -1.0;
    for (int32_t i = 0; i < N; i++) {
        cx[0] += x[i] * x[i];
        cdx[i] = 2.0 * x[i];
        for (int32_t j = 0; j < N; j++) {
            int32_t location = i * N + j;
            if (i == j) cddx[location] = 2.0;
            else        cddx[location] = 0.0;
        }
    }
}

int main() {
    srand(time(NULL));

    int32_t dim = 10;
    double * x = new double[dim];
    double norm;

    std::cout << "Steepest descent" << std::endl;
    for (int32_t i = 0; i < dim; i++) x[i] = (double)rand() / (double)RAND_MAX;
    Foptim::steepest_descent(f, f_fd, x, dim);
    norm = 0.0;
    for (int32_t i = 0; i < dim; i++) norm += x[i] * x[i];
    std::cout << sqrt(norm) << '\n' << std::endl;

    std::cout << "Steepest descent verbose" << std::endl;
    for (int32_t i = 0; i < dim; i++) x[i] = (double)rand() / (double)RAND_MAX;
    Foptim::steepest_descent_verbose(f, f_fd, x, dim);
    norm = 0.0;
    for (int32_t i = 0; i < dim; i++) norm += x[i] * x[i];
    std::cout << sqrt(norm) << '\n' << std::endl;

    std::cout << "Dai-Yuan conjugate gradient" << std::endl;
    for (int32_t i = 0; i < dim; i++) x[i] = (double)rand() / (double)RAND_MAX;
    Foptim::CGDY(f, f_fd, x, dim);
    norm = 0.0;
    for (int32_t i = 0; i < dim; i++) norm += x[i] * x[i];
    std::cout << sqrt(norm) << '\n' << std::endl;

    std::cout << "Dai-Yuan conjugate gradient" << std::endl;
    for (int32_t i = 0; i < dim; i++) x[i] = (double)rand() / (double)RAND_MAX;
    Foptim::CGDY_verbose(f, f_fd, x, dim);
    norm = 0.0;
    for (int32_t i = 0; i < dim; i++) norm += x[i] * x[i];
    std::cout << sqrt(norm) << '\n' << std::endl;

    std::cout << "Polak-Ribiere+ conjugate gradient" << std::endl;
    for (int32_t i = 0; i < dim; i++) x[i] = (double)rand() / (double)RAND_MAX;
    Foptim::CGPR(f, f_fd, x, dim);
    norm = 0.0;
    for (int32_t i = 0; i < dim; i++) norm += x[i] * x[i];
    std::cout << sqrt(norm) << '\n' << std::endl;

    std::cout << "Polak-Ribiere+ conjugate gradient" << std::endl;
    for (int32_t i = 0; i < dim; i++) x[i] = (double)rand() / (double)RAND_MAX;
    Foptim::CGPR_verbose(f, f_fd, x, dim);
    norm = 0.0;
    for (int32_t i = 0; i < dim; i++) norm += x[i] * x[i];
    std::cout << sqrt(norm) << '\n' << std::endl;

    std::cout << "Newton-Raphson" << std::endl;
    for (int32_t i = 0; i < dim; i++) x[i] = (double)rand() / (double)RAND_MAX;
    Foptim::NewtonRaphson(f, f_fd, fdd, x, dim);
    norm = 0.0;
    for (int32_t i = 0; i < dim; i++) norm += x[i] * x[i];
    std::cout << sqrt(norm) << '\n' << std::endl;

    std::cout << "BFGS" << std::endl;
    for (int32_t i = 0; i < dim; i++) x[i] = (double)rand() / (double)RAND_MAX;
    Foptim::BFGS(f, f_fd, fdd, x, dim);
    norm = 0.0;
    for (int32_t i = 0; i < dim; i++) norm += x[i] * x[i];
    std::cout << sqrt(norm) << '\n' << std::endl;

    std::cout << "LBFGS" << std::endl;
    for (int32_t i = 0; i < dim; i++) x[i] = (double)rand() / (double)RAND_MAX;
    Foptim::LBFGS(f, f_fd, x, dim);
    norm = 0.0;
    for (int32_t i = 0; i < dim; i++) norm += x[i] * x[i];
    std::cout << sqrt(norm) << '\n' << std::endl;

    std::cout << "LBFGS verbose" << std::endl;
    for (int32_t i = 0; i < dim; i++) x[i] = (double)rand() / (double)RAND_MAX;
    Foptim::LBFGS_verbose(f, f_fd, x, dim);
    norm = 0.0;
    for (int32_t i = 0; i < dim; i++) norm += x[i] * x[i];
    std::cout << sqrt(norm) << '\n' << std::endl;

    std::cout << "Trust region" << std::endl;
    for (int32_t i = 0; i < dim; i++) x[i] = (double)rand() / (double)RAND_MAX;
    Foptim::trust_region(fd_residue, fdd_Jacobian, x, dim, dim);
    norm = 0.0;
    for (int32_t i = 0; i < dim; i++) norm += x[i] * x[i];
    std::cout << sqrt(norm) << '\n' << std::endl;

    std::cout << "Trust region verbose" << std::endl;
    for (int32_t i = 0; i < dim; i++) x[i] = (double)rand() / (double)RAND_MAX;
    Foptim::trust_region_verbose(fd_residue, fdd_Jacobian, x, dim, dim);
    norm = 0.0;
    for (int32_t i = 0; i < dim; i++) norm += x[i] * x[i];
    std::cout << sqrt(norm) << '\n' << std::endl;

    std::cout << "Gauss-BFGS" << std::endl;
    for (int32_t i = 0; i < dim; i++) x[i] = (double)rand() / (double)RAND_MAX;
    Foptim::Gauss_BFGS(fd_residue, fdd_Jacobian, x, dim, dim);
    norm = 0.0;
    for (int32_t i = 0; i < dim; i++) norm += x[i] * x[i];
    std::cout << sqrt(norm) << '\n' << std::endl;

    std::cout << "Augmented Lagrangian based on Newton-Raphson" << std::endl;
    for (int32_t i = 0; i < dim; i++) x[i] = (double)rand() / (double)RAND_MAX;
    Foptim::ALagrangian_NewtonRaphson(f, f_fd, fdd, c, c_cd, c_cd_cdd, x, 10, 1);
    norm = 0.0;
    for (int32_t i = 0; i < dim; i++) norm += x[i] * x[i];
    std::cout << sqrt(norm) - 1.0 << '\n' << std::endl;

    std::cout << "Augmented Lagrangian based on BFGS" << std::endl;
    for (int32_t i = 0; i < dim; i++) x[i] = (double)rand() / (double)RAND_MAX;
    Foptim::ALagrangian_BFGS(f, f_fd, fdd, c, c_cd, c_cd_cdd, x, 10, 1);
    norm = 0.0;
    for (int32_t i = 0; i < dim; i++) norm += x[i] * x[i];
    std::cout << sqrt(norm) - 1.0 << '\n' << std::endl;
}