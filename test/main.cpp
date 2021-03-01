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

void fd_tr(double * fdx, const double * x, const int32_t & M, const int32_t & N) {
    for (int32_t i = 0; i < M; i++) fdx[i] = 4.0 * x[i] * x[i] * x[i];
}

void fdd_tr(double * fddx, const double * x, const int32_t & M, const int32_t & N) {
    for (int32_t i = 0; i < M; i++)
    for (int32_t j = 0; j < N; j++) {
        int32_t location = i * M + j;
        if (i == j) fddx[location] = 12.0 * x[i] * x[i];
        else        fddx[location] = 0.0;
    }
}

int main() {
    srand(time(NULL));

    int32_t dim = 10;
    double * x = new double[dim];
    double norm;

    std::cout << "BFGS" << std::endl;
    for (int32_t i = 0; i < dim; i++) x[i] = (double)rand() / (double)RAND_MAX; 
    Foptim::BFGS(f, f_fd, fdd, x, dim);
    norm = 0.0;
    for (int32_t i = 0; i < dim; i++) norm += x[i] * x[i];
    std::cout << sqrt(norm) << '\n' << std::endl;

    std::cout << "Trust region" << std::endl;
    for (int32_t i = 0; i < dim; i++) x[i] = (double)rand() / (double)RAND_MAX; 
    Foptim::trust_region(fd_tr, fdd_tr, x, dim, dim);
    norm = 0.0;
    for (int32_t i = 0; i < dim; i++) norm += x[i] * x[i];
    std::cout << sqrt(norm) << '\n' << std::endl;
}