# Foptim
FORTRAN nonlinear optimization library with c++ frontend

As [Fortran-Library](https://github.com/YifanShenSZ/Fortran-Library) supports more and more functionalities, it becomes too heavy and clumsy to add and debug c++ features. This library extracts only the nonlinear optimization routines and provides c++ frontend

This library aims at serving c++, so compared to [Fortran-Library](https://github.com/YifanShenSZ/Fortran-Library) the FORTRAN-only tricks are removed, e.g. the presence of optional arguments

See `trust_region.md` and `line_search.md` for details of available nonlinear optimization solvers

## Reference
> 1. J. Nocedal, S. J. Wright, *Numerical Optimization 2nd edition* (Springer, 2006)