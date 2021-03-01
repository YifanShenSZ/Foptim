# Line search solver
Available solvers:
* steepest descent
* conjugate gradient
* limited-memory Broyden–Fletcher–Goldfarb–Shanno
* Broyden–Fletcher–Goldfarb–Shanno
* Newton-Raphson

Suggestion:
* If dimensionality is low, adopt quasi-Newton (or Newton if Hessian is cheap and initial guess is close)
* If dimensionality is so high that O(dim^2) memory is unaffordable, adopt conjugate gradient or L-BFGS

Nomenclature:
* f = the target function to be minimized
* a = the line search step length
* p = the line search direction
* phi(a) = f( x + a * p ), so phi'(a) = f'( x + a * p ) . p

External procedure:
* subroutine f(f(x), x, dim)
* subroutine fd(f'(x), x, dim)
* subroutine f_fd(f(x), f'(x), x, dim)
* subroutine fdd(f''(x), x, dim)
* dim dimensional vector x & f'(x), dim order matrix f''(x)
    
Required argument:
* subroutine f & fd
* dim dimensional vector x
* integer dim

Common optional argument:
* max_iteration: (default = 100) max number of iterations to perform
* precision: (default = 1d-15) convergence considered when || f'(x) ||_2 < precision
* min_StepLength: (default = 1d-15) terminate if search step < min_StepLength before || f'(x) ||_2 converges

On input x is an initial guess, on exit x is a local minimum of f(x)