# Trust region solver
MKL trust-region nonlinear least square problem (trnlsp) solver wrapper

Minimize r(x) . r(x) through trust-region method with model function m(p) = [ r(x) + J(x) . p ]^2, where J(x) is the Jacobian

External procedure:
* subroutine residue(r(x), x, M, N)
* subroutine Jacobian(J(x), x, M, N)
* M dimensional vector r(x), N dimensional vector x, M x N matrix J(x), M >= N

Required arguments:
* subroutine residue
* subroutine Jacobian
* N dimensional vector x
* integer*4 M & N

Optional arguments:
* max_iteration: (default = 100) max number of iterations to perform
* max_StepIteration: (default = 100) max number of iterations for determining each step length
* precision: (default = 1d-15) convergence considered when || r(x) ||_2 < precision
* min_StepLength: (default = 1d-15) terminate if search step < min_StepLength before || r(x) ||_2 converges

On input x is an initial guess, on exit x is a local minimum of r(x) . r(x)