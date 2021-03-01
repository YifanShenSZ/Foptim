# Line search with equality constraint
Lagrangian multiplier method is a classical way to treat equality constraint:

    L = f - lambda . c

where f is the target function to be minimized, lambda is Lagrangian multiplier, c is equality constraint c(x) = 0

Textbook is wrong:
1. it claims Lagrangian multiplier method transforms constrained optimization into unconstrained one
2. However, L has no lower bound, since lambda . c can approach infinity when c != 0 and lambda diverges
3. Lagrangian multiplier method actually turns a minimization problem into a saddle point problem, which cannot necessarily be solved through decreasing L, deteriorating all unconstrained minimizers

Lagrangian multiplier method (minimizing || L'(x) ||) is numerically feasible only when at least 1 of the following statements is true:
* L has unique saddle point
* The initial guess is sufficiently close to the exact solution

In general case, we have to turn to the augmented Lagrangian method:

    Augmented Lagrangian = f - lambda . c + miu / 2 * c . c

where miu is constraint violation penalty strength, who should >= 1 because c ~ ( lambda - lambda_true ) / miu

Available solvers:
* Lagrangian
* augmented Lagrangian

Suggestion:
* Augmented Lagrangian has ill conditioned Hessian when miu is too large, deteriorating performance of line searchers
* so do not take too much iterations nor push accuracy to double precision limit

External procedure:
* subroutine f(f(x), x, dim)
* subroutine fd(f'(x), x, dim)
* subroutine f_fd(f(x), f'(x), x, dim)
* subroutine fdd(f''(x), x, dim)
* subroutine c(c(x), x, M, N)
* subroutine c_cd(c'(x), c'(x), x, M, N)
* subroutine c_cd_cdd(c(x), c'(x), c''(x), x, M, N)
* N dimensional vector x & f'(x), N order matrix f''(x)
* M dimensional vector c, N x M matrix c'(x), N x N x M 3rd-order tensor c''(x)

Required arguments:
* subroutine f & fd & c & cd
* N dimensional vector x
* integer N & M

On input x is an initial guess, on exit x is a local minimum of f(x) subject to constraint c(x)