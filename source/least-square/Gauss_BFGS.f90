!Gauss-Newton method for least square problem,
!replacing Newton-Raphson with BFGS
subroutine Gauss_BFGS(residue, Jacobian, x, M, N, &
max_iteration, precision, min_StepLength)

use linalg
implicit none

external::residue, Jacobian
interface 
    subroutine residue(r, x, M, N)
        integer*4, intent(in)::M, N
        real*8, dimension(M), intent(out)::r
        real*8, dimension(N), intent(in)::x
    end subroutine residue

    subroutine Jacobian(J, x, M, N)
        integer*4, intent(in)::M, N
        real*8, dimension(M, N), intent(out)::J
        real*8, dimension(N), intent(in)::x
    end subroutine Jacobian
end interface

integer*4, intent(in)::M, N
real*8, dimension(N), intent(inout)::x
integer*4, intent(in)::max_iteration
real*8 , intent(in)::precision, min_StepLength

integer*4::iIteration, po, i
real*8::precision_square, min_StepLength_square, a, fnew, phidnew, rho
real*8, dimension(M)::r
real*8, dimension(M, N)::J
real*8, dimension(N)::p, fdnew, s, y
real*8, dimension(N, N)::U, Hinv

!initialize
precision_square = 0.5d0 * precision * precision
min_StepLength_square = min_StepLength * min_StepLength
call residue(r, x, M, N)
call Jacobian(J, x, M, N)
fnew = 0.5d0 * dot_product(r, r)
fdnew = matmul(transpose(J), r)
!initial approximate inverse Hessian & direction & step length
Hinv = matmul(transpose(J), J)
p = -fdnew
po = My_dpotri(Hinv, N)
call dsyL2U(Hinv, N)
p = -matmul(Hinv, fdnew)
phidnew = dot_product(fdnew, p)
a = 1d0

!main loop
do iIteration = 1, max_iteration
    !prepare
    s = x
    y = fdnew
    !line search
    call strong_Wolfe_2nd(merit, merit_gradient, x, a, p, fnew, phidnew, fdnew, N)
    !check convergence
    if (fnew < precision_square) return
    if (dot_product(p, p) * a * a < min_StepLength_square) then
        write(*,*)"Gauss-BFGS warning: step length has converged, but gradient norm has not met accuracy goal"
        call residue(r, x, M, N)
        write(*,*)"Final residual = ", norm2(r)
        return
    end if
    !determine new direction and step length, update approximate inverse Hessian
    s = x - s
    y = fdnew - y
    rho = 1d0 / dot_product(y, s)
    U = -rho * vector_direct_product(y, s, N, N)
    forall(i = 1 : N)
        U(i, i) = U(i, i) + 1d0
    end forall
    Hinv = matmul(transpose(U), matmul(Hinv, U)) &
         + rho * vector_direct_product(s, s, N, N)
    p = -matmul(Hinv, fdnew)
    phidnew = dot_product(fdnew, p)
    a = 1d0
end do

!warn
if (iIteration > max_iteration) then
    write(*,*)"Failed Gauss-BFGS: max iteration exceeded!"
    call residue(r, x, M, N)
    write(*,*)"Final residual = ", norm2(r)
end if

contains
subroutine merit(fx, x, N)
    integer*4, intent(in)::N
    real*8, intent(out)::fx
    real*8, dimension(N), intent(in)::x
    real*8, dimension(M)::r
    call residue(r, x, M, N)
    fx = 0.5d0 * dot_product(r, r)
end subroutine merit

subroutine merit_gradient(fx, fdx, x, N)
    integer*4, intent(in)::N
    real*8, intent(out)::fx
    real*8, dimension(N), intent(out)::fdx
    real*8, dimension(N), intent(in)::x
    real*8, dimension(M)::r
    real*8, dimension(M, N)::J
    call residue(r, x, M, N)
    call Jacobian(J, x, M, N)
    fx = 0.5d0 * dot_product(r, r)
    fdx = matmul(transpose(J), r)
end subroutine merit_gradient

end subroutine Gauss_BFGS