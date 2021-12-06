!Newton-Raphson method
subroutine NewtonRaphson(f, f_fd, fdd, x, dim, &
max_iteration, precision, min_StepLength)

use linalg
implicit none

external::f, f_fd, fdd
interface
    subroutine f(fx, x, dim)
        integer*4, intent(in)::dim
        real*8, intent(out)::fx
        real*8, dimension(dim), intent(in)::x
    end subroutine f

    subroutine f_fd(fx, fdx, x, dim)
        integer*4, intent(in)::dim
        real*8, intent(out)::fx
        real*8, dimension(dim), intent(out)::fdx
        real*8, dimension(dim), intent(in)::x
    end subroutine f_fd

    subroutine fdd(fddx, x, dim)
        integer*4, intent(in)::dim
        real*8, dimension(dim, dim), intent(out)::fddx
        real*8, dimension(dim), intent(in)::x
    end subroutine fdd
end interface

integer*4, intent(in)::dim
real*8, dimension(dim), intent(inout)::x
integer*4, intent(in)::max_iteration
real*8 , intent(in)::precision, min_StepLength

integer*4::iIteration, po
real*8::precision_square, min_StepLength_square, a, fnew, phidnew
real*8, dimension(dim)::p, fdnew
real*8, dimension(dim, dim)::Hinv

!Initialize
precision_square = precision * precision
min_StepLength_square = min_StepLength * min_StepLength
call f_fd(fnew, fdnew, x, dim)
call fdd(Hinv, x, dim)
po = My_dpotri(Hinv, dim)
if (po /= 0) then
    write(*,*)"Failed Newton-Raphson: Hessian is not positive-definite"
    write(*,*)"Final ||gradient|| = ", norm2(fdnew)
    return
end if
p = -matmul(Hinv, fdnew)
phidnew = dot_product(fdnew, p)
a = 1d0

!Main loop
do iIteration = 1, max_iteration
    !Line search
    call strong_Wolfe_2nd(f, f_fd, x, a, p, fnew, phidnew, fdnew, dim)
    !Check convergence
    phidnew = dot_product(fdnew, fdnew)
    if (phidnew < precision_square) return
    if (dot_product(p, p) * a * a < min_StepLength_square) then
        write(*,*)"Newton-Raphson warning: step length has converged, but gradient norm has not met accuracy goal"
        write(*,*)"Euclidean norm of gradient =", dSqrt(phidnew)
        return
    end if
    !Update search direction
    call fdd(Hinv, x, dim)
    po = My_dpotri(Hinv, dim)
    if (po /= 0) then
        write(*,*)"Failed Newton-Raphson: Hessian is not positive-definite"
        write(*,*)"Final ||gradient|| = ", norm2(fdnew)
        return
    end if
    p = -matmul(Hinv, fdnew)
    phidnew = dot_product(fdnew, p)
    a = 1d0
end do

!Warn
if (iIteration > max_iteration) then
    write(*,*)"Failed Newton-Raphson: max iteration exceeded!"
    write(*,*)"Final ||gradient|| = ", norm2(fdnew)
end if

end subroutine NewtonRaphson