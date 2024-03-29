!Broyden–Fletcher–Goldfarb–Shanno (BFGS) quasi-Newton method
subroutine BFGS(f, f_fd, fdd, x, dim, &
Hessian_step, max_iteration, precision, min_StepLength)

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
integer*4, intent(in)::Hessian_step, max_iteration
real*8, intent(in)::precision, min_StepLength

integer*4::iIteration, po, i, Hessian_time
real*8::precision_square, min_StepLength_square, a, fnew, phidnew, rho, ys
real*8, dimension(dim)::p, fdnew, s, y
real*8, dimension(dim, dim)::U, Hinv

!initialize
precision_square = precision * precision
min_StepLength_square = min_StepLength * min_StepLength
call f_fd(fnew, fdnew, x, dim)
!initial approximate inverse Hessian & direction & step length
call fdd(Hinv, x, dim)
p = -fdnew
po = My_dpotri(Hinv, dim)
if (po == 0) then
    call dsyL2U(Hinv, dim)
    p = -matmul(Hinv, fdnew)
    phidnew = dot_product(fdnew, p)
    a = 1d0
!exact Hessian is not positive definite, let initial approximate inverse Hessian = a
else
    !perform a steepest descent step to find a
    p = -fdnew
    phidnew = -dot_product(fdnew, fdnew)
    if (-phidnew < precision_square) return
    !although initial learning rate may need some trial and error,
    !since this is just a workaround for bad initial Hessian,
    !user should try a 1st-order optimizer and tune its initial learning rate to walk out of the bad Hessian region
    !instead of wasting time here with BFGS
    a = 1d0
    s = x
    y = fdnew
    call strong_Wolfe_1st(f, f_fd, x, a, p, fnew, phidnew, fdnew, dim)
    phidnew = dot_product(fdnew,fdnew)
    if (phidnew < precision_square) return
    if (dot_product(p, p) * a * a < min_StepLength_square) then
        write(*,*)"BFGS warning: step length has converged, but gradient norm has not met accuracy goal"
        write(*,*)"Euclidean norm of gradient =", dSqrt(phidnew)
        return
    end if
    !update inverse Hessian
    s = x - s
    y = fdnew - y
    ys = dot_product(y, s)
    if (ys < 1d-24) then
        write(*,*)"Failed BFGS: curvature condition is violated"
        write(*,*)"y . s =", ys
        return
    end if
    rho = 1d0 / ys
    U = -rho * vector_direct_product(y, s, dim, dim)
    forall (i = 1 : dim)
        U(i, i) = U(i, i) + 1d0
    end forall
    Hinv = matmul(transpose(U), a * U) &
         + rho * vector_direct_product(s, s, dim, dim)
    !set direction & step length
    p = -matmul(Hinv, fdnew)
    phidnew = dot_product(fdnew, p)
    a = 1d0
end if

!main loop
do iIteration = 1, max_iteration
    !prepare
    s = x
    y = fdnew
    !line search
    call strong_Wolfe_2nd(f, f_fd, x, a, p, fnew, phidnew, fdnew, dim)
    !check convergence
    phidnew = dot_product(fdnew, fdnew)
    if (phidnew < precision_square) return
    if (dot_product(p, p) * a * a < min_StepLength_square) then
        write(*,*)"BFGS warning: step length has converged, but gradient norm has not met accuracy goal"
        write(*,*)"Euclidean norm of gradient =", dSqrt(phidnew)
        return
    end if
    !determine new direction and step length, update approximate inverse Hessian
    Hessian_time = mod(iIteration, Hessian_step)
    !every `Hessian_step` steps compute exact Hessian
    if (Hessian_time == 0) then
        call fdd(U, x, dim)
        po = My_dpotri(U, dim)
        !use exact Hessian if positive definite
        if (po == 0) then
            Hinv = U
            call dsyL2U(Hinv, dim)
            p = -matmul(Hinv, fdnew)
            phidnew = dot_product(fdnew, p)
            a = 1d0
        end if
    end if
    !exact Hessian is either uncomputed or not positive definite, update approximate Hessian
    if (Hessian_time /= 0 .or. po /= 0) then
        s = x - s
        y = fdnew - y
        ys = dot_product(y, s)
        if (ys < 1d-24) then
            write(*,*)"Failed BFGS: curvature condition is violated"
            write(*,*)"y . s =", ys
            return
        end if
        rho = 1d0 / ys
        U = -rho * vector_direct_product(y, s, dim, dim)
        forall(i = 1 : dim)
            U(i, i) = U(i, i) + 1d0
        end forall
        Hinv = matmul(transpose(U), matmul(Hinv, U)) &
             + rho * vector_direct_product(s, s, dim, dim)
        p = -matmul(Hinv, fdnew)
        phidnew = dot_product(fdnew, p)
        a = 1d0
    end if
end do

!warn
if (iIteration > max_iteration) then
    write(*,*)"Failed BFGS: max iteration exceeded!"
    write(*,*)"Final ||gradient|| = ", norm2(fdnew)
end if

end subroutine BFGS