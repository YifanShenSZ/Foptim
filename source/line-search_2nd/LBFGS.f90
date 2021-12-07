!Limited-memory Broyden–Fletcher–Goldfarb–Shanno (L-BFGS) quasi-Newton method
subroutine LBFGS(f, f_fd, x, dim, &
initial_learning_rate, memory, max_iteration, precision, min_StepLength)

use linalg
implicit none

external::f, f_fd
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
end interface

integer*4, intent(in)::dim
real*8, dimension(dim), intent(inout)::x
real*8, intent(in)::initial_learning_rate
integer*4, intent(in)::memory, max_iteration
real*8, intent(in)::precision, min_StepLength

integer*4::recent, iIteration, i
real*8::precision_square, min_StepLength_square, a, fnew, phidnew
real*8, dimension(dim)::xold, p, fdold, fdnew
real*8, dimension(0 : memory - 1)::rho, alpha
real*8, dimension(dim, 0 : memory - 1)::s, y

!initialize
precision_square = precision * precision
min_StepLength_square = min_StepLength * min_StepLength
call f_fd(fnew, fdnew, x, dim)
!initial iteration history
p = -fdnew
phidnew = -dot_product(fdnew, fdnew)
if (-phidnew < precision_square) return
a = initial_learning_rate
!initial approximate inverse Hessian = a,
!perform a steepest descent step to find a
xold = x
fdold = fdnew
call strong_Wolfe_1st(f, f_fd, x, a, p, fnew, phidnew, fdnew, dim)
phidnew = dot_product(fdnew,fdnew)
if (phidnew < precision_square) return
if (dot_product(p, p) * a * a < min_StepLength_square) then
    write(*,*)"LBFGS warning: step length has converged, but gradient norm has not met accuracy goal"
    write(*,*)"Euclidean norm of gradient =", dSqrt(phidnew)
    return
end if
!preiterate to get enough history
recent = 0
s(:, 0) = x - xold
y(:, 0) = fdnew - fdold
rho(0) = 1d0 / dot_product(y(:, 0), s(:, 0))
do iIteration = 1, memory - 1
    !prepare
    xold = x
    fdold = fdnew
    !determine new direction
    p = fdnew
    do i = recent, 0, -1
        alpha(i) = rho(i) * dot_product(s(:, i), p)
        p = p - alpha(i) * y(:, i)
    end do
    p = p / rho(recent) / dot_product(y(:, recent), y(:, recent))
    do i = 0, recent
        phidnew = rho(i) * dot_product(y(:, i), p)
        p = p + (alpha(i) -phidnew) * s(:, i)
    end do
    p = -p
    phidnew = dot_product(fdnew, p)
    a = 1d0
    !line search
    call strong_Wolfe_1st(f, f_fd, x, a, p, fnew, phidnew, fdnew, dim)
    !check convergence
    phidnew = dot_product(fdnew, fdnew)
    if (phidnew < precision_square) return
    if (dot_product(p, p) * a * a < min_StepLength_square) then
        write(*,*)"LBFGS warning: step length has converged, but gradient norm has not met accuracy goal"
        write(*,*)"Euclidean norm of gradient =", dSqrt(phidnew)
        return
    end if
    !memorize history
    recent = recent + 1
    s(:, recent) = x - xold
    y(:, recent) = fdnew - fdold
    rho(recent) = 1d0 / dot_product(y(:, recent), s(:, recent))
end do

!main loop
do iIteration = 1, max_iteration
    !prepare
    xold = x
    fdold = fdnew
    !determine new direction
    p = fdnew
    do i = recent, 0, -1
        alpha(i) = rho(i) * dot_product(s(:, i), p)
        p = p - alpha(i) * y(:, i)
    end do
    do i = memory - 1, recent + 1, -1
        alpha(i) = rho(i) * dot_product(s(:, i), p)
        p = p - alpha(i) * y(:, i)
    end do
    p = p / rho(recent) /dot_product(y(:, recent), y(:, recent))
    do i = recent + 1, memory - 1
        phidnew = rho(i) * dot_product(y(:, i), p)
        p = p + (alpha(i) -phidnew) * s(:, i)
    end do
    do i = 0, recent
        phidnew = rho(i) * dot_product(y(:, i), p)
        p = p + (alpha(i) - phidnew) * s(:, i)
    end do
    p = -p
    phidnew = dot_product(fdnew, p)
    a = 1d0
    !line search
    call strong_Wolfe_2nd(f, f_fd, x, a, p, fnew, phidnew, fdnew, dim)
    !check convergence
    phidnew = dot_product(fdnew, fdnew)
    if (phidnew < precision_square) return
    if (dot_product(p, p) * a * a < min_StepLength_square) then
        write(*,*)"LBFGS warning: step length has converged, but gradient norm has not met accuracy goal"
        write(*,*)"Euclidean norm of gradient =", dSqrt(phidnew)
        return
    end if
    !replace earliest history with latest
    recent = mod(recent + 1, memory)
    s(:, recent) = x - xold
    y(:, recent) = fdnew - fdold
    rho(recent) = 1d0 / dot_product(y(:, recent), s(:, recent))
end do

!warn
if (iIteration > max_iteration) then
    write(*,*)"Failed LBFGS: max iteration exceeded!"
    write(*,*)"Final ||gradient|| = ", norm2(fdnew)
end if

end subroutine LBFGS