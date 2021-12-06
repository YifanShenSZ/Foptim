!Dai-Yuan conjugate gradient
subroutine CGDY_verbose(f, f_fd, x, dim, &
initial_learning_rate, max_iteration, precision, min_StepLength)

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
end interface

integer*4, intent(in)::dim
real*8, dimension(dim), intent(inout)::x
real*8, intent(in)::initial_learning_rate
integer*4, intent(in)::max_iteration
real*8, intent(in)::precision, min_StepLength

integer*4::iIteration
real*8::precision_square, min_StepLength_square, a, fnew, fdnewfdnew, phidold, phidnew
real*8, dimension(dim)::p, fdold, fdnew

!initialize
precision_square = precision * precision
min_StepLength_square = min_StepLength * min_StepLength
call f_fd(fnew, fdnew, x, dim)
!initial direction and step length
p = -fdnew
phidnew = -dot_product(fdnew, fdnew)
if (-phidnew < precision_square) return
a = initial_learning_rate

!main loop
do iIteration = 1, max_iteration
    !prepare
    fdold = fdnew
    phidold = phidnew
    !line search
    call strong_Wolfe_1st(f, f_fd, x, a, p, fnew, phidnew, fdnew, dim)
    !check convergence
    fdnewfdnew = dot_product(fdnew, fdnew)
    if (fdnewfdnew < precision_square) return
    if (dot_product(p, p) * a * a < min_StepLength_square) then
        write(*,*)"Dai-Yuan conjugate gradient warning: step length has converged, but gradient norm has not met accuracy goal"
        write(*,*)"Euclidean norm of gradient =", dSqrt(fdnewfdnew)
        return
    end if
    !determine new direction
    p = -fdnew + fdnewfdnew / dot_product(fdnew - fdold, p) * p
    phidnew = dot_product(fdnew, p)
    !ascent direction, reset to steepest descent direction
    if (phidnew > 0d0) then
        p = -fdnew
        phidnew = -fdnewfdnew
    end if
    !determine new step length
    a = a * phidold / phidnew
    !verbosity
    write(*,*)
    call show_time()
    write(*,*)"Iteration", iIteration, ":"
    write(*,*)"* loss", fnew
    write(*,*)"* ||gradient||", dSqrt(-phidnew)
    write(*,*)"* learning rate ", a
end do

!warn
if (iIteration > max_iteration) then
    write(*,*)"Failed Dai-Yuan conjugate gradient: max iteration exceeded!"
    write(*,*)"Final ||gradient|| = ", norm2(fdnew)
end if

end subroutine CGDY_verbose