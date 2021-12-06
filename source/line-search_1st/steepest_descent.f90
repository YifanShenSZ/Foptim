subroutine steepest_descent(f, f_fd, x, dim, &
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
real*8::precision_square, min_StepLength_square, a, fnew, phidold, phidnew
real*8, dimension(dim)::p, fdnew

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
    phidold = phidnew
    !line search
    call strong_Wolfe_1st(f, f_fd, x, a, p, fnew, phidnew, fdnew, dim)
    !check convergence
    phidnew = dot_product(fdnew, fdnew)
    if (phidnew < precision_square) return
    if (dot_product(p, p) * a * a < min_StepLength_square) then
        write(*,*)"Steepest descent warning: step length has converged, but gradient norm has not met accuracy goal"
        write(*,*)"Euclidean norm of gradient =", dSqrt(phidnew)
        return
    end if
    !determine new direction and step length
    p = -fdnew
    phidnew = -phidnew
    a = a * phidold / phidnew
end do

!warn
if (iIteration > max_iteration) then
    write(*,*)"Failed steepest descent: max iteration exceeded!"
    write(*,*)"Final ||gradient|| = ", norm2(fdnew)
end if

end subroutine steepest_descent