!Line search for a step length satisfying Wolfe condition
!This routine is designed to minimize gradient computation
!Input:  x is current x
!        a is initial guess of a
!        p is current p
!        fx = f(x)
!        phid0 = phi'(0)
!Output: a harvests the step length satisfying certain condition
!        x = x + a * p
!        fx = f(x)
!        fdx = f'(x)
subroutine Wolfe(f, fd, x, a, p, fx, phid0, fdx, dim)

implicit none

!c1 and c2 are Wolfe parameters, who should satisfy
!* 0 < c1 < c2 <  1  for Newton & quasi-Newton
!* 0 < c1 < c2 < 0.5 for steepest descent & conjugate gradient
!This c2 is best for steepest descent & conjugate gradient
real*8, parameter::c1 = 1d-4, c2 = 0.45d0
!In each trial a is updated by a *= (/=) increment
real*8, parameter::increment = 1.05d0

external::f, fd
interface
    subroutine f(fx, x, dim)
        integer*4, intent(in)::dim
        real*8, intent(out)::fx
        real*8, dimension(dim), intent(in)::x
    end subroutine f

    subroutine fd(fdx, x, dim)
        integer*4, intent(in)::dim
        real*8, dimension(dim), intent(out)::fdx
        real*8, dimension(dim), intent(in)::x
    end subroutine fd
end interface

integer*4, intent(in)::dim
real*8, dimension(dim), intent(inout)::x
real*8, intent(inout)::a
real*8, dimension(dim), intent(in)::p
real*8, intent(inout)::fx
real*8, intent(in)::phid0
real*8, dimension(dim), intent(out)::fdx

real*8::c2_m_phid0, fx0, ftemp, atemp, aold, fold, phidx
real*8,dimension(dim)::x0

!Initialize
x0 = x
fx0 = fx
c2_m_phid0 = c2 * phid0

!Check whether initial guess satisfies the sufficient decrease condition
x = x0 + a * p
call f(fx, x, dim)
!Satisfied, search for a larger a
if (fx <= fx0 + c1 * a * phid0) then
    do
        aold = a
        fold = fx
        a = aold * increment
        x = x0 + a * p
        call f(fx, x, dim)
        !Found an upper limit who violates the sufficient decrease condition
        if (fx > fx0 + c1 * a * phid0) then
            x = x0 + aold * p
            call fd(fdx, x, dim)
            phidx = dot_product(fdx, p)
            !aold satisfies the curvature condition
            if(phidx >= c2_m_phid0) then
                a = aold
                fx = fold
            !Search for an a sastisfying the curvature condition within (aold, a)
            else
                atemp = a
                ftemp = fx
                call zoom(aold, atemp, fold, ftemp, phidx)
            end if
            return
        end if
    end do
!Violated, search for a smaller a
else
    do
        aold = a
        fold = fx
        a = aold / increment
        x = x0 + a * p
        call f(fx, x, dim)
        !Found a lower limit who satisfies the sufficient decrease condition
        if (fx <= fx0 + c1 * a * phid0) then
            call fd(fdx, x, dim)
            phidx = dot_product(fdx, p)
            !Search for an a sastisfying the curvature condition within (a, aold)
            if(phidx < c2_m_phid0) then
                atemp = a
                ftemp = fx
                call zoom(atemp, aold, ftemp, fold, phidx)
            end if
            !else a already satisfies the curvature condition
            return
        end if
        !Vanising step length, stop searching
        if (a < 1d-12) then
            call fd(fdx, x, dim)
            return
        end if
    end do
end if

contains
!low & up must satisfy:
!* low < up
!* low satisfies the sufficient decrease condition, but up violates
!* phi'(low) < 0
subroutine zoom(low, up, flow, fup, phidlow)
    real*8,intent(inout)::low, up, flow, fup, phidlow
    real*8::phidnew, phidlow_m_a
    !Initialize
    phidlow_m_a = phidlow * a
    !Main loop
    do
        !Vanising range, stop searching
        if(up - low < 1d-12 .or. (up - low) / max(dAbs(low), dAbs(up)) < 1d-12) then
            call fd(fdx, x, dim)
            return
        end if
        !Updata a by quadratic interpolation
        a = phidlow_m_a * a /2d0 / (flow + phidlow_m_a - fup)
        !Fail safe
        if (.not.( a > low .and. a < up)) a = (low + up) / 2d0
        !Check whether current a satisfies the sufficient decrease condition
        x = x0 + a * p
        call f(fx, x, dim)
        !Violated, current a becomes the upper limit
        if(fx > fx0 + c1 * a * phid0) then
            up = a
            fup = fx
        !Satisfied, further check the curvature condition
        else
            call fd(fdx, x, dim)
            phidnew = dot_product(fdx, p)
            !Satisfied, mission success
            if(phidnew >= c2_m_phid0) return
            !Violated, current a becomes the lower limit
            low = a
            flow = fx
            phidlow = phidnew
            phidlow_m_a = phidlow*a
        end if
    end do
end subroutine zoom

end subroutine Wolfe