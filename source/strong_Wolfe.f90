!Line search for a step length satisfying strong Wolfe condition
!Input:  x is current x
!        a is initial guess of a
!        p is current p
!        fx = f(x)
!        phid0 = phi'(0)
!Output: a harvests the step length satisfying certain condition
!        x = x + a * p
!        fx = f(x)
!        fdx = f'(x)
subroutine strong_Wolfe(f, f_fd, x, a, p, fx, phid0, fdx, dim)

implicit none

!c1 and c2 are Wolfe parameters, who should satisfy
!* 0 < c1 < c2 <  1  for Newton & quasi-Newton
!* 0 < c1 < c2 < 0.5 for steepest descent & conjugate gradient
!This c2 is best for Newton & quasi-Newton
real*8, parameter::c1 = 1d-4, c2 = 0.9d0
!In each trial a is updated by a *= (/=) increment
real*8, parameter::increment = 1.05d0

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
real*8, intent(inout)::a
real*8, dimension(dim), intent(in)::p
real*8, intent(inout)::fx
real*8, intent(in)::phid0
real*8, dimension(dim), intent(out)::fdx

real*8::c2_m_AbsPhid0, fx0, ftemp, atemp, aold, fold, phidnew, phidold
real*8, dimension(dim)::x0

!Initialize
x0 = x
fx0 = fx
c2_m_AbsPhid0 = c2 * dAbs(phid0)

!Check whether initial guess satisfies the sufficient decrease condition
x = x0 + a * p
call f_fd(fx, fdx, x, dim)
!Satisfied
if (fx <= fx0 + c1 * a * phid0) then
    phidnew = dot_product(fdx, p)
    !Curve is heading up, search for a smaller a
    if (phidnew > 0d0) then
        !a also satisfies the strong curvature condition
        if (dAbs(phidnew) <= c2_m_AbsPhid0) return
        !else search for a smaller a that phi(a) >= phi(aold) or phid(a) <= 0
        do
            aold = a
            fold = fx
            phidold = phidnew
            a = aold / increment
            x = x0 + a * p
            call f_fd(fx, fdx, x, dim)
            phidnew = dot_product(fdx, p)
            !Found such an a
            if (fx >= fold .or. phidnew <= 0d0) then
                atemp = a
                ftemp = fx
                call zoom(aold, atemp, fold, ftemp, phidold, phidnew)
                return
            end if
            !Vanising step length, stop searching
            if (a < 1d-15) return
        end do
    !Curve still heads down, search for a larger a
    else
        do
            aold = a
            fold = fx
            phidold = phidnew
            a = aold * increment
            x = x0 + a * p
            call f_fd(fx, fdx, x, dim)
            phidnew = dot_product(fdx, p)
            !Current a breaks the sufficient decrease condition or phi(a) >= phiold(a)
            if (fx > fx0 + c1 * a * phid0 .or. fx >= fold) then
                atemp = a
                ftemp = fx
                call zoom(aold, atemp, fold, ftemp, phidold, phidnew)
                return
            end if
            !Curve is heading up
            if (phidnew > 0d0) then
                !Current a satisfies the strong curvature condition
                if (dAbs(phidnew) <= c2_m_AbsPhid0) return
                atemp = a
                ftemp = fx
                call zoom(atemp, aold, ftemp, fold, phidnew, phidold)
                return
            end if
        end do
    end if
!Violated, search for a smaller a
else
    do
        aold = a
        fold = fx
        a = aold / increment
        x = x0 + a * p
        call f(fx, x, dim)
        !Current a satisfies the sufficient decrease condition, then look at slope
        if (fx <= fx0 + c1 * a *phid0) then
            call f_fd(fx, fdx, x, dim)
            phidnew = dot_product(fdx, p)
            !a also satisfies the strong curvature condition
            if (dAbs(phidnew) <= c2_m_AbsPhid0) return
            !Curve is heading down, the true a lies within (a, aold)
            if (phidnew < 0d0) then
                x = x0 + aold * p
                call f_fd(fx, fdx, x, dim)
                phidold = dot_product(fdx, p)
                atemp = a
                ftemp = fx
                call zoom(atemp, aold, ftemp, fold, phidnew, phidold)
                return
            !Curve is heading up, search for a smaller a that phi(a) >= phi(aold) or phid(a) <= 0
            else
                do
                    aold = a
                    fold = fx
                    phidold = phidnew
                    a = aold / increment
                    x = x0 + a * p
                    call f_fd(fx, fdx, x, dim)
                    phidnew = dot_product(fdx, p)
                    !Found such an a
                    if (fx >= fold .or. phidnew <= 0d0) then
                        atemp = a
                        ftemp = fx
                        call zoom(aold, atemp, fold, ftemp, phidold, phidnew)
                        return
                    end if
                    !Vanising step length, stop searching
                    if(a<1d-15) return
                end do
            end if
        end if
        !Vanising step length, stop searching
        if (a < 1d-15) then
            call f_fd(fx, fdx, x, dim)
            return
        end if
    end do
end if

contains
!low & up must satisfy:
!* low satisfies sufficient decrease condition
!* (up - low) * phi'(low) < 0
!* (low, up) (or (up, low)) contains a step length satisfying strong Wolfe condition,
!  which means at least 1 of 3 following statements is true:
!  * up violates the sufficient decrease condition
!  * phi(up) >= phi(low)
!  * up < low & phi'(up) <= 0
subroutine zoom(low, up, flow, fup, phidlow, phidup)
    real*8, intent(inout)::low, up, flow, fup, phidlow, phidup
    real*8::phidnew, d1, d2
    do
        !Updata a by cubic interpolation
        d1 = phidlow + phidup - 3d0 * (flow - fup) / (low - up)
        if (up > low) then
            d2 =  dSqrt(d1 * d1 - phidlow * phidup)
        else
            d2 = -dSqrt(d1 * d1 - phidlow * phidup)
        end if
        a = up - (up - low) * (phidup + d2 - d1) / (phidup - phidlow + 2d0 * d2)
        !Fail safe
        if (.not.(a > min(low, up) .and. a < max(low, up))) a = (low + up) / 2d0
        !Check new a
        x = x0 + a * p
        call f_fd(fx, fdx, x, dim)
        phidnew = dot_product(fdx, p)
        !Current a breaks the sufficient decrease condition or phi(a) >= phi(low)
        if (fx > fx0 + c1 * a * phid0 .or. fx >= flow) then
            up = a
            fup = fx
            phidup = phidnew
        !Current a satisfies the sufficient decrease condition
        else
            !Current a satisfies the curvature condition
            if (dAbs(phidnew) <= c2_m_AbsPhid0) return
            !Exchange low and up
            if (phidnew * (up - low) >= 0d0) then
                up = low
                fup = flow
                phidup = phidlow
            end if
            low = a
            flow = fx
            phidlow = phidnew
        end if
        !Vanising range, stop searching
        if (dAbs(up - low) < 1d-15 .or. dAbs(up - low) / max(dAbs(low), dAbs(up)) < 1d-15) return
    end do
end subroutine zoom

end subroutine strong_Wolfe