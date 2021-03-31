! Augmented Lagrangian multiplier method for equality constraint
! The underlying unconstrained solver is Newton-Raphson
subroutine ALagrangian_Newton_Raphson(f, f_fd, fdd, c, c_cd, c_cd_cdd, x, N, M, &
lambda0, miu0, &
max_iteration, max_StepIteration, precision, min_StepLength)

implicit none

!In each trial miu is updated by miu *= increment
real*8, parameter::increment = 1.05d0

external::f, f_fd, fdd, c, c_cd, c_cd_cdd
interface
    subroutine f(fx, x, N)
        integer*4, intent(in)::N
        real*8, intent(out)::fx
        real*8, dimension(N), intent(in)::x
    end subroutine f

    subroutine f_fd(fx, fdx, x, N)
        integer*4, intent(in)::N
        real*8, intent(out)::fx
        real*8, dimension(N), intent(out)::fdx
        real*8, dimension(N), intent(in)::x
    end subroutine f_fd

    subroutine fdd(fddx, x, N)
        integer*4, intent(in)::N
        real*8, dimension(N, N), intent(out)::fddx
        real*8, dimension(N), intent(in)::x
    end subroutine fdd

    subroutine c(cx, x, M, N)
        integer*4, intent(in)::M, N
        real*8, dimension(M), intent(out)::cx
        real*8, dimension(N), intent(in)::x
    end subroutine c

    subroutine c_cd(cx, cdx, x, M, N)
        integer*4, intent(in)::M, N
        real*8, dimension(M), intent(out)::cx
        real*8, dimension(N, M), intent(out)::cdx
        real*8, dimension(N), intent(in)::x
    end subroutine c_cd

    subroutine c_cd_cdd(cx, cdx, cddx, x, M, N)
        integer*4, intent(in)::M, N
        real*8, dimension(M), intent(out)::cx
        real*8, dimension(N, M), intent(out)::cdx
        real*8, dimension(N, N, M), intent(out)::cddx
        real*8, dimension(N), intent(in)::x
    end subroutine c_cd_cdd
end interface

integer*4, intent(in)::N, M
real*8, dimension(N), intent(inout)::x
real*8, dimension(M), intent(in)::lambda0
real*8, intent(in)::miu0
integer*4, intent(in)::max_iteration
integer*4, intent(in)::max_StepIteration
real*8 , intent(in)::precision, min_StepLength

real*8, dimension(M)::lambda
real*8::miu

integer*4::iIteration
real*8::precision_square
real*8, dimension(M)::cx

!Initialize
lambda = lambda0
miu = miu0
precision_square = precision * precision

!Main loop
do iIteration = 1, max_iteration
    call Newton_Raphson(L, L_Ld, Ldd, x, N, &
              max_StepIteration, precision, min_StepLength)
    call c(cx, x, M, N)
    if (dot_product(cx, cx) < precision_square) exit
    lambda = lambda - miu * cx
    miu = miu * increment
end do

!Warn
if (iIteration > max_iteration) then
    write(*,*)"Failed augmented Lagrangian: max iteration exceeded!"
    write(*,*)"Euclidean norm of constraint violation =", norm2(cx)
end if

contains
subroutine L(Lx, x, N)
    integer*4, intent(in)::N
    real*8, dimension(N), intent(in)::x
    real*8, intent(out)::Lx
    real*8, dimension(M)::cx
    call f(Lx, x, N)
    call c(cx, x, M, N)
    Lx = Lx - dot_product(lambda, cx) + miu / 2d0 * dot_product(cx, cx)
end subroutine L

subroutine L_Ld(Lx, Ldx, x, N)
    integer*4, intent(in)::N
    real*8, intent(out)::Lx
    real*8, dimension(N), intent(out)::Ldx
    real*8, dimension(N), intent(in)::x
    real*8, dimension(M)::cx
    real*8, dimension(N, M)::cdx
    call f_fd(Lx, Ldx, x, N)
    call c_cd(cx, cdx, x, M, N)
    Lx = Lx - dot_product(lambda, cx) + miu / 2d0 * dot_product(cx, cx)
    Ldx = Ldx + matmul(cdx, miu * cx - lambda)
end subroutine L_Ld

subroutine Ldd(Lddx, x, N)
    integer*4, intent(in)::N
    real*8, dimension(N), intent(in)::x
    real*8, dimension(N, N), intent(out)::Lddx
    real*8, dimension(M)::cx
    real*8, dimension(N, M)::cdx
    real*8, dimension(N, N, M)::cddx
    real*8, dimension(N, N)::cddx_term
    integer*4::i
    call fdd(Lddx, x, N)
    call c_cd_cdd(cx, cdx, cddx, x, M, N)
    cx = miu * cx - lambda
    forall (i = 1 : N)
        cddx_term(:, i) = matmul(cddx(i, :, :), cx)
    end forall
    Lddx = Lddx + cddx_term + miu * matmul(cdx, transpose(cdx))
end subroutine Ldd

end subroutine ALagrangian_Newton_Raphson