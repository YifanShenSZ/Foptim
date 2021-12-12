#include "mkl_rci.f90"

subroutine trust_region(residue, Jacobian, x, M, N, &
max_iteration, max_StepIteration, precision, min_StepLength)

use mkl_rci_type
use mkl_rci
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
integer*4, intent(in)::max_iteration, max_StepIteration
real*8 , intent(in)::precision, min_StepLength

!Reverse communication interface (RCI)
integer*4::RCI_request, & !Recieve job request
           brief          !Brief of initialization and input parameter checking
integer*4, dimension(6)::info !Results of input parameter checking
type(handle_tr)::handle !Trust-region solver handle
!Job control
!tol(1:5) contains the stopping criteria for solving f'(x) = 0:
!    1, trust region radius < tol(1)
!    2, || f'(x) ||_2 < tol(2)
!    3, || Jacobian ||_1 < tol(3) 
!    4, || s ||_2 < tol(4), where s is the trial step
!    5, || f'(x) ||_2 - || f'(x) - Jacobian . s ||_2 < tol(5)
!tol(6) is the precision of s calculation
real*8, dimension(6)::tol
!total_iteration harvests the solver stops after how many iterations
!stop_reason harvests why the solver has stopped:
!    1,   max iteration exceeded
!    2-6, tol(stop_reason - 1) is met
integer*4::total_iteration, stop_reason
real*8::init_residue, final_residue, step_bound
real*8, dimension(M)::fdx
real*8, dimension(M, N)::J

!Initialize
tol = [min_StepLength, precision, 1d-12, min_StepLength, min_StepLength, 1d-12]
call residue(fdx, x, M, N)
call Jacobian(J, x, M, N)
step_bound = 100d0
RCI_request = 0
brief = dtrnlsp_init(handle, N, M, x, tol, max_iteration, max_StepIteration, step_bound)
if (brief /= TR_SUCCESS) then
    write(*,*)"Trust region abort: invalid initialization"
    if (brief == TR_INVALID_OPTION) then
        write(*,*)"There was an error in the input parameters"
    else if (brief == TR_OUT_OF_MEMORY) then
        write(*,*)"There was a memory error"
    else
        write(*,*)"There was an unexpected crash"
    end if
    call mkl_free_buffers
    return
end if
brief = dtrnlsp_check(handle, N, M, J, fdx, tol, info)
if (brief /= TR_SUCCESS) then
    write(*,*)"Trust region abort: check failed"
    call mkl_free_buffers
    return
else if(info(1) /= 0 .or. info(2) /= 0 .or. info(3) /= 0 .or. info(4) /= 0) then
    write(*,*)"Trust region abort: check was not passed, the information is:"
    write(*,*)info
    call mkl_free_buffers
    return
end if

!Main loop
do
    if (dtrnlsp_solve(handle, fdx, J, RCI_request) /= TR_SUCCESS) then
        call mkl_free_buffers
        return
    end if
    select case (RCI_request)
    case (-1,-2,-3,-4,-5,-6); exit
    case (1); call residue(fdx, x, M, N)
    case (2); call Jacobian(J, x, M, N)
    end select
end do

!Clean up
if (dtrnlsp_get(handle, total_iteration, stop_reason, init_residue, final_residue) /= TR_SUCCESS) then
    call mkl_free_buffers
    return
end if
if (dtrnlsp_delete(handle) /= TR_SUCCESS) then
    call mkl_free_buffers
    return
end if
call mkl_free_buffers

!Warn
if (stop_reason /= 3) then
    select case(stop_reason)
    case(1); write(*,*)"Failed trust region: max iteration exceeded!"
    case(4); write(*,*)"Failed trust region: singular Jacobian encountered!"
    case default; write(*,*)"Trust region warning: step length has converged, but residual has not met accuracy goal"
    end select
    write(*,*)"Final residual = ", final_residue
end if
end subroutine trust_region