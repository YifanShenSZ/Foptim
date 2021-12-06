program main

implicit none

real*8, dimension(10)::x
real*8, dimension(1)::lambda0

call random_seed()

write(*,*)"Steepest descent"
call random_number(x)
call steepest_descent(f, f_fd, x, 10, &
     1d-3, 100, 1d-12, 1d-12)
write(*,*)norm2(x)
write(*,*)

write(*,*)"Steepest descent verbose"
call random_number(x)
call steepest_descent_verbose(f, f_fd, x, 10, &
     1d-3, 100, 1d-12, 1d-12)
write(*,*)norm2(x)
write(*,*)

write(*,*)"Dai-Yuan conjugate gradient"
call random_number(x)
call CGDY(f, f_fd, x, 10, &
     1d-3, 100, 1d-12, 1d-12)
write(*,*)norm2(x)
write(*,*)

write(*,*)"Dai-Yuan conjugate gradient"
call random_number(x)
call CGDY_verbose(f, f_fd, x, 10, &
     1d-3, 100, 1d-12, 1d-12)
write(*,*)norm2(x)
write(*,*)

write(*,*)"Polak-Ribiere+ conjugate gradient"
call random_number(x)
call CGPR(f, f_fd, x, 10, &
     1d-3, 100, 1d-12, 1d-12)
write(*,*)norm2(x)
write(*,*)

write(*,*)"Polak-Ribiere+ conjugate gradient"
call random_number(x)
call CGPR_verbose(f, f_fd, x, 10, &
     1d-3, 100, 1d-12, 1d-12)
write(*,*)norm2(x)
write(*,*)

write(*,*)"Newton-Raphson"
call random_number(x)
call NewtonRaphson(f, f_fd, fdd, x, 10, &
     100, 1d-12, 1d-12)
write(*,*)norm2(x)
write(*,*)

write(*,*)"BFGS"
call random_number(x)
call BFGS(f, f_fd, fdd, x, 10, &
     10, 100, 1d-12, 1d-12)
write(*,*)norm2(x)
write(*,*)

write(*,*)"Trust region"
call random_number(x)
call trust_region(fd_residue, fdd_Jacobian, x, 10, 10, &
     100, 100, 1d-12, 1d-12)
write(*,*)norm2(x)
write(*,*)

write(*,*)"Trust region verbose"
call random_number(x)
call trust_region_verbose(fd_residue, fdd_Jacobian, x, 10, 10, &
     100, 100, 1d-12, 1d-12)
write(*,*)norm2(x)
write(*,*)

write(*,*)"Gauss-BFGS"
call random_number(x)
call Gauss_BFGS(fd_residue, fdd_Jacobian, x, 10, 10, &
                1000, 1d-12, 1d-12)
write(*,*)norm2(x)
write(*,*)

write(*,*)"Augmented Lagrangian based on Newton-Raphson"
call random_number(x)
lambda0 = 0d0
call ALagrangian_NewtonRaphson( &
    f, f_fd, fdd, c, c_cd, c_cd_cdd, x, 10, 1, &
    lambda0, 1d0, &
    100, 10, 100, 1d-12, 1d-12)
write(*,*)norm2(x) - 1d0
write(*,*)

write(*,*)"Augmented Lagrangian based on BFGS"
call random_number(x)
lambda0 = 0d0
call ALagrangian_BFGS( &
    f, f_fd, fdd, c, c_cd, c_cd_cdd, x, 10, 1, &
    lambda0, 1d0, &
    100, 10, 100, 1d-12, 1d-12)
write(*,*)norm2(x) - 1d0
write(*,*)

contains
subroutine f(fx,x,dim)
    real*8,intent(out)::fx
    integer,intent(in)::dim
    real*8,dimension(dim),intent(in)::x
    integer::i
    fx=0d0
    do i=1,dim
        fx=fx+x(i)**4
    end do
end subroutine f
    
subroutine fd(fdx,x,dim)
    integer,intent(in)::dim
    real*8,dimension(dim),intent(out)::fdx
    real*8,dimension(dim),intent(in)::x
    integer::i
    forall(i=1:dim)
        fdx(i)=4d0*x(i)**3
    end forall
end subroutine fd

subroutine f_fd(fx,fdx,x,dim)
    integer,intent(in)::dim
    real*8,intent(out)::fx
    real*8,dimension(dim),intent(out)::fdx
    real*8,dimension(dim),intent(in)::x
    integer::i
    fx=0d0
    do i=1,dim
        fx=fx+x(i)**4
        fdx(i)=4d0*x(i)**3
    end do
end subroutine f_fd

subroutine fdd(fddx,x,dim)
    integer,intent(in)::dim
    real*8,dimension(dim,dim),intent(out)::fddx
    real*8,dimension(dim),intent(in)::x
    integer::i
    fddx=0d0
    forall(i=1:dim)
        fddx(i,i)=12d0*x(i)*x(i)
    end forall
end subroutine fdd

subroutine fd_residue(fdx,x,M,N)
    integer,intent(in)::M,N
    real*8,dimension(M),intent(out)::fdx
    real*8,dimension(N),intent(in)::x
    integer::i
    forall(i=1:M)
        fdx(i)=4d0*x(i)**3
    end forall
end subroutine fd_residue

subroutine fdd_Jacobian(fddx,x,M,N)
    integer,intent(in)::M,N
    real*8,dimension(M,N),intent(out)::fddx
    real*8,dimension(N),intent(in)::x
    integer::i
    fddx=0d0
    forall(i=1:M)
        fddx(i,i)=12d0*x(i)*x(i)
    end forall
end subroutine fdd_Jacobian

subroutine c(cx, x, M, N)
    integer, intent(in)::M, N
    real*8, dimension(M), intent(out)::cx
    real*8, dimension(N), intent(in)::x
    cx(1) = dot_product(x, x) - 1d0
end subroutine c

subroutine c_cd(cx, cdx, x, M, N)
    integer, intent(in)::M, N
    real*8, dimension(M), intent(out)::cx
    real*8, dimension(N, M), intent(out)::cdx
    real*8, dimension(N), intent(in)::x
    cx(1) = dot_product(x, x) - 1d0
    cdx(:, 1) = 2d0 * x
end subroutine c_cd

subroutine c_cd_cdd(cx, cdx, cddx, x, M, N)
    integer, intent(in)::M, N
    real*8, dimension(M), intent(out)::cx
    real*8, dimension(N, M), intent(out)::cdx
    real*8, dimension(N, N, M),intent(out)::cddx
    real*8, dimension(N), intent(in)::x
    integer::i
    cx(1) = dot_product(x, x) - 1d0
    cdx(:, 1) = 2d0 * x
    cddx = 0d0
    forall (i = 1 : N)
        cddx(i, i, 1) = 2d0
    end forall
end subroutine c_cd_cdd

end program main