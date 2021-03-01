program main

implicit none

integer::dim = 10, M = 10, N = 10
real*8,dimension(10)::x

call random_seed()
    
write(*,*)"BFGS"
call random_number(x)
call BFGS(f, f_fd, fdd, x, dim, &
          10, 100, 1d-15, 1d-15)
write(*,*)norm2(x)
write(*,*)

write(*,*)"Trust region"
call random_number(x)
call trust_region(fd_tr, fdd_tr, x, M, N, &
                  100, 100, 1d-15, 1d-15)
write(*,*)norm2(x)
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

subroutine fdd(fddx,x,N)
    integer,intent(in)::N
    real*8,dimension(N,N),intent(out)::fddx
    real*8,dimension(N),intent(in)::x
    integer::i
    fddx=0d0
    forall(i=1:dim)
        fddx(i,i)=12d0*x(i)*x(i)
    end forall
end subroutine fdd

subroutine fd_tr(fdx,x,M,N)
    integer,intent(in)::M,N
    real*8,dimension(M),intent(out)::fdx
    real*8,dimension(N),intent(in)::x
    integer::i
    forall(i=1:M)
        fdx(i)=4d0*x(i)**3
    end forall
end subroutine fd_tr

subroutine fdd_tr(fddx,x,M,N)
    integer,intent(in)::M,N
    real*8,dimension(M,N),intent(out)::fddx
    real*8,dimension(N),intent(in)::x
    integer::i
    fddx=0d0
    forall(i=1:M)
        fddx(i,i)=12d0*x(i)*x(i)
    end forall
end subroutine fdd_tr

subroutine constraint(cx,x,M,N)
    integer,intent(in)::M,N
    real*8,dimension(N),intent(in)::x
    real*8,dimension(M),intent(out)::cx
    cx(1)=dot_product(x,x)-1d0
end subroutine constraint

subroutine constraintd(cdx,x,M,N)
    integer,intent(in)::M,N
    real*8,dimension(N),intent(in)::x
    real*8,dimension(N,M),intent(out)::cdx
    cdx(:,1)=2d0*x
end subroutine constraintd

subroutine constraintdd(cddx,x,M,N)
    integer,intent(in)::M,N
    real*8,dimension(N),intent(in)::x
    real*8,dimension(N,N,M),intent(out)::cddx
    integer::i
    cddx = 0d0
    forall (i = 1 : N); cddx(i, i, 1) = 2d0; end forall
end subroutine constraintdd

end program main