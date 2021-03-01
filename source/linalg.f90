!Naming convention (following LAPACK):
!    d = real*8, z = complex*16
!    ge = general matrix
!    sy = real symmetric matrix
!    he = Hermitian matrix
!    po = real symmetric or Hermitian positive definite matrix
!Only use lower triangle of sy & he & po
module linalg
    implicit none

contains
!--------------- Vector ---------------
    !M dimensional vector a, N dimensional vector b, vector_direct_product(a,b) = a b
    function vector_direct_product(a, b, M, N)
        integer, intent(in)::M, N
        real*8, dimension(M), intent(in)::a
        real*8, dimension(N), intent(in)::b
        real*8, dimension(M, N)::vector_direct_product
        integer::i
        forall (i = 1 : N)
            vector_direct_product(:, i) = a * b(i)
        end forall
    end function vector_direct_product
!---------------- End -----------------

!--------------- Matrix ---------------
    !N order matrix A, strictly upper triangle is blank,
    !copy strictly lower triangle elements to strictly upper triangle
    subroutine dsyL2U(A, N)
        integer, intent(in)::N
        real*8, dimension(N, N), intent(inout)::A
        integer::i, j
        forall(i = 1 : N - 1, j = 2 : N, i < j)
            A(i, j) = A(j, i)
        end forall
    end subroutine dsyL2U
!---------------- End -----------------

!-------------- Inverse ---------------
    !Input : N x N positive definite matrix A
    !Output: A harvests its inverse
    !        If A is indeed po, info returns 0
    !        else, the info-th leading minor of A <= 0 so failed
    !A will be overwritten even if failed
    subroutine My_dpotri(A, N, info)
        integer, intent(in)::N
        real*8, dimension(N,N), intent(inout)::A
        integer, intent(out)::info
        call dpotrf('L', N, A, N, info)
        if (info /= 0) return
        call dpotri('L', N, A, N, info)
    end subroutine My_dpotri
!---------------- End -----------------

end module linalg