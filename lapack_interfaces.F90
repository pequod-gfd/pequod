#include "config.h"

! Simple interfaces to LAPACK.
module lapack_interfaces

  use qg_parallel
  use qg_parameters, only : complex_8, real_8

  implicit none
  
  private
  
  public :: eigendecompose, eigendecompose_symmetric, eigenrecompose, &
    & eigenrecompose_symmetric, determinant, invert, inverse, solve
    
  ! API documentation for LAPACK: http://www.netlib.org/lapack/double/
    
  interface dgeev
    subroutine dgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, &
    & work, lwork, info)
      use qg_parameters, only : real_8
      implicit none
      integer, intent(in) :: n
      integer, intent(in) :: lda
      integer, intent(in) :: ldvl
      integer, intent(in) :: ldvr
      integer, intent(in) :: lwork
      
      character(len = 1), intent(in) :: jobvl
      character(len = 1), intent(in) :: jobvr
      real(kind = real_8), dimension(lda, n), intent(inout) :: a
      real(kind = real_8), dimension(n), intent(out) :: wr
      real(kind = real_8), dimension(n), intent(out) :: wi
      real(kind = real_8), dimension(ldvl, n), intent(out) :: vl
      real(kind = real_8), dimension(ldvr, n), intent(out) :: vr
      real(kind = real_8), dimension(lwork), intent(out) :: work
      integer, intent(out) :: info
    end subroutine dgeev
  end interface dgeev
  
  interface dgesv
    subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
      use qg_parameters, only : real_8
      implicit none
      integer, intent(in) :: n
      integer, intent(in) :: nrhs
      integer, intent(in) :: lda
      integer, intent(in) :: ldb
      
      real(kind = real_8), dimension(lda, n), intent(inout) :: a
      integer, dimension(n), intent(out) :: ipiv
      real(kind = real_8), dimension(ldb, nrhs), intent(inout) :: b
      integer, intent(out) :: info
    end subroutine dgesv
  end interface dgesv
  
  interface dgetrf
    subroutine dgetrf(m, n, a, lda, ipiv, info)
      use qg_parameters, only : real_8
      implicit none
      integer, intent(in) :: m
      integer, intent(in) :: n
      integer, intent(in) :: lda
      
      real(kind = real_8), dimension(lda, n), intent(inout) :: a
      integer, dimension(min(m, n)), intent(out) :: ipiv
      integer, intent(out) :: info
    end subroutine dgetrf
  end interface dgetrf
  
  interface dgetri
    subroutine dgetri(n, a, lda, ipiv, work, lwork, info)
      use qg_parameters, only : real_8
      implicit none
      integer, intent(in) :: n
      integer, intent(in) :: lda
      integer, intent(in) :: lwork
      
      real(kind = real_8), dimension(lda, n), intent(inout) :: a
      integer, dimension(n), intent(in) :: ipiv
      real(kind = real_8), dimension(lwork), intent(out) :: work
      integer, intent(out) :: info
    end subroutine dgetri
  end interface dgetri
  
  interface dspev
    subroutine dspev(jobz, uplo, n, ap, w, z, ldz, work, info)
      use qg_parameters, only : real_8
      implicit none
      integer, intent(in) :: n
      integer, intent(in) :: ldz
      
      character(len = 1), intent(in) :: jobz
      character(len = 1), intent(in) :: uplo
      real(kind = real_8), dimension(n * (n + 1) / 2), intent(inout) :: ap
      real(kind = real_8), dimension(n), intent(out) :: w
      real(kind = real_8), dimension(ldz, n), intent(out) :: z
      real(kind = real_8), dimension(3 * n), intent(out) :: work
      integer, intent(out) :: info
    end subroutine dspev
  end interface dspev
  
  interface eigendecompose
    module procedure eigendecompose_complex, eigendecompose_real
  end interface eigendecompose
  
  interface eigenrecompose
    module procedure eigenrecompose_complex, eigenrecompose_real
  end interface eigenrecompose

contains

  ! Eigendecompose a general real square matrix.
  subroutine eigendecompose_complex(matrix, eigenvalues, eigenvectors)
    real(kind = real_8), dimension(:, :), intent(in) :: matrix
    complex(kind = complex_8), dimension(size(matrix, 1)), intent(out) :: eigenvalues
    real(kind = real_8), dimension(size(matrix, 1), size(matrix, 1)), intent(out) :: eigenvectors
    
    integer :: i
    
    integer :: n
    integer :: lda
    integer :: ldvl
    integer :: ldvr
    integer :: lwork
    
    character(len = 1), save :: jobvl = "N"
    character(len = 1), save :: jobvr = "V"
    real(kind = real_8), dimension(size(matrix, 1), size(matrix, 1)) :: a
    real(kind = real_8), dimension(size(matrix, 1)) :: wr
    real(kind = real_8), dimension(size(matrix, 1)) :: wi
    real(kind = real_8), dimension(1, size(matrix, 1)) :: vl
    !real(kind = real_8), dimension(size(matrix, 1), size(matrix, 1)) :: vr
    real(kind = real_8), dimension(4 * size(matrix, 1)) :: work
    integer :: info
    
    assert(size(matrix, 1) == size(matrix, 2))
    assert(size(matrix, 1) > 0)
    
    n = size(matrix, 1)
    lda = size(matrix, 1)
    ldvl = 1
    ldvr = size(matrix, 1)
    lwork = size(work)
    
    a = matrix
    
    call dgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, eigenvectors, ldvr, &
      & work, lwork, info)
    assert(info == 0)
    
    do i = 1, size(matrix, 1)
      eigenvalues(i) = cmplx(wr(i), wi(i), kind = complex_8)
    end do
    
  end subroutine eigendecompose_complex
  
  ! Eigendecompose a general real square matrix, assuming real eigenvalues.
  subroutine eigendecompose_real(matrix, eigenvalues, eigenvectors)
    real(kind = real_8), dimension(:, :), intent(in) :: matrix
    real(kind = real_8), dimension(size(matrix, 1)), intent(out) :: eigenvalues
    real(kind = real_8), dimension(size(matrix, 1), size(matrix, 1)), intent(out) :: eigenvectors
        
    integer :: n
    integer :: lda
    integer :: ldvl
    integer :: ldvr
    integer :: lwork
    
    character(len = 1), save :: jobvl = "N"
    character(len = 1), save :: jobvr = "V"
    real(kind = real_8), dimension(size(matrix, 1), size(matrix, 1)) :: a
    !real(kind = real_8), dimension(size(matrix, 1)) :: wr
    real(kind = real_8), dimension(size(matrix, 1)) :: wi
    real(kind = real_8), dimension(1, size(matrix, 1)) :: vl
    !real(kind = real_8), dimension(size(matrix, 1), size(matrix, 1)) :: vr
    real(kind = real_8), dimension(4 * size(matrix, 1)) :: work
    integer :: info
    
    assert(size(matrix, 1) == size(matrix, 2))
    assert(size(matrix, 1) > 0)
    
    n = size(matrix, 1)
    lda = size(matrix, 1)
    ldvl = 1
    ldvr = size(matrix, 1)
    lwork = size(work)
    
    a = matrix
    
    call dgeev(jobvl, jobvr, n, a, lda, eigenvalues, wi, vl, ldvl, eigenvectors, ldvr, &
      & work, lwork, info)
    assert(info == 0)
    
  end subroutine eigendecompose_real
  
  ! Eigendecompose a square symmetric matrix.
  subroutine eigendecompose_symmetric(matrix, eigenvalues, eigenvectors)
    real(kind = real_8), dimension(:, :), intent(in) :: matrix
    real(kind = real_8), dimension(size(matrix, 1)), intent(out) :: eigenvalues
    real(kind = real_8), dimension(size(matrix, 1), size(matrix, 1)), intent(out) :: eigenvectors
        
    integer :: i, j
    
    integer :: n
    integer :: ldz
    
    character(len = 1), save :: jobz = "V"
    character(len = 1), save :: uplo = "U"
    real(kind = real_8), dimension(size(matrix, 1) * (size(matrix, 1) + 1) / 2) :: ap
    real(kind = real_8), dimension(3 * size(matrix, 1)) :: work
    integer :: info
    
    assert(size(matrix, 1) == size(matrix, 2))
    assert(size(matrix, 1) > 0)
    
    n = size(matrix, 1)
    ldz = size(matrix, 1)
    
    do j = 1, size(matrix, 1)
      do i = 1, j
        ap(i + (j - 1) * j / 2) = matrix(i, j)
      end do
    end do
    
    call dspev(jobz, uplo, n, ap, eigenvalues, eigenvectors, ldz, work, info)
    assert(info == 0)
    
  end subroutine eigendecompose_symmetric
  
  ! Recompose a general square matrix from its eigenvalues and
  ! eigenvectors.
  subroutine eigenrecompose_complex(matrix, eigenvalues, eigenvectors)
    real(kind = real_8), dimension(:, :), intent(out) :: matrix
    complex(kind = complex_8), dimension(size(matrix, 1)), intent(in) :: eigenvalues
    real(kind = real_8), dimension(size(matrix, 1), size(matrix, 1)), intent(in) :: eigenvectors
        
    assert(size(matrix, 1) == size(matrix, 2))
    
    matrix = real(matmul(eigenvectors * spread(eigenvalues, 1, size(matrix, 1)), inverse(eigenvectors)), kind = real_8)
    
  end subroutine eigenrecompose_complex
  
  ! Recompose a general real square matrix from its eigenvalues and
  ! eigenvectors.
  subroutine eigenrecompose_real(matrix, eigenvalues, eigenvectors)
    real(kind = real_8), dimension(:, :), intent(out) :: matrix
    real(kind = real_8), dimension(size(matrix, 1)), intent(in) :: eigenvalues
    real(kind = real_8), dimension(size(matrix, 1), size(matrix, 1)), intent(in) :: eigenvectors
        
    assert(size(matrix, 1) == size(matrix, 2))
    
    matrix = matmul(eigenvectors * spread(eigenvalues, 1, size(matrix, 1)), inverse(eigenvectors))
  
  end subroutine eigenrecompose_real
  
  ! Recompose a symmetric square matrix from its eigenvalues and eigenvectors.
  subroutine eigenrecompose_symmetric(matrix, eigenvalues, eigenvectors)
    real(kind = real_8), dimension(:, :), intent(out) :: matrix
    real(kind = real_8), dimension(size(matrix, 1)), intent(in) :: eigenvalues
    real(kind = real_8), dimension(size(matrix, 1), size(matrix, 1)), intent(in) :: eigenvectors
        
    assert(size(matrix, 1) == size(matrix, 2))
    
    matrix = matmul(eigenvectors * spread(eigenvalues, 1, size(matrix, 1)), transpose(eigenvectors))
  
  end subroutine eigenrecompose_symmetric
    
  ! Compute the determinant of a square matrix.
  function determinant(matrix)
    real(kind = real_8), dimension(:, :), intent(in) :: matrix
    
    real(kind = real_8) :: determinant
    
    integer :: i
    complex(kind = complex_8) :: determinant_c
    
    integer :: n
    integer :: lda
    integer :: ldvl
    integer :: ldvr
    integer :: lwork
    
    character(len = 1), save :: jobvl = "N"
    character(len = 1), save :: jobvr = "N"
    real(kind = real_8), dimension(size(matrix, 1), size(matrix, 1)) :: a
    real(kind = real_8), dimension(size(matrix, 1)) :: wr
    real(kind = real_8), dimension(size(matrix, 1)) :: wi
    real(kind = real_8), dimension(1, size(matrix, 1)) :: vl
    real(kind = real_8), dimension(1, size(matrix, 1)) :: vr
    real(kind = real_8), dimension(4 * size(matrix, 1)) :: work
    integer :: info
    
    assert(size(matrix, 1) == size(matrix, 2))
    assert(size(matrix, 1) > 0)
    
    n = size(matrix, 1)
    lda = size(matrix, 1)
    ldvl = 1
    ldvr = 1
    lwork = size(work)
    
    a = matrix
    
    call dgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, &
      & work, lwork, info)
    assert(info == 0)
    
    determinant_c = 1.0D0
    do i = 1, size(matrix, 1)
      determinant_c = determinant_c * cmplx(wr(i), wi(i), kind = complex_8)
    end do
    
    determinant = real(determinant_c, kind = real_8)
    
  end function determinant
    
  ! Invert a square matrix.
  subroutine invert(matrix)
    real(kind = real_8), dimension(:, :), intent(inout) :: matrix
    
    integer :: m
    integer :: n
    integer :: lda
    
    integer, dimension(size(matrix, 1)) :: ipiv
    integer :: info
    
    integer :: lwork
    
    real(kind = real_8), dimension(size(matrix, 1)) :: work
        
    assert(size(matrix, 1) == size(matrix, 2))
    assert(size(matrix, 1) > 0)
    
    m = size(matrix, 1)
    n = size(matrix, 1)
    lda = size(matrix, 1)
    
    call dgetrf(m, n, matrix, lda, ipiv, info)
    assert(info == 0)
    
    lwork = size(work)
    
    call dgetri(n, matrix, lda, ipiv, work, lwork, info)
    assert(info == 0)
        
  end subroutine invert
  
  ! Return the inverse of a square matrix.
  function inverse(matrix)
    real(kind = real_8), dimension(:, :), intent(in) :: matrix
    
    real(kind = real_8), dimension(size(matrix, 1), size(matrix, 1)) :: inverse
    
    inverse = matrix
    call invert(inverse)
    
  end function inverse
  
  ! Solve a system of linear equations.
  subroutine solve(matrix, rhs, soln)
    real(kind = real_8), dimension(:, :), intent(in) :: matrix
    real(kind = real_8), dimension(size(matrix, 1)), intent(in) :: rhs
    real(kind = real_8), dimension(size(matrix, 1)), intent(out) :: soln
    
    integer :: info, lda, ldb, n, nrhs
    integer, dimension(size(matrix, 1)) :: ipiv
    real(kind = real_8), dimension(size(matrix, 1), 1) :: x
    real(kind = real_8), dimension(size(matrix, 1), size(matrix, 1)) :: a
    
    assert(size(matrix, 1) == size(matrix, 2))
    assert(size(matrix, 1) > 0)
    
    n = size(matrix, 1)
    nrhs = 1
    a = matrix
    lda = size(matrix, 1)
    x(:, 1) = rhs
    ldb = size(matrix, 1)
    
    call dgesv(n, nrhs, a, lda, ipiv, x, ldb, info)
    assert(info == 0)
    
    soln = x(:, 1)
    
  end subroutine solve

end module lapack_interfaces
