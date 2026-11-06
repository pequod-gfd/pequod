#include "config.h"

module test_tools

  use qg_parameters, only : complex_8, real_8
  
  implicit none
  
  include "mpif.h"
  
  private
  
  public :: report_test, preport_test, fe, fne
  
  interface fe
    module procedure fe_rank0, fe_rank0_tol, fe_rank1, fe_rank1_tol
  end interface fe
  
  interface fne
    module procedure fne_rank0, fne_rank0_tol, fne_rank1, fne_rank1_tol, &
      & fne_rank2, fne_rank2_tol, fne_rank3, fne_rank3_tol, fne_rank0_complex, &
      & fne_rank1_complex
  end interface fne
  
  ! Tolerance factor for comparison of floating point values (when no overriding
  ! tolerance is supplied).
  real(kind = real_8), parameter :: fne_factor = 10.0D0
  
contains

  ! Print the result of a test to standard output.
  subroutine report_test(name, fail)
    character(len = *), intent(in) :: name
    logical, intent(in) :: fail
    
#ifdef DDEBUG
    integer :: ierr, init, rank

    call mpi_initialized(init, ierr)
    if(ierr /= MPI_SUCCESS) then
      write(0, "(a)") "Failed to determine MPI initialisation status"
      stop 1
    else if(init == 1) then
      call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
      if(ierr /= MPI_SUCCESS) then
        write(0, "(a)") "Failed to determine rank"
        call mpi_abort(MPI_COMM_WORLD, MPI_ERR_OTHER, ierr)
        stop 1
      else if(rank /= 0) then
        write(0, "(a)") "report_test called by process with rank > 0"
        call mpi_abort(MPI_COMM_WORLD, MPI_ERR_OTHER, ierr)
        stop 1
      end if
    end if
#endif
    
    if(fail) then
      write(6, "(a)") "[" // trim(name) // "] : FAIL"
    else
      write(6, "(a)") "[" // trim(name) // "] : PASS"
    end if
  
  end subroutine report_test
  
  ! Print the result of a test to standard output, performing an allor for the
  ! value of fail.
  subroutine preport_test(name, fail)
    character(len = *), intent(in) :: name
    logical, intent(in) :: fail
    
    integer :: p, ierr, lfail, gfail
    
    if(fail) then
      lfail = 1
    else
      lfail = 0
    end if    
    call mpi_allreduce(lfail, gfail, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    
    call mpi_comm_rank(MPI_COMM_WORLD, p, ierr)
    if(p == 0) call report_test(name, gfail == 1)
    
  end subroutine preport_test
  
  ! Return whether real variables are equal.
  pure function fe_rank0(arr1, arr2) result(fe)
    real(kind = real_8), intent(in) :: arr1
    real(kind = real_8), intent(in) :: arr2
    
    logical :: fe

    fe = dabs(arr1 - arr2) <= fne_factor * max(spacing(arr1), epsilon(arr1))
    
  end function fe_rank0   
  
  ! Return whether real variables are equal with a specified tolerance.
  pure function fe_rank0_tol(arr1, arr2, tol) result(fe)
    real(kind = real_8), intent(in) :: arr1
    real(kind = real_8), intent(in) :: arr2
    real(kind = real_8), intent(in) :: tol
    
    logical :: fe
    
    fe = dabs(arr1 - arr2) <= tol
    
  end function fe_rank0_tol
  
  ! Return whether real variables are equal.
  pure function fe_rank1(arr1, arr2)
    real(kind = real_8), dimension(:), intent(in) :: arr1
    real(kind = real_8), dimension(size(arr1)), intent(in) :: arr2
    
    logical :: fe_rank1
    
    integer :: i
    
    fe_rank1 = .true.
    do i = 1, size(arr1)
      if(.not. fe(arr1(i), arr2(i))) then
        fe_rank1 = .false.
        exit
      end if
    end do
    
  end function fe_rank1  
  
  ! Return whether real variables are equal with a specified tolerance.
  pure function fe_rank1_tol(arr1, arr2, tol)
    real(kind = real_8), dimension(:), intent(in) :: arr1
    real(kind = real_8), dimension(size(arr1)), intent(in) :: arr2
    real(kind = real_8), intent(in) :: tol
    
    logical :: fe_rank1_tol
    
    integer :: i
    
    fe_rank1_tol = .true.
    do i = 1, size(arr1)
      if(.not. fe(arr1(i), arr2(i), tol = tol)) then
        fe_rank1_tol = .false.
        exit
      end if
    end do
    
  end function fe_rank1_tol
  
  ! Return whether real variables are not equal.
  pure function fne_rank0(arr1, arr2) result(fne)
    real(kind = real_8), intent(in) :: arr1
    real(kind = real_8), intent(in) :: arr2
    
    logical :: fne
    
    fne = dabs(arr1 - arr2) > fne_factor * max(spacing(arr1), epsilon(arr1))
#ifdef DDEBUG
    if(arr1 /= arr1 .or. arr2 /= arr2) fne = .true.
#endif
    
  end function fne_rank0
  
  ! Return whether real variables are not equal with a specified tolerance.
  pure function fne_rank0_tol(arr1, arr2, tol) result(fne)
    real(kind = real_8), intent(in) :: arr1
    real(kind = real_8), intent(in) :: arr2
    real(kind = real_8), intent(in) :: tol
    
    logical :: fne
    
    fne = dabs(arr1 - arr2) > tol
#ifdef DDEBUG
    if(arr1 /= arr1 .or. arr2 /= arr2) fne = .true.
#endif
    
  end function fne_rank0_tol
  
  ! Return whether real variables are not equal.
  pure function fne_rank1(arr1, arr2)
    real(kind = real_8), dimension(:), intent(in) :: arr1
    real(kind = real_8), dimension(size(arr1)), intent(in) :: arr2
    
    logical :: fne_rank1
    
    integer :: i
    
    fne_rank1 = .false.
    do i = 1, size(arr1)
      if(fne(arr1(i), arr2(i))) then
        fne_rank1 = .true.
        exit
      end if
    end do
    
  end function fne_rank1
  
  ! Return whether real variables are not equal with a specified tolerance.
  pure function fne_rank1_tol(arr1, arr2, tol) result(fne_rank1)
    real(kind = real_8), dimension(:), intent(in) :: arr1
    real(kind = real_8), dimension(size(arr1)), intent(in) :: arr2
    real(kind = real_8), intent(in) :: tol
    
    logical :: fne_rank1
    
    integer :: i
    
    fne_rank1 = .false.
    do i = 1, size(arr1)
      if(fne(arr1(i), arr2(i), tol = tol)) then
        fne_rank1 = .true.
        exit
      end if
    end do
    
  end function fne_rank1_tol
  
  ! Return whether real variables are not equal.
  pure function fne_rank2(arr1, arr2)
    real(kind = real_8), dimension(:, :), intent(in) :: arr1
    real(kind = real_8), dimension(size(arr1, 1), size(arr1, 2)), intent(in) :: arr2
    
    logical :: fne_rank2
    
    integer :: i, j
    
    fne_rank2 = .false.
    i_loop: do i = 1, size(arr1, 1)
      do j = 1, size(arr1, 2)
        if(fne(arr1(i, j), arr2(i, j))) then
          fne_rank2 = .true.
          exit i_loop
        end if
      end do
    end do i_loop
    
  end function fne_rank2
  
  ! Return whether real variables are not equal with a specified tolerance.
  pure function fne_rank2_tol(arr1, arr2, tol) result(fne_rank2)
    real(kind = real_8), dimension(:, :), intent(in) :: arr1
    real(kind = real_8), dimension(size(arr1, 1), size(arr1, 2)), intent(in) :: arr2
    real(kind = real_8), intent(in) :: tol
    
    logical :: fne_rank2
    
    integer :: i, j
    
    fne_rank2 = .false.
    i_loop: do i = 1, size(arr1, 1)
      do j = 1, size(arr1, 2)
        if(fne(arr1(i, j), arr2(i, j), tol = tol)) then
          fne_rank2 = .true.
          exit i_loop
        end if
      end do
    end do i_loop
    
  end function fne_rank2_tol
  
  ! Return whether real variables are not equal.
  pure function fne_rank3(arr1, arr2)
    real(kind = real_8), dimension(:, :, :), intent(in) :: arr1
    real(kind = real_8), dimension(size(arr1, 1), size(arr1, 2), size(arr1, 3)), intent(in) :: arr2
    
    logical :: fne_rank3
    
    integer :: i, j, k
    
    fne_rank3 = .false.
    i_loop: do i = 1, size(arr1, 1)
      do j = 1, size(arr1, 2)
        do k = 1, size(arr1, 3)
          if(fne(arr1(i, j, k), arr2(i, j, k))) then
            fne_rank3 = .true.
            exit i_loop
          end if
        end do
      end do
    end do i_loop
    
  end function fne_rank3
  
  ! Return whether real variables are not equal with a specified tolerance.
  pure function fne_rank3_tol(arr1, arr2, tol) result(fne_rank3)
    real(kind = real_8), dimension(:, :, :), intent(in) :: arr1
    real(kind = real_8), dimension(size(arr1, 1), size(arr1, 2), size(arr1, 3)), intent(in) :: arr2
    real(kind = real_8), intent(in) :: tol
    
    logical :: fne_rank3
    
    integer :: i, j, k
    
    fne_rank3 = .false.
    i_loop: do i = 1, size(arr1, 1)
      do j = 1, size(arr1, 2)
        do k = 1, size(arr1, 3)
          if(fne(arr1(i, j, k), arr2(i, j, k), tol = tol)) then
            fne_rank3 = .true.
            exit i_loop
          end if
        end do
      end do
    end do i_loop
    
  end function fne_rank3_tol
  
  ! Return whether complex variables are not equal.
  pure function fne_rank0_complex(arr1, arr2) result(fne)
    complex(kind = complex_8), intent(in) :: arr1
    complex(kind = complex_8), intent(in) :: arr2
    
    logical :: fne
    
    fne = dabs(real(arr1, kind = real_8) - real(arr2, kind = real_8))   > fne_factor * max(spacing(real(arr1, kind = real_8)), epsilon(real(arr1, kind = real_8))) &
   & .or. dabs(aimag(arr1) - aimag(arr2)) > fne_factor * max(spacing(aimag(arr1)), epsilon(real(arr1, kind = real_8)))
#ifdef DDEBUG
    if(arr1 /= arr1 .or. arr2 /= arr2) fne = .true.
#endif
    
  end function fne_rank0_complex
  
  ! Return whether complex variables are not equal.
  pure function fne_rank1_complex(arr1, arr2)
    complex(kind = complex_8), dimension(:), intent(in) :: arr1
    complex(kind = complex_8), dimension(size(arr1)), intent(in) :: arr2
    
    logical :: fne_rank1_complex
    
    integer :: i
    
    fne_rank1_complex = .false.
    do i = 1, size(arr1)
      if(fne(arr1(i), arr2(i))) then
        fne_rank1_complex = .true.
        exit
      end if
    end do
    
  end function fne_rank1_complex

end module test_tools
