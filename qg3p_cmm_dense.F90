#include "config.h"
#ifdef PETSC
#include "petscconf.h"
#endif

! Optional distributed capacitance factors. PETSc owns block-cyclic storage.
module qg_cmm_dense
#if defined(PETSC) && defined(PETSC_HAVE_SCALAPACK)
#include "petsc/finclude/petscmat.h"
  use petscmat
  use qg_parallel, only: pabort
  use qg_parameters, only: real_8
#endif
  implicit none
  private
  public :: has_scalapack

#if defined(PETSC) && defined(PETSC_HAVE_SCALAPACK)
  logical, parameter :: has_scalapack = .true.
  public :: dense_factor, create_factor, set_factor_column, assemble_factor, factor_dense, solve_dense, destroy_factor

  type dense_factor
    private
    Mat :: matrix
    Vec :: rhs,solution,gathered
    VecScatter :: gather
    PetscInt :: lo,hi,physical_size,padding
  end type dense_factor

contains

  subroutine check_petsc(ierr)
    PetscErrorCode, intent(in) :: ierr

    if (ierr /= 0) call pabort('PETSc distributed CMM factor failed',__FILE__,__LINE__)
  end subroutine check_petsc

  subroutine create_factor(factor, comm, ranks, size_c)
    type(dense_factor), intent(out) :: factor
    integer, intent(in) :: comm,ranks,size_c
    PetscErrorCode :: ierr
    PetscInt :: n,block,row
    PetscBool :: supported

    ! PETSc's ScaLAPACK Mat/Vec bridge requires equal local row lengths.
    ! Append an uncoupled identity block to complete the last set of rows.
    factor%physical_size = size_c
    factor%padding = modulo(-size_c,ranks)
    n = size_c+factor%padding
    block = max(1,min(32,n/8))
    call MatCreateScaLAPACK(comm,block,block,n,n,0,0,factor%matrix,ierr)
    call check_petsc(ierr)
    call MatGetFactorAvailable(factor%matrix,MATSOLVERSCALAPACK,MAT_FACTOR_CHOLESKY,supported,ierr)
    call check_petsc(ierr)
    if (.not. supported) call pabort('ScaLAPACK Cholesky is unavailable',__FILE__,__LINE__)
    call MatCreateVecs(factor%matrix,factor%solution,factor%rhs,ierr)
    call check_petsc(ierr)
    call VecGetOwnershipRange(factor%rhs,factor%lo,factor%hi,ierr)
    call check_petsc(ierr)
    call VecScatterCreateToAll(factor%solution,factor%gather,factor%gathered,ierr)
    call check_petsc(ierr)
    do row = max(factor%lo,factor%physical_size),factor%hi-1
      call MatSetValue(factor%matrix,row,row,1.0D0,INSERT_VALUES,ierr)
      call check_petsc(ierr)
    end do
  end subroutine create_factor

  subroutine set_factor_column(factor, column, owned_rows, values)
    type(dense_factor), intent(inout) :: factor
    integer, intent(in) :: column,owned_rows(:) ! Fortran one-based global indices
    real(real_8), intent(in) :: values(:)
    PetscInt :: n,rows(size(owned_rows)),cols(1)
    PetscErrorCode :: ierr

    n = size(owned_rows)
    if (n == 0) return
    rows = owned_rows-1
    cols = [column-1]
    call MatSetValues(factor%matrix,n,rows,1,cols,values(owned_rows),INSERT_VALUES,ierr)
    call check_petsc(ierr)
  end subroutine set_factor_column

  subroutine assemble_factor(factor, last)
    type(dense_factor), intent(inout) :: factor
    logical, intent(in) :: last
    MatAssemblyType :: stage
    PetscErrorCode :: ierr

    stage = MAT_FLUSH_ASSEMBLY
    if (last) stage = MAT_FINAL_ASSEMBLY
    call MatAssemblyBegin(factor%matrix,stage,ierr)
    call check_petsc(ierr)
    call MatAssemblyEnd(factor%matrix,stage,ierr)
    call check_petsc(ierr)
  end subroutine assemble_factor

  ! Compare the assembled matrix with its transpose before Cholesky.
  ! The condition bound comes from the background rectangle eigenvalues.
  subroutine factor_dense(factor, condition_bound, assembly_error, asymmetry, bytes, error)
    type(dense_factor), intent(inout) :: factor
    real(real_8), intent(in) :: condition_bound,assembly_error
    real(real_8), intent(out) :: asymmetry
    integer(selected_int_kind(18)), intent(out) :: bytes
    character(*), intent(out) :: error
    PetscErrorCode :: ierr
    Mat :: transpose
    MatFactorInfo :: info
    MatInfo :: memory
    PetscReal :: norm,difference

    call MatNorm(factor%matrix,NORM_FROBENIUS,norm,ierr)
    call check_petsc(ierr)
    call MatTranspose(factor%matrix,MAT_INITIAL_MATRIX,transpose,ierr)
    call check_petsc(ierr)
    call MatAXPY(transpose,-1.0D0,factor%matrix,SAME_NONZERO_PATTERN,ierr)
    call check_petsc(ierr)
    call MatNorm(transpose,NORM_FROBENIUS,difference,ierr)
    call check_petsc(ierr)
    call MatDestroy(transpose,ierr)
    call check_petsc(ierr)
    norm = sqrt(max(tiny(norm),norm**2-real(factor%padding,real_8)))
    asymmetry = difference/norm
    error = ''
    if (asymmetry > 5.0D-12 .or. condition_bound*1.0D-8 < max(asymmetry,assembly_error,epsilon(norm))) then
      error = 'distributed capacitance assembly accuracy is insufficient for its condition bound'
      return
    end if
    call MatGetInfo(factor%matrix,MAT_LOCAL,memory,ierr)
    call check_petsc(ierr)
    bytes = int(memory%nz_allocated,kind(bytes))*(storage_size(0.0D0)/8)
    call MatFactorInfoInitialize(info,ierr)
    call check_petsc(ierr)
    call MatCholeskyFactor(factor%matrix,PETSC_NULL_IS,info,ierr)
    call check_petsc(ierr)
  end subroutine factor_dense

  subroutine solve_dense(factor, values)
    type(dense_factor), intent(inout) :: factor
    real(real_8), intent(inout) :: values(:)
    PetscInt :: count_rows,row,physical_end
    PetscErrorCode :: ierr
    real(real_8) :: local(factor%hi-factor%lo)

    count_rows = factor%hi-factor%lo
    if (count_rows > 0) then
      physical_end = min(factor%hi,factor%physical_size)
      local = 0.0D0
      if (physical_end > factor%lo) local(:physical_end-factor%lo) = values(factor%lo+1:physical_end)
      call VecSetValues(factor%rhs,count_rows,[(row,row=factor%lo,factor%hi-1)],local,INSERT_VALUES,ierr)
      call check_petsc(ierr)
    end if
    call VecAssemblyBegin(factor%rhs,ierr)
    call check_petsc(ierr)
    call VecAssemblyEnd(factor%rhs,ierr)
    call check_petsc(ierr)
    call MatSolve(factor%matrix,factor%rhs,factor%solution,ierr)
    call check_petsc(ierr)
    call VecScatterBegin(factor%gather,factor%solution,factor%gathered,INSERT_VALUES,SCATTER_FORWARD,ierr)
    call check_petsc(ierr)
    call VecScatterEnd(factor%gather,factor%solution,factor%gathered,INSERT_VALUES,SCATTER_FORWARD,ierr)
    call check_petsc(ierr)
    count_rows = size(values)
    call VecGetValues(factor%gathered,count_rows,[(row,row=0,count_rows-1)],values,ierr)
    call check_petsc(ierr)
  end subroutine solve_dense

  subroutine destroy_factor(factor)
    type(dense_factor), intent(inout) :: factor
    PetscErrorCode :: ierr

    call VecScatterDestroy(factor%gather,ierr)
    call check_petsc(ierr)
    call VecDestroy(factor%gathered,ierr)
    call check_petsc(ierr)
    call VecDestroy(factor%rhs,ierr)
    call check_petsc(ierr)
    call VecDestroy(factor%solution,ierr)
    call check_petsc(ierr)
    call MatDestroy(factor%matrix,ierr)
    call check_petsc(ierr)
  end subroutine destroy_factor

#else
  logical, parameter :: has_scalapack = .false.
#endif
end module qg_cmm_dense
