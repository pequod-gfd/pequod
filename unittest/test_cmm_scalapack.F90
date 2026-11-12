#include "config.h"
#ifdef PETSC
#include "petscconf.h"
#endif

program test_cmm_scalapack
#if defined(PETSC) && defined(PETSC_HAVE_SCALAPACK)
#include "petsc/finclude/petscksp.h"
  use petscksp
  use qg_parallel, only: parallel_init, parallel_finalise, pabort
  implicit none

  Mat :: matrix
  Vec :: rhs,solution,exact
  KSP :: solver
  PC :: pc
  PetscInt, parameter :: n=96
  PetscInt :: lo,hi,row
  PetscErrorCode :: ierr
  PetscReal :: error
  PetscBool :: available

  call parallel_init()
  call MatCreate(PETSC_COMM_WORLD,matrix,ierr)
  call check(ierr)
  call MatSetSizes(matrix,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr)
  call check(ierr)
  call MatSetType(matrix,MATSCALAPACK,ierr)
  call check(ierr)
  call MatSetUp(matrix,ierr)
  call check(ierr)
  call MatGetFactorAvailable(matrix,MATSOLVERSCALAPACK,MAT_FACTOR_CHOLESKY,available,ierr)
  call check(ierr)
  if (.not. available) call pabort('installed ScaLAPACK Cholesky factor unavailable',__FILE__,__LINE__)
  call MatGetOwnershipRange(matrix,lo,hi,ierr)
  call check(ierr)
  do row = lo,hi-1
    call MatSetValue(matrix,row,row,4.0D0,INSERT_VALUES,ierr)
    call check(ierr)
    if (row > 0) then
      call MatSetValue(matrix,row,row-1,-1.0D0,INSERT_VALUES,ierr)
      call check(ierr)
    end if
    if (row < n-1) then
      call MatSetValue(matrix,row,row+1,-1.0D0,INSERT_VALUES,ierr)
      call check(ierr)
    end if
  end do
  call MatAssemblyBegin(matrix,MAT_FINAL_ASSEMBLY,ierr)
  call check(ierr)
  call MatAssemblyEnd(matrix,MAT_FINAL_ASSEMBLY,ierr)
  call check(ierr)
  call MatCreateVecs(matrix,solution,rhs,ierr)
  call check(ierr)
  call VecDuplicate(solution,exact,ierr)
  call check(ierr)
  call KSPCreate(PETSC_COMM_WORLD,solver,ierr)
  call check(ierr)
  call KSPSetOperators(solver,matrix,matrix,ierr)
  call check(ierr)
  call KSPSetType(solver,KSPPREONLY,ierr)
  call check(ierr)
  call KSPGetPC(solver,pc,ierr)
  call check(ierr)
  call PCSetType(pc,PCCHOLESKY,ierr)
  call check(ierr)
  call PCFactorSetMatSolverType(pc,MATSOLVERSCALAPACK,ierr)
  call check(ierr)
  call KSPSetUp(solver,ierr)
  call check(ierr)

  call VecSet(exact,1.25D0,ierr)
  call check(ierr)
  call MatMult(matrix,exact,rhs,ierr)
  call check(ierr)
  call KSPSolve(solver,rhs,solution,ierr)
  call check(ierr)
  call VecAXPY(solution,-1.0D0,exact,ierr)
  call check(ierr)
  call VecNorm(solution,NORM_INFINITY,error,ierr)
  call check(ierr)
  if (error>5e-13) call pabort('distributed Cholesky disagrees with exact solution',__FILE__,__LINE__)
  call VecScale(exact,-0.375D0,ierr)
  call check(ierr)
  call MatMult(matrix,exact,rhs,ierr)
  call check(ierr)
  call KSPSolve(solver,rhs,solution,ierr)
  call check(ierr)
  call VecAXPY(solution,-1.0D0,exact,ierr)
  call check(ierr)
  call VecNorm(solution,NORM_INFINITY,error,ierr)
  call check(ierr)
  if (error>5e-13) call pabort('reused distributed Cholesky disagrees with exact solution',__FILE__,__LINE__)
  print '(a,es14.6,a)', 'ScaLAPACK repeated Cholesky solve error: ',error,' PASS'

  call KSPDestroy(solver,ierr)
  call check(ierr)
  call VecDestroy(exact,ierr)
  call check(ierr)
  call VecDestroy(solution,ierr)
  call check(ierr)
  call VecDestroy(rhs,ierr)
  call check(ierr)
  call MatDestroy(matrix,ierr)
  call check(ierr)
  call parallel_finalise()

contains

  subroutine check(status)
    PetscErrorCode, intent(in) :: status

    if (status/=0) call pabort('PETSc error in ScaLAPACK capability fixture',__FILE__,__LINE__)
  end subroutine check
#else
  implicit none
  print '(a)', 'SKIP ScaLAPACK Cholesky unavailable in this build'
#endif
end program test_cmm_scalapack
