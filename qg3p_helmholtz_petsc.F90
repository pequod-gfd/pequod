#include "config.h"

#ifdef PETSC

! Helmholtz solver routines for the 3D QG model.
module qg_helmholtz_petsc

#include "petsc/finclude/petscksp.h"
  use petscksp

  use qg_laplace_bc
  use qg_parallel
  use qg_parameters, only : real_8, real_8_chr

  implicit none
  
#include "petscconf.h"

  private
  
  public :: solve_helmholtz
  public :: helmholtz_matrices
  public :: allocate, deallocate
  
  ! PETSc Helmholtz solver cache.
  type helmholtz_matrices
    ! Decomposition
    type(decomposition), pointer :: decomp  
    ! Helmholtz operator coefficeints
    real(kind = real_8), dimension(:), pointer :: ss
    
    ! LHS matrices
    Mat, dimension(:), pointer :: matrix
    ! RHS vector
    Vec :: rhs
    ! Global solution vector
    Vec :: soln
    
    ! Preconditioners
    PC, dimension(:), pointer :: pc
    ! KSP
    KSP, dimension(:), pointer :: ksp
    
    ! Global and local indices in PETSc numbering
    integer, dimension(:, :), pointer :: gix, lix
    
    ! Local solution vector
    Vec :: lsoln
    ! Global to local index sets
    IS :: gis, lis
    ! Global to local scatter
    VecScatter :: scatter
  end type helmholtz_matrices
  
  ! Default maximum number of iterations for iterative methods.
  integer, parameter :: default_max_its = 1000
  ! Default absolute tolerance for iterative methods.
  real(kind = real_8), parameter :: default_atol = 1.0D-16
  ! Default relative tolerance for iterative methods.
  real(kind = real_8), parameter :: default_rtol = 1.0D-7
  
  interface allocate
    module procedure allocate_helmholtz_matrices
  end interface allocate
  
  interface deallocate
    module procedure deallocate_helmholtz_matrices
  end interface deallocate
  
contains

  ! Allocate PETSc Helmholtz solver cache. This assembles the matrices, sets up
  ! the KSPs and preconditioners, and configures parallel communication objects.
  ! The equations solved are:
  !   \left( \partial_{xx} + \partial_{yy} - ss_i \right) \psi_i = \zeta_i
  ! subject to homogeneous Dirichlet boundary conditions.
  subroutine allocate_helmholtz_matrices(hmats, decomp, ss, ksp_type, pc_type, hypre_type, atol, rtol, max_its)
    type(helmholtz_matrices), intent(out) :: hmats
    type(decomposition), target, intent(in) :: decomp    
    real(kind = real_8), dimension(decomp%layers), intent(in) :: ss
    character(len = *), intent(in) :: ksp_type
    character(len = *), intent(in) :: pc_type
    character(len = *), optional, intent(in) :: hypre_type
    real(kind = real_8), optional, intent(in) :: atol
    real(kind = real_8), optional, intent(in) :: rtol
    integer, optional, intent(in) :: max_its
       
    integer :: g_m_ii, i, ierr, j, index, layer, lmax_its, m_ii
    integer, dimension(decomp%ii2 * decomp%jj2) :: gix, lix, nnz
    real(kind = real_8) :: latol, lrtol
    
    if(present(atol)) then
      latol = atol
    else
      latol = default_atol
    end if
    if(present(rtol)) then
      lrtol = rtol
    else
      lrtol = default_rtol
    end if
    if(present(max_its)) then
      lmax_its = max_its
    else
      lmax_its = default_max_its
    end if
    
    hmats%decomp => decomp
    
    allocate(hmats%ss(decomp%layers))
    hmats%ss = ss
    
    ! Assemble the LHS matrix. Note that the matrix is multiplied by -1 to
    ! guarantee positive definiteness.
    
    m_ii = decomp%ii2 * decomp%jj2
    g_m_ii = decomp%g_ii2 * decomp%g_jj2
    nnz = 5
    allocate(hmats%matrix(decomp%layers))
    allocate(hmats%ksp(decomp%layers))
    allocate(hmats%pc(decomp%layers))
    do layer = 1, decomp%layers
      ! \left( \partial_{xx} + \partial_{yy} - ss_i \right), discretised using
      ! finite differencing with dx = dy = 1 and with a boundary condition of
      ! psi = 0 on all boundaries
      call MatCreateAIJ(PETSC_COMM_WORLD, m_ii, m_ii, g_m_ii, g_m_ii, nnz(1), nnz, nnz(1), nnz, hmats%matrix(layer), ierr);  assert(ierr == 0)
      call MatSetUp(hmats%matrix(layer), ierr);  assert(ierr == 0)
      do i = 1, decomp%ii2
        do j = 1, decomp%jj2
          if(decomp%p_b >= 0 .or. i > 1)          call MatSetValue(hmats%matrix(layer), grid_ix(decomp, i, j), grid_ix(decomp, i - 1, j),     -1.0D0, ADD_VALUES, ierr);  assert(ierr == 0)
          if(decomp%p_c >= 0 .or. j > 1)          call MatSetValue(hmats%matrix(layer), grid_ix(decomp, i, j), grid_ix(decomp, i,     j - 1), -1.0D0, ADD_VALUES, ierr);  assert(ierr == 0)
                                                  call MatSetValue(hmats%matrix(layer), grid_ix(decomp, i, j), grid_ix(decomp, i,     j),      4.0D0 + ss(layer), ADD_VALUES, ierr);  assert(ierr == 0)
          if(decomp%p_d >= 0 .or. i < decomp%ii2) call MatSetValue(hmats%matrix(layer), grid_ix(decomp, i, j), grid_ix(decomp, i + 1, j),     -1.0D0, ADD_VALUES, ierr);  assert(ierr == 0)
          if(decomp%p_a >= 0 .or. j < decomp%jj2) call MatSetValue(hmats%matrix(layer), grid_ix(decomp, i, j), grid_ix(decomp, i,     j + 1), -1.0D0, ADD_VALUES, ierr);  assert(ierr == 0)
        end do
      end do      
      call MatAssemblyBegin(hmats%matrix(layer), MAT_FINAL_ASSEMBLY, ierr);  assert(ierr == 0)
      call MatAssemblyEnd(hmats%matrix(layer), MAT_FINAL_ASSEMBLY, ierr);  assert(ierr == 0)
    
      ! Configure the linear solver      
      call KSPCreate(PETSC_COMM_WORLD, hmats%ksp(layer), ierr);  assert(ierr == 0)
      call KSPSetOperators(hmats%ksp(layer), hmats%matrix(layer), hmats%matrix(layer), ierr);  assert(ierr == 0)
      call KSPGetPC(hmats%ksp(layer), hmats%pc(layer), ierr);  assert(ierr == 0)
      call PCSetType(hmats%pc(layer), pc_type, ierr);  assert(ierr == 0)
      if(pc_type == "hypre") then
#ifdef PETSC_HAVE_HYPRE
        if(.not. present(hypre_type)) then
          labort("hypre_type not supplied")
        end if
        call PCHYPRESetType(hmats%pc(layer), hypre_type, ierr);  assert(ierr == 0)
#else
        labort("No Hypre support")
#endif
      end if
      call PCSetFromOptions(hmats%pc(layer), ierr);  assert(ierr == 0)
      call KSPSetType(hmats%ksp(layer), ksp_type, ierr);  assert(ierr == 0)
      call KSPSetFromOptions(hmats%ksp(layer), ierr);  assert(ierr == 0)
      call KSPCGSetType(hmats%ksp(layer), KSP_CG_SYMMETRIC, ierr);  assert(ierr == 0)
      if(pc_type /= "lu") call KSPSetInitialGuessNonzero(hmats%ksp(layer), PETSC_TRUE, ierr) ;  assert(ierr == 0) 
      call KSPSetTolerances(hmats%ksp(layer), lrtol, latol, huge(0.0D0), lmax_its, ierr);  assert(ierr == 0)
    end do
      
    ! Allocate the RHS and (global) solution vectors
    
    call VecCreateMPI(PETSC_COMM_WORLD, m_ii, g_m_ii, hmats%rhs, ierr);  assert(ierr == 0)
    call VecCreateMPI(PETSC_COMM_WORLD, m_ii, g_m_ii, hmats%soln, ierr);  assert(ierr == 0)
    call VecZeroEntries(hmats%soln, ierr);  assert(ierr == 0)
    call VecAssemblyBegin(hmats%soln, ierr);  assert(ierr == 0)
    call VecAssemblyEnd(hmats%soln, ierr);  assert(ierr == 0)
    
    ! Allocate the local solution vector
    
    call VecCreateMPI(PETSC_COMM_SELF, m_ii, m_ii, hmats%lsoln, ierr);  assert(ierr == 0)
    
    ! Configure the global to local solution vector VecScatter
         
    allocate(hmats%gix(decomp%ii2, decomp%jj2)) 
    do j = 1, decomp%jj2 
      do i = 1, decomp%ii2
        hmats%gix(i, j) = grid_ix(decomp, i, j)
      end do
    end do
    allocate(hmats%lix(decomp%ii2, decomp%jj2)) 
    index = 0
    do j = 1, decomp%jj2 
      do i = 1, decomp%ii2
        hmats%lix(i, j) = index
        index = index + 1
      end do
    end do
    
    index = 0
    do j = 1, decomp%jj2 
      do i = 1, decomp%ii2
        index = index + 1
        gix(index) = hmats%gix(i, j)
        lix(index) = index - 1
      end do
    end do
    call ISCreateGeneral(PETSC_COMM_SELF, m_ii, gix, PETSC_COPY_VALUES, hmats%gis, ierr);  assert(ierr == 0)
    call ISCreateGeneral(PETSC_COMM_SELF, m_ii, lix, PETSC_COPY_VALUES, hmats%lis, ierr);  assert(ierr == 0)
    call VecScatterCreate(hmats%soln, hmats%gis, hmats%lsoln, hmats%lis, hmats%scatter, ierr);  assert(ierr == 0)
          
  end subroutine allocate_helmholtz_matrices
    
  ! Return the PETSc index associated with a given i and j.
  pure function grid_ix(decomp, i, j) result(ix)
    type(decomposition), intent(in) :: decomp
    integer, intent(in) :: i
    integer, intent(in) :: j
    
    integer :: ix
    
    ix = decomp%ii2_0 + i - 1 + modulo(decomp%jj2_0 + j - 1, decomp%g_jj2) * decomp%g_ii2
    
  end function grid_ix
  
  ! Deallocate a PETSc Helmholtz cache.
  subroutine deallocate_helmholtz_matrices(hmats)
    type(helmholtz_matrices), intent(inout) :: hmats
    
    integer :: ierr, layer
    
    do layer = 1, hmats%decomp%layers
      call MatDestroy(hmats%matrix(layer), ierr);  assert(ierr == 0)
      call KSPDestroy(hmats%ksp(layer), ierr);  assert(ierr == 0)
    end do
    call VecDestroy(hmats%rhs, ierr);  assert(ierr == 0)
    call VecDestroy(hmats%soln, ierr);  assert(ierr == 0)
    call VecDestroy(hmats%lsoln, ierr);  assert(ierr == 0)
    call ISDestroy(hmats%gis, ierr);  assert(ierr == 0)
    call ISDestroy(hmats%lis, ierr);  assert(ierr == 0)
    call VecScatterDestroy(hmats%scatter, ierr);  assert(ierr == 0)
    
    deallocate(hmats%ss)
    deallocate(hmats%matrix)
    deallocate(hmats%ksp)
    deallocate(hmats%pc)
    deallocate(hmats%gix)
    deallocate(hmats%lix)

  end subroutine deallocate_helmholtz_matrices
  
  ! Perform Helmholtz solves using PETSc. This does not apply the boundary
  ! condition to the solution. On input zeta is the RHS, on output zeta is the
  ! solution.
  subroutine solve_helmholtz(hmats, zeta)
    type(helmholtz_matrices), target, intent(inout) :: hmats
    real(kind = real_8), dimension(:, :, :), intent(inout) :: zeta
    
    integer :: ierr, its, j, layer
    KSPConvergedReason :: reason
    type(decomposition), pointer :: decomp
    
    assert(size(zeta, 1) == hmats%decomp%ii)
    assert(size(zeta, 2) == hmats%decomp%jj)
    assert(size(zeta, 3) == hmats%decomp%layers)
    
    decomp => hmats%decomp
      
    ! Copy the RHS into PETSc format
    do layer = 1, decomp%layers
      do j = 1, decomp%jj2 
        call VecSetValues(hmats%rhs, decomp%ii2, hmats%gix(:, j), zeta(decomp%l_i_min:decomp%l_i_max, decomp%l_j_min - 1 + j, layer), INSERT_VALUES, ierr);  assert(ierr == 0)
      end do
      call VecAssemblyBegin(hmats%rhs, ierr);  assert(ierr == 0)
      call VecAssemblyEnd(hmats%rhs, ierr);  assert(ierr == 0)
      ! The LHS matrix has been scaled by -1 to guarantee positive definiteness.
      ! Hence we scale the RHS.
      call VecScale(hmats%rhs, -1.0D0, ierr);  assert(ierr == 0) 
      
      ! Linear solve
      write(6, "(a)") "Performing linear solve using PETSc"
      call KSPSolve(hmats%ksp(layer), hmats%rhs, hmats%soln, ierr);  assert(ierr == 0)  
      call KSPGetConvergedReason(hmats%ksp(layer), reason, ierr);  assert(ierr == 0)
      call KSPGetIterationNumber(hmats%ksp(layer), its, ierr);  assert(ierr == 0)
      write(6, "(a,"//real_8_chr//")") "For ss = ", hmats%ss(layer)
      write(6, "(a,i0)") "Convergence reason: ", reason%v
      write(6, "(a,i0)") "Iterations: ", its
      if(reason%v <= 0) then
        labort("Convergence failure")
      end if
       
      ! Communicate the solution
      call VecScatterBegin(hmats%scatter, hmats%soln, hmats%lsoln, INSERT_VALUES, SCATTER_FORWARD, ierr);  assert(ierr == 0)
      call VecScatterEnd(hmats%scatter, hmats%soln, hmats%lsoln, INSERT_VALUES, SCATTER_FORWARD, ierr);  assert(ierr == 0)
          
      ! Copy the solution from PETSc format 
      do j = 1, decomp%jj2 
        call VecGetValues(hmats%lsoln, decomp%ii2, hmats%lix(:, j), zeta(decomp%l_i_min:decomp%l_i_max, decomp%l_j_min - 1 + j, layer), ierr);  assert(ierr == 0)
      end do 
    end do
    
  end subroutine solve_helmholtz

end module qg_helmholtz_petsc

#endif
