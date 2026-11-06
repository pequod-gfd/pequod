#include "config.h"
#include "leapfrog.h"

! Auxilliary routines for use by the qg3 main loop. These routines may directly
! access parameters in qg_parameters (e.g. for static memory allocation).
module qg_aux

  use lapack_interfaces
  use qg_diagnostics
#ifdef PETSC
  use qg_helmholtz_petsc, only : helmholtz_solver_cache => helmholtz_matrices
  use qg_helmholtz_petsc, only : allocate, deallocate, solve_helmholtz
#else
  use qg_helmholtz, only : solve_helmholtz => helm, &
    & helmholtz_solver_cache => helm_cache
  use qg_helmholtz, only : allocate, deallocate
#endif
  use qg_laplace_bc
  use qg_parallel
  use qg_parameters, only : INFINITY, real_8, real_8_chr
#ifndef CABARET
  use qg_parameters, only : ra_nu
#ifdef RA_DOUBLE_FILTER
  use qg_parameters, only : ra_alpha
#endif
#endif
  use test_tools

  implicit none
  
  private
 
  public :: check_decomp, pv_inversion, rel_from_psi, zeta_from_rel_and_psi, &
    & update_zeta_bc, rel_from_zeta_and_psi, apply_rel_bc
  public :: allocate, deallocate
  public :: helmholtz_cache
#ifdef CABARET
  public :: cabaret_adapt_dt, cabaret_cfl, cabaret_predictor, cabaret_beta, &
    & cabaret_forcing, cabaret_cell_boundaries, cabaret_extrapolator, &
    & cabaret_flux, cabaret_corrector
#else
  public :: zeta_forward_euler, zeta_centred, robert_asselin_filter
#endif
#ifdef BACKGROUND_HSHEAR
  public :: initialise_background
#endif

#ifndef CABARET
    ! Robert-Asselin-Williams filter weights
#ifdef RA_DOUBLE_FILTER
                                      ! Second filter weights
    real(kind = real_8), parameter :: ra_weight21 = 0.5D0 * ra_nu * ra_alpha, &
                                    & ra_weight22 = 1.0D0 - 2.0D0 * ra_weight21, &
                                    !&  ra_weight13 = ra_weight11, &
                                      ! First filter weights
                                    & ra_weight11 = -0.5D0 * ra_nu * (1.0D0 - ra_alpha), &
                                    & ra_weight12 = -2.0D0 * ra_weight11, &
                                    & ra_weight13 = 1.0D0 + ra_weight11
#else
    real(kind = real_8), parameter :: ra_weight1 = 0.5D0 * ra_nu, &
                                    & ra_weight2 = 1.0D0 - 2.0D0 * ra_weight1!, &
                                    !&  ra_weight3 = ra_weight1
#endif
#endif
 
  ! Caches data associated with the Helmholtz solve.
  type helmholtz_cache
    ! Cache used by solve_helmholtz / helm
    type(helmholtz_solver_cache) :: helm_c
    ! Normalised boundary contributions to the stream function
    real(kind = real_8), dimension(:, :, :, :), pointer :: corr
    ! Inverse of LHS matrix for the channel constraint
    real(kind = real_8), dimension(:, :, :), pointer :: inv_corr_a
    ! Working memory used within pv_inversion
    real(kind = real_8), dimension(:, :, :), pointer :: phi
  end type helmholtz_cache

  interface allocate
    module procedure allocate_helmholtz_cache
  end interface allocate
  
  interface deallocate
    module procedure deallocate_helmholtz_cache
  end interface deallocate
  
#ifdef CABARET
  interface cabaret_cell_boundaries
    module procedure cabaret_cell_boundaries_velocity, &
      & cabaret_cell_boundaries_time_integrated
  end interface cabaret_cell_boundaries
#else
  interface robert_asselin_filter
    module procedure robert_asselin_filter_rank1, robert_asselin_filter_rank2, &
      & robert_asselin_filter_rank3, robert_asselin_filter_pdata
  end interface robert_asselin_filter
#endif

contains

  ! Check the decomposition parameters.
  subroutine check_decomp()
    use qg_parameters, only : iis, ii_m, jjs, jj_m, layers, periodic
  
    integer :: g_ii, g_jj, n
    logical :: warn
    
    warn = .false.
    
    n = nprocs()
    g_ii = sum(iis)
    g_jj = sum(jjs)
  
    if(size(iis) * size(jjs) /= n) then
      write(6, "(a,i0,a,i0)") "Decomposition topology: ", size(iis), " x ", size(jjs)
      write(6, "(a,i0)") "Number of processes: ", n
      labort("check_decomp: Invalid number of processes in decomposition")
    end if
    
#ifdef CABARET
    if(minval(iis) <= 2 .or. minval(jjs) <= 2) then
#else
    if(minval(iis) <= 1 .or. minval(jjs) <= 1) then
#endif
      labort("check_decomp: Insufficient partition size")
    end if
    if(maxval(iis) /= ii_m .or. maxval(jjs) /= jj_m) then
      labort("check_decomp: Incorrect decomposition size upper bounds")
    end if
    
    if(2 ** nint(dlog(g_ii - 1.0D0) / dlog(2.0D0)) /= g_ii - 1) then
      write(6, "(a,i0,a,i0)") "Domain size: ", g_ii, " x ", g_jj
      labort("check_decomp: Latitudinal decomposition must equal 2^n + 1 for some natural number n")
    end if
    if(.not. periodic .and. g_ii /= g_jj) then
      write(0, "(a,i0,a,i0)") "Domain size: ", g_ii, " x ", g_jj
      write(0, "(a)") "check_decomp: WARNING - Domain not square"
      warn = .true.
    end if
    
    if(layers <= 0) then
      write(6, "(a,i0)") "Number of layers: ", layers
      labort("check_decomp: Zero or negative number of layers")
    end if
    
    if(.not. warn) then
      write(6, "(a)") "check_decomp: success"
    end if
  
  end subroutine check_decomp

  ! Allocate the Helmholtz solver cache.
  subroutine allocate_helmholtz_cache(hcache, decomp, strat, vbc)
    type(helmholtz_cache), intent(out) :: hcache
    type(decomposition), intent(in) :: decomp
    type(stratification), intent(in) :: strat
    type(velocity_boundary_condition), intent(in) :: vbc
    
    ! Allocate the correction functions
    if(decomp%periodic_j) then
      allocate(hcache%corr(decomp%ii, decomp%jj, decomp%layers, 2))
      allocate(hcache%inv_corr_a(2, 2, decomp%layers - 1))
    else
      allocate(hcache%corr(decomp%ii, decomp%jj, decomp%layers, 1))
      nullify(hcache%inv_corr_a)
    end if
    
    ! Allocate and initialise the Helmholtz solvers
#ifdef PETSC
    if(decomp%n > 1) then
#ifdef SUPERLU_DIST
      call allocate(hcache%helm_c, decomp, strat%ss, ksp_type = "preonly", pc_type = "lu")
#else
      call allocate(hcache%helm_c, decomp, strat%ss, ksp_type = "cg", pc_type = "sor", atol = 1.0D-16, rtol = 1.0D-7)
#endif
    else
      call allocate(hcache%helm_c, decomp, strat%ss, ksp_type = "preonly", pc_type = "lu")
    end if
#else
    call allocate(hcache%helm_c, decomp, strat%ss)
#endif

    ! Compute the correction functions
    if(decomp%periodic_j) then
      call correction_function_momentum(hcache%helm_c, strat, vbc, hcache%corr, hcache%inv_corr_a)
    else
      call correction_function_mass(hcache%helm_c, hcache%corr(:, :, :, 1))
    end if
    
    ! Allocate the working memory
    allocate(hcache%phi(decomp%ii, decomp%jj, decomp%layers))
    hcache%phi = 0.0D0
    
  end subroutine allocate_helmholtz_cache

  ! Compute the correction functions used to impose mass conservation.
  subroutine correction_function_mass(helm_c, corr)
    type(helmholtz_solver_cache), intent(inout) :: helm_c
    real(kind = real_8), dimension(:, :, :), intent(out) :: corr
    
    integer :: layer
    real(kind = real_8), dimension(helm_c%decomp%layers) :: mean_corr
    type(decomposition), pointer :: decomp
       
    assert(.not. helm_c%decomp%periodic_j)
    assert(size(corr, 1) == helm_c%decomp%ii)
    assert(size(corr, 2) == helm_c%decomp%jj)
    assert(size(corr, 3) == helm_c%decomp%layers)
    
    decomp => helm_c%decomp

    ! Solve:
    !   (\partial_xx + \partial_yy + S) \psi = 1
    ! subject to a Dirichlet boundary condition \psi = 0 on all boundaries
    corr = 1.0D0
    call solve_helmholtz(helm_c, corr)
    call apply_dirichlet_bc(decomp, 0.0D0, corr)

    ! Scale the correction function. The mass conservation correction then
    ! equals -\int psi_0 times the correction function, where \psi_0 is the
    ! solution to the problem with homogeneous boundary conditions.
    do layer = 1, helm_c%decomp%layers
      mean_corr(layer) = sum(corr(decomp%l_i_min:decomp%l_i_max, decomp%l_j_min:decomp%l_j_max, layer))
    end do
    call allsum(mean_corr)
    mean_corr = mean_corr / real(decomp%g_ii_int * decomp%g_jj_int, kind = real_8)
    do layer = 1, helm_c%decomp%layers
      corr(:, :, layer) = (1.0D0 + helm_c%ss(layer) * corr(:, :, layer)) / (1.0D0 + helm_c%ss(layer) * mean_corr(layer))    
    end do
    
  end subroutine correction_function_mass
  
  ! Compute the correction functions used for the channel.
  subroutine correction_function_momentum(helm_c, strat, vbc, corr, inv_corr_a)
    type(helmholtz_solver_cache), intent(inout) :: helm_c
    type(stratification), intent(in) :: strat
    type(velocity_boundary_condition), intent(in) :: vbc
    real(kind = real_8), dimension(:, :, :, :), intent(out) :: corr
    real(kind = real_8), dimension(:, :, :), intent(out) :: inv_corr_a
    
    integer :: mode, layers
    real(kind = real_8), dimension(helm_c%decomp%layers, 2, 2) :: z_b
    type(decomposition), pointer :: decomp
       
    assert(helm_c%decomp%periodic_j)
    assert(size(corr, 1) == helm_c%decomp%ii)
    assert(size(corr, 2) == helm_c%decomp%jj)
    assert(size(corr, 3) == helm_c%decomp%layers)
    assert(size(corr, 4) == 2)
    assert(size(inv_corr_a, 1) == 2)
    assert(size(inv_corr_a, 2) == 2)
    assert(size(inv_corr_a, 3) == helm_c%decomp%layers - 1)
    
    decomp => helm_c%decomp
    layers = decomp%layers
    
    corr = 0.0D0
    
    ! Solve Helmholtz equations with zero RHS subject to bc:
    !   psi = 1 on upper boundary
    !   psi = 0 on lower boundary
    call add_laplace_bc_rhs(decomp, bc_b = 1.0D0, bc_d = 0.0D0, rhs = corr(:, :, :, 1))
    call solve_helmholtz(helm_c, corr(:, :, :, 1))
    call apply_dirichlet_bc(decomp, bc_b = 1.0D0, bc_d = 0.0D0, zeta = corr(:, :, :, 1))
    
    ! Solve Helmholtz equations with zero RHS subject to bc:
    !   psi = 0 on upper boundary
    !   psi = 1 on lower boundary
    call add_laplace_bc_rhs(decomp, bc_b = 0.0D0, bc_d = 1.0D0, rhs = corr(:, :, :, 2))
    call solve_helmholtz(helm_c, corr(:, :, :, 2))
    call apply_dirichlet_bc(decomp, bc_b = 0.0D0, bc_d = 1.0D0, zeta = corr(:, :, :, 2))

    ! Set up the LHS of the linear solves for the correction function
    ! coefficients    
#ifdef DDEBUG
    inv_corr_a = INFINITY
#endif

    call zeta_boundary(decomp, strat, vbc, corr(:, :, :, 1), z_b(:, 1, 1), z_b(:, 2, 1))
    call zeta_boundary(decomp, strat, vbc, corr(:, :, :, 2), z_b(:, 1, 2), z_b(:, 2, 2))

    do mode = 2, layers
      inv_corr_a(:, :, mode - 1) = inverse(z_b(mode, :, :))
    end do
    
  end subroutine correction_function_momentum  

  ! Compute the mean relative potential vorticity on the boundaries. Input and
  ! output are modal variables.
  subroutine zeta_boundary(decomp, strat, vbc, psi, z_b, z_d)
    type(decomposition), intent(in) :: decomp
    type(stratification), intent(in) :: strat
    type(velocity_boundary_condition), intent(in) :: vbc
    real(kind = real_8), dimension(:, :, :), intent(in) :: psi
    real(kind = real_8), dimension(decomp%layers), intent(out) :: z_b
    real(kind = real_8), dimension(decomp%layers), intent(out) :: z_d
    
    integer :: j, layer, layers
    real(kind = real_8) :: cff
    
    assert(strat%layers == decomp%layers)
    assert(decomp%periodic_j)
    assert(size(psi, 1) == decomp%ii)
    assert(size(psi, 2) == decomp%jj)
    assert(size(psi, 3) == decomp%layers)
    
    layers = decomp%layers
    
    z_b = 0.0D0
    z_d = 0.0D0
    
    ! Relative vorticity
    
    select case(vbc%bc_type)
      case(BC_PARTIAL_SLIP)
        cff = vbc%alpha / (vbc%alpha + 3.0D0)
        if(decomp%p_b == -1) then
          do layer = 1, layers
            do j = decomp%l_j_min, decomp%l_j_max
              z_b(layer) = z_b(layer) + cff * (-3.5D0 * psi(1, j, layer) &
                                            & + 4.0D0 * psi(2, j, layer) &
                                            & - 0.5D0 * psi(3, j, layer))
            end do
          end do
        end if
        if(decomp%p_d == -1) then
          do layer = 1, layers
            do j = decomp%l_j_min, decomp%l_j_max
              z_d(layer) = z_d(layer) + cff * (-3.5D0 * psi(decomp%ii, j, layer) &
                                            & + 4.0D0 * psi(decomp%ii - 1, j, layer) &
                                            & - 0.5D0 * psi(decomp%ii - 2, j, layer))
            end do
          end do
        end if          
      case(BC_FREE_SLIP)
        !z_b = z_b + 0.0D0
        !z_d = z_d + 0.0D0
      case(BC_NO_SLIP)
        if(decomp%p_b == -1) then
          do layer = 1, layers
            do j = decomp%l_j_min, decomp%l_j_max
              z_b(layer) = z_b(layer) + 2.0D0 * (psi(2, j, layer) - psi(1, j, layer))
            end do
          end do
        end if
        if(decomp%p_d == -1) then
          do layer = 1, layers
            do j = decomp%l_j_min, decomp%l_j_max
              z_d(layer) = z_d(layer) + 2.0D0 * (psi(decomp%ii - 1, j, layer) - psi(decomp%ii, j, layer))
            end do
          end do
        end if
      case default
        write(0, "(a,i0)") "For boundary condition type: ", vbc%bc_type
        labort("Unrecognised boundary condition type")
    end select
    
    ! Thickness
    
    if(decomp%p_b == -1) then
      do layer = 1, layers
        do j = decomp%l_j_min, decomp%l_j_max
          z_b(layer) = z_b(layer) - strat%ss(layer) * psi(1, j, layer)
        end do
      end do
    end if
    if(decomp%p_d == -1) then
      do layer = 1, layers
        do j = decomp%l_j_min, decomp%l_j_max
          z_d(layer) = z_d(layer) - strat%ss(layer) * psi(decomp%ii, j, layer)
        end do
      end do
    end if
    
    call allsum(z_b)
    call allsum(z_d)    
    z_b = z_b / real(decomp%g_jj_int, kind = real_8)
    z_d = z_d / real(decomp%g_jj_int, kind = real_8)
    
  end subroutine zeta_boundary

  ! Deallocate the Helmholtz solver cache.
  subroutine deallocate_helmholtz_cache(hcache)
    type(helmholtz_cache), intent(inout) :: hcache
    
    call deallocate(hcache%helm_c)
    deallocate(hcache%corr)
    if(associated(hcache%inv_corr_a)) deallocate(hcache%inv_corr_a)
    deallocate(hcache%phi)
    
  end subroutine deallocate_helmholtz_cache

#ifdef BACKGROUND_HSHEAR
  ! Initialise the background stream function. This sets the halo region for
  ! psi_b but *not* for rel_b.
  subroutine initialise_background(psi_b, rel_b)
    use qg_parameters, only : u_background_h, u_scale
  
    type(parallel_data), intent(inout) :: psi_b
    type(parallel_data), intent(inout) :: rel_b
    
    integer :: gi, i, j
    real(kind = real_8) :: y, yi
    real(kind = real_8), dimension(:), allocatable :: psi, u
    type(decomposition), pointer :: decomp
    type(velocity_boundary_condition) :: vbc
    
    decomp => psi_b%decomp
    
    if(decomp%layers > 1) then
      write(0, "(a,i0)") "Number of layers: ", decomp%layers
      labort("Background flow should only be used with a barotropic system")
    end if
    
    ! Compute the (global) background velocity between nodes
    allocate(u(decomp%g_ii  - 1))
    do i = 1, decomp%g_ii - 1
      yi = real(decomp%g_ii - i, kind = real_8) - 0.5D0
      y  = yi / real(decomp%g_ii - 1, kind = real_8)
      u(i) = u_background_h(y) / u_scale
    end do
    ! Compute the (global) background stream function
    allocate(psi(decomp%g_ii))
    psi(decomp%g_ii) = 0.0D0
    do i = 2, decomp%g_ii
      psi(decomp%g_ii - i + 1) = psi(decomp%g_ii - i + 2) - u(decomp%g_ii - i + 1)
    end do
    
    ! Set the (local) background stream function
    do j = decomp%fa_j_min, decomp%fa_j_max
      do i = decomp%fa_i_min, decomp%fa_i_max
        gi = i - 1 + decomp%ii_0
        psi_b%data(i, j, 1) = psi(gi)
      end do
    end do
    deallocate(psi)
    
    ! Set the (local) relative vorticity
    do j = decomp%f_j_min, decomp%f_j_max
      do i = decomp%f_i_min, decomp%f_i_max
        rel_b%data(i, j, 1) = psi_b%data(i + 1, j, 1) + psi_b%data(i - 1, j, 1) - 2.0D0 * psi_b%data(i, j, 1)
      end do
    end do    
    ! Apply free slip boundary conditions to the background flow profile.
    call initialise(vbc, BC_FREE_SLIP)
    call apply_rel_bc(vbc, psi_b, rel_b)
      
  end subroutine initialise_background
#endif
  
  ! Invert the relative PV for the stream function. This does *not* update the
  ! halo region.
  subroutine pv_inversion(strat, vbc, zeta, hcache, psi)
    type(stratification), intent(in) :: strat
    type(velocity_boundary_condition), intent(in) :: vbc
    type(parallel_data), intent(in) :: zeta
    type(helmholtz_cache), intent(inout) :: hcache
    type(parallel_data), intent(inout) :: psi
        
    type(decomposition), pointer :: decomp

    assert(associated(zeta%decomp, psi%decomp))
    
    decomp => psi%decomp
    
    ! Convert from layers to modes
    call z_to_modes(strat, zeta%sdata, hcache%phi)
     
    ! Solve the problem with homogeneous boundary conditions
    call solve_helmholtz(hcache%helm_c, hcache%phi)     
    call apply_dirichlet_bc(decomp, 0.0D0, hcache%phi)
    
    ! Convert from modes to layers, applying corrections
    if(decomp%periodic_j) then
      call modes_to_z_correct_momentum(decomp, strat, vbc, hcache%phi, psi%sdata, hcache%corr, hcache%inv_corr_a)
    else
      call modes_to_z_correct_mass(decomp, strat, hcache%phi, psi%sdata, hcache%corr(:, :, :, 1))
    end if
        
  end subroutine pv_inversion
  
  ! Convert from modes to layers, applying mass correction.
  subroutine modes_to_z_correct_mass(decomp, strat, phi, psi, corr)
    type(decomposition), intent(in) :: decomp
    type(stratification), intent(in) :: strat
    real(kind = real_8), dimension(:, :, :), intent(in) :: phi
    real(kind = real_8), dimension(:, :, :), intent(out) :: psi
    real(kind = real_8), dimension(:, :, :), intent(in) :: corr
    
    integer :: layer
    real(kind = real_8), dimension(decomp%layers) :: psiC
    
    assert(size(phi, 1) == decomp%ii)
    assert(size(phi, 2) == decomp%jj)
    assert(size(phi, 3) == decomp%layers)
    assert(size(psi, 1) == decomp%ii)
    assert(size(psi, 2) == decomp%jj)
    assert(size(psi, 3) == decomp%layers)
    assert(size(corr, 1) == decomp%ii)
    assert(size(corr, 2) == decomp%jj)
    assert(size(corr, 3) == decomp%layers)
    
    psiC = mean(decomp, phi)
      
    if(any(psiC /= psiC)) then
      write(6, "(a)") "psiC"
      do layer = 1, decomp%layers
        write(6, "(a,i0,a,"//real_8_chr//")") "Layer ", layer, ": ", psiC(layer)
      end do
      labort("NaN psiC encountered")
    else if(any(dabs(psiC) > 1.0D100)) then
      write(6, "(a)") "psiC"
      do layer = 1, decomp%layers
        write(6, "(a,i0,a,"//real_8_chr//")") "Layer ", layer, ": ", psiC(layer)
      end do
      write(6, "(a)") "WARNING: Large psiC encountered"
      write(6, "(a)") "Simulation may be blowing up"
    end if
  
    call modes_to_z_btcu(strat, phi, psi, psiC, corr) 
      
  end subroutine modes_to_z_correct_mass
  
  ! Convert from modes to layers, applying appropriate corrections.
  subroutine modes_to_z_correct_momentum(decomp, strat, vbc, phi, psi, corr, inv_corr_a)
    type(decomposition), intent(in) :: decomp
    type(stratification), intent(in) :: strat
    type(velocity_boundary_condition), intent(in) :: vbc
    real(kind = real_8), dimension(:, :, :), intent(inout) :: phi
    real(kind = real_8), dimension(:, :, :), intent(out) :: psi
    real(kind = real_8), dimension(:, :, :, :), intent(in) :: corr
    real(kind = real_8), dimension(:, :, :), intent(in) :: inv_corr_a
    
    integer :: i, j, layer
    real(kind = real_8), dimension(2) :: psiC
    real(kind = real_8), dimension(decomp%layers, 2) :: corr_b
    
    assert(size(phi, 1) == decomp%ii)
    assert(size(phi, 2) == decomp%jj)
    assert(size(phi, 3) == decomp%layers)
    assert(size(psi, 1) == decomp%ii)
    assert(size(psi, 2) == decomp%jj)
    assert(size(psi, 3) == decomp%layers)
    assert(size(corr, 1) == decomp%ii)
    assert(size(corr, 2) == decomp%jj)
    assert(size(corr, 3) == decomp%layers)
    assert(size(corr, 4) == 2)
    assert(size(inv_corr_a, 1) == 2)
    assert(size(inv_corr_a, 2) == 2)
    assert(size(inv_corr_a, 3) == decomp%layers - 1)
    
    call zeta_boundary(decomp, strat, vbc, phi, corr_b(:, 1), corr_b(:, 2))
    do layer = 2, decomp%layers
      psiC = matmul(inv_corr_a(:, :, layer - 1), corr_b(layer, :))
    
      if(any(psiC /= psiC)) then
        write(6, "(a)") "psiC"
        write(6, "(a,i0,a,"//real_8_chr//",a,"//real_8_chr//")") "Layer ", layer, ": ", psiC(1), ", ", psiC(2)
        labort("NaN psiC encountered")
      else if(any(dabs(psiC) > 1.0D100)) then
        write(6, "(a)") "psiC"
        write(6, "(a,i0,a,"//real_8_chr//",a,"//real_8_chr//")") "Layer ", layer, ": ", psiC(1), ", ", psiC(2)
        write(6, "(a)") "WARNING: Large psiC encountered"
        write(6, "(a)") "Simulation may be blowing up"
      end if

      do j = 1, decomp%jj
        do i = 1, decomp%ii
          phi(i, j, layer) = phi(i, j, layer) - psiC(1) * corr(i, j, layer, 1) - psiC(2) * corr(i, j, layer, 2)
        end do
      end do
    end do
    
    call modes_to_z(strat, phi, psi)
    
  end subroutine modes_to_z_correct_momentum
  
  ! Compute the relative vorticity from the stream function. This does *not*
  ! update the halo region.
  subroutine rel_from_psi(vbc, psi, rel)
    type(velocity_boundary_condition), intent(in) :: vbc
    type(parallel_data), intent(in) :: psi
    type(parallel_data), intent(inout) :: rel
    
    integer :: i, j, layer
    type(decomposition), pointer :: decomp
    
    decomp => psi%decomp 
    
    do layer = 1, decomp%layers
      do j = decomp%f_j_min, decomp%f_j_max
        do i = decomp%f_i_min, decomp%f_i_max
          rel%data(i, j, layer) = -4.0D0 * psi%data(i, j, layer) &
                                & +        psi%data(i - 1, j, layer) &
                                & +        psi%data(i, j + 1, layer) &
                                & +        psi%data(i + 1, j, layer) &
                                & +        psi%data(i, j - 1, layer)
        end do
      end do
    end do
    
    call apply_rel_bc(vbc, psi, rel)
    
  end subroutine rel_from_psi
  
  ! Compute the relative potential vorticity from the relative vorticity and
  ! stream function. This does *not* update the halo region
  subroutine zeta_from_rel_and_psi(strat, rel, psi, zeta)
    type(stratification), intent(in) :: strat
    type(parallel_data), intent(in) :: rel
    type(parallel_data), intent(in) :: psi
    type(parallel_data), intent(inout) :: zeta
    
    integer :: i, j, layer, layers
    type(decomposition), pointer :: decomp
    
    layers = psi%decomp%layers
    if(layers == 1) then
      call set(zeta, rel)
      return
    end if
    
    decomp => psi%decomp
    
    do j = 1, decomp%jj
      do i = 1, decomp%ii
        zeta%sdata(i, j, 1) = rel%sdata(i, j, 1) - strat%s2(1) * (psi%sdata(i, j, 1) - psi%sdata(i, j, 2))
        zeta%sdata(i, j, layers) = rel%sdata(i, j, layers) - strat%s1(layers) * (psi%sdata(i, j, layers) - psi%sdata(i, j, layers - 1))
      end do
    end do
    do layer = 2, layers - 1
      do j = 1, decomp%jj
        do i = 1, decomp%ii
          zeta%sdata(i, j, layer) = rel%sdata(i, j, layer) - strat%s1(layer) * (psi%sdata(i, j, layer) - psi%sdata(i, j, layer - 1)) &
                                                         & - strat%s2(layer) * (psi%sdata(i, j, layer) - psi%sdata(i, j, layer + 1))
        end do
      end do
    end do
    
  end subroutine zeta_from_rel_and_psi
  
  ! Update the boundary values of the relative potential vorticity from the
  ! relative vorticity and stream function. This updates the halo region.
  subroutine update_zeta_bc(strat, rel, psi, zeta)
    type(stratification), intent(in) :: strat
    type(parallel_data), intent(in) :: rel
    type(parallel_data), intent(in) :: psi
    type(parallel_data), intent(inout) :: zeta
    
    integer :: i, j, layer, layers
    type(decomposition), pointer :: decomp
    
    decomp => psi%decomp
    layers = decomp%layers
    if(layers == 1) then
      ! Left boundary
      if(decomp%p_c == -1) then
        do i = decomp%fa_i_min, decomp%fa_i_max
          zeta%data(i, 2, 1) = rel%data(i, 2, 1)
        end do
      end if
      ! Right boundary
      if(decomp%p_a == -1) then
        do i = decomp%fa_i_min, decomp%fa_i_max
          zeta%data(i, decomp%jj + 1, 1) = rel%data(i, decomp%jj + 1, 1)
        end do
      end if
      ! Top boundary
      if(decomp%p_b == -1) then
        do j = decomp%fa_j_min, decomp%fa_j_max
          zeta%data(2, j, 1) = rel%data(2, j, 1)
        end do
      end if
      ! Bottom boundary
      if(decomp%p_d == -1) then
        do j = decomp%fa_j_min, decomp%fa_j_max
          zeta%data(decomp%ii + 1, j, 1) = rel%data(decomp%ii + 1, j, 1)
        end do
      end if
    else
      ! Left boundary
      if(decomp%p_c == -1) then
        do i = decomp%fa_i_min, decomp%fa_i_max
          zeta%data(i, 2, 1) = rel%data(i, 2, 1) - strat%s2(1) * (psi%data(i, 2, 1) - psi%data(i, 2, 2))
          zeta%data(i, 2, layers) = rel%data(i, 2, layers) - strat%s1(layers) * (psi%data(i, 2, layers) - psi%data(i, 2, layers - 1))
        end do
        do layer = 2, layers - 1
          do i = decomp%fa_i_min, decomp%fa_i_max
            zeta%data(i, 2, layer) = rel%data(i, 2, layer) - strat%s1(layer) * (psi%data(i, 2, layer) - psi%data(i, 2, layer - 1)) &
                                                         & - strat%s2(layer) * (psi%data(i, 2, layer) - psi%data(i, 2, layer + 1))
          end do
        end do
      end if
      ! Right boundary
      if(decomp%p_a == -1) then
        do i = decomp%fa_i_min, decomp%fa_i_max
          zeta%data(i, decomp%jj + 1, 1) = rel%data(i, decomp%jj + 1, 1) - strat%s2(1) * (psi%data(i, decomp%jj + 1, 1) - psi%data(i, decomp%jj + 1, 2))
          zeta%data(i, decomp%jj + 1, layers) = rel%data(i, decomp%jj + 1, layers) - strat%s1(layers) * (psi%data(i, decomp%jj + 1, layers) - psi%data(i, decomp%jj + 1, layers - 1))
        end do
        do layer = 2, layers - 1
          do i = decomp%fa_i_min, decomp%fa_i_max
            zeta%data(i, decomp%jj + 1, layer) = rel%data(i, decomp%jj + 1, layer) - strat%s1(layer) * (psi%data(i, decomp%jj + 1, layer) - psi%data(i, decomp%jj + 1, layer - 1)) &
                                                                                 & - strat%s2(layer) * (psi%data(i, decomp%jj + 1, layer) - psi%data(i, decomp%jj + 1, layer + 1))
          end do
        end do
      end if
      ! Top boundary
      if(decomp%p_b == -1) then
        do j = decomp%fa_j_min, decomp%fa_j_max
          zeta%data(2, j, 1) = rel%data(2, j, 1) - strat%s2(1) * (psi%data(2, j, 1) - psi%data(2, j, 2))
          zeta%data(2, j, layers) = rel%data(2, j, layers) - strat%s1(layers) * (psi%data(2, j, layers) - psi%data(2, j, layers - 1))
        end do
        do layer = 2, layers - 1
          do j = decomp%fa_j_min, decomp%fa_j_max
            zeta%data(2, j, layer) = rel%data(2, j, layer) - strat%s1(layer) * (psi%data(2, j, layer) - psi%data(2, j, layer - 1)) &
                                                         & - strat%s2(layer) * (psi%data(2, j, layer) - psi%data(2, j, layer + 1))
          end do
        end do
      end if
      ! Bottom boundary
      if(decomp%p_d == -1) then
        do j = decomp%fa_j_min, decomp%fa_j_max
          zeta%data(decomp%ii + 1, j, 1) = rel%data(decomp%ii + 1, j, 1) - strat%s2(1) * (psi%data(decomp%ii + 1, j, 1) - psi%data(decomp%ii + 1, j, 2))
          zeta%data(decomp%ii + 1, j, layers) = rel%data(decomp%ii + 1, j, layers) - strat%s1(layers) * (psi%data(decomp%ii + 1, j, layers) - psi%data(decomp%ii + 1, j, layers - 1))
        end do
        do layer = 2, layers - 1
          do j = decomp%fa_j_min, decomp%fa_j_max
            zeta%data(decomp%ii + 1, j, layer) = rel%data(decomp%ii + 1, j, layer) - strat%s1(layer) * (psi%data(decomp%ii + 1, j, layer) - psi%data(decomp%ii + 1, j, layer - 1)) &
                                                                                 & - strat%s2(layer) * (psi%data(decomp%ii + 1, j, layer) - psi%data(decomp%ii + 1, j, layer + 1))
          end do
        end do
      end if
    end if
    
  end subroutine update_zeta_bc
  
  ! Compute the relative vorticity from the relative potential vorticity and
  ! stream function. This updates the halo region.
  subroutine rel_from_zeta_and_psi(strat, vbc, zeta, psi, rel)
    type(stratification), intent(in) :: strat
    type(velocity_boundary_condition), intent(in) :: vbc
    type(parallel_data), intent(in) :: zeta
    type(parallel_data), intent(in) :: psi
    type(parallel_data), intent(inout) :: rel
    
    integer :: i, j, layer, layers
    type(decomposition), pointer :: decomp
    
    decomp => psi%decomp  
    layers = decomp%layers
    if(layers == 1) then
      call set(rel, zeta)
    else
      do j = decomp%fa_j_min, decomp%fa_j_max
        do i = decomp%fa_i_min, decomp%fa_i_max
          rel%data(i, j, 1) = zeta%data(i, j, 1) + strat%s2(1) * (psi%data(i, j, 1) - psi%data(i, j, 2))
        end do
      end do
      do layer = 2, layers - 1
        do j = decomp%fa_j_min, decomp%fa_j_max
          do i = decomp%fa_i_min, decomp%fa_i_max
            rel%data(i, j, layer) = zeta%data(i, j, layer) + strat%s1(layer) * (psi%data(i, j, layer) - psi%data(i, j, layer - 1)) &
                                                         & + strat%s2(layer) * (psi%data(i, j, layer) - psi%data(i, j, layer + 1))
          end do
        end do
      end do
      do j = decomp%fa_j_min, decomp%fa_j_max
        do i = decomp%fa_i_min, decomp%fa_i_max
          rel%data(i, j, layers) = zeta%data(i, j, layers) + strat%s1(layers) * (psi%data(i, j, layers) - psi%data(i, j, layers - 1))
        end do
      end do
    end if
    
    call apply_rel_bc(vbc, psi, rel)
    
  end subroutine rel_from_zeta_and_psi
  
  ! Apply boundary conditions to the relative vorticity. This updates the halo
  ! region.
  subroutine apply_rel_bc(vbc, psi, rel)
    type(velocity_boundary_condition), intent(in) :: vbc
    type(parallel_data), intent(in) :: psi
    type(parallel_data), intent(inout) :: rel
    
    integer :: i, j, layer, layers
    real(kind = real_8) :: cff
    type(decomposition), pointer :: decomp
    
    decomp => psi%decomp
    layers = decomp%layers
    
    select case(vbc%bc_type)
      case(BC_PARTIAL_SLIP)
        cff = vbc%alpha / (vbc%alpha + 3.0D0)
      
        ! Top boundary condition
        if(decomp%p_b == -1) then
          do layer = 1, layers
            do j = decomp%fa_j_min, decomp%fa_j_max
              rel%data(2, j, layer) = cff * (-3.5D0 * psi%data(2, j, layer) &
                                          & + 4.0D0 * psi%data(3, j, layer) &
                                          & - 0.5D0 * psi%data(4, j, layer))
            end do
          end do
        end if
        
        ! Bottom boundary condition
        if(decomp%p_d == -1) then
          do layer = 1, layers
            do j = decomp%fa_j_min, decomp%fa_j_max
              rel%data(decomp%ii + 1, j, layer) = cff * (-3.5D0 * psi%data(decomp%ii + 1, j, layer) &
                                                      & + 4.0D0 * psi%data(decomp%ii, j, layer) &
                                                      & - 0.5D0 * psi%data(decomp%ii - 1, j, layer))
            end do
          end do
        end if
        
        ! Left boundary condition
        if(decomp%p_c == -1) then
          do layer = 1, layers
            do i = decomp%fa_i_min, decomp%fa_i_max
              rel%data(i, 2, layer) = cff * (-3.5D0 * psi%data(i, 2, layer) &
                                          & + 4.0D0 * psi%data(i, 3, layer) &
                                          & - 0.5D0 * psi%data(i, 4, layer))
            end do
          end do
        end if
        
        ! Right boundary condition
        if(decomp%p_a == -1) then
          do layer = 1, layers
            do i = decomp%fa_i_min, decomp%fa_i_max
               rel%data(i, decomp%jj + 1, layer) = cff * (-3.5D0 * psi%data(i, decomp%jj + 1, layer) &
                                                       & + 4.0D0 * psi%data(i, decomp%jj, layer) &
                                                       & - 0.5D0 * psi%data(i, decomp%jj - 1, layer))
            end do
          end do
        end if
      case(BC_NO_SLIP)
        ! Top boundary condition
        if(decomp%p_b == -1) then
          do layer = 1, layers
            do j = decomp%fa_j_min, decomp%fa_j_max
              rel%data(2, j, layer) = 2.0D0 * (psi%data(3, j, layer) - psi%data(2, j, layer))
            end do
          end do
        end if
        
        ! Bottom boundary condition
        if(decomp%p_d == -1) then
          do layer = 1, layers
            do j = decomp%fa_j_min, decomp%fa_j_max
              rel%data(decomp%ii + 1, j, layer) = 2.0D0 * (psi%data(decomp%ii, j, layer) - psi%data(decomp%ii + 1, j, layer))
            end do
          end do
        end if
        
        ! Left boundary condition
        if(decomp%p_c == -1) then
          do layer = 1, layers
            do i = decomp%fa_i_min, decomp%fa_i_max
              rel%data(i, 2, layer) = 2.0D0 * (psi%data(i, 3, layer)  - psi%data(i, 2, layer))
            end do
          end do
        end if
        
        ! Right boundary condition
        if(decomp%p_a == -1) then
          do layer = 1, layers
            do i = decomp%fa_i_min, decomp%fa_i_max
               rel%data(i, decomp%jj + 1, layer) = 2.0D0 * (psi%data(i, decomp%jj, layer) - psi%data(i, decomp%jj + 1, layer))
            end do
          end do
        end if
      case(BC_FREE_SLIP)
        ! Top boundary condition
        if(decomp%p_b == -1) then
          do layer = 1, layers
            do j = decomp%fa_j_min, decomp%fa_j_max
              rel%data(2, j, layer) = 0.0D0
            end do
          end do
        end if
        
        ! Bottom boundary condition
        if(decomp%p_d == -1) then
          do layer = 1, layers
            do j = decomp%fa_j_min, decomp%fa_j_max
              rel%data(decomp%ii + 1, j, layer) = 0.0D0
            end do
          end do
        end if
        
        ! Left boundary condition
        if(decomp%p_c == -1) then
          do layer = 1, layers
            do i = decomp%fa_i_min, decomp%fa_i_max
              rel%data(i, 2, layer) = 0.0D0
            end do
          end do
        end if
        
        ! Right boundary condition
        if(decomp%p_a == -1) then
          do layer = 1, layers
            do i = decomp%fa_i_min, decomp%fa_i_max
               rel%data(i, decomp%jj + 1, layer) = 0.0D0
            end do
          end do
        end if
      case default
        write(0, "(a,i0)") "For boundary condition type: ", vbc%bc_type
        labort("Unrecognised boundary condition type")
    end select
    
  end subroutine apply_rel_bc
  
#ifdef CABARET
  ! Compute the CABARET Courant number. This performs a collective max
  ! operation.
  function cabaret_cfl(dt, ui, uj) result(cfl)
    real(kind = real_8), intent(in) :: dt
    type(parallel_flux_data), intent(in) :: ui
    type(parallel_flux_data), intent(in) :: uj
  
    real(kind = real_8) :: cfl
    
    type(decomposition), pointer :: decomp
    
    assert(dt > 0.0D0)
    assert(ui%flux_dir == FLUX_II)
    assert(associated(uj%decomp, ui%decomp))
    assert(uj%flux_dir == FLUX_JJ)
    
    decomp => ui%decomp
  
    cfl = max( &
      &   maxval(dabs(ui%data(decomp%fa_i_min + 1:decomp%fa_i_max, decomp%fa_j_min + 1:decomp%fa_j_max - 1, :))), &
      &   maxval(dabs(uj%data(decomp%fa_i_min + 1:decomp%fa_i_max - 1, decomp%fa_j_min + 1:decomp%fa_j_max, :))) &
      & ) * dt    
    call allmax(cfl)
  
  end function cabaret_cfl
  
  ! Adapt the timestep according to the CABARET Courant number. This performs a
  ! collective max operation.
  subroutine cabaret_adapt_dt(target_cfl, max_dt, max_dt_increase, ui, uj, dt, tweight_a, tweight_u)  
    real(kind = real_8), intent(in) :: target_cfl
    real(kind = real_8), intent(in) :: max_dt
    real(kind = real_8), intent(in) :: max_dt_increase
    type(parallel_flux_data), intent(in) :: ui
    type(parallel_flux_data), intent(in) :: uj
    real(kind = real_8), intent(inout) :: dt
    real(kind = real_8), intent(out) :: tweight_a
    real(kind = real_8), intent(out) :: tweight_u
    
    real(kind = real_8) :: old_dt, max_u
    type(decomposition), pointer :: decomp
    
    assert(ui%flux_dir == FLUX_II)
    assert(associated(uj%decomp, ui%decomp))
    assert(uj%flux_dir == FLUX_JJ)
    
    decomp => ui%decomp
    
    max_u = max( &
      &   maxval(dabs(ui%data(decomp%fa_i_min + 1:decomp%fa_i_max, decomp%fa_j_min + 1:decomp%fa_j_max - 1, :))), &
      &   maxval(dabs(uj%data(decomp%fa_i_min + 1:decomp%fa_i_max - 1, decomp%fa_j_min + 1:decomp%fa_j_max, :)))&
      & )
    call allmax(max_u)
    
    old_dt = dt
    if(fne(max_u, 0.0D0)) then  ! Do not adapt the time step if the max Courant
                                ! no. is zero
      dt = target_cfl / max_u
      dt = min(dt, old_dt * max_dt_increase, max_dt)
    end if
    
    tweight_a = 0.5D0 * dt / old_dt
    tweight_u = dt / (old_dt + dt)
    
    assert(dt > 0.0D0)
    if(dt < 1.0D-6) then
      write(6, "(a,"//real_8_chr//")") "dt = ", dt
      write(6, "(a)") "WARNING: Small time step encountered"
      write(6, "(a)") "Simulation may be blowing up"
    end if
    assert(tweight_a >= 0.0D0)
    assert(tweight_u >= 0.0D0 .and. tweight_u <= 1.0D0)
    
  end subroutine cabaret_adapt_dt
  
  ! CABARET predictor step. This does *not* update the halo region.
  subroutine cabaret_predictor(dt, zeta_old, res, zeta)
    real(kind = real_8), intent(in) :: dt
    real(kind = real_8), dimension(:, :, :), intent(in) :: zeta_old
    real(kind = real_8), dimension(:, :, :), intent(in) :: res
    type(parallel_data), intent(inout) :: zeta
  
    integer :: i, j, layer
    type(decomposition), pointer :: decomp
    
    assert(dt > 0.0D0)
    assert(size(zeta_old, 1) >= zeta%decomp%ii + 2)
    assert(size(zeta_old, 2) >= zeta%decomp%jj + 2)
    assert(size(zeta_old, 3) == zeta%decomp%layers)
    assert(size(res, 1) >= zeta%decomp%ii)
    assert(size(res, 2) >= zeta%decomp%jj)
    assert(size(res, 3) == zeta%decomp%layers)
    
    decomp => zeta%decomp
  
    do layer = 1, decomp%layers
      do j = decomp%f_j_min, decomp%f_j_max
        do i = decomp%f_i_min, decomp%f_i_max
          zeta%data(i, j, layer) = zeta_old(i, j, layer) + 0.5D0 * dt * res(i - 1, j - 1, layer)
        end do
      end do
    end do
    
  end subroutine cabaret_predictor

  ! Time integration of the beta term. This does *not* update the halo region.
  subroutine cabaret_beta(dt, tweight_a, beta, ui, zeta, res_beta_o)    
    real(kind = real_8), intent(in) :: dt
    real(kind = real_8), intent(in) :: tweight_a
    type(parallel_data), intent(inout) :: zeta
    real(kind = real_8), dimension(zeta%decomp%layers), intent(in) :: beta
    type(parallel_flux_data), intent(in) :: ui
    real(kind = real_8), dimension(:, :, :), intent(inout) :: res_beta_o
  
    integer :: i, j, layer, layers
    real(kind = real_8) :: res_beta
    type(decomposition), pointer :: decomp    
    
    assert(dt > 0.0D0)
    assert(tweight_a >= 0.0D0)
    assert(associated(ui%decomp, zeta%decomp))
    assert(ui%flux_dir == FLUX_II)
    assert(size(res_beta_o, 1) >= zeta%decomp%ii)
    assert(size(res_beta_o, 2) >= zeta%decomp%jj)
    assert(size(res_beta_o, 3) == zeta%decomp%layers)
    
    decomp => zeta%decomp
    layers = decomp%layers
    
    ! Second-order in time integration, which supresses spurious beta_y-waves;
    ! NB: Interstingly, the simple operator
    ! zeta(n + 1) - zeta(n + 1/2) = L(zeta(n + 1/2)), as suggested by VM, also
    ! leads to a stable solution but it is contaminated by spurious
    ! large-scale beta_y-waves.  
    do layer = 1, layers
      do j = decomp%l_j_min, decomp%l_j_max
        do i = decomp%l_i_min, decomp%l_i_max
          res_beta = dt * (0.5D0 * (ui%data(i + 2, j + 1, layer) + ui%data(i + 1, j + 1, layer)) * beta(layer))
          zeta%sdata(i, j, layer) = zeta%sdata(i, j, layer) + (1.0D0 + tweight_a) * res_beta - tweight_a * res_beta_o(i, j, layer)
          res_beta_o(i, j, layer) = res_beta
        end do
      end do
    end do
    
  end subroutine cabaret_beta
  
  ! Time integration of all other terms (viscosity, forcing and bottom
  ! friction). This does *not* update the halo region.
  subroutine cabaret_forcing(dt, &
#ifdef WIND_FORCING
    & force_wind, &
#endif
    & viscosity, visc_bot, rel_p, zeta)
    use qg_parameters, only : ii_m, jj_m, layers
  
    real(kind = real_8), intent(in) :: dt
    real(kind = real_8), intent(in) :: visc_bot
#ifdef WIND_FORCING
    real(kind = real_8), dimension(:, :), intent(in) :: force_wind
#endif
    real(kind = real_8), intent(in) :: viscosity
    type(parallel_data), intent(in) :: rel_p
    type(parallel_data), intent(inout) :: zeta
  
    integer :: i, j, layer!, layers
    type(decomposition), pointer :: decomp
    
    real(kind = real_8), dimension(ii_m, jj_m, layers), save :: res_visc = 0.0D0
  
    assert(dt > 0.0D0)
#ifdef WIND_FORCING
    assert(size(force_wind, 1) >= zeta%decomp%ii)
    assert(size(force_wind, 2) >= zeta%decomp%jj)
#endif
    assert(associated(rel_p%decomp, zeta%decomp))
  
#ifdef DDEBUG
    ! Initialise the CABARET variables with junk to aid invalid access
    ! debugging
    res_visc = INFINITY
#endif
  
    decomp => zeta%decomp
    !layers = decomp%layers    
    if(layers == 1) then
      do j = decomp%f_j_min, decomp%f_j_max
        do i = decomp%f_i_min, decomp%f_i_max
          res_visc(i - 1, j - 1, 1) = dt * &
            & ( &
#ifdef WIND_FORCING
            &   force_wind(i - 1, j - 1) &
#endif
            & - visc_bot * rel_p%data(i, j, 1) &
            & + viscosity * ( &
            &     -4.0D0 * rel_p%data(i, j, 1) + rel_p%data(i - 1, j, 1) + rel_p%data(i + 1, j, 1) + rel_p%data(i, j - 1, 1) + rel_p%data(i, j + 1, 1) &
            &   ) &
            & )
        end do
      end do
    else  
      ! Upper layer
      do j = decomp%f_j_min, decomp%f_j_max
        do i = decomp%f_i_min, decomp%f_i_max
          res_visc(i - 1, j - 1, 1) = dt * &
            & ( &
#ifdef WIND_FORCING
            &   force_wind(i - 1, j - 1) &
#endif
            & + viscosity * ( &
            &     -4.0D0 * rel_p%data(i, j, 1) + rel_p%data(i - 1, j, 1) + rel_p%data(i + 1, j, 1) + rel_p%data(i, j - 1, 1) + rel_p%data(i, j + 1, 1) &
            &   ) &
            & )
        end do
      end do
      ! Middle layers
      do layer = 2, layers - 1
        do j = decomp%f_j_min, decomp%f_j_max
          do i = decomp%f_i_min, decomp%f_i_max
            res_visc(i - 1, j - 1, layer) = dt * &
              &   viscosity * ( &
              &     -4.0D0 * rel_p%data(i, j, layer) + rel_p%data(i - 1, j, layer) + rel_p%data(i + 1, j, layer) + rel_p%data(i, j - 1, layer) + rel_p%data(i, j + 1, layer) &
              &   )
          end do
        end do
      end do
      ! Lower layer
      do j = decomp%f_j_min, decomp%f_j_max
        do i = decomp%f_i_min, decomp%f_i_max
          res_visc(i - 1, j - 1, layers) = dt * &
            & ( &
            & - visc_bot * rel_p%data(i, j, layers) &
            & + viscosity * ( &
            &     -4.0D0 * rel_p%data(i, j, layers) + rel_p%data(i - 1, j, layers) + rel_p%data(i + 1, j, layers) + rel_p%data(i, j - 1, layers) + rel_p%data(i, j + 1, layers) &
            &   ) &
            & )
        end do
      end do
    end if
    
    do layer = 1, layers
      do j = decomp%f_j_min, decomp%f_j_max
        do i = decomp%f_i_min, decomp%f_i_max
          zeta%data(i, j, layer) = zeta%data(i, j, layer) + res_visc(i - 1, j - 1, layer)
        end do
      end do
    end do
        
  end subroutine cabaret_forcing
  
  ! Compute CABARET variables stored on the edges and centres of cells. This
  ! updates the first order halo region only.
  subroutine cabaret_cell_boundaries_velocity(psi, ui, uj &
#ifdef BACKGROUND_VSHEAR
    & , u_background_v &
#endif
    & )
    use qg_parameters, only : ii_m, jj_m, layers
  
    type(parallel_data), intent(in) :: psi
    type(parallel_flux_data), intent(inout) :: ui
    type(parallel_flux_data), intent(inout) :: uj
#ifdef BACKGROUND_VSHEAR
    real(kind = real_8), dimension(:), intent(in) :: u_background_v
#endif
    
    integer :: i, j, layer!, layers
    type(decomposition), pointer :: decomp
    
    real(kind = real_8), dimension(ii_m + 2, jj_m + 2, layers), save :: psi_edge = 0.0D0

    assert(associated(ui%decomp, psi%decomp))
    assert(ui%flux_dir == FLUX_II)
    assert(associated(uj%decomp, psi%decomp))
    assert(uj%flux_dir == FLUX_JJ)
#ifdef BACKGROUND_VSHEAR
    assert(size(u_background_v) == psi%decomp%layers)
#endif
    
    !layers = psi%layers
    decomp => psi%decomp
    
#ifdef DDEBUG
    ! Initialise the CABARET variables with junk to aid invalid access
    ! debugging
    psi_edge = INFINITY
#endif

    do layer = 1, layers
      do j = decomp%fa_j_min + 1, decomp%fa_j_max
        do i = decomp%fa_i_min + 1, decomp%fa_i_max
          psi_edge(i, j, layer) = 0.25D0 * (psi%data(i, j, layer) + psi%data(i, j - 1, layer) + psi%data(i - 1, j, layer) + psi%data(i - 1, j - 1, layer))
        end do
      end do
    end do
    
    do layer = 1, layers
      do j = decomp%fa_j_min + 1, decomp%fa_j_max - 1
        do i = decomp%fa_i_min + 1, decomp%fa_i_max
          ui%data(i, j, layer) = -psi_edge(i, j + 1, layer) + psi_edge(i, j, layer)
        end do
      end do
    end do
    do layer = 1, layers
      do j = decomp%fa_j_min + 1, decomp%fa_j_max
        do i = decomp%fa_i_min + 1, decomp%fa_i_max - 1
          uj%data(i, j, layer) = &
#ifdef BACKGROUND_VSHEAR
            & u_background_v(layer) &
#endif
            & + psi_edge(i + 1, j, layer) - psi_edge(i, j, layer)
        end do
      end do
    end do

  end subroutine cabaret_cell_boundaries_velocity
  
  ! Compute CABARET variables stored on the edges and centres of cells. This
  ! updates the first order halo region only.
  subroutine cabaret_cell_boundaries_time_integrated(tweight_u, psi, ui, uj, uimo, ujmo &
#ifdef BACKGROUND_VSHEAR
    & , u_background_v &
#endif
    & )
    use qg_parameters, only : ii_m, jj_m, layers
  
    real(kind = real_8), intent(in) :: tweight_u
    type(parallel_data), intent(in) :: psi
    type(parallel_flux_data), intent(inout) :: ui
    type(parallel_flux_data), intent(inout) :: uj
    type(parallel_flux_data), intent(inout) :: uimo
    type(parallel_flux_data), intent(inout) :: ujmo
#ifdef BACKGROUND_VSHEAR
    real(kind = real_8), dimension(:), intent(in) :: u_background_v
#endif
    
    integer :: i, j, layer!, layers
    real(kind = real_8) :: uim, ujm
    type(decomposition), pointer :: decomp
    
    real(kind = real_8), dimension(ii_m + 2, jj_m + 2, layers), save :: psi_edge = 0.0D0

    assert(tweight_u >= 0.0D0 .and. tweight_u <= 1.0D0)
    assert(associated(ui%decomp, psi%decomp))
    assert(ui%flux_dir == FLUX_II)
    assert(associated(uj%decomp, psi%decomp))
    assert(uj%flux_dir == FLUX_JJ)
    assert(associated(uimo%decomp, psi%decomp))
    assert(uimo%flux_dir == FLUX_II)
    assert(associated(ujmo%decomp, psi%decomp))
    assert(ujmo%flux_dir == FLUX_JJ)
#ifdef BACKGROUND_VSHEAR
    assert(size(u_background_v) == psi%decomp%layers)
#endif
    
    !layers = psi%layers
    decomp => psi%decomp
    
#ifdef DDEBUG
    ! Initialise the CABARET variables with junk to aid invalid access
    ! debugging
    psi_edge = INFINITY
    call set(ui, INFINITY)
    call set(uj, INFINITY)
#endif

    do layer = 1, layers
      do j = decomp%fa_j_min + 1, decomp%fa_j_max
        do i = decomp%fa_i_min + 1, decomp%fa_i_max
          psi_edge(i, j, layer) = 0.25D0 * (psi%data(i, j, layer) + psi%data(i, j - 1, layer) + psi%data(i - 1, j, layer) + psi%data(i - 1, j - 1, layer))
        end do
      end do
    end do
    
    do layer = 1, layers
      do j = decomp%fa_j_min + 1, decomp%fa_j_max - 1
        do i = decomp%fa_i_min + 1, decomp%fa_i_max
          uim = -psi_edge(i, j + 1, layer) + psi_edge(i, j, layer)
          ui%data(i, j, layer) = (1.0D0 + tweight_u) * uim - tweight_u * uimo%data(i, j, layer)
          uimo%data(i, j, layer) = uim
        end do
      end do
    end do
    do layer = 1, layers
      do j = decomp%fa_j_min + 1, decomp%fa_j_max
        do i = decomp%fa_i_min + 1, decomp%fa_i_max - 1
          ujm = &
#ifdef BACKGROUND_VSHEAR
            & u_background_v(layer) &
#endif
            & + psi_edge(i + 1, j, layer) - psi_edge(i, j, layer)
          uj%data(i, j, layer) = (1.0D0 + tweight_u) * ujm - tweight_u * ujmo%data(i, j, layer)
          ujmo%data(i, j, layer) = ujm
        end do
      end do
    end do

  end subroutine cabaret_cell_boundaries_time_integrated
    
  ! CABARET extrapolator step, including the flux limiter. This updates the
  ! first order halo region only.
  subroutine cabaret_extrapolator(dt, zeta, zeta_old, ui, uj, zeta_fi, zeta_fj)
    use qg_parameters, only : ii_m, jj_m, layers
  
    real(kind = real_8), intent(in) :: dt
    type(parallel_data), intent(in) :: zeta
    real(kind = real_8), dimension(:, :, :), intent(in) :: zeta_old
    type(parallel_flux_data), intent(in) :: ui
    type(parallel_flux_data), intent(in) :: uj
    type(parallel_flux_data), intent(inout) :: zeta_fi
    type(parallel_flux_data), intent(inout) :: zeta_fj
    
    integer :: i, j, layer!, layers
    type(decomposition), pointer :: decomp
    
    real(kind = real_8) :: fi2, psi_n, u1
#ifdef CABARET_LIMITER
    real(kind = real_8) :: d_f, fi1, fmax, fmin, psi_o, u2, v0
#endif
    
    real(kind = real_8), dimension(ii_m + 2, jj_m + 2, layers), save :: zeta_fin = 0.0D0
    real(kind = real_8), dimension(ii_m + 2, jj_m + 2, layers), save :: zeta_fjn = 0.0D0

    assert(dt > 0.0D0)
    assert(size(zeta_old, 1) >= zeta%decomp%ii + 2)
    assert(size(zeta_old, 2) >= zeta%decomp%jj + 2)
    assert(size(zeta_old, 3) == layers)
    assert(associated(ui%decomp, zeta%decomp))
    assert(ui%flux_dir == FLUX_II)
    assert(associated(uj%decomp, zeta%decomp))
    assert(uj%flux_dir == FLUX_JJ)
    assert(associated(zeta_fi%decomp, zeta%decomp))
    assert(zeta_fi%flux_dir == FLUX_II)
    assert(associated(zeta_fj%decomp, zeta%decomp))
    assert(zeta_fj%flux_dir == FLUX_JJ)

    !layers = zeta%layers
    decomp => zeta%decomp
    
#ifdef DDEBUG
    ! Initialise the CABARET variables with junk to aid invalid access
    ! debugging
    zeta_fin = INFINITY
    zeta_fjn = INFINITY
#endif
    
    do layer = 1, layers
      do j = decomp%fa_j_min + 1, decomp%fa_j_max - 1
        do i = decomp%fa_i_min + 1, decomp%fa_i_max   
          u1 = ui%data(i, j, layer)
          
#ifdef CABARET_LIMITER   
          fi1 = zeta_fi%data(i, j, layer)
          
          if(u1 >= 0.0D0) then
            ! Positive speed
            if(decomp%p_b >= 0 .or. i > decomp%fb_i_min + 1) then
              u2 = ui%data(i - 1, j, layer)
              fi2 = zeta_fi%data(i - 1, j, layer)
              v0 = (u1 + u2) / 2.0D0
              psi_n = zeta%data(i - 1, j, layer)
              psi_o = zeta_old(i - 1, j, layer)
              d_f = -(psi_n - psi_o) * 2.0D0 - dt * v0 * (fi1 - fi2)
            else
              ! Boundary conditions
              fi2 = zeta%data(i - 1, j, layer)
              v0 = 0.0D0
              psi_n = zeta%data(i - 1, j, layer)
              psi_o = zeta_old(i - 1, j, layer)
              d_f = -(psi_n - psi_o)
            end if
          else
            ! Negative speed
            if(decomp%p_d >= 0 .or. i < decomp%fb_i_max) then
              u2 = ui%data(i + 1, j, layer)
              fi2 = zeta_fi%data(i + 1, j, layer)
              v0 = (u1 + u2) / 2.0D0
              psi_n = zeta%data(i, j, layer)
              psi_o = zeta_old(i, j, layer)
              d_f = -(psi_n - psi_o) * 2.0D0 - dt * v0 * (fi2 - fi1)
            else
              ! Boundary conditions
              fi2 = zeta%data(i, j, layer)
              v0 = 0.0D0
              psi_n = zeta%data(i, j, layer)
              psi_o = zeta_old(i, j, layer)
              d_f = -(psi_n - psi_o)
            end if
          end if
#else
          if(u1 >= 0.0D0) then
            ! Positive speed
            if(decomp%p_b >= 0 .or. i > decomp%fb_i_min + 1) then
              fi2 = zeta_fi%data(i - 1, j, layer)
              psi_n = zeta%data(i - 1, j, layer)
            else
              ! Boundary conditions
              fi2 = zeta%data(i - 1, j, layer)
              psi_n = zeta%data(i - 1, j, layer)
            end if
          else
            ! Negative speed
            if(decomp%p_d >= 0 .or. i < decomp%fb_i_max) then
              fi2 = zeta_fi%data(i + 1, j, layer)
              psi_n = zeta%data(i, j, layer)
            else
              ! Boundary conditions
              fi2 = zeta%data(i, j, layer)
              psi_n = zeta%data(i, j, layer)
            end if
          end if
#endif
            
          zeta_fin(i, j, layer) = 2.0D0 * psi_n - fi2
            
#ifdef CABARET_LIMITER
          ! Non-linear correction to enforce the maximum principle
          fmax = max(fi1, fi2, psi_n) - d_f
          fmin = min(fi1, fi2, psi_n) - d_f
          if(zeta_fin(i, j, layer) > fmax) zeta_fin(i, j, layer) = fmax
          if(zeta_fin(i, j, layer) < fmin) zeta_fin(i, j, layer) = fmin
#endif
        end do
      end do
    end do
    
    do layer = 1, layers
      do j = decomp%fa_j_min + 1, decomp%fa_j_max
        do i = decomp%fa_i_min + 1, decomp%fa_i_max - 1  
          u1 = uj%data(i, j, layer)
          
#ifdef CABARET_LIMITER        
          fi1 = zeta_fj%data(i, j, layer)
          
          if(u1 >= 0.0D0) then
            ! Positive speed
            if(decomp%p_c >= 0 .or. j > decomp%fb_j_min + 1) then
              u2 = uj%data(i, j - 1, layer)
              fi2 = zeta_fj%data(i, j - 1, layer)
              v0 = (u1 + u2) / 2.0D0
              psi_n = zeta%data(i, j - 1, layer)
              psi_o = zeta_old(i, j - 1, layer)
              d_f = -(psi_n - psi_o) * 2.0D0 - dt * v0 * (fi1 - fi2)
            else
              ! Boundary conditions
              fi2 = zeta%data(i, j - 1, layer)
              v0 = 0.0D0
              psi_n = zeta%data(i, j - 1, layer)
              psi_o = zeta_old(i, j - 1, layer)
              d_f = -(psi_n - psi_o)
            end if
          else
            ! Negative speed
            if(decomp%p_a >= 0 .or. j < decomp%fb_j_max) then
              u2 = uj%data(i, j + 1, layer)
              fi2 = zeta_fj%data(i, j + 1, layer)
              v0 = (u1 + u2) / 2.0D0
              psi_n = zeta%data(i, j, layer)
              psi_o = zeta_old(i, j, layer)
              d_f = -(psi_n - psi_o) * 2.0D0 - dt * v0 * (fi2 - fi1)
            else
              ! Boundary conditions
              fi2 = zeta%data(i, j, layer)
              v0 = 0.0D0
              psi_n = zeta%data(i, j, layer)
              psi_o = zeta_old(i, j, layer)
              d_f = -(psi_n - psi_o)
            end if             
          end if
#else
          if(u1 >= 0.0D0) then
            ! Positive speed
            if(decomp%p_c >= 0 .or. j > decomp%fb_j_min + 1) then
              fi2 = zeta_fj%data(i, j - 1, layer)
              psi_n = zeta%data(i, j - 1, layer)
            else
              ! Boundary conditions
              fi2 = zeta%data(i, j - 1, layer)
              psi_n = zeta%data(i, j - 1, layer)
            end if
          else
            ! Negative speed
            if(decomp%p_a >= 0 .or. j < decomp%fb_j_max) then
              fi2 = zeta_fj%data(i, j + 1, layer)
              psi_n = zeta%data(i, j, layer)
            else
              ! Boundary conditions
              fi2 = zeta%data(i, j, layer)
              psi_n = zeta%data(i, j, layer)
            end if             
          end if
#endif
          
          zeta_fjn(i, j, layer) = 2.0D0 * psi_n - fi2
          
#ifdef CABARET_LIMITER
          ! Non-linear correction to enforce the maximum principle
          fmax = max(fi1, fi2, psi_n) - d_f
          fmin = min(fi1, fi2, psi_n) - d_f
          if(zeta_fjn(i, j, layer) > fmax) zeta_fjn(i, j, layer) = fmax
          if(zeta_fjn(i, j, layer) < fmin) zeta_fjn(i, j, layer) = fmin
#endif          
        end do
      end do
    end do

#ifdef DDEBUG
    ! Initialise the CABARET variables with junk to aid invalid access
    ! debugging
    call set(zeta_fi, INFINITY)
    call set(zeta_fj, INFINITY)
#endif
    
    do layer = 1, layers
      do j = decomp%fa_j_min + 1, decomp%fa_j_max - 1
        do i = decomp%fa_i_min + 1, decomp%fa_i_max 
          zeta_fi%data(i, j, layer) = zeta_fin(i, j, layer)
        end do
      end do
    end do
    do layer = 1, layers
      do j = decomp%fa_j_min + 1, decomp%fa_j_max
        do i = decomp%fa_i_min + 1, decomp%fa_i_max - 1   
          zeta_fj%data(i, j, layer) = zeta_fjn(i, j, layer)
        end do
      end do
    end do
    
  end subroutine cabaret_extrapolator
  
  ! Compute the CABARET flux. This updates the first order halo region only.
  subroutine cabaret_flux(ui, uj, zeta_fi, zeta_fj, flux_i, flux_j)
    type(parallel_flux_data), intent(in) :: ui
    type(parallel_flux_data), intent(in) :: uj
    type(parallel_flux_data), intent(in) :: zeta_fi
    type(parallel_flux_data), intent(in) :: zeta_fj
    type(parallel_flux_data), intent(inout) :: flux_i
    type(parallel_flux_data), intent(inout) :: flux_j
    
    integer :: i, j, layer
    type(decomposition), pointer :: decomp
    
    assert(ui%flux_dir == FLUX_II)
    assert(associated(uj%decomp, ui%decomp))
    assert(uj%flux_dir == FLUX_JJ)
    assert(associated(zeta_fi%decomp, ui%decomp))
    assert(zeta_fi%flux_dir == FLUX_II)
    assert(associated(zeta_fj%decomp, ui%decomp))
    assert(zeta_fj%flux_dir == FLUX_JJ)
    assert(associated(flux_i%decomp, ui%decomp))
    assert(flux_i%flux_dir == FLUX_II)
    assert(associated(flux_j%decomp, ui%decomp))
    assert(flux_j%flux_dir == FLUX_JJ)
    
    decomp => ui%decomp
    
#ifdef DDEBUG
    ! Initialise the CABARET variables with junk to aid invalid access
    ! debugging
    call set(flux_i, INFINITY)
    call set(flux_j, INFINITY)
#endif
    
    do layer = 1, decomp%layers
      do j = decomp%fa_j_min + 1, decomp%fa_j_max - 1
        do i = decomp%fa_i_min, decomp%fa_i_max + 1
          flux_i%data(i, j, layer) = ui%data(i, j, layer) * zeta_fi%data(i, j, layer)
        end do
      end do
    end do
    do layer = 1, decomp%layers
      do j = decomp%fa_j_min, decomp%fa_j_max + 1
        do i = decomp%fa_i_min + 1, decomp%fa_i_max - 1   
          flux_j%data(i, j, layer) = uj%data(i, j, layer) * zeta_fj%data(i, j, layer)
        end do
      end do
    end do
    
  end subroutine cabaret_flux
  
  ! CABARET corrector step. This does *not* update the halo region.
  subroutine cabaret_corrector(dt, flux_i, flux_j, res, zeta)
    real(kind = real_8), intent(in) :: dt
    type(parallel_data), intent(inout) :: zeta
    type(parallel_flux_data), intent(in) :: flux_i
    type(parallel_flux_data), intent(in) :: flux_j
    real(kind = real_8), dimension(:, :, :), intent(inout) :: res
    
    integer :: i, j, layer
    type(decomposition), pointer :: decomp
    
    assert(dt > 0.0D0)
    assert(associated(flux_i%decomp, zeta%decomp))
    assert(flux_i%flux_dir == FLUX_II)
    assert(associated(flux_j%decomp, zeta%decomp))
    assert(flux_j%flux_dir == FLUX_JJ)
    assert(size(res, 1) >= zeta%decomp%ii)
    assert(size(res, 2) >= zeta%decomp%jj)
    assert(size(res, 3) == zeta%decomp%layers)
    
    decomp => zeta%decomp
    
#ifdef DDEBUG
    ! Initialise the CABARET variables with junk to aid invalid access
    ! debugging
    res = INFINITY
#endif

    do layer = 1, decomp%layers
      do j = decomp%f_j_min, decomp%f_j_max
        do i = decomp%f_i_min, decomp%f_i_max
          res(i - 1, j - 1, layer) = -(flux_i%data(i + 1, j, layer) - flux_i%data(i, j, layer) + flux_j%data(i, j + 1, layer) - flux_j%data(i, j, layer))
#ifdef DDEBUG
          if(dabs(res(i - 1, j - 1, layer)) > 1.0D100 .or. res(i - 1, j - 1, layer) /= res(i - 1, j - 1, layer)) then
            write(6, "(a,i0,a,i0)") "For i = ", i, ", j = ", j
            write(6, "(a,i0,a,"//real_8_chr//")") "Layer ", layer, " residual: ", res(i - 1, j - 1, layer)
            write(6, "(a,"//real_8_chr//",a,"//real_8_chr//",a,"//real_8_chr//",a,"//real_8_chr//")") "Fluxes: ", flux_i%data(i + 1, j, layer), ", ", flux_i%data(i, j, layer), ", ", flux_j%data(i, j + 1, layer), ", ", flux_j%data(i, j, layer)
            labort("Large or NaN residual")
          end if
#endif
          
          zeta%data(i, j, layer) = zeta%data(i, j, layer) + 0.5D0 * dt * res(i - 1, j - 1, layer)
        end do
      end do
    end do
    
  end subroutine cabaret_corrector
  
#else

  ! Apply the Robert-Asselin-Williams time filter.
  subroutine robert_asselin_filter_rank1(phi_l, phi_m, phi_old_l, phi_old_m)
    real(kind = real_8), dimension(:), intent(inout) :: phi_l
    real(kind = real_8), dimension(:), intent(inout) :: phi_m
    real(kind = real_8), dimension(:), intent(in) :: phi_old_l
    real(kind = real_8), dimension(:), intent(out) :: phi_old_m
    
    integer :: i
#ifdef RA_DOUBLE_FILTER
    real(kind = real_8) :: lphi_l
#endif
    
    assert(size(phi_m) == size(phi_l))
    assert(size(phi_old_l) == size(phi_l))
    assert(size(phi_old_m) == size(phi_l))
    
    do i = 1, size(phi_l)
      ! phi_old_l  \phi_{n - 1}
      ! phi_m      \phi_n
      ! phi_l      \phi_{n + 1}
#ifdef RA_DOUBLE_FILTER       
      lphi_l = phi_l(i)
      phi_l(i) = ra_weight11 * phi_old_l(i) + ra_weight12 * phi_m(i) + ra_weight13 * lphi_l
      phi_m(i) = ra_weight21 * (phi_old_l(i) + lphi_l) + ra_weight22 * phi_m(i)
#else
      phi_m(i) = ra_weight1 * (phi_old_l(i) + phi_l(i)) + ra_weight2 * phi_m(i)
#endif
      phi_old_m = phi_m(i)
    end do
    
  end subroutine robert_asselin_filter_rank1
    
  ! Apply the Robert-Asselin-Williams time filter.
  subroutine robert_asselin_filter_rank2(phi_l, phi_m, phi_old_l, phi_old_m)
    real(kind = real_8), dimension(:, :), intent(inout) :: phi_l
    real(kind = real_8), dimension(:, :), intent(inout) :: phi_m
    real(kind = real_8), dimension(:, :), intent(in) :: phi_old_l
    real(kind = real_8), dimension(:, :), intent(out) :: phi_old_m
    
    integer :: i, j
#ifdef RA_DOUBLE_FILTER
    real(kind = real_8) :: lphi_l
#endif
    
    assert(size(phi_m, 1) == size(phi_l, 1))
    assert(size(phi_m, 2) == size(phi_l, 2))
    assert(size(phi_old_l, 1) == size(phi_l, 1))
    assert(size(phi_old_l, 2) == size(phi_l, 2))
    assert(size(phi_old_m, 1) == size(phi_l, 1))
    assert(size(phi_old_m, 2) == size(phi_l, 2))
    
    do j = 1, size(phi_l, 2)
      do i = 1, size(phi_l, 1)
        ! phi_old_l  \phi_{n - 1}
        ! phi_m      \phi_n
        ! phi_l      \phi_{n + 1}
#ifdef RA_DOUBLE_FILTER       
        lphi_l = phi_l(i, j)
        phi_l(i, j) = ra_weight11 * phi_old_l(i, j) + ra_weight12 * phi_m(i, j) + ra_weight13 * lphi_l
        phi_m(i, j) = ra_weight21 * (phi_old_l(i, j) + lphi_l) + ra_weight22 * phi_m(i, j)
#else
        phi_m(i, j) = ra_weight1 * (phi_old_l(i, j) + phi_l(i, j)) + ra_weight2 * phi_m(i, j)
#endif
        phi_old_m = phi_m(i, j)
      end do
    end do
    
  end subroutine robert_asselin_filter_rank2
    
  ! Apply the Robert-Asselin-Williams time filter.
  subroutine robert_asselin_filter_rank3(phi_l, phi_m, phi_old_l, phi_old_m)
    real(kind = real_8), dimension(:, :, :), intent(inout) :: phi_l
    real(kind = real_8), dimension(:, :, :), intent(inout) :: phi_m
    real(kind = real_8), dimension(:, :, :), intent(in) :: phi_old_l
    real(kind = real_8), dimension(:, :, :), intent(out) :: phi_old_m
    
    integer :: i, j, k
#ifdef RA_DOUBLE_FILTER
    real(kind = real_8) :: lphi_l
#endif

    assert(size(phi_m, 1) == size(phi_l, 1))
    assert(size(phi_m, 2) == size(phi_l, 2))
    assert(size(phi_m, 3) == size(phi_l, 3))
    assert(size(phi_old_l, 1) == size(phi_l, 1))
    assert(size(phi_old_l, 2) == size(phi_l, 2))
    assert(size(phi_old_l, 3) == size(phi_l, 3))
    assert(size(phi_old_m, 1) == size(phi_l, 1))
    assert(size(phi_old_m, 2) == size(phi_l, 2))
    assert(size(phi_old_m, 3) == size(phi_l, 3))
    
    do k = 1, size(phi_l, 3)
      do j = 1, size(phi_l, 2)
        do i = 1, size(phi_l, 1)
          ! phi_old_l  \phi_{n - 1}
          ! phi_m      \phi_n
          ! phi_l      \phi_{n + 1}
#ifdef RA_DOUBLE_FILTER
          lphi_l = phi_l(i, j, k)
          phi_l(i, j, k) = ra_weight11 * phi_old_l(i, j, k) + ra_weight12 * phi_m(i, j, k) + ra_weight13 * lphi_l
          phi_m(i, j, k) = ra_weight21 * (phi_old_l(i, j, k) + lphi_l) + ra_weight22 * phi_m(i, j, k)
#else
          phi_m(i, j, k) = ra_weight1 * (phi_old_l(i, j, k) + phi_l(i, j, k)) + ra_weight2 * phi_m(i, j, k)
#endif
          phi_old_m(i, j, k) = phi_m(i, j, k)
        end do
      end do
    end do
    
  end subroutine robert_asselin_filter_rank3

  ! Apply the Robert-Asselin-Williams time filter. This updates the halo region.
  subroutine robert_asselin_filter_pdata(leapfrog, phi, phi_old)
    use qg_parameters, only : ra_nu
#ifdef RA_DOUBLE_FILTER
    use qg_parameters, only : ra_alpha
#endif
    
    integer, intent(in) :: leapfrog
    type(parallel_data), dimension(2), intent(inout) :: phi
    real(kind = real_8), dimension(:, :, :, :), intent(inout) :: phi_old
    
    integer :: m
    type(decomposition), pointer :: decomp
    
    assert(associated(phi(1)%decomp, phi(2)%decomp))
    assert(size(phi_old, 1) >= phi(2)%decomp%ii + 2)
    assert(size(phi_old, 2) >= phi(2)%decomp%jj + 2)
    assert(size(phi_old, 3) == phi(2)%decomp%layers)
    assert(size(phi_old, 4) == 2)
    
    decomp => phi(2)%decomp
    
    m = mod(leapfrog, 2) + 1
    
    call robert_asselin_filter(phi(leapfrog)%data(decomp%fa_i_min:decomp%fa_i_max, decomp%fa_j_min:decomp%fa_j_max, :), &
                             & phi(m)%data(decomp%fa_i_min:decomp%fa_i_max, decomp%fa_j_min:decomp%fa_j_max, :), &
                             & phi_old(decomp%fa_i_min:decomp%fa_i_max, decomp%fa_j_min:decomp%fa_j_max, :, leapfrog), &
                             & phi_old(decomp%fa_i_min:decomp%fa_i_max, decomp%fa_j_min:decomp%fa_j_max, :, m))
    
  end subroutine robert_asselin_filter_pdata
  
  ! Compute the relative potential vorticity tendency. This does *not* update
  ! the halo region.
  subroutine zeta_tendency(dt, strat, beta_y, &
#ifdef WIND_FORCING
    & force_wind, &
#endif
#ifdef BACKGROUND_VSHEAR
    & u_background_v, &
#endif
    & viscosity, visc_bot, psi, rel, rel_p, zeta_t)
    real(kind = real_8), intent(in) :: dt
    type(stratification), intent(in) :: strat
    type(parallel_data), intent(in) :: psi
    real(kind = real_8), dimension(:, :), intent(in) :: beta_y
#ifdef WIND_FORCING
    real(kind = real_8), dimension(:, :), intent(in) :: force_wind
#endif
#ifdef BACKGROUND_VSHEAR
    real(kind = real_8), dimension(psi%decomp%layers), intent(in) :: u_background_v
#endif
    real(kind = real_8), intent(in) :: viscosity
    real(kind = real_8), intent(in) :: visc_bot
    type(parallel_data), intent(in) :: rel
    type(parallel_data), intent(in) :: rel_p
    type(parallel_data), intent(inout) :: zeta_t
    
    integer :: i, j, layer, layers
    type(decomposition), pointer :: decomp
    
    real(kind = real_8), dimension(:, :, :), pointer :: lpsi, lrel, lrel_p, &
      & lzeta_t
      
    assert(associated(rel%decomp, psi%decomp))
    assert(associated(rel_p%decomp, psi%decomp))
    assert(associated(zeta_t%decomp, psi%decomp))
    
    decomp => psi%decomp
    layers = decomp%layers
    
    lpsi => psi%data
    lrel => rel%data
    lrel_p => rel_p%data
    lzeta_t => zeta_t%data
    
    if(layers == 1) then
      do j = decomp%f_j_min, decomp%f_j_max
        do i = decomp%f_i_min, decomp%f_i_max
          lzeta_t(i, j, 1) = &
            & dt * ( &
            &   0.25D0 * ( &
#ifdef BACKGROUND_VSHEAR
            &       (beta_y(i, 1) + lrel(i, j + 1, 1)) * (lpsi(i - 1, j + 1, 1) - lpsi(i + 1, j + 1, 1) - 2.0D0 * u_background_v(1)) &
            &     - (beta_y(i, 1) + lrel(i, j - 1, 1)) * (lpsi(i - 1, j - 1, 1) - lpsi(i + 1, j - 1, 1) - 2.0D0 * u_background_v(1)) &
#else
            &       (beta_y(i, 1) + lrel(i, j + 1, 1)) * (lpsi(i - 1, j + 1, 1) - lpsi(i + 1, j + 1, 1)) &
            &     - (beta_y(i, 1) + lrel(i, j - 1, 1)) * (lpsi(i - 1, j - 1, 1) - lpsi(i + 1, j - 1, 1)) &
#endif
            &     - (beta_y(i - 1, 1) + lrel(i - 1, j, 1)) * (lpsi(i - 1, j + 1, 1) - lpsi(i - 1, j - 1, 1)) &
            &     + (beta_y(i + 1, 1) + lrel(i + 1, j, 1)) * (lpsi(i + 1, j + 1, 1) - lpsi(i + 1, j - 1, 1)) &
            &   ) &
#ifdef WIND_FORCING
            &   + force_wind(i - 1, j - 1) &
#endif
            &   + viscosity * ( &
            &     -4.0D0 * lrel_p(i, j, 1) + lrel_p(i - 1, j, 1) + lrel_p(i + 1, j, 1) + lrel_p(i, j - 1, 1) + lrel_p(i, j + 1, 1) &
            &   - visc_bot * lrel_p(i, j, layers) &
            &   ) &
            & )
        end do
      end do
    else
      ! Upper layer
      do j = decomp%f_j_min, decomp%f_j_max
        do i = decomp%f_i_min, decomp%f_i_max
          lzeta_t(i, j, 1) = &
            & dt * ( &
            &   0.25D0 * ( &
#ifdef BACKGROUND_VSHEAR
            &       (beta_y(i, 1) + lrel(i, j + 1, 1) + strat%s2(1) * lpsi(i, j + 1, 2)) * (lpsi(i - 1, j + 1, 1) - lpsi(i + 1, j + 1, 1) - 2.0D0 * u_background_v(1)) &
            &     - (beta_y(i, 1) + lrel(i, j - 1, 1) + strat%s2(1) * lpsi(i, j - 1, 2)) * (lpsi(i - 1, j - 1, 1) - lpsi(i + 1, j - 1, 1) - 2.0D0 * u_background_v(1)) &
            &     + (strat%s2(1) * lpsi(i, j + 1, 1)) * (2.0D0 * u_background_v(1)) &
            &     - (strat%s2(1) * lpsi(i, j - 1, 1)) * (2.0D0 * u_background_v(1)) &
#else
            &       (beta_y(i, 1) + lrel(i, j + 1, 1) + strat%s2(1) * lpsi(i, j + 1, 2)) * (lpsi(i - 1, j + 1, 1) - lpsi(i + 1, j + 1, 1)) &
            &     - (beta_y(i, 1) + lrel(i, j - 1, 1) + strat%s2(1) * lpsi(i, j - 1, 2)) * (lpsi(i - 1, j - 1, 1) - lpsi(i + 1, j - 1, 1)) &
#endif
            &     - (beta_y(i - 1, 1) + lrel(i - 1, j, 1) + strat%s2(1) * lpsi(i - 1, j, 2)) * (lpsi(i - 1, j + 1, 1) - lpsi(i - 1, j - 1, 1)) &
            &     + (beta_y(i + 1, 1) + lrel(i + 1, j, 1) + strat%s2(1) * lpsi(i + 1, j, 2)) * (lpsi(i + 1, j + 1, 1) - lpsi(i + 1, j - 1, 1)) &
            &   ) &
#ifdef WIND_FORCING
            &   + force_wind(i - 1, j - 1) &
#endif
            &   + viscosity * ( &
            &     -4.0D0 * lrel_p(i, j, 1) + lrel_p(i - 1, j, 1) + lrel_p(i + 1, j, 1) + lrel_p(i, j - 1, 1) + lrel_p(i, j + 1, 1) &
            &   ) &
            & )
        end do
      end do
      ! Middle layers
      do layer = 2, layers - 1
        do j = decomp%f_j_min, decomp%f_j_max
          do i = decomp%f_i_min, decomp%f_i_max
            lzeta_t(i, j, layer) = &
              & dt * ( &
              &   0.25D0 * ( &
#ifdef BACKGROUND_VSHEAR
              &       (beta_y(i, layer) + lrel(i, j + 1, layer) + strat%s1(layer) * lpsi(i, j + 1, layer - 1) + strat%s2(layer) * lpsi(i, j + 1, layer + 1)) * (lpsi(i - 1, j + 1, layer) - lpsi(i + 1, j + 1, layer) - 2.0D0 * u_background_v(layer)) &
              &     - (beta_y(i, layer) + lrel(i, j - 1, layer) + strat%s1(layer) * lpsi(i, j - 1, layer - 1) +  strat%s2(layer) * lpsi(i, j - 1, layer + 1)) * (lpsi(i - 1, j - 1, layer) - lpsi(i + 1, j - 1, layer) - 2.0D0 * u_background_v(layer)) &
              &     + (strat%s1(layer) * lpsi(i, j + 1, layer) + strat%s2(layer) * lpsi(i, j + 1, layer)) * (2.0D0 * u_background_v(layer)) &
              &     - (strat%s1(layer) * lpsi(i, j - 1, layer) +  strat%s2(layer) * lpsi(i, j - 1, layer)) * (2.0D0 * u_background_v(layer)) &
#else
              &       (beta_y(i, layer) + lrel(i, j + 1, layer) + strat%s1(layer) * lpsi(i, j + 1, layer - 1) + strat%s2(layer) * lpsi(i, j + 1, layer + 1)) * (lpsi(i - 1, j + 1, layer) - lpsi(i + 1, j + 1, layer)) &
              &     - (beta_y(i, layer) + lrel(i, j - 1, layer) + strat%s1(layer) * lpsi(i, j - 1, layer - 1) +  strat%s2(layer) * lpsi(i, j - 1, layer + 1)) * (lpsi(i - 1, j - 1, layer) - lpsi(i + 1, j - 1, layer)) &
#endif
              &     - (beta_y(i - 1, layer) + lrel(i - 1, j, layer) + strat%s1(layer) * lpsi(i - 1, j, layer - 1) + strat%s2(layer) * lpsi(i - 1, j, layer + 1)) * (lpsi(i - 1, j + 1, layer) - lpsi(i - 1, j - 1, layer)) &
              &     + (beta_y(i + 1, layer) + lrel(i + 1, j, layer) + strat%s1(layer) * lpsi(i + 1, j, layer - 1) + strat%s2(layer) * lpsi(i + 1, j, layer + 1)) * (lpsi(i + 1, j + 1, layer) - lpsi(i + 1, j - 1, layer)) &
              &   ) &
              &   + viscosity * ( &
              &     -4.0D0 * lrel_p(i, j, layer) + lrel_p(i - 1, j, layer) + lrel_p(i + 1, j, layer) + lrel_p(i, j - 1, layer) + lrel_p(i, j + 1, layer) &
              &   ) &
              & )
          end do
        end do
      end do
      ! Lower layer
      do j = decomp%f_j_min, decomp%f_j_max
        do i = decomp%f_i_min, decomp%f_i_max      
          lzeta_t(i, j, layers) = &
            & dt * ( &
            &   0.25D0 * ( &
#ifdef BACKGROUND_VSHEAR
            &       (beta_y(i, layers) + lrel(i, j + 1, layers) + strat%s1(layers) * lpsi(i, j + 1, layers - 1)) * (lpsi(i - 1, j + 1, layers) - lpsi(i + 1, j + 1, layers) - 2.0D0 * u_background_v(layers)) &
            &     - (beta_y(i, layers) + lrel(i, j - 1, layers) + strat%s1(layers) * lpsi(i, j - 1, layers - 1)) * (lpsi(i - 1, j - 1, layers) - lpsi(i + 1, j - 1, layers) - 2.0D0 * u_background_v(layers)) &
            &     + (strat%s1(layers) * lpsi(i, j + 1, layers)) * (2.0D0 * u_background_v(layers)) &
            &     - (strat%s1(layers) * lpsi(i, j - 1, layers)) * (2.0D0 * u_background_v(layers)) &
#else
            &       (beta_y(i, layers) + lrel(i, j + 1, layers) + strat%s1(layers) * lpsi(i, j + 1, layers - 1)) * (lpsi(i - 1, j + 1, layers) - lpsi(i + 1, j + 1, layers)) &
            &     - (beta_y(i, layers) + lrel(i, j - 1, layers) + strat%s1(layers) * lpsi(i, j - 1, layers - 1)) * (lpsi(i - 1, j - 1, layers) - lpsi(i + 1, j - 1, layers)) &
#endif
            &     - (beta_y(i - 1, layers) + lrel(i - 1, j, layers) + strat%s1(layers) * lpsi(i - 1, j, layers - 1)) * (lpsi(i - 1, j + 1, layers) - lpsi(i - 1, j - 1, layers)) &
            &     + (beta_y(i + 1, layers) + lrel(i + 1, j, layers) + strat%s1(layers) * lpsi(i + 1, j, layers - 1)) * (lpsi(i + 1, j + 1, layers) - lpsi(i + 1, j - 1, layers)) &
            &   ) &
            &   + viscosity * ( &
            &     -4.0D0 * lrel_p(i, j, layers) + lrel_p(i - 1, j, layers) + lrel_p(i + 1, j, layers) + lrel_p(i, j - 1, layers) + lrel_p(i, j + 1, layers) &
            &   ) &
            &   - visc_bot * lrel_p(i, j, layers) &
            & )
        end do
      end do
    end if
    
  end subroutine zeta_tendency
  
  ! Timestep the relative potential vorticity using a forward Euler scheme. This
  ! does *not* update the halo region.
  subroutine zeta_forward_euler(dt, strat, beta_y, &
#ifdef WIND_FORCING
    & force_wind, &
#endif
#ifdef BACKGROUND_VSHEAR
    & u_background_v, &
#endif
    & viscosity, visc_bot, psi, rel, rel_p, zeta, zeta_np1, zeta_t)
    real(kind = real_8), intent(in) :: dt
    type(stratification), intent(in) :: strat
    type(parallel_data), intent(in) :: psi
    real(kind = real_8), dimension(:, :), intent(in) :: beta_y
#ifdef WIND_FORCING
    real(kind = real_8), dimension(:, :), intent(in) :: force_wind
#endif
#ifdef BACKGROUND_VSHEAR
    real(kind = real_8), dimension(psi%decomp%layers), intent(in) :: u_background_v
#endif
    real(kind = real_8), intent(in) :: viscosity
    real(kind = real_8), intent(in) :: visc_bot
    type(parallel_data), intent(in) :: rel
    type(parallel_data), intent(in) :: rel_p
    type(parallel_data), target, intent(inout) :: zeta
    type(parallel_data), target, optional, intent(inout) :: zeta_np1
    type(parallel_data), target, optional, intent(inout) :: zeta_t
    
    integer :: i, j, layer
    type(decomposition), pointer :: decomp
    ! The following variable is used (rather than allocating lzeta_t if
    ! required) in order to work around an ICE in PGI 11.6.0
    type(parallel_data), target :: lzeta_t_t
    type(parallel_data), pointer :: lzeta_np1, lzeta_t
    
#ifdef DDEBUG
    assert(associated(rel%decomp, psi%decomp))
    assert(associated(rel_p%decomp, psi%decomp))
    assert(associated(zeta%decomp, psi%decomp))
    if(present(zeta_np1)) then
      assert(associated(zeta_np1%decomp, psi%decomp))
    end if
    if(present(zeta_t)) then
      assert(associated(zeta_t%decomp, psi%decomp))
    end if
#endif
    
    decomp => psi%decomp
    
    if(present(zeta_np1)) then
      lzeta_np1 => zeta_np1
    else
      lzeta_np1 => zeta
    end if
    if(present(zeta_t)) then
      lzeta_t => zeta_t
    else
      call allocate(lzeta_t_t, decomp)
      lzeta_t => lzeta_t_t
    end if
    
    call zeta_tendency(dt, strat, beta_y, &
#ifdef WIND_FORCING
      & force_wind, &
#endif
#ifdef BACKGROUND_VSHEAR
      & u_background_v, &
#endif
      & viscosity, visc_bot, psi, rel, rel_p, lzeta_t)
    
    do layer = 1, decomp%layers
      do j = decomp%l_j_min, decomp%l_j_max
        do i = decomp%l_i_min, decomp%l_i_max
          lzeta_np1%sdata(i, j, layer) = zeta%sdata(i, j, layer) + lzeta_t%sdata(i, j, layer)
        end do
      end do
    end do
    
    if(.not. present(zeta_t)) then
      call deallocate(lzeta_t_t)
    end if
    
  end subroutine zeta_forward_euler   
  
  ! Timestep the relative potential vorticity using the centred scheme. This
  ! does *not* update the halo region.
  subroutine zeta_centred(leapfrog, T, strat, beta_y, &
#ifdef WIND_FORCING
    & force_wind, &
#endif
#ifdef BACKGROUND_VSHEAR
    & u_background_v, &
#endif
    & viscosity, visc_bot, psi, rel, rel_p, zeta)    
    integer, intent(in) :: leapfrog
    real(kind = real_8), intent(in) :: T
    type(stratification), intent(in) :: strat
    type(parallel_data), dimension(2), intent(in) :: psi
    real(kind = real_8), dimension(:, :), intent(in) :: beta_y
#ifdef WIND_FORCING
    real(kind = real_8), dimension(:, :), intent(in) :: force_wind
#endif
#ifdef BACKGROUND_VSHEAR
    real(kind = real_8), dimension(psi(2)%decomp%layers), intent(in) :: u_background_v
#endif
    real(kind = real_8), intent(in) :: viscosity
    real(kind = real_8), intent(in) :: visc_bot
    type(parallel_data), dimension(2), intent(in) :: rel
    type(parallel_data), dimension(2), intent(in) :: rel_p
    type(parallel_data), dimension(2), intent(inout) :: zeta
    
    integer :: i, j, layer, layers, m
    type(decomposition), pointer :: decomp
    
    real(kind = real_8), dimension(:, :, :), pointer :: psi_m
    real(kind = real_8), dimension(:, :, :), pointer :: rel_l, rel_m, rel_p_l, &
      & rel_p_m
    
    real(kind = real_8), dimension(:, :, :), pointer :: zeta_l
        
    assert(associated(psi(1)%decomp, psi(2)%decomp))
    assert(size(beta_y, 1) >= psi(2)%decomp%ii + 2)
    assert(size(beta_y, 2) == psi(2)%decomp%layers)
#ifdef WIND_FORCING
    assert(size(force_wind, 1) >= psi(2)%decomp%ii)
    assert(size(force_wind, 2) >= psi(2)%decomp%jj)
#endif
    assert(associated(rel(1)%decomp, psi(2)%decomp))
    assert(associated(rel(2)%decomp, psi(2)%decomp))
    assert(associated(rel_p(1)%decomp, psi(2)%decomp))
    assert(associated(rel_p(2)%decomp, psi(2)%decomp))
    assert(associated(zeta(1)%decomp, psi(2)%decomp))
    assert(associated(zeta(2)%decomp, psi(2)%decomp))
    
    decomp => psi(2)%decomp
    layers = decomp%layers
    
    m = mod(leapfrog, 2) + 1
        
    psi_m => psi(m)%data
    rel_l => rel(leapfrog)%data
    rel_m => rel(m)%data
    rel_p_l => rel_p(leapfrog)%data
    rel_p_m => rel_p(m)%data
    zeta_l => zeta(leapfrog)%data
    
    if(layers == 1) then
      do j = decomp%f_j_min, decomp%f_j_max
        do i = decomp%f_i_min, decomp%f_i_max
          zeta_l(i, j, 1) = zeta_l(i, j, 1) + &
            & T * ( &
            &   0.25D0 * ( &
#ifdef BACKGROUND_VSHEAR
            &       (beta_y(i, 1) + rel_m(i, j + 1, 1)) * (psi_m(i - 1, j + 1, 1) - psi_m(i + 1, j + 1, 1) - 2.0D0 * u_background_v(1)) &
            &     - (beta_y(i, 1) + rel_m(i, j - 1, 1)) * (psi_m(i - 1, j - 1, 1) - psi_m(i + 1, j - 1, 1) - 2.0D0 * u_background_v(1)) &
#else
            &       (beta_y(i, 1) + rel_m(i, j + 1, 1)) * (psi_m(i - 1, j + 1, 1) - psi_m(i + 1, j + 1, 1)) &
            &     - (beta_y(i, 1) + rel_m(i, j - 1, 1)) * (psi_m(i - 1, j - 1, 1) - psi_m(i + 1, j - 1, 1)) &
#endif
            &     - (beta_y(i - 1, 1) + rel_m(i - 1, j, 1)) * (psi_m(i - 1, j + 1, 1) - psi_m(i - 1, j - 1, 1)) &
            &     + (beta_y(i + 1, 1) + rel_m(i + 1, j, 1)) * (psi_m(i + 1, j + 1, 1) - psi_m(i + 1, j - 1, 1)) &
            &   ) &
#ifdef WIND_FORCING
            &   + force_wind(i - 1, j - 1) &
#endif
            &   + viscosity * ( &
            &     -4.0D0 * rel_p_l(i, j, 1) + rel_p_l(i - 1, j, 1) + rel_p_l(i + 1, j, 1) + rel_p_l(i, j - 1, 1) + rel_p_l(i, j + 1, 1) &
#ifdef CENTRED_BOTTOM_FRICTION
            &   - visc_bot * rel_p_m(i, j, 1) &
#else
            &   - visc_bot * rel_p_l(i, j, 1) &
#endif
            &   ) &
            & )
        end do
      end do
    else
      ! Upper layer
      do j = decomp%f_j_min, decomp%f_j_max
        do i = decomp%f_i_min, decomp%f_i_max
          zeta_l(i, j, 1) = zeta_l(i, j, 1) + &
            & T * ( &
            &   0.25D0 * ( &
#ifdef BACKGROUND_VSHEAR
            &       (beta_y(i, 1) + rel_m(i, j + 1, 1) + strat%s2(1) * psi_m(i, j + 1, 2)) * (psi_m(i - 1, j + 1, 1) - psi_m(i + 1, j + 1, 1) - 2.0D0 * u_background_v(1)) &
            &     - (beta_y(i, 1) + rel_m(i, j - 1, 1) + strat%s2(1) * psi_m(i, j - 1, 2)) * (psi_m(i - 1, j - 1, 1) - psi_m(i + 1, j - 1, 1) - 2.0D0 * u_background_v(1)) &
            &     + (strat%s2(1) * psi_m(i, j + 1, 1)) * (2.0D0 * u_background_v(1)) &
            &     - (strat%s2(1) * psi_m(i, j - 1, 1)) * (2.0D0 * u_background_v(1)) &
#else
            &       (beta_y(i, 1) + rel_m(i, j + 1, 1) + strat%s2(1) * psi_m(i, j + 1, 2)) * (psi_m(i - 1, j + 1, 1) - psi_m(i + 1, j + 1, 1)) &
            &     - (beta_y(i, 1) + rel_m(i, j - 1, 1) + strat%s2(1) * psi_m(i, j - 1, 2)) * (psi_m(i - 1, j - 1, 1) - psi_m(i + 1, j - 1, 1)) &
#endif
            &     - (beta_y(i - 1, 1) + rel_m(i - 1, j, 1) + strat%s2(1) * psi_m(i - 1, j, 2)) * (psi_m(i - 1, j + 1, 1) - psi_m(i - 1, j - 1, 1)) &
            &     + (beta_y(i + 1, 1) + rel_m(i + 1, j, 1) + strat%s2(1) * psi_m(i + 1, j, 2)) * (psi_m(i + 1, j + 1, 1) - psi_m(i + 1, j - 1, 1)) &
            &   ) &
#ifdef WIND_FORCING
            &   + force_wind(i - 1, j - 1) &
#endif
            &   + viscosity * ( &
            &     -4.0D0 * rel_p_l(i, j, 1) + rel_p_l(i - 1, j, 1) + rel_p_l(i + 1, j, 1) + rel_p_l(i, j - 1, 1) + rel_p_l(i, j + 1, 1) &
            &   ) &
            & )
        end do
      end do
      ! Middle layers
      do layer = 2, layers - 1
        do j = decomp%f_j_min, decomp%f_j_max
          do i = decomp%f_i_min, decomp%f_i_max
            zeta_l(i, j, layer) = zeta_l(i, j, layer) + &
              & T * ( &
              &   0.25D0 * ( &
#ifdef BACKGROUND_VSHEAR
              &       (beta_y(i, layer) + rel_m(i, j + 1, layer) + strat%s1(layer) * psi_m(i, j + 1, layer - 1) + strat%s2(layer) * psi_m(i, j + 1, layer + 1)) * (psi_m(i - 1, j + 1, layer) - psi_m(i + 1, j + 1, layer) - 2.0D0 * u_background_v(layer)) &
              &     - (beta_y(i, layer) + rel_m(i, j - 1, layer) + strat%s1(layer) * psi_m(i, j - 1, layer - 1) +  strat%s2(layer) * psi_m(i, j - 1, layer + 1)) * (psi_m(i - 1, j - 1, layer) - psi_m(i + 1, j - 1, layer) - 2.0D0 * u_background_v(layer)) &
              &     + (strat%s1(layer) * psi_m(i, j + 1, layer) + strat%s2(layer) * psi_m(i, j + 1, layer)) * (2.0D0 * u_background_v(layer)) &
              &     - (strat%s1(layer) * psi_m(i, j - 1, layer) + strat%s2(layer) * psi_m(i, j - 1, layer)) * (2.0D0 * u_background_v(layer)) &
#else
              &       (beta_y(i, layer) + rel_m(i, j + 1, layer) + strat%s1(layer) * psi_m(i, j + 1, layer - 1) + strat%s2(layer) * psi_m(i, j + 1, layer + 1)) * (psi_m(i - 1, j + 1, layer) - psi_m(i + 1, j + 1, layer)) &
              &     - (beta_y(i, layer) + rel_m(i, j - 1, layer) + strat%s1(layer) * psi_m(i, j - 1, layer - 1) +  strat%s2(layer) * psi_m(i, j - 1, layer + 1)) * (psi_m(i - 1, j - 1, layer) - psi_m(i + 1, j - 1, layer)) &
#endif
              &     - (beta_y(i - 1, layer) + rel_m(i - 1, j, layer) + strat%s1(layer) * psi_m(i - 1, j, layer - 1) + strat%s2(layer) * psi_m(i - 1, j, layer + 1)) * (psi_m(i - 1, j + 1, layer) - psi_m(i - 1, j - 1, layer)) &
              &     + (beta_y(i + 1, layer) + rel_m(i + 1, j, layer) + strat%s1(layer) * psi_m(i + 1, j, layer - 1) + strat%s2(layer) * psi_m(i + 1, j, layer + 1)) * (psi_m(i + 1, j + 1, layer) - psi_m(i + 1, j - 1, layer)) &
              &   ) &
              &   + viscosity * ( &
              &     -4.0D0 * rel_p_l(i, j, layer) + rel_p_l(i - 1, j, layer) + rel_p_l(i + 1, j, layer) + rel_p_l(i, j - 1, layer) + rel_p_l(i, j + 1, layer) &
              &   ) &
              & )
          end do
        end do
      end do
      ! Lower layer
      do j = decomp%f_j_min, decomp%f_j_max
        do i = decomp%f_i_min, decomp%f_i_max      
          zeta_l(i, j, layers) = zeta_l(i, j, layers) + &
            & T * ( &
            &   0.25D0 * ( &
#ifdef BACKGROUND_VSHEAR
            &       (beta_y(i, layers) + rel_m(i, j + 1, layers) + strat%s1(layers) * psi_m(i, j + 1, layers - 1)) * (psi_m(i - 1, j + 1, layers) - psi_m(i + 1, j + 1, layers) - 2.0D0 * u_background_v(layers)) &
            &     - (beta_y(i, layers) + rel_m(i, j - 1, layers) + strat%s1(layers) * psi_m(i, j - 1, layers - 1)) * (psi_m(i - 1, j - 1, layers) - psi_m(i + 1, j - 1, layers) - 2.0D0 * u_background_v(layers)) &
            &     + (strat%s1(layers) * psi_m(i, j + 1, layers)) * (2.0D0 * u_background_v(layers)) &
            &     - (strat%s1(layers) * psi_m(i, j - 1, layers)) * (2.0D0 * u_background_v(layers)) &
#else
            &       (beta_y(i, layers) + rel_m(i, j + 1, layers) + strat%s1(layers) * psi_m(i, j + 1, layers - 1)) * (psi_m(i - 1, j + 1, layers) - psi_m(i + 1, j + 1, layers)) &
            &     - (beta_y(i, layers) + rel_m(i, j - 1, layers) + strat%s1(layers) * psi_m(i, j - 1, layers - 1)) * (psi_m(i - 1, j - 1, layers) - psi_m(i + 1, j - 1, layers)) &
#endif
            &     - (beta_y(i - 1, layers) + rel_m(i - 1, j, layers) + strat%s1(layers) * psi_m(i - 1, j, layers - 1)) * (psi_m(i - 1, j + 1, layers) - psi_m(i - 1, j - 1, layers)) &
            &     + (beta_y(i + 1, layers) + rel_m(i + 1, j, layers) + strat%s1(layers) * psi_m(i + 1, j, layers - 1)) * (psi_m(i + 1, j + 1, layers) - psi_m(i + 1, j - 1, layers)) &
            &   ) &
            &   + viscosity * ( &
            &     -4.0D0 * rel_p_l(i, j, layers) + rel_p_l(i - 1, j, layers) + rel_p_l(i + 1, j, layers) + rel_p_l(i, j - 1, layers) + rel_p_l(i, j + 1, layers) &
            &   ) &
#ifdef CENTRED_BOTTOM_FRICTION
            &   - visc_bot * rel_p_m(i, j, layers) &
#else
            &   - visc_bot * rel_p_l(i, j, layers) &
#endif
            & )
        end do
      end do
    end if
    
  end subroutine zeta_centred
#endif

end module qg_aux
