#include "config.h"

! Diagnostics for the 3D QG model, including configuration of the
! stratification.
module qg_diagnostics

  use lapack_interfaces
  use qg_parallel
  use qg_parameters, only : real_8, real_8_chr
  use test_tools
  
  implicit none
  
  private
  
  public :: mean, psum, energy
  public :: allocate, initialise, deallocate, zero, add, print
  public :: stratification, s_matrix, z_to_modes, modes_to_z, modes_to_z_btcu
  public :: velocity_boundary_condition, zonal_velocity, &
    & meridional_velocity, buoyancy, interface_velocity, thickness, &
    & bf_pv_flux, add_scaled_buoyancy_corr
  public :: eddy_covariance
  
  ! Stratification data, including z-layer / modal conversion.
  type stratification
    ! Number of layers
    integer :: layers
    ! Layer thicknesses
    real(kind = real_8), dimension(:), pointer :: h
    ! Stratification parameters associated with layers above a given layer
    real(kind = real_8), dimension(:), pointer :: s1
    ! Stratification parameters associated with layers below a given layer
    real(kind = real_8), dimension(:), pointer :: s2
    ! Non-dimensionalisation
    real(kind = real_8) :: l_scale
    
    ! Stratification matrix
    real(kind = real_8), dimension(:, :), pointer :: s_mat
    ! Stratification parameters in modes (-ve eigenvalues of s_mat)
    real(kind = real_8), dimension(:), pointer :: ss
    ! Conversion matrix from modes to z-layers
    real(kind = real_8), dimension(:, :), pointer :: modes_to_z
    ! Conversion matrix from z-layers to modes
    real(kind = real_8), dimension(:, :), pointer :: z_to_modes
    ! Modal advection tensor
    real(kind = real_8), dimension(:, :, :), pointer :: adv
  end type stratification
  
  ! Velocity boundary condition.
  integer, parameter, public :: BC_PARTIAL_SLIP = 1, BC_NO_SLIP = 2, &
    & BC_FREE_SLIP = 3
  type velocity_boundary_condition
    integer :: bc_type            ! Boundary condition type
    real(kind = real_8) :: alpha  ! Partial slip parameter
  end type velocity_boundary_condition
  
  ! Eddy diagnostic.
  integer, parameter, public :: EC_LAYERS = 1, EC_INTERFACES = 2
  type eddy_covariance
    integer :: ec_type  ! Eddy covariance diagnostic type (layers or interfaces)
    type(decomposition), pointer :: decomp                    ! Decomposition
    real(kind = real_8), dimension(:, :, :), pointer :: data  ! Time-integrated
                                                              ! eddy covariance
  end type eddy_covariance
  
  interface mean
    module procedure mean_decomp_layer, mean_decomp_layers, mean_pdata
  end interface mean
  
  interface psum
    module procedure psum_decomp, psum_pdata
  end interface psum
  
  interface bf_pv_flux
    module procedure bf_pv_flux_periodic_data, &
      & bf_pv_flux_periodic_pdata
  end interface bf_pv_flux
  
  interface allocate
    module procedure allocate_stratification_l1, allocate_stratification_l2, &
      & allocate_stratification_l3, allocate_stratification_parameters, &
      & allocate_eddy_covariance
  end interface allocate
  
  interface initialise
    module procedure initialise_vbc
  end interface initialise
  
  interface deallocate
    module procedure deallocate_stratification, deallocate_eddy_covariance
  end interface deallocate
  
  interface zero
    module procedure zero_eddy_covariance
  end interface zero
  
  interface add
    module procedure add_eddy_covariance_pdata, &
      & add_eddy_covariance_data, add_eddy_covariance_pdata_pdata, &
      & add_eddy_covariance_data_pdata, add_eddy_covariance_data_data, &
      & add_eddy_covariance_pdata_scale, add_eddy_covariance_data_scale, &
      & add_eddy_covariance_pdata_pdata_scale, &
      & add_eddy_covariance_data_pdata_scale, &
      & add_eddy_covariance_data_data_scale
  end interface add
  
  interface print
    module procedure print_stratification, print_energy
  end interface print
  
  interface z_to_modes
    module procedure z_to_modes_rank0, z_to_modes_rank1, z_to_modes_rank2
  end interface z_to_modes
  
  interface modes_to_z
    module procedure modes_to_z_rank0, modes_to_z_rank2, modes_to_z_correction
  end interface modes_to_z
  
  interface modes_to_z_btcu
    module procedure modes_to_z_correction_btcu
  end interface modes_to_z_btcu

contains
  
  ! Allocate a (dummy) 1-layer stratification object.
  subroutine allocate_stratification_l1(strat, h, l_scale)
    type(stratification), intent(out) :: strat
    real(kind = real_8), intent(in) :: h
    real(kind = real_8), optional, intent(in) :: l_scale
    
    strat%layers = 1
    allocate(strat%h(1))    
    allocate(strat%s1(1))
    allocate(strat%s2(1))
    
    if(present(l_scale)) then
      strat%l_scale = l_scale
      strat%h(1) = h / strat%l_scale
    else
      strat%l_scale = 1.0D0
      strat%h(1) = h
    end if    
    
    call add_s_matrices(strat)
    
  end subroutine allocate_stratification_l1
  
  ! Allocate a 2-layer stratification object.
  subroutine allocate_stratification_l2(strat, rd_bcl, h1, h2, l_scale)
    type(stratification), intent(out) :: strat
    real(kind = real_8), intent(in) :: rd_bcl
    real(kind = real_8), intent(in) :: h1
    real(kind = real_8), intent(in) :: h2
    real(kind = real_8), optional, intent(in) :: l_scale
    
    strat%layers = 2
    allocate(strat%h(2))
    allocate(strat%s1(2))
    allocate(strat%s2(2))
    
    call strat_l2(rd_bcl, h1, h2, strat%s2(1), strat%s1(2))
    
    if(present(l_scale)) then
      strat%l_scale = l_scale  
      strat%h = (/h1, h2/) / l_scale
      strat%s2(1) = strat%s2(1) * (l_scale ** 2)
      strat%s1(2) = strat%s1(2) * (l_scale ** 2)
    else
      strat%l_scale = 1.0D0
      strat%h = (/h1, h2/)
    end if
    
    call add_s_matrices(strat)
    
  end subroutine allocate_stratification_l2
  
  ! Allocate a 3-layer stratification object.
  subroutine allocate_stratification_l3(strat, rd_bcl1, rd_bcl2, h1, h2, h3, l_scale)
    type(stratification), intent(out) :: strat
    real(kind = real_8), intent(in) :: rd_bcl1
    real(kind = real_8), intent(in) :: rd_bcl2
    real(kind = real_8), intent(in) :: h1
    real(kind = real_8), intent(in) :: h2
    real(kind = real_8), intent(in) :: h3
    real(kind = real_8), optional, intent(in) :: l_scale
    
    strat%layers = 3
    allocate(strat%h(3))
    allocate(strat%s1(3))
    allocate(strat%s2(3))
    
    call strat_l3(rd_bcl1, rd_bcl2, h1, h2, h3, strat%s2(1), strat%s1(2), strat%s2(2), strat%s1(3))
    
    if(present(l_scale)) then
      strat%l_scale = l_scale  
      strat%h = (/h1, h2, h3/) / l_scale
      strat%s2(1:2) = strat%s2(1:2) * (l_scale ** 2)
      strat%s1(2:3) = strat%s1(2:3) * (l_scale ** 2)
    else
      strat%l_scale = 1.0D0
      strat%h = (/h1, h2, h3/)
    end if
    
    call add_s_matrices(strat)
    
  end subroutine allocate_stratification_l3
  
  ! Allocate a stratification object based upon the supplied stratification
  ! parameters.
  subroutine allocate_stratification_parameters(strat, h, s1, s2, l_scale)
    type(stratification), intent(out) :: strat
    real(kind = real_8), dimension(:), intent(in) :: h
    real(kind = real_8), dimension(size(h) - 1), intent(in) :: s1
    real(kind = real_8), dimension(size(h) - 1), intent(in) :: s2
    real(kind = real_8), optional, intent(in) :: l_scale
        
    strat%layers = size(h)
    allocate(strat%h(size(h)))
    allocate(strat%s1(size(h)))
    allocate(strat%s2(size(h)))
    
    if(present(l_scale)) then
      strat%l_scale = l_scale
      strat%h = h / l_scale
      strat%s1(2:)        = (l_scale ** 2) * s1
      strat%s2(:size(s2)) = (l_scale ** 2) * s2
    else
      strat%l_scale = 1.0D0
      strat%h = h
      strat%s1(2:)        = s1
      strat%s2(:size(s2)) = s2
    end if
    
    call add_s_matrices(strat)
    
  end subroutine allocate_stratification_parameters
  
  ! Add stratification matrix and conversion matrices to a stratification
  ! object.
  subroutine add_s_matrices(strat)
    type(stratification), intent(inout) :: strat
    
    allocate(strat%s_mat(strat%layers, strat%layers))    
    strat%s_mat = s_matrix(strat)

    allocate(strat%ss(strat%layers))
    allocate(strat%modes_to_z(strat%layers, strat%layers))
    allocate(strat%z_to_modes(strat%layers, strat%layers))
    call s_conversion_matrices(strat%s_mat, strat%modes_to_z, strat%z_to_modes, strat%ss)
    
    allocate(strat%adv(strat%layers, strat%layers, strat%layers))
    call s_advection_tensor(strat%modes_to_z, strat%z_to_modes, strat%adv)
    
  end subroutine add_s_matrices
  
  ! Deallocate a stratification object.
  subroutine deallocate_stratification(strat)
    type(stratification), intent(inout) :: strat
    
    deallocate(strat%h)
    deallocate(strat%s1)
    deallocate(strat%s2)
    
    deallocate(strat%s_mat)
    deallocate(strat%ss)
    deallocate(strat%modes_to_z)
    deallocate(strat%z_to_modes)
    deallocate(strat%adv)
  
  end subroutine deallocate_stratification
  
  ! Compute 2-layer stratification parmeters given layer thicknesses and
  ! baroclinic radius of deformation. This inverts the eigendecomposition
  ! problem by hand.
  subroutine strat_l2(rd_bcl, h1, h2, s1, s2)
    real(kind = real_8), intent(in) :: rd_bcl
    real(kind = real_8), intent(in) :: h1
    real(kind = real_8), intent(in) :: h2
    real(kind = real_8), intent(out) :: s1
    real(kind = real_8), intent(out) :: s2
    
    real(kind = real_8) :: fsq_over_nsq

    fsq_over_nsq = (h1 * h2) / (2.0D0 * (rd_bcl ** 2))
    s1 = fsq_over_nsq * 2.0D0 / (h1 * (h1 + h2))
    s2 = fsq_over_nsq * 2.0D0 / (h2 * (h1 + h2))
    
  end subroutine strat_l2
  
  ! Compute 3-layer stratification parmeters given layer thicknesses and
  ! baroclinic radius of deformation. This inverts the eigendecomposition
  ! problem by hand.
  subroutine strat_l3(rd_bcl1, rd_bcl2, h1, h2, h3, s1, s21, s22, s3)
    real(kind = real_8), intent(in) :: rd_bcl1
    real(kind = real_8), intent(in) :: rd_bcl2
    real(kind = real_8), intent(in) :: h1
    real(kind = real_8), intent(in) :: h2
    real(kind = real_8), intent(in) :: h3
    real(kind = real_8), intent(out) :: s1
    real(kind = real_8), intent(out) :: s21
    real(kind = real_8), intent(out) :: s22
    real(kind = real_8), intent(out) :: s3  
    
    real(kind = real_8) :: a2, b1, b2, c1, c2, l
    real(kind = real_8), dimension(2) :: fsq_over_nsq, lambda
    real(kind = real_8), dimension(4) :: h
    
    lambda(1) = -1.0D0 / (rd_bcl1 ** 2)
    lambda(2) = -1.0D0 / (rd_bcl2 ** 2)
    b1 = -(lambda(1) + lambda(2))
    c1 = lambda(1) * lambda(2)
    
    h(1) = 2.0D0 / (h1 * (h1 + h2))
    h(2) = 2.0D0 / (h2 * (h1 + h2))
    h(3) = 2.0D0 / (h2 * (h2 + h3))
    h(4) = 2.0D0 / (h3 * (h2 + h3))
    
    l = h(1) * h(3) + h(1) * h(4) + h(2) * h(4)
    a2 = l * (h(3) + h(4))
    b2 = -b1 * l
    c2 = c1 * (h(1) + h(2))
    
                           ! Ambiguity in sign here. Choose -ve sign for
                           ! consistency with Pavel's code.
    fsq_over_nsq(2) = (-b2 - dsqrt(b2 ** 2 - 4.0D0 * a2 * c2)) / (2.0D0 * a2)
    fsq_over_nsq(1) = (b1 - fsq_over_nsq(2) * (h(3) + h(4))) / (h(1) + h(2))
    
    s1  = fsq_over_nsq(1) * h(1)
    s21 = fsq_over_nsq(1) * h(2)
    s22 = fsq_over_nsq(2) * h(3)
    s3 =  fsq_over_nsq(2) * h(4)

#ifdef DDEBUG    
    b1 = fsq_over_nsq(1) * (h(1) + h(2)) + fsq_over_nsq(2) * (h(3) + h(4)) 
    c1 = fsq_over_nsq(1) * fsq_over_nsq(2) * l
    assert(fe(lambda ** 2 + b1 * lambda + c1, (/0.0D0, 0.0D0/)))
#endif
    
  end subroutine strat_l3
  
  ! Form the stratification matrix.
  function s_matrix(strat)
    type(stratification), intent(in) :: strat
    
    real(kind = real_8), dimension(strat%layers, strat%layers) :: s_matrix
    
    integer :: i
    
    s_matrix = 0.0D0
    
    do i = 1, size(strat%s1)
      if(i > 1) then
        s_matrix(i, i - 1) = strat%s1(i)
        if(i < size(strat%s1)) then
          s_matrix(i, i) = -(strat%s1(i) + strat%s2(i))
          s_matrix(i, i + 1) = strat%s2(i)
        else
          s_matrix(i, i) = -strat%s1(i)
        end if
      else if(i < size(strat%s1)) then
        s_matrix(i, i) = -strat%s2(i)
        s_matrix(i, i + 1) = strat%s2(i)
      end if
    end do
    
  end function s_matrix
  
  ! Compute conversion matrices between z-layers and baroclinic modes. Also
  ! compute the modal stratification parameters. The barotropic mode is
  ! guaranteed to be associated with ss(1), and the barotropic mode eigenvector
  ! is guaranteed to have all components equal to unity.
  subroutine s_conversion_matrices(s_matrix, modes_to_z, z_to_modes, ss)
    real(kind = real_8), dimension(:, :), intent(in) :: s_matrix
    real(kind = real_8), dimension(size(s_matrix, 1), size(s_matrix, 1)), intent(out) :: modes_to_z
    real(kind = real_8), dimension(size(s_matrix, 1), size(s_matrix, 1)), intent(out) :: z_to_modes
    real(kind = real_8), dimension(size(s_matrix, 1)), intent(out) :: ss
                
    assert(size(s_matrix, 1) == size(s_matrix, 2))
    
    call eigendecompose(s_matrix, ss, modes_to_z)
    ss = -ss
    
    call sort_barotropic(ss, modes_to_z)
        
    z_to_modes = inverse(modes_to_z)
  
  contains
  
    subroutine sort_barotropic(ss, modes_to_z)
      real(kind = real_8), dimension(:), intent(inout) :: ss
      real(kind = real_8), dimension(size(ss), size(ss)), intent(inout) :: modes_to_z
      
      integer :: bt_index, i
      real(kind = real_8) :: min_s
      
      ! Find the index of the barotropic mode
      bt_index = 1
      min_s = ss(1)
      do i = 2, size(ss)
        if(ss(i) < min_s) then
          bt_index = i
          min_s = ss(i)
        end if
      end do
      write(6, "(a,"//real_8_chr//")") "Found barotropic mode with ss = ", ss(bt_index)
      if(size(ss) > 1) then
        assert(fe(modes_to_z(2:, bt_index), spread(modes_to_z(1, bt_index), 1, size(modes_to_z, 1) - 1)))
      end if
      
      ! Move the barotropic mode to the first entry
      ss(2:bt_index) = ss(:bt_index - 1)
      ss(1) = 0.0D0
      
      ! Move the barotropic mode to the first entry, and set the barotropic
      ! mode eigenvalue components to unity.
      modes_to_z(:, 2:bt_index) = modes_to_z(:, :bt_index - 1) / modes_to_z(1, 1)
      modes_to_z(:, bt_index + 1:) = modes_to_z(:, bt_index + 1:) / modes_to_z(1, 1)
      modes_to_z(:, 1) = 1.0D0
      
    end subroutine sort_barotropic
    
  end subroutine s_conversion_matrices
  
  ! Compute the modal advection tensor.
  subroutine s_advection_tensor(modes_to_z, z_to_modes, adv)
    real(kind = real_8), dimension(:, :), intent(in) :: modes_to_z
    real(kind = real_8), dimension(size(modes_to_z, 1), size(modes_to_z, 1)), intent(in) :: z_to_modes
    real(kind = real_8), dimension(size(modes_to_z, 1), size(modes_to_z, 1), size(modes_to_z, 1)), intent(out) :: adv

    integer :: i, j, k, l, layers
    
    assert(size(modes_to_z, 1) == size(modes_to_z, 2))
    
    layers = size(modes_to_z, 1)
    assert(layers >= 1)
    
    do i = 1, layers
      do k = 1, layers
        do l = 1, layers
          adv(i, k, l) = 0.0D0
          do j = 1, layers
            adv(i, k, l) = adv(i, k, l) + z_to_modes(i, j) * modes_to_z(j, k) * modes_to_z(j, l)
          end do
          if(fe(adv(i, k, l), 0.0D0, tol = 1.0D-10)) adv(i, k, l) = 0.0D0
        end do
      end do
    end do
            
  end subroutine s_advection_tensor
  
  ! Print information about a stratification object to standard output.
  subroutine print_stratification(strat)
    type(stratification), intent(in) :: strat
    
    integer :: i, j
    
    write(6, "(a,i0,a)") "Stratification, ", strat%layers, " layers:"
    write(6, "(a)") "Non-dimensional stratification parameters:"
    if(strat%layers > 1) then
      write(6, "(a,"//real_8_chr//")") "Layer 1: ", strat%s2(1)
      do i = 2, strat%layers - 1
        write(6, "(a,i0,a,"//real_8_chr//",a,"//real_8_chr//")") "Layer ", i, ": ", strat%s1(i), ", ", strat%s2(i)
      end do
      write(6, "(a,i0,a,"//real_8_chr//")") "Layer ", strat%layers, ": ", strat%s1(strat%layers)
    end if
    write(6, "(a)") "Dimensional stratification parameters:"
    if(strat%layers > 1) then
      write(6, "(a,"//real_8_chr//")") "Layer 1: ", strat%s2(1) / (strat%l_scale ** 2)
      do i = 2, strat%layers - 1
        write(6, "(a,i0,a,"//real_8_chr//",a,"//real_8_chr//")") "Layer ", i, ": ", strat%s1(i) / (strat%l_scale ** 2), ", ", strat%s2(i) / (strat%l_scale ** 2)
      end do
      write(6, "(a,i0,a,"//real_8_chr//")") "Layer ", strat%layers, ": ", strat%s1(strat%layers) / (strat%l_scale ** 2)
    end if
    write(6, "(a)") "Non-dimensional deformation radii:"
    do i = 1, strat%layers
      write(6, "("//real_8_chr//",a)", advance = "no") dsqrt(1.0D0 / strat%ss(i))
    end do
    write(6, "(a)") ""
    write(6, "(a)") "Dimensional deformation radii:"
    do i = 1, strat%layers
      write(6, "("//real_8_chr//",a)", advance = "no") dsqrt(1.0D0 / strat%ss(i)) * strat%l_scale
    end do
    write(6, "(a)") ""
    write(6, "(a)") "Layers to modes projection matrix:"
    do i = 1, strat%layers
      do j = 1, strat%layers - 1
        write(6, "("//real_8_chr//",a)", advance = "no") strat%z_to_modes(i, j), " "
      end do
      write(6, "("//real_8_chr//")") strat%z_to_modes(i, strat%layers)
    end do
    write(6, "(a)") "Modes to layers projection matrix:"
    do i = 1, strat%layers
      do j = 1, strat%layers - 1
        write(6, "("//real_8_chr//",a)", advance = "no") strat%modes_to_z(i, j), " "
      end do
      write(6, "("//real_8_chr//")") strat%modes_to_z(i, strat%layers)
    end do
    
  end subroutine print_stratification
  
  ! Convert a local real array from z-layers to modes.
  subroutine z_to_modes_rank0(strat, z, modes)
    type(stratification), intent(in) :: strat
    real(kind = real_8), dimension(strat%layers), intent(in) :: z
    real(kind = real_8), dimension(strat%layers), intent(out) :: modes
    
    integer :: layer
    
    do layer = 1, strat%layers
      modes(layer) = dot_product(strat%z_to_modes(layer, :), z)
    end do
  
  end subroutine z_to_modes_rank0
  
  ! Convert a local real array from z-layers to modes.
  subroutine z_to_modes_rank1(strat, z, modes)
    type(stratification), intent(in) :: strat
    real(kind = real_8), dimension(:, :), intent(in) :: z
    real(kind = real_8), dimension(size(z, 1), size(z, 2)), intent(out) :: modes
    
    integer :: i, layer
    
    assert(size(z, 2) == strat%layers)
    
    do i = 1, size(z, 1)
      do layer = 1, size(z, 2)
        modes(i, layer) = dot_product(strat%z_to_modes(layer, :), z(i, :))
      end do
    end do
    
  end subroutine z_to_modes_rank1
  
  ! Convert a local real array from z-layers to modes.
  subroutine z_to_modes_rank2(strat, z, modes)
    type(stratification), intent(in) :: strat
    real(kind = real_8), dimension(:, :, :), intent(in) :: z
    real(kind = real_8), dimension(size(z, 1), size(z, 2), size(z, 3)), intent(out) :: modes
    
    integer :: i, j, layer
    
    assert(size(z, 3) == strat%layers)
    
    do j = 1, size(z, 2)
      do i = 1, size(z, 1)
        ! dot_product version seems to be slightly faster
!        modes(i, j, :) = matmul(strat%z_to_modes, z(i, j, :))
        do layer = 1, size(z, 3)
          modes(i, j, layer) = dot_product(strat%z_to_modes(layer, :), z(i, j, :))
        end do
      end do
    end do
    
  end subroutine z_to_modes_rank2
  
  ! Convert a local real array from modes to z-layers.
  subroutine modes_to_z_rank0(strat, modes, z)
    type(stratification), intent(in) :: strat
    real(kind = real_8), dimension(strat%layers), intent(in) :: modes
    real(kind = real_8), dimension(strat%layers), intent(out) :: z
    
    integer :: layer
    
    do layer = 1, strat%layers
      z(layer) = modes(1) + dot_product(strat%modes_to_z(layer, 2:), modes(2:))
    end do
    
  end subroutine modes_to_z_rank0
  
  ! Convert a local real array from modes to z-layers.
  subroutine modes_to_z_rank2(strat, modes, z)
    type(stratification), intent(in) :: strat
    real(kind = real_8), dimension(:, :, :), intent(in) :: modes
    real(kind = real_8), dimension(size(modes, 1), size(modes, 2), size(modes, 3)), intent(out) :: z
    
    integer :: i, j, layer
    
    assert(size(modes, 3) == strat%layers)
    
    do j = 1, size(modes, 2)
      do i = 1, size(modes, 1)
        do layer = 1, size(modes, 3)
          z(i, j, layer) = modes(i, j, 1) + dot_product(strat%modes_to_z(layer, 2:), modes(i, j, 2:))
        end do
      end do
    end do
    
  end subroutine modes_to_z_rank2
  
  ! Convert a local real array from modes to z-layers, applying boundary
  ! correction in the conversion.
  subroutine modes_to_z_correction(strat, modes, z, psiC, corr)
    type(stratification), intent(in) :: strat
    real(kind = real_8), dimension(:, :, :), intent(in) :: modes
    real(kind = real_8), dimension(size(modes, 1), size(modes, 2), size(modes, 3)), intent(out) :: z
    real(kind = real_8), dimension(size(modes, 3)), intent(in) :: psiC
    real(kind = real_8), dimension(size(modes, 1), size(modes, 2), size(modes, 3)), intent(in) :: corr
    
    integer :: i, j, layer
    
    assert(size(modes, 3) == strat%layers)
    
    do j = 1, size(modes, 2)
      do i = 1, size(modes, 1)
        do layer = 1, size(modes, 3)
          z(i, j, layer) = modes(i, j, 1) + dot_product(strat%modes_to_z(layer, :), modes(i, j, :) - psiC * corr(i, j, :))
        end do
      end do
    end do
    
  end subroutine modes_to_z_correction
  
  ! Convert a local real array from modes to z-layers applying boundary
  ! correction in the conversion. This assumes that the barotropic correction
  ! function is equal to unity.
  subroutine modes_to_z_correction_btcu(strat, modes, z, psiC, corr)
    type(stratification), intent(in) :: strat
    real(kind = real_8), dimension(:, :, :), intent(in) :: modes
    real(kind = real_8), dimension(size(modes, 1), size(modes, 2), size(modes, 3)), intent(out) :: z
    real(kind = real_8), dimension(size(modes, 3)), intent(in) :: psiC
    real(kind = real_8), dimension(size(modes, 1), size(modes, 2), size(modes, 3)), intent(in) :: corr
    
    integer :: i, j, layer
    
    assert(size(modes, 3) == strat%layers)
    assert(fe(maxval(corr(:, :, 1)), 1.0D0) .and. fe(minval(corr(:, :, 1)), 1.0D0))
    
    do j = 1, size(modes, 2)
      do i = 1, size(modes, 1)
        do layer = 1, size(modes, 3)
          z(i, j, layer) = modes(i, j, 1) - psiC(1) + dot_product(strat%modes_to_z(layer, 2:), modes(i, j, 2:) - psiC(2:) * corr(i, j, 2:))
        end do
      end do
    end do
    
  end subroutine modes_to_z_correction_btcu
  
  ! Initialise velocity boundary conditions.
  subroutine initialise_vbc(vbc, bc_type, alpha)
    type(velocity_boundary_condition), intent(out) :: vbc
    integer, intent(in) :: bc_type
    real(kind = real_8), optional, intent(in) :: alpha
    
    vbc%bc_type = bc_type
    select case(bc_type)
      case(BC_PARTIAL_SLIP)
        if(.not. present(alpha)) then
          labort("alpha required for initialise_bc with BC_PARTIAL_SLIP")
        end if
        assert(alpha >= 0.0D0)
        vbc%alpha = alpha
      case(BC_NO_SLIP)
        if(present(alpha)) then
          write(6, "(a)") "WARNING: alpha supplied to initialise_bc with BC_NO_SLIP"
        end if
        vbc%alpha = huge(0.0D0)
      case(BC_FREE_SLIP)
        if(present(alpha)) then
          write(6, "(a)") "WARNING: alpha supplied to initialise_bc with BC_FREE_SLIP"
        end if
        vbc%alpha = 0.0D0
      case default
        write(0, "(a,i0)") "For boundary condition type: ", bc_type
        labort("Unrecognised boundary condition type")
    end select
  
  end subroutine initialise_vbc

  ! Compute the local integral of decomposed data.
  function local_integral(decomp, psi) result (integral)
    type(decomposition), intent(in) :: decomp
    real(kind = real_8), dimension(:, :), intent(in) :: psi
    
    real(kind = real_8) :: integral
    
    ! Interior
    integral = sum(psi(decomp%l_i_min:decomp%l_i_max, &
                     & decomp%l_j_min:decomp%l_j_max))
                  
    ! Top boundary
    if(decomp%p_b == -1) then
      integral = integral + 0.5D0 * sum(psi(1, decomp%l_j_min:decomp%l_j_max))
      ! Top left corner
      if(decomp%p_c == -1) integral = integral + 0.25D0 * psi(1, 1)
      ! Top right corner
      if(decomp%p_a == -1) integral = integral + 0.25D0 * psi(1, decomp%jj)
    end if
    ! Bottom boundary
    if(decomp%p_d == -1) then
      integral = integral + 0.5D0 * sum(psi(decomp%ii, decomp%l_j_min:decomp%l_j_max))
      ! Bottom left corner
      if(decomp%p_c == -1) integral = integral + 0.25D0 * psi(decomp%ii, 1)
      ! Bottom right corner
      if(decomp%p_a == -1) integral = integral + 0.25D0 * psi(decomp%ii, decomp%jj)
    end if 
    ! Left boundary
    if(decomp%p_c == -1) then
      integral = integral + 0.5D0 * sum(psi(decomp%l_i_min:decomp%l_i_max, 1))
    end if
    ! Right boundary
    if(decomp%p_a == -1) then
      integral = integral + 0.5D0 * sum(psi(decomp%l_i_min:decomp%l_i_max, decomp%jj))
    end if
    
  end function local_integral

  ! Compute the mean of decomposed data.
  function mean_decomp_layer(decomp, psi) result(psiC)
    type(decomposition), intent(in) :: decomp
    real(kind = real_8), dimension(:, :), intent(in) :: psi
    
    real(kind = real_8) :: psiC
    
    assert(size(psi, 1) == decomp%ii)
    assert(size(psi, 2) == decomp%jj)
    
    psiC = local_integral(decomp, psi)    
    call allsum(psiC)
    psiC = psiC / real(decomp%g_ii_int * decomp%g_jj_int, kind = real_8)
    
  end function mean_decomp_layer

  ! Compute the mean of decomposed data.
  function mean_decomp_layers(decomp, psi) result(psiC)
    type(decomposition), intent(in) :: decomp
    real(kind = real_8), dimension(:, :, :), intent(in) :: psi
    
    real(kind = real_8), dimension(decomp%layers) :: psiC
    
    integer :: layer
    
    assert(size(psi, 1) == decomp%ii)
    assert(size(psi, 2) == decomp%jj)
    assert(size(psi, 3) == decomp%layers)
    
    do layer = 1, decomp%layers
      psiC(layer) = local_integral(decomp, psi(:, :, layer))
    end do
    call allsum(psiC)
    psiC = psiC / real(decomp%g_ii_int * decomp%g_jj_int, kind = real_8)
    
  end function mean_decomp_layers
  
  ! Compute the mean of decomposed data.
  function mean_pdata(psi, layer) result(psiC)
    type(parallel_data), intent(in) :: psi
    integer, intent(in) :: layer
    
    real(kind = real_8) :: psiC
    
    psiC = mean(psi%decomp, psi%sdata(:, :, layer))
  
  end function mean_pdata
  
  ! Sum decomposed data.
  function psum_decomp(decomp, psi) result(psiC)
    type(decomposition), intent(in) :: decomp
    real(kind = real_8), dimension(:, :), intent(in) :: psi
    
    real(kind = real_8) :: psiC
    
    psiC = sum(psi(:decomp%ii, :decomp%jj))
    call allsum(psiC)
  
  end function psum_decomp    
  
  ! Sum decomposed data.
  function psum_pdata(psi, layer) result(psiC)
    type(parallel_data), intent(in) :: psi
    integer, intent(in) :: layer
    
    real(kind = real_8) :: psiC
    
    psiC = psum(psi%decomp, psi%sdata(:, :, layer))
  
  end function psum_pdata

  ! Compute energy diagnostics.
  subroutine energy(strat, psi, ekin, epot)
    type(parallel_data), intent(in) :: psi
    type(stratification), intent(in) :: strat
    real(kind = real_8), dimension(psi%decomp%layers), intent(out) :: ekin
    real(kind = real_8), intent(out) :: epot
    
    integer :: i, j, k, layers
    real(kind = real_8) :: cff, h_total
    real(kind = real_8), dimension(psi%decomp%layers - 1) :: fac
    type(decomposition), pointer :: decomp
    
    decomp => psi%decomp
    layers = decomp%layers
    
    h_total = sum(strat%h)
    do i = 1, layers - 1
      fac(i) = 0.5D0 * strat%s2(i) * strat%h(i) / h_total
    end do
    
    epot = 0.0D0    
    do k = 1, layers - 1
      ! Interior
      do j = decomp%l_j_min, decomp%l_j_max
        do i = decomp%l_i_min, decomp%l_i_max
          epot = epot + fac(k) * (psi%sdata(i, j, k) - psi%sdata(i, j, k + 1)) ** 2
        end do
      end do
      
      ! Left boundary
      if(decomp%p_c == -1) then
        do i = decomp%l_i_min, decomp%l_i_max
          epot = epot + 0.5D0 * fac(k) * (psi%sdata(i, 1, k) - psi%sdata(i, 1, k + 1)) ** 2
        end do
        ! Top left corner
        if(decomp%p_b == -1) epot = epot + 0.25D0 * fac(k) * (psi%sdata(1, 1, k) - psi%sdata(1, 1, k + 1)) ** 2
        ! Bottom left corner
        if(decomp%p_d == -1) epot = epot + 0.25D0 * fac(k) * (psi%sdata(decomp%ii, 1, k) - psi%sdata(decomp%ii, 1, k + 1)) ** 2
      end if
      ! Right boundary
      if(decomp%p_a == -1) then
        do i = decomp%l_i_min, decomp%l_i_max
          epot = epot + 0.5D0 * fac(k) * (psi%sdata(i, decomp%jj, k) - psi%sdata(i, decomp%jj, k + 1)) ** 2
        end do
        ! Top right corner
        if(decomp%p_b == -1) epot = epot + 0.25D0 * fac(k) * (psi%sdata(1, decomp%jj, k) - psi%sdata(1, decomp%jj, k + 1)) ** 2
        ! Bottom right corner
        if(decomp%p_d == -1) epot = epot + 0.25D0 * fac(k) * (psi%sdata(decomp%ii, decomp%jj, k) - psi%sdata(decomp%ii, decomp%jj, k + 1)) ** 2
      end if
      ! Top boundary
      if(decomp%p_b == -1) then
        do j = decomp%l_j_min, decomp%l_j_max
          epot = epot + 0.5D0 * fac(k) * (psi%sdata(1, j, k) - psi%sdata(1, j, k + 1)) ** 2
        end do
      end if
      ! Bottom boundary
      if(decomp%p_d == -1) then
        do j = decomp%l_j_min, decomp%l_j_max
          epot = epot + 0.5D0 * fac(k) * (psi%sdata(decomp%ii, j, k) - psi%sdata(decomp%ii, j, k + 1)) ** 2
        end do
      end if
    end do
    
    ekin = 0.0D0
    do k = 1, layers
      do j = decomp%f_j_min, decomp%fb_j_max
        do i = decomp%f_i_min, decomp%fb_i_max
          ekin(k) = ekin(k) + 0.25D0 * ( &
            &     (psi%data(i - 1, j, k) + psi%data(i, j, k) - psi%data(i - 1, j - 1, k) - psi%data(i, j - 1, k)) ** 2 &
            &   + (psi%data(i, j - 1, k) + psi%data(i, j, k) - psi%data(i - 1, j - 1, k) - psi%data(i - 1, j, k)) ** 2 &
            & )
        end do
      end do
    end do
    
    call allsum(epot)
    call allsum(ekin)
    cff = 1.0D0 / real(decomp%g_ii_int * decomp%g_jj_int, kind = real_8)
    epot = cff * epot
    ekin = 0.5D0 * cff * (strat%h / h_total) * ekin
    
    if(epot /= epot .or. any(ekin /= ekin)) then
      call print(ekin, epot)
      labort("NaN value encountered")
    end if
    
  end subroutine energy
  
  ! Print energy diagnostics to standard output.
  subroutine print_energy(ekin, epot)
    real(kind = real_8), dimension(:), intent(in) :: ekin
    real(kind = real_8), intent(in) :: epot
    
    integer :: layer
    
    write(6, "(a,"//real_8_chr//",a)", advance = "no") "epot = ", epot, "; "
    do layer = 1, size(ekin) - 1
      write(6, "(a,i0,a,"//real_8_chr//",a)", advance = "no") "ekin", layer, " = ", ekin(layer), "; "
    end do
    write(6, "(a,i0,a,"//real_8_chr//")") "ekin", size(ekin), " = ", ekin(size(ekin))
    
  end subroutine print_energy
  
  ! Compute modal bottom friction flux through East-West boundaries.
  subroutine bf_pv_flux_periodic_data(decomp, strat, psi, f_b, f_d, visc_bot)
    type(decomposition), intent(in) :: decomp
    type(stratification), intent(in) :: strat
    real(kind = real_8), dimension(:, :, :), intent(in) :: psi
    real(kind = real_8), dimension(decomp%layers), intent(out) :: f_b
    real(kind = real_8), dimension(decomp%layers), intent(out) :: f_d
    real(kind = real_8), optional, intent(in) :: visc_bot
       
    integer :: j, layer, layer_i, layers
       
    assert(decomp%ii >= 2)
    assert(decomp%periodic_j)
    assert(size(psi, 1) == decomp%ii)
    assert(size(psi, 2) == decomp%jj)
    assert(size(psi, 3) == decomp%layers)
    
    layers = decomp%layers
    
    f_b = 0.0D0
    f_d = 0.0D0
    if(decomp%p_b == -1) then
      do layer = 1, layers
        do j = 1, decomp%jj
          do layer_i = 1, layers
            f_b(layer) = f_b(layer) + strat%z_to_modes(layer, layers) * strat%modes_to_z(layers, layer_i) * (psi(2, j, layer_i) - psi(1, j, layer_i))
          end do
        end do
      end do
    end if
    if(decomp%p_d == -1) then
      do layer = 1, layers
        do j = 1, decomp%jj
          do layer_i = 1, layers
            f_d(layer) = f_d(layer) + strat%z_to_modes(layer, layers) * strat%modes_to_z(layers, layer_i) * (psi(decomp%ii - 1, j, layer_i) - psi(decomp%ii, j, layer_i))
          end do
        end do
      end do
    end if
    call allsum(f_b)
    call allsum(f_d)
    if(present(visc_bot)) then
      f_b = f_b * visc_bot
      f_d = f_d * visc_bot
    end if
    
  end subroutine bf_pv_flux_periodic_data
  
  ! Compute the bottom friction flux through East-West boundaries.
  subroutine bf_pv_flux_periodic_pdata(strat, psi, f_b, f_d, visc_bot)
    type(stratification), intent(in) :: strat
    type(parallel_data), intent(in) :: psi
    real(kind = real_8), dimension(psi%decomp%layers), intent(out) :: f_b
    real(kind = real_8), dimension(psi%decomp%layers), intent(out) :: f_d
    real(kind = real_8), optional, intent(in) :: visc_bot
    
    call bf_pv_flux(psi%decomp, strat, psi%sdata, f_b, f_d, visc_bot = visc_bot)
    
  end subroutine bf_pv_flux_periodic_pdata
  
  ! Compute the zonal velocity.
  subroutine zonal_velocity(vbc, psi, rel, u)
    type(velocity_boundary_condition), intent(in) :: vbc
    type(parallel_data), intent(in) :: psi
    type(parallel_data), intent(in) :: rel
    real(kind = real_8), dimension(:, :, :), intent(out) :: u
    
    integer :: i, j, layer
    type(decomposition), pointer :: decomp
    
    assert(associated(rel%decomp, psi%decomp))
    assert(size(u, 1) == psi%decomp%ii)
    assert(size(u, 2) == psi%decomp%jj)
    assert(size(u, 3) == psi%decomp%layers)
  
    decomp => psi%decomp
  
    ! Interior
    do layer = 1, decomp%layers
      do j = decomp%f_j_min, decomp%f_j_max
        do i = decomp%f_i_min, decomp%f_i_max
          u(i - 1, j - 1, layer) = 0.5D0 * (psi%data(i + 1, j, layer) - psi%data(i - 1, j, layer))
        end do
      end do
    end do
    
    select case(vbc%bc_type)
      case(BC_PARTIAL_SLIP)
        ! Upper boundary condition
        if(decomp%p_b == -1) then
          do layer = 1, decomp%layers
            do j = decomp%f_j_min, decomp%f_j_max
              u(1, j - 1, layer) = rel%sdata(1, j - 1, layer) / vbc%alpha
            end do
          end do
        end if
        ! Lower boundary condition
        if(decomp%p_d == -1) then
          do layer = 1, decomp%layers
            do j = decomp%f_j_min, decomp%f_j_max
              u(decomp%ii, j - 1, layer) = rel%sdata(decomp%ii, j - 1, layer) / vbc%alpha
            end do
          end do
        end if
      case(BC_NO_SLIP)        
        ! Upper boundary condition
        if(decomp%p_b == -1) then
          do layer = 1, decomp%layers
            do j = decomp%f_j_min, decomp%f_j_max
              u(1, j - 1, layer) = 0.0D0
            end do
          end do
        end if
        ! Lower boundary condition
        if(decomp%p_d == -1) then
          do layer = 1, decomp%layers
            do j = decomp%f_j_min, decomp%f_j_max
              u(decomp%ii, j - 1, layer) = 0.0D0
            end do
          end do
        end if      
      case(BC_FREE_SLIP)        
        ! Upper boundary condition
        if(decomp%p_b == -1) then
          do layer = 1, decomp%layers
            do j = decomp%f_j_min, decomp%f_j_max
              u(1, j - 1, layer) = nan()
            end do
          end do
        end if
        ! Lower boundary condition
        if(decomp%p_d == -1) then
          do layer = 1, decomp%layers
            do j = decomp%f_j_min, decomp%f_j_max
              u(decomp%ii, j - 1, layer) = nan()
            end do
          end do
        end if      
      case default
        write(0, "(a,i0)") "For boundary condition type: ", vbc%bc_type
        labort("Unrecognised boundary condition type")
    end select
    
  end subroutine zonal_velocity
  
  ! Compute the meridional velocity.
  subroutine meridional_velocity(vbc, psi, rel, v)
    type(velocity_boundary_condition), intent(in) :: vbc
    type(parallel_data), intent(in) :: psi
    type(parallel_data), intent(in) :: rel
    real(kind = real_8), dimension(:, :, :), intent(out) :: v
    
    integer :: i, j, layer
    type(decomposition), pointer :: decomp
    
    assert(associated(rel%decomp, psi%decomp))
    assert(size(v, 1) == psi%decomp%ii)
    assert(size(v, 2) == psi%decomp%jj)
    assert(size(v, 3) == psi%decomp%layers)
  
    decomp => psi%decomp
  
    ! Interior
    do layer = 1, decomp%layers
      do j = decomp%f_j_min, decomp%f_j_max
        do i = decomp%f_i_min, decomp%f_i_max
            v(i - 1, j - 1, layer) = 0.5D0 * (psi%data(i, j + 1, layer) - psi%data(i, j - 1, layer))
          end do
        end do
    end do
    select case(vbc%bc_type)
      case(BC_PARTIAL_SLIP)
        ! Left boundary
        if(decomp%p_c == -1) then
          do layer = 1, decomp%layers
            do i = decomp%f_i_min, decomp%f_i_max
              v(i - 1, 1, layer) = rel%sdata(i - 1, 1, layer) / vbc%alpha
            end do
          end do
        end if
        ! Right boundary
        if(decomp%p_a == -1) then
          do layer = 1, decomp%layers
           do i = decomp%f_i_min, decomp%f_i_max
              v(i - 1, decomp%jj, layer) = rel%sdata(i - 1, decomp%jj, layer) / vbc%alpha
            end do
          end do
        end if
      case(BC_NO_SLIP)
        ! Left boundary
        if(decomp%p_c == -1) then
          do layer = 1, decomp%layers
            do i = decomp%f_i_min, decomp%f_i_max
              v(i - 1, 1, layer) = 0.0D0
            end do
          end do
        end if
        ! Right boundary
        if(decomp%p_a == -1) then
          do layer = 1, decomp%layers
            do i = decomp%f_i_min, decomp%f_i_max
              v(i - 1, decomp%jj, layer) = 0.0D0
            end do
          end do
        end if  
      case(BC_FREE_SLIP)
        ! Left boundary
        if(decomp%p_c == -1) then
          do layer = 1, decomp%layers
            do i = decomp%f_i_min, decomp%f_i_max
              v(i - 1, 1, layer) = nan()
            end do
          end do
        end if
        ! Right boundary
        if(decomp%p_a == -1) then
          do layer = 1, decomp%layers
            do i = decomp%f_i_min, decomp%f_i_max
              v(i - 1, decomp%jj, layer) = nan()
            end do
          end do
        end if           
      case default
        write(0, "(a,i0)") "For boundary condition type: ", vbc%bc_type
        labort("Unrecognised boundary condition type")
    end select
  
  end subroutine meridional_velocity
  
  ! Compute the buoyancy, interpolated onto layers.
  subroutine buoyancy(strat, psi, b)
    type(stratification), intent(in) :: strat
    type(parallel_data), intent(in) :: psi
    real(kind = real_8), dimension(:, :, :), intent(out) :: b
    
    integer :: i, j, layer, layers
    type(decomposition), pointer :: decomp
    
    assert(size(b, 1) == psi%decomp%ii)
    assert(size(b, 2) == psi%decomp%jj)
    assert(size(b, 3) == psi%decomp%layers)
    
    decomp => psi%decomp
    layers = decomp%layers
    
    select case(layers)
      case(3)
        do j = 1, decomp%jj
          do i = 1, decomp%ii
            b(i, j, 1) = (psi%sdata(i, j, 1) - psi%sdata(i, j, 2)) / (strat%h(1) + strat%h(2))
            b(i, j, 3) = (psi%sdata(i, j, 2) - psi%sdata(i, j, 3)) / (strat%h(2) + strat%h(3))
            b(i, j, 2) = b(i, j, 1) + b(i, j, 3)
          end do
        end do
      case(2)
        do j = 1, decomp%jj
          do i = 1, decomp%ii
            b(i, j, 1) = (psi%sdata(i, j, 1) - psi%sdata(i, j, 2)) / (strat%h(1) + strat%h(2))
            b(i, j, 2) = b(i, j, 1)
          end do
        end do
      case(1)
        b = 0.0D0
      case default
        do j = 1, decomp%jj
          do i = 1, decomp%ii
            ! Upper layer
            b(i, j, 1) = (psi%sdata(i, j, 1) - psi%sdata(i, j, 2)) / (strat%h(1) + strat%h(2))
            ! Middle layers
            do layer = 2, layers - 1
              b(i, j, layer) = &
                &  ((psi%sdata(i, j, layer - 1) - psi%sdata(i, j, layer)) / (strat%h(layer - 1) + strat%h(layer))) + &
                &  ((psi%sdata(i, j, layer) - psi%sdata(i, j, layer + 1)) / (strat%h(layer) + strat%h(layer + 1)))
            end do
            ! Lower layer
            b(i, j, layers) = (psi%sdata(i, j, layers - 1) - psi%sdata(i, j, layers)) / (strat%h(layers - 1) + strat%h(layers))
          end do
        end do
    end select
    
  end subroutine buoyancy
  
  ! Compute the velocity, interpolated onto interfaces.
  subroutine interface_velocity(decomp, u, v, u_i, v_i)
    type(decomposition), intent(in) :: decomp
    real(kind = real_8), dimension(:, :, :), intent(in) :: u
    real(kind = real_8), dimension(:, :, :), intent(in) :: v
    real(kind = real_8), dimension(:, :, :), intent(out) :: u_i
    real(kind = real_8), dimension(:, :, :), intent(out) :: v_i
    
    integer :: i, j, layer
    
    assert(size(u, 1) == decomp%ii)
    assert(size(u, 2) == decomp%jj)
    assert(size(u, 3) == decomp%layers)
    assert(size(v, 1) == decomp%ii)
    assert(size(v, 2) == decomp%jj)
    assert(size(v, 3) == decomp%layers)
    assert(size(u_i, 1) == decomp%ii)
    assert(size(u_i, 2) == decomp%jj)
    assert(size(u_i, 3) == decomp%layers - 1)
    assert(size(v_i, 1) == decomp%ii)
    assert(size(v_i, 2) == decomp%jj)
    assert(size(v_i, 3) == decomp%layers - 1)
    
    do layer = 1, decomp%layers - 1
      do j = 1, decomp%jj
        do i = 1, decomp%ii
          u_i(i, j, layer) = 0.5D0 * (u(i, j, layer) + u(i, j, layer + 1))
          v_i(i, j, layer) = 0.5D0 * (v(i, j, layer) + v(i, j, layer + 1))
        end do
      end do
    end do
    
  end subroutine interface_velocity

  ! Compute the thickness part of the potential vorticity.
  subroutine thickness(psi, strat, t)
    type(parallel_data), intent(in) :: psi
    type(stratification), intent(in) :: strat
    real(kind = real_8), dimension(:, :, :), intent(out) :: t

    integer :: i, j, layer, layers
    type(decomposition), pointer :: decomp

    assert(size(t, 1) == psi%decomp%ii)
    assert(size(t, 2) == psi%decomp%jj)
    assert(size(t, 3) == psi%decomp%layers)

    decomp => psi%decomp
    layers = decomp%layers

    do j = 1, decomp%jj
      do i = 1, decomp%ii
        t(i, j, 1) = - strat%s2(1) * (psi%sdata(i, j, 1) - psi%sdata(i, j, 2))
        t(i, j, layers) = - strat%s1(layers) * (psi%sdata(i, j, layers) - psi%sdata(i, j, layers - 1))
      end do
    end do
    do layer = 2, layers - 1
      do j = 1, decomp%jj
        do i = 1, decomp%ii
          t(i, j, layer) = - strat%s1(layer) * (psi%sdata(i, j, layer) - psi%sdata(i, j, layer - 1)) &
                         & - strat%s2(layer) * (psi%sdata(i, j, layer) - psi%sdata(i, j, layer + 1))
        end do
      end do
    end do

  end subroutine thickness
  
  ! Allocate an eddy diagnostic type.
  subroutine allocate_eddy_covariance(corr, decomp, ec_type)
    type(eddy_covariance), intent(out) :: corr
    type(decomposition), target, intent(in) :: decomp
    integer, optional, intent(in) :: ec_type
    
    integer :: lec_type
    
    if(present(ec_type)) then
      lec_type = ec_type
    else
      lec_type = EC_LAYERS
    end if
    
    corr%ec_type = lec_type
    corr%decomp => decomp
    select case(lec_type)
      case(EC_LAYERS)
        allocate(corr%data(decomp%ii, decomp%jj, decomp%layers))
      case(EC_INTERFACES)
        allocate(corr%data(decomp%ii, decomp%jj, decomp%layers - 1))
      case default
        write(6, "(a,i0)") "For eddy covariance type: ", lec_type
        labort("Unrecognised eddy covariance type")
    end select
  
  end subroutine allocate_eddy_covariance
  
  ! Deallocate an eddy diagnostic type.
  subroutine deallocate_eddy_covariance(corr)
    type(eddy_covariance), intent(inout) :: corr
    
    deallocate(corr%data)
  
  end subroutine deallocate_eddy_covariance
  
  ! Zero an eddy diagnostic type.
  subroutine zero_eddy_covariance(corr)
    type(eddy_covariance), intent(inout) :: corr
    
    corr%data = 0.0D0
  
  end subroutine zero_eddy_covariance
  
  ! Add to an eddy diagnostic type.
  subroutine add_eddy_covariance_pdata(corr, data1)
    type(eddy_covariance), intent(inout) :: corr
    type(parallel_data), intent(in) :: data1
    
    assert(corr%ec_type == EC_LAYERS)
    assert(associated(data1%decomp, corr%decomp))
    
    call add(corr, data1%sdata)
  
  end subroutine add_eddy_covariance_pdata
  
  ! Add to an eddy diagnostic type.
  subroutine add_eddy_covariance_data(corr, data1)
    type(eddy_covariance), intent(inout) :: corr
    real(kind = real_8), dimension(:, :, :), intent(in) :: data1
    
    integer :: i, j, k
    
    assert(size(data1, 1) == corr%decomp%ii)
    assert(size(data1, 2) == corr%decomp%jj)
    assert(size(data1, 3) == size(corr%data, 3))
    
    do k = 1, size(corr%data, 3)
      do j = 1, corr%decomp%jj
        do i = 1, corr%decomp%ii
          corr%data(i, j, k) = corr%data(i, j, k) + data1(i, j, k) ** 2
        end do
      end do
    end do
  
  end subroutine add_eddy_covariance_data
  
  ! Add to an eddy diagnostic type.
  subroutine add_eddy_covariance_pdata_pdata(corr, data1, data2)
    type(eddy_covariance), intent(inout) :: corr
    type(parallel_data), intent(in) :: data1
    type(parallel_data), intent(in) :: data2
    
    assert(corr%ec_type == EC_LAYERS)
    assert(associated(data1%decomp, corr%decomp))
    assert(associated(data2%decomp, corr%decomp))
    
    call add(corr, data1%sdata, data2%sdata)
  
  end subroutine add_eddy_covariance_pdata_pdata
  
  ! Add to an eddy diagnostic type.
  subroutine add_eddy_covariance_data_pdata(corr, data1, data2)
    type(eddy_covariance), intent(inout) :: corr
    real(kind = real_8), dimension(:, :, :), intent(in) :: data1
    type(parallel_data), intent(in) :: data2
    
    assert(corr%ec_type == EC_LAYERS)
    assert(associated(data2%decomp, corr%decomp))
    
    call add(corr, data1, data2%sdata)
  
  end subroutine add_eddy_covariance_data_pdata
  
  ! Add to an eddy diagnostic type.
  subroutine add_eddy_covariance_data_data(corr, data1, data2)
    type(eddy_covariance), intent(inout) :: corr
    real(kind = real_8), dimension(:, :, :), intent(in) :: data1
    real(kind = real_8), dimension(:, :, :), intent(in) :: data2
    
    integer :: i, j, k
    
    assert(size(data1, 1) == corr%decomp%ii)
    assert(size(data1, 2) == corr%decomp%jj)
    assert(size(data1, 3) == size(corr%data, 3))
    assert(size(data2, 1) == corr%decomp%ii)
    assert(size(data2, 2) == corr%decomp%jj)
    assert(size(data2, 3) == size(corr%data, 3))
    
    do k = 1, size(corr%data, 3)
      do j = 1, corr%decomp%jj
        do i = 1, corr%decomp%ii
          corr%data(i, j, k) = corr%data(i, j, k) + data1(i, j, k) * data2(i, j, k)
        end do
      end do
    end do
  
  end subroutine add_eddy_covariance_data_data
  
  ! Add to an eddy diagnostic type.
  subroutine add_eddy_covariance_pdata_scale(corr, data1, dt)
    type(eddy_covariance), intent(inout) :: corr
    type(parallel_data), intent(in) :: data1
    real(kind = real_8), intent(in) :: dt
    
    assert(corr%ec_type == EC_LAYERS)
    assert(associated(data1%decomp, corr%decomp))
    
    call add(corr, data1%sdata, dt = dt)
  
  end subroutine add_eddy_covariance_pdata_scale
  
  ! Add to an eddy diagnostic type.
  subroutine add_eddy_covariance_data_scale(corr, data1, dt)
    type(eddy_covariance), intent(inout) :: corr
    real(kind = real_8), dimension(:, :, :), intent(in) :: data1
    real(kind = real_8), intent(in) :: dt
    
    integer :: i, j, k
    
    assert(size(data1, 1) == corr%decomp%ii)
    assert(size(data1, 2) == corr%decomp%jj)
    assert(size(data1, 3) == size(corr%data, 3))
    
    do k = 1, size(corr%data, 3)
      do j = 1, corr%decomp%jj
        do i = 1, corr%decomp%ii
          corr%data(i, j, k) = corr%data(i, j, k) + dt * data1(i, j, k) ** 2
        end do
      end do
    end do
  
  end subroutine add_eddy_covariance_data_scale
  
  ! Add to an eddy diagnostic type.
  subroutine add_eddy_covariance_pdata_pdata_scale(corr, data1, data2, dt)
    type(eddy_covariance), intent(inout) :: corr
    type(parallel_data), intent(in) :: data1
    type(parallel_data), intent(in) :: data2
    real(kind = real_8), intent(in) :: dt
    
    assert(corr%ec_type == EC_LAYERS)
    assert(associated(data1%decomp, corr%decomp))
    assert(associated(data2%decomp, corr%decomp))
    
    call add(corr, data1%sdata, data2%sdata, dt = dt)
  
  end subroutine add_eddy_covariance_pdata_pdata_scale
  
  ! Add to an eddy diagnostic type.
  subroutine add_eddy_covariance_data_pdata_scale(corr, data1, data2, dt)
    type(eddy_covariance), intent(inout) :: corr
    real(kind = real_8), dimension(:, :, :), intent(in) :: data1
    type(parallel_data), intent(in) :: data2
    real(kind = real_8), intent(in) :: dt
    
    assert(corr%ec_type == EC_LAYERS)
    assert(associated(data2%decomp, corr%decomp))
    
    call add(corr, data1, data2%sdata, dt = dt)
  
  end subroutine add_eddy_covariance_data_pdata_scale
  
  ! Add to an eddy diagnostic type.
  subroutine add_eddy_covariance_data_data_scale(corr, data1, data2, dt)
    type(eddy_covariance), intent(inout) :: corr
    real(kind = real_8), dimension(:, :, :), intent(in) :: data1
    real(kind = real_8), dimension(:, :, :), intent(in) :: data2
    real(kind = real_8), intent(in) :: dt
    
    integer :: i, j, k
    
    assert(size(data1, 1) == corr%decomp%ii)
    assert(size(data1, 2) == corr%decomp%jj)
    assert(size(data1, 3) == size(corr%data, 3))
    assert(size(data2, 1) == corr%decomp%ii)
    assert(size(data2, 2) == corr%decomp%jj)
    assert(size(data2, 3) == size(corr%data, 3))
    
    do k = 1, size(corr%data, 3)
      do j = 1, corr%decomp%jj
        do i = 1, corr%decomp%ii
          corr%data(i, j, k) = corr%data(i, j, k) + dt * data1(i, j, k) * data2(i, j, k)
        end do
      end do
    end do
  
  end subroutine add_eddy_covariance_data_data_scale
  
  ! Add scaled buoyancy covariances.
  subroutine add_scaled_buoyancy_corr(sc_ub_av, sc_vb_av, sc_bb_av, sc_bb_i_av, strat, u, v, psi, dt)
    type(eddy_covariance), intent(inout) :: sc_ub_av
    type(eddy_covariance), intent(inout) :: sc_vb_av
    type(eddy_covariance), intent(inout) :: sc_bb_av
    type(eddy_covariance), intent(inout) :: sc_bb_i_av
    type(stratification), intent(in) :: strat
    real(kind = real_8), dimension(:, :, :), intent(in) :: u
    real(kind = real_8), dimension(:, :, :), intent(in) :: v
    type(parallel_data), intent(in) :: psi
    real(kind = real_8), intent(in) :: dt
    
    integer :: i, j, layer, layers
    real(kind = real_8) :: sc_b
    type(decomposition), pointer :: decomp
    
    assert(sc_ub_av%ec_type == EC_INTERFACES)
    assert(sc_vb_av%ec_type == EC_INTERFACES)
    assert(sc_bb_av%ec_type == EC_LAYERS)
    assert(sc_bb_i_av%ec_type == EC_INTERFACES)
    assert(associated(sc_ub_av%decomp, psi%decomp))
    assert(associated(sc_vb_av%decomp, psi%decomp))
    assert(associated(sc_bb_av%decomp, psi%decomp))
    assert(associated(sc_bb_i_av%decomp, psi%decomp))
    assert(size(u, 1) == psi%decomp%ii)
    assert(size(u, 2) == psi%decomp%jj)
    assert(size(u, 3) == psi%decomp%layers)
    assert(size(v, 1) == psi%decomp%ii)
    assert(size(v, 2) == psi%decomp%jj)
    assert(size(v, 3) == psi%decomp%layers)
    
    decomp => psi%decomp
    layers = decomp%layers
    
    if(layers == 1) return
    
    do layer = 1, layers - 1
      do j = 1, decomp%jj
        do i = 1, decomp%ii
          sc_b = strat%s2(layer) * strat%h(layer) * (psi%sdata(i, j, layer) - psi%sdata(i, j, layer + 1))
          sc_ub_av%data(i, j, layer) = sc_ub_av%data(i, j, layer) + dt * &
            & ( &
            &   0.5D0 * (u(i, j, layer + 1) + u(i, j, layer)) * sc_b &
            & )
          sc_vb_av%data(i, j, layer) = sc_vb_av%data(i, j, layer) + dt * &
            & ( &
            &   0.5D0 * (v(i, j, layer + 1) + v(i, j, layer)) * sc_b &
            & )
            
          sc_bb_i_av%data(i, j, layer) = sc_bb_i_av%data(i, j, layer) + dt * &
            & ( &
            &   strat%s2(layer) * strat%h(layer) * ((psi%sdata(i, j, layer) - psi%sdata(i, j, layer + 1)) ** 2) / (strat%h(layer) + strat%h(layer + 1)) &
            & )
        end do
      end do
    end do  
    
    do j = 1, decomp%jj
      do i = 1, decomp%ii
        sc_bb_av%data(i, j, 1) = sc_bb_av%data(i, j, 1) + dt * &
          & ( &
          &   0.25D0 * (strat%s2(1) * ((psi%sdata(i, j, 2) - psi%sdata(i, j, 1)) ** 2)) &
          & )
      end do
    end do
      
    do layer = 2, layers - 1
      do j = 1, decomp%jj
        do i = 1, decomp%ii
          sc_bb_av%data(i, j, layer) = sc_bb_av%data(i, j, layer) + dt * &
            & ( &
            &   0.25D0 * (strat%s1(layer) * ((psi%sdata(i, j, layer - 1) - psi%sdata(i, j, layer)) ** 2) &
            &           + strat%s2(layer) * ((psi%sdata(i, j, layer + 1) - psi%sdata(i, j, layer)) ** 2)) &
            & )
        end do
      end do
    end do
    
    do j = 1, decomp%jj
      do i = 1, decomp%ii
        sc_bb_av%data(i, j, layers) = sc_bb_av%data(i, j, layers) + dt * &
          & ( &
          &   0.25D0 * (strat%s1(layers) * ((psi%sdata(i, j, layers - 1) - psi%sdata(i, j, layers)) ** 2)) &
          & )
      end do
    end do
    
  end subroutine add_scaled_buoyancy_corr
              
end module qg_diagnostics
