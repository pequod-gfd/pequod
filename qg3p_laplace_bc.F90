#include "config.h"

! Routines for the application of Dirichlet boundary conditions to a discrete
! Laplace / Helmholtz operator.
module qg_laplace_bc
  
  use qg_parameters, only : real_8
  use qg_parallel
  
  implicit none
  
  private
  
  public :: laplace_bc_rhs, add_laplace_bc_rhs, apply_dirichlet_bc
  
  interface add_laplace_bc_rhs
    module procedure add_laplace_bc_rhs_layer, &
      & add_laplace_bc_rhs_layers_single, add_laplace_bc_rhs_layers_multiple, &
      & add_laplace_bc_rhs_periodic_layer, &
      & add_laplace_bc_rhs_periodic_layers_single, &
      & add_laplace_bc_rhs_periodic_layers_multiple
  end interface add_laplace_bc_rhs
  
  interface apply_dirichlet_bc
    module procedure apply_dirichlet_bc_layer, &
      & apply_dirichlet_bc_layers_single_data, &
      & apply_dirichlet_bc_layers_single_pdata, &
      & apply_dirichlet_bc_layers_multiple_data, &
      & apply_dirichlet_bc_layers_multiple_pdata, &
      & apply_dirichlet_bc_periodic_layer_data, &
      & apply_dirichlet_bc_periodic_layer_pdata, &
      & apply_dirichlet_bc_periodic_layers_single_data, &
      & apply_dirichlet_bc_periodic_layers_single_pdata, &
      & apply_dirichlet_bc_periodic_layers_multiple_data, &
      & apply_dirichlet_bc_periodic_layers_multiple_pdata
  end interface apply_dirichlet_bc
  
contains

  ! Return the RHS associated with the Dirichlet boundary condition for a
  ! Poisson / Helmholtz system.
  function laplace_bc_rhs(decomp, bc) result(rhs)
    type(decomposition), intent(in) :: decomp
    real(kind = real_8), intent(in) :: bc
    
    real(kind = real_8), dimension(decomp%ii, decomp%jj) :: rhs
    
    rhs = 0.0D0
    call add_laplace_bc_rhs(decomp, bc, rhs)
    
  end function laplace_bc_rhs
  
  ! Add the RHS associated with the Dirichlet boundary condition for a Poisson
  ! / Helmholtz system.
  subroutine add_laplace_bc_rhs_layer(decomp, bc, rhs)
    type(decomposition), intent(in) :: decomp
    real(kind = real_8), intent(in) :: bc
    real(kind = real_8), dimension(:, :), intent(inout) :: rhs
    
    integer :: i, j
    
    assert(size(rhs, 1) == decomp%ii)
    assert(size(rhs, 2) == decomp%jj)
    
    if(decomp%p_a == -1) then
      do i = decomp%l_i_min, decomp%l_i_max
        rhs(i, decomp%l_j_max) = rhs(i, decomp%l_j_max) - bc
      end do
    end if
    if(decomp%p_c == -1) then
      do i = decomp%l_i_min, decomp%l_i_max
        rhs(i, decomp%l_j_min) = rhs(i, decomp%l_j_min) - bc
      end do
    end if
    if(decomp%p_b == -1) then
      do j = decomp%l_j_min, decomp%l_j_max
        rhs(decomp%l_i_min, j) = rhs(decomp%l_i_min, j) - bc
      end do
    end if
    if(decomp%p_d == -1) then
      do j = decomp%l_j_min, decomp%l_j_max
        rhs(decomp%l_i_max, j) = rhs(decomp%l_i_max, j) - bc
      end do
    end if
    
  end subroutine add_laplace_bc_rhs_layer
  
  ! Add the RHS associated with the Dirichlet boundary condition for a Poisson
  ! / Helmholtz system.
  subroutine add_laplace_bc_rhs_layers_single(decomp, bc, rhs)
    type(decomposition), intent(in) :: decomp
    real(kind = real_8), intent(in) :: bc
    real(kind = real_8), dimension(:, :, :), intent(inout) :: rhs
    
    integer :: layer
    
    assert(size(rhs, 1) == decomp%ii)
    assert(size(rhs, 2) == decomp%jj)
    assert(size(rhs, 3) == decomp%layers)
    
    do layer = 1, decomp%layers
      call add_laplace_bc_rhs(decomp, bc, rhs(:, :, layer))
    end do
    
  end subroutine add_laplace_bc_rhs_layers_single
  
  ! Add the RHS associated with the Dirichlet boundary condition for a Poisson
  ! / Helmholtz system.
  subroutine add_laplace_bc_rhs_layers_multiple(decomp, bc, rhs)
    type(decomposition), intent(in) :: decomp
    real(kind = real_8), dimension(decomp%layers), intent(in) :: bc
    real(kind = real_8), dimension(:, :, :), intent(inout) :: rhs
    
    integer :: layer
    
    assert(size(rhs, 1) == decomp%ii)
    assert(size(rhs, 2) == decomp%jj)
    assert(size(rhs, 3) == decomp%layers)
    
    do layer = 1, decomp%layers
      call add_laplace_bc_rhs(decomp, bc(layer), rhs(:, :, layer))
    end do
    
  end subroutine add_laplace_bc_rhs_layers_multiple
  
  ! Add the RHS associated with the Dirichlet boundary condition for a Poisson
  ! / Helmholtz system.
  subroutine add_laplace_bc_rhs_periodic_layer(decomp, bc_b, bc_d, rhs)
    type(decomposition), intent(in) :: decomp
    real(kind = real_8), intent(in) :: bc_b
    real(kind = real_8), intent(in) :: bc_d
    real(kind = real_8), dimension(:, :), intent(inout) :: rhs
    
    integer :: j
    
    assert(size(rhs, 1) == decomp%ii)
    assert(size(rhs, 2) == decomp%jj)
    
    if(decomp%p_b == -1) then
      do j = decomp%l_j_min, decomp%l_j_max
        rhs(decomp%l_i_min, j) = rhs(decomp%l_i_min, j) - bc_b
      end do
    end if
    if(decomp%p_d == -1) then
      do j = decomp%l_j_min, decomp%l_j_max
        rhs(decomp%l_i_max, j) = rhs(decomp%l_i_max, j) - bc_d
      end do
    end if
    
  end subroutine add_laplace_bc_rhs_periodic_layer
  
  ! Add the RHS associated with the Dirichlet boundary condition for a Poisson
  ! / Helmholtz system.
  subroutine add_laplace_bc_rhs_periodic_layers_single(decomp, bc_b, bc_d, rhs)
    type(decomposition), intent(in) :: decomp
    real(kind = real_8), intent(in) :: bc_b
    real(kind = real_8), intent(in) :: bc_d
    real(kind = real_8), dimension(:, :, :), intent(inout) :: rhs
    
    integer :: layer
    
    assert(size(rhs, 1) == decomp%ii)
    assert(size(rhs, 2) == decomp%jj)
    assert(size(rhs, 3) == decomp%layers)
    
    do layer = 1, decomp%layers
      call add_laplace_bc_rhs(decomp, bc_b, bc_d, rhs(:, :, layer))
    end do
    
  end subroutine add_laplace_bc_rhs_periodic_layers_single
  
  ! Add the RHS associated with the Dirichlet boundary condition for a Poisson
  ! / Helmholtz system.
  subroutine add_laplace_bc_rhs_periodic_layers_multiple(decomp, bc_b, bc_d, rhs)
    type(decomposition), intent(in) :: decomp
    real(kind = real_8), dimension(decomp%layers), intent(in) :: bc_b
    real(kind = real_8), dimension(decomp%layers), intent(in) :: bc_d
    real(kind = real_8), dimension(:, :, :), intent(inout) :: rhs
    
    integer :: layer
    
    assert(size(rhs, 1) == decomp%ii)
    assert(size(rhs, 2) == decomp%jj)
    assert(size(rhs, 2) == decomp%layers)
    
    do layer = 1, decomp%layers
      call add_laplace_bc_rhs(decomp, bc_b(layer), bc_d(layer), rhs(:, :, layer))
    end do
    
  end subroutine add_laplace_bc_rhs_periodic_layers_multiple
  
  ! Apply Dirichlet boundary conditions.
  subroutine apply_dirichlet_bc_layer(decomp, bc, zeta)
    type(decomposition), intent(in) :: decomp
    real(kind = real_8), intent(in) :: bc
    real(kind = real_8), dimension(:, :), intent(inout) :: zeta
    
    integer :: i, j
    
    assert(size(zeta, 1) == decomp%ii)
    assert(size(zeta, 2) == decomp%jj)
    
    if(decomp%p_a == -1) then
      do i = 1, decomp%ii
        zeta(i, decomp%jj) = bc
      end do
    end if
    if(decomp%p_c == -1) then
      do i = 1, decomp%ii
        zeta(i, 1) = bc
      end do
    end if
    if(decomp%p_b == -1) then
      do j = 1, decomp%jj
        zeta(1, j) = bc
      end do
    end if
    if(decomp%p_d == -1) then
      do j = 1, decomp%jj
        zeta(decomp%ii, j) = bc
      end do
    end if
    
  end subroutine apply_dirichlet_bc_layer
  
  ! Apply Dirichlet boundary conditions.
  subroutine apply_dirichlet_bc_layers_single_data(decomp, bc, zeta)
    type(decomposition), intent(in) :: decomp
    real(kind = real_8), intent(in) :: bc
    real(kind = real_8), dimension(:, :, :), intent(inout) :: zeta
    
    integer :: layer
    
    assert(size(zeta, 1) == decomp%ii)
    assert(size(zeta, 2) == decomp%jj)
    assert(size(zeta, 3) == decomp%layers)
  
    do layer = 1, decomp%layers
      call apply_dirichlet_bc(decomp, bc, zeta(:, :, layer))
    end do
    
  end subroutine apply_dirichlet_bc_layers_single_data
  
  ! Apply Dirichlet boundary conditions. This does *not* update the halo region.
  subroutine apply_dirichlet_bc_layers_single_pdata(bc, zeta)
    real(kind = real_8), intent(in) :: bc
    type(parallel_data), intent(inout) :: zeta
    
    call apply_dirichlet_bc(zeta%decomp, bc, zeta%sdata)
    
  end subroutine apply_dirichlet_bc_layers_single_pdata
  
  ! Apply Dirichlet boundary conditions.
  subroutine apply_dirichlet_bc_layers_multiple_data(decomp, bc, zeta)
    type(decomposition), intent(in) :: decomp
    real(kind = real_8), dimension(decomp%layers), intent(in) :: bc
    real(kind = real_8), dimension(:, :, :), intent(inout) :: zeta
    
    integer :: layer
    
    assert(size(zeta, 1) == decomp%ii)
    assert(size(zeta, 2) == decomp%jj)
    assert(size(zeta, 3) == decomp%layers)
  
    do layer = 1, decomp%layers
      call apply_dirichlet_bc(decomp, bc(layer), zeta(:, :, layer))
    end do
    
  end subroutine apply_dirichlet_bc_layers_multiple_data
  
  ! Apply Dirichlet boundary conditions. This does *not* update the halo region.
  subroutine apply_dirichlet_bc_layers_multiple_pdata(bc, zeta)
    type(parallel_data), intent(inout) :: zeta
    real(kind = real_8), dimension(zeta%decomp%layers), intent(in) :: bc
    
    call apply_dirichlet_bc(zeta%decomp, bc, zeta%sdata)
    
  end subroutine apply_dirichlet_bc_layers_multiple_pdata
  
  ! Apply Dirichlet boundary conditions.
  subroutine apply_dirichlet_bc_periodic_layer_data(decomp, bc_b, bc_d, zeta)
    type(decomposition), intent(in) :: decomp
    real(kind = real_8), intent(in) :: bc_b
    real(kind = real_8), intent(in) :: bc_d
    real(kind = real_8), dimension(:, :), intent(inout) :: zeta
    
    integer :: j
    
    assert(size(zeta, 1) == decomp%ii)
    assert(size(zeta, 2) == decomp%jj)
    
    if(decomp%p_b == -1) then
      do j = 1, decomp%jj
        zeta(1, j) = bc_b
      end do
    end if
    if(decomp%p_d == -1) then
      do j = 1, decomp%jj
        zeta(decomp%ii, j) = bc_d
      end do
    end if
    
  end subroutine apply_dirichlet_bc_periodic_layer_data
  
  ! Apply Dirichlet boundary conditions. This does *not* update the halo region.
  subroutine apply_dirichlet_bc_periodic_layer_pdata(bc_b, bc_d, zeta, layer)
    real(kind = real_8), intent(in) :: bc_b
    real(kind = real_8), intent(in) :: bc_d
    type(parallel_data), intent(inout) :: zeta
    integer, intent(in) :: layer
    
    assert(layer >= 1 .and. layer <= zeta%decomp%layers)
    
    call apply_dirichlet_bc(zeta%decomp, bc_b, bc_d, zeta%sdata(:, :, layer))
    
  end subroutine apply_dirichlet_bc_periodic_layer_pdata
  
  ! Apply Dirichlet boundary conditions.
  subroutine apply_dirichlet_bc_periodic_layers_single_data(decomp, bc_b, bc_d, zeta)
    type(decomposition), intent(in) :: decomp
    real(kind = real_8), intent(in) :: bc_b
    real(kind = real_8), intent(in) :: bc_d
    real(kind = real_8), dimension(:, :, :), intent(inout) :: zeta
    
    integer :: layer
    
    assert(size(zeta, 1) == decomp%ii)
    assert(size(zeta, 2) == decomp%jj)
    assert(size(zeta, 3) == decomp%layers)
  
    do layer = 1, decomp%layers
      call apply_dirichlet_bc(decomp, bc_b, bc_d, zeta(:, :, layer))
    end do
    
  end subroutine apply_dirichlet_bc_periodic_layers_single_data
  
  ! Apply Dirichlet boundary conditions. This does *not* update the halo region.
  subroutine apply_dirichlet_bc_periodic_layers_single_pdata(bc_b, bc_d, zeta)
    real(kind = real_8), intent(in) :: bc_b
    real(kind = real_8), intent(in) :: bc_d
    type(parallel_data), intent(inout) :: zeta
    
    call apply_dirichlet_bc(zeta%decomp, bc_b, bc_d, zeta%sdata)
    
  end subroutine apply_dirichlet_bc_periodic_layers_single_pdata
  
  ! Apply Dirichlet boundary conditions.
  subroutine apply_dirichlet_bc_periodic_layers_multiple_data(decomp, bc_b, bc_d, zeta)
    type(decomposition), intent(in) :: decomp
    real(kind = real_8), dimension(decomp%layers), intent(in) :: bc_b
    real(kind = real_8), dimension(decomp%layers), intent(in) :: bc_d
    real(kind = real_8), dimension(:, :, :), intent(inout) :: zeta
    
    integer :: layer
    
    assert(size(zeta, 1) == decomp%ii)
    assert(size(zeta, 2) == decomp%jj)
    assert(size(zeta, 3) == decomp%layers)
  
    do layer = 1, decomp%layers
      call apply_dirichlet_bc(decomp, bc_b(layer), bc_d(layer), zeta(:, :, layer))
    end do
    
  end subroutine apply_dirichlet_bc_periodic_layers_multiple_data
  
  ! Apply Dirichlet boundary conditions. This does *not* update the halo region.
  subroutine apply_dirichlet_bc_periodic_layers_multiple_pdata(bc_b, bc_d, zeta)
    type(parallel_data), intent(inout) :: zeta
    real(kind = real_8), dimension(zeta%decomp%layers), intent(in) :: bc_b
    real(kind = real_8), dimension(zeta%decomp%layers), intent(in) :: bc_d
    
    call apply_dirichlet_bc(zeta%decomp, bc_b, bc_d, zeta%sdata)
    
  end subroutine apply_dirichlet_bc_periodic_layers_multiple_pdata
  
end module qg_laplace_bc
