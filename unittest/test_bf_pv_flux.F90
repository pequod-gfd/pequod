#include "config.h"

program test_bf_pv_flux

  use qg_diagnostics
  use qg_parallel
  use qg_parameters, only : real_8
  use test_tools

  implicit none
  
  integer :: i, j, layer
  real(kind = real_8), parameter :: Rd_bcl1 = 20.0D5, Rd_bcl2 = 10.0D5, &
    & visc_bot = 42.0D0
  real(kind = real_8), dimension(3), parameter :: H = (/1.0D5, 2.0D5, 3.0D5/)
  real(kind = real_8), dimension(:, :), allocatable :: f_m, f_z
  real(kind = real_8), dimension(:, :, :), allocatable :: psi_m, psi_z
  type(decomposition) :: decomp
  type(stratification) :: strat
  
  call parallel_init()
  
  call allocate(decomp, iis = (/4/), jjs = (/128/), layers = 3, periodic_j = .true.)
  call allocate(strat, Rd_bcl1, Rd_bcl2, H(1), H(2), H(3))
  allocate(psi_z(decomp%ii, decomp%jj, decomp%layers))
  allocate(psi_m(decomp%ii, decomp%jj, decomp%layers))
  allocate(f_z(decomp%layers, 2))
  allocate(f_m(decomp%layers, 2))
  
  do layer = 1, decomp%layers
    do j = 1, decomp%jj
      do i = 1, decomp%ii
        psi_z(i, j, layer) = real(decomp%ii - i + 1, kind = real_8)
      end do
    end do
  end do
  call z_to_modes(strat, psi_z, psi_m)
  
  call bf_pv_flux(decomp, strat, psi_m, f_m(:, 1), f_m(:, 2), visc_bot = visc_bot)
  call modes_to_z(strat, f_m(:, 1), f_z(:, 1))
  call modes_to_z(strat, f_m(:, 2), f_z(:, 2))
  
  call report_test("bf_pv_flux: f_b", fne(f_z(:, 1), (/0.0D0, 0.0D0, -visc_bot * real(decomp%jj, kind = real_8)/), tol = 1.0D-10))
  call report_test("bf_pv_flux: f_d", fne(f_z(:, 2), (/0.0D0, 0.0D0,  visc_bot * real(decomp%jj, kind = real_8)/), tol = 1.0D-10))
    
  deallocate(f_m)
  deallocate(f_z)
  deallocate(psi_m)
  deallocate(psi_z)
  call deallocate(strat)
  call deallocate(decomp)
  
  call parallel_finalise()

end program test_bf_pv_flux
