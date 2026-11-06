#include "config.h"

program test_laplace_bc

  use qg_helmholtz
  use qg_laplace_bc
  use qg_parallel
  use qg_parameters, only : real_8
  use test_tools

  implicit none
  
  real(kind = real_8), dimension(:, :, :), allocatable :: phi, phi0, phi1
  type(decomposition) :: decomp
  type(helm_cache) :: hcache
  integer :: i, j
  real(kind = real_8) :: m
  
  call parallel_init()
  
  call allocate(decomp, iis = (/9/), jjs = (/9/), layers = 1)
  allocate(phi(decomp%ii, decomp%jj, decomp%layers))
  allocate(phi0(decomp%ii, decomp%jj, decomp%layers))
  allocate(phi1(decomp%ii, decomp%jj, decomp%layers))
  call allocate(hcache, decomp, ss = (/0.0D0/))
  
  phi(:, :, 1) = laplace_bc_rhs(decomp, 1.0D0)
  call helm(hcache, phi)
  call apply_dirichlet_bc(decomp, 1.0D0, phi(:, :, 1))
  
  call report_test("laplace_bc_rhs", fne(minval(phi), 1.0D0) .or. fne(maxval(phi), 1.0D0))
  
  phi0 = 1.0D0
  !call add_laplace_bc(decomp, 0.0D0, phi0(:, :, 1))
  call helm(hcache, phi0)
  call apply_dirichlet_bc(decomp, 0.0D0, phi0(:, :, 1))
  
  phi1 = 1.0D0
  call add_laplace_bc_rhs(decomp, 1.0D0, phi1(:, :, 1))
  call helm(hcache, phi1)
  call apply_dirichlet_bc(decomp, 1.0D0, phi1(:, :, 1))
  
  phi = phi1 - phi0
  call report_test("add_laplace_bc_rhs", fne(minval(phi), 1.0D0) .or. fne(maxval(phi), 1.0D0))
  
  call deallocate(hcache)
  deallocate(phi)
  deallocate(phi0)
  deallocate(phi1)
  call deallocate(decomp)

  call allocate(decomp, iis = (/9/), jjs = (/9/), layers = 1, periodic_j = .true.)
  allocate(phi(decomp%ii, decomp%jj, decomp%layers))
  call allocate(hcache, decomp, ss = (/0.0D0/))
  
  phi = 0.0D0
  call add_laplace_bc_rhs(decomp, 1.0D0, 0.0D0, phi(:, :, 1))
  call helm(hcache, phi)
  call apply_dirichlet_bc(decomp, 1.0D0, 0.0D0, phi(:, :, 1))
  
  m = 0.0D0
  do j = 1, decomp%jj
    do i = 1, decomp%ii
      m = max(m, dabs(phi(i, j, 1) - real(decomp%ii - i, kind = real_8) / real(decomp%ii - 1, kind = real_8)))
    end do
  end do
  call report_test("add_laplace_bc_rhs: periodic_j", fne(m, 0.0D0))
  
  phi = 0.0D0
  call add_laplace_bc_rhs(decomp, 0.0D0, 1.0D0, phi(:, :, 1))
  call helm(hcache, phi)
  call apply_dirichlet_bc(decomp, 0.0D0, 1.0D0, phi(:, :, 1))
  
  m = 0.0D0
  do j = 1, decomp%jj
    do i = 1, decomp%ii
      m = max(m, dabs(phi(i, j, 1) - real(i - 1, kind = real_8) / real(decomp%ii - 1, kind = real_8)))
    end do
  end do
  call report_test("add_laplace_bc_rhs: periodic_j", fne(m, 0.0D0))
  
  call deallocate(hcache)
  deallocate(phi)
  call deallocate(decomp)

  call parallel_finalise()
  
end program test_laplace_bc
