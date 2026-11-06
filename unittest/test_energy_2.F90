!!NP 2

program test_energy_2
  
  use qg_diagnostics
  use qg_parallel
  use qg_parameters, only : real_8
  use test_tools
  
  implicit none
  
  integer :: i, n
  real(kind = real_8), dimension(2) :: ekin
  real(kind = real_8) :: epot
  type(decomposition) :: decomp
  type(parallel_data) :: psi
  type(stratification) :: strat
  
  call parallel_init()
  
  n = nprocs()
  call preport_test("nprocs == 2", n /= 2)
  
  call allocate(strat, dsqrt(0.5D0), 1.0D0, 1.0D0)
  call allocate(decomp, iis = (/5/), jjs = (/4, 5/), layers = 2)
  call allocate(psi, decomp)
  
  call zero(psi)
  call energy(strat, psi, ekin, epot)
  call preport_test("energy: zero valued psi", fne(epot, 0.0D0) .and. fne(ekin, (/0.0D0, 0.0D0/)))
  
  call set(psi, 1.0D0)
  call energy(strat, psi, ekin, epot)
  call preport_test("energy: unity valued psi", fne(epot, 0.0D0) .and. fne(ekin, (/0.0D0, 0.0D0/)))
  
  call set(psi, 1.0D0, layer = 1)
  call set(psi, 0.0D0, layer = 2)
  call energy(strat, psi, ekin, epot)
  call preport_test("energy: psi constant in layers", fne(epot, 0.25D0) .and. fne(ekin, (/0.0D0, 0.0D0/)))
  
  call set(psi, 0.0D0, layer = 1)
  call set(psi, 1.0D0, layer = 2)
  call energy(strat, psi, ekin, epot)
  call preport_test("energy: psi constant in layers", fne(epot, 0.25D0) .and. fne(ekin, (/0.0D0, 0.0D0/)))
  
  call set(psi, 2.0D0, layer = 1)
  call set(psi, 0.0D0, layer = 2)
  call energy(strat, psi, ekin, epot)
  call preport_test("energy: psi constant in layers", fne(epot, 1.0D0) .and. fne(ekin, (/0.0D0, 0.0D0/)))
  
  call zero(psi)
  do i = decomp%fa_i_min, decomp%fa_i_max
    psi%data(i, :, 1) = real(decomp%ii_0 + i - 1, kind = real_8)
    psi%data(i, :, 2) = real(decomp%ii_0 + i - 1, kind = real_8) + 1.0D0
  end do
  call energy(strat, psi, ekin, epot)
  call preport_test("energy: linearly varying psi in layers", fne(epot, 0.25D0) .and. fne(ekin, (/0.25D0, 0.25D0/)))
  
  call zero(psi)
  do i = decomp%fa_i_min, decomp%fa_i_max
    psi%data(i, :, 1) = 2.0D0 * real(decomp%ii_0 + i - 1, kind = real_8)
    psi%data(i, :, 2) = 2.0D0 * real(decomp%ii_0 + i - 1, kind = real_8) + 2.0D0
  end do
  call energy(strat, psi, ekin, epot)
  call preport_test("energy: linearly varying psi in layers", fne(epot, 1.0D0) .and. fne(ekin, (/1.0D0, 1.0D0/)))
  
  call deallocate(psi)
  call deallocate(decomp)
  call deallocate(strat)
  
  call parallel_finalise()
  
end program test_energy_2
