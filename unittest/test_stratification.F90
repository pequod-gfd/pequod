program test_stratification
  
  use lapack_interfaces
  use qg_diagnostics
  use qg_parallel
  use qg_parameters, only : real_8
  use test_tools

  implicit none
  
  real(kind = real_8), dimension(:), allocatable :: m, z
  real(kind = real_8), parameter :: rd_bcl1 = 20.0D5, rd_bcl2 = 10.0D5
  real(kind = real_8), dimension(3), parameter :: h = (/1.0D5, 2.0D5, 3.0D5/)
  type(stratification) :: strat
  
  call parallel_init()
  
  ! 1 layer
  call allocate(strat, h(1))
  
  call report_test("s_matrix: 1 layer", fne(strat%ss(1), 0.0D0))
  
  call deallocate(strat)
  
  ! 2 layer
  call allocate(strat, rd_bcl1, h(1), h(2))
  
  call report_test("s_matrix: 2 layer", fne(1.0D0 / dsqrt(strat%ss(2)), rd_bcl1))
  call report_test("s_matrix: 2 layer", fne(strat%ss(1), 0.0D0))
    
  allocate(m(2), z(2))
  z = 1.0D0
  call z_to_modes(strat, z, m)
  call report_test("z_to_modes: 2 layer", fne(m, (/1.0D0, 0.0D0/)))
  call modes_to_z(strat, m, z)
  call report_test("modes_to_z: 2 layer", fne(z, (/1.0D0, 1.0D0/)))
  deallocate(m, z)
  
  call deallocate(strat)
  
  ! 3 layer
  call allocate(strat, rd_bcl1, rd_bcl2, h(1), h(2), h(3))
  
  call report_test("s_matrix: 3 layer", fne(1.0D0 / dsqrt(strat%ss(2:3)), (/rd_bcl2, rd_bcl1/)))
  call report_test("s_matrix: 3 layer", fne(strat%ss(1), 0.0D0))
    
  allocate(m(3), z(3))
  z = 1.0D0
  call z_to_modes(strat, z, m)
  call report_test("z_to_modes: 3 layer", fne(m, (/1.0D0, 0.0D0, 0.0D0/)))
  call modes_to_z(strat, m, z)
  call report_test("modes_to_z: 3 layer", fne(z, (/1.0D0, 1.0D0, 1.0D0/)))
  deallocate(m, z)
  
  call deallocate(strat)
  
  call parallel_finalise()

end program test_stratification
