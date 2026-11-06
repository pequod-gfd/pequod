program test_mean

  use qg_diagnostics
  use qg_parallel
  use qg_parameters, only : real_8
  use test_tools

  implicit none
  
  real(kind = real_8) :: m
  real(kind = real_8), dimension(3, 3) :: arr_3
  real(kind = real_8), dimension(4, 4) :: arr_4
  real(kind = real_8), dimension(5, 5) :: arr_5
  type(decomposition) :: decomp
  
  call parallel_init()
  
  call allocate(decomp, iis = (/3/), jjs = (/3/), layers = 1)
  
  arr_3 = 1.0D0
  m = mean(decomp, arr_3)
  call report_test("mean: 3 x 3 unity valued field", fne(m, 1.0D0))
  
  arr_3(2, 2) = 0.0D0
  m = mean(decomp, arr_3)
  call report_test("mean: 3 x 3", fne(m, 0.75D0))
  
  call deallocate(decomp)
  
  call allocate(decomp, iis = (/4/), jjs = (/4/), layers = 1)
  
  arr_4 = 1.0D0
  m = mean(decomp, arr_4)
  call report_test("mean: 4 x 4 unity valued field", fne(m, 1.0D0))
  
  arr_4(2:3, 2:3) = 0.0D0
  m = mean(decomp, arr_4)
  call report_test("mean: 4 x 4", fne(m, 5.0D0 / 9.0D0))
  
  call deallocate(decomp)
  
  call allocate(decomp, iis = (/5/), jjs = (/5/), layers = 1)
  
  arr_5 = 1.0D0
  m = mean(decomp, arr_5)
  call report_test("mean: 5 x 5 unity valued field", fne(m, 1.0D0))
  
  arr_5(2:4, 2:4) = 0.0D0
  m = mean(decomp, arr_5)
  call report_test("mean: 5 x 5", fne(m, 7.0D0 / 16.0D0))
  
  arr_5 = 1.0D0
  arr_5(decomp%l_i_min:decomp%l_i_max, decomp%l_j_min:decomp%l_j_max) = 0.0D0
  m = mean(decomp, arr_5)
  call report_test("mean: 5 x 5", fne(m, 7.0D0 / 16.0D0))
  
  call deallocate(decomp)
  
  call parallel_finalise()

end program test_mean
