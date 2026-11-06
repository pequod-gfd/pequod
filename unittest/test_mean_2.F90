!!NP 2

program test_mean_2

  use qg_diagnostics
  use qg_parallel
  use qg_parameters, only : real_8
  use test_tools

  implicit none
  
  integer :: i, j, n
  real(kind = real_8) :: m
  real(kind = real_8), dimension(5, 5) :: arr
  type(decomposition) :: decomp
  type(parallel_data) :: pdata
  
  call parallel_init()
  
  n = nprocs()
  call preport_test("nprocs == 2", n /= 2)
  
  call allocate(decomp, iis = (/5/), jjs = (/2, 3/), layers = 1)
  
  arr = 1.0D0
  m = mean(decomp, arr(:decomp%ii, :decomp%jj))
  call preport_test("mean: 5 x 5 unity valued field", fne(m, 1.0D0))
  
  arr(decomp%l_i_min:decomp%l_i_max, decomp%l_j_min:decomp%l_j_max) = 0.0D0
  m = mean(decomp, arr(:decomp%ii, :decomp%jj))
  call preport_test("mean: 5 x 5", fne(m, 7.0D0 / 16.0D0))
  
  call allocate(pdata, decomp)
  
  call set(pdata, 1.0D0)
  m = mean(pdata, 1)
  call preport_test("mean: 5 x 5 unity valued field", fne(m, 1.0D0))
  
  pdata%sdata(decomp%l_i_min:decomp%l_i_max, decomp%l_j_min:decomp%l_j_max, 1) = 0.0D0
  m = mean(pdata, 1)
  call preport_test("mean: 5 x 5", fne(m, 7.0D0 / 16.0D0))
  
  call deallocate(pdata)
  
  call deallocate(decomp)
  
  call allocate(decomp, iis = (/2, 3/), jjs = (/5/), layers = 1)
  
  arr = 1.0D0
  m = mean(decomp, arr(:decomp%ii, :decomp%jj))
  call preport_test("mean: 5 x 5 unity valued field", fne(m, 1.0D0))
  
  arr(decomp%l_i_min:decomp%l_i_max, decomp%l_j_min:decomp%l_j_max) = 0.0D0
  m = mean(decomp, arr(:decomp%ii, :decomp%jj))
  call preport_test("mean: 5 x 5", fne(m, 7.0D0 / 16.0D0))
  
  do i = 1, decomp%ii
    arr(i, :) = real(decomp%ii_0 + i - 1, kind = real_8) / real(decomp%g_ii_int, kind = real_8)
  end do
  m = mean(decomp, arr(:decomp%ii, :decomp%jj))
  call preport_test("mean: 5 x 5", fne(m, 0.5D0))
  
  do j = 1, decomp%jj
    arr(:, j) = real(decomp%jj_0 + j - 1, kind = real_8) / real(decomp%g_jj_int, kind = real_8)
  end do
  m = mean(decomp, arr(:decomp%ii, :decomp%jj))
  call preport_test("mean: 5 x 5", fne(m, 0.5D0))
  
  call allocate(pdata, decomp)
  
  call set(pdata, 1.0D0)
  m = mean(pdata, 1)
  call preport_test("mean: 5 x 5 unity valued field", fne(m, 1.0D0))
  
  pdata%sdata(decomp%l_i_min:decomp%l_i_max, decomp%l_j_min:decomp%l_j_max, 1) = 0.0D0
  m = mean(pdata, 1)
  call preport_test("mean: 5 x 5", fne(m, 7.0D0 / 16.0D0))
  
  call zero(pdata)
  do i = decomp%fa_i_min, decomp%fa_i_max
    pdata%data(i, :, 1) = real(decomp%ii_0 + i - 2, kind = real_8) / real(decomp%g_ii_int, kind = real_8)
  end do
  m = mean(pdata, 1)  
  call preport_test("mean: 5 x 5", fne(m, 0.5D0))
  
  call zero(pdata)
  do j = decomp%fa_j_min, decomp%fa_j_max
    pdata%data(:, j, 1) = real(decomp%jj_0 + j - 2, kind = real_8) / real(decomp%g_jj_int, kind = real_8)
  end do
  m = mean(pdata, 1)  
  call preport_test("mean: 5 x 5", fne(m, 0.5D0))
  
  call deallocate(pdata)
  
  call deallocate(decomp)
  
  call parallel_finalise()

end program test_mean_2
