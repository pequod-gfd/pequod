!!NP 2

program test_update_flux_i_2

  use qg_parallel
  use qg_parameters, only : real_8
  use test_tools

  implicit none
  
  integer :: j, n
  real(kind = real_8) :: m
  type(decomposition) :: decomp
  type(parallel_flux_data) :: fpdata
  
  call parallel_init()
  
  n = nprocs()
  call preport_test("nprocs == 2", n /= 2)
  
  call allocate(decomp, iis = (/4, 4/), jjs = (/6/), layers = 1)
  
  call allocate(fpdata, decomp, FLUX_II)
  call zero(fpdata)
  
  do j = decomp%fa_j_min + 1, decomp%fa_j_max - 1
    fpdata%data(decomp%fa_i_min + 2, j, 1) = real(j, kind = real_8)
    fpdata%data(decomp%fa_i_max - 1, j, 1) = -real(j, kind = real_8)
  end do
  
  call update(fpdata)
  call flush()
  
  m = 0.0D0
  if(decomp%p_d >= 0) then
    do j = decomp%fa_j_min + 1, decomp%fa_j_max - 1
      m = max(m, dabs(fpdata%data(decomp%fa_i_max + 1, j, 1) - real(j, kind = real_8)))
    end do
  end if
  if(decomp%p_b >= 0) then
    do j = decomp%fa_j_min + 1, decomp%fa_j_max - 1
      m = max(m, dabs(fpdata%data(decomp%fa_i_min, j, 1) + real(j, kind = real_8)))
    end do
  end if
  call preport_test("update_flux_i", fne(m, 0.0D0))
  
  call zero(fpdata)
  
  do j = decomp%fa_j_min + 1, decomp%fa_j_max - 1
    fpdata%data(decomp%fa_i_min + 1, j, 1) = real(j + decomp%g_jj, kind = real_8)
    fpdata%data(decomp%fa_i_min + 2, j, 1) = real(j, kind = real_8)
    fpdata%data(decomp%fa_i_max - 1, j, 1) = -real(j, kind = real_8)
  end do
  
  call update(fpdata, level_1 = .true.)
  call flush()
  
  m = 0.0D0
  if(decomp%p_d >= 0) then
    do j = decomp%fa_j_min + 1, decomp%fa_j_max - 1
      m = max(m, dabs(fpdata%data(decomp%fa_i_max, j, 1) - real(j + decomp%g_jj, kind = real_8)))
      m = max(m, dabs(fpdata%data(decomp%fa_i_max + 1, j, 1) - real(j, kind = real_8)))
    end do
  end if
  if(decomp%p_b >= 0) then
    do j = decomp%fa_j_min + 1, decomp%fa_j_max - 1
      m = max(m, dabs(fpdata%data(decomp%fa_i_min, j, 1) + real(j, kind = real_8)))
    end do
  end if
  call preport_test("update_flux_i", fne(m, 0.0D0))
  
  call deallocate(fpdata)
  
  call deallocate(decomp)
  
  call parallel_finalise()

end program test_update_flux_i_2
