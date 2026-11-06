!!NP 2

program test_update_flux_j_2

  use qg_parallel
  use qg_parameters, only : real_8
  use test_tools

  implicit none
  
  integer :: i, n
  real(kind = real_8) :: m
  type(decomposition) :: decomp
  type(parallel_flux_data) :: fpdata
  
  call parallel_init()
  
  n = nprocs()
  call preport_test("nprocs == 2", n /= 2)
  
  call allocate(decomp, iis = (/6/), jjs = (/4, 4/), layers = 1)
  
  call allocate(fpdata, decomp, FLUX_JJ)
  call zero(fpdata)
  
  do i = decomp%fa_i_min + 1, decomp%fa_i_max - 1
    fpdata%data(i, decomp%fa_j_min + 2, 1) = real(i, kind = real_8)
    fpdata%data(i, decomp%fa_j_max - 1, 1) = -real(i, kind = real_8)
  end do
  
  call update(fpdata)
  call flush()
  
  m = 0.0D0
  if(decomp%p_c >= 0) then
    do i = decomp%fa_i_min + 1, decomp%fa_i_max - 1
      m = max(m, dabs(fpdata%data(i, decomp%fa_j_min, 1) + real(i, kind = real_8)))
    end do
  end if
  if(decomp%p_a >= 0) then
    do i = decomp%fa_i_min + 1, decomp%fa_i_max - 1
      m = max(m, dabs(fpdata%data(i, decomp%fa_j_max + 1, 1) - real(i, kind = real_8)))
    end do
  end if
  call preport_test("update_flux_j", fne(m, 0.0D0))
  
  call zero(fpdata)
  
  do i = decomp%fa_i_min + 1, decomp%fa_i_max - 1
    fpdata%data(i, decomp%fa_j_min + 1, 1) = real(i + decomp%g_ii, kind = real_8)
    fpdata%data(i, decomp%fa_j_min + 2, 1) = real(i, kind = real_8)
    fpdata%data(i, decomp%fa_j_max - 1, 1) = -real(i, kind = real_8)
  end do
  
  call update(fpdata, level_1 = .true.)
  call flush()
  
  m = 0.0D0
  if(decomp%p_c >= 0) then
    do i = decomp%fa_i_min + 1, decomp%fa_i_max - 1
      m = max(m, dabs(fpdata%data(i, decomp%fa_j_min, 1) + real(i, kind = real_8)))
    end do
  end if
  if(decomp%p_a >= 0) then
    do i = decomp%fa_i_min + 1, decomp%fa_i_max - 1
      m = max(m, dabs(fpdata%data(i, decomp%fa_j_max, 1) - real(i + decomp%g_ii, kind = real_8)))
      m = max(m, dabs(fpdata%data(i, decomp%fa_j_max + 1, 1) - real(i, kind = real_8)))
    end do
  end if
  call preport_test("update_flux_j", fne(m, 0.0D0))
  
  call deallocate(fpdata)
  
  call deallocate(decomp)
  
  call parallel_finalise()

end program test_update_flux_j_2
