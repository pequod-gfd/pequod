!!NP 4

#include "config.h"

program test_transpose_pdata_i_4

  use qg_parallel
  use qg_parameters, only : real_8
  use test_tools

  implicit none
  
  integer :: i, j, n
  real(kind = real_8) :: m
  type(decomposition) :: decomp
  type(parallel_data) :: donor, target
  
  call parallel_init()
  
  n = nprocs()
  call preport_test("nprocs == 4", n /= 4)
  
  call allocate(decomp, iis = (/2, 2, 2, 3/), jjs = (/4/), layers = 1)
  call allocate(donor, decomp)
  call allocate(target, decomp%trans_i_decomp)
  
  call set(donor, -1.0D0)
  do j = 1, decomp%jj
    do i = 1, decomp%ii
      donor%sdata(i, j, 1) = real(j + decomp%jj_0 + (i + decomp%ii_0 - 1) * decomp%g_jj, kind = real_8)
    end do
  end do
  call zero(target)
  
  call ptranspose(donor, target, TRANSPOSE_II)
  call flush(donor, TRANSPOSE_II)
  
  m = 0.0D0
  do j = 1, target%decomp%jj
    do i = 1, target%decomp%ii
      m = max(m, dabs(target%sdata(i, j, 1) - (j + target%decomp%jj_0 + (i + target%decomp%ii_0 - 1) * target%decomp%g_jj)))
    end do
  end do
  call preport_test("transpose_i: 4 x 1", fne(m, 0.0D0))
  
  call deallocate(donor)
  call deallocate(target)
  call deallocate(decomp)

  call allocate(decomp, iis = (/128, 128, 128, 129/), jjs = (/513/), layers = 1)
  call allocate(donor, decomp)
  call allocate(target, decomp%trans_i_decomp)
  
  call set(donor, -1.0D0)
  do j = 1, decomp%jj
    do i = 1, decomp%ii
      donor%sdata(i, j, 1) = real(j + decomp%jj_0 + (i + decomp%ii_0 - 1) * decomp%g_jj, kind = real_8)
    end do
  end do
  call zero(target)
  
  call ptranspose(donor, target, TRANSPOSE_II)
  call flush(donor, TRANSPOSE_II)
  
  m = 0.0D0
  do j = 1, target%decomp%jj
    do i = 1, target%decomp%ii
      m = max(m, dabs(target%sdata(i, j, 1) - (j + target%decomp%jj_0 + (i + target%decomp%ii_0 - 1) * target%decomp%g_jj)))
    end do
  end do
  call preport_test("transpose_i: 4 x 1", fne(m, 0.0D0))
  
  call deallocate(donor)
  call deallocate(target)
  call deallocate(decomp)
  
  call allocate(decomp, iis = (/3, 4/), jjs = (/3, 4/), layers = 1)
  call allocate(donor, decomp)
  call allocate(target, decomp%trans_i_decomp)
  
  call set(donor, -1.0D0)
  do j = 1, decomp%jj
    do i = 1, decomp%ii
      donor%sdata(i, j, 1) = real(j + decomp%jj_0 + (i + decomp%ii_0 - 1) * decomp%g_jj, kind = real_8)
    end do
  end do
  call zero(target)
  
  call ptranspose(donor, target, TRANSPOSE_II)
  call flush(donor, TRANSPOSE_II)
  
  m = 0.0D0
  do j = 1, target%decomp%jj
    do i = 1, target%decomp%ii
      m = max(m, dabs(target%sdata(i, j, 1) - (j + target%decomp%jj_0 + (i + target%decomp%ii_0 - 1) * target%decomp%g_jj)))
    end do
  end do
  call preport_test("transpose_i: 2 x 2", fne(m, 0.0D0))
  
  call deallocate(donor)
  call deallocate(target)
  call deallocate(decomp)
  
  call allocate(decomp, iis = (/256, 257/), jjs = (/256, 257/), layers = 1)
  call allocate(donor, decomp)
  call allocate(target, decomp%trans_i_decomp)
  
  call set(donor, -1.0D0)
  do j = 1, decomp%jj
    do i = 1, decomp%ii
      donor%sdata(i, j, 1) = real(j + decomp%jj_0 + (i + decomp%ii_0 - 1) * decomp%g_jj, kind = real_8)
    end do
  end do
  call zero(target)
  
  call ptranspose(donor, target, TRANSPOSE_II)
  call flush(donor, TRANSPOSE_II)
  
  m = 0.0D0
  do j = 1, target%decomp%jj
    do i = 1, target%decomp%ii
      m = max(m, dabs(target%sdata(i, j, 1) - (j + target%decomp%jj_0 + (i + target%decomp%ii_0 - 1) * target%decomp%g_jj)))
    end do
  end do
  call preport_test("transpose_i: 2 x 2", fne(m, 0.0D0))
  
  call deallocate(donor)
  call deallocate(target)
  call deallocate(decomp)

  call parallel_finalise()

end program test_transpose_pdata_i_4
