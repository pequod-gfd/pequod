!!NP 2

#include "config.h"

program test_inverse_transpose_pdata_j_2

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
  call preport_test("nprocs == 2", n /= 2)
  
  call allocate(decomp, iis = (/2/), jjs = (/1, 2/), layers = 1)
  call allocate(target, decomp)  
  call allocate(donor, decomp%trans_j_decomp)
  
  call set(donor, -1.0D0)
  do j = 1, donor%decomp%jj
    do i = 1, donor%decomp%ii
      donor%sdata(i, j, 1) = real(j + donor%decomp%jj_0 + (i + donor%decomp%ii_0 - 1) * donor%decomp%g_jj, kind = real_8)
    end do
  end do
  call zero(target)
  
  call pinverse_transpose(donor, target, TRANSPOSE_JJ)
  call flush(target, TRANSPOSE_JJ)
  
  m = 0.0D0
  do j = 1, target%decomp%jj
    do i = 1, target%decomp%ii
      m = max(m, dabs(target%sdata(i, j, 1) - (j + target%decomp%jj_0 + (i + target%decomp%ii_0 - 1) * target%decomp%g_jj)))
    end do
  end do
  call preport_test("inverse_transpose_j: 1 x 2", fne(m, 0.0D0))
  
  call deallocate(donor)
  call deallocate(target)
  call deallocate(decomp)
  
  call allocate(decomp, iis = (/4/), jjs = (/3, 4/), layers = 1)
  call allocate(target, decomp)  
  call allocate(donor, decomp%trans_j_decomp)
  
  call set(donor, -1.0D0)
  do j = 1, donor%decomp%jj
    do i = 1, donor%decomp%ii
      donor%sdata(i, j, 1) = real(j + donor%decomp%jj_0 + (i + donor%decomp%ii_0 - 1) * donor%decomp%g_jj, kind = real_8)
    end do
  end do
  call zero(target)
  
  call pinverse_transpose(donor, target, TRANSPOSE_JJ)
  call flush(target, TRANSPOSE_JJ)
  
  m = 0.0D0
  do j = 1, target%decomp%jj
    do i = 1, target%decomp%ii
      m = max(m, dabs(target%sdata(i, j, 1) - (j + target%decomp%jj_0 + (i + target%decomp%ii_0 - 1) * target%decomp%g_jj)))
    end do
  end do
  call preport_test("inverse_transpose_j: 1 x 2", fne(m, 0.0D0))
  
  call deallocate(donor)
  call deallocate(target)
  call deallocate(decomp)
  
  call allocate(decomp, iis = (/513/), jjs = (/256, 257/), layers = 1)
  call allocate(target, decomp)  
  call allocate(donor, decomp%trans_j_decomp)
  
  call set(donor, -1.0D0)
  do j = 1, donor%decomp%jj
    do i = 1, donor%decomp%ii
      donor%sdata(i, j, 1) = real(j + donor%decomp%jj_0 + (i + donor%decomp%ii_0 - 1) * donor%decomp%g_jj, kind = real_8)
    end do
  end do
  call zero(target)
  
  call pinverse_transpose(donor, target, TRANSPOSE_JJ)
  call flush(target, TRANSPOSE_JJ)
  
  m = 0.0D0
  do j = 1, target%decomp%jj
    do i = 1, target%decomp%ii
      m = max(m, dabs(target%sdata(i, j, 1) - (j + target%decomp%jj_0 + (i + target%decomp%ii_0 - 1) * target%decomp%g_jj)))
    end do
  end do
  call preport_test("inverse_transpose_j: 1 x 2", fne(m, 0.0D0))
  
  call deallocate(donor)
  call deallocate(target)
  call deallocate(decomp)
  
  call parallel_finalise()

end program test_inverse_transpose_pdata_j_2
