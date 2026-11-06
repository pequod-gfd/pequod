!!NP 4

#include "config.h"

program test_gather_pdata_4

  use qg_parallel
  use qg_parameters, only : real_8
  use test_tools
  
  implicit none

  integer :: i, j, n, root
  real(kind = real_8) :: m
  real(kind = real_8), dimension(:, :, :), allocatable :: data
  type(decomposition) :: decomp
  type(parallel_data) :: pdata
  
  call parallel_init()
  
  n = nprocs()
  call preport_test("nprocs == 4", n /= 4)
  
  call allocate(decomp, iis = (/3, 4/), jjs = (/5, 6/), layers = 1)
  call allocate(pdata, decomp)  
  allocate(data(decomp%g_ii, decomp%g_jj, 1))
  
  do root = 1, n
    call set(pdata, -1.0D0)
    do j = 1, decomp%jj
      do i = 1, decomp%ii
        pdata%sdata(i, j, 1) = real(j + decomp%jj_0 + (i + decomp%ii_0 - 1) * decomp%g_jj, kind = real_8)
      end do
    end do
  
    m = 0.0D0
    if(decomp%p + 1 == root) then
      data = 0.0D0
      call gather(pdata, data)
      
      do j = 1, decomp%g_jj
        do i = 1, decomp%g_ii
          m = max(m, dabs(data(i, j, 1) - real(j + (i - 1) * decomp%g_jj, kind = real_8)))
        end do
      end do
    else
      call gather(pdata, root)
    end if
    call preport_test("gather: root " // int2str(root), fne(m, 0.0D0))
  end do

  deallocate(data)
  call deallocate(pdata)
  call deallocate(decomp)
  
  call parallel_finalise()

end program test_gather_pdata_4
