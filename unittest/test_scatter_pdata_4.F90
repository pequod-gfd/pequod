!!NP 4

#include "config.h"

program test_scatter_pdata_4

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
    call zero(pdata)
    if(decomp%p + 1 == root) then
      do j = 1, decomp%g_jj
        do i = 1, decomp%g_ii
          data(i, j, 1) = i + decomp%g_ii * (j - 1)
        end do
      end do
      
      call scatter(pdata, data)
    else
      call scatter(pdata, root)
    end if
  
    m = 0.0D0
    do j = decomp%fa_j_min, decomp%fa_j_max
      do i = decomp%fa_i_min, decomp%fa_i_max
        m = max(m, dabs(pdata%data(i, j, 1) - real((i - 1 + decomp%ii_0) + decomp%g_ii * ((j - 1 + decomp%jj_0) - 1), kind = real_8)))
      end do
    end do
    call preport_test("scatter: root " // int2str(root), fne(m, 0.0D0))
  end do

  deallocate(data)
  call deallocate(pdata)
  call deallocate(decomp)
  
  call parallel_finalise()

end program test_scatter_pdata_4
