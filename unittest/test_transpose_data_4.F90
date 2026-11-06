!!NP 4

#include "config.h"

program test_transpose_data_4

  use qg_parallel
  use qg_parameters, only : real_8
  use test_tools

  implicit none
  
  integer :: i, j, n
  real(kind = real_8) :: m
  real(kind = real_8), dimension(:, :, :), allocatable :: donor, target
  type(decomposition) :: donor_decomp
  type(decomposition), pointer :: target_decomp
  
  call parallel_init()
  
  n = nprocs()
  call preport_test("nprocs == 4", n /= 4)
  
  call allocate(donor_decomp, iis = (/128, 128, 128, 129/), jjs = (/513/), layers = 1)  
  allocate(donor(donor_decomp%ii + 2, donor_decomp%jj + 2, 1))
  
  target_decomp => donor_decomp%trans_i_decomp
  allocate(target(target_decomp%ii + 2, target_decomp%jj + 2, 1))
  
  donor = -1.0D0
  do j = 1, donor_decomp%jj
    do i = 1, donor_decomp%ii
      donor(i + 1, j + 1, 1) = real(j + donor_decomp%jj_0 + (i + donor_decomp%ii_0 - 1) * donor_decomp%g_jj, kind = real_8)
    end do
  end do
  target = 0.0D0
  
  call ptranspose(donor_decomp, donor, target, TRANSPOSE_II)
  call flush(donor_decomp, TRANSPOSE_II)
  
  m = 0.0D0
  do j = 1, target_decomp%jj
    do i = 1, target_decomp%ii
      m = max(m, dabs(target(i + 1, j + 1, 1) - (j + target_decomp%jj_0 + (i + target_decomp%ii_0 - 1) * target_decomp%g_jj)))
    end do
  end do
  call preport_test("transpose_i: 4 x 1", fne(m, 0.0D0))
  
  donor = 0.0D0
  
  call pinverse_transpose(donor_decomp, target, donor, TRANSPOSE_II)
  call flush(donor_decomp, TRANSPOSE_II)
  
  m = 0.0D0
  do j = 1, donor_decomp%jj
    do i = 1, donor_decomp%ii
      m = max(m, dabs(donor(i + 1, j + 1, 1) - (j + donor_decomp%jj_0 + (i + donor_decomp%ii_0 - 1) * donor_decomp%g_jj)))
    end do
  end do
  call preport_test("inverse_transpose_i: 4 x 1", fne(m, 0.0D0))
  
  deallocate(donor)
  deallocate(target)
  
  call deallocate(donor_decomp)
  
  call allocate(donor_decomp, iis = (/513/), jjs = (/128, 128, 128, 129/), layers = 1)  
  allocate(donor(donor_decomp%ii + 2, donor_decomp%jj + 2, 1))
  
  target_decomp => donor_decomp%trans_j_decomp  
  allocate(target(target_decomp%ii + 2, target_decomp%jj + 2, 1))
  
  donor = -1.0D0
  do j = 1, donor_decomp%jj
    do i = 1, donor_decomp%ii
      donor(i + 1, j + 1, 1) = real(j + donor_decomp%jj_0 + (i + donor_decomp%ii_0 - 1) * donor_decomp%g_jj, kind = real_8)
    end do
  end do
  target = 0.0D0
  
  call ptranspose(donor_decomp, donor, target, TRANSPOSE_JJ)
  call flush(donor_decomp, TRANSPOSE_JJ)
  
  m = 0.0D0
  do j = 1, target_decomp%jj
    do i = 1, target_decomp%ii
      m = max(m, dabs(target(i + 1, j + 1, 1) - (j + target_decomp%jj_0 + (i + target_decomp%ii_0 - 1) * target_decomp%g_jj)))
    end do
  end do
  call preport_test("transpose_j: 1 x 4", fne(m, 0.0D0))
  
  donor = 0.0D0
  
  call pinverse_transpose(donor_decomp, target, donor, TRANSPOSE_JJ)
  call flush(donor_decomp, TRANSPOSE_JJ)
  
  m = 0.0D0
  do j = 1, donor_decomp%jj
    do i = 1, donor_decomp%ii
      m = max(m, dabs(donor(i + 1, j + 1, 1) - (j + donor_decomp%jj_0 + (i + donor_decomp%ii_0 - 1) * donor_decomp%g_jj)))
    end do
  end do
  call preport_test("inverse_transpose_j: 1 x 4", fne(m, 0.0D0))
  
  deallocate(donor)
  deallocate(target)
  
  call deallocate(donor_decomp)
  
  call allocate(donor_decomp, iis = (/256, 257/), jjs = (/256, 257/), layers = 1)  
  allocate(donor(donor_decomp%ii + 2, donor_decomp%jj + 2, 1))
  
  target_decomp => donor_decomp%trans_i_decomp
  allocate(target(target_decomp%ii + 2, target_decomp%jj + 2, 1))
  
  donor = -1.0D0
  do j = 1, donor_decomp%jj
    do i = 1, donor_decomp%ii
      donor(i + 1, j + 1, 1) = real(j + donor_decomp%jj_0 + (i + donor_decomp%ii_0 - 1) * donor_decomp%g_jj, kind = real_8)
    end do
  end do
  target = 0.0D0
  
  call ptranspose(donor_decomp, donor, target, TRANSPOSE_II)
  call flush(donor_decomp, TRANSPOSE_II)
  
  m = 0.0D0
  do j = 1, target_decomp%jj
    do i = 1, target_decomp%ii
      m = max(m, dabs(target(i + 1, j + 1, 1) - (j + target_decomp%jj_0 + (i + target_decomp%ii_0 - 1) * target_decomp%g_jj)))
    end do
  end do
  call preport_test("transpose_i: 2 x 2", fne(m, 0.0D0))
  
  donor = 0.0D0
  
  call pinverse_transpose(donor_decomp, target, donor, TRANSPOSE_II)
  call flush(donor_decomp, TRANSPOSE_II)
  
  m = 0.0D0
  do j = 1, donor_decomp%jj
    do i = 1, donor_decomp%ii
      m = max(m, dabs(donor(i + 1, j + 1, 1) - (j + donor_decomp%jj_0 + (i + donor_decomp%ii_0 - 1) * donor_decomp%g_jj)))
    end do
  end do
  call preport_test("inverse_transpose_i: 2 x 2", fne(m, 0.0D0))
  
  deallocate(target)
  
  target_decomp => donor_decomp%trans_j_decomp
  allocate(target(target_decomp%ii + 2, target_decomp%jj + 2, 1))
  
  donor = -1.0D0
  do j = 1, donor_decomp%jj
    do i = 1, donor_decomp%ii
      donor(i + 1, j + 1, 1) = real(j + donor_decomp%jj_0 + (i + donor_decomp%ii_0 - 1) * donor_decomp%g_jj, kind = real_8)
    end do
  end do
  target = 0.0D0
  
  call ptranspose(donor_decomp, donor, target, TRANSPOSE_JJ)
  call flush(donor_decomp, TRANSPOSE_JJ)
  
  m = 0.0D0
  do j = 1, target_decomp%jj
    do i = 1, target_decomp%ii
      m = max(m, dabs(target(i + 1, j + 1, 1) - (j + target_decomp%jj_0 + (i + target_decomp%ii_0 - 1) * target_decomp%g_jj)))
    end do
  end do
  call preport_test("transpose_j: 2 x 2", fne(m, 0.0D0))
  
  donor = 0.0D0
  
  call pinverse_transpose(donor_decomp, target, donor, TRANSPOSE_JJ)
  call flush(donor_decomp, TRANSPOSE_JJ)
  
  m = 0.0D0
  do j = 1, donor_decomp%jj
    do i = 1, donor_decomp%ii
      m = max(m, dabs(donor(i + 1, j + 1, 1) - (j + donor_decomp%jj_0 + (i + donor_decomp%ii_0 - 1) * donor_decomp%g_jj)))
    end do
  end do
  call preport_test("inverse_transpose_j: 2 x 2", fne(m, 0.0D0))
  
  deallocate(donor)
  deallocate(target)
  
  call deallocate(donor_decomp)
  
  call parallel_finalise()
  
end program test_transpose_data_4
