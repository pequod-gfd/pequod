!!NP 2

#include "config.h"

program test_decomposition_2
  
  use qg_parallel
  use test_tools
 
  implicit none
  
  integer :: n, p
  logical :: fail
  type(decomposition) :: decomp
  
  call parallel_init()
  
  n = nprocs()
  p = procno()
  call preport_test("nprocs == 2", n /= 2)
  
  if(p == 1) write(6, "(a)") "Decomposed in j"
  call allocate(decomp, iis = (/5/), jjs = (/8, 9/), layers = 1)
  
  if(decomp%p_j == 1) then
    fail = decomp%p_a /= mod(decomp%p + 1, 2)
  else
    fail = decomp%p_c /= mod(decomp%p + 1, 2)
  end if
  call preport_test("decomposition: neighbours", fail) 
  
  fail = decomp%layers /= 1
  call preport_test("decomposition: layers", fail) 
  
  fail = decomp%ii /= 5
  call preport_test("decomposition: ii", fail)
  
  if(decomp%p_j == 1) then
    fail = decomp%jj /= 8
  else
    fail = decomp%jj /= 9
  end if
  call preport_test("decomposition: jj", fail)
  
  fail = decomp%ii2 /= 3
  call preport_test("decomposition: ii2", fail)
  
  if(decomp%p_j == 1) then
    fail = decomp%jj2 /= 7
  else
    fail = decomp%jj2 /= 8
  end if
  call preport_test("decomposition: jj2", fail)
  
  call preport_test("decomposition: ii2s", any(decomp%ii2s /= (/3/)))
  
  call preport_test("decomposition: jj2s", any(decomp%jj2s /= (/7, 8/)))
  
  fail = decomp%ii_0 /= 0
  call preport_test("decomposition: ii_0", fail)
  
  if(decomp%p_j == 1) then
    fail = decomp%jj_0 /= 0
  else
    fail = decomp%jj_0 /= 8
  end if
  call preport_test("decomposition: jj_0", fail)
  
  fail = decomp%ii2_0 /= 0
  call preport_test("decomposition: ii2_0", fail)
  
  if(decomp%p_j == 1) then
    fail = decomp%jj2_0 /= 0
  else
    fail = decomp%jj2_0 /= 7
  end if
  call preport_test("decomposition: jj2_0", fail)
  
  fail = decomp%fa_ii /= 5
  call preport_test("decomposition: fa_ii", fail)
  
  if(decomp%p_j == 1) then
    fail = decomp%fa_jj /= 9
  else
    fail = decomp%fa_jj /= 10
  end if
  call preport_test("decomposition: fa_jj", fail)
  
  fail = decomp%g_ii /= 5
  call preport_test("decomposition: g_ii", fail)
  
  fail = decomp%g_jj /= 17
  call preport_test("decomposition: g_jj", fail)
  
  fail = decomp%g_ii2 /= 3
  call preport_test("decomposition: g_ii2", fail)
  
  fail = decomp%g_jj2 /= 15
  call preport_test("decomposition: g_jj2", fail)
  
  fail = decomp%g_ii_int /= 4
  call preport_test("decomposition: g_ii_int", fail)
  
  fail = decomp%g_jj_int /= 16
  call preport_test("decomposition: g_jj_int", fail)
  
  fail = decomp%l_i_min /= 2
  call preport_test("decomposition: l_i_min", fail)
  
  if(decomp%p_j == 1) then
    fail = decomp%l_j_min /= 2
  else
    fail = decomp%l_j_min /= 1
  end if
  call preport_test("decomposition: l_j_min", fail)
  
  fail = decomp%l_i_max /= 4
  call preport_test("decomposition: l_i_max", fail)
  
  fail = decomp%l_j_max /= 8
  call preport_test("decomposition: l_j_max", fail)
  
  fail = decomp%f_i_min /= 3
  call preport_test("decomposition: f_i_min", fail)
  
  if(decomp%p_j == 1) then
    fail = decomp%f_j_min /= 3
  else
    fail = decomp%f_j_min /= 2
  end if
  call preport_test("decomposition: f_j_min", fail)
  
  fail = decomp%f_i_max /= 5
  call preport_test("decomposition: f_i_max", fail)
  
  fail = decomp%f_j_max /= 9
  call preport_test("decomposition: f_j_max", fail)
  
  fail = decomp%fa_i_min /= 2
  call preport_test("decomposition: fa_i_min", fail)
  
  if(decomp%p_j == 1) then
    fail = decomp%fa_j_min /= 2
  else
    fail = decomp%fa_j_min /= 1
  end if
  call preport_test("decomposition: fa_j_min", fail)
  
  fail = decomp%fa_i_max /= 6
  call preport_test("decomposition: fa_i_max", fail)
  
  fail = decomp%fa_j_max /= 10
  call preport_test("decomposition: fa_j_max", fail)
  
  fail = decomp%fb_i_min /= 2
  call preport_test("decomposition: fb_i_min", fail)
  
  fail = decomp%fb_j_min /= 2
  call preport_test("decomposition: fb_j_min", fail)
  
  fail = decomp%fb_i_max /= 6
  call preport_test("decomposition: fb_i_max", fail)
  
  if(decomp%p_j == 1) then
    fail = decomp%fb_j_max /= 9
  else
    fail = decomp%fb_j_max /= 10
  end if
  call preport_test("decomposition: fb_j_max", fail)
  
  call deallocate(decomp)
  
  if(p == 1) write(6, "(a)") "Decomposed in i"
  call allocate(decomp, iis = (/4, 5/), jjs = (/9/), layers = 1)
  
  if(decomp%p_i == 1) then
    fail = decomp%p_d /= mod(decomp%p + 1, 2)
  else
    fail = decomp%p_b /= mod(decomp%p + 1, 2)
  end if
  call preport_test("decomposition: neighbours", fail)  
  
  fail = decomp%layers /= 1
  call preport_test("decomposition: layers", fail)
  
  if(decomp%p_i == 1) then
    fail = decomp%ii /= 4
  else
    fail = decomp%ii /= 5
  end if
  call preport_test("decomposition: ii", fail)
  
  fail = decomp%jj /= 9
  call preport_test("decomposition: jj", fail)
  
  if(decomp%p_i == 1) then
    fail = decomp%ii2 /= 3
  else
    fail = decomp%ii2 /= 4
  end if
  call preport_test("decomposition: ii2", fail)
  
  fail = decomp%jj2 /= 7
  call preport_test("decomposition: jj2", fail)
  
  call preport_test("decomposition: ii2s", any(decomp%ii2s /= (/3, 4/)))
  
  call preport_test("decomposition: jj2s", any(decomp%jj2s /= (/7/)))

  if(decomp%p_i == 1) then
    fail = decomp%ii_0 /= 0
  else
    fail = decomp%ii_0 /= 4
  end if
  call preport_test("decomposition: ii_0", fail)
  
  fail = decomp%jj_0 /= 0
  call preport_test("decomposition: jj_0", fail)
  
  if(decomp%p_i == 1) then
    fail = decomp%ii2_0 /= 0
  else
    fail = decomp%ii2_0 /= 3
  end if
  call preport_test("decomposition: ii2_0", fail)
  
  fail = decomp%jj2_0 /= 0
  call preport_test("decomposition: jj2_0", fail)
  
  if(decomp%p_i == 1) then
    fail = decomp%fa_ii /= 5
  else
    fail = decomp%fa_ii /= 6
  end if
  call preport_test("decomposition: fa_ii", fail)
  
  fail = decomp%fa_jj /= 9
  call preport_test("decomposition: fa_jj", fail)
  
  fail = decomp%g_ii /= 9
  call preport_test("decomposition: g_ii", fail)
  
  fail = decomp%g_jj /= 9
  call preport_test("decomposition: g_jj", fail)
  
  fail = decomp%g_ii2 /= 7
  call preport_test("decomposition: g_ii2", fail)
  
  fail = decomp%g_jj2 /= 7
  call preport_test("decomposition: g_jj2", fail)
  
  fail = decomp%g_ii_int /= 8
  call preport_test("decomposition: g_ii_int", fail)
  
  fail = decomp%g_jj_int /= 8
  call preport_test("decomposition: g_jj_int", fail)
  
  if(decomp%p_i == 1) then
    fail = decomp%l_i_min /= 2
  else
    fail = decomp%l_i_min /= 1
  end if
  call preport_test("decomposition: l_i_min", fail)
  
  fail = decomp%l_j_min /= 2
  call preport_test("decomposition: l_j_min", fail)
  
  fail = decomp%l_i_max /= 4
  call preport_test("decomposition: l_i_max", fail)
  
  fail = decomp%l_j_max /= 8
  call preport_test("decomposition: l_j_max", fail)
  
  if(decomp%p_i == 1) then
    fail = decomp%f_i_min /= 3
  else
    fail = decomp%f_i_min /= 2
  end if
  call preport_test("decomposition: f_i_min", fail)
  
  fail = decomp%f_j_min /= 3
  call preport_test("decomposition: f_j_min", fail)
  
  fail = decomp%f_i_max /= 5
  call preport_test("decomposition: f_i_max", fail)
  
  fail = decomp%f_j_max /= 9
  call preport_test("decomposition: f_j_max", fail)
  
  if(decomp%p_i == 1) then
    fail = decomp%fa_i_min /= 2
  else
    fail = decomp%fa_i_min /= 1
  end if
  call preport_test("decomposition: fa_i_min", fail)
  
  fail = decomp%fa_j_min /= 2
  call preport_test("decomposition: fa_j_min", fail)
  
  fail = decomp%fa_i_max /= 6
  call preport_test("decomposition: fa_i_max", fail)
  
  fail = decomp%fa_j_max /= 10
  call preport_test("decomposition: fa_j_max", fail)
  
  fail = decomp%fb_i_min /= 2
  call preport_test("decomposition: fb_i_min", fail)
  
  fail = decomp%fb_j_min /= 2
  call preport_test("decomposition: fb_j_min", fail)
  
  if(decomp%p_i == 1) then
    fail = decomp%fb_i_max /= 5
  else
    fail = decomp%fb_i_max /= 6
  end if
  call preport_test("decomposition: fb_i_max", fail)
  
  fail = decomp%fb_j_max /= 10
  call preport_test("decomposition: fb_j_max", fail)
  
  call deallocate(decomp)
  
  if(p == 1) write(6, "(a)") "Decomposed in j, periodic in j"
  call allocate(decomp, iis = (/5/), jjs = (/8, 8/), layers = 1, periodic_j = .true.)
  
  fail = decomp%p_a /= mod(decomp%p + 1, 2) .and. &
       & decomp%p_c /= mod(decomp%p + 1, 2)
  call preport_test("decomposition: neighbours", fail)  
  
  fail = decomp%layers /= 1
  call preport_test("decomposition: layers", fail)
  
  fail = decomp%ii /= 5
  call preport_test("decomposition: ii", fail)
  
  fail = decomp%jj /= 8
  call preport_test("decomposition: jj", fail)
  
  fail = decomp%ii2 /= 3
  call preport_test("decomposition: ii2", fail)
  
  fail = decomp%jj2 /= 8
  call preport_test("decomposition: jj2", fail)
  
  call preport_test("decomposition: ii2s", any(decomp%ii2s /= (/3/)))
  
  call preport_test("decomposition: jj2s", any(decomp%jj2s /= (/8, 8/)))
  
  fail = decomp%ii_0 /= 0
  call preport_test("decomposition: ii_0", fail)
  
  if(decomp%p_j == 1) then
    fail = decomp%jj_0 /= 0
  else
    fail = decomp%jj_0 /= 8
  end if
  call preport_test("decomposition: jj_0", fail)
  
  fail = decomp%ii2_0 /= 0
  call preport_test("decomposition: ii2_0", fail)
  
  if(decomp%p_j == 1) then
    fail = decomp%jj2_0 /= 0
  else
    fail = decomp%jj2_0 /= 8
  end if
  call preport_test("decomposition: jj2_0", fail)
  
  fail = decomp%fa_ii /= 5
  call preport_test("decomposition: fa_ii", fail)

  fail = decomp%fa_jj /= 10
  call preport_test("decomposition: fa_jj", fail)
  
  fail = decomp%g_ii /= 5
  call preport_test("decomposition: g_ii", fail)
  
  fail = decomp%g_jj /= 16
  call preport_test("decomposition: g_jj", fail)
  
  fail = decomp%g_ii2 /= 3
  call preport_test("decomposition: g_ii2", fail)
  
  fail = decomp%g_jj2 /= 16
  call preport_test("decomposition: g_jj2", fail)
  
  fail = decomp%g_ii_int /= 4
  call preport_test("decomposition: g_ii_int", fail)
  
  fail = decomp%g_jj_int /= 16
  call preport_test("decomposition: g_jj_int", fail)
  
  fail = decomp%l_i_min /= 2
  call preport_test("decomposition: l_i_min", fail)
  
  fail = decomp%l_j_min /= 1
  call preport_test("decomposition: l_j_min", fail)
  
  fail = decomp%l_i_max /= 4
  call preport_test("decomposition: l_i_max", fail)
  
  fail = decomp%l_j_max /= 8
  call preport_test("decomposition: l_j_max", fail)
  
  fail = decomp%f_i_min /= 3
  call preport_test("decomposition: f_i_min", fail)
  
  fail = decomp%f_j_min /= 2
  call preport_test("decomposition: f_j_min", fail)
  
  fail = decomp%f_i_max /= 5
  call preport_test("decomposition: f_i_max", fail)
  
  fail = decomp%f_j_max /= 9
  call preport_test("decomposition: f_j_max", fail)
  
  fail = decomp%fa_i_min /= 2
  call preport_test("decomposition: fa_i_min", fail)
  
  fail = decomp%fa_j_min /= 1
  call preport_test("decomposition: fa_j_min", fail)
  
  fail = decomp%fa_i_max /= 6
  call preport_test("decomposition: fa_i_max", fail)
  
  fail = decomp%fa_j_max /= 10
  call preport_test("decomposition: fa_j_max", fail)
  
  fail = decomp%fb_i_min /= 2
  call preport_test("decomposition: fb_i_min", fail)
  
  fail = decomp%fb_j_min /= 2
  call preport_test("decomposition: fb_j_min", fail)
  
  fail = decomp%fb_i_max /= 6
  call preport_test("decomposition: fb_i_max", fail)
  
  fail = decomp%fb_j_max /= 9
  call preport_test("decomposition: fb_j_max", fail)
  
  call deallocate(decomp)
    
  call parallel_finalise()
  
end program test_decomposition_2
