#include "config.h"

program test_cabaret_cfl

#ifdef CABARET
  use qg_aux
  use qg_parallel
  use test_tools
  
  implicit none
  
  type(decomposition) :: decomp
  type(parallel_flux_data) :: u_i, u_j
  
  call parallel_init()
  
  call allocate(decomp, iis = (/65/), jjs = (/65/), layers = 1)
  call allocate(u_i, decomp, FLUX_II)
  call allocate(u_j, decomp, FLUX_JJ)
  
  call set(u_i, 1.0D0)
  call zero(u_j)
  call report_test("cabaret_cfl", fne(cabaret_cfl(1.0D0, u_i, u_j), 1.0D0))
  call report_test("cabaret_cfl", fne(cabaret_cfl(2.0D0, u_i, u_j), 2.0D0))
  
  call set(u_j, 2.0D0)
  call report_test("cabaret_cfl", fne(cabaret_cfl(1.0D0, u_i, u_j), 2.0D0))
  call report_test("cabaret_cfl", fne(cabaret_cfl(2.0D0, u_i, u_j), 4.0D0))

  u_i%data(10, 10, 1) = 10.0D0
  call report_test("cabaret_cfl", fne(cabaret_cfl(1.0D0, u_i, u_j), 10.0D0))
  call report_test("cabaret_cfl", fne(cabaret_cfl(2.0D0, u_i, u_j), 20.0D0))

  u_j%data(10, 10, 1) = 20.0D0
  call report_test("cabaret_cfl", fne(cabaret_cfl(1.0D0, u_i, u_j), 20.0D0))
  call report_test("cabaret_cfl", fne(cabaret_cfl(2.0D0, u_i, u_j), 40.0D0))
  
  call deallocate(u_i)
  call deallocate(u_j)
  call deallocate(decomp)

  call parallel_finalise()

#else
  implicit none
#endif

end program test_cabaret_cfl
