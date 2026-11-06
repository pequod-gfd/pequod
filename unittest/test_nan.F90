#include "config.h"

program test_nan

#ifdef DDEBUG
  use qg_parallel
  use qg_parameters, only : real_8
  use test_tools
  
  implicit none
  
  real(kind = real_8) :: r
  
  r = nan()
  call report_test("nan == nan", r == r)
  call report_test("nan > nan", r > r)
  call report_test("nan < nan", r < r)
#else
  implicit none
#endif
  
end program test_nan
