!!NP 4

program test_pdata_4

  use qg_parallel
  use test_tools

  implicit none
  
  integer :: n, stat
  type(decomposition) :: decomp
  
  call parallel_init()
  
  n = nprocs()
  call preport_test("nprocs == 4", n /= 4)
  
  call allocate(decomp, iis = (/5/), jjs = (/4, 4, 4, 5/), layers = 1)
  
  call test_pdata(decomp, stat = stat)
  call preport_test("test_pdata: 1 x 4", stat /= 0)
  
  call deallocate(decomp)
  
  call allocate(decomp, iis = (/4, 4, 4, 5/), jjs = (/5/), layers = 1)
  
  call test_pdata(decomp, stat = stat)
  call preport_test("test_pdata: 4 x 1", stat /= 0)
  
  call deallocate(decomp)
  
  call allocate(decomp, iis = (/4, 5/), jjs = (/4, 5/), layers = 1)
  
  call test_pdata(decomp, stat = stat)
  call preport_test("test_pdata: 2 x 2", stat /= 0)
  
  call deallocate(decomp)
  
  call parallel_finalise()

end program test_pdata_4
