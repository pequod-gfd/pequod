program test_solve

  use lapack_interfaces
  use qg_parameters, only : real_8
  use test_tools
  
  implicit none
  
  real(kind = real_8), dimension(3, 3) :: matrix, matrix_inverse
  real(kind = real_8), dimension(3) :: rhs, soln1, soln2

  matrix(1, 1) = 1.0D0;  matrix(1, 2) = 0.0D0;  matrix(1, 3) = 0.0D0
  matrix(2, 1) = 0.0D0;  matrix(2, 2) = 1.0D0;  matrix(2, 3) = 0.0D0
  matrix(3, 1) = 0.0D0;  matrix(3, 2) = 0.0D0;  matrix(3, 3) = 1.0D0
  
  rhs(1) = 2.0D0;  rhs(2) = 3.0D0;  rhs(3) = 4.0D0
  call solve(matrix, rhs, soln1)
  
  call report_test("solve", fne(soln1, rhs))

  matrix(1, 1) = 1.0D0;  matrix(1, 2) = 2.0D0;   matrix(1, 3) = 3.0D0
  matrix(2, 1) = 0.5D0;  matrix(2, 2) = 0.25D0;  matrix(2, 3) = 0.0D0
  matrix(3, 1) = 4.1D0;  matrix(3, 2) = 0.1D0;   matrix(3, 3) = 4.0D0
  
  matrix_inverse = inverse(matrix)
  soln1 = matmul(matrix_inverse, rhs)
  
  call solve(matrix, rhs, soln2)
  
  call report_test("solve", fne(soln1, soln2))

end program test_solve
