program test_inverse
  
  use lapack_interfaces
  use qg_parameters, only : real_8
  use test_tools
  
  implicit none
  
  real(kind = real_8), dimension(2, 2) :: matrix, matrix2
  
  matrix(1, 1) = 1.0D0;  matrix(1, 2) = 0.0D0
  matrix(2, 1) = 0.0D0;  matrix(2, 2) = 1.0D0
  call report_test("determinant: correct determinant", fne(determinant(matrix), 1.0D0))
  matrix2 = inverse(matrix)
  call report_test("inverse: correct inverse", fne(matrix2, reshape((/1.0D0, 0.0D0, 0.0D0, 1.0D0/), (/2, 2/))))
  
  matrix(1, 1) = 1.0D0;  matrix(1, 2) = 0.0D0
  matrix(2, 1) = 0.0D0;  matrix(2, 2) = 2.0D0
  call report_test("determinant: correct determinant", fne(determinant(matrix), 2.0D0))
  matrix2 = inverse(matrix)
  call report_test("inverse: correct inverse", fne(matrix2, reshape((/1.0D0, 0.0D0, 0.0D0, 0.5D0/), (/2, 2/))))
  
  matrix(1, 1) = 1.0D0;  matrix(1, 2) = 1.0D0
  matrix(2, 1) = 0.0D0;  matrix(2, 2) = 1.0D0
  call report_test("determinant: correct determinant", fne(determinant(matrix), 1.0D0))
  matrix2 = inverse(matrix)
  call report_test("inverse: correct inverse", fne(matrix2, reshape((/1.0D0, 0.0D0, -1.0D0, 1.0D0/), (/2, 2/))))
  
  matrix(1, 1) = 1.0D0;  matrix(1, 2) = -4.0D0
  matrix(2, 1) = 1.0D0;  matrix(2, 2) = 1.0D0
  call report_test("determinant: correct determinant", fne(determinant(matrix), 5.0D0))
  matrix2 = inverse(matrix)
  call report_test("inverse: correct inverse", fne(matrix2, reshape((/0.2D0, -0.2D0, 0.8D0, 0.2D0/), (/2, 2/))))

end program test_inverse
