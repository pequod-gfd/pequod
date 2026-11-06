program test_eigendecompose
  
  use lapack_interfaces
  use qg_diagnostics
  use qg_parallel
  use qg_parameters, only : complex_8, real_8
  use test_tools
  
  implicit none
  
  real(kind = real_8), dimension(2, 2) :: matrix, matrix2
  real(kind = real_8), dimension(2) :: eigenvalues_r
  complex(kind = complex_8), dimension(2) :: eigenvalues_c
  real(kind = real_8), dimension(2, 2) :: eigenvectors
  
  type(stratification) :: strat
  real(kind = real_8), parameter :: rd_bcl1 = 20.0D5, rd_bcl2 = 10.0D5
  real(kind = real_8), dimension(3), parameter :: h = (/1.0D5, 2.0D5, 3.0D5/)
  real(kind = real_8), dimension(3, 3) :: matrix_3
  real(kind = real_8), dimension(3) :: eigenvalues_3
  real(kind = real_8), dimension(3, 3) :: eigenvectors_3
  
  call parallel_init()
  
  matrix(1, 1) = 1.0D0;  matrix(1, 2) = 0.0D0
  matrix(2, 1) = 0.0D0;  matrix(2, 2) = 1.0D0
  call eigendecompose(matrix, eigenvalues_r, eigenvectors)
  call report_test("eigendecompose_real: correct eigenvalues", fne(eigenvalues_r, (/1.0D0, 1.0D0/)))
  
  matrix(1, 1) = 1.0D0;  matrix(1, 2) = 0.0D0
  matrix(2, 1) = 0.0D0;  matrix(2, 2) = 2.0D0
  call eigendecompose(matrix, eigenvalues_r, eigenvectors)
  call report_test("eigendecompose_real: correct eigenvalues", fne(eigenvalues_r, (/1.0D0, 2.0D0/)))
  
  matrix(1, 1) = 1.0D0;  matrix(1, 2) = 0.0D0
  matrix(2, 1) = 0.0D0;  matrix(2, 2) = 1.0D0
  call eigendecompose(matrix, eigenvalues_r, eigenvectors)
  call eigenrecompose(matrix2, eigenvalues_r, eigenvectors)
  call report_test("eigerecompose_real: correct recomposition", fne(matrix, matrix2))  
  
  matrix(1, 1) = 4.0D0;  matrix(1, 2) = 5.0D0
  matrix(2, 1) = 1.0D0;  matrix(2, 2) = 18.0D0
  call eigendecompose(matrix, eigenvalues_r, eigenvectors)
  call eigenrecompose(matrix2, eigenvalues_r, eigenvectors)
  call report_test("eigerecompose_real: correct recomposition", fne(matrix, matrix2))  
  
  matrix(1, 1) = 1.0D0;  matrix(1, 2) = 0.0D0
  matrix(2, 1) = 0.0D0;  matrix(2, 2) = 1.0D0
  call eigendecompose(matrix, eigenvalues_c, eigenvectors)
  call report_test("eigendecompose_complex: correct eigenvalues", fne(eigenvalues_c, (/cmplx(1.0D0, 0.0D0, kind = complex_8), cmplx(1.0D0, 0.0D0, kind = complex_8)/)))
  
  matrix(1, 1) = 1.0D0;  matrix(1, 2) = -4.0D0
  matrix(2, 1) = 1.0D0;  matrix(2, 2) = 1.0D0
  call eigendecompose(matrix, eigenvalues_c, eigenvectors)
  call report_test("eigendecompose_complex: correct eigenvalues", fne(eigenvalues_c, (/cmplx(1.0D0, 2.0D0, kind = complex_8), cmplx(1.0D0, -2.0D0, kind = complex_8)/)))
  
  matrix(1, 1) = 1.0D0;  matrix(1, 2) = 0.0D0
  matrix(2, 1) = 0.0D0;  matrix(2, 2) = 1.0D0
  call eigendecompose_symmetric(matrix, eigenvalues_r, eigenvectors)
  call report_test("eigendecompose_symmetric: correct eigenvalues", fne(eigenvalues_r, (/1.0D0, 1.0D0/)))
  
  matrix(1, 1) = 1.0D0;  matrix(1, 2) = 0.0D0
  matrix(2, 1) = 0.0D0;  matrix(2, 2) = 2.0D0
  call eigendecompose_symmetric(matrix, eigenvalues_r, eigenvectors)
  call report_test("eigendecompose_symmetric: correct eigenvalues", fne(eigenvalues_r, (/1.0D0, 2.0D0/)))
  
  matrix(1, 1) = 1.0D0;  matrix(1, 2) = 0.0D0
  matrix(2, 1) = 0.0D0;  matrix(2, 2) = 1.0D0
  call eigendecompose_symmetric(matrix, eigenvalues_r, eigenvectors)
  call eigenrecompose_symmetric(matrix2, eigenvalues_r, eigenvectors)
  call report_test("eigerecompose_symmetric: correct recomposition", fne(matrix, matrix2))  
  
  matrix(1, 1) = 4.0D0;  matrix(1, 2) = 5.0D0
  matrix(2, 1) = 5.0D0;  matrix(2, 2) = 18.0D0
  call eigendecompose_symmetric(matrix, eigenvalues_r, eigenvectors)
  call eigenrecompose_symmetric(matrix2, eigenvalues_r, eigenvectors)
  call report_test("eigerecompose_symmetric: correct recomposition", fne(matrix, matrix2)) 
  
  call allocate(strat, rd_bcl1, rd_bcl2, h(1), h(2), h(3))
  call report_test("eigendecompose_symmetric: correct eigenvalues", fne((1.0D0 / dsqrt(strat%ss(2:3))), (/Rd_bcl2, Rd_bcl1/)))
  
  call eigendecompose(strat%s_mat, eigenvalues_3, eigenvectors_3)
  call report_test("eigendecompose_symmetric: correct eigenvalues", fne((1.0D0 / dsqrt(-eigenvalues_3(1:2))), (/Rd_bcl2, Rd_bcl1/)))
    
  matrix_3(1, 1) = -strat%s2(1);  matrix_3(1, 2) = strat%s2(1);                   matrix_3(1, 3) = 0.0D0
  matrix_3(2, 1) = strat%s1(2);   matrix_3(2, 2) = -(strat%s1(2) + strat%s2(2));  matrix_3(2, 3) = strat%s1(3)
  matrix_3(3, 1) = 0.0D0;         matrix_3(3, 2) = strat%s2(2);                   matrix_3(3, 3) = -strat%s1(3)
  call report_test("eigendecompose_symmetric: correct eigenvalues", fne((1.0D0 / dsqrt(-eigenvalues_3(1:2))), (/Rd_bcl2, Rd_bcl1/)))
  
  call deallocate(strat)
  
  call parallel_finalise()
     
end program test_eigendecompose
