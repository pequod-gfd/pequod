program test_irregular_moments
  use qg_irregular_parallel, only: make_sector_moments, cell_mean_moments
  use qg_moments
  use qg_parameters, only: real_8
  implicit none
  real(real_8), parameter :: gauss(4) = [0.06943184420297371D0,0.33000947820757187D0, &
                                         0.66999052179242813D0,0.93056815579702629D0]
  real(real_8), parameter :: weight(4) = [0.17392742256872693D0,0.32607257743127307D0, &
                                          0.32607257743127307D0,0.17392742256872693D0]
  real(real_8), parameter :: boole(0:4) = [7,32,12,32,7]/90.0D0
  real(real_8), parameter :: ray(2,0:7) = 0.5D0*reshape([1,0, 1,1, 0,1, -1,1, -1,0, -1,-1, 0,-1, 1,-1],[2,8])
  real(real_8), parameter :: origin(2) = [0.375D0,-0.75D0]
  real(real_8) :: sector(n_moments,0:7), quadrature(n_moments,0:7), predicted(n_moments), expected(n_moments)
  real(real_8) :: x(2), p(2), q(2), defect
  integer :: r,s,i,j,m,c,mask

  ! Independent Duffy transform of a square, with tensor Gauss weights
  ! including the Jacobian. Normalization gives mean rather than integral.
  sector = make_sector_moments()
  defect = 0.0D0
  do r = 0,7
    s = modulo(r+1,8)
    p = ray(:,r)
    q = ray(:,s)
    quadrature(:,r) = 0.0D0
    expected = 0.0D0
    do i = 1,4
      do j = 1,4
        x = gauss(i)*(p+gauss(j)*(q-p))+origin
        quadrature(:,r) = quadrature(:,r)+2*gauss(i)*weight(i)*weight(j)*independent_monomials(x)
      end do
    end do
    do i = 0,4
      expected = expected+boole(i)*independent_monomials(p+real(i,real_8)/4*(q-p)+origin)
    end do
    defect = max(defect,maxval(abs(expected-segment_moments(p+origin,q+origin))))
  end do
  do c = 1,8
    do r = 0,7
      mask = 0
      expected = 0.0D0
      do i = 0,c-1
        s = modulo(r+i,8)
        mask = ibset(mask,s)
        expected = expected+quadrature(:,s)/c
      end do
      predicted = translate_moments(cell_mean_moments(mask,sector),origin)
      defect = max(defect,maxval(abs(predicted-expected)))
    end do
  end do
  if (defect > 4e-15) stop 1
  print '(a,es14.6,1x,a)', 'All rotated degree-four cell and segment moments, max defect:', defect, 'PASS'

contains
  pure function independent_monomials(x) result(row)
    real(real_8), intent(in) :: x(2)
    real(real_8) :: row(15)
    integer :: d,b,k

    k = 0
    do d = 0,4
      do b = 0,d
        k = k+1
        row(k) = x(1)**(d-b)*x(2)**b
      end do
    end do
  end function independent_monomials
end program test_irregular_moments
