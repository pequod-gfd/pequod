! Monomial averages in grid units. Use descending x power within each degree.
module qg_moments
  use qg_parameters, only: real_8
  implicit none
  private
  public :: n_moments, x_power, y_power, monomial_index, point_moments
  public :: triangle_moments, segment_moments, translate_moments

  integer, parameter :: n_moments = 15
  integer, parameter :: x_power(n_moments) = [0,1,0,2,1,0,3,2,1,0,4,3,2,1,0]
  integer, parameter :: y_power(n_moments) = [0,0,1,0,1,2,0,1,2,3,0,1,2,3,4]
  integer, parameter :: factorial(0:6) = [1,1,2,6,24,120,720]

contains

  pure integer function monomial_index(a,b) result(index)
    integer, intent(in) :: a,b
    index = (a+b)*(a+b+1)/2+b+1
  end function monomial_index

  pure integer function binomial(a,b) result(coefficient)
    integer, intent(in) :: a,b
    coefficient = factorial(a)/(factorial(b)*factorial(a-b))
  end function binomial

  pure function point_moments(x) result(row)
    real(real_8), intent(in) :: x(2)
    real(real_8) :: row(n_moments)
    row = x(1)**x_power * x(2)**y_power
  end function point_moments

  ! Triangle vertices (0,p,q). Barycentric u,v have mean
  ! u^r v^s = 2*r!*s!/(r+s+2)! independently of triangle size.
  pure function triangle_moments(p,q) result(row)
    real(real_8), intent(in) :: p(2),q(2)
    real(real_8) :: row(n_moments), weight
    integer :: m,a,b,r,s,k,l

    row = 0.0D0
    do m = 1,n_moments
      a = x_power(m)
      b = y_power(m)
      do r = 0,a
        do s = 0,b
          k = r+s
          l = a+b-k
          weight = real(2*factorial(k)*factorial(l),real_8)/factorial(a+b+2)
          row(m) = row(m)+weight*binomial(a,r)*binomial(b,s)*p(1)**r*q(1)**(a-r)*p(2)**s*q(2)**(b-s)
        end do
      end do
    end do
  end function triangle_moments

  ! Segment mean with the physical/native length supplied separately.
  pure function segment_moments(p,q) result(row)
    real(real_8), intent(in) :: p(2),q(2)
    real(real_8) :: row(n_moments), midpoint(2), delta(2)
    midpoint = 0.5D0*(p+q)
    delta = 0.5D0*sqrt(0.6D0)*(q-p)
    row = (5*point_moments(midpoint-delta)+8*point_moments(midpoint)+5*point_moments(midpoint+delta))/18
  end function segment_moments

  ! Evaluate moments of (x+offset_x, y+offset_y).
  pure function translate_moments(local, offset) result(row)
    real(real_8), intent(in) :: local(n_moments), offset(2)
    real(real_8) :: row(n_moments)
    integer :: m,a,b,r,s

    row = 0.0D0
    do m = 1,n_moments
      a = x_power(m)
      b = y_power(m)
      do r = 0,a
        do s = 0,b
          row(m) = row(m)+binomial(a,r)*binomial(b,s)*offset(1)**(a-r)*offset(2)**(b-s)*local(monomial_index(r,s))
        end do
      end do
    end do
  end function translate_moments

end module qg_moments
