! Compatible cut-cell curl and component torque from a parent momentum force.
module qg_irregular_forcing
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use qg_irregular_parallel, only: distributed_geometry, half_coordinates, share_face_flux, broadcast_error
  use qg_irregular_operators, only: sum_open_flux
  use qg_parameters, only: real_8, pi
  implicit none
  include 'mpif.h'
  private
  public :: wind_curl_torque, bottom_drag

  real(real_8), parameter :: gauss(4) = [0.1834346424956498D0,0.5255324099163290D0, &
                                        0.7966664774136267D0,0.9602898564975363D0]
  real(real_8), parameter :: weight(4) = [0.3626837833783620D0,0.3137066458778873D0, &
                                         0.2223810344533745D0,0.1012285362903763D0]

contains

  ! The legacy piecewise-sinusoidal curl is the negative y derivative
  ! of this continuous zonal acceleration. The coordinates are native.
  ! strength is the layer acceleration used in the legacy PV curl.
  subroutine wind_curl_torque(geom, strength, shift, tilt, asymmetry, curl, torque, error)
    type(distributed_geometry), intent(in) :: geom
    real(real_8), intent(in) :: strength,shift,tilt,asymmetry
    real(real_8), intent(out) :: curl(:,:,:),torque(0:,:)
    character(*), intent(out) :: error
    real(real_8), allocatable :: open(:,:)
    real(real_8) :: local(0:ubound(torque,1),size(torque,2)),params(4),p(2),q(2),r(2)
    real(real_8) :: xscale,yscale,d
    integer :: f,k,i,j,g_i,edge,previous,comp,ierr,nl,counts(2)

    error = ''
    nl = size(curl,3)
    if (any(shape(curl) /= [geom%ni,geom%nj,nl]) .or. &
        any(shape(torque) /= [size(geom%global%loops),nl]) .or. nl < 1) error = 'incompatible wind arrays'
    if (.not. all(ieee_is_finite([strength,shift,tilt,asymmetry])) .or. &
        abs(shift)+0.5D0*abs(tilt) >= 0.5D0 .or. asymmetry <= 0.0D0) error = 'invalid piecewise wind profile'
    call broadcast_error(geom%comm,geom%rank,geom%ranks,error)
    if (len_trim(error) /= 0) return
    call mpi_allreduce([nl,-nl],counts,2,MPI_INTEGER,MPI_MAX,geom%comm,ierr)
    if (counts(1) /= -counts(2)) then
      error = 'wind layer counts differ across ranks'
      return
    end if
    params = [strength,shift,tilt,asymmetry]
    call mpi_bcast(params,size(params),MPI_DOUBLE_PRECISION,0,geom%comm,ierr)
    if (any(params /= [strength,shift,tilt,asymmetry])) error = 'wind parameters differ across ranks'
    call broadcast_error(geom%comm,geom%rank,geom%ranks,error)
    if (len_trim(error) /= 0) return

    g_i = size(geom%owner,1)
    xscale = real(size(geom%owner,2)-1,real_8)
    yscale = real(g_i-1,real_8)
    allocate(open(1,size(geom%faces)))
    open = 0.0D0
    do f = 1,size(geom%faces)
      if (geom%faces(f)%owner /= geom%rank) cycle
      p = 0.5D0*geom%faces(f)%endpoints(:,1)
      q = 0.5D0*geom%faces(f)%endpoints(:,2)
      open(1,f) = wind_line(p,q)
    end do
    call share_face_flux(geom,open,error)
    if (len_trim(error) /= 0) return
    call sum_open_flux(geom,open,curl(:,:,1:1))
    curl(:,:,2:) = 0.0D0

    local = 0.0D0
    do k = 1,size(geom%walls)
      comp = geom%walls(k)%component
      edge = geom%walls(k)%edge
      associate(loop => geom%global%loops(comp)%ij)
        previous = modulo(edge-2,size(loop,2))+1
        r = 0.5D0*half_coordinates(loop(1,previous),loop(2,previous),g_i)
        i = loop(1,edge)-geom%i0
        j = loop(2,edge)-geom%j0
      end associate
      p = 0.5D0*geom%walls(k)%endpoints(:,1)
      q = 0.5D0*geom%walls(k)%endpoints(:,2)
      d = wind_line(0.5D0*(r+p),p)+wind_line(p,0.5D0*(p+q))
      curl(i,j,1) = curl(i,j,1)+d
      local(comp,1) = local(comp,1)+d
    end do
    call mpi_allreduce(local,torque,size(local),MPI_DOUBLE_PRECISION,MPI_SUM,geom%comm,ierr)

  contains

    function wind_line(first,second) result(integral)
      real(real_8), intent(in) :: first(2),second(2)
      real(real_8) :: integral,delta(2),line(0:2),a,b,p(2),x,y,mid,width,value,rates(2)
      integer :: n,j,k,s

      integral = 0.0D0
      delta = second-first
      if (delta(1) == 0.0D0) return
      rates = [first(2),second(2)]/yscale-0.5D0-shift-tilt*([first(1),second(1)]/xscale-0.5D0)
      line = [0.0D0,1.0D0,1.0D0]
      n = 1
      if (minval(rates) < 0.0D0 .and. maxval(rates) > 0.0D0) then
        line(1) = rates(1)/(rates(1)-rates(2))
        n = 2
      end if
      do j = 1,n
        a = 0.5D0*(line(j-1)+line(j))
        b = 0.5D0*(line(j)-line(j-1))
        do k = 1,size(gauss)
          do s = -1,1,2
            p = first+(a+s*b*gauss(k))*delta
            x = p(1)/xscale-0.5D0
            y = p(2)/yscale-0.5D0
            mid = shift+tilt*x
            if (y < mid) then
              width = 0.5D0+mid
              value = -asymmetry*width*(cos(pi*(0.5D0+y)/width)+1.0D0)
            else
              width = 0.5D0-mid
              value = width*(cos(pi*(y-mid)/width)-1.0D0)/asymmetry
            end if
            integral = integral+b*weight(k)*value
          end do
        end do
      end do
      integral = integral*delta(1)*2.0D0*yscale*strength
    end function wind_line

  end subroutine wind_curl_torque

  ! Add linear bottom momentum drag and its integrated curl in place.
  ! omega and gamma must be from the same compatible Poisson stage.
  subroutine bottom_drag(geom, coefficient, omega, gamma, curl, torque)
    type(distributed_geometry), intent(in) :: geom
    real(real_8), intent(in) :: coefficient,omega(:,:,:),gamma(0:,:)
    real(real_8), intent(inout) :: curl(:,:,:),torque(0:,:)
    integer :: i,j,bottom

    bottom = size(omega,3)
    do j = 1,geom%nj
      do i = 1,geom%ni
        if (geom%area(i,j) == 0.0D0) cycle
        curl(i,j,bottom) = curl(i,j,bottom)-coefficient*geom%area(i,j)/geom%h**2*omega(i,j,bottom)
      end do
    end do
    torque(:,bottom) = torque(:,bottom)-coefficient*gamma(:,bottom)
  end subroutine bottom_drag

end module qg_irregular_forcing
