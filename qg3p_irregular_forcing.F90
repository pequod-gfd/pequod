! Compatible cut-cell curl and component torque from a parent momentum force.
module qg_irregular_forcing
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use qg_irregular_parallel, only: distributed_geometry, half_coordinates, share_face_flux, broadcast_error
  use qg_irregular_operators, only: sum_open_flux
  use qg_parameters, only: real_8, pi
  implicit none
  include 'mpif.h'
  private
  public :: wind_curl_torque, cape_wind_curl_torque, bottom_drag

  real(real_8), parameter :: gauss(4) = [0.1834346424956498D0,0.5255324099163290D0, &
                                        0.7966664774136267D0,0.9602898564975363D0]
  real(real_8), parameter :: weight(4) = [0.3626837833783620D0,0.3137066458778873D0, &
                                         0.2223810344533745D0,0.1012285362903763D0]

  abstract interface
    pure function line_integral(first, second, params) result(value)
      import real_8
      real(real_8), intent(in) :: first(2),second(2),params(:)
      real(real_8) :: value
    end function line_integral

    pure function zonal_stress(point, params) result(value)
      import real_8
      real(real_8), intent(in) :: point(2),params(:)
      real(real_8) :: value
    end function zonal_stress
  end interface

contains

  ! The legacy piecewise-sinusoidal curl is the negative y derivative
  ! of this continuous zonal acceleration. The coordinates are native.
  ! strength is the layer acceleration used in the legacy PV curl.
  subroutine wind_curl_torque(geom, strength, shift, tilt, asymmetry, curl, torque, error)
    type(distributed_geometry), intent(in) :: geom
    real(real_8), intent(in) :: strength,shift,tilt,asymmetry
    real(real_8), intent(out) :: curl(:,:,:),torque(0:,:)
    character(*), intent(out) :: error

    error = ''
    if (.not. all(ieee_is_finite([strength,shift,tilt,asymmetry])) .or. &
        abs(shift)+0.5D0*abs(tilt) >= 0.5D0 .or. asymmetry <= 0.0D0) error = 'invalid piecewise wind profile'
    call broadcast_error(geom%comm,geom%rank,geom%ranks,error)
    if (len_trim(error) /= 0) return
    call accumulate_wind_curl(geom,double_gyre_line,[strength,shift,tilt,asymmetry, &
                             real(shape(geom%owner)-1,real_8)],curl,torque,error)
  end subroutine wind_curl_torque

  ! The gyre occupies [west,east wall] x [south,north wall] in native coordinates.
  subroutine cape_wind_curl_torque(geom, amplitude, west, south, curl, torque, error)
    type(distributed_geometry), intent(in) :: geom
    real(real_8), intent(in) :: amplitude,west,south
    real(real_8), intent(out) :: curl(:,:,:),torque(0:,:)
    character(*), intent(out) :: error
    real(real_8) :: east,north

    east = size(geom%owner,2)-1
    north = size(geom%owner,1)-1
    error = ''
    if (.not. all(ieee_is_finite([amplitude,west,south])) .or. &
        west < 0.0D0 .or. west >= east .or. south < 0.0D0 .or. south >= north) error = 'invalid cape wind patch'
    call broadcast_error(geom%comm,geom%rank,geom%ranks,error)
    if (len_trim(error) /= 0) return
    call accumulate_wind_curl(geom,cape_gyre_line,[amplitude,west,east,south,north],curl,torque,error)
  end subroutine cape_wind_curl_torque

  ! Use the same oriented line integral for shared open faces and the coast.
  subroutine accumulate_wind_curl(geom, integrate, params, curl, torque, error)
    type(distributed_geometry), intent(in) :: geom
    procedure(line_integral) :: integrate
    real(real_8), intent(in) :: params(:)
    real(real_8), intent(out) :: curl(:,:,:),torque(0:,:)
    character(*), intent(out) :: error
    real(real_8), allocatable :: open(:,:)
    real(real_8) :: local(0:ubound(torque,1),size(torque,2)),root_params(size(params)),p(2),q(2),r(2),d
    integer :: f,k,i,j,g_i,edge,previous,comp,ierr,nl,counts(4)

    error = ''
    nl = size(curl,3)
    if (any(shape(curl) /= [geom%ni,geom%nj,nl]) .or. &
        any(shape(torque) /= [size(geom%global%loops),nl]) .or. nl < 1) error = 'incompatible wind arrays'
    call broadcast_error(geom%comm,geom%rank,geom%ranks,error)
    if (len_trim(error) /= 0) return
    call mpi_allreduce([nl,-nl,size(params),-size(params)],counts,size(counts),MPI_INTEGER,MPI_MAX,geom%comm,ierr)
    if (counts(1) /= -counts(2)) then
      error = 'wind layer counts differ across ranks'
      return
    end if
    if (counts(3) /= -counts(4)) then
      error = 'wind profiles differ across ranks'
      return
    end if
    root_params = params
    call mpi_bcast(root_params,size(root_params),MPI_DOUBLE_PRECISION,0,geom%comm,ierr)
    if (any(root_params /= params)) error = 'wind parameters differ across ranks'
    call broadcast_error(geom%comm,geom%rank,geom%ranks,error)
    if (len_trim(error) /= 0) return

    g_i = size(geom%owner,1)
    allocate(open(1,size(geom%faces)))
    open = 0.0D0
    do f = 1,size(geom%faces)
      if (geom%faces(f)%owner /= geom%rank) cycle
      p = 0.5D0*geom%faces(f)%endpoints(:,1)
      q = 0.5D0*geom%faces(f)%endpoints(:,2)
      open(1,f) = integrate(p,q,params)
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
      d = integrate(0.5D0*(r+p),p,params)+integrate(p,0.5D0*(p+q),params)
      curl(i,j,1) = curl(i,j,1)+d
      local(comp,1) = local(comp,1)+d
    end do
    call mpi_allreduce(local,torque,size(local),MPI_DOUBLE_PRECISION,MPI_SUM,geom%comm,ierr)
  end subroutine accumulate_wind_curl

  ! Integrate zonal stress along the given parameter interval of a face.
  pure function integrate_zonal_stress(first, delta, lo, hi, params, stress) result(integral)
    real(real_8), intent(in) :: first(2),delta(2),lo,hi,params(:)
    procedure(zonal_stress) :: stress
    real(real_8) :: integral,mid,half,point(2)
    integer :: k,s

    mid = 0.5D0*(lo+hi)
    half = 0.5D0*(hi-lo)
    integral = 0.0D0
    do k = 1,size(gauss)
      do s = -1,1,2
        point = first+(mid+s*half*gauss(k))*delta
        integral = integral+half*weight(k)*stress(point,params)
      end do
    end do
    integral = integral*delta(1)
  end function integrate_zonal_stress

  pure function double_gyre_line(first, second, params) result(integral)
    real(real_8), intent(in) :: first(2),second(2),params(:)
    real(real_8) :: integral,delta(2),rates(2),split

    delta = second-first
    integral = 0.0D0
    if (delta(1) == 0.0D0) return
    rates = [first(2),second(2)]/params(5)-0.5D0-params(2)-params(3)*([first(1),second(1)]/params(6)-0.5D0)
    if (minval(rates) < 0.0D0 .and. maxval(rates) > 0.0D0) then
      split = rates(1)/(rates(1)-rates(2))
      integral = integrate_zonal_stress(first,delta,0.0D0,split,params,double_gyre_stress) &
                 +integrate_zonal_stress(first,delta,split,1.0D0,params,double_gyre_stress)
    else
      integral = integrate_zonal_stress(first,delta,0.0D0,1.0D0,params,double_gyre_stress)
    end if
    integral = integral*2.0D0*params(5)*params(1)
  end function double_gyre_line

  pure function double_gyre_stress(point, params) result(value)
    real(real_8), intent(in) :: point(2),params(:)
    real(real_8) :: value,x,y,mid,width

    x = point(1)/params(6)-0.5D0
    y = point(2)/params(5)-0.5D0
    mid = params(2)+params(3)*x
    if (y < mid) then
      width = 0.5D0+mid
      value = -params(4)*width*(cos(pi*(0.5D0+y)/width)+1.0D0)
    else
      width = 0.5D0-mid
      value = width*(cos(pi*(y-mid)/width)-1.0D0)/params(4)
    end if
  end function double_gyre_stress

  pure function cape_gyre_line(first, second, params) result(integral)
    real(real_8), intent(in) :: first(2),second(2),params(:)
    real(real_8) :: integral,delta(2),lo,hi,bounds(2),crossings(2)
    integer :: axis

    integral = 0.0D0
    delta = second-first
    if (delta(1) == 0.0D0) return
    lo = 0.0D0
    hi = 1.0D0
    do axis = 1,2
      bounds = params(2*axis:2*axis+1)
      if (delta(axis) == 0.0D0) then
        if (first(axis) < bounds(1) .or. first(axis) > bounds(2)) return
      else
        crossings = (bounds-first(axis))/delta(axis)
        lo = max(lo,minval(crossings))
        hi = min(hi,maxval(crossings))
      end if
    end do
    if (hi > lo) integral = integrate_zonal_stress(first,delta,lo,hi,params,cape_gyre_stress)
  end function cape_gyre_line

  pure function cape_gyre_stress(point, params) result(value)
    real(real_8), intent(in) :: point(2),params(:)
    real(real_8) :: value,x,y

    ! params holds peak amplitude, west, east, south and north.
    x = (point(1)-params(2))/(params(3)-params(2))
    y = (point(2)-params(4))/(params(5)-params(4))
    value = -params(1)*(sin(pi*x)*sin(0.5D0*pi*y))**2
  end function cape_gyre_stress

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
