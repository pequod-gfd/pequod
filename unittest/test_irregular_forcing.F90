program test_irregular_forcing
  use qg_irregular_parallel, only: distributed_geometry, load_distributed_geometry
  use qg_geometry, only: write_rectangular_cape
  use qg_irregular_forcing, only: wind_curl_torque, cape_wind_curl_torque, bottom_drag
  use qg_irregular_transport, only: source_momentum_budget
  use qg_parallel, only: parallel_init, parallel_finalise, pabort
  use qg_parameters, only: real_8, pi
  implicit none
  include 'mpif.h'

  integer, parameter :: n=33,layers=2
  integer, parameter :: cape_left=5,cape_width=8,cape_height=18
  real(real_8), parameter :: strength=0.7D0,asym=0.8D0,drag=0.125D0
  real(real_8), parameter :: rays(2,0:7) = reshape([1,0,1,1,0,1,-1,1,-1,0,-1,-1,0,-1,1,-1],[2,8])*0.5D0
  real(real_8), parameter :: gauss(6) = [0.03376524289842399D0,0.1693953067668677D0,0.3806904069584016D0, &
                                         0.6193095930415985D0,0.8306046932331323D0,0.966234757101576D0]
  real(real_8), parameter :: weight(6) = [0.08566224618958517D0,0.1803807865240693D0,0.2339569672863455D0, &
                                          0.2339569672863455D0,0.1803807865240693D0,0.08566224618958517D0]
  type(distributed_geometry) :: geom
  real(real_8), allocatable :: curl(:,:,:),torque(:,:),base(:,:),omega(:,:,:),gamma(:,:)
  real(real_8) :: p(2),q(2),a(2),b(2),u,y,expected,local_error,global_error,residual(layers)
  real(real_8) :: first,second,theta,delta,width,amplitude,direct
  integer :: comm,rank,ranks,ierr,grid(2),box(4),row,column,i,j,sector,ii,jj,mask,component,k,next
  character(1024) :: path,error
  logical :: cape_case

  if (command_argument_count() /= 1) then
    print '(a)', 'Irregular forcing fixture requires a .dat path or --cape'
    stop 2
  end if
  call get_command_argument(1,path)
  cape_case = path=='--cape'
  call parallel_init()
  call mpi_comm_rank(MPI_COMM_WORLD,rank,ierr)
  call mpi_comm_size(MPI_COMM_WORLD,ranks,ierr)
  if (cape_case) then
    path = 'cape-wind.dat'
    if (rank==0) then
      call write_rectangular_cape(path,n,n,cape_left,cape_width,cape_height,2.5D0,error)
      call require(len_trim(error)==0,error)
    end if
    call mpi_barrier(MPI_COMM_WORLD,ierr)
  end if
  call mpi_comm_split(MPI_COMM_WORLD,0,ranks-1-rank,comm,ierr)
  grid = 0
  call mpi_dims_create(ranks,2,grid,ierr)
  row = rank/grid(2)
  column = mod(rank,grid(2))
  box = [row*n/grid(1),column*n/grid(2),(row+1)*n/grid(1)-row*n/grid(1),(column+1)*n/grid(2)-column*n/grid(2)]
  call load_distributed_geometry(path,n,n,box,comm,2.5D0,geom,error)
  call require(len_trim(error)==0,error)
  allocate(curl(geom%ni,geom%nj,layers),omega(geom%ni,geom%nj,layers),base(geom%ni,geom%nj))
  allocate(torque(0:ubound(geom%global%loops,1),layers),gamma(0:ubound(geom%global%loops,1),layers))
  if (cape_case) then
    call check_cape_wind()
    call mpi_comm_free(comm,ierr)
    call parallel_finalise()
    stop
  end if

  call wind_curl_torque(geom,strength,0.0D0,0.0D0,asym,curl,torque,error)
  call require(len_trim(error)==0,error)
  call require(maxval(abs(curl(:,:,2)))==0.0D0 .and. maxval(abs(torque(:,2)))==0.0D0,'wind entered lower layers')
  local_error = 0.0D0
  do j = 1,geom%nj
    do i = 1,geom%ni
      mask = geom%global%wedge(geom%i0+i,geom%j0+j)
      p = [real(geom%j0+j-1,real_8),real(n-geom%i0-i,real_8)]
      expected = 0.0D0
      do sector = 0,7
        if (.not. btest(mask,sector)) cycle
        a = rays(:,sector)
        b = rays(:,modulo(sector+1,8))
        do ii = 1,size(gauss)
          u = gauss(ii)
          do jj = 1,size(gauss)
            q = p+u*(a+gauss(jj)*(b-a))
            y = q(2)/(n-1)-0.5D0
            if (y < 0.0D0) then
              direct = -asym*sin(pi*(1+2*y))
            else
              direct = sin(2*pi*y)/asym
            end if
            expected = expected+0.25D0*u*weight(ii)*weight(jj)*strength*2*pi*direct
          end do
        end do
      end do
      local_error = max(local_error,abs(curl(i,j,1)-expected))
    end do
  end do
  call mpi_allreduce(local_error,global_error,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
  call require(global_error<6.0D-13,'wind curl disagrees with independent volume quadrature')

  do component = 0,ubound(torque,1)
    expected = 0.0D0
    associate(loop=>geom%global%loops(component)%ij)
      do k = 1,size(loop,2)
        next = mod(k,size(loop,2))+1
        first = real(n-loop(1,k),real_8)/(n-1)-0.5D0
        second = real(n-loop(1,next),real_8)/(n-1)-0.5D0
        width = 0.5D0
        theta = pi*first/width
        delta = pi*(second-first)/width
        direct = cos(theta)
        if (delta /= 0) direct = (sin(theta+delta)-sin(theta))/delta
        amplitude = merge(asym,1/asym,first+second<0.0D0)
        expected = expected+2*strength*(n-1)*width*amplitude*(direct-1)*(loop(2,next)-loop(2,k))
      end do
    end associate
    call require(abs(expected-torque(component,1))<1.0D-10,'component wind torque disagrees with analytic line integral')
  end do

  base = curl(:,:,1)
  omega = 0.0D0
  where (geom%area>0.0D0) omega(:,:,layers) = base*geom%h**2/max(geom%area,tiny(0.0D0))
  gamma = 0.0D0
  gamma(:,layers) = torque(:,1)
  call bottom_drag(geom,drag,omega,gamma,curl,torque)
  call require(maxval(abs(curl(:,:,layers)+drag*base))<1.0D-13,'bottom drag did not act on the cell-mean omega')
  call source_momentum_budget(geom,curl,torque,residual,error)
  call require(len_trim(error)==0,error)
  call wind_curl_torque(geom,strength,0.13D0,0.29D0,asym,curl,torque,error)
  call require(len_trim(error)==0,error)
  call source_momentum_budget(geom,curl,torque,residual,error)
  call require(len_trim(error)==0,error)
  if (rank==0) print '(a,2es14.6,a)', 'irregular wind and drag errors: ', global_error,maxval(abs(residual)),' PASS'

  call mpi_comm_free(comm,ierr)
  call parallel_finalise()

contains

  subroutine check_cape_wind()
    real(real_8), parameter :: west=cape_left+cape_width-1, east=n-1, south=n-1-cape_height, north=n-1
    real(real_8) :: x,y,xlo,xhi,ylo,yhi,zonal,meridional,expected,local_error,global_error
    integer :: i,j

    call cape_wind_curl_torque(geom,strength,west,south,curl,torque,error)
    call require(len_trim(error)==0,error)
    call require(all(curl(:,:,2)==0.0D0) .and. all(torque(:,2)==0.0D0),'cape wind entered lower layers')
    local_error = 0.0D0
    do j = 1,geom%nj
      x = geom%j0+j-1
      xlo = max(west,x-0.5D0)
      xhi = min(east,x+0.5D0)
      do i = 1,geom%ni
        y = n-geom%i0-i
        ylo = max(south,y-0.5D0)
        yhi = min(north,y+0.5D0)
        expected = 0.0D0
        if (xlo<xhi .and. ylo<yhi) then
          ! The entire patch is water. Integrate the analytic curl on the
          ! intersection with the full square even at partially wet cells.
          zonal = (xhi-xlo)/2-(east-west)/(4*pi)* &
                  (sin(2*pi*(xhi-west)/(east-west))-sin(2*pi*(xlo-west)/(east-west)))
          meridional = (cos(pi*(ylo-south)/(north-south))-cos(pi*(yhi-south)/(north-south)))/2
          expected = strength*zonal*meridional
        else
          call require(curl(i,j,1)==0.0D0,'cape wind curl leaked outside the patch')
        end if
        local_error = max(local_error,abs(curl(i,j,1)-expected))
      end do
    end do
    call mpi_allreduce(local_error,global_error,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call require(global_error<5.0D-14,'cape wind disagrees with exact cell curl integral')
    call require(abs(torque(0,1)-strength*(east-west)/2)<5.0D-14,'cape wind has incorrect north-wall torque')
    call source_momentum_budget(geom,curl,torque,residual,error)
    call require(len_trim(error)==0,error)
    base = curl(:,:,1)
    call cape_wind_curl_torque(geom,-strength,west,south,curl,torque,error)
    call require(len_trim(error)==0,error)
    call require(all(curl(:,:,1)==-base),'negative stress did not reverse the cape wind curl')
    if (rank==0) print '(a,2es14.6,a)', 'cape wind errors: ', global_error,maxval(abs(residual)),' PASS'
  end subroutine check_cape_wind

  subroutine require(ok,message)
    logical, intent(in) :: ok
    character(*), intent(in) :: message
    if (.not. ok) call pabort(message,__FILE__,__LINE__)
  end subroutine require

end program test_irregular_forcing
