program test_irregular_vertex_order
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  use qg_cmm, only: cmm_cache, setup_cmm, cmm_solve, release_cmm
  use qg_geometry, only: fluid_node
  use qg_irregular_bc, only: vertex_plan, setup_vertex_plan, instantaneous_transport, extrapolated_transport
  use qg_irregular_operators, only: sum_open_flux
  use qg_irregular_parallel, only: distributed_geometry, load_distributed_geometry
  use qg_parallel, only: decomposition, allocate, deallocate, parallel_init, parallel_finalise, pabort
  use qg_parameters, only: real_8
  implicit none
  include 'mpif.h'

  type(decomposition), target :: decomp
  type(distributed_geometry), target :: geom
  type(vertex_plan) :: plan
  type(cmm_cache) :: inverse
  integer :: n, i,j,gi,gj,f,k, ranks,grid(2),ierr, component
  integer, allocatable :: iis(:),jjs(:)
  real(real_8), allocatable :: psi(:,:,:), wall(:,:), current(:,:), past(:,:), flux(:,:), divergence(:,:,:), rhs(:,:,:), velocity(:,:,:)
  real(real_8) :: h, x,y,nan,length, reference, error_local(4), error_global(4), p(2,2),normal(2),analytic
  real(real_8), parameter :: weight = 0.3D0, previous = 0.625D0
  character(1024) :: path,error
  character(30) :: fixture,arg
  logical :: use_inverse

  if (command_argument_count() /= 4) then
    print '(a)', 'Vertex order fixture requires .dat, dimension, case and direct/inverted'
    stop 2
  end if
  call get_command_argument(1,path)
  call get_command_argument(2,arg)
  read(arg,*) n
  call get_command_argument(3,fixture)
  call get_command_argument(4,arg)
  use_inverse = arg == 'inverted'
  call parallel_init()
  call mpi_comm_size(MPI_COMM_WORLD,ranks,ierr)
  grid = 0
  call mpi_dims_create(ranks,2,grid,ierr)
  iis = [(i*n/grid(1)-(i-1)*n/grid(1),i=1,grid(1))]
  jjs = [(j*n/grid(2)-(j-1)*n/grid(2),j=1,grid(2))]
  call allocate(decomp,iis=iis,jjs=jjs,layers=1)
  h = 1.0D0/(n-1)
  call load_distributed_geometry(path,n,n,[decomp%ii_0,decomp%jj_0,decomp%ii,decomp%jj],decomp%cart_comm,h,geom,error)
  call require(len_trim(error) == 0,error)
  call setup_vertex_plan(geom,plan,error,normal_profile=fixture=="cubic")
  call require(len_trim(error) == 0,error)
  allocate(psi(geom%ni,geom%nj,1),divergence(geom%ni,geom%nj,1),wall(0:ubound(geom%global%loops,1),1), &
           current(1,size(geom%faces)),past(1,size(geom%faces)),flux(1,size(geom%faces)))
  do component = 0,ubound(wall,1)
    wall(component,1) = 0.375D0+0.25D0*component
  end do
  nan = ieee_value(0.0D0,ieee_quiet_nan)
  psi = nan
  do j = 1,geom%nj
    gj = geom%j0+j
    do i = 1,geom%ni
      gi = geom%i0+i
      if (geom%global%kind(gi,gj) /= fluid_node) cycle
      psi(i,j,1) = exact(h*(gj-1),h*(n-gi))
    end do
  end do
  if (use_inverse) then
    allocate(rhs(geom%ni,geom%nj,1))
    rhs = 0.0D0
    do j = 1,geom%nj
      gj = geom%j0+j
      do i = 1,geom%ni
        gi = geom%i0+i
        if (geom%global%kind(gi,gj) /= fluid_node) cycle
        x = h*(gj-1)
        y = h*(n-gi)
        rhs(i,j,1) = exact(x-h,y)+exact(x+h,y)+exact(x,y-h)+exact(x,y+h)-4*psi(i,j,1)
      end do
    end do
    call setup_cmm(inverse,decomp,geom,[0.0D0],error)
    call require(len_trim(error) == 0,error)
    call cmm_solve(inverse,rhs,wall,psi,error)
    call require(len_trim(error) == 0,error)
    call release_cmm(inverse)
  end if

  if (fixture == 'cubic') then
    allocate(velocity(2,1,size(geom%faces)))
    call instantaneous_transport(geom,plan,psi,wall,current,error,normal_velocity=velocity)
    call require(len_trim(error) == 0,error)
    do f = 1,size(geom%faces)
      if (geom%faces(f)%owner /= geom%rank) cycle
      p = 0.5D0*h*geom%faces(f)%endpoints
      length = sum(abs(p(:,2)-p(:,1)))
      normal = (p(:,2)-p(:,1))/length
      do k = 1,2
        x = p(1,k)
        y = p(2,k)
        analytic = -dot_product(normal,[y*(1-2*x-y),x*(1-x-2*y)])
        call require(abs(velocity(k,1,f)/h-analytic)<2e-11,'endpoint cubic velocity reproduction')
      end do
      call require(abs(sum(velocity(:,1,f))*length/(2*h)-current(1,f))<2e-15,'profile does not match vertex transport')
    end do
  else
    call instantaneous_transport(geom,plan,psi,wall,current,error)
    call require(len_trim(error) == 0,error)
  end if
  call instantaneous_transport(geom,plan,previous*psi,previous*wall,past,error)
  call require(len_trim(error) == 0,error)
  call extrapolated_transport(geom,plan,psi,wall,weight,past,flux,error)
  call require(len_trim(error) == 0,error)
  call require(all(ieee_is_finite(current)).and.all(ieee_is_finite(flux)),'vertex read poisoned values')
  call require(all(past == current),'extrapolated vertex did not advance instantaneous history')
  error_local = 0.0D0
  do f = 1,size(geom%faces)
    p = 0.5D0*h*geom%faces(f)%endpoints
    length = sum(abs(p(:,1)-p(:,2)))
    reference = exact(p(1,1),p(2,1))-exact(p(1,2),p(2,2))
    error_local(1) = max(error_local(1),abs(current(1,f)-reference)/length)
    error_local(2) = max(error_local(2),abs(flux(1,f)-(1+weight-weight*previous)*reference)/length)
  end do
  call sum_open_flux(geom,current,divergence)
  error_local(3) = maxval(abs(divergence))
  call sum_open_flux(geom,flux,divergence)
  error_local(4) = maxval(abs(divergence))
  call mpi_allreduce(error_local,error_global,4,MPI_DOUBLE_PRECISION,MPI_MAX,geom%comm,ierr)
  call require(maxval(error_global(3:)) < 5e-14,'shared vertex transport failed cell telescoping')
  if (fixture == 'cubic') call require(maxval(error_global(:2)) < 1e-10,'compatible cubic velocity reproduction')
  if (geom%rank == 0) print '(a,4(1x,es23.15),1x,a)', 'vertex-errors', error_global,'PASS'
  call deallocate(decomp)
  call parallel_finalise()

contains
  real(real_8) function exact(x,y) result(value)
    real(real_8), intent(in) :: x,y
    real(real_8) :: r,t, inner,outer,cx,cy
    integer :: island,k

    value = 0.375D0
    select case(fixture)
      case('rectangle')
        value = value+20*(x*(1-x)*y*(1-y))**2
      case('triangle')
        value = value+20*(x*y*(1-x-y))**2
      case('cubic')
        value = value+x*y*(1-x-y)
      case('diamond')
        value = value+100*( (0.375D0-x-y+1)*(0.375D0+x+y-1)*(0.375D0-x+y)*(0.375D0+x-y) )**2
      case('concave')
        value = value+100*(x*y*(1-x)*(1-y)*(x-0.5D0)*(y-0.5D0))**2
      case('islands','near-islands')
        inner = 0.13D0
        outer = 0.19D0
        if (fixture == 'near-islands') outer = 0.15D0
        do island = 1,2
          cx = merge(0.34375D0,0.75D0,island==1)
          cy = merge(0.4375D0,0.65625D0,island==1)
          if (fixture == 'near-islands') then
            if (island == 2) cx = 0.65625D0
            cy = 0.5D0
          end if
          r = hypot(x-cx,y-cy)
          t = max(0.0D0,min(1.0D0,(r-inner)/(outer-inner)))
          ! Component 1 is the first (northmost) loop in the generated data.
          k = 3-island
          if (fixture == 'near-islands') k = island
          value = value+0.25D0*k*(1-10*t**3+15*t**4-6*t**5)
        end do
      case default
        call pabort('unknown vertex geometry',__FILE__,__LINE__)
    end select
  end function exact

  subroutine require(ok,message)
    logical, intent(in) :: ok
    character(*), intent(in) :: message
    if (.not. ok) call pabort(message,__FILE__,__LINE__)
  end subroutine require
end program test_irregular_vertex_order
