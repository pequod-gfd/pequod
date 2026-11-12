program test_weak_diffusion
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use qg_irregular_parallel, only: distributed_geometry, load_distributed_geometry
  use qg_weak_diffusion, only: weak_diffusion_cache, setup_weak_diffusion, advance_weak_diffusion, release_weak_diffusion
  use qg_parallel, only: decomposition, allocate, deallocate, parallel_init, parallel_finalise, pabort
  use qg_parameters, only: real_8
  implicit none
  include 'mpif.h'

  real(real_8), parameter :: finish=0.03125D0,physical_viscosity=0.01D0
  real(real_8), parameter :: nodes(4) = [0.06943184420297371D0,0.3300094782075719D0, &
                                          0.6699905217924281D0,0.9305681557970262D0]
  real(real_8), parameter :: weights(4) = [0.1739274225687269D0,0.3260725774312731D0, &
                                            0.3260725774312731D0,0.1739274225687269D0]
  real(real_8), parameter :: island_nodes(12) = [0.0092196828766403782D0,0.047941371814762601D0,0.11504866290284765D0, &
    0.20634102285669126D0,0.31608425050090994D0,0.43738329574426554D0, &
    0.5626167042557344D0,0.68391574949909006D0,0.79365897714330869D0, &
    0.88495133709715235D0,0.95205862818523745D0,0.99078031712335957D0]
  real(real_8), parameter :: island_weights(12) = [0.023587668193255706D0,0.053469662997659533D0,0.080039164271673208D0, &
    0.10158371336153287D0,0.11674626826917731D0,0.12457352290670134D0, &
    0.12457352290670134D0,0.11674626826917731D0,0.10158371336153287D0, &
    0.080039164271673208D0,0.053469662997659533D0,0.023587668193255706D0]
  real(real_8), parameter :: centres(2)=[0.25D0,0.75D0],levels(2)=[0.6D0,-0.35D0]
  real(real_8), parameter :: amplitudes(2)=[2.0D0,-1.25D0],radius=0.0625D0,normal_width=0.25D0
  real(real_8), parameter :: ramp(0:7)=[0,0,0,0,35,-84,70,-20]*1.0D0
  real(real_8), parameter :: tangent(0:10)=[1,0,-5,0,10,0,-10,0,5,0,-1]*1.0D0
  real(real_8), parameter :: rays(2,0:7) = 0.5D0*reshape([1,0,1,1,0,1,-1,1,-1,0,-1,-1,0,-1,1,-1],[2,8])
  type(decomposition), target :: d
  type(distributed_geometry), target :: geom
  type(weak_diffusion_cache) :: cache
  real(real_8), allocatable :: psi(:,:,:),old(:,:,:),source(:,:,:),base(:,:,:),omega(:,:,:),wall(:,:,:)
  integer, allocatable :: iis(:),jjs(:)
  real(real_8) :: coeff(0:4,0:4),normal(0:13),x,y,h,dt,t,factor,local(4),global(4),energy,expected,contour
  real(real_8), allocatable :: torque(:),quadrature(:),quad_weights(:)
  integer :: grid(2),ranks,ierr,n,steps,i,j,k,gi,gj,its,max_iterations
  character(1024) :: path,error
  character(20) :: fixture,arg

  if (command_argument_count() /= 4) then
    print '(a)', 'Weak diffusion fixture requires geometry, dimension, steps and rectangle/triangle'
    stop 2
  end if
  call get_command_argument(1,path)
  call get_command_argument(2,arg)
  read(arg,*) n
  call get_command_argument(3,arg)
  read(arg,*) steps
  call get_command_argument(4,fixture)

  coeff = 0.0D0
  select case(trim(fixture))
    case('islands')
      quadrature = island_nodes
      quad_weights = island_weights
      normal = 0.0D0
      normal(2:12:2) = tangent(0:10:2)
      normal(3:13:2) = tangent(0:10:2)
    case('rectangle')
      do j = 2,4
        do i = 2,4
          coeff(i,j) = merge(-2,1,i==3)*merge(-2,1,j==3)
        end do
      end do
    case('triangle')
      coeff(2,2) = 1.0D0
      coeff(4,2) = 1.0D0
      coeff(2,4) = 1.0D0
      coeff(3,2) = -2.0D0
      coeff(2,3) = -2.0D0
      coeff(3,3) = 2.0D0
    case default
      error stop 'unknown weak diffusion fixture'
  end select
  if (.not. allocated(quadrature)) then
    quadrature = nodes
    quad_weights = weights
  end if
  h = 1.0D0/(n-1)
  call parallel_init()
  call mpi_comm_size(MPI_COMM_WORLD,ranks,ierr)
  grid = 0
  call mpi_dims_create(ranks,2,grid,ierr)
  iis = [(i*n/grid(1)-(i-1)*n/grid(1),i=1,grid(1))]
  jjs = [(j*n/grid(2)-(j-1)*n/grid(2),j=1,grid(2))]
  call allocate(d,iis=iis,jjs=jjs,layers=1)
  call load_distributed_geometry(path,n,n,[d%ii_0,d%jj_0,d%ii,d%jj],d%cart_comm,h,geom,error)
  call check(error)
  if (fixture=='islands' .and. size(geom%global%loops)/=3) call check('expected two square islands')
  call setup_weak_diffusion(cache,d,geom,error)
  call check(error)
  allocate(psi(geom%ni,geom%nj,1),old(geom%ni,geom%nj,1),source(geom%ni,geom%nj,1),base(geom%ni,geom%nj,1), &
           omega(geom%ni,geom%nj,1),wall(geom%ni,geom%nj,1))
  allocate(torque(0:ubound(geom%global%loops,1)))
  contour = 2*radius
  do k = 1,5
    contour = contour*(2*k)/real(2*k+1,real_8)
  end do
  psi = ieee_value(0.0D0,ieee_quiet_nan)
  base = psi
  do j = 1,geom%nj
    gj = geom%j0+j
    x = (gj-1)*h
    do i = 1,geom%ni
      gi = geom%i0+i
      y = (n-gi)*h
      if (geom%area(i,j)>0.0D0) psi(i,j,1) = derivative(x,y,0,0)
      if (geom%area(i,j)>0.0D0) base(i,j,1) = curl_integral([x,y],geom%global%wedge(gi,gj))
    end do
  end do
  max_iterations = 0
  t = 0.0D0
  do k = 1,steps
    dt = finish/steps
    if (mod(steps,2)==0) dt = dt*merge(0.75D0,1.25D0,mod(k,2)==1)
    factor = (exp(-t)-exp(-t-dt))/dt
    old = psi
    source = factor*base
    torque = 0.0D0
    if (fixture=='islands') torque(1:) = physical_viscosity*factor*amplitudes*6*contour/normal_width
    call advance_weak_diffusion(cache,geom%comm,old,source,dt,physical_viscosity/h**2,psi,omega,wall,its,error,torque)
    call check(error)
    max_iterations = max(max_iterations,its)
    t = t+dt
    energy = sum(geom%area*omega(:,:,1))/h**2
    call mpi_allreduce(energy,expected,1,MPI_DOUBLE_PRECISION,MPI_SUM,geom%comm,ierr)
    if (abs(expected)>1.0D-12) call check('no-slip weak Stokes total is nonzero')
  end do

  local = 0.0D0
  do j = 1,geom%nj
    gj = geom%j0+j
    x = (gj-1)*h
    do i = 1,geom%ni
      gi = geom%i0+i
      y = (n-gi)*h
      if (geom%area(i,j)==0.0D0) cycle
      local(1) = max(local(1),abs(psi(i,j,1)-exp(-finish)*derivative(x,y,0,0)))
      local(2) = local(2)+source(i,j,1)+physical_viscosity/h**2*wall(i,j,1)
      local(3) = local(3)+psi(i,j,1)**2
      if (gi==n .and. abs(x-0.5D0)<0.15D0) then
        expected = -factor*(derivative(x,y,2,1)+derivative(x,y,0,3))
        local(4) = max(local(4),abs(wall(i,j,1)/h**3-expected))
      end if
    end do
  end do
  call mpi_allreduce(local(1),global(1),1,MPI_DOUBLE_PRECISION,MPI_MAX,geom%comm,ierr)
  call mpi_allreduce(local(2:3),global(2:3),2,MPI_DOUBLE_PRECISION,MPI_SUM,geom%comm,ierr)
  call mpi_allreduce(local(4),global(4),1,MPI_DOUBLE_PRECISION,MPI_MAX,geom%comm,ierr)
  if (abs(global(2))>2.0D-10) call check('coupled no-slip source budget did not close')
  if (geom%rank==0) print '(a,4(1x,es23.15),1x,i0,1x,a)', 'weak-diffusion-errors', global,max_iterations,'PASS'
  call release_weak_diffusion(cache)
  call deallocate(d)
  call parallel_finalise()

contains

  function derivative(x,y,dx,dy) result(value)
    real(real_8), intent(in) :: x,y
    integer, intent(in) :: dx,dy
    real(real_8) :: value,term
    integer :: a,b,k

    value = 0.0D0
    if (fixture=='islands') then
      do a = 1,2
        value = value+levels(a)*plateau(x-centres(a),radius,0.25D0,dx)*plateau(y-0.5D0,radius,0.5D0,dy)
        if (abs(x-centres(a))<radius .and. y>=0.5D0+radius .and. y<0.5D0+radius+normal_width) then
          value = value+amplitudes(a)*polynomial(tangent,(x-centres(a))/radius,dx)/radius**dx* &
                  polynomial(normal,(y-0.5D0-radius)/normal_width,dy)*normal_width**(2-dy)
        end if
      end do
      return
    end if
    do b = dy,4
      do a = dx,4
        term = coeff(a,b)
        do k = 0,dx-1
          term = term*(a-k)
        end do
        do k = 0,dy-1
          term = term*(b-k)
        end do
        value = value+term*x**(a-dx)*y**(b-dy)
      end do
    end do
  end function derivative

  function plateau(point,inner,outer,order) result(value)
    real(real_8), intent(in) :: point,inner,outer
    integer, intent(in) :: order
    real(real_8) :: value,t,scale

    value = 0.0D0
    if (abs(point)<=inner) then
      value = merge(1.0D0,0.0D0,order==0)
    else if (abs(point)<outer) then
      scale = outer-inner
      t = (outer-abs(point))/scale
      value = polynomial(ramp,t,order)*(-sign(1.0D0,point)/scale)**order
    end if
  end function plateau

  function polynomial(coefficients,t,order) result(value)
    real(real_8), intent(in) :: coefficients(0:),t
    integer, intent(in) :: order
    real(real_8) :: value,term
    integer :: j,k

    value = 0.0D0
    do j = ubound(coefficients,1),order,-1
      term = coefficients(j)
      do k = 0,order-1
        term = term*(j-k)
      end do
      value = t*value+term
    end do
  end function polynomial

  function curl_integral(origin,mask) result(value)
    real(real_8), intent(in) :: origin(2)
    integer, intent(in) :: mask
    real(real_8) :: value,p(2),a(2),b(2),x,y,omega,delta_omega
    integer :: k,i,j

    value = 0.0D0
    do k = 0,7
      if (.not. btest(mask,k)) cycle
      a = rays(:,k)
      b = rays(:,modulo(k+1,8))
      do i = 1,size(quadrature)
        do j = 1,size(quadrature)
          p = origin+h*quadrature(i)*(a+quadrature(j)*(b-a))
          x = p(1)
          y = p(2)
          omega = derivative(x,y,2,0)+derivative(x,y,0,2)
          delta_omega = derivative(x,y,4,0)+2*derivative(x,y,2,2)+derivative(x,y,0,4)
          value = value-0.25D0*h**2*quadrature(i)*quad_weights(i)*quad_weights(j)*(omega+physical_viscosity*delta_omega)
        end do
      end do
    end do
  end function curl_integral

  subroutine check(message)
    character(*), intent(in) :: message

    if (len_trim(message)>0) call pabort(message,__FILE__,__LINE__)
  end subroutine check
end program test_weak_diffusion
