program test_irregular_island_constraints
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  use qg_diagnostics, only: stratification, allocate, deallocate
  use qg_geometry, only: fluid_node
  use qg_irregular_parallel, only: distributed_geometry, load_distributed_geometry
  use qg_island_constraints
  use qg_parallel, only: decomposition, allocate, deallocate, parallel_init, parallel_finalise, pabort
  use qg_parameters, only: real_8
  implicit none
  include 'mpif.h'

  real(real_8), parameter :: outer(2) = [0.2D0,-0.2D0]
  real(real_8), parameter :: expected(2,2) = reshape([0.4D0,-0.2D0,-0.3D0,0.6D0],[2,2])
  real(real_8), parameter :: amplitude(2,2) = reshape([0.7D0,-0.4D0,-0.55D0,1.0D0],[2,2])
  real(real_8), parameter :: centres(2) = [0.25D0,0.75D0], y0 = 0.5D0, edge = 0.5625D0
  real(real_8), parameter :: tangent_width = 0.05D0, normal_width = 0.08D0
  real(real_8), parameter :: rays(2,0:7) = reshape([1,0,1,1,0,1,-1,1,-1,0,-1,-1,0,-1,1,-1],[2,8])
  real(real_8), parameter :: qpoint(6) = [0.033765242898423986D0,0.16939530676686776D0,0.38069040695840156D0, &
                                         0.61930959304159849D0,0.83060469323313224D0,0.96623475710157603D0]
  real(real_8), parameter :: qweight(6) = [0.08566224618958517D0,0.18038078652406930D0,0.23395696728634552D0, &
                                          0.23395696728634552D0,0.18038078652406930D0,0.08566224618958517D0]
  type(decomposition), target :: d
  type(distributed_geometry), target :: g
  type(stratification), target :: s
  type(weak_island_cache) :: cache
  integer :: i,j,gi,gj,n,ranks,grid(2),ierr,kind,k,offset
  integer, allocatable :: iis(:),jjs(:)
  integer :: have(2,4),sample_count(2,4)
  real(real_8) :: h,x,y,psi0(2),vort(2),p(2),gamma(2,2),wrong(2,2),mass(1),local_mass,contour_scale
  real(real_8) :: samples(2,2,4),near_wall(2,2,4),slip_error,mean_psi(2)
  logical :: physical_means
  real(real_8) :: wall(0:2,2),shifted(0:2,2),error_max,local_error,backward_error
  real(real_8), allocatable :: q(:,:,:),psi(:,:,:)
  character(1024) :: path,error
  character(30) :: arg

  if (command_argument_count() < 2 .or. command_argument_count() > 3) then
    print '(a)', 'Weak island fixture requires a .dat path and dimension'
    stop 2
  end if
  call get_command_argument(1,path)
  call get_command_argument(2,arg)
  read(arg,*) n
  physical_means = .false.
  if (command_argument_count() == 3) then
    call get_command_argument(3,arg)
    if (arg /= 'averages') error stop 'unknown island fixture input representation'
    physical_means = .true.
  end if
  h = 1.0D0/(n-1)
  call parallel_init()
  call mpi_comm_size(MPI_COMM_WORLD,ranks,ierr)
  grid = 0
  call mpi_dims_create(ranks,2,grid,ierr)
  iis = [(i*n/grid(1)-(i-1)*n/grid(1),i=1,grid(1))]
  jjs = [(j*n/grid(2)-(j-1)*n/grid(2),j=1,grid(2))]
  call allocate(d,iis=iis,jjs=jjs,layers=2)
  call load_distributed_geometry(path,n,n,[d%ii_0,d%jj_0,d%ii,d%jj],d%cart_comm,h,g,error)
  call require(len_trim(error)==0,error)
  call require(size(g%global%loops)==3,'fixture requires two distinct islands')
  call allocate(s,[1.0D0,1.0D0],[0.5D0*h*h],[0.5D0*h*h])
  call setup_weak_islands(cache,d,g,s,error)
  call require(len_trim(error)==0,error)
  allocate(q(g%ni,g%nj,2),psi(g%ni,g%nj,2))
  q = ieee_value(0.0D0,ieee_quiet_nan)
  local_mass = 0.0D0
  do j = 1,g%nj
    gj = g%j0+j
    x = h*(gj-1)
    do i = 1,g%ni
      gi = g%i0+i
      if (g%area(i,j)==0) cycle
      y = h*(n-gi)
      call exact(x,y,psi0,vort)
      mean_psi = psi0
      kind = g%global%kind(gi,gj)
      if (kind == fluid_node .and. .not. physical_means) then
        q(i,j,:) = -4*psi0+matmul(s%s_mat,psi0)
        call exact(x-h,y,p,vort)
        q(i,j,:) = q(i,j,:)+p
        call exact(x+h,y,p,vort)
        q(i,j,:) = q(i,j,:)+p
        call exact(x,y-h,p,vort)
        q(i,j,:) = q(i,j,:)+p
        call exact(x,y+h,p,vort)
        q(i,j,:) = q(i,j,:)+p
      else
        call mean_pv([x,y],g%global%wedge(gi,gj),q(i,j,:),p)
        if (physical_means) mean_psi = p
      end if
      local_mass = local_mass+g%area(i,j)*dot_product(s%z_to_modes(2,:),mean_psi)
    end do
  end do
  call mpi_allreduce(local_mass,mass,1,MPI_DOUBLE_PRECISION,MPI_SUM,g%comm,ierr)
  contour_scale = 2*tangent_width
  do k = 1,5
    contour_scale = contour_scale*(2*k)/real(2*k+1,real_8)
  end do
  gamma = -contour_scale*transpose(amplitude)
  if (ranks > 1) then
    wrong = gamma
    if (g%rank == 1) wrong(1,2) = wrong(1,2)+0.1D0
    call invert_weak_islands(cache,g%comm,q,wrong,mass,psi,wall,backward_error,error)
    call require(error=='weak island targets differ across ranks','rank-dependent island targets were accepted')
  end if
  call invert_weak_islands(cache,g%comm,q,gamma,mass,psi,wall,backward_error,error)
  call require(len_trim(error)==0,error)
  call require(all(ieee_is_finite(wall)),'nonfinite island constants')
  local_error = 0.0D0
  samples = 0.0D0
  have = 0
  do j = 1,g%nj
    x = h*(g%j0+j-1)
    do i = 1,g%ni
      if (g%area(i,j)==0) cycle
      gi = g%i0+i
      y = h*(n-gi)
      call exact(x,y,psi0,vort)
      local_error = max(local_error,maxval(abs(psi(i,j,:)-psi0)))
      offset = n-gi-nint(edge/h)
      if (offset < 1 .or. offset > 4) cycle
      do k = 1,2
        if (g%j0+j-1 /= nint(centres(k)/h)) cycle
        samples(k,:,offset) = psi(i,j,:)
        have(k,offset) = 1
      end do
    end do
  end do
  call mpi_allreduce(local_error,error_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,g%comm,ierr)
  call mpi_allreduce(samples,near_wall,size(samples),MPI_DOUBLE_PRECISION,MPI_SUM,g%comm,ierr)
  call mpi_allreduce(have,sample_count,size(have),MPI_INTEGER,MPI_SUM,g%comm,ierr)
  call require(all(sample_count==1),'wall slip fixture samples are missing or duplicated')
  slip_error = maxval(abs((35*wall(1:,:)-104*near_wall(:,:,1)+114*near_wall(:,:,2) &
                         -56*near_wall(:,:,3)+11*near_wall(:,:,4))/(12*h*h)))
  if (g%rank==0) print '(a,1x,es23.15)', 'independent-island-slip',slip_error
  if (g%rank==0) print '(a,3(1x,es23.15))', 'island-errors', error_max, &
                  max(maxval(abs(wall(0,:)-outer)),maxval(abs(wall(1:,:)-transpose(expected)))),backward_error

  ! Reuse the same component factor at independent circulation targets.
  wrong = gamma
  wrong(1,1) = wrong(1,1)+0.01D0
  call invert_weak_islands(cache,g%comm,q,wrong,mass,psi,shifted,backward_error,error)
  call require(len_trim(error)==0,error)
  call require(maxval(abs(shifted-wall))>1e-5,'changing one island target did not change the constants')
  if (g%rank==0) print '(a)', 'Weak island constraints and reused influence factor PASS'
  call release_weak_islands(cache)
  call deallocate(s)
  call deallocate(d)
  call parallel_finalise()

contains
  subroutine exact(x,y,value,lap)
    real(real_8), intent(in) :: x,y
    real(real_8), intent(out) :: value(2),lap(2)
    real(real_8) :: r,t,cut,cut_lap,shape,d1,d2,normal_poly,normal_derivative,normal_second
    real(real_8) :: tangent_poly,tangent_first,tangent_second,dx,dy
    integer :: k

    value = outer
    lap = 0.0D0
    do k = 1,2
      dx = x-centres(k)
      dy = y-edge
      r = hypot(dx,y-y0)
      cut = 0.0D0
      cut_lap = 0.0D0
      if (r <= 0.1D0) then
        cut = 1.0D0
      else if (r < 0.19D0) then
        t = (r-0.1D0)/0.09D0
        cut = 1-10*t**3+15*t**4-6*t**5
        d1 = (-30*t**2+60*t**3-30*t**4)/0.09D0
        d2 = (-60*t+180*t**2-120*t**3)/0.09D0**2
        cut_lap = d2+d1/r
      end if
      value = value+(expected(:,k)-outer)*cut
      lap = lap+(expected(:,k)-outer)*cut_lap
      call bump(dx,tangent_width,tangent_poly,tangent_first,tangent_second)
      call bump(dy,normal_width,normal_poly,normal_derivative,normal_second)
      shape = dy*normal_poly*tangent_poly
      value = value+amplitude(:,k)*shape
      lap = lap+amplitude(:,k)*(dy*normal_poly*tangent_second+(2*normal_derivative+dy*normal_second)*tangent_poly)
    end do
  end subroutine exact

  subroutine bump(x,width,value,derivative,second)
    real(real_8), intent(in) :: x,width
    real(real_8), intent(out) :: value,derivative,second
    real(real_8) :: t,u

    t = x/width
    u = max(0.0D0,1-t*t)
    value = u**5
    derivative = -10*t*u**4/width
    second = (80*t*t*u**3-10*u**4)/width**2
  end subroutine bump

  subroutine mean_pv(point,mask,value,mean_psi)
    real(real_8), intent(in) :: point(2)
    integer, intent(in) :: mask
    real(real_8), intent(out) :: value(2),mean_psi(2)
    real(real_8) :: p(2),phi(2),vort(2),a(2),b(2),weight
    integer :: r,i,j

    value = 0.0D0
    mean_psi = 0.0D0
    do r = 0,7
      if (.not.btest(mask,r)) cycle
      a = rays(:,r)
      b = rays(:,modulo(r+1,8))
      do j = 1,6
        do i = 1,6
          p = point+0.5D0*h*qpoint(i)*(a+qpoint(j)*(b-a))
          call exact(p(1),p(2),phi,vort)
          weight = 2*qpoint(i)*qweight(i)*qweight(j)
          value = value+weight*(h*h*vort+matmul(s%s_mat,phi))
          mean_psi = mean_psi+weight*phi
        end do
      end do
    end do
    value = value/popcnt(mask)
    mean_psi = mean_psi/popcnt(mask)
  end subroutine mean_pv

  subroutine require(ok,message)
    logical, intent(in) :: ok
    character(*), intent(in) :: message
    if (.not.ok) call pabort(message,__FILE__,__LINE__)
  end subroutine require
end program test_irregular_island_constraints
