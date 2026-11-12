program test_cape_cabaret
  use mpi
  use qg_parameters, only: real_8, pi
  use qg_geometry, only: write_rectangular_cape
  use qg_parallel, only: parallel_init, parallel_finalise, pabort
  use qg_irregular_parallel, only: distributed_geometry, load_distributed_geometry
  use qg_irregular_transport, only: cut_cabaret_plan, setup_cut_cabaret, cabaret_face_traces, cabaret_cell_velocity, &
                                    transport_fluxes, wet_cell_rate
  implicit none

  type(distributed_geometry) :: geom
  type(cut_cabaret_plan) :: plan
  real(real_8), parameter :: finish=0.18D0
  real(real_8), parameter :: abscissa(3)=[-0.7745966692414834D0,0.0D0,0.7745966692414834D0]
  real(real_8), parameter :: weights(3)=[5.0D0,8.0D0,5.0D0]/18.0D0
  real(real_8), allocatable :: q(:,:,:),initial(:,:,:),exact(:,:,:),half(:,:,:),source(:,:,:),rate(:,:,:), &
                               area(:,:,:),incoming(:,:,:),outgoing(:,:,:),s0(:,:,:),s1(:,:,:),s2(:,:,:), &
                               s3(:,:,:),s4(:,:,:),speed(:,:,:),trace(:,:,:,:),next(:,:,:,:),velocity(:,:,:,:),flux(:,:)
  logical, allocatable :: transition(:,:,:)
  real(real_8) :: dx,dt,t,base_dt,peak,cfl,local(5),global(5),endpoints(2,2),a,b,m,length
  integer :: n,steps,rank,ranks,ierr,grid(2),box(4),row,column,k,i,j,f,side,neighbor,failures
  logical :: smooth
  character(24) :: arg
  character(1024) :: error

  if (command_argument_count()/=2) then
    print '(a)', 'Cape CABARET fixture requires a grid size and affine or smooth'
    stop 2
  end if
  call get_command_argument(1,arg)
  read(arg,*) n
  call get_command_argument(2,arg)
  smooth=arg=='smooth'
  if ((.not. smooth .and. arg/='affine') .or. n<9 .or. mod(n-1,4)/=0) stop 2
  dx=1.0D0/(n-1)
  steps=2*(n-1)
  base_dt=finish/steps
  call parallel_init()
  call mpi_comm_rank(MPI_COMM_WORLD,rank,ierr)
  call mpi_comm_size(MPI_COMM_WORLD,ranks,ierr)
  if (rank==0) then
    call write_rectangular_cape('transport-cape.dat',n,n,1+(n-1)/4,(n-1)/4,(n-1)/2,dx,error)
    call check(error)
  end if
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  grid=0
  call mpi_dims_create(ranks,2,grid,ierr)
  row=rank/grid(2)
  column=mod(rank,grid(2))
  box=[row*n/grid(1),column*n/grid(2),(row+1)*n/grid(1)-row*n/grid(1),(column+1)*n/grid(2)-column*n/grid(2)]
  call load_distributed_geometry('transport-cape.dat',n,n,box,MPI_COMM_WORLD,dx,geom,error)
  call check(error)
  call setup_cut_cabaret(geom,plan,error,include_full=.true.)
  call check(error)
  allocate(q(box(3),box(4),1),transition(box(3),box(4),1))
  allocate(initial,exact,half,source,rate,area,incoming,outgoing,s0,s1,s2,s3,s4,mold=q)
  f=size(geom%faces)
  allocate(speed(2,1,f),trace(2,2,1,f),next(2,2,1,f),flux(1,f),velocity(2,box(3),box(4),1))
  source=0.0D0
  area(:,:,1)=geom%area/dx**2
  transition=area>0 .and. area<1
  do j=1,geom%nj
    do i=1,geom%ni
      do side=1,4
        f=geom%face_index(side,i,j)
        if (f==0) cycle
        neighbor=merge(2,1,side<=2)
        associate(ij=>geom%faces(f)%cells(:,neighbor))
          transition(i,j,1)=transition(i,j,1) .or. geom%global%wedge(ij(1),ij(2))/=255
        end associate
      end do
    end do
  end do
  do f=1,size(geom%faces)
    endpoints=0.5D0*dx*geom%faces(f)%endpoints
    a=streamfunction(endpoints(:,1))
    b=streamfunction(endpoints(:,2))
    m=streamfunction(0.5D0*sum(endpoints,dim=2))
    length=sum(abs(endpoints(:,2)-endpoints(:,1)))
    speed(:,1,f)=[3*(a-m)-(m-b),3*(m-b)-(a-m)]/(dx*length)
    do k=1,2
      trace(k,:,1,f)=scalar(endpoints(:,k),0.0D0)
    end do
  end do
  call cabaret_cell_velocity(geom,speed,velocity)
  t=0.0D0
  call sample_cells(t,initial,s0)
  q=initial
  peak=0.0D0
  do k=1,steps
    dt=base_dt*merge(0.75D0,1.25D0,mod(k,2)==1)
    call sample_cells(t+dt/4,exact,s1)
    call sample_cells(t+dt/2,exact,s2)
    call sample_cells(t+3*dt/4,exact,s3)
    call sample_cells(t+dt,exact,s4)
    incoming=-dt/12*(s0+4*s1+s2)
    outgoing=dt/12*(s2+4*s3+s4)
    call transport_fluxes(geom,speed,trace,flux,cfl,error)
    call check(error)
    peak=max(peak,dt*cfl)
    call wet_cell_rate(geom,flux,source,rate,error)
    call check(error)
    half=q+dt/2*rate-incoming
    call cabaret_face_traces(geom,plan,half,trace,speed,velocity,incoming,outgoing,dt,.true.,next,failures,error)
    call check(error)
    if (failures/=0) call check('rank-deficient stencil on axial cape')
    call transport_fluxes(geom,speed,next,flux,cfl,error)
    call check(error)
    call wet_cell_rate(geom,flux,source,rate,error)
    call check(error)
    q=half+dt/2*rate+outgoing
    trace=next
    initial=initial+outgoing-incoming
    s0=s4
    t=t+dt
  end do
  local=[sum(area*(q-exact)**2),sum(area),sum(area*(q-exact)**2,mask=transition), &
         sum(area,mask=transition),sum(area*(q-initial))]
  call mpi_allreduce(local,global,5,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  if (peak>0.45D0 .or. abs(global(5))>1.0D-11*global(2)) call check('cape transport CFL or integrated budget failed')
  if (.not. smooth .and. maxval(global([1,3]))>1.0D-23*global(2)) call check('affine cape trajectory lost accuracy')
  if (rank==0) print '(a,4(1x,es24.16),1x,a)', 'cape CABARET errors:', sqrt(global(1)/global(2)), &
                   sqrt(global(3)/global(4)),global(5),peak,'PASS'
  call parallel_finalise()

contains

  function streamfunction(point) result(value)
    real(real_8), intent(in) :: point(2)
    real(real_8) :: value,x,y

    x=point(1)
    y=point(2)
    value=80*x*(1-x)*(x-0.25D0)*(x-0.5D0)*y*(1-y)*(y-0.5D0)
  end function streamfunction

  function scalar(point,time) result(value)
    real(real_8), intent(in) :: point(2),time
    real(real_8) :: value

    value=0.3D0+0.5D0*point(1)-0.25D0*point(2)+0.7D0*time
    if (smooth) value=value+0.03D0*sin(2*pi*point(1)+0.7D0*time)*sin(pi*point(2)-0.2D0*time)
  end function scalar

  subroutine sample_cells(time,mean,forcing)
    real(real_8), intent(in) :: time
    real(real_8), intent(out) :: mean(:,:,:),forcing(:,:,:)
    integer, parameter :: sx(0:3)=[1,-1,-1,1],sy(0:3)=[1,1,-1,-1]
    real(real_8) :: point(2),centre(2),x,y,phase_x,phase_y,fx,g,df,dg,u,v,w
    integer :: ii,jj,mask,quadrant,p,r

    mean=0.0D0
    forcing=0.0D0
    do jj=1,geom%nj
      do ii=1,geom%ni
        if (area(ii,jj,1)==0) cycle
        centre=dx*[geom%j0+jj-1,n-geom%i0-ii]
        mask=geom%global%wedge(geom%i0+ii,geom%j0+jj)
        do quadrant=0,3
          if (.not. btest(mask,2*quadrant)) cycle
          do r=1,3
            do p=1,3
              point=centre+0.25D0*dx*[sx(quadrant)*(1+abscissa(p)),sy(quadrant)*(1+abscissa(r))]
              x=point(1)
              y=point(2)
              w=0.25D0*weights(p)*weights(r)/area(ii,jj,1)
              mean(ii,jj,1)=mean(ii,jj,1)+w*scalar(point,time)
              if (.not. smooth) cycle
              fx=x*(1-x)*(x-0.25D0)*(x-0.5D0)
              df=(1-2*x)*(x-0.25D0)*(x-0.5D0)+x*(1-x)*(2*x-0.75D0)
              g=y*(1-y)*(y-0.5D0)
              dg=(1-2*y)*(y-0.5D0)+y*(1-y)
              u=-80*fx*dg
              v=80*df*g
              phase_x=2*pi*x+0.7D0*time
              phase_y=pi*y-0.2D0*time
              forcing(ii,jj,1)=forcing(ii,jj,1)+w*(0.7D0+0.5D0*u-0.25D0*v+0.03D0*( &
                (0.7D0+2*pi*u)*cos(phase_x)*sin(phase_y)+(-0.2D0+pi*v)*sin(phase_x)*cos(phase_y)))
            end do
          end do
        end do
        if (.not. smooth) forcing(ii,jj,1)=0.7D0+dx*(0.5D0*velocity(1,ii,jj,1)-0.25D0*velocity(2,ii,jj,1))
      end do
    end do
  end subroutine sample_cells

  subroutine check(message)
    character(*), intent(in) :: message
    if (len_trim(message)/=0) call pabort(trim(message),__FILE__,__LINE__)
  end subroutine check
end program test_cape_cabaret
