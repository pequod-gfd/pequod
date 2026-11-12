program test_irregular_operators
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
  use qg_geometry, only: solid_node, fluid_node
  use qg_irregular_parallel
  use qg_irregular_operators
  use qg_parallel, only: parallel_init, parallel_finalise, pabort
  use qg_parameters, only: real_8
  implicit none
  include 'mpif.h'

  integer, parameter :: n = 33
  type(distributed_geometry) :: geom
  type(viscous_face_plan) :: viscous
  type(wall_trace_plan) :: trace
  integer :: rank,ranks,comm,ierr,grid(2),row,column,box(4),i,j,gi,gj,f,k,component,next,end_node,prev,side
  integer :: endpoints(2,2), ij(2)
  real(real_8) :: nan, x,y,step(2), normal(2), midpoint(2), points(2,2), value(2), defect, maximum
  real(real_8) :: integrals(2), total(2), strict_defect, trace_max, vertex(2), boundary(2,2), sample(2)
  real(real_8), allocatable :: reference(:,:)
  real(real_8), allocatable :: nodal(:,:,:), open_flux(:,:), wall_flux(:,:,:), tendency(:,:,:), q(:,:,:), gamma(:,:)
  character(1024) :: path,error

  if (command_argument_count() /= 1) then
    print '(a)', 'Irregular operator fixture requires a .dat path'
    stop 2
  end if
  call get_command_argument(1,path)
  call parallel_init()
  call mpi_comm_rank(MPI_COMM_WORLD,rank,ierr)
  call mpi_comm_size(MPI_COMM_WORLD,ranks,ierr)
  call mpi_comm_split(MPI_COMM_WORLD,0,ranks-1-rank,comm,ierr)
  grid = 0
  call mpi_dims_create(ranks,2,grid,ierr)
  row = rank/grid(2)
  column = mod(rank,grid(2))
  box = [row*n/grid(1),column*n/grid(2), (row+1)*n/grid(1)-row*n/grid(1), (column+1)*n/grid(2)-column*n/grid(2)]
  call load_distributed_geometry(path,n,n,box,comm,2.5D0,geom,error)
  call require(len_trim(error) == 0,error)
  allocate(nodal(geom%ni+2,geom%nj+2,2),open_flux(2,size(geom%faces)),wall_flux(geom%ni,geom%nj,2))
  allocate(tendency(geom%ni,geom%nj,2),q(geom%ni,geom%nj,2),gamma(0:ubound(geom%global%loops,1),2))
  allocate(reference(0:ubound(geom%global%loops,1),2))

  nan = ieee_value(0.0D0,ieee_quiet_nan)
  nodal = nan
  do j = 1,geom%nj+2
    gj = geom%j0+j-1
    do i = 1,geom%ni+2
      gi = geom%i0+i-1
      if (min(gi,gj) < 1 .or. max(gi,gj) > n) cycle
      if (geom%global%kind(gi,gj) == solid_node) cycle
      x = real(gj-1,real_8)
      y = real(n-gi,real_8)
      nodal(i,j,:) = [3*x-2*y+0.5D0,-x+4*y-0.25D0]
    end do
  end do
  call stiffness_flux(geom,nodal,open_flux,error)
  call require(len_trim(error) == 0,error)
  defect = 0.0D0
  do f = 1,size(geom%faces)
    step = (geom%faces(f)%endpoints(:,2)-geom%faces(f)%endpoints(:,1))/2.0D0
    normal = [step(2),-step(1)]
    value = [3*normal(1)-2*normal(2),-normal(1)+4*normal(2)]
    call require(all(ieee_is_finite(open_flux(:,f))),'stiffness face read poisoned land')
    defect = max(defect,maxval(abs(open_flux(:,f)-value)))
  end do
  call require(defect == 0.0D0,'affine Poisson stiffness flux is not exact')
  q = 0.0D0
  call weak_circulation(geom,q,open_flux,gamma,strict_defect,tendency,error)
  call require(len_trim(error) == 0,error)
  call require(maxval(abs(gamma)) == 0.0D0 .and. strict_defect == 0.0D0,'affine circulation balance')

  ! Independently integrate the actual outward flux of two quadratic omega
  ! fields through every open and closed segment. Only owners write values.
  open_flux = nan
  do f = 1,size(geom%faces)
    if (geom%faces(f)%owner /= geom%rank) cycle
    points = geom%faces(f)%endpoints/2.0D0
    open_flux(:,f) = exact_flux(points)
  end do
  call share_face_flux(geom,open_flux,error)
  call require(len_trim(error) == 0,error)
  wall_flux = nan
  where (spread(geom%global%component(geom%i0+1:geom%i0+geom%ni,geom%j0+1:geom%j0+geom%nj) >= 0,3,2)) wall_flux = 0.0D0
  reference = 0.0D0
  do component = 0,ubound(geom%global%loops,1)
    associate(loop => geom%global%loops(component)%ij)
      do k = 1,size(loop,2)
        next = modulo(k,size(loop,2))+1
        endpoints(:,1) = [loop(2,k)-1,n-loop(1,k)]
        endpoints(:,2) = [loop(2,next)-1,n-loop(1,next)]
        points = real(endpoints,real_8)
        reference(component,:) = reference(component,:)+exact_flux(points)
        midpoint = (points(:,1)+points(:,2))/2
        do end_node = 1,2
          ij = loop(:,merge(k,next,end_node==1))
          if (geom%owner(ij(1),ij(2)) /= geom%rank) cycle
          i = ij(1)-geom%i0
          j = ij(2)-geom%j0
          if (end_node == 1) then
            points(:,2) = midpoint
          else
            points(:,1) = midpoint
            points(:,2) = endpoints(:,2)
          end if
          wall_flux(i,j,:) = wall_flux(i,j,:)+exact_flux(points)
        end do
      end do
    end associate
  end do
  call diffusion_balance(geom,open_flux,wall_flux,tendency,gamma,error)
  call require(len_trim(error) == 0,error)
  defect = 0.0D0
  integrals = 0.0D0
  do j = 1,geom%nj
    do i = 1,geom%ni
      if (geom%area(i,j) == 0) cycle
      call require(all(ieee_is_finite(tendency(i,j,:))),'open or wall source poisoned a wet cell')
      defect = max(defect,maxval(abs(tendency(i,j,:)-[6.0D0,-2.0D0])))
      integrals = integrals+tendency(i,j,:)*geom%area(i,j)/geom%h**2
    end do
  end do
  call mpi_allreduce(defect,maximum,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
  call mpi_allreduce(integrals,total,2,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  call require(maximum < 1e-13,'cut-cell diffusion of exact quadratic flux is incorrect')
  call require(maxval(abs(gamma-reference)) < 1e-12,'component torque disagrees with contour integral')
  call require(maxval(abs(sum(gamma,dim=1)-total)) < 1e-12,'diffusion global budget failed')

  ! Now reconstruct cubic face fluxes from independent cell means.
  ! No point value of omega is passed to the numerical face plan.
  q = nan
  do j = 1,geom%nj
    gj = geom%j0+j
    do i = 1,geom%ni
      gi = geom%i0+i
      if (geom%area(i,j) == 0) cycle
      q(i,j,:) = cubic_cell_mean([real(gj-1,real_8),real(n-gi,real_8)],geom%global%wedge(gi,gj))
    end do
  end do
  call setup_viscous_faces(geom,viscous,error)
  call require(len_trim(error) == 0,error)
  call viscous_open_flux(geom,viscous,q,open_flux,error)
  call require(len_trim(error) == 0,error)
  defect = 0.0D0
  do f = 1,size(geom%faces)
    points = geom%faces(f)%endpoints/2.0D0
    call require(all(ieee_is_finite(open_flux(:,f))),'viscous reconstruction read poisoned land')
    defect = max(defect,maxval(abs(open_flux(:,f)-cubic_face_integral(points))))
  end do
  call mpi_allreduce(defect,maximum,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
  call require(maximum < 2e-10,'cubic shared viscous face flux')

  call setup_wall_trace(geom,trace,error)
  call require(len_trim(error) == 0,error)
  call reconstruct_wall_trace(trace,q,wall_flux)
  defect = 0.0D0
  do component = 0,ubound(geom%global%loops,1)
    associate(loop => geom%global%loops(component)%ij)
      do k = 1,size(loop,2)
        i = loop(1,k)-geom%i0
        j = loop(2,k)-geom%j0
        if (min(i,j)<1 .or. i>geom%ni .or. j>geom%nj) cycle
        prev = modulo(k-2,size(loop,2))+1
        next = mod(k,size(loop,2))+1
        vertex = [real(loop(2,k)-1,real_8),real(n-loop(1,k),real_8)]
        boundary(:,1) = [real(loop(2,prev)-1,real_8),real(n-loop(1,prev),real_8)]
        boundary(:,2) = [real(loop(2,next)-1,real_8),real(n-loop(1,next),real_8)]
        value = 0.0D0
        do side = 1,2
          step = 0.5D0*(boundary(:,side)-vertex)
          do f = -1,1,2
            sample = vertex+(0.5D0+f/(2*sqrt(3.0D0)))*step
            x = sample(1)
            y = sample(2)
            value = value+0.5D0*norm2(step)*[x*x+x*y+2*y*y+x**3/32+x*y*y/64, &
                                             -2*x*x+x*y+y*y+y**3/32+x*x*y/64]
          end do
        end do
        defect = max(defect,maxval(abs(value-wall_flux(i,j,:))))
      end do
    end associate
  end do
  call mpi_allreduce(defect,trace_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
  call require(trace_max < 2e-10,'cubic wall trace from cut means failed half-segment integration')
  if (rank == 0) print '(a,es14.6,1x,a)', 'Cut wall half-segment trace, max cubic defect:', trace_max,'PASS'

  ! Use quadratic nodal data to check each strict five-point row. Its exact
  ! Laplacian is independent of the implementation of the stiffness sum.
  do j = 1,geom%nj+2
    gj = geom%j0+j-1
    do i = 1,geom%ni+2
      gi = geom%i0+i-1
      if (min(gi,gj) < 1 .or. max(gi,gj) > n) cycle
      if (geom%global%kind(gi,gj) == solid_node) cycle
      x = real(gj-1,real_8)
      y = real(n-gi,real_8)
      nodal(i,j,:) = [x*x+x*y+2*y*y,-2*x*x+x*y+y*y]
    end do
  end do
  call stiffness_flux(geom,nodal,open_flux,error)
  call require(len_trim(error) == 0,error)
  q(:,:,1) = 6.0D0
  q(:,:,2) = -2.0D0
  call weak_circulation(geom,q,open_flux,gamma,strict_defect,tendency,error)
  call require(len_trim(error) == 0,error)
  call require(strict_defect == 0.0D0,'weak Poisson does not match the strict five-point stencil')
  call require(maxval(abs(sum(gamma,dim=1)-total)) < 1e-12,'weak Poisson global Stokes sum')
  call mean_to_pv_nodes(geom,q,tendency)
  do j = 1,geom%nj
    do i = 1,geom%ni
      if (geom%global%kind(geom%i0+i,geom%j0+j) == fluid_node) then
        call require(all(tendency(i,j,:) == q(i,j,:)),'mean-to-node full cell')
      else
        call require(all(tendency(i,j,:) == 0.0D0),'wall/solid mean leaked into an inversion row')
      end if
    end do
  end do
  if (rank == 0) print '(a,es14.6,1x,a)', 'Cut diffusion and shared Stokes/torque, max cubic face defect:', maximum,'PASS'
  call mpi_comm_free(comm,ierr)
  call parallel_finalise()

contains
  pure function exact_flux(p) result(result_flux)
    real(real_8), intent(in) :: p(2,2)
    real(real_8) :: result_flux(2), t(2), m(2), normal(2)
    t = p(:,2)-p(:,1)
    m = (p(:,1)+p(:,2))/2
    normal = [t(2),-t(1)]
    result_flux(1) = dot_product([2*m(1)+m(2),m(1)+4*m(2)],normal)
    result_flux(2) = dot_product([-4*m(1)+m(2),m(1)+2*m(2)],normal)
  end function exact_flux

  pure function cubic_cell_mean(origin,mask) result(mean)
    real(real_8), intent(in) :: origin(2)
    integer, intent(in) :: mask
    real(real_8), parameter :: nodes(3) = [0.1127016653792583D0,0.5D0,0.8872983346207417D0]
    real(real_8), parameter :: weights(3) = [5,8,5]/18.0D0
    real(real_8), parameter :: rays(2,0:7) = 0.5D0*reshape([1,0,1,1,0,1,-1,1,-1,0,-1,-1,0,-1,1,-1],[2,8])
    real(real_8) :: mean(2), p(2), a(2), b(2), x, y
    integer :: ray, i, j

    mean = 0.0D0
    do ray = 0,7
      if (.not. btest(mask,ray)) cycle
      a = rays(:,ray)
      b = rays(:,modulo(ray+1,8))
      do i = 1,3
        do j = 1,3
          p = origin+nodes(i)*((1-nodes(j))*a+nodes(j)*b)
          x = p(1)
          y = p(2)
          mean = mean+2*nodes(i)*weights(i)*weights(j)*[x*x+x*y+2*y*y+x**3/32+x*y*y/64, &
                                                       -2*x*x+x*y+y*y+y**3/32+x*x*y/64]
        end do
      end do
    end do
    mean = mean/popcnt(mask)
  end function cubic_cell_mean

  pure function cubic_face_integral(p) result(result_flux)
    real(real_8), intent(in) :: p(2,2)
    real(real_8) :: result_flux(2), t(2), normal(2), x, y, q(2)
    integer :: i

    t = p(:,2)-p(:,1)
    normal = [t(2),-t(1)]
    result_flux = exact_flux(p)
    do i = -1,1,2
      q = (p(:,1)+p(:,2))/2+i*t/(2*sqrt(3.0D0))
      x = q(1)
      y = q(2)
      result_flux(1) = result_flux(1)+0.5D0*dot_product([3*x*x/32+y*y/64,x*y/32],normal)
      result_flux(2) = result_flux(2)+0.5D0*dot_product([x*y/32,3*y*y/32+x*x/64],normal)
    end do
  end function cubic_face_integral

  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message
    if (.not. condition) call pabort(message,__FILE__,__LINE__)
  end subroutine require
end program test_irregular_operators
