program test_irregular_smooth_corners
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
  use qg_diagnostics, only: BC_FREE_SLIP, BC_PARTIAL_SLIP, BC_NO_SLIP
  use qg_geometry, only: fluid_node, wall_node
  use qg_irregular_parallel, only: distributed_geometry, load_distributed_geometry
  use qg_parallel, only: parallel_init, parallel_finalise, pabort
  use qg_smooth_corners
  use qg_parameters, only: real_8
  implicit none
  include 'mpif.h'

  real(real_8), parameter :: rays(2,0:7) = reshape([1,0,1,1,0,1,-1,1,-1,0,-1,-1,0,-1,1,-1],[2,8])
  real(real_8), parameter :: nodes(3) = [0.1127016653792583D0,0.5D0,0.8872983346207417D0]
  real(real_8), parameter :: quadrature(3) = [5,8,5]/18.0D0, pi = acos(-1.0D0)
  integer, parameter :: n = 33
  real(real_8), parameter :: c0 = 0.375D0, ell = 0.7D0
  type(distributed_geometry) :: geom
  type(smooth_corner_plan) :: plan
  real(real_8) :: points(2,81), weights(81,5), data(81), oracle(5), calculated(5), derivative(2)
  real(real_8) :: angle, theta, length, tangent(2,2), inward(2,2), exponent_vector(2), p(2), h
  real(real_8) :: value, omega, error_max, error_local, norm, nan, distance, closest
  real(real_8), allocatable :: psi(:,:,:), diagnostics(:,:,:,:), wall(:,:)
  integer :: sectors, start, mode, mask, count, power, i,j,k,l,r,a,b,last,rank,ranks,comm,grid(2),row,column,box(4)
  integer :: gi,gj,corner(2),candidate(2),ierr,ncorners,total_corners
  logical :: accepted
  character(1024) :: arg,error

  call get_command_argument(1,arg)
  if (command_argument_count() /= 1) then
    print '(a)', 'Smooth corner fixture requires a .dat path or --scalar'
    stop 2
  end if
  error_max = 0.0D0
  if (arg == '--scalar') then
    do sectors = 1,7
      do start = 0,7
        last = modulo(start+sectors,8)
        mask = 0
        do r = 0,sectors-1
          mask = ibset(mask,modulo(start+r,8))
        end do
        tangent = rays(:,[start,last])
        inward(:,1) = [-tangent(2,1),tangent(1,1)]/norm2(tangent(:,1))
        inward(:,2) = [tangent(2,2),-tangent(1,2)]/norm2(tangent(:,2))
        norm = inward(1,1)*inward(2,2)-inward(1,2)*inward(2,1)
        if (sectors == 4) then
          exponent_vector = inward(:,1)/(2*ell)
        else
          exponent_vector = [inward(2,2)-inward(2,1),inward(1,1)-inward(1,2)]/(2*ell*norm)
        end if
        power = 1
        do while (mod(power*sectors,4) /= 0)
          power = power+1
        end do
        theta = start*pi/4
        count = 0
        do j = -4,4
          do i = -4,4
            if (i == 0 .and. j == 0) cycle
            angle = modulo(atan2(real(j,real_8),real(i,real_8))-theta,2*pi)
            if (angle <= 1e-12 .or. angle >= sectors*pi/4-1e-12) cycle
            count = count+1
            points(:,count) = [i,j]
          end do
        end do
        do mode = 1,3
          call smooth_corner_weights(mask,mode,ell,points(:,:count),weights(:count,:),accepted)
          if (.not. accepted) stop 1
          do i = 1,count
            call exact(points(:,i),value,derivative,omega)
            data(i) = c0+value
          end do
          calculated = matmul(data(:count)-c0,weights(:count,:))
          oracle = 0.0D0
          call exact([0.0D0,0.0D0],value,derivative,oracle(corner_omega_point))
          do r = 0,sectors-1
            a = modulo(start+r,8)
            b = modulo(a+1,8)
            do k = 1,3
              do l = 1,3
                p = 0.5D0*nodes(k)*(rays(:,a)+nodes(l)*(rays(:,b)-rays(:,a)))
                call exact(p,value,derivative,omega)
                h = 2*nodes(k)*quadrature(k)*quadrature(l)/sectors
                oracle(corner_psi_mean) = oracle(corner_psi_mean)+h*value
                oracle(corner_omega_mean) = oracle(corner_omega_mean)+h*omega
              end do
            end do
          end do
          do k = 1,2
            length = 0.5D0*norm2(tangent(:,k))
            do l = 1,3
              p = 0.5D0*nodes(l)*tangent(:,k)
              call exact(p,value,derivative,omega)
              oracle(3+k) = oracle(3+k)-length*quadrature(l)*dot_product(inward(:,k),derivative)
            end do
          end do
          error_max = max(error_max,maxval(abs(calculated-oracle)))
        end do
      end do
    end do
    if (error_max > 2e-11) stop 1
    print '(a,es14.6,1x,a)', '168 oriented smooth corner polynomial checks, max error:', error_max,'PASS'
    stop
  end if

  call parallel_init()
  call mpi_comm_rank(MPI_COMM_WORLD,rank,ierr)
  call mpi_comm_size(MPI_COMM_WORLD,ranks,ierr)
  call mpi_comm_split(MPI_COMM_WORLD,0,ranks-1-rank,comm,ierr)
  grid = 0
  call mpi_dims_create(ranks,2,grid,ierr)
  row = rank/grid(2)
  column = mod(rank,grid(2))
  box = [row*n/grid(1),column*n/grid(2), (row+1)*n/grid(1)-row*n/grid(1), (column+1)*n/grid(2)-column*n/grid(2)]
  call load_distributed_geometry(arg,n,n,box,comm,1.0D0,geom,error)
  call require(len_trim(error)==0,error)
  call setup_smooth_corners(geom,BC_NO_SLIP,0.0D0,plan,error)
  call require(len_trim(error)==0,error)
  allocate(psi(geom%ni,geom%nj,1),diagnostics(geom%ni,geom%nj,1,5),wall(0:0,1))
  nan = ieee_value(0.0D0,ieee_quiet_nan)
  psi = nan
  diagnostics = nan
  wall = c0
  do j = 1,geom%nj
    gj = geom%j0+j
    do i = 1,geom%ni
      gi = geom%i0+i
      if (geom%global%kind(gi,gj) /= fluid_node) cycle
      closest = huge(0.0D0)
      do k = 1,size(geom%global%loops(0)%ij,2)
        candidate = geom%global%loops(0)%ij(:,k)
        mask = geom%global%wedge(candidate(1),candidate(2))
        if (popcnt(mask) == 4) cycle
        distance = sum(([gi,gj]-candidate)**2)
        if (distance >= closest) cycle
        closest = distance
        corner = candidate
      end do
      psi(i,j,1) = c0+real(gi-corner(1),real_8)**2*real(gj-corner(2),real_8)**2
    end do
  end do
  call reconstruct_smooth_corners(plan,psi,wall,diagnostics)
  error_local = 0.0D0
  ncorners = 0
  do j = 1,geom%nj
    gj = geom%j0+j
    do i = 1,geom%ni
      gi = geom%i0+i
      if (geom%global%kind(gi,gj) /= wall_node .or. popcnt(geom%global%wedge(gi,gj)) == 4) then
        call require(.not. any(ieee_is_finite(diagnostics(i,j,1,:))),'non-corner output changed')
        cycle
      end if
      ncorners = ncorners+1
      call require(all(ieee_is_finite(diagnostics(i,j,1,:))),'corner reconstruction read poisoned data')
      oracle = [c0+1.0D0/144,0.0D0,1.0D0/3,0.0D0,0.0D0]
      error_local = max(error_local,maxval(abs(diagnostics(i,j,1,:)-oracle)))
    end do
  end do
  call mpi_allreduce(error_local,error_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
  call mpi_allreduce(ncorners,total_corners,1,MPI_INTEGER,MPI_SUM,comm,ierr)
  call require(total_corners == 6,'expected L-shaped geometry with six corners')
  call require(error_max < 2e-11,'MPI smooth corner polynomial accuracy')
  if (rank == 0) print '(a,es14.6,1x,a)', 'MPI smooth corner max error:', error_max,'PASS'
  call mpi_comm_free(comm,ierr)
  call parallel_finalise()

contains
  subroutine exact(q,value,gradient,laplacian)
    real(real_8), intent(in) :: q(2)
    real(real_8), intent(out) :: value,gradient(2),laplacian
    real(real_8) :: a,b,d, poly(0:4), s, product_value, product_gradient(2)
    complex(real_8) :: z, phase, first
    integer :: k

    if (mode == BC_NO_SLIP) then
      a = dot_product(q,inward(:,1))
      b = dot_product(q,inward(:,2))
      d = dot_product(inward(:,1),inward(:,2))
      value = a*a*b*b
      gradient = 2*a*b*(b*inward(:,1)+a*inward(:,2))
      laplacian = 2*(a*a+b*b)+8*a*b*d
    else
      phase = cmplx(cos(theta),-sin(theta),kind=real_8)
      z = cmplx(q(1),q(2),kind=real_8)*phase
      product_value = aimag(z**power)
      first = power*z**(power-1)*phase
      product_gradient = [aimag(first),real(first,real_8)]
      poly = 0.0D0
      poly(0) = 1.0D0
      if (mode == BC_PARTIAL_SLIP) then
        s = dot_product(q,exponent_vector)
        do k = 1,4-power
          poly(k) = poly(k-1)*s/k
        end do
      end if
      value = product_value*sum(poly)
      gradient = product_gradient*sum(poly)
      laplacian = 0.0D0
      if (mode == BC_PARTIAL_SLIP .and. power < 4) then
        gradient = gradient+product_value*exponent_vector*sum(poly(:3-power))
        laplacian = 2*dot_product(exponent_vector,product_gradient)*sum(poly(:3-power))
        if (power <= 2) laplacian = laplacian+sum(exponent_vector**2)*product_value*sum(poly(:2-power))
      end if
    end if
  end subroutine exact

  subroutine require(ok,message)
    logical, intent(in) :: ok
    character(*), intent(in) :: message
    if (.not. ok) call pabort(message,__FILE__,__LINE__)
  end subroutine require
end program test_irregular_smooth_corners
