#include "config.h"
program test_irregular_mpi_geometry
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
  use qg_geometry, only: fluid_node, wall_node
  use qg_irregular_parallel
  use qg_parallel, only: parallel_init, parallel_finalise, pabort
  use qg_parameters, only: real_8
  implicit none
  include 'mpif.h'

  integer, parameter :: n = 33
  type(distributed_geometry) :: geom, rejected
  type(donor_plan) :: donors
  real(real_8), allocatable :: field(:,:,:), samples(:,:)
  integer, allocatable :: points(:,:), owner_face(:,:), face_records(:,:), incidence(:), total_incidence(:)
  integer :: world_rank, ranks, comm, ierr, row, column, grid(2), box(4), k, i,j,gi,gj,side, f, id
  integer :: req, di, dj, ti, tj, remote, remote_total, dry, dry_total, walls, total_walls
  integer :: v(2), count_faces, mode, status, bad_point(2,1)
  real(real_8) :: area, perimeter, totals(2), expected_perimeter, x, y, dx, dy
  real(real_8) :: raw_moments(5), moment_totals(5), expected_moments(5), wall_integral, global_wall_integral
  character(1024) :: path, error
  character(30) :: fixture

  ! The external fixture runner passes the geometry path and runs 1/2/4/12 ranks.
  if (command_argument_count() == 0) then
    print '(a)', 'Geometry MPI fixture requires a .dat path and case name'
    stop
  end if
  call get_command_argument(1,path)
  call get_command_argument(2,fixture)
  call parallel_init()
  call mpi_comm_rank(MPI_COMM_WORLD,world_rank,ierr)
  call mpi_comm_size(MPI_COMM_WORLD,ranks,ierr)
  grid = 0
  call mpi_dims_create(ranks,2,grid,ierr)
  row = world_rank/grid(2)
  column = mod(world_rank,grid(2))
  box = [row*n/grid(1),column*n/grid(2), (row+1)*n/grid(1)-row*n/grid(1), &
         (column+1)*n/grid(2)-column*n/grid(2)]
  ! Reverse rank numbers while keeping each rank's original physical patch.
  call mpi_comm_split(MPI_COMM_WORLD,0,ranks-1-world_rank,comm,ierr)
  call load_distributed_geometry(trim(path)//'.missing',n,n,box,comm,1.0D0,rejected,error)
  if (len_trim(error) == 0) call pabort('missing geometry accepted',__FILE__,__LINE__)
  call load_distributed_geometry(trim(path),n,n,box,comm,1.0D0,geom,error)
  if (len_trim(error) /= 0) call pabort(error,__FILE__,__LINE__)
  if (geom%rank /= ranks-1-world_rank) call pabort('rank order was not reversed',__FILE__,__LINE__)
  if (any(geom%owner(box(1)+1:box(1)+box(3),box(2)+1:box(2)+box(4)) /= geom%rank)) &
    call pabort('partition ownership is incorrect',__FILE__,__LINE__)

  area = sum(geom%area)
  perimeter = 0.0D0
  wall_integral = 0.0D0
  walls = size(geom%walls)
  do k = 1,walls
    perimeter = perimeter+geom%walls(k)%length
    v = geom%walls(k)%endpoints(:,2)-geom%walls(k)%endpoints(:,1)
    if (maxval(abs(v)) /= 2) call pabort('invalid contour segment',__FILE__,__LINE__)
    if (abs(sum(geom%walls(k)%outward**2)-1.0D0) > 4e-15) call pabort('invalid wall normal',__FILE__,__LINE__)
    wall_integral = wall_integral + geom%walls(k)%length * &
      dot_product(0.25D0*sum(geom%walls(k)%endpoints,dim=2),geom%walls(k)%outward)
  end do
  call mpi_allreduce([area,perimeter],totals,2,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  call mpi_allreduce(walls,total_walls,1,MPI_INTEGER,MPI_SUM,comm,ierr)
  select case(fixture)
    case('triangle')
      expected_perimeter = 64 + 32*sqrt(2.0D0)
      expected_moments = [32.0D0/3,32.0D0/3,1024.0D0/6,1024.0D0/12,1024.0D0/6]
      if (totals(1) /= 512 .or. total_walls /= 96) call pabort('wrong triangle area or wall count',__FILE__,__LINE__)
    case('dry')
      expected_perimeter = 16*sqrt(2.0D0)
      expected_moments = [6.0D0,6.0D0,36.0D0+16.0D0/6,36.0D0,36.0D0+16.0D0/6]
      if (totals(1) /= 32 .or. total_walls /= 16) call pabort('wrong diamond area or wall count',__FILE__,__LINE__)
    case default
      call pabort('unknown test fixture',__FILE__,__LINE__)
  end select
  if (abs(totals(2)-expected_perimeter) > 2e-13) call pabort('wrong contour length',__FILE__,__LINE__)
  call mpi_allreduce(wall_integral,global_wall_integral,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  if (abs(global_wall_integral-2*totals(1)) > 1e-12) call pabort('contour normal integral mismatch',__FILE__,__LINE__)
  raw_moments = 0.0D0
  do j = 1,geom%nj
    x = real(geom%j0+j-1,real_8)
    do i = 1,geom%ni
      y = real(n-geom%i0-i,real_8)
      dx = geom%moments(2,i,j)
      dy = geom%moments(3,i,j)
      raw_moments = raw_moments + geom%area(i,j) * &
        [x+dx,y+dy,x*x+2*x*dx+geom%moments(4,i,j), &
         x*y+x*dy+y*dx+geom%moments(5,i,j),y*y+2*y*dy+geom%moments(6,i,j)]
    end do
  end do
  call mpi_allreduce(raw_moments,moment_totals,5,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  if (maxval(abs(moment_totals/totals(1)-expected_moments)) > 2e-12) &
    call pabort('cut-cell polynomial moments mismatch',__FILE__,__LINE__)
  dry = merge(1,0,area == 0.0D0)
  call mpi_allreduce(dry,dry_total,1,MPI_INTEGER,MPI_SUM,comm,ierr)
  if (ranks > 1 .and. fixture == 'dry' .and. dry_total == 0) call pabort('fixture must have dry ranks',__FILE__,__LINE__)

  allocate(owner_face(5,2*n*n),face_records(5,2*n*n),incidence(2*n*n),total_incidence(2*n*n))
  owner_face = 0
  incidence = 0
  do k = 1,size(geom%faces)
    associate(face => geom%faces(k))
      if (face%owner /= geom%rank) cycle
      if (owner_face(1,face%id) /= 0) call pabort('duplicate local face owner',__FILE__,__LINE__)
      owner_face(:,face%id) = [1,reshape(face%endpoints,[4])]
    end associate
  end do
  call mpi_allreduce(owner_face,face_records,size(owner_face),MPI_INTEGER,MPI_SUM,comm,ierr)
  do k = 1,size(geom%faces)
    associate(face => geom%faces(k))
      if (any(face_records(:,face%id) /= [1,reshape(face%endpoints,[4])])) &
        call pabort('face records disagree across ranks',__FILE__,__LINE__)
    end associate
  end do
  do j = 1,geom%nj
    do i = 1,geom%ni
      do side = 1,4
        f = geom%face_index(side,i,j)
        if (f == 0) cycle
        id = geom%faces(f)%id
        incidence(id) = incidence(id)+1
      end do
    end do
  end do
  call mpi_allreduce(incidence,total_incidence,size(incidence),MPI_INTEGER,MPI_SUM,comm,ierr)
  if (any(total_incidence /= 2*face_records(1,:))) call pabort('shared-face incidence mismatch',__FILE__,__LINE__)
  count_faces = sum(face_records(1,:))

  allocate(field(0:geom%ni+1,0:geom%nj+1,2))
  field = ieee_value(0.0D0,ieee_quiet_nan)
  do j = 1,geom%nj
    gj = geom%j0+j
    do i = 1,geom%ni
      gi = geom%i0+i
      if (geom%global%kind(gi,gj) /= fluid_node) cycle
      field(i,j,:) = [encoded(gi,gj,1),encoded(gi,gj,2)]
      if (abs(geom%moments(4,i,j)-1.0D0/12) > 1e-15) call pabort('full cell moment mismatch',__FILE__,__LINE__)
    end do
  end do
  ! Include remote donors at Chebyshev distance two and repeated donor IDs.
  allocate(points(2,24*walls))
  req = 0
  remote = 0
  do k = 1,walls
    gi = n-geom%walls(k)%endpoints(2,1)/2
    gj = 1+geom%walls(k)%endpoints(1,1)/2
    if (geom%global%kind(gi,gj) /= wall_node) call pabort('wall owner is not a wall node',__FILE__,__LINE__)
    do dj = -2,2
      do di = -2,2
        if (di == 0 .and. dj == 0) cycle
        ti = gi+di
        tj = gj+dj
        if (min(ti,tj) < 1 .or. max(ti,tj) > n) cycle
        if (geom%global%kind(ti,tj) /= fluid_node) cycle
        req = req+1
        points(:,req) = [ti,tj]
        if (ti < geom%i0 .or. ti > geom%i0+geom%ni+1 .or. tj < geom%j0 .or. tj > geom%j0+geom%nj+1) remote = remote+1
      end do
    end do
  end do
  call mpi_allreduce(remote,remote_total,1,MPI_INTEGER,MPI_SUM,comm,ierr)
  if (ranks > 1 .and. fixture == 'triangle' .and. remote_total == 0) &
    call pabort('fixture does not exercise donors beyond one halo ring',__FILE__,__LINE__)
  call build_donor_plan(geom,points(:,:req),donors,error)
  if (len_trim(error) /= 0) call pabort(error,__FILE__,__LINE__)
  allocate(samples(2,req))
  do mode = 1,2
    call gather_donors(donors,field(1:geom%ni,1:geom%nj,:),samples)
    do k = 1,req
      if (.not. all(ieee_is_finite(samples(:,k)))) call pabort('poisoned donor data',__FILE__,__LINE__)
      if (any(samples(:,k) /= mode*[encoded(points(1,k),points(2,k),1),encoded(points(1,k),points(2,k),2)])) &
        call pabort('donor values or result order incorrect',__FILE__,__LINE__)
    end do
    field = 2*field
  end do
  status = 0
  if (geom%rank == 0) status = 1
  bad_point = 1
  call build_donor_plan(geom,bad_point(:,:status),donors,error)
  if (len_trim(error) == 0) call pabort('solid or wall donor accepted',__FILE__,__LINE__)

  if (world_rank == 0) then
    print '(a,1x,f12.6,1x,f12.6,2(1x,i0))', geom%global%digest, totals, total_walls, count_faces
    print '(a,3(1x,i0),1x,a)', 'ranks, dry ranks and remote donors:', ranks, dry_total, remote_total, 'PASS'
  end if
  call mpi_comm_free(comm,ierr)
  call parallel_finalise()

contains

  real(real_8) function encoded(i,j,layer)
    integer, intent(in) :: i,j,layer
    encoded = 10000*layer+100*i+j
  end function encoded

end program test_irregular_mpi_geometry
