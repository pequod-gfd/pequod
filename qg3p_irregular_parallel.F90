#include "config.h"

module qg_irregular_parallel
  use qg_geometry, only: irregular_geometry, read_geometry, fluid_node, solid_node
  use qg_moments, only: n_moments, triangle_moments
  use qg_parallel, only: pabort
  use qg_parameters, only: real_8
  implicit none
  include 'mpif.h'
  private
  public :: distributed_geometry, open_face, wall_segment, donor_plan
  public :: load_distributed_geometry, build_donor_plan, gather_donors
  public :: half_coordinates, cell_mean_moments, make_sector_moments, broadcast_error
  public :: share_face_flux, gather_second_face_trace

  integer, parameter :: face_di(4) = [1,0,-1,0], face_dj(4) = [0,1,0,-1]
  integer, parameter :: opposite(4) = [3,4,1,2]
  integer, parameter :: face_halves(2,4) = reshape([5,6, 7,0, 2,1, 4,3],[2,4])
  integer, parameter :: ray_x(0:7) = [1,1,0,-1,-1,-1,0,1]
  integer, parameter :: ray_y(0:7) = [0,1,1,1,0,-1,-1,-1]

  type open_face
    integer :: id, owner
    integer :: cells(2,2) ! global (i,j), with north/west cell first
    integer :: endpoints(2,2) ! half-grid (x,y), increasing x or y
  end type open_face

  type wall_segment
    integer :: component, edge
    integer :: endpoints(2,2) ! half-grid (x,y), fluid on the left
    real(real_8) :: length, outward(2)
  end type wall_segment

  type face_exchange
    integer, allocatable :: send_count(:), receive_count(:), send_offset(:), receive_offset(:)
    integer, allocatable :: send_face(:), receive_face(:)
  end type face_exchange

  type distributed_geometry
    type(irregular_geometry) :: global
    integer :: comm, rank, ranks ! communicator borrowed from the caller
    integer :: i0, j0, ni, nj
    real(real_8) :: h
    integer, allocatable :: owner(:,:) ! rank in comm, indexed by global (i,j)
    real(real_8), allocatable :: area(:,:), moments(:,:,:)
    type(open_face), allocatable :: faces(:)
    integer, allocatable :: face_index(:,:,:) ! south, east, north, west; 0 for closed
    type(wall_segment), allocatable :: walls(:) ! owned only
    type(face_exchange), private :: face_comm
  end type distributed_geometry

  type donor_plan
    private
    integer :: comm
    integer :: ni, nj
    integer, allocatable :: request_count(:), request_offset(:)
    integer, allocatable :: send_count(:), send_offset(:), send_local(:,:), result_slot(:)
  end type donor_plan

contains

  subroutine load_distributed_geometry(path, ni, nj, box, comm, physical_h, geom, error)
    character(*), intent(in) :: path
    integer, intent(in) :: ni, nj, box(4), comm ! local [i0,j0,ni,nj] (offsets zero-based)
    real(real_8), intent(in) :: physical_h
    type(distributed_geometry), intent(out) :: geom
    character(*), intent(out) :: error
    integer, allocatable :: boxes(:,:)
    integer :: ierr, rank, ranks, count, n, component, a, b, c, d, dimensions(2)
    real(real_8) :: metadata(8), grid_step

    call mpi_comm_rank(comm, rank, ierr)
    call mpi_comm_size(comm, ranks, ierr)
    error = ''
    dimensions = [ni,nj]
    grid_step = physical_h
    call mpi_bcast(dimensions,2,MPI_INTEGER,0,comm,ierr)
    call mpi_bcast(grid_step,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    if (any(dimensions /= [ni,nj]) .or. grid_step /= physical_h) error = 'geometry dimensions or grid step differ across ranks'
    call broadcast_error(comm,rank,ranks,error)
    if (len_trim(error) /= 0) return
    if (rank == 0) then
      if (int(ni,selected_int_kind(18))*nj > huge(ni)/2) then
        error = 'embedding is too large for version-1 face IDs'
      else
        call read_geometry(path,ni,nj,physical_h,geom%global,error)
      end if
    end if
    call mpi_bcast(error,len(error),MPI_CHARACTER,0,comm,ierr)
    if (len_trim(error) /= 0) return
    geom%rank = rank
    geom%ranks = ranks
    geom%comm = comm
    geom%i0 = box(1)
    geom%j0 = box(2)
    geom%ni = box(3)
    geom%nj = box(4)
    allocate(boxes(4,0:ranks-1),geom%owner(ni,nj))
    call mpi_allgather(box,4,MPI_INTEGER,boxes,4,MPI_INTEGER,comm,ierr)
    geom%owner = -1
    do n = 0,ranks-1
      if (any(boxes(:2,n) < 0) .or. any(boxes(3:,n) < 1)) then
        error = 'invalid geometry partition'
        return
      end if
      if (boxes(3,n) > ni-boxes(1,n) .or. boxes(4,n) > nj-boxes(2,n)) then
        error = 'geometry partition extends outside the embedding'
        return
      end if
      a = boxes(1,n)+1
      b = boxes(2,n)+1
      c = boxes(1,n)+boxes(3,n)
      d = boxes(2,n)+boxes(4,n)
      if (any(geom%owner(a:c,b:d) >= 0)) then
        error = 'overlapping geometry partitions'
        return
      end if
      geom%owner(a:c,b:d) = n
    end do
    if (any(geom%owner < 0)) then
      error = 'geometry partitions do not cover the embedding'
      return
    end if

    if (rank == 0) then
      count = size(geom%global%loops)
      metadata = [geom%global%bounds, geom%global%spacing, geom%global%distance_tolerance, &
                  geom%global%max_snap_distance, physical_h]
    end if
    call mpi_bcast(count,1,MPI_INTEGER,0,comm,ierr)
    call mpi_bcast(metadata,size(metadata),MPI_DOUBLE_PRECISION,0,comm,ierr)
    geom%h = metadata(8)
    if (rank /= 0) then
      allocate(geom%global%kind(ni,nj),geom%global%wedge(ni,nj),geom%global%component(ni,nj))
      allocate(geom%global%loops(0:count-1))
      geom%global%bounds = metadata(:4)
      geom%global%spacing = metadata(5)
      geom%global%distance_tolerance = metadata(6)
      geom%global%max_snap_distance = metadata(7)
    end if
    call mpi_bcast(geom%global%max_depth,1,MPI_INTEGER,0,comm,ierr)
    call mpi_bcast(geom%global%coordinate_units,len(geom%global%coordinate_units),MPI_CHARACTER,0,comm,ierr)
    call mpi_bcast(geom%global%digest,len(geom%global%digest),MPI_CHARACTER,0,comm,ierr)
    call mpi_bcast(geom%global%kind,ni*nj,MPI_INTEGER,0,comm,ierr)
    call mpi_bcast(geom%global%component,ni*nj,MPI_INTEGER,0,comm,ierr)
    call mpi_bcast(geom%global%wedge,ni*nj,MPI_INTEGER,0,comm,ierr)
    do component = 0,count-1
      if (rank == 0) n = size(geom%global%loops(component)%ij,2)
      call mpi_bcast(n,1,MPI_INTEGER,0,comm,ierr)
      if (rank /= 0) allocate(geom%global%loops(component)%ij(2,n))
      call mpi_bcast(geom%global%loops(component)%ij,2*n,MPI_INTEGER,0,comm,ierr)
    end do
    call derive_owned_geometry(geom, error)
    call broadcast_error(comm,rank,ranks,error)
    if (len_trim(error) /= 0) return
    call setup_face_exchange(geom, geom%face_comm, error)
    call broadcast_error(comm,rank,ranks,error)
  end subroutine load_distributed_geometry

  subroutine broadcast_error(comm, rank, ranks, error)
    integer, intent(in) :: comm, rank, ranks
    character(*), intent(inout) :: error
    integer :: failing_rank, first, ierr

    failing_rank = merge(rank,ranks,len_trim(error) /= 0)
    call mpi_allreduce(failing_rank,first,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (first < ranks) call mpi_bcast(error,len(error),MPI_CHARACTER,first,comm,ierr)
  end subroutine broadcast_error

  pure function half_coordinates(i,j,ni) result(q)
    integer, intent(in) :: i,j,ni
    integer :: q(2)

    q = [2*(j-1),2*(ni-i)]
  end function half_coordinates

  pure function make_sector_moments() result(sector)
    real(real_8) :: sector(n_moments,0:7)
    integer :: ray, next

    do ray = 0,7
      next = modulo(ray+1,8)
      sector(:,ray) = triangle_moments(0.5D0*[ray_x(ray),ray_y(ray)],0.5D0*[ray_x(next),ray_y(next)])
    end do
  end function make_sector_moments

  pure function cell_mean_moments(mask, sector) result(moments)
    integer, intent(in) :: mask
    real(real_8), intent(in) :: sector(n_moments,0:7)
    real(real_8) :: moments(n_moments)
    integer :: ray, count

    count = popcnt(mask)
    moments = 0.0D0
    if (count == 0) return
    do ray = 0,7
      if (btest(mask,ray)) moments = moments + sector(:,ray)
    end do
    moments = moments/count
  end function cell_mean_moments

  subroutine derive_owned_geometry(geom, error)
    type(distributed_geometry), intent(inout) :: geom
    character(*), intent(out) :: error
    type(open_face), allocatable :: faces(:)
    type(wall_segment), allocatable :: walls(:)
    integer :: i,j,gi,gj,ti,tj,ni,nj,mask,neighbor_mask,side,nface,nwall,n,k,component
    integer :: a(2),b(2),q(2),tangent(2),axis,offset(2)
    logical :: open_half(2), neighbor_half(2)
    real(real_8) :: delta(2), sector(n_moments,0:7), all_moments(n_moments,0:255)

    error = ''
    ni = size(geom%owner,1)
    nj = size(geom%owner,2)
    allocate(geom%area(geom%ni,geom%nj),geom%moments(n_moments,geom%ni,geom%nj))
    sector = make_sector_moments()
    do mask = 0,255
      all_moments(:,mask) = cell_mean_moments(mask,sector)
    end do
    allocate(geom%face_index(4,geom%ni,geom%nj),faces(4*geom%ni*geom%nj))
    geom%face_index = 0
    nface = 0
    do i = 1,geom%ni
      gi = geom%i0+i
      do j = 1,geom%nj
        gj = geom%j0+j
        mask = geom%global%wedge(gi,gj)
        geom%area(i,j) = (geom%h**2/8)*popcnt(mask)
        geom%moments(:,i,j) = all_moments(:,mask)
        do side = 1,4
          ti = gi+face_di(side)
          tj = gj+face_dj(side)
          neighbor_mask = 0
          if (ti >= 1 .and. ti <= ni .and. tj >= 1 .and. tj <= nj) neighbor_mask = geom%global%wedge(ti,tj)
          open_half = btest(mask,face_halves(:,side))
          neighbor_half = btest(neighbor_mask,face_halves(:,opposite(side)))
          if (any(open_half .neqv. neighbor_half)) then
            write(error,'(a,3(1x,i0))') 'inconsistent open-face geometry at i, j, side', gi,gj,side
            return
          end if
          if (.not. any(open_half)) cycle
          if (side > 2 .and. geom%owner(ti,tj) == geom%rank) then
            geom%face_index(side,i,j) = geom%face_index(opposite(side),ti-geom%i0,tj-geom%j0)
            cycle
          end if
          a = [min(gi,ti),min(gj,tj)]
          b = [max(gi,ti),max(gj,tj)]
          axis = merge(1,2,gi /= ti)
          q = half_coordinates(gi,gj,ni) + [face_dj(side),-face_di(side)]
          tangent = [2-axis,axis-1]
          nface = nface+1
          faces(nface)%id = 2*((a(1)-1)*nj+a(2)-1)+axis
          faces(nface)%owner = geom%owner(a(1),a(2))
          faces(nface)%cells(:,1) = a
          faces(nface)%cells(:,2) = b
          faces(nface)%endpoints(:,1) = q-merge(1,0,open_half(1))*tangent
          faces(nface)%endpoints(:,2) = q+merge(1,0,open_half(2))*tangent
          geom%face_index(side,i,j) = nface
        end do
      end do
    end do
    geom%faces = faces(:nface)

    allocate(walls(geom%ni*geom%nj))
    nwall = 0
    do component = 0,ubound(geom%global%loops,1)
      associate(ij => geom%global%loops(component)%ij)
        n = size(ij,2)
        do k = 1,n
          a = ij(:,k)
          if (geom%owner(a(1),a(2)) /= geom%rank) cycle
          b = ij(:,mod(k,n)+1)
          nwall = nwall+1
          walls(nwall)%component = component
          walls(nwall)%edge = k
          walls(nwall)%endpoints(:,1) = half_coordinates(a(1),a(2),ni)
          walls(nwall)%endpoints(:,2) = half_coordinates(b(1),b(2),ni)
          offset = walls(nwall)%endpoints(:,2)-walls(nwall)%endpoints(:,1)
          delta = offset*(geom%h/2)
          walls(nwall)%length = norm2(delta)
          walls(nwall)%outward = [delta(2),-delta(1)]/walls(nwall)%length
        end do
      end associate
    end do
    geom%walls = walls(:nwall)
  end subroutine derive_owned_geometry

  subroutine setup_face_exchange(geom, plan, error)
    type(distributed_geometry), intent(in) :: geom
    type(face_exchange), intent(out) :: plan
    character(*), intent(out) :: error
    integer :: n, f, r, k, id, axis, i, j, nj, ierr
    integer, allocatable :: cursor(:), send_ids(:), receive_ids(:)

    n = geom%ranks-1
    allocate(plan%send_count(0:n),plan%receive_count(0:n),plan%send_offset(0:n),plan%receive_offset(0:n),cursor(0:n))
    plan%send_count = 0
    do f = 1,size(geom%faces)
      if (geom%faces(f)%owner /= geom%rank) cycle
      associate(c => geom%faces(f)%cells(:,2))
        r = geom%owner(c(1),c(2))
      end associate
      if (r /= geom%rank) plan%send_count(r) = plan%send_count(r)+1
    end do
    call mpi_alltoall(plan%send_count,1,MPI_INTEGER,plan%receive_count,1,MPI_INTEGER,geom%comm,ierr)
    call offsets(plan%send_count,plan%send_offset)
    call offsets(plan%receive_count,plan%receive_offset)
    allocate(send_ids(sum(plan%send_count)),receive_ids(sum(plan%receive_count)))
    allocate(plan%send_face(size(send_ids)),plan%receive_face(size(receive_ids)))
    cursor = plan%send_offset
    do f = 1,size(geom%faces)
      if (geom%faces(f)%owner /= geom%rank) cycle
      associate(c => geom%faces(f)%cells(:,2))
        r = geom%owner(c(1),c(2))
      end associate
      if (r == geom%rank) cycle
      cursor(r) = cursor(r)+1
      send_ids(cursor(r)) = geom%faces(f)%id
      plan%send_face(cursor(r)) = f
    end do
    call mpi_alltoallv(send_ids,plan%send_count,plan%send_offset,MPI_INTEGER, &
                       receive_ids,plan%receive_count,plan%receive_offset,MPI_INTEGER,geom%comm,ierr)
    error = ''
    nj = size(geom%owner,2)
    do k = 1,size(receive_ids)
      id = receive_ids(k)
      axis = modulo(id-1,2)+1
      i = (id-1)/2/nj+1-geom%i0+face_di(axis)
      j = modulo((id-1)/2,nj)+1-geom%j0+face_dj(axis)
      if (i < 1 .or. i > geom%ni .or. j < 1 .or. j > geom%nj) then
        error = 'shared face ID was sent to the wrong owner'
        return
      end if
      f = geom%face_index(opposite(axis),i,j)
      if (f == 0) then
        error = 'shared face is missing on its receiving rank'
        return
      end if
      if (geom%faces(f)%id /= id) then
        error = 'shared face ID differs across its two incidences'
        return
      end if
      plan%receive_face(k) = f
    end do
  end subroutine setup_face_exchange

  ! Each owner supplies a single value. The adjacent rank receives it in its
  ! local face ordering. Fields of any kind can share this incidence plan.
  subroutine share_face_flux(geom, flux, error)
    type(distributed_geometry), intent(in) :: geom
    real(real_8), intent(inout) :: flux(:,:) ! (fields, local face IDs)
    character(*), intent(out) :: error
    real(real_8), allocatable :: send(:,:), receive(:,:)
    integer :: fields, ierr, counts(2)

    error = ''
    fields = size(flux,1)
    if (size(flux,2) /= size(geom%faces)) error = 'shared face flux has incompatible shape'
    call broadcast_error(geom%comm,geom%rank,geom%ranks,error)
    if (len_trim(error) /= 0) return
    call mpi_allreduce([fields,-fields],counts,2,MPI_INTEGER,MPI_MAX,geom%comm,ierr)
    if (counts(1) /= -counts(2)) then
      error = 'shared face field counts differ across ranks'
      return
    end if
    associate(plan => geom%face_comm)
      send = flux(:,plan%send_face)
      allocate(receive(fields,size(plan%receive_face)))
      call mpi_alltoallv(send,fields*plan%send_count,fields*plan%send_offset,MPI_DOUBLE_PRECISION, &
                         receive,fields*plan%receive_count,fields*plan%receive_offset,MPI_DOUBLE_PRECISION,geom%comm,ierr)
      flux(:,plan%receive_face) = receive
    end associate
  end subroutine share_face_flux

  ! Reverse the shared face exchange for a trace computed by cell 2.
  ! Each rank fills the slots for cells it owns. After this exchange the
  ! north/west face owner has both endpoint traces, including at MPI seams.
  subroutine gather_second_face_trace(geom,trace,error)
    type(distributed_geometry), intent(in) :: geom
    real(real_8), intent(inout) :: trace(:,:,:,:) ! (endpoint, cell, layer, face)
    character(*), intent(out) :: error
    real(real_8), allocatable :: send(:,:,:),receive(:,:,:)
    integer :: nl,fields,counts(2),ierr,k

    error = ''
    nl = size(trace,3)
    if (nl < 1 .or. any(shape(trace) /= [2,2,nl,size(geom%faces)])) error = 'face trace has incompatible shape'
    call broadcast_error(geom%comm,geom%rank,geom%ranks,error)
    if (len_trim(error) /= 0) return
    call mpi_allreduce([nl,-nl],counts,2,MPI_INTEGER,MPI_MAX,geom%comm,ierr)
    if (counts(1) /= -counts(2)) then
      error = 'face trace layer counts differ across ranks'
      return
    end if
    fields = 2*nl
    associate(plan => geom%face_comm)
      allocate(send(2,nl,size(plan%receive_face)),receive(2,nl,size(plan%send_face)))
      do k = 1,size(plan%receive_face)
        send(:,:,k) = trace(:,2,:,plan%receive_face(k))
      end do
      call mpi_alltoallv(send,fields*plan%receive_count,fields*plan%receive_offset,MPI_DOUBLE_PRECISION, &
                         receive,fields*plan%send_count,fields*plan%send_offset,MPI_DOUBLE_PRECISION,geom%comm,ierr)
      do k = 1,size(plan%send_face)
        trace(:,2,:,plan%send_face(k)) = receive(:,:,k)
      end do
    end associate
  end subroutine gather_second_face_trace

  subroutine offsets(counts, displacement)
    integer, intent(in) :: counts(0:)
    integer, intent(out) :: displacement(0:)
    integer :: r

    displacement(0) = 0
    do r = 1,ubound(counts,1)
      displacement(r) = displacement(r-1)+counts(r-1)
    end do
  end subroutine offsets

  ! points(:,k) is the global (i,j) of strict-fluid donor k.
  ! Sort requests by owner once; preserve the caller's result ordering.
  subroutine build_donor_plan(geom, points, plan, error, cell_means)
    type(distributed_geometry), intent(in) :: geom
    integer, intent(in) :: points(:,:)
    type(donor_plan), intent(out) :: plan
    character(*), intent(out) :: error
    logical, intent(in), optional :: cell_means
    integer, allocatable :: requests(:,:), cursor(:)
    integer :: k, i, j, r, n, nr, ierr
    logical :: allow_cut

    error = ''
    allow_cut = .false.
    if (present(cell_means)) allow_cut = cell_means
    n = size(points,2)
    nr = geom%ranks-1
    if (size(points,1) /= 2) error = 'donor coordinates must have shape (2,n)'
    do k = 1,n
      if (len_trim(error) /= 0) exit
      i = points(1,k)
      j = points(2,k)
      if (min(i,j) < 1 .or. i > size(geom%owner,1) .or. j > size(geom%owner,2)) then
        error = 'donor lies outside the embedding'
      else if (geom%global%kind(i,j) == solid_node .and. allow_cut) then
        write(error,'(a,2(1x,i0))') 'cell mean donor must be wet at', i,j
      else if (.not. allow_cut .and. geom%global%kind(i,j) /= fluid_node) then
        write(error,'(a,2(1x,i0))') 'donor must be a strict-fluid node at', i,j
      end if
    end do
    call broadcast_error(geom%comm,geom%rank,geom%ranks,error)
    if (len_trim(error) /= 0) return

    plan%comm = geom%comm
    plan%ni = geom%ni
    plan%nj = geom%nj
    allocate(plan%request_count(0:nr),plan%request_offset(0:nr),plan%send_count(0:nr),plan%send_offset(0:nr))
    allocate(plan%result_slot(n),requests(2,n),cursor(0:nr))
    plan%request_count = 0
    do k = 1,n
      r = geom%owner(points(1,k),points(2,k))
      plan%request_count(r) = plan%request_count(r)+1
    end do
    call offsets(plan%request_count,plan%request_offset)
    cursor = plan%request_offset
    do k = 1,n
      r = geom%owner(points(1,k),points(2,k))
      cursor(r) = cursor(r)+1
      requests(:,cursor(r)) = points(:,k)
      plan%result_slot(k) = cursor(r)
    end do
    call mpi_alltoall(plan%request_count,1,MPI_INTEGER,plan%send_count,1,MPI_INTEGER,plan%comm,ierr)
    call offsets(plan%send_count,plan%send_offset)
    allocate(plan%send_local(2,sum(plan%send_count)))
    call mpi_alltoallv(requests,2*plan%request_count,2*plan%request_offset,MPI_INTEGER, &
                       plan%send_local,2*plan%send_count,2*plan%send_offset,MPI_INTEGER,plan%comm,ierr)
    plan%send_local(1,:) = plan%send_local(1,:)-geom%i0
    plan%send_local(2,:) = plan%send_local(2,:)-geom%j0
  end subroutine build_donor_plan

  subroutine gather_donors(plan, field, samples)
    type(donor_plan), intent(in) :: plan
    real(real_8), intent(in) :: field(:,:,:) ! owned nodes only
    real(real_8), intent(out) :: samples(:,:) ! (layer, requested point)
    real(real_8), allocatable :: send(:,:), receive(:,:)
    integer :: nl, n, k, ierr, counts(2)

    nl = size(field,3)
    n = size(plan%result_slot)
    if (any(shape(field) /= [plan%ni,plan%nj,nl]) .or. any(shape(samples) /= [nl,n])) then
      labort('donor samples have incompatible shape')
    end if
    call mpi_allreduce([nl,-nl],counts,2,MPI_INTEGER,MPI_MAX,plan%comm,ierr)
    if (counts(1) /= -counts(2)) then
      labort('donor field counts differ across ranks')
    end if
    allocate(send(nl,size(plan%send_local,2)),receive(nl,n))
    do k = 1,size(send,2)
      send(:,k) = field(plan%send_local(1,k),plan%send_local(2,k),:)
    end do
    call mpi_alltoallv(send,nl*plan%send_count,nl*plan%send_offset,MPI_DOUBLE_PRECISION, &
                      receive,nl*plan%request_count,nl*plan%request_offset,MPI_DOUBLE_PRECISION,plan%comm,ierr)
    samples = receive(:,plan%result_slot)
  end subroutine gather_donors

end module qg_irregular_parallel
