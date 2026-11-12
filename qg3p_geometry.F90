! Version-1 geometry on a closed square-spacing embedding grid.
module qg_geometry
  use, intrinsic :: iso_fortran_env, only: iostat_end, iostat_eor
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use qg_parameters, only: real_8
  use qg_sha256, only: sha256, hash_line, hash_hex
  implicit none
  private
  public :: irregular_geometry, boundary_loop, read_geometry, solid_node, fluid_node, wall_node
  public :: write_rectangular_cape

  integer, parameter :: wide = selected_int_kind(18)
  integer, parameter :: solid_node = 0, fluid_node = 1, wall_node = 2
  character(*), parameter :: signature = 'PEQUOD_IRREGULAR 1'
  character(*), parameter :: index_convention = 'GLOBAL_1_BASED_I_SOUTH_J_EAST'
  integer, parameter :: di(0:7) = [0,-1,-1,-1,0,1,1,1]
  integer, parameter :: dj(0:7) = [1,1,0,-1,-1,-1,0,1]

  type boundary_loop
    integer, allocatable :: ij(:,:) ! (2,n), without a duplicate closing vertex
  end type boundary_loop

  type irregular_geometry
    real(real_8) :: bounds(4), spacing, distance_tolerance, max_snap_distance
    integer :: max_depth
    character(10) :: coordinate_units
    character(64) :: digest
    integer, allocatable :: kind(:,:), component(:,:), wedge(:,:)
    type(boundary_loop), allocatable :: loops(:) ! Lower bound 0, outer coast first
  end type irregular_geometry

contains

  ! The cape removes whole primal squares. Its boundary is one fluid-left
  ! loop. left is a one-based node column, width and depth count grid cells.
  subroutine write_rectangular_cape(path, ni, nj, left, width, depth, h, error)
    character(*), intent(in) :: path
    integer, intent(in) :: ni,nj,left,width,depth
    real(real_8), intent(in) :: h
    character(*), intent(out) :: error
    type(irregular_geometry) :: geom
    type(sha256) :: hash
    integer, parameter :: row_offset(4) = [-1,-1,0,0], col_offset(4) = [0,-1,-1,0]
    integer :: right,tip,corners(2,9),start(2),step(2),length,k,n,i,j,a,b,q,mask,unit,status
    character(16) :: hex(4)
    character(512) :: line

    error = ''
    if (min(ni,nj)<3 .or. int(ni,wide)*nj>huge(ni)/2 .or. .not. ieee_is_finite(h) .or. h<=0.0D0) then
      error = 'invalid cape embedding dimensions or spacing'
      return
    end if
    if (left<2 .or. left>=nj-1 .or. width<1 .or. width>nj-left-1 .or. depth<1 .or. depth>ni-3) then
      error = 'cape must lie between the side walls and leave a strict-fluid row south of its tip'
      return
    end if
    right = left+width
    tip = 1+depth
    geom%bounds = [0.0D0,(nj-1)*h,0.0D0,(ni-1)*h]
    geom%spacing = h
    geom%coordinate_units = 'cm'
    geom%distance_tolerance = h*1.0D-8
    geom%max_snap_distance = h
    geom%max_depth = 1
    if (.not. all(ieee_is_finite(geom%bounds)) .or. geom%distance_tolerance<=0.0D0) then
      error = 'cape bounds or geometry tolerance are outside the supported numeric range'
      return
    end if
    allocate(geom%kind(ni,nj),geom%component(ni,nj),geom%wedge(ni,nj),geom%loops(0:0))
    allocate(geom%loops(0)%ij(2,2*(ni+nj-2+depth)))
    corners = reshape([1,1, ni,1, ni,nj, 1,nj, 1,right, tip,right, tip,left, 1,left, 1,1],[2,9])
    n = 0
    do k = 1,size(corners,2)-1
      start = corners(:,k)
      step = corners(:,k+1)-start
      length = sum(abs(step))
      step = step/length
      do q = 0,length-1
        n = n+1
        geom%loops(0)%ij(:,n) = start+q*step
      end do
    end do
    do i = 1,ni
      do j = 1,nj
        mask = 0
        do q = 1,4
          a = i+row_offset(q)
          b = j+col_offset(q)
          if (a<1 .or. a>=ni .or. b<1 .or. b>=nj) cycle
          if (a<=depth .and. b>=left .and. b<right) cycle
          mask = ibset(ibset(mask,2*q-2),2*q-1)
        end do
        geom%wedge(i,j) = mask
        geom%component(i,j) = -1
        select case(mask)
          case(0)
            geom%kind(i,j) = solid_node
          case(255)
            geom%kind(i,j) = fluid_node
          case default
            geom%kind(i,j) = wall_node
            geom%component(i,j) = 0
        end select
      end do
    end do
    call verify_topology(geom,ni,nj,error)
    if (len_trim(error)/=0) return

    call records(.false.)
    geom%digest = hash_hex(hash)
    open(newunit=unit,file=path,status='new',action='write',iostat=status,iomsg=error)
    if (status/=0) return
    call records(.true.)
    if (status/=0) then
      close(unit)
    else
      close(unit,iostat=status,iomsg=error)
    end if

  contains

    subroutine records(writing)
      logical, intent(in) :: writing
      integer :: row,col,vertex

      call record(signature,writing)
      write(line,'(a,2(1x,i0))') 'GRID',ni,nj
      call record(line,writing)
      if (writing) then
        write(line,'(a,4(1x,es26.17e3))') 'BOUNDS',geom%bounds
      else
        do vertex = 1,4
          write(hex(vertex),'(z16.16)') transfer(geom%bounds(vertex),0_wide)
          call lower_hex(hex(vertex))
        end do
        write(line,'(a,4(1x,a))') 'BOUNDS',hex
      end if
      call record(line,writing)
      call record('COORD_UNITS cm',writing)
      call record('INDEX '//index_convention,writing)
      if (writing) then
        write(line,'(a,1x,es26.17e3,1x,i0,1x,es26.17e3)') &
          'SAMPLING',geom%distance_tolerance,geom%max_depth,geom%max_snap_distance
        call record(line,writing)
        call record('HASH SHA256 '//geom%digest,writing)
      end if
      write(line,'(a,1x,i0)') 'NODES',ni*nj
      call record(line,writing)
      do row = 1,ni
        do col = 1,nj
          write(line,'(i0,4(1x,i0))') row,col,geom%kind(row,col),geom%component(row,col),geom%wedge(row,col)
          call record(line,writing)
        end do
      end do
      call record('COMPONENTS 1',writing)
      write(line,'(a,1x,i0)') '0 OUTER',n
      call record(line,writing)
      do vertex = 1,n
        write(line,'(i0,1x,i0)') geom%loops(0)%ij(:,vertex)
        call record(line,writing)
      end do
      call record('END',writing)
    end subroutine records

    subroutine record(text,writing)
      character(*), intent(in) :: text
      logical, intent(in) :: writing

      if (writing) then
        if (status==0) write(unit,'(a)',iostat=status,iomsg=error) trim(text)
      else
        call hash_line(hash,text)
      end if
    end subroutine record
  end subroutine write_rectangular_cape

  ! Return an error to the caller so an MPI root can communicate one decision.
  subroutine read_geometry(path, ni, nj, physical_h, geom, error)
    character(*), intent(in) :: path
    integer, intent(in) :: ni, nj
    real(real_8), intent(in) :: physical_h
    type(irregular_geometry), intent(out) :: geom
    character(*), intent(out) :: error
    integer :: unit, status

    error = ''
    if (min(ni,nj) < 3 .or. int(ni,wide)*int(nj,wide) >= huge(ni)) then
      error = 'invalid embedding dimensions'
      return
    end if
    if (.not. ieee_is_finite(physical_h) .or. physical_h <= 0.0D0) then
      error = 'invalid physical grid spacing'
      return
    end if
    open(newunit=unit, file=path, status='old', action='read', iostat=status, iomsg=error)
    if (status /= 0) return
    call read_records(unit, ni, nj, physical_h, geom, error)
    close(unit)
    if (len_trim(error) == 0) call verify_topology(geom, ni, nj, error)
  end subroutine read_geometry

  subroutine read_records(unit, ni, nj, physical_h, geom, error)
    integer, intent(in) :: unit, ni, nj
    real(real_8), intent(in) :: physical_h
    type(irregular_geometry), intent(out) :: geom
    character(*), intent(out) :: error
    type(sha256) :: hash
    character(128) :: field(5)
    character(512) :: line
    character(16) :: hex(4)
    character(6) :: role
    integer :: nline, i, j, k, component, count, length, left, status
    integer :: values(5)
    real(real_8) :: dx, dy
    logical :: valid

    error = ''
    nline = 0
    call record(2, 'PEQUOD_IRREGULAR')
    if (.not. valid .or. field(2) /= '1') goto 900
    call hash_line(hash, signature)
    call record(3, 'GRID')
    call integers(2, values(:2))
    if (.not. valid) goto 900
    if (any(values(:2) /= [ni,nj])) then
      error = 'geometry dimensions differ from the decomposition'
      return
    end if
    write(line,'(a,2(1x,i0))') 'GRID', ni, nj
    call hash_line(hash, line)

    call record(5, 'BOUNDS')
    do k = 1,4
      geom%bounds(k) = decimal(k+1)
    end do
    if (.not. valid) goto 900
    dx = (geom%bounds(2) - geom%bounds(1))/(nj-1)
    dy = (geom%bounds(4) - geom%bounds(3))/(ni-1)
    if (.not. ieee_is_finite(dx) .or. .not. ieee_is_finite(dy) .or. min(dx,dy) <= 0.0D0) goto 900
    if (abs(dx-dy) > 1.0D-12*max(dx,dy)) then
      error = 'geometry requires square coordinate spacing'
      return
    end if
    geom%spacing = dx
    do k = 1,4
      write(hex(k),'(z16.16)') transfer(geom%bounds(k),0_wide)
      call lower_hex(hex(k))
    end do
    write(line,'(a,4(1x,a))') 'BOUNDS', hex
    call hash_line(hash, line)
    call record(2, 'COORD_UNITS')
    if (.not. valid .or. (field(2) /= 'normalized' .and. field(2) /= 'cm')) goto 900
    geom%coordinate_units = field(2)(:len(geom%coordinate_units))
    if (field(2) == 'cm' .and. abs(dx-physical_h) > 1.0D-12*max(dx,physical_h)) then
      error = 'geometry cm spacing differs from the model physical spacing'
      return
    end if
    call hash_line(hash, 'COORD_UNITS '//trim(geom%coordinate_units))
    call record(2, 'INDEX')
    if (.not. valid .or. field(2) /= index_convention) goto 900
    call hash_line(hash, 'INDEX '//index_convention)
    call record(4, 'SAMPLING')
    geom%distance_tolerance = decimal(2)
    geom%max_snap_distance = decimal(4)
    call integers(3, values(:1))
    if (.not. valid .or. min(geom%distance_tolerance,geom%max_snap_distance) <= 0.0D0) goto 900
    geom%max_depth = values(1)
    if (geom%max_depth < 1 .or. geom%max_depth > 10) goto 900
    call record(3, 'HASH')
    if (.not. valid .or. field(2) /= 'SHA256' .or. len_trim(field(3)) /= 64) goto 900
    if (verify(trim(field(3)), '0123456789abcdef') /= 0) goto 900
    geom%digest = field(3)(:len(geom%digest))
    call record(2, 'NODES')
    call integers(2, values(:1))
    if (.not. valid .or. values(1) /= ni*nj) goto 900
    call hash_ints('NODES', values(:1))
    allocate(geom%kind(ni,nj), geom%component(ni,nj), geom%wedge(ni,nj))
    do i = 1,ni
      do j = 1,nj
        call record(5)
        call integers(1, values)
        if (.not. valid .or. any(values(:2) /= [i,j])) goto 900
        geom%kind(i,j) = values(3)
        geom%component(i,j) = values(4)
        geom%wedge(i,j) = values(5)
        call hash_ints('', values)
      end do
    end do

    call record(2, 'COMPONENTS')
    call integers(2, values(:1))
    if (.not. valid) goto 900
    count = values(1)
    if (count < 1 .or. count > ni*nj/3) goto 900
    call hash_ints('COMPONENTS', values(:1))
    allocate(geom%loops(0:count-1))
    left = ni*nj
    do component = 0,count-1
      role = 'OUTER'
      if (component > 0) role = 'ISLAND'
      call record(3)
      call integers(1, values(:1))
      call integers(3, values(2:2))
      if (.not. valid .or. values(1) /= component .or. field(2) /= role) goto 900
      length = values(2)
      if (length < 3 .or. length > left) goto 900
      left = left - length
      write(line,'(i0,1x,a,1x,i0)') component, trim(role), length
      call hash_line(hash, line)
      allocate(geom%loops(component)%ij(2,length))
      do k = 1,length
        call record(2)
        call integers(1, values(:2))
        if (.not. valid) goto 900
        geom%loops(component)%ij(:,k) = values(:2)
        call hash_ints('', values(:2))
      end do
    end do
    call record(1, 'END')
    if (.not. valid) goto 900
    call hash_line(hash, 'END')
    do
      read(unit,'(a)',advance='no',iostat=status) line
      if (status == iostat_end) exit
      nline = nline + 1
      if (status /= iostat_eor .or. verify(trim(line), ' '//achar(9)//achar(13)) /= 0) goto 900
    end do
    if (hash_hex(hash) /= geom%digest) error = 'geometry SHA256 mismatch'
    return

900 write(error,'(a,i0)') 'invalid geometry record at line ', nline

  contains

    subroutine record(nfields, first)
      integer, intent(in) :: nfields
      character(*), optional, intent(in) :: first
      integer :: start, finish, l, n

      nline = nline + 1
      field = ''
      valid = .false.
      read(unit,'(a)',advance='no',iostat=status) line
      if (status /= iostat_eor) return
      l = len_trim(line)
      start = 1
      do n = 1,nfields
        do while (start <= l)
          if (.not. separator(line(start:start))) exit
          start = start + 1
        end do
        if (start > l) return
        finish = start
        do while (finish <= l)
          if (separator(line(finish:finish))) exit
          finish = finish + 1
        end do
        if (finish-start > len(field(n))) return
        field(n) = line(start:finish-1)
        start = finish
      end do
      if (start <= l) then
        if (verify(line(start:l), ' '//achar(9)//achar(13)) /= 0) return
      end if
      if (present(first)) then
        if (field(1) /= first) return
      end if
      valid = .true.
    end subroutine record

    logical function separator(ch)
      character, intent(in) :: ch
      separator = index(' '//achar(9)//achar(13), ch) /= 0
    end function separator

    subroutine integers(first, data)
      integer, intent(in) :: first
      integer, intent(out) :: data(:)
      integer :: n, stat

      data = 0
      if (.not. valid) return
      do n = 1,size(data)
        if (verify(trim(field(first+n-1)), '0123456789+-') /= 0) then
          valid = .false.
          return
        end if
        read(field(first+n-1),*,iostat=stat) data(n)
        if (stat /= 0) then
          valid = .false.
          return
        end if
      end do
    end subroutine integers

    real(real_8) function decimal(n) result(value)
      integer, intent(in) :: n
      integer :: stat

      value = 0.0D0
      if (.not. valid) return
      if (verify(trim(field(n)), '0123456789+-.eE') /= 0) then
        valid = .false.
        return
      end if
      read(field(n),*,iostat=stat) value
      valid = stat == 0 .and. ieee_is_finite(value)
    end function decimal

    subroutine hash_ints(prefix, data)
      character(*), intent(in) :: prefix
      integer, intent(in) :: data(:)

      write(line,'(a,*(1x,i0))') prefix, data
      call hash_line(hash, adjustl(line))
    end subroutine hash_ints
  end subroutine read_records

  subroutine lower_hex(hex)
    character(*), intent(inout) :: hex
    integer :: k

    do k = 1,len(hex)
      if (hex(k:k) >= 'A' .and. hex(k:k) <= 'F') hex(k:k) = achar(iachar(hex(k:k)) + 32)
    end do
  end subroutine lower_hex

  logical function before(a,b)
    integer, intent(in) :: a(2), b(2)

    before = a(1) < b(1) .or. (a(1) == b(1) .and. a(2) < b(2))
  end function before

  integer function ray(offset)
    integer, intent(in) :: offset(2)
    integer :: r

    ray = -1
    do r = 0,7
      if (offset(1) == di(r) .and. offset(2) == dj(r)) ray = r
    end do
  end function ray

  logical function inside(ij, point)
    integer, intent(in) :: ij(:,:), point(2)
    integer :: k, q, north

    inside = .false.
    do k = 1,size(ij,2)
      q = mod(k,size(ij,2))+1
      if (ij(1,k) == ij(1,q)) cycle
      north = k
      if (ij(1,q) < ij(1,k)) north = q
      if (ij(1,north) == point(1) .and. ij(2,north) > point(2)) inside = .not. inside
    end do
  end function inside

  subroutine verify_topology(geom, ni, nj, error)
    type(irregular_geometry), intent(in) :: geom
    integer, intent(in) :: ni,nj
    character(*), intent(out) :: error
    integer, allocatable :: expected_wall(:,:), masks(:,:), queue(:,:)
    logical, allocatable :: crosses(:,:), diagonal(:,:), visited(:,:)
    integer :: component, other, k, r, n, i, j, ii, jj, next, prev, outgoing, backward
    integer :: sector_count, fluid_count
    integer :: a(2), b(2), offset(2), cell(2), head, tail, kind
    integer(wide) :: area2, total2, eighths
    logical :: fluid

    error = ''
    allocate(expected_wall(ni,nj), masks(ni,nj), crosses(ni,nj), diagonal(ni-1,nj-1))
    expected_wall = -1
    masks = 0
    crosses = .false.
    diagonal = .false.
    total2 = 0
    do component = 0,ubound(geom%loops,1)
      associate(ij => geom%loops(component)%ij)
        n = size(ij,2)
        if (any(ij(1,:) < 1) .or. any(ij(1,:) > ni) .or. any(ij(2,:) < 1) .or. any(ij(2,:) > nj)) then
          error = 'boundary node is outside the embedding'
          return
        end if
        if (component > 1) then
          if (.not. before(geom%loops(component-1)%ij(:,1),ij(:,1))) then
            error = 'island components are out of order'
            return
          end if
        end if
        area2 = 0
        do k = 1,n
          a = ij(:,k)
          i = a(1)
          j = a(2)
          if (before(a,ij(:,1)) .or. expected_wall(i,j) >= 0) then
            write(error,'(a,2(1x,i0))') 'noncanonical, repeated or shared boundary node', a
            return
          end if
          expected_wall(i,j) = component
          next = mod(k,n)+1
          prev = modulo(k-2,n)+1
          b = ij(:,next)
          offset = b-a
          outgoing = ray(offset)
          backward = ray(ij(:,prev)-a)
          if (min(outgoing,backward) < 0 .or. outgoing == backward) then
            write(error,'(a,2(1x,i0))') 'nonadjacent or backtracking boundary step at', a
            return
          end if
          sector_count = modulo(backward-outgoing,8)
          do r = 0,sector_count-1
            masks(i,j) = ibset(masks(i,j),modulo(outgoing+r,8))
          end do
          if (all(offset /= 0)) then
            cell = min(a,b)
            if (diagonal(cell(1),cell(2))) then
              write(error,'(a,2(1x,i0))') 'crossing diagonals in primal cell', cell
              return
            end if
            diagonal(cell(1),cell(2)) = .true.
          end if
          if (offset(1) /= 0) then
            if (offset(1) > 0) b = a
            crosses(b(1),b(2)) = .not. crosses(b(1),b(2))
          end if
          area2 = area2 + int(i,wide)*ij(2,next) - int(j,wide)*ij(1,next)
        end do
        if ((component == 0 .and. area2 <= 0) .or. (component > 0 .and. area2 >= 0)) then
          error = 'invalid component orientation or zero area'
          return
        end if
        total2 = total2 + area2
      end associate
    end do

    do component = 1,ubound(geom%loops,1)
      a = geom%loops(component)%ij(:,1)
      if (.not. inside(geom%loops(0)%ij,a)) then
        error = 'island is outside the outer coastline'
        return
      end if
      do other = 1,ubound(geom%loops,1)
        if (other == component) cycle
        if (inside(geom%loops(other)%ij,a)) then
          error = 'nested island would enclose a separate basin'
          return
        end if
      end do
    end do

    eighths = 0
    do i = 1,ni
      fluid = .false.
      do j = nj,1,-1
        kind = solid_node
        if (expected_wall(i,j) >= 0) then
          kind = wall_node
        else if (fluid) then
          kind = fluid_node
          masks(i,j) = 255
        end if
        if (geom%kind(i,j) /= kind .or. geom%wedge(i,j) /= masks(i,j) .or. &
            geom%component(i,j) /= expected_wall(i,j)) then
          write(error,'(a,2(1x,i0))') 'node flags disagree with reconstructed polygon at', i,j
          return
        end if
        eighths = eighths + popcnt(masks(i,j))
        if (crosses(i,j)) fluid = .not. fluid
      end do
    end do
    if (eighths /= 4*total2) then
      error = 'cut-cell eighths do not sum to polygon area'
      return
    end if

    fluid_count = count(geom%kind == fluid_node)
    if (fluid_count == 0) then
      error = 'the snapped domain has no strict-fluid nodes'
      return
    end if
    allocate(queue(2,fluid_count), visited(ni,nj))
    visited = .false.
    queue(:,1) = maxloc(merge(1,0,geom%kind == fluid_node))
    head = 1
    tail = 1
    visited(queue(1,1),queue(2,1)) = .true.
    do while (head <= tail)
      i = queue(1,head)
      j = queue(2,head)
      head = head+1
      do r = 0,7,2
        ii = i + di(r)
        jj = j + dj(r)
        if (ii < 1 .or. ii > ni .or. jj < 1 .or. jj > nj) then
          write(error,'(a,2(1x,i0))') 'five-point stencil leaves the embedding at', i,j
          return
        end if
        if (geom%kind(ii,jj) == solid_node) then
          write(error,'(a,2(1x,i0))') 'five-point stencil leaks into land at', ii,jj
          return
        end if
        if (geom%kind(ii,jj) /= fluid_node .or. visited(ii,jj)) cycle
        tail = tail+1
        queue(:,tail) = [ii,jj]
        visited(ii,jj) = .true.
      end do
    end do
    if (tail /= fluid_count) then
      error = 'disconnected strict-fluid graph'
      return
    end if
    do j = 1,nj
      do i = 1,ni
        if (geom%kind(i,j) /= wall_node) cycle
        do r = 0,2,2
          ii = i + di(r)
          jj = j + dj(r)
          if (ii < 1 .or. jj > nj) cycle
          if (geom%kind(ii,jj) == wall_node .and. geom%component(i,j) /= geom%component(ii,jj)) then
            write(error,'(a,2(1x,i0))') 'unresolved one-spacing gap between components at', i,j
            return
          end if
        end do
      end do
    end do
  end subroutine verify_topology

end module qg_geometry
