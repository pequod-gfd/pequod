#include "config.h"

module qg_irregular_bc
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use qg_diagnostics, only: BC_FREE_SLIP, BC_PARTIAL_SLIP, BC_NO_SLIP
  use qg_geometry, only: wall_node, solid_node
  use qg_irregular_parallel, only: distributed_geometry, donor_plan, build_donor_plan, gather_donors, &
                                   share_face_flux, broadcast_error
  use qg_reconstruction, only: reconstruction_plan, cubic_terms, build_reconstruction, evaluate_reconstruction
  use qg_parallel, only: pabort
  use qg_parameters, only: real_8
  implicit none
  include 'mpif.h'
  private
  public :: wall_slip, straight_wall_plan, initialise_wall_slip, straight_vorticity
  public :: setup_straight_walls, reconstruct_straight_walls, vertex_streamfunction, face_transports
  public :: vertex_plan, setup_vertex_plan, instantaneous_transport, extrapolated_transport

  integer, parameter :: di(0:7) = [0,-1,-1,-1,0,1,1,1]
  integer, parameter :: dj(0:7) = [1,1,0,-1,-1,-1,0,1]

  type wall_slip
    private
    integer :: kind = BC_FREE_SLIP
    real(real_8) :: length = 0.0D0 ! in the same units as the normal distance
  end type wall_slip

  type straight_wall_plan
    private
    type(donor_plan) :: donors
    integer, allocatable :: site(:,:) ! (owned i, owned j, component) in stencil order
    real(real_8), allocatable :: normal_step(:)
  end type straight_wall_plan

  type vertex_plan
    private
    type(reconstruction_plan) :: reconstruction
    integer, allocatable :: owned_face(:), endpoint(:,:), midpoint(:)
    integer :: open_vertices
  end type vertex_plan

contains

  subroutine initialise_wall_slip(slip, kind, error, length)
    type(wall_slip), intent(out) :: slip
    integer, intent(in) :: kind
    character(*), intent(out) :: error
    real(real_8), optional, intent(in) :: length

    error = ''
    select case(kind)
      case(BC_FREE_SLIP,BC_NO_SLIP)
        slip%kind = kind
      case(BC_PARTIAL_SLIP)
        if (.not. present(length)) then
          error = 'partial slip requires a positive finite slip length'
        else if (length <= 0.0D0 .or. .not. ieee_is_finite(length)) then
          error = 'partial slip requires a positive finite slip length'
        else
          slip%kind = kind
          slip%length = length
        end if
      case default
        error = 'invalid irregular wall slip selector'
    end select
  end subroutine initialise_wall_slip

  elemental real(real_8) function straight_vorticity(slip, wall, near, far, distance) result(omega)
    type(wall_slip), intent(in) :: slip
    real(real_8), intent(in) :: wall, near, far, distance
    real(real_8) :: scale, difference

    omega = 0.0D0
    if (slip%kind == BC_FREE_SLIP) return
    ! Eliminate the cubic term using samples at s and 2s. In strong
    ! no-slip the length is zero, preserving the same second-order limit.
    scale = max(distance,slip%length)
    difference = 8.0D0*(near-wall)-(far-wall)
    omega = (difference/scale)/(distance*(2.0D0*(distance/scale)+6.0D0*(slip%length/scale)))
  end function straight_vorticity

  subroutine setup_straight_walls(geom, plan, error)
    type(distributed_geometry), intent(in) :: geom
    type(straight_wall_plan), intent(out) :: plan
    character(*), intent(out) :: error
    integer, allocatable :: points(:,:), site(:,:)
    real(real_8), allocatable :: step(:)
    integer :: edge, qx, qy, gi,gj,i,j,mask,ray,normal,n,offset(2)

    allocate(points(2,2*size(geom%walls)),site(3,size(geom%walls)),step(size(geom%walls)))
    n = 0
    do edge = 1,size(geom%walls)
      qx = geom%walls(edge)%endpoints(1,1)
      qy = geom%walls(edge)%endpoints(2,1)
      gi = size(geom%owner,1)-qy/2
      gj = 1+qx/2
      mask = geom%global%wedge(gi,gj)
      if (popcnt(mask) /= 4) cycle
      do ray = 0,7
        if (btest(mask,ray) .and. .not. btest(mask,modulo(ray-1,8))) exit
      end do
      normal = modulo(ray+2,8)
      offset = [di(normal),dj(normal)]
      n = n+1
      i = gi-geom%i0
      j = gj-geom%j0
      site(:,n) = [i,j,geom%global%component(gi,gj)]
      points(:,2*n-1) = [gi,gj]+offset
      points(:,2*n) = [gi,gj]+2*offset
      step(n) = sqrt(real(sum(offset**2),real_8))
    end do
    ! A unit normal-ray step stays in one primal cell or edge. With the
    ! straight start and strict-fluid endpoints, it cannot cross a wall
    ! once and finish on the same side, or cross two distinct diagonals.
    call build_donor_plan(geom,points(:,:2*n),plan%donors,error)
    if (len_trim(error) /= 0) then
      error = 'straight-wall normal support: '//trim(error)
      return
    end if
    plan%site = site(:,:n)
    plan%normal_step = step(:n)
  end subroutine setup_straight_walls

  ! Only writes straight wall point values. Corner closure is separate.
  subroutine reconstruct_straight_walls(plan, slip, psi, wall_values, omega)
    type(straight_wall_plan), intent(in) :: plan
    type(wall_slip), intent(in) :: slip
    real(real_8), intent(in) :: psi(:,:,:), wall_values(0:,:)
    real(real_8), intent(inout) :: omega(:,:,:)
    real(real_8), allocatable :: values(:,:)
    integer :: k,n,i,j,component

    n = size(plan%site,2)
    allocate(values(size(psi,3),2*n))
    call gather_donors(plan%donors,psi,values)
    do k = 1,n
      i = plan%site(1,k)
      j = plan%site(2,k)
      component = plan%site(3,k)
      omega(i,j,:) = straight_vorticity(slip,wall_values(component,:),values(:,2*k-1),values(:,2*k),plan%normal_step(k))
    end do
  end subroutine reconstruct_straight_walls

  integer function boundary_ray_component(geom,i,j,ray) result(component)
    type(distributed_geometry), intent(in) :: geom
    integer, intent(in) :: i,j,ray
    integer :: mask

    component = -1
    if (geom%global%kind(i,j) /= wall_node) return
    mask = geom%global%wedge(i,j)
    if (btest(mask,ray) .neqv. btest(mask,modulo(ray-1,8))) component = geom%global%component(i,j)
  end function boundary_ray_component

  integer function vertex_component(geom,q) result(component)
    type(distributed_geometry), intent(in) :: geom
    integer, intent(in) :: q(2) ! half-grid physical (x,y)
    integer :: i,j
    logical :: odd_x, odd_y

    i = size(geom%owner,1)-q(2)/2
    j = 1+q(1)/2
    odd_x = mod(q(1),2) /= 0
    odd_y = mod(q(2),2) /= 0
    if (odd_x .and. odd_y) then
      component = boundary_ray_component(geom,i,j,1) ! SW to NE
      if (component < 0) component = boundary_ray_component(geom,i-1,j,7) ! NW to SE
    else if (odd_x) then
      component = boundary_ray_component(geom,i,j,0)
    else if (odd_y) then
      component = boundary_ray_component(geom,i,j,2)
    else
      component = geom%global%component(i,j)
    end if
  end function vertex_component

  ! Four-centre comparison path retained for the R2 regression.
  subroutine vertex_streamfunction(geom, q, psi, wall_values, value)
    type(distributed_geometry), intent(in) :: geom
    integer, intent(in) :: q(2)
    real(real_8), intent(in) :: psi(:,:,:), wall_values(0:,:)
    real(real_8), intent(out) :: value(:)
    integer :: i,j,ii,jj,component

    if (any(q < 0) .or. q(1) > 2*(size(geom%owner,2)-1) .or. q(2) > 2*(size(geom%owner,1)-1)) then
      labort('vertex lies outside the embedding')
    end if
    component = vertex_component(geom,q)
    if (component >= 0) then
      value = wall_values(component,:)
      return
    end if

    if (any(mod(q,2) == 0)) then
      labort('open endpoint must be a dual corner or an actual wall crossing')
    end if
    i = size(geom%owner,1)-q(2)/2
    j = 1+q(1)/2
    if (any(geom%global%kind(i-1:i,j:j+1) == solid_node)) then
      labort('interior vertex interpolation would read a solid node')
    end if
    ii = i-geom%i0+1
    jj = j-geom%j0+1
    if (ii < 2 .or. ii > size(psi,1) .or. jj < 1 .or. jj >= size(psi,2)) then
      labort('open vertex donors lie outside the available halo')
    end if
    value = 0.25D0*(psi(ii-1,jj,:)+psi(ii,jj,:)+psi(ii-1,jj+1,:)+psi(ii,jj+1,:))
  end subroutine vertex_streamfunction

  subroutine face_transports(geom, psi, wall_values, transport)
    type(distributed_geometry), intent(in) :: geom
    real(real_8), intent(in) :: psi(:,:,:),wall_values(0:,:)
    real(real_8), intent(out) :: transport(:,:)
    real(real_8) :: first(size(psi,3)),second(size(psi,3))
    integer :: f

    do f = 1,size(geom%faces)
      call vertex_streamfunction(geom,geom%faces(f)%endpoints(:,1),psi,wall_values,first)
      call vertex_streamfunction(geom,geom%faces(f)%endpoints(:,2),psi,wall_values,second)
      ! Endpoints run east on a south face and north on an east face.
      ! For u = k cross grad(psi), the right-normal flux is psi(a)-psi(b).
      transport(f,:) = first-second
    end do
  end subroutine face_transports

  ! Cubic point reproduction at all open vertices avoids an order switch
  ! between wall-adjacent and interior faces. Wall sites use the exact C.
  subroutine setup_vertex_plan(geom,plan,error,normal_profile)
    type(distributed_geometry), intent(in) :: geom
    type(vertex_plan), intent(out) :: plan
    character(*), intent(out) :: error
    integer, allocatable :: lookup(:,:), origins(:,:)
    real(real_8), allocatable :: target(:,:)
    logical, optional, intent(in) :: normal_profile
    integer :: n,f,k,j,lo(2),hi(2),p(2),component,sites
    logical :: profile

    n = count(geom%faces%owner == geom%rank)
    profile = .false.
    if (present(normal_profile)) profile = normal_profile
    sites = merge(3,2,profile)*n
    allocate(plan%owned_face(n),plan%endpoint(2,n),origins(2,sites),target(cubic_terms,sites))
    if (profile) allocate(plan%midpoint(n))
    lo = huge(0)
    hi = -huge(0)
    if (n == 0) then
      lo = 1
      hi = 0
    end if
    k = 0
    do f = 1,size(geom%faces)
      if (geom%faces(f)%owner /= geom%rank) cycle
      k = k+1
      plan%owned_face(k) = f
      lo = min(lo,minval(geom%faces(f)%endpoints,dim=2))
      hi = max(hi,maxval(geom%faces(f)%endpoints,dim=2))
    end do
    allocate(lookup(lo(1):hi(1),lo(2):hi(2)))
    lookup = 0
    plan%open_vertices = 0
    do k = 1,n
      f = plan%owned_face(k)
      do j = 1,2
        p = geom%faces(f)%endpoints(:,j)
        if (lookup(p(1),p(2)) == 0) then
          component = vertex_component(geom,p)
          if (component >= 0) then
            lookup(p(1),p(2)) = -component-1
          else
            plan%open_vertices = plan%open_vertices+1
            lookup(p(1),p(2)) = plan%open_vertices
            origins(:,plan%open_vertices) = 2*p
          end if
        end if
        plan%endpoint(j,k) = lookup(p(1),p(2))
      end do
    end do
    target = 0.0D0
    target(1,:) = 1.0D0
    sites = plan%open_vertices
    if (profile) then
      do k = 1,n
        sites = sites+1
        plan%midpoint(k) = sites
        origins(:,sites) = sum(geom%faces(plan%owned_face(k))%endpoints,dim=2)
      end do
    end if
    call build_reconstruction(geom,origins(:,:sites),target(:,:sites),.false.,plan%reconstruction,error)
  end subroutine setup_vertex_plan

  subroutine instantaneous_transport(geom,plan,psi,wall_values,transport,error,normal_velocity)
    type(distributed_geometry), intent(in) :: geom
    type(vertex_plan), intent(in) :: plan
    real(real_8), intent(in) :: psi(:,:,:), wall_values(0:,:) ! owned nodal psi
    real(real_8), intent(out) :: transport(:,:)
    character(*), intent(out) :: error
    real(real_8), optional, intent(out) :: normal_velocity(:,:,:) ! (endpoint, layer, face), owner only
    real(real_8), allocatable :: value(:,:), replicated(:,:)
    real(real_8) :: length,a(size(psi,3)),b(size(psi,3)),middle(size(psi,3))
    integer :: k,f,nl,ierr,counts(2),sites

    error = ''
    nl = size(psi,3)
    if (any(shape(psi) /= [geom%ni,geom%nj,nl]) .or. any(shape(transport) /= [nl,size(geom%faces)]) .or. &
        any(shape(wall_values) /= [size(geom%global%loops),nl])) then
      error = 'vertex transport has incompatible shapes'
    else if (.not. all(ieee_is_finite(wall_values))) then
      error = 'vertex wall streamfunction must be finite'
    end if
    if (present(normal_velocity)) then
      if (.not. allocated(plan%midpoint) .or. any(shape(normal_velocity) /= [2,nl,size(geom%faces)])) &
        error = 'normal velocity requires a midpoint plan and compatible shape'
    end if
    call broadcast_error(geom%comm,geom%rank,geom%ranks,error)
    if (len_trim(error) /= 0) return
    call mpi_allreduce([nl,-nl],counts,2,MPI_INTEGER,MPI_MAX,geom%comm,ierr)
    if (counts(1) /= -counts(2)) then
      error = 'vertex layer counts differ across ranks'
      return
    end if
    replicated = wall_values
    call mpi_bcast(replicated,size(replicated),MPI_DOUBLE_PRECISION,0,geom%comm,ierr)
    if (any(wall_values /= replicated)) error = 'vertex wall streamfunction differs across ranks'
    call broadcast_error(geom%comm,geom%rank,geom%ranks,error)
    if (len_trim(error) /= 0) return

    sites = plan%open_vertices
    if (allocated(plan%midpoint)) sites = sites+size(plan%midpoint)
    allocate(value(nl,sites))
    call evaluate_reconstruction(plan%reconstruction,psi,value,wall_values)
    transport = 0.0D0
    if (present(normal_velocity)) normal_velocity = 0.0D0
    do k = 1,size(plan%owned_face)
      f = plan%owned_face(k)
      a = at_vertex(plan%endpoint(1,k))
      b = at_vertex(plan%endpoint(2,k))
      transport(:,f) = a-b
      if (present(normal_velocity)) then
        middle = value(:,plan%midpoint(k))
        length = 0.5D0*sum(abs(geom%faces(f)%endpoints(:,2)-geom%faces(f)%endpoints(:,1)))
        ! Negative tangential derivative of the endpoint/midpoint quadratic.
        ! The integral of this linear normal velocity is exactly a-b.
        normal_velocity(1,:,f) = (3*(a-middle)-(middle-b))/length
        normal_velocity(2,:,f) = (3*(middle-b)-(a-middle))/length
      end if
    end do
    call share_face_flux(geom,transport,error)
  contains
    function at_vertex(index) result(psi_at)
      integer, intent(in) :: index
      real(real_8) :: psi_at(nl)

      if (index > 0) then
        psi_at = value(:,index)
      else
        psi_at = wall_values(-index-1,:)
      end if
    end function at_vertex
  end subroutine instantaneous_transport

  subroutine extrapolated_transport(geom,plan,psi,wall_values,weight,history,transport,error)
    type(distributed_geometry), intent(in) :: geom
    type(vertex_plan), intent(in) :: plan
    real(real_8), intent(in) :: psi(:,:,:), wall_values(0:,:), weight
    real(real_8), intent(inout) :: history(:,:) ! preceding instantaneous transport
    real(real_8), intent(out) :: transport(:,:)
    character(*), intent(out) :: error
    real(real_8), allocatable :: instantaneous(:,:)
    real(real_8) :: extrema(2)
    integer :: ierr

    error = ''
    if (any(shape(history) /= shape(transport))) error = 'vertex transport history has incompatible shape'
    if (.not. ieee_is_finite(weight) .or. weight < 0.0D0 .or. weight > 1.0D0) error = 'invalid vertex extrapolation weight'
    call broadcast_error(geom%comm,geom%rank,geom%ranks,error)
    if (len_trim(error) /= 0) return
    call mpi_allreduce([weight,-weight],extrema,2,MPI_DOUBLE_PRECISION,MPI_MAX,geom%comm,ierr)
    if (extrema(1) /= -extrema(2)) then
      error = 'vertex extrapolation weight differs across ranks'
      return
    end if
    allocate(instantaneous(size(transport,1),size(transport,2)))
    call instantaneous_transport(geom,plan,psi,wall_values,instantaneous,error)
    if (len_trim(error) /= 0) return
    transport = instantaneous
    if (weight > 0.0D0) transport = (1+weight)*instantaneous-weight*history
    call share_face_flux(geom,transport,error)
    if (len_trim(error) == 0) history = instantaneous
  end subroutine extrapolated_transport

end module qg_irregular_bc
