#include "config.h"

module qg_irregular_bc
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use qg_diagnostics, only: BC_FREE_SLIP, BC_PARTIAL_SLIP, BC_NO_SLIP
  use qg_geometry, only: wall_node, solid_node
  use qg_irregular_parallel, only: distributed_geometry, donor_plan, build_donor_plan, gather_donors
  use qg_parallel, only: pabort
  use qg_parameters, only: real_8
  implicit none
  private
  public :: wall_slip, straight_wall_plan, initialise_wall_slip, straight_vorticity
  public :: setup_straight_walls, reconstruct_straight_walls, vertex_streamfunction, face_transports

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

  ! One lookup for both CABARET velocity time levels. An endpoint on an
  ! actual oriented wall edge takes that component's value exactly.
  subroutine vertex_streamfunction(geom, q, psi, wall_values, value)
    type(distributed_geometry), intent(in) :: geom
    integer, intent(in) :: q(2) ! half-grid physical (x,y) from the southwest embedding node
    real(real_8), intent(in) :: psi(:,:,:), wall_values(0:,:) ! psi includes one halo ring
    real(real_8), intent(out) :: value(:)
    integer :: i,j,ii,jj,component
    logical :: odd_x, odd_y

    if (any(q < 0) .or. q(1) > 2*(size(geom%owner,2)-1) .or. q(2) > 2*(size(geom%owner,1)-1)) then
      labort('vertex lies outside the embedding')
    end if
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
    if (component >= 0) then
      value = wall_values(component,:)
      return
    end if

    if (.not. (odd_x .and. odd_y)) then
      labort('open endpoint must be a dual corner or an actual wall crossing')
    end if
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

end module qg_irregular_bc
