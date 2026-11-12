#include "config.h"
program test_irregular_straight_wall
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  use qg_diagnostics, only: BC_NO_SLIP, BC_FREE_SLIP, BC_PARTIAL_SLIP
  use qg_geometry, only: fluid_node, wall_node, solid_node
  use qg_irregular_bc
  use qg_irregular_parallel, only: distributed_geometry, load_distributed_geometry
  use qg_parallel, only: parallel_init, parallel_finalise, pabort
  use qg_parameters, only: real_8
  implicit none
  include 'mpif.h'

  integer, parameter :: n = 33
  integer, parameter :: ox(0:7) = [1,1,0,-1,-1,-1,0,1], oy(0:7) = [0,1,1,1,0,-1,-1,-1]
  real(real_8), parameter :: wall = 0.37D0, curvature = 1.3D0, ell = 0.4D0
  type(distributed_geometry) :: geom
  type(straight_wall_plan) :: plan
  type(wall_slip) :: partial, no_slip, free, tiny_slip, large_slip
  real(real_8) :: nan_value, h, distance, errors(3,3), rates(2,3), maximum, local
  real(real_8) :: vertex(2), local_flux(2), reference(2), step, radius
  real(real_8), allocatable :: psi(:,:,:), omega(:,:,:), values(:,:), halo(:,:,:), transports(:,:), owned(:,:), shared(:,:)
  integer :: refinement, angle, step_start, rank, ranks, ierr, grid(2), row,column,box(4),i,j,gi,gj
  integer :: f,side,component,ray,q(2),wall_count,global_count,k
  character(1024) :: path,error
  character(30) :: fixture
  logical :: scalar_only

  call get_command_argument(1,path)
  scalar_only = command_argument_count() == 1 .and. path == '--scalar'
  if (.not. scalar_only .and. command_argument_count() /= 2) then
    print '(a)', 'Straight-wall fixture requires a .dat path and case name, or --scalar'
    stop 2
  end if

  call initialise_wall_slip(no_slip,BC_NO_SLIP,error)
  call require(len_trim(error) == 0,error)
  call initialise_wall_slip(free,BC_FREE_SLIP,error)
  call require(len_trim(error) == 0,error)
  call initialise_wall_slip(partial,BC_PARTIAL_SLIP,error,length=0.0D0)
  call require(len_trim(error) > 0,'accepted zero length as partial slip')
  call initialise_wall_slip(partial,BC_PARTIAL_SLIP,error,length=ell)
  call require(len_trim(error) == 0,error)
  call initialise_wall_slip(tiny_slip,BC_PARTIAL_SLIP,error,length=1e-200_real_8)
  call require(len_trim(error) == 0,error)
  call initialise_wall_slip(large_slip,BC_PARTIAL_SLIP,error,length=1e308_real_8)
  call require(len_trim(error) == 0,error)
  nan_value = ieee_value(0.0D0,ieee_quiet_nan)
  call require(straight_vorticity(free,wall,nan_value,nan_value,1.0D0) == 0.0D0,'free slip must be exact')
  call require(abs(straight_vorticity(tiny_slip,wall,1.0D0,2.0D0,1.0D0)- &
    straight_vorticity(no_slip,wall,1.0D0,2.0D0,1.0D0)) < 1e-15,'small-slip limit')
  call require(abs(straight_vorticity(large_slip,wall,1.0D0,2.0D0,1.0D0)) < 1e-300_real_8,'large-slip limit')
  do k = 1,2
    do refinement = 1,3
      h = 0.04D0/2**(refinement-1)
      distance = h*sqrt(real(k,real_8))
      errors(refinement,1) = abs(straight_vorticity(no_slip,wall,profile(distance,0.0D0), &
                                                     profile(2*distance,0.0D0),distance)-curvature)
      errors(refinement,2) = abs(straight_vorticity(partial,wall,profile(distance,ell), &
                                                     profile(2*distance,ell),distance)-curvature)
      ! Independent one-sided estimate of inward derivative, with no slip
      ! condition in the weights. For no slip the exact value is zero.
      errors(refinement,3) = abs((4*(profile(distance,0.0D0)-wall)-(profile(2*distance,0.0D0)-wall))/(2*distance))
    end do
    rates = log(errors(1:2,:)/errors(2:3,:))/log(2.0D0)
    call require(all(rates > 1.8D0),'normal Taylor or independent no-slip velocity order')
  end do
  if (scalar_only) then
    print '(a)', '[straight slip selectors and normal polynomial convergence] : PASS'
    stop
  end if

  call get_command_argument(1,path)
  call get_command_argument(2,fixture)
  call parallel_init()
  call mpi_comm_rank(MPI_COMM_WORLD,rank,ierr)
  call mpi_comm_size(MPI_COMM_WORLD,ranks,ierr)
  grid = 0
  call mpi_dims_create(ranks,2,grid,ierr)
  row = rank/grid(2)
  column = mod(rank,grid(2))
  box = [row*n/grid(1),column*n/grid(2), (row+1)*n/grid(1)-row*n/grid(1), (column+1)*n/grid(2)-column*n/grid(2)]
  call load_distributed_geometry(trim(path),n,n,box,MPI_COMM_WORLD,1.0D0,geom,error)
  call require(len_trim(error) == 0,error)

  allocate(psi(geom%ni,geom%nj,2),omega(geom%ni,geom%nj,2),halo(geom%ni+2,geom%nj+2,2))
  allocate(values(0:ubound(geom%global%loops,1),2),transports(size(geom%faces),2),owned(2,2*n*n),shared(2,2*n*n))
  if (fixture == 'rectangle' .or. fixture == 'diamond') then
    call setup_straight_walls(geom,plan,error)
    call require(len_trim(error) == 0,error)
    step_start = 0
    radius = 16
    if (fixture == 'diamond') then
      step_start = 1
      radius = 12
    end if
    do angle = step_start,7,2
      step = sqrt(real(ox(angle)**2+oy(angle)**2,real_8))
      psi = nan_value
      do j = 1,geom%nj
        gj = geom%j0+j
        do i = 1,geom%ni
          gi = geom%i0+i
          if (geom%global%kind(gi,gj) /= fluid_node) cycle
          distance = (radius-ox(angle)*(gj-17)-oy(angle)*(17-gi))/step
          psi(i,j,1) = wall + curvature*(ell*distance+0.5D0*distance**2) + 0.03D0*distance**3
          psi(i,j,2) = wall + curvature*(ell*distance+0.5D0*distance**2) - 0.01D0*distance**3
        end do
      end do
      values = wall
      omega = nan_value
      call reconstruct_straight_walls(plan,partial,psi,values,omega)
      wall_count = 0
      local = 0.0D0
      do j = 1,geom%nj
        gj = geom%j0+j
        do i = 1,geom%ni
          gi = geom%i0+i
          if (geom%global%kind(gi,gj) /= wall_node .or. popcnt(geom%global%wedge(gi,gj)) /= 4) then
            call require(all(.not.ieee_is_finite(omega(i,j,:))),'straight closure overwrote a non-straight value')
            cycle
          end if
          if (abs(radius-ox(angle)*(gj-17)-oy(angle)*(17-gi)) > 1e-15) cycle
          call require(all(ieee_is_finite(omega(i,j,:))),'straight closure sampled poisoned data')
          wall_count = wall_count+1
          local = max(local,maxval(abs(omega(i,j,:)-curvature)))
        end do
      end do
      call mpi_allreduce(local,maximum,1,MPI_DOUBLE_PRECISION,MPI_MAX,geom%comm,ierr)
      call mpi_allreduce(wall_count,global_count,1,MPI_INTEGER,MPI_SUM,geom%comm,ierr)
      call require(global_count > 0 .and. maximum < 1e-13,'axial or diagonal donor slip error')
    end do
  end if

  do component = 0,ubound(values,1)
    values(component,:) = [0.37D0+component,0.19D0-0.5D0*component]
  end do
  halo = nan_value
  do j = 1,geom%nj+2
    gj = geom%j0+j-1
    do i = 1,geom%ni+2
      gi = geom%i0+i-1
      if (min(gi,gj) < 1 .or. max(gi,gj) > n) cycle
      if (geom%global%kind(gi,gj) == solid_node) cycle
      halo(i,j,:) = [real(gj-1,real_8),real(n-gi,real_8)]
    end do
  end do
  call face_transports(geom,halo,values,transports)
  do f = 1,size(geom%faces)
    ! Affine sign/orientation check on open vertices well away from all walls.
    q = geom%faces(f)%endpoints(:,1)
    gi = n-q(2)/2
    gj = q(1)/2+1
    if (fixture == 'rectangle' .and. min(gi,gj) > 4 .and. max(gi,gj) < n-3) then
      reference = -0.5D0*(geom%faces(f)%endpoints(:,2)-q)
      call require(all(transports(f,:) == reference),'psi=x or psi=y transport sign')
    end if
  end do
  do j = 1,geom%nj+2
    gj = geom%j0+j-1
    do i = 1,geom%ni+2
      gi = geom%i0+i-1
      if (min(gi,gj) < 1 .or. max(gi,gj) > n) cycle
      component = geom%global%component(gi,gj)
      if (component >= 0) halo(i,j,:) = values(component,:)
    end do
  end do
  do k = 1,size(geom%walls)
    component = geom%walls(k)%component
    q = sum(geom%walls(k)%endpoints,dim=2)/2
    call vertex_streamfunction(geom,q,halo,values,vertex)
    call require(all(vertex == values(component,:)),'wall midpoint did not take component streamfunction')
  end do

  call face_transports(geom,halo,values,transports)
  owned = 0.0D0
  do f = 1,size(geom%faces)
    if (geom%faces(f)%owner == geom%rank) owned(:,geom%faces(f)%id) = transports(f,:)
  end do
  call mpi_allreduce(owned,shared,size(owned),MPI_DOUBLE_PRECISION,MPI_SUM,geom%comm,ierr)
  do f = 1,size(geom%faces)
    call require(all(transports(f,:) == shared(:,geom%faces(f)%id)),'rank-shared transport mismatch')
  end do
  local = 0.0D0
  do j = 1,geom%nj
    do i = 1,geom%ni
      local_flux = 0.0D0
      do side = 1,4
        f = geom%face_index(side,i,j)
        if (f > 0) local_flux = local_flux+merge(1,-1,side<=2)*transports(f,:)
      end do
      local = max(local,maxval(abs(local_flux)))
    end do
  end do
  call mpi_allreduce(local,maximum,1,MPI_DOUBLE_PRECISION,MPI_MAX,geom%comm,ierr)
  call require(maximum < 1e-13,'face transports failed polygon telescoping')
  if (rank == 0) print '(a,1x,i0,1x,a,1x,es12.4,1x,a)', trim(fixture), ranks, 'divergence:', maximum, 'PASS'
  call parallel_finalise()

contains
  real(real_8) function profile(distance,length)
    real(real_8), intent(in) :: distance,length
    profile = wall + curvature*(length*distance+0.5D0*distance**2) + 0.2D0*distance**3 + 0.4D0*distance**4
  end function profile

  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message
    if (.not. condition) then
      if (scalar_only) then
        print '(a)', trim(message)
        stop 1
      end if
      call pabort(message,__FILE__,__LINE__)
    end if
  end subroutine require
end program test_irregular_straight_wall
