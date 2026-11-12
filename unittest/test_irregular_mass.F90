#include "config.h"
program test_irregular_mass
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  use qg_aux, only: helmholtz_cache, allocate, deallocate, pv_inversion
  use qg_constraints
  use qg_diagnostics, only: stratification, velocity_boundary_condition, allocate, deallocate, &
                            initialise, BC_FREE_SLIP, z_to_modes, modes_to_z, mean
  use qg_geometry, only: fluid_node, wall_node, solid_node
  use qg_irregular_diagnostics
  use qg_irregular_parallel, only: distributed_geometry, load_distributed_geometry
  use qg_parallel, only: decomposition, parallel_data, allocate, deallocate, parallel_init, parallel_finalise, pabort
  use qg_parameters, only: real_8
  implicit none
  include 'mpif.h'

  integer, parameter :: n = 33, nl = 3
  real(real_8), parameter :: gauge(nl) = [0.0D0,0.7D0,-0.3D0], h = 0.125D0
  type(decomposition), target :: d
  type(distributed_geometry), target :: g
  type(stratification), target :: s
  type(velocity_boundary_condition) :: bc
  type(helmholtz_cache) :: legacy
  type(parallel_data) :: old_psi, old_zeta
  type(basin_mass_cache) :: cache
  real(real_8) :: exact(n,n,nl), target(nl), sums(nl), weights(n,n), modal(nl), area, potential, kinetic(nl)
  real(real_8) :: local_error, error_max, gauge_error, x,y, expected_area, expected_potential, layer_values(nl), recovered(nl)
  real(real_8), allocatable :: pv(:,:,:), psi(:,:,:), halo(:,:,:), omega(:,:,:), wall_values(:,:)
  integer :: i,j,gi,gj,k,ranks,grid(2),ierr,kind
  integer, allocatable :: iis(:),jjs(:)
  character(1024) :: path, error
  character(30) :: fixture

  if (command_argument_count() == 0) then
    print '(a)', 'Mass fixture requires a .dat path and case name'
    stop
  end if
  call get_command_argument(1,path)
  call get_command_argument(2,fixture)
  call parallel_init()
  call mpi_comm_size(MPI_COMM_WORLD,ranks,ierr)
  grid = 0
  call mpi_dims_create(ranks,2,grid,ierr)
  iis = [(i*n/grid(1)-(i-1)*n/grid(1),i=1,grid(1))]
  jjs = [(j*n/grid(2)-(j-1)*n/grid(2),j=1,grid(2))]
  call allocate(d,iis=iis,jjs=jjs,layers=nl)
  call load_distributed_geometry(trim(path),n,n,[d%ii_0,d%jj_0,d%ii,d%jj],d%cart_comm,h,g,error)
  call require(len_trim(error) == 0,error)
  call allocate(s,[1.0D0,2.0D0,4.0D0],[0.17D0,0.08D0],[0.34D0,0.16D0])
  allocate(pv(g%ni,g%nj,nl),psi(g%ni,g%nj,nl),omega(g%ni,g%nj,nl),halo(g%ni+2,g%nj+2,nl))
  call setup_basin_mass(cache,d,g,s,error)
  if (fixture == 'islands') then
    call require(index(error,'requires one coastline') > 0,'basin constraint accepted an island case')
    expected_area = 777*h*h
  else
    call require(len_trim(error) == 0,error)
    do j = 1,n
      x = real(j-1,real_8)/(n-1)
      do i = 1,n
        y = real(n-i,real_8)/(n-1)
        if (fixture == 'triangle') then
          exact(i,j,:) = gauge + [0.6D0,0.8D0,1.3D0]*x*y*(1-x-y)
          weights(i,j) = triangle_weight(j-1,n-i)
        else
          exact(i,j,:) = gauge + [0.6D0,0.8D0,1.3D0]*x*y*(1-x)*(1-y)
          weights(i,j) = (1.0D0-0.5D0*merge(1,0,j==1.or.j==n))*(1.0D0-0.5D0*merge(1,0,i==1.or.i==n))
        end if
      end do
    end do
    weights = weights*h*h
    expected_area = sum(weights)
    do k = 1,nl
      target(k) = sum(exact(:,:,k)*weights)
    end do

    pv = ieee_value(0.0D0,ieee_quiet_nan)
    do j = 1,g%nj
      gj = g%j0+j
      do i = 1,g%ni
        gi = g%i0+i
        if (g%global%kind(gi,gj) /= fluid_node) cycle
        modal = exact(gi-1,gj,:)+exact(gi+1,gj,:)+exact(gi,gj-1,:)+exact(gi,gj+1,:) &
                -(4.0D0+s%ss)*exact(gi,gj,:)
        call modes_to_z(s,modal,layer_values)
        pv(i,j,:) = layer_values
      end do
    end do
    call invert_basin_mass(cache,pv,target(2:),psi,recovered,error)
    call require(len_trim(error) == 0,error)
    call z_to_modes(s,psi,omega)
    call fluid_integrals(g,omega,sums)
    call require(maxval(abs(sums(2:)-target(2:))) < 1e-11,'nonzero mass targets were not preserved')
    local_error = 0.0D0
    do j = 1,g%nj
      do i = 1,g%ni
        if (g%global%kind(g%i0+i,g%j0+j) == solid_node) cycle
        local_error = max(local_error,maxval(abs(omega(i,j,:)-exact(g%i0+i,g%j0+j,:))))
      end do
    end do
    call mpi_allreduce(local_error,error_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,g%comm,ierr)
    call require(error_max < 1e-10,'manufactured constrained solution error')
    call z_to_modes(s,recovered,modal)
    call require(maxval(abs(modal-gauge)) < 1e-10,'incorrect modal wall constants')
    if (g%rank == 0) print '(a,es14.6)', 'Max constrained modal solution error: ', error_max

    allocate(wall_values(0:0,nl))
    wall_values(0,:) = recovered
    omega = 1.2D0
    psi = -42.0D0
    call set_wall_pv(g,s,wall_values,omega,psi)
    layer_values = 1.2D0 + [s%s2(1)*(recovered(2)-recovered(1)), &
                           s%s1(2)*(recovered(1)-recovered(2))+s%s2(2)*(recovered(3)-recovered(2)), &
                           s%s1(3)*(recovered(2)-recovered(3))]
    do j = 1,g%nj
      do i = 1,g%ni
        kind = g%global%kind(g%i0+i,g%j0+j)
        if (kind == wall_node) then
          call require(maxval(abs(psi(i,j,:)-layer_values)) < 1e-14,'incorrect wall PV stretching')
        else
          call require(all(psi(i,j,:) == -42.0D0),'wall PV overwrote a non-wall node')
        end if
      end do
    end do

    if (fixture == 'rectangle') then
      call allocate(old_psi,d)
      call allocate(old_zeta,d)
      old_zeta%sdata = pv
      where (.not. ieee_is_finite(old_zeta%sdata)) old_zeta%sdata = 0.0D0
      call initialise(bc,BC_FREE_SLIP)
      call allocate(legacy,d,s,bc &
#ifdef PETSC
                    ,petsc_rtol=5.0D-14 &
#endif
                    )
      call pv_inversion(s,bc,old_zeta,legacy,old_psi)
      call invert_basin_mass(cache,pv,[0.0D0,0.0D0],psi,recovered,error)
      call require(len_trim(error) == 0,error)
      call z_to_modes(s,psi,omega)
      gauge_error = mean(d,omega(:,:,1))
      local_error = maxval(abs(psi-gauge_error-old_psi%sdata))
      call mpi_allreduce(local_error,error_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,g%comm,ierr)
      call require(error_max < 1e-11,'legacy rectangle mismatch after gauge alignment at matched solver accuracy')
      if (g%rank == 0) print '(a,es14.6)', 'Max gauge-aligned rectangle difference: ', error_max
      call deallocate(legacy)
      call deallocate(old_psi)
      call deallocate(old_zeta)
    end if
    call release_basin_mass(cache)
  end if

  halo = ieee_value(0.0D0,ieee_quiet_nan)
  do j = 1,g%nj+2
    gj = g%j0+j-1
    do i = 1,g%ni+2
      gi = g%i0+i-1
      if (min(gi,gj) < 1 .or. max(gi,gj) > n) cycle
      if (g%global%kind(gi,gj) == solid_node) cycle
      halo(i,j,:) = 3.0D0*(gj-1)-2.0D0*(n-gi)+[2.0D0,4.0D0,6.0D0]
    end do
  end do
  call irregular_energy_totals(g,s,halo,area,kinetic,potential)
  expected_potential = 2*area*sum(s%h(:nl-1)*s%s2(:nl-1))/sum(s%h)
  call require(ieee_is_finite(potential).and.all(ieee_is_finite(kinetic)),'energy read a nonfluid sample')
  call require(abs(area-expected_area) < 1e-14,'unexpected diagnostic fluid area')
  call require(maxval(abs(kinetic-6.5D0*area*s%h/sum(s%h))) < 1e-12,'affine kinetic energy mismatch')
  call require(abs(potential-expected_potential) < 1e-12,'potential energy mismatch')
  if (g%rank == 0) print '(a,5(1x,es14.6),1x,a)', 'Area, kinetic totals and potential total:', area,kinetic,potential,'PASS'
  call deallocate(s)
  call deallocate(d)
  call parallel_finalise()

contains
  real(real_8) function triangle_weight(ix,iy)
    integer, intent(in) :: ix,iy
    integer :: q

    q = ix+iy
    triangle_weight = 0.0D0
    if (q > n-1) return
    if (q == 0) then
      triangle_weight = 0.25D0
    else if (q == n-1 .and. min(ix,iy) == 0) then
      triangle_weight = 0.125D0
    else if (q == n-1 .or. min(ix,iy) == 0) then
      triangle_weight = 0.5D0
    else
      triangle_weight = 1.0D0
    end if
  end function triangle_weight

  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message
    if (.not. condition) call pabort(message,__FILE__,__LINE__)
  end subroutine require
end program test_irregular_mass
