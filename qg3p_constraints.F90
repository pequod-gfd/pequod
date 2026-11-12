! Simply connected modal mass constraints for the irregular Dirichlet inverse.
module qg_constraints
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use qg_cmm, only: cmm_cache, setup_cmm, cmm_solve, release_cmm
  use qg_diagnostics, only: stratification, z_to_modes, modes_to_z
  use qg_geometry, only: fluid_node, wall_node
  use qg_irregular_parallel, only: distributed_geometry
  use qg_parallel, only: decomposition
  use qg_parameters, only: real_8
  implicit none
  include 'mpif.h'
  private
  public :: basin_mass_cache, setup_basin_mass, invert_basin_mass, release_basin_mass, fluid_integrals

  type basin_mass_cache
    private
    type(distributed_geometry), pointer :: geom => null()
    type(stratification), pointer :: strat => null()
    type(cmm_cache) :: cmm
    real(real_8), allocatable :: harmonic(:,:,:), modal_rhs(:,:,:), particular(:,:,:), response(:)
    real(real_8) :: area
  end type basin_mass_cache

contains

  ! Owned-node lumped area quadrature; ignore fictitious values even if NaN.
  subroutine fluid_integrals(geom, field, integrals)
    type(distributed_geometry), intent(in) :: geom
    real(real_8), intent(in) :: field(:,:,:)
    real(real_8), intent(out) :: integrals(size(field,3))
    real(real_8) :: local(size(field,3))
    integer :: i, j, ierr

    local = 0.0D0
    do j = 1,geom%nj
      do i = 1,geom%ni
        if (geom%area(i,j) > 0.0D0) local = local+geom%area(i,j)*field(i,j,:)
      end do
    end do
    call mpi_allreduce(local,integrals,size(local),MPI_DOUBLE_PRECISION,MPI_SUM,geom%comm,ierr)
  end subroutine fluid_integrals

  subroutine setup_basin_mass(cache, decomp, geom, strat, error)
    type(basin_mass_cache), intent(inout) :: cache
    type(decomposition), target, intent(in) :: decomp
    type(distributed_geometry), target, intent(in) :: geom
    type(stratification), target, intent(in) :: strat
    character(*), intent(out) :: error
    integer :: nl, ierr, invalid, global_invalid
    real(real_8) :: local_area

    ! geom and strat are successfully constructed collective inputs.
    invalid = 0
    if (size(geom%global%loops) /= 1) invalid = 2
    if (strat%ss(1) /= 0.0D0) invalid = max(invalid,1)
    call mpi_allreduce(invalid,global_invalid,1,MPI_INTEGER,MPI_MAX,geom%comm,ierr)
    select case(global_invalid)
      case(2)
        error = 'simply connected mass inverter requires one coastline'
        return
      case(1)
        error = 'mass inverter requires a barotropic first mode'
        return
    end select
    call setup_cmm(cache%cmm,decomp,geom,strat%ss,error)
    if (len_trim(error) /= 0) return
    cache%geom => geom
    cache%strat => strat
    nl = strat%layers
    allocate(cache%harmonic(geom%ni,geom%nj,nl),cache%modal_rhs(geom%ni,geom%nj,nl), &
             cache%particular(geom%ni,geom%nj,nl),cache%response(nl))
    cache%modal_rhs = 0.0D0
    call cmm_solve(cache%cmm,cache%modal_rhs,spread([1.0D0],2,nl),cache%harmonic,error)
    if (len_trim(error) /= 0) then
      call release_basin_mass(cache)
      return
    end if
    local_area = sum(geom%area)
    call mpi_allreduce(local_area,cache%area,1,MPI_DOUBLE_PRECISION,MPI_SUM,geom%comm,ierr)
    call fluid_integrals(geom,cache%harmonic,cache%response)
    if (cache%area <= 0.0D0 .or. any(cache%response <= 0.0D0) .or. .not. all(ieee_is_finite(cache%response))) then
      error = 'invalid mass response of a unit-wall harmonic field'
      call release_basin_mass(cache)
    end if
  end subroutine setup_basin_mass

  ! All ranks call on a common intracommunicator with the same group as the cache.
  ! Pass comm explicitly so even a rank without a live cache can participate
  ! in validation before any rank touches the cache's solver communicator.
  subroutine invert_basin_mass(cache, comm, zeta, target_mass, psi, wall_values, error)
    type(basin_mass_cache), intent(inout) :: cache
    integer, intent(in) :: comm
    real(real_8), intent(in) :: zeta(:,:,:), target_mass(:) ! modal totals, modes 2..L
    real(real_8), intent(out) :: psi(:,:,:), wall_values(:) ! layer streamfunction
    character(*), intent(out) :: error
    real(real_8) :: modal_constants(size(zeta,3)), mass(size(zeta,3)), scale
    real(real_8) :: extremes(2*size(target_mass))
    integer :: i, j, k, invalid, global_invalid, comparison, ierr, counts(2)

    error = ''
    invalid = 0
    if (.not. associated(cache%geom)) then
      invalid = 4
    else
      call mpi_comm_compare(comm,cache%geom%comm,comparison,ierr)
      if (comparison == MPI_UNEQUAL) then
        invalid = 3
      else if (any(shape(zeta) /= shape(cache%modal_rhs)) .or. any(shape(psi) /= shape(cache%modal_rhs)) .or. &
               size(target_mass) /= cache%strat%layers-1 .or. size(wall_values) /= cache%strat%layers) then
        invalid = 2
      else if (.not. all(ieee_is_finite(target_mass))) then
        invalid = 1
      end if
    end if
    call mpi_allreduce(invalid,global_invalid,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    select case(global_invalid)
      case(4)
        error = 'basin mass cache is not initialized'
      case(3)
        error = 'basin mass cache belongs to a different communicator group'
      case(2)
        error = 'incompatible simply connected mass inversion dimensions'
      case(1)
        error = 'baroclinic mass targets must be finite'
    end select
    if (global_invalid /= 0) return

    k = size(target_mass)
    call mpi_allreduce([k,-k],counts,2,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (counts(1) /= -counts(2)) then
      error = 'incompatible simply connected mass inversion dimensions'
      return
    end if
    if (k > 0) then
      call mpi_allreduce([target_mass,-target_mass],extremes,2*k,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      if (any(extremes(:k) /= -extremes(k+1:))) then
        error = 'baroclinic mass targets differ across ranks'
        return
      end if
    end if
    associate(g => cache%geom, s => cache%strat)
      cache%particular = 0.0D0
      do j = 1,g%nj
        do i = 1,g%ni
          if (g%global%kind(g%i0+i,g%j0+j) == fluid_node) cache%particular(i,j,:) = zeta(i,j,:)
        end do
      end do
      call z_to_modes(s,cache%particular,cache%modal_rhs)
      modal_constants = 0.0D0
      call cmm_solve(cache%cmm,cache%modal_rhs,reshape(modal_constants,[1,s%layers]),cache%particular,error)
      if (len_trim(error) /= 0) return

      call fluid_integrals(g,cache%particular,mass)
      modal_constants(2:) = (target_mass-mass(2:))/cache%response(2:)
      if (.not. all(ieee_is_finite(modal_constants))) then
        error = 'nonfinite baroclinic wall value from mass inversion'
        return
      end if
      do k = 2,s%layers
        cache%particular(:,:,k) = cache%particular(:,:,k)+modal_constants(k)*cache%harmonic(:,:,k)
      end do
      do j = 1,g%nj
        do i = 1,g%ni
          if (g%global%kind(g%i0+i,g%j0+j) == wall_node) cache%particular(i,j,:) = modal_constants
        end do
      end do
      call fluid_integrals(g,cache%particular,mass)
      do k = 2,s%layers
        scale = max(g%h**2,abs(target_mass(k-1)),abs(modal_constants(k))*cache%response(k))
        if (.not. ieee_is_finite(mass(k)) .or. abs(mass(k)-target_mass(k-1)) > 1.0D-11*scale) then
          error = 'baroclinic mass residual exceeds the inversion tolerance'
          return
        end if
      end do
      call modes_to_z(s,cache%particular,psi)
      call modes_to_z(s,modal_constants,wall_values)
    end associate
  end subroutine invert_basin_mass

  subroutine release_basin_mass(cache)
    type(basin_mass_cache), intent(inout) :: cache

    if (.not. associated(cache%geom)) return
    call release_cmm(cache%cmm)
    deallocate(cache%harmonic,cache%modal_rhs,cache%particular,cache%response)
    nullify(cache%geom,cache%strat)
  end subroutine release_basin_mass

end module qg_constraints
