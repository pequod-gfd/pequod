! Experimental static circulation/mass system. The weak circulation measure
! uses prescribed stage means. Forced/no-slip wall evolution is not here.
module qg_island_constraints
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use lapack_interfaces, only: dgetrf, dgetrs, dgecon
  use qg_cmm, only: cmm_cache, setup_cmm, cmm_solve, release_cmm
  use qg_constraints, only: fluid_integrals
  use qg_diagnostics, only: stratification, z_to_modes, modes_to_z
  use qg_geometry, only: wall_node
  use qg_irregular_parallel, only: distributed_geometry, broadcast_error
  use qg_irregular_operators, only: wall_mean_plan, setup_wall_means, reconstruct_wall_means, &
                                    stiffness_flux, sum_open_flux, mean_to_pv_nodes
  use qg_parallel, only: decomposition, parallel_data, request_cache, allocate, deallocate, update, flush
  use qg_parameters, only: real_8
  implicit none
  include 'mpif.h'
  private
  public :: weak_island_cache, setup_weak_islands, invert_weak_islands, release_weak_islands

  type weak_island_cache
    private
    type(distributed_geometry), pointer :: geom => null()
    type(stratification), pointer :: strat => null()
    type(cmm_cache) :: cmm
    type(wall_mean_plan) :: mean
    type(parallel_data) :: halo
    type(request_cache), pointer :: requests => null()
    real(real_8), allocatable :: harmonic(:,:,:,:), q(:,:,:), rhs(:,:,:), phi(:,:,:), means(:,:,:)
    real(real_8), allocatable :: open(:,:), integrated(:,:,:), matrix(:,:), factor(:,:), row_scale(:)
    integer, allocatable :: pivots(:)
    real(real_8) :: rcond
  end type weak_island_cache

contains

  subroutine setup_weak_islands(cache,decomp,geom,strat,error)
    type(weak_island_cache), intent(inout) :: cache
    type(decomposition), target, intent(in) :: decomp
    type(distributed_geometry), target, intent(in) :: geom
    type(stratification), target, intent(in) :: strat
    character(*), intent(out) :: error
    real(real_8), allocatable :: wall(:,:), col(:), work(:)
    integer, allocatable :: iwork(:)
    integer :: nk,nl,n,i,j,k,mode,col_index,info,ierr
    real(real_8) :: matrix_norm

    error = ''
    if (associated(cache%geom)) error = 'weak island cache is already initialized'
    if (strat%layers /= decomp%layers .or. strat%ss(1) /= 0.0D0) error = 'invalid weak island modal stratification'
    call broadcast_error(geom%comm,geom%rank,geom%ranks,error)
    if (len_trim(error) /= 0) return
    call setup_wall_means(geom,cache%mean,error)
    if (len_trim(error) /= 0) return
    call setup_cmm(cache%cmm,decomp,geom,strat%ss,error)
    if (len_trim(error) /= 0) return
    cache%geom => geom
    cache%strat => strat
    allocate(cache%requests)
    call allocate(cache%requests,comm=decomp%cart_comm)
    call allocate(cache%halo,decomp,requests=cache%requests)
    nk = ubound(geom%global%loops,1)
    nl = strat%layers
    n = (nk+1)*nl
    i = geom%ni
    j = geom%nj
    allocate(cache%harmonic(i,j,nl,0:nk),cache%q(i,j,nl),cache%rhs(i,j,nl),cache%phi(i,j,nl), &
             cache%means(i,j,nl),cache%open(nl,size(geom%faces)),cache%integrated(i,j,nl))
    allocate(wall(0:nk,nl),col(n))
    if (geom%rank == 0) allocate(cache%matrix(n,n),cache%factor(n,n),cache%row_scale(n),cache%pivots(n),work(4*n),iwork(n))

    cache%rhs = 0.0D0
    do k = 0,nk
      wall = 0.0D0
      wall(k,:) = 1.0D0
      call cmm_solve(cache%cmm,cache%rhs,wall,cache%harmonic(:,:,:,k),error)
      if (len_trim(error) /= 0) exit
      do mode = 1,nl
        wall = 0.0D0
        wall(k,mode) = 1.0D0
        cache%phi = 0.0D0
        cache%phi(:,:,mode) = cache%harmonic(:,:,mode,k)
        call measure_response(cache,cache%phi,wall,col,error)
        if (len_trim(error) /= 0) exit
        col_index = k*nl+mode
        if (geom%rank == 0) cache%matrix(:,col_index) = col
      end do
      if (len_trim(error) /= 0) exit
    end do

    if (geom%rank == 0 .and. len_trim(error) == 0) then
      cache%row_scale = maxval(abs(cache%matrix),dim=2)
      if (any(cache%row_scale <= tiny(1.0D0)) .or. .not. all(ieee_is_finite(cache%matrix))) then
        error = 'weak island influence has a zero or nonfinite constraint row'
      else
        cache%factor = cache%matrix/spread(cache%row_scale,2,n)
        matrix_norm = maxval(sum(abs(cache%factor),dim=2))
        call dgetrf(n,n,cache%factor,n,cache%pivots,info)
        if (info == 0) call dgecon('I',n,cache%factor,n,matrix_norm,cache%rcond,work,iwork,info)
        if (info /= 0) then
          error = 'weak island constraint matrix is rank deficient'
        else if (cache%rcond < 1.0D-12) then
          error = 'weak island constraint matrix is ill-conditioned'
        end if
      end if
    end if
    call mpi_bcast(error,len(error),MPI_CHARACTER,0,geom%comm,ierr)
    if (len_trim(error) /= 0) call release_weak_islands(cache)
  end subroutine setup_weak_islands

  ! The affine circulation from boundary Q is added separately. With
  ! omega^a = Q - (S psi)^a, the psi response is ss*A*psi_mean - B psi.
  subroutine measure_response(cache,phi,wall,response,error)
    type(weak_island_cache), intent(inout) :: cache
    real(real_8), intent(in) :: phi(:,:,:),wall(0:,:)
    real(real_8), intent(out) :: response(:)
    character(*), intent(out) :: error
    real(real_8) :: mass(size(phi,3)), local(size(response)), area
    integer :: i,j,m,k,nl,nk,ierr

    associate(g => cache%geom)
      nl = size(phi,3)
      nk = ubound(wall,1)
      cache%halo%sdata = phi
      ! Impose the prescribed trace before evaluating the physical flux.
      do j = 1,g%nj
        do i = 1,g%ni
          k = g%global%component(g%i0+i,g%j0+j)
          if (k >= 0) cache%halo%sdata(i,j,:) = wall(k,:)
        end do
      end do
      call update(cache%halo)
      call flush(cache%halo)
      call stiffness_flux(g,cache%halo%data,cache%open,error)
      if (len_trim(error) /= 0) return
      call sum_open_flux(g,cache%open,cache%integrated)
      call reconstruct_wall_means(cache%mean,cache%halo%sdata,wall,cache%means)
      local = 0.0D0
      do j = 1,g%nj
        do i = 1,g%ni
          k = g%global%component(g%i0+i,g%j0+j)
          if (k <= 0) cycle
          area = g%area(i,j)/g%h**2
          do m = 1,nl
            local((k-1)*nl+m) = local((k-1)*nl+m)+area*cache%strat%ss(m)*cache%means(i,j,m)-cache%integrated(i,j,m)
          end do
        end do
      end do
      call mpi_allreduce(local,response,size(local),MPI_DOUBLE_PRECISION,MPI_SUM,g%comm,ierr)
      call fluid_integrals(g,cache%halo%sdata,mass)
      response(nk*nl+1:size(response)-1) = mass(2:)
      response(size(response)) = wall(0,1)
    end associate
  end subroutine measure_response

  subroutine invert_weak_islands(cache,comm,pv_means,gamma_target,mass_target,psi,wall_layer,backward_error,error)
    type(weak_island_cache), intent(inout) :: cache
    integer, intent(in) :: comm
    real(real_8), intent(in) :: pv_means(:,:,:), gamma_target(:,:), mass_target(:)
    real(real_8), intent(out) :: psi(:,:,:),wall_layer(0:,:),backward_error
    character(*), intent(out) :: error
    real(real_8), allocatable :: wall(:,:), target(:), replicated(:), affine(:), local(:), response(:), coefficients(:,:)
    real(real_8) :: row_error,scale,area
    integer :: i,j,k,m,n,nl,nk,ierr,invalid,global_invalid,info,comparison,counts(2)

    error = ''
    invalid = 0
    if (.not. associated(cache%geom)) then
      invalid = 1
    else
      call mpi_comm_compare(comm,cache%geom%comm,comparison,ierr)
      nl = cache%strat%layers
      nk = ubound(cache%geom%global%loops,1)
      if (comparison == MPI_UNEQUAL .or. any(shape(pv_means) /= shape(cache%q)) .or. &
          any(shape(psi) /= shape(cache%q)) .or. any(shape(gamma_target) /= [nk,nl]) .or. &
          size(mass_target) /= nl-1 .or. any(shape(wall_layer) /= [nk+1,nl])) invalid = 1
      if (.not. all(ieee_is_finite(gamma_target)) .or. .not. all(ieee_is_finite(mass_target))) invalid = 1
    end if
    call mpi_allreduce(invalid,global_invalid,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (global_invalid /= 0) then
      error = 'weak island inversion requires a live cache, matching arrays and finite targets'
      return
    end if
    n = (nk+1)*nl
    call mpi_allreduce([n,-n],counts,2,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (counts(1) /= -counts(2)) then
      error = 'weak island constraint counts differ across ranks'
      return
    end if
    allocate(wall(0:nk,nl),target(n),replicated(n),affine(n),local(n),response(n),coefficients(n,1))
    target = [reshape(gamma_target,[nk*nl]),mass_target,0.0D0]
    replicated = target
    call mpi_bcast(replicated,n,MPI_DOUBLE_PRECISION,0,comm,ierr)
    invalid = merge(1,0,any(replicated /= target))
    call mpi_allreduce(invalid,global_invalid,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (global_invalid /= 0) then
      error = 'weak island targets differ across ranks'
      return
    end if

    associate(g => cache%geom, s => cache%strat)
      do k = 1,nk
        call z_to_modes(s,gamma_target(k,:),target((k-1)*nl+1:k*nl))
      end do
      cache%rhs = 0.0D0
      do j = 1,g%nj
        do i = 1,g%ni
          if (g%area(i,j) > 0.0D0) cache%rhs(i,j,:) = pv_means(i,j,:)
        end do
      end do
      call z_to_modes(s,cache%rhs,cache%q)
      local = 0.0D0
      do j = 1,g%nj
        do i = 1,g%ni
          k = g%global%component(g%i0+i,g%j0+j)
          if (k <= 0) cycle
          area = g%area(i,j)/g%h**2
          local((k-1)*nl+1:k*nl) = local((k-1)*nl+1:k*nl)+area*cache%q(i,j,:)
        end do
      end do
      call mpi_allreduce(local,affine,n,MPI_DOUBLE_PRECISION,MPI_SUM,g%comm,ierr)
      if (.not. all(ieee_is_finite(affine))) then
        error = 'nonfinite stage PV mean on an island cut cell'
        return
      end if
      call mean_to_pv_nodes(g,cache%q,cache%rhs)
      wall = 0.0D0
      call cmm_solve(cache%cmm,cache%rhs,wall,cache%phi,error)
      if (len_trim(error) /= 0) return
      call measure_response(cache,cache%phi,wall,response,error)
      if (len_trim(error) /= 0) return
      if (g%rank == 0) then
        coefficients(:,1) = (target-affine-response)/cache%row_scale
        call dgetrs('N',n,1,cache%factor,n,cache%pivots,coefficients,n,info)
        if (info /= 0 .or. .not. all(ieee_is_finite(coefficients))) error = 'weak island LU solve failed'
      end if
      call mpi_bcast(error,len(error),MPI_CHARACTER,0,g%comm,ierr)
      if (len_trim(error) /= 0) return
      call mpi_bcast(coefficients,n,MPI_DOUBLE_PRECISION,0,g%comm,ierr)
      do k = 0,nk
        wall(k,:) = coefficients(k*nl+1:(k+1)*nl,1)
        do m = 1,nl
          cache%phi(:,:,m) = cache%phi(:,:,m)+wall(k,m)*cache%harmonic(:,:,m,k)
        end do
      end do
      call measure_response(cache,cache%phi,wall,response,error)
      if (len_trim(error) /= 0) return
      backward_error = 0.0D0
      do k = 1,n
        row_error = abs(response(k)+affine(k)-target(k))
        scale = max(1.0D0,abs(response(k))+abs(affine(k))+abs(target(k)))
        backward_error = max(backward_error,row_error/scale)
      end do
      if (g%rank == 0) then
        if (backward_error > min(1.0D-10,1.0D-8*cache%rcond)) error = 'weak island constraint residual exceeds tolerance'
      end if
      call mpi_bcast(error,len(error),MPI_CHARACTER,0,g%comm,ierr)
      if (len_trim(error) /= 0) return
      ! measure_response stamped the exact component constants in halo.
      call modes_to_z(s,cache%halo%sdata,psi)
      do k = 0,nk
        call modes_to_z(s,wall(k,:),wall_layer(k,:))
      end do
    end associate
  end subroutine invert_weak_islands

  subroutine release_weak_islands(cache)
    type(weak_island_cache), intent(inout) :: cache

    if (.not. associated(cache%geom)) return
    call release_cmm(cache%cmm)
    call deallocate(cache%halo)
    call deallocate(cache%requests)
    deallocate(cache%requests)
    deallocate(cache%harmonic,cache%q,cache%rhs,cache%phi,cache%means,cache%open,cache%integrated)
    if (allocated(cache%matrix)) deallocate(cache%matrix,cache%factor,cache%row_scale,cache%pivots)
    nullify(cache%geom,cache%strat)
  end subroutine release_weak_islands

end module qg_island_constraints
