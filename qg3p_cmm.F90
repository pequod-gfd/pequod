#include "config.h"
#ifdef PETSC
#include "petscconf.h"
#endif

! Collocated Dirichlet capacitance correction for the native L = Delta - ss.
module qg_cmm
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use lapack_interfaces, only: dpotrf, dpotrs, dpocon
  use qg_geometry, only: wall_node, fluid_node
  use qg_irregular_parallel, only: distributed_geometry, broadcast_error
  use qg_laplace_bc, only: add_laplace_bc_rhs, apply_dirichlet_bc
  use qg_parallel, only: decomposition, request_cache, parallel_data, allocate, deallocate, update, flush
  use qg_parameters, only: real_8, pi
  use qg_cmm_dense, only: has_scalapack
#if defined(PETSC) && defined(PETSC_HAVE_SCALAPACK)
  use qg_cmm_dense, only: dense_factor, create_factor, set_factor_column, assemble_factor, factor_dense, &
                          solve_dense, destroy_factor
#endif
#ifdef PETSC
  use qg_helmholtz_petsc, only: rectangle_cache => helmholtz_matrices, allocate, deallocate, solve_helmholtz
#else
  use qg_helmholtz, only: rectangle_cache => helm_cache, allocate, deallocate, helm
#endif
  implicit none
  include 'mpif.h'
  private
  public :: cmm_cache, cmm_statistics, setup_cmm, release_cmm, cmm_solve, cmm_get_statistics, cmm_copy_matrix

  integer, parameter :: wide = selected_int_kind(18)
  real(real_8), parameter :: rectangle_tolerance = 5.0D-13

  type cmm_statistics
    integer :: constraints = 0, factors = 0, background_batches = 0
    logical :: distributed_factor = .false.
    integer(wide) :: factor_bytes = 0, local_green_bytes = 0
    integer(wide) :: explicit_array_bytes = 0, dense_setup_peak_bytes = 0
    real(real_8) :: max_rectangle_backward_error = 0.0D0
    real(real_8) :: relative_asymmetry = 0.0D0, min_rcond = 1.0D0
  end type cmm_statistics

  type cmm_cache
    private
    type(rectangle_cache) :: rectangle
    type(decomposition), pointer :: decomp => null()
    type(distributed_geometry), pointer :: geom => null()
    type(parallel_data) :: halo
    type(request_cache), pointer :: requests => null()
    type(cmm_statistics) :: stats
#if defined(PETSC) && defined(PETSC_HAVE_SCALAPACK)
    type(dense_factor), allocatable :: distributed(:)
#endif
    real(real_8), allocatable :: shifts(:), factor(:,:,:), green(:,:,:,:)
    real(real_8), allocatable :: work(:,:,:), lifted_rhs(:,:,:), trace(:,:), local_trace(:,:)
    integer, allocatable :: local_ij(:,:), wall_component(:), factor_of_mode(:)
  end type cmm_cache

contains

  ! decomp and geom are immutable over this cache's lifetime.
  subroutine setup_cmm(cache, decomp, geom, shifts, error, retain_green, max_dense_bytes, distributed_factor)
    type(cmm_cache), intent(inout) :: cache
    type(decomposition), target, intent(in) :: decomp
    type(distributed_geometry), target, intent(in) :: geom
    real(real_8), intent(in) :: shifts(:)
    character(*), intent(out) :: error
    logical, optional, intent(in) :: retain_green,distributed_factor
    integer(wide), optional, intent(in) :: max_dense_bytes
    integer :: m, nf, l, k, i, j, n, mode, ierr, compare, invalid, any_invalid
    integer, allocatable :: representative(:), iwork(:)
    real(real_8), allocatable :: condition_work(:)
    real(real_8) :: root_shifts(size(shifts))
    integer(wide) :: budget,root_budget,dense_peak,green_bytes,local_points
    integer(wide), parameter :: real_bytes = storage_size(0.0D0)/8, integer_bytes = storage_size(0)/8
    integer :: settings(3), factor_map(size(shifts))
    logical :: keep_fields,use_distributed
#if defined(PETSC) && defined(PETSC_HAVE_SCALAPACK)
    integer, allocatable :: owned_rows(:)
    integer(wide) :: local_factor_bytes
    real(real_8) :: base_eigenvalue,condition_bound,asymmetry
#endif

    error = ''
    invalid = 0
    if (associated(cache%decomp)) invalid = 1
    if (decomp%periodic_j .or. decomp%layers /= size(shifts) .or. size(shifts) < 1) invalid = 1
    if (any(shifts < 0.0D0) .or. .not. all(ieee_is_finite(shifts))) invalid = 1
    if (any([decomp%ii_0,decomp%jj_0,decomp%ii,decomp%jj] /= [geom%i0,geom%j0,geom%ni,geom%nj])) invalid = 1
    if (any(shape(geom%owner) /= [decomp%g_ii,decomp%g_jj])) invalid = 1
    call mpi_comm_compare(geom%comm,decomp%cart_comm,compare,ierr)
    if (compare == MPI_UNEQUAL) invalid = 1
    call mpi_allreduce(invalid,any_invalid,1,MPI_INTEGER,MPI_MAX,geom%comm,ierr)
    if (any_invalid /= 0) then
      error = 'CMM requires a new cache, closed matching embedding and finite nonnegative modal shifts'
      return
    end if

    keep_fields = .false.
    if (present(retain_green)) keep_fields = retain_green
    use_distributed = .false.
    if (present(distributed_factor)) use_distributed = distributed_factor
    budget = huge(budget)
    if (present(max_dense_bytes)) budget = max_dense_bytes
    settings = [size(shifts),merge(1,0,keep_fields),merge(1,0,use_distributed)]
    call mpi_bcast(settings,size(settings),MPI_INTEGER,0,geom%comm,ierr)
    if (any(settings /= [size(shifts),merge(1,0,keep_fields),merge(1,0,use_distributed)])) &
      error = 'CMM settings differ across ranks'
    if (use_distributed .and. .not. has_scalapack) error = 'distributed CMM factors require PETSc ScaLAPACK'
    call broadcast_error(geom%comm,geom%rank,geom%ranks,error)
    if (len_trim(error) /= 0) return
    root_budget = budget
    root_shifts = shifts
    call mpi_bcast(root_budget,1,MPI_INTEGER8,0,geom%comm,ierr)
    call mpi_bcast(root_shifts,size(shifts),MPI_DOUBLE_PRECISION,0,geom%comm,ierr)
    if (root_budget /= budget .or. any(root_shifts /= shifts)) error = 'CMM settings differ across ranks'
    if (budget < 0) error = 'CMM dense memory budget must be nonnegative'
    call broadcast_error(geom%comm,geom%rank,geom%ranks,error)
    if (len_trim(error) /= 0) return

    n = size(shifts)
    allocate(representative(n))
    nf = 0
    do l = 1,n
      do k = 1,nf
        if (shifts(l) == shifts(representative(k))) exit
      end do
      if (k > nf) then
        nf = nf+1
        representative(nf) = l
      end if
      factor_map(l) = k
    end do
    m = count(geom%global%kind(2:decomp%g_ii-1,2:decomp%g_jj-1) == wall_node)
    local_points = int(geom%ni,wide)*geom%nj
    dense_peak = 0
    green_bytes = 0
    if (m > 0 .and. geom%rank == 0 .and. .not. use_distributed) then
      if (int(m,wide) > (huge(dense_peak)-int(m,wide)*(3*real_bytes+integer_bytes))/real_bytes/nf/m) then
        error = 'CMM dense setup byte count overflows'
      else
        dense_peak = real_bytes*int(m,wide)*m*nf+int(m,wide)*(3*real_bytes+integer_bytes)
        if (dense_peak > budget) then
          write(error,'(a,i0,a,i0)') 'CMM dense setup requires ',dense_peak,' bytes, budget ',budget
        end if
      end if
    end if
    if (keep_fields .and. m > 0) then
      if (local_points > (huge(green_bytes)/2)/real_bytes/m/nf) then
        error = 'CMM retained field byte count overflows'
      else
        green_bytes = local_points*m*nf*real_bytes
      end if
    end if
    call broadcast_error(geom%comm,geom%rank,geom%ranks,error)
    if (len_trim(error) /= 0) return

    cache%geom => geom
    cache%decomp => decomp
    cache%shifts = shifts
    cache%factor_of_mode = factor_map
    cache%stats%distributed_factor = use_distributed
    cache%stats%dense_setup_peak_bytes = dense_peak
    cache%stats%explicit_array_bytes = real_bytes*(n+2*local_points*n+int(geom%ni+2,wide)*(geom%nj+2)*n+ &
      2*int(m,wide)*n)+integer_bytes*(int(n,wide)+3*int(m,wide))+green_bytes
    allocate(cache%local_ij(2,m),cache%wall_component(m))
    cache%local_ij = 0
    cache%stats%constraints = m
    k = 0
    do i = 2,decomp%g_ii-1
      do j = 2,decomp%g_jj-1
        if (geom%global%kind(i,j) /= wall_node) cycle
        k = k+1
        cache%wall_component(k) = geom%global%component(i,j)
        if (geom%owner(i,j) == geom%rank) cache%local_ij(:,k) = [i-geom%i0,j-geom%j0]
      end do
    end do

    allocate(cache%requests)
    call allocate(cache%requests,comm=decomp%cart_comm)
    call allocate(cache%halo,decomp,requests=cache%requests)
#ifdef PETSC
    if (decomp%n == 1) then
      call allocate(cache%rectangle,decomp,shifts,ksp_type='preonly',pc_type='lu',quiet=.true.)
    else
      call allocate(cache%rectangle,decomp,shifts,ksp_type='cg',pc_type='sor',atol=1.0D-50,rtol=5.0D-14, &
                    max_its=5000,zero_initial_guess=.true.,quiet=.true.)
    end if
#else
    call allocate(cache%rectangle,decomp,shifts)
#endif
    allocate(cache%work(geom%ni,geom%nj,decomp%layers),cache%lifted_rhs(geom%ni,geom%nj,decomp%layers))
    if (m == 0) return

    allocate(cache%trace(m,n),cache%local_trace(m,n))
#if defined(PETSC) && defined(PETSC_HAVE_SCALAPACK)
    if (use_distributed) then
      allocate(cache%distributed(nf))
      owned_rows = pack([(k,k=1,m)],cache%local_ij(1,:) > 0)
      do l = 1,nf
        call create_factor(cache%distributed(l),geom%comm,geom%ranks,m)
      end do
    end if
#endif
    if (geom%rank == 0 .and. .not. use_distributed) then
      allocate(cache%factor(m,m,nf),condition_work(3*m),iwork(m))
      cache%stats%factor_bytes = int(m,wide)*m*nf*real_bytes
      cache%stats%explicit_array_bytes = cache%stats%explicit_array_bytes+cache%stats%factor_bytes
    end if
    if (keep_fields) then
      allocate(cache%green(geom%ni,geom%nj,m,nf))
      cache%stats%local_green_bytes = green_bytes
    end if

    do k = 1,m
      cache%work = 0.0D0
      i = cache%local_ij(1,k)
      j = cache%local_ij(2,k)
      if (i > 0) cache%work(i,j,representative(:nf)) = 1.0D0
      call rectangle_inverse(cache,error)
      if (len_trim(error) /= 0) exit
      call sample_wall(cache,cache%work)
      if (.not. use_distributed) call reduce_wall(cache)
      do l = 1,nf
        mode = representative(l)
#if defined(PETSC) && defined(PETSC_HAVE_SCALAPACK)
        if (use_distributed) then
          call set_factor_column(cache%distributed(l),k,owned_rows,-cache%local_trace(:,mode))
          if (mod(k,32)==0 .or. k==m) call assemble_factor(cache%distributed(l),k==m)
        end if
#endif
        if (geom%rank == 0 .and. .not. use_distributed) cache%factor(:,k,l) = -cache%trace(:,mode)
        if (keep_fields) cache%green(:,:,k,l) = cache%work(:,:,mode)
      end do
    end do

#if defined(PETSC) && defined(PETSC_HAVE_SCALAPACK)
    if (use_distributed .and. len_trim(error) == 0) then
      base_eigenvalue = 4*(sin(pi/(2*(decomp%g_ii-1)))**2+sin(pi/(2*(decomp%g_jj-1)))**2)
      do l = 1,nf
        mode = representative(l)
        condition_bound = min(1.0D0,base_eigenvalue+shifts(mode))/(8+shifts(mode))
        call factor_dense(cache%distributed(l),condition_bound,cache%stats%max_rectangle_backward_error, &
                          asymmetry,local_factor_bytes,error)
        if (len_trim(error) /= 0) exit
        cache%stats%min_rcond = min(cache%stats%min_rcond,condition_bound)
        cache%stats%relative_asymmetry = max(cache%stats%relative_asymmetry,asymmetry)
        cache%stats%factor_bytes = cache%stats%factor_bytes+local_factor_bytes
      end do
      cache%stats%explicit_array_bytes = cache%stats%explicit_array_bytes+cache%stats%factor_bytes
    end if
#endif
    if (geom%rank == 0 .and. .not. use_distributed .and. len_trim(error) == 0) then
      do l = 1,nf
        call factor_capacitance(cache%factor(:,:,l), cache%stats, condition_work, iwork, error)
        if (len_trim(error) /= 0) exit
      end do
    end if
    call mpi_bcast(error,len(error),MPI_CHARACTER,0,geom%comm,ierr)
    if (len_trim(error) /= 0) then
      call release_cmm(cache)
      return
    end if
    cache%stats%factors = nf
    call mpi_bcast(cache%stats%min_rcond,1,MPI_DOUBLE_PRECISION,0,geom%comm,ierr)
    call mpi_bcast(cache%stats%relative_asymmetry,1,MPI_DOUBLE_PRECISION,0,geom%comm,ierr)
  end subroutine setup_cmm

  subroutine factor_capacitance(c, stats, work, iwork, error)
    real(real_8), intent(inout) :: c(:,:)
    type(cmm_statistics), intent(inout) :: stats
    real(real_8), intent(out) :: work(:)
    integer, intent(out) :: iwork(:)
    character(*), intent(out) :: error
    integer :: m, j, info
    real(real_8) :: anorm, asymmetry, rcond, assembly_error

    error = ''
    m = size(c,1)
    anorm = 0.0D0
    asymmetry = 0.0D0
    do j = 1,m
      anorm = max(anorm,sum(abs(c(:,j))))
      asymmetry = max(asymmetry,maxval(abs(c(:,j)-c(j,:))))
    end do
    asymmetry = asymmetry/max(tiny(anorm),maxval(abs(c)))
    stats%relative_asymmetry = max(stats%relative_asymmetry,asymmetry)
    if (asymmetry > 5.0D-12) then
      write(error,'(a,es12.4)') 'capacitance assembly is asymmetric, relative defect ', asymmetry
      return
    end if
    call dpotrf('L',m,c,m,info)
    if (info /= 0) then
      write(error,'(a,i0)') 'capacitance Cholesky failed at index ', info
      return
    end if
    call dpocon('L',m,c,m,anorm,rcond,work,iwork,info)
    assembly_error = max(asymmetry,stats%max_rectangle_backward_error,epsilon(rcond))
    if (info /= 0 .or. rcond*1.0D-8 < assembly_error) then
      write(error,'(a,es12.4)') 'capacitance assembly accuracy is insufficient for reciprocal condition ', rcond
      return
    end if
    stats%min_rcond = min(stats%min_rcond,rcond)
  end subroutine factor_capacitance

  subroutine rectangle_inverse(cache, error)
    type(cmm_cache), intent(inout) :: cache
    character(*), intent(out) :: error
    real(real_8) :: residual
    integer :: layer

    cache%lifted_rhs = cache%work
#ifdef PETSC
    call solve_helmholtz(cache%rectangle,cache%work)
#else
    call helm(cache%rectangle,cache%work)
#endif
    do layer = 1,size(cache%work,3)
      call apply_dirichlet_bc(cache%decomp,0.0D0,cache%work(:,:,layer))
    end do
    cache%stats%background_batches = cache%stats%background_batches+1
    call measure_residual(cache,cache%work,cache%lifted_rhs,.false.,residual)
    cache%stats%max_rectangle_backward_error = max(cache%stats%max_rectangle_backward_error,residual)
    error = ''
    if (residual > rectangle_tolerance) write(error,'(a,es12.4)') 'rectangle solve backward error ', residual
  end subroutine rectangle_inverse

  subroutine measure_residual(cache, solution, rhs, physical_only, residual)
    type(cmm_cache), intent(inout) :: cache
    real(real_8), intent(in) :: solution(:,:,:), rhs(:,:,:)
    logical, intent(in) :: physical_only
    real(real_8), intent(out) :: residual
    integer :: i, j, l, i0, j0, ierr
    real(real_8) :: local(3,size(rhs,3)), global(3,size(rhs,3)), value, defect

    cache%halo%sdata = solution
    call update(cache%halo)
    call flush(cache%halo)
    local = 0.0D0
    associate(d => cache%decomp, h => cache%halo%data)
      i0 = d%ii_0
      j0 = d%jj_0
      do l = 1,size(rhs,3)
        do j = d%f_j_min,d%f_j_max
          do i = d%f_i_min,d%f_i_max
            if (physical_only) then
              if (cache%geom%global%kind(i+i0-1,j+j0-1) /= fluid_node) cycle
            end if
            value = rhs(i-1,j-1,l)
            defect = h(i-1,j,l)+h(i+1,j,l)+h(i,j-1,l)+h(i,j+1,l)-(4.0D0+cache%shifts(l))*h(i,j,l)-value
            if (.not. ieee_is_finite(defect)) then
              local(1,l) = huge(0.0D0)
            else
              local(1,l) = max(local(1,l),abs(defect))
            end if
            local(2,l) = max(local(2,l),abs(value))
            local(3,l) = max(local(3,l),abs(h(i,j,l)),abs(h(i-1,j,l)),abs(h(i+1,j,l)),abs(h(i,j-1,l)),abs(h(i,j+1,l)))
          end do
        end do
      end do
    end associate
    call mpi_allreduce(local,global,size(local),MPI_DOUBLE_PRECISION,MPI_MAX,cache%geom%comm,ierr)
    residual = maxval(global(1,:)/max(tiny(0.0D0),global(2,:)+(8.0D0+cache%shifts)*global(3,:)))
    if (.not. ieee_is_finite(residual)) residual = huge(residual)
  end subroutine measure_residual

  subroutine sample_wall(cache, field)
    type(cmm_cache), intent(inout) :: cache
    real(real_8), intent(in) :: field(:,:,:)
    integer :: k, i, j

    cache%local_trace = 0.0D0
    do k = 1,size(cache%local_ij,2)
      i = cache%local_ij(1,k)
      j = cache%local_ij(2,k)
      if (i > 0) cache%local_trace(k,:) = field(i,j,:)
    end do
  end subroutine sample_wall

  subroutine reduce_wall(cache)
    type(cmm_cache), intent(inout) :: cache
    integer :: ierr

    if (cache%stats%distributed_factor) then
      call mpi_allreduce(cache%local_trace,cache%trace,size(cache%trace),MPI_DOUBLE_PRECISION,MPI_SUM,cache%geom%comm,ierr)
    else
      call mpi_reduce(cache%local_trace,cache%trace,size(cache%trace),MPI_DOUBLE_PRECISION,MPI_SUM,0,cache%geom%comm,ierr)
    end if
  end subroutine reduce_wall

  ! Solve the physical strict-fluid equations with the supplied modal component
  ! constants. The returned fictitious values are part of the embedding extension.
  subroutine cmm_solve(cache, rhs, constants, solution, error)
    type(cmm_cache), intent(inout) :: cache
    real(real_8), intent(in) :: rhs(:,:,:), constants(0:,:)
    real(real_8), intent(out) :: solution(:,:,:)
    character(*), intent(out) :: error
    integer :: l, k, i, j, m, ierr, info, invalid, any_invalid
    real(real_8) :: residual, wall_error, global_wall_error, root_values(size(constants,1),size(constants,2))

    error = ''
    if (.not. associated(cache%decomp)) then
      error = 'CMM cache is not initialized'
      return
    end if
    invalid = 0
    if (any(shape(rhs) /= shape(cache%work)) .or. any(shape(solution) /= shape(cache%work))) invalid = 1
    if (any(shape(constants) /= [size(cache%geom%global%loops),size(cache%shifts)])) invalid = 1
    if (.not. all(ieee_is_finite(constants)) .or. .not. all(ieee_is_finite(rhs))) invalid = 1
    call mpi_allreduce(invalid,any_invalid,1,MPI_INTEGER,MPI_MAX,cache%geom%comm,ierr)
    if (any_invalid /= 0) then
      error = 'CMM requires finite RHS and wall values with the cache dimensions'
      return
    end if
    root_values = constants
    call mpi_bcast(root_values,size(root_values),MPI_DOUBLE_PRECISION,0,cache%geom%comm,ierr)
    if (any(root_values /= constants)) error = 'CMM wall values differ across ranks'
    call broadcast_error(cache%geom%comm,cache%geom%rank,cache%geom%ranks,error)
    if (len_trim(error) /= 0) return

    cache%work = rhs
    do l = 1,size(cache%shifts)
      call add_laplace_bc_rhs(cache%decomp,constants(0,l),cache%work(:,:,l))
    end do
    call rectangle_inverse(cache,error)
    if (len_trim(error) /= 0) return
    solution = cache%work
    do l = 1,size(cache%shifts)
      call apply_dirichlet_bc(cache%decomp,constants(0,l),solution(:,:,l))
    end do
    m = cache%stats%constraints
    if (m == 0) return

    call sample_wall(cache,solution)
    call reduce_wall(cache)
#if defined(PETSC) && defined(PETSC_HAVE_SCALAPACK)
    if (cache%stats%distributed_factor) then
      do l = 1,size(cache%shifts)
        cache%trace(:,l) = cache%trace(:,l)-constants(cache%wall_component,l)
        call solve_dense(cache%distributed(cache%factor_of_mode(l)),cache%trace(:,l))
      end do
    end if
#endif
    if (cache%geom%rank == 0 .and. .not. cache%stats%distributed_factor) then
      do l = 1,size(cache%shifts)
        cache%trace(:,l) = cache%trace(:,l)-constants(cache%wall_component,l)
        call dpotrs('L',m,1,cache%factor(:,:,cache%factor_of_mode(l)),m,cache%trace(:,l:l),m,info)
        if (info /= 0) then
          write(error,'(a,i0)') 'capacitance factor solve failed at index ', info
          exit
        end if
      end do
    end if
    call mpi_bcast(error,len(error),MPI_CHARACTER,0,cache%geom%comm,ierr)
    if (len_trim(error) /= 0) return
    if (.not. cache%stats%distributed_factor) &
      call mpi_bcast(cache%trace,size(cache%trace),MPI_DOUBLE_PRECISION,0,cache%geom%comm,ierr)

    if (allocated(cache%green)) then
      do l = 1,size(cache%shifts)
        do k = 1,m
          solution(:,:,l) = solution(:,:,l)+cache%trace(k,l)*cache%green(:,:,k,cache%factor_of_mode(l))
        end do
      end do
    else
      cache%work = 0.0D0
      do k = 1,m
        i = cache%local_ij(1,k)
        j = cache%local_ij(2,k)
        if (i > 0) cache%work(i,j,:) = cache%trace(k,:)
      end do
      call rectangle_inverse(cache,error)
      if (len_trim(error) /= 0) return
      solution = solution+cache%work
    end if

    wall_error = 0.0D0
    do k = 1,m
      i = cache%local_ij(1,k)
      j = cache%local_ij(2,k)
      if (i == 0) cycle
      wall_error = max(wall_error,maxval(abs(solution(i,j,:)-constants(cache%wall_component(k),:))))
    end do
    call mpi_allreduce(wall_error,global_wall_error,1,MPI_DOUBLE_PRECISION,MPI_MAX,cache%geom%comm,ierr)
    call measure_residual(cache,solution,rhs,.true.,residual)
    if (residual > 1.0D-11 .or. global_wall_error > 1.0D-11*max(1.0D0,maxval(abs(constants)))) then
      write(error,'(a,2(1x,es12.4))') 'CMM physical PDE backward error and absolute wall trace error:', residual, global_wall_error
    end if
  end subroutine cmm_solve

  function cmm_get_statistics(cache) result(stats)
    type(cmm_cache), intent(in) :: cache
    type(cmm_statistics) :: stats

    stats = cache%stats
  end function cmm_get_statistics

  ! Root-only diagnostic reconstruction from the stored triangular factor.
  ! The runtime solve never materializes this second full matrix.
  subroutine cmm_copy_matrix(cache, mode, matrix)
    type(cmm_cache), intent(in) :: cache
    integer, intent(in) :: mode
    real(real_8), allocatable, intent(out) :: matrix(:,:)
    integer :: m, i, j, k

    m = cache%stats%constraints
    allocate(matrix(0,0))
    if (.not. allocated(cache%factor)) return
    deallocate(matrix)
    allocate(matrix(m,m))
    k = cache%factor_of_mode(mode)
    do j = 1,m
      do i = j,m
        matrix(i,j) = dot_product(cache%factor(i,:j,k),cache%factor(j,:j,k))
        matrix(j,i) = matrix(i,j)
      end do
    end do
  end subroutine cmm_copy_matrix

  subroutine release_cmm(cache)
    type(cmm_cache), intent(inout) :: cache
#if defined(PETSC) && defined(PETSC_HAVE_SCALAPACK)
    integer :: l
#endif

    if (.not. associated(cache%decomp)) return
    call deallocate(cache%rectangle)
    call deallocate(cache%halo)
    call deallocate(cache%requests)
    deallocate(cache%requests)
    deallocate(cache%work,cache%lifted_rhs,cache%shifts,cache%factor_of_mode,cache%local_ij,cache%wall_component)
    if (allocated(cache%factor)) deallocate(cache%factor)
#if defined(PETSC) && defined(PETSC_HAVE_SCALAPACK)
    if (allocated(cache%distributed)) then
      do l = 1,size(cache%distributed)
        call destroy_factor(cache%distributed(l))
      end do
      deallocate(cache%distributed)
    end if
#endif
    if (allocated(cache%green)) deallocate(cache%green)
    if (allocated(cache%trace)) deallocate(cache%trace,cache%local_trace)
    nullify(cache%decomp,cache%geom)
    cache%stats = cmm_statistics()
  end subroutine release_cmm

end module qg_cmm
