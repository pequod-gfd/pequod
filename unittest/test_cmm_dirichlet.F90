#include "config.h"
program test_cmm_dirichlet
  use qg_cmm
  use qg_cmm_dense, only: has_scalapack
  use qg_geometry, only: solid_node, wall_node
  use qg_irregular_parallel, only: distributed_geometry, load_distributed_geometry
  use qg_parallel, only: decomposition, parallel_init, parallel_finalise, allocate, deallocate, pabort
  use qg_parameters, only: real_8
  implicit none
  include 'mpif.h'

  integer, parameter :: wide = selected_int_kind(18)
  type(decomposition), target :: decomp
  type(distributed_geometry), target :: geom
  type(cmm_cache) :: cache, full_cache
  type(cmm_statistics) :: stats, full_stats
  integer, allocatable :: iis(:),jjs(:)
  real(real_8), allocatable :: rhs(:,:,:), solution(:,:,:), full_solution(:,:,:), constants(:,:), matrix(:,:)
  real(real_8), allocatable :: send(:,:,:), receive(:,:,:), first(:,:,:)
  integer :: n, rank, ranks, ierr, grid(2), i,j,gi,gj,mode,trial,component,unit,mode_count,comm
  real(real_8) :: shifts(3),wrong(3), comparison, maximum, step, x, y, p, px, py
  integer(wide) :: storage_totals(3),local_bytes(3)
  character(1024) :: path, output, error
  character(30) :: arg, fixture
  logical :: reordered,distributed
  real(real_8) :: setup_time(2),solve_time(2),max_times(4),tick

  if (command_argument_count() < 4) then
    print '(a)', 'CMM fixture requires geometry, output, dimension and case name'
    stop 2
  end if
  call get_command_argument(1,path)
  call get_command_argument(2,output)
  call get_command_argument(3,arg)
  read(arg,*) n
  call get_command_argument(4,fixture)
  call parallel_init()
  call mpi_comm_rank(MPI_COMM_WORLD,rank,ierr)
  call mpi_comm_size(MPI_COMM_WORLD,ranks,ierr)
  grid = 0
  call mpi_dims_create(ranks,2,grid,ierr)
  iis = [(i*n/grid(1)-(i-1)*n/grid(1),i=1,grid(1))]
  jjs = [(j*n/grid(2)-(j-1)*n/grid(2),j=1,grid(2))]
  mode_count = size(shifts)
  call allocate(decomp,iis=iis,jjs=jjs,layers=mode_count)
  comm = decomp%cart_comm
  reordered = .false.
  distributed = .false.
  if (command_argument_count() >= 5) then
    call get_command_argument(5,arg)
    reordered = arg == 'reordered' .or. arg == 'distributed-reordered'
    distributed = arg == 'distributed' .or. arg == 'distributed-reordered'
  end if
  if (reordered) call mpi_comm_split(decomp%cart_comm,0,ranks-1-rank,comm,ierr)
  call load_distributed_geometry(trim(path),n,n,[decomp%ii_0,decomp%jj_0,decomp%ii,decomp%jj], &
                                 comm,1.0D0,geom,error)
  call assert_ok(error)
  shifts = [0.0D0,0.3D0,0.3D0]
  step = 1.0D0/(n-1)
  if (fixture == 'smooth') shifts = [0.0D0,0.7D0,0.7D0]*step**2
  if (any(geom%global%kind(2:n-1,2:n-1)==wall_node)) then
    call setup_cmm(cache,decomp,geom,shifts,error,max_dense_bytes=0_wide)
    if (index(error,'CMM dense setup requires ') /= 1) call pabort('dense memory cap was ignored',__FILE__,__LINE__)
  end if
  if (ranks > 1) then
    wrong = shifts
    if (geom%rank == 1) wrong(2) = wrong(2)+0.1D0
    call setup_cmm(cache,decomp,geom,wrong,error)
    if (error /= 'CMM settings differ across ranks') call pabort('inconsistent shifts were accepted',__FILE__,__LINE__)
  end if
  tick = mpi_wtime()
  call setup_cmm(cache,decomp,geom,shifts,error,distributed_factor=distributed)
  if (distributed .and. .not. has_scalapack) then
    if (error /= 'distributed CMM factors require PETSc ScaLAPACK') call assert_ok('unavailable backend not rejected')
    if (geom%rank == 0) print '(a)', 'SKIP distributed CMM factors require PETSc ScaLAPACK'
    if (reordered) call mpi_comm_free(comm,ierr)
    call deallocate(decomp)
    call parallel_finalise()
    stop
  end if
  setup_time(1) = mpi_wtime()-tick
  call assert_ok(error)
  tick = mpi_wtime()
  call setup_cmm(full_cache,decomp,geom,shifts,error,retain_green=.true.)
  setup_time(2) = mpi_wtime()-tick
  call assert_ok(error)
  solve_time = 0.0D0

  allocate(rhs(geom%ni,geom%nj,mode_count),solution(geom%ni,geom%nj,mode_count), &
           full_solution(geom%ni,geom%nj,mode_count),first(geom%ni,geom%nj,mode_count), &
           constants(0:ubound(geom%global%loops,1),mode_count))
  if (ranks > 1) then
    constants = 0.0D0
    if (geom%rank == 1) constants(0,1) = 1.0D0
    rhs = 0.0D0
    call cmm_solve(cache,rhs,constants,solution,error)
    if (error /= 'CMM wall values differ across ranks') call pabort('inconsistent traces were accepted',__FILE__,__LINE__)
  end if
  allocate(send(n,n,mode_count),receive(n,n,mode_count))
  if (geom%rank == 0) open(newunit=unit,file=trim(output),form='unformatted',access='stream',status='replace')
  do trial = 1,3
    do mode = 1,mode_count
      do component = 0,ubound(constants,1)
        constants(component,mode) = 0.25D0*(1+2*component-mode)
      end do
      do j = 1,geom%nj
        gj = geom%j0+j
        do i = 1,geom%ni
          gi = geom%i0+i
          rhs(i,j,mode) = (mod(23*gi+17*gj+gi*gj*13+5*mode,43)-21)/31.0D0
          if (trial == 2 .and. geom%global%kind(gi,gj) == solid_node) rhs(i,j,mode) = rhs(i,j,mode)+0.75D0
          if (fixture == 'smooth') then
            x = step*(gj-1)
            y = step*(n-gi)
            p = x*y*(1-x-y)
            px = y*(1-2*x-y)
            py = x*(1-x-2*y)
            rhs(i,j,mode) = 2*step**2*(px*px+py*py-2*p*(x+y))-shifts(mode)*p*p
          end if
        end do
      end do
    end do
    if (fixture == 'smooth' .or. trial == 3) constants = 0.0D0
    if (trial == 3) rhs = 0.0D0
    tick = mpi_wtime()
    call cmm_solve(cache,rhs,constants,solution,error)
    solve_time(1) = solve_time(1)+mpi_wtime()-tick
    call assert_ok(error)
    tick = mpi_wtime()
    call cmm_solve(full_cache,rhs,constants,full_solution,error)
    solve_time(2) = solve_time(2)+mpi_wtime()-tick
    call assert_ok(error)
    comparison = maxval(abs(solution-full_solution))
    if (trial == 1) first = solution
    if (trial == 2 .and. fixture /= 'smooth') then
      do j = 1,geom%nj
        do i = 1,geom%ni
          if (geom%global%kind(geom%i0+i,geom%j0+j) == solid_node) cycle
          comparison = max(comparison,maxval(abs(solution(i,j,:)-first(i,j,:))))
        end do
      end do
    end if
    if (trial == 3) comparison = max(comparison,maxval(abs(solution)))
    call mpi_allreduce(comparison,maximum,1,MPI_DOUBLE_PRECISION,MPI_MAX,geom%comm,ierr)
    if (maximum > 1.0D-10) call pabort('CMM Green cache, land extension, or zero-RHS check failed',__FILE__,__LINE__)
    send = 0.0D0
    send(geom%i0+1:geom%i0+geom%ni,geom%j0+1:geom%j0+geom%nj,:) = solution
    call mpi_reduce(send,receive,size(send),MPI_DOUBLE_PRECISION,MPI_SUM,0,geom%comm,ierr)
    if (geom%rank == 0) write(unit) receive
  end do
  stats = cmm_get_statistics(cache)
  full_stats = cmm_get_statistics(full_cache)
  if (stats%constraints == 0) then
    if (stats%factors /= 0 .or. stats%background_batches /= 3) call pabort('zero-constraint rectangle failed',__FILE__,__LINE__)
  else
    if (stats%factors /= 2 .or. stats%background_batches /= stats%constraints+6 .or. &
        full_stats%background_batches /= stats%constraints+3) call pabort('factor/rectangle reuse count failed',__FILE__,__LINE__)
    if (stats%local_green_bytes /= 0 .or. full_stats%local_green_bytes == 0) &
      call pabort('CMM storage accounting failed',__FILE__,__LINE__)
  end if
  if (full_stats%explicit_array_bytes-stats%explicit_array_bytes /= &
      full_stats%local_green_bytes+full_stats%factor_bytes-stats%factor_bytes) &
    call pabort('CMM explicit storage accounting failed',__FILE__,__LINE__)
  if (geom%rank == 0 .and. stats%constraints > 0 .and. .not. distributed) then
    if (stats%dense_setup_peak_bytes <= stats%factor_bytes) call pabort('dense workspace missing',__FILE__,__LINE__)
  end if
  local_bytes = [stats%explicit_array_bytes,full_stats%explicit_array_bytes,stats%dense_setup_peak_bytes]
  call mpi_reduce(local_bytes,storage_totals,3,MPI_INTEGER8,MPI_SUM,0,geom%comm,ierr)
  call mpi_reduce([setup_time,solve_time/3],max_times,4,MPI_DOUBLE_PRECISION,MPI_MAX,0,geom%comm,ierr)
  if (geom%rank == 0) then
    print '(a,3(1x,i0))', 'aggregate explicit two-solve, Green, dense setup bytes:', storage_totals
    print '(a,4(1x,es14.6))', 'two-solve and Green max setup and mean solve seconds:', max_times
    do mode = 1,mode_count
      call cmm_copy_matrix(full_cache,mode,matrix)
      write(unit) matrix
    end do
    close(unit)
    print '(a,4(1x,i0),3(1x,es14.6))', 'constraints, factors, root factor bytes, root cached-Green bytes, residual, asymmetry, rcond:', &
      stats%constraints,stats%factors,stats%factor_bytes,full_stats%local_green_bytes, &
      stats%max_rectangle_backward_error,stats%relative_asymmetry,stats%min_rcond
  end if
  call release_cmm(full_cache)
  call release_cmm(cache)
  if (reordered) call mpi_comm_free(comm,ierr)
  call deallocate(decomp)
  call parallel_finalise()

contains
  subroutine assert_ok(error)
    character(*), intent(in) :: error
    if (len_trim(error) /= 0) call pabort(error,__FILE__,__LINE__)
  end subroutine assert_ok
end program test_cmm_dirichlet
