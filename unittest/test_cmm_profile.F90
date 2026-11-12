program test_cmm_profile
  use qg_cmm
  use qg_cmm_dense, only: has_scalapack
  use qg_irregular_parallel, only: distributed_geometry, load_distributed_geometry
  use qg_parallel, only: decomposition, allocate, deallocate, parallel_init, parallel_finalise, pabort
  use qg_parameters, only: real_8
  implicit none
  include 'mpif.h'

  integer, parameter :: modes=3,trials=8,wide=selected_int_kind(18)
  type(decomposition), target :: decomp
  type(distributed_geometry), target :: geometry
  type(cmm_cache) :: inverse
  type(cmm_statistics) :: stats
  integer, allocatable :: iis(:),jjs(:)
  real(real_8), allocatable :: rhs(:,:,:),solution(:,:,:),constant(:,:)
  real(real_8) :: time0,times(3),elapsed(3),norm,local_norm
  integer(wide) :: local(3),aggregate(3),maximum(3)
  integer :: n,i,j,k,grid(2),ranks,ierr
  logical :: green,distributed
  character(1024) :: path,error
  character(32) :: arg,method

  if (command_argument_count() /= 3) then
    print '(a)', 'CMM profile requires geometry, dimension and two-solve/green/distributed'
    stop 2
  end if
  call get_command_argument(1,path)
  call get_command_argument(2,arg)
  read(arg,*) n
  call get_command_argument(3,method)
  green = method == 'green'
  distributed = method == 'distributed'
  if (.not.(green .or. distributed .or. method=='two-solve')) error stop 'unknown CMM profile method'

  call parallel_init()
  call mpi_comm_size(MPI_COMM_WORLD,ranks,ierr)
  grid = 0
  call mpi_dims_create(ranks,2,grid,ierr)
  iis = [(i*n/grid(1)-(i-1)*n/grid(1),i=1,grid(1))]
  jjs = [(j*n/grid(2)-(j-1)*n/grid(2),j=1,grid(2))]
  call allocate(decomp,iis=iis,jjs=jjs,layers=modes)
  time0 = mpi_wtime()
  call load_distributed_geometry(path,n,n,[decomp%ii_0,decomp%jj_0,decomp%ii,decomp%jj], &
                                decomp%cart_comm,1.0D0,geometry,error)
  call check(error)
  times(1) = mpi_wtime()-time0
  time0 = mpi_wtime()
  call setup_cmm(inverse,decomp,geometry,[0.0D0,0.3D0,0.3D0],error, &
                 retain_green=green,distributed_factor=distributed)
  if (distributed .and. .not. has_scalapack) then
    if (geometry%rank == 0) print '(a)', 'SKIP distributed CMM requires PETSc ScaLAPACK'
    call deallocate(decomp)
    call parallel_finalise()
    stop
  end if
  call check(error)
  times(2) = mpi_wtime()-time0
  allocate(rhs(geometry%ni,geometry%nj,modes),solution(geometry%ni,geometry%nj,modes), &
           constant(0:ubound(geometry%global%loops,1),modes))
  do j = 1,geometry%nj
    do i = 1,geometry%ni
      rhs(i,j,:) = sin(0.01D0*(geometry%i0+i))*cos(0.02D0*(geometry%j0+j))*[1.0D0,0.5D0,-2.0D0]
    end do
  end do
  constant = 0.0D0
  time0 = mpi_wtime()
  do k = 1,trials
    call cmm_solve(inverse,rhs,constant,solution,error)
    call check(error)
  end do
  times(3) = (mpi_wtime()-time0)/trials
  stats = cmm_get_statistics(inverse)
  if (stats%constraints>0) then
    if (stats%factors/=2 .or. stats%background_batches/=stats%constraints+merge(trials,2*trials,green)) &
      call pabort('profile did not reuse capacitance factors',__FILE__,__LINE__)
  end if
  local = [stats%explicit_array_bytes,stats%factor_bytes,stats%local_green_bytes]
  call mpi_reduce(local,aggregate,3,MPI_INTEGER8,MPI_SUM,0,geometry%comm,ierr)
  call mpi_reduce(local,maximum,3,MPI_INTEGER8,MPI_MAX,0,geometry%comm,ierr)
  call mpi_reduce(times,elapsed,3,MPI_DOUBLE_PRECISION,MPI_MAX,0,geometry%comm,ierr)
  local_norm = sum(solution**2)
  call mpi_reduce(local_norm,norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,geometry%comm,ierr)
  if (geometry%rank==0) then
    print '(a,3(1x,i0))','aggregate explicit solver/factor/Green bytes:',aggregate
    print '(a,3(1x,i0))','max-rank explicit solver/factor/Green bytes:',maximum
    print '(a,4(1x,es15.7))','max geometry/setup/mean solve seconds and solution norm²:',elapsed,norm
    print '(a)','CMM profile PASS'
  end if
  call release_cmm(inverse)
  call deallocate(decomp)
  call parallel_finalise()

contains

  subroutine check(message)
    character(*), intent(in) :: message

    if (len_trim(message)>0) call pabort(message,__FILE__,__LINE__)
  end subroutine check

end program test_cmm_profile
