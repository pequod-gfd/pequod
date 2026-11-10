#include "config.h"

! Base parallisation routines for the QG model.
module qg_parallel

  use qg_parameters, only : real_8
  use test_tools
#ifdef PETSC
  use petscsys, only : PetscInitialize, PetscFinalize, PETSC_NULL_CHARACTER
#endif

  implicit none

  include "mpif.h"
  
  private
  
  public :: nan, rseed, pabort, parallel_init, parallel_finalise, nprocs, &
    & procno, alland, allor, allmax, allmin, allsum, parallel_filename_len, &
    & parallel_filename, int2str, int2str_len
  public :: allocate, deallocate, zero, set, add, subtract, init_rand, &
    & scatter, gather, print, update, ptranspose, pinverse_transpose, flush
  public :: decomposition
  public :: transpose_cache, request_cache, rcache, increment_tag, next_request
  public :: parallel_data, test_pdata
  public :: parallel_flux_data
    
  interface allmax
    module procedure allmax_integer, allmax_real
  end interface allmax
    
  interface allmin
    module procedure allmin_integer, allmin_real
  end interface allmin
    
  interface allsum
    module procedure allsum_integer, allsum_real, allsum_real_rank1, &
      & allsum_real_rank2
  end interface allsum    
  
  interface parallel_filename_len
    module procedure parallel_filename_decomp_len, parallel_filename_procno_len
  end interface parallel_filename_len
  
  interface parallel_filename
    module procedure parallel_filename_decomp, parallel_filename_procno
  end interface parallel_filename
  
  ! Non-blocking send / receive requests cache.
  type request_cache
    integer :: comm                              ! MPI communicator
    integer, dimension(:), pointer :: requests   ! MPI requests
    integer, dimension(:), pointer :: statuses   ! MPI statuses from the last waitall
    integer :: max_nrequests                     ! The maximum number of requests
    integer :: nrequests                         ! Active requests
    integer :: tag_base                          ! Current tag base
  end type request_cache    
  ! Global request_cache. Used by default in relevant routines if no custom
  ! request_cache is supplied.  
  type(request_cache), target, save :: rcache
    
  ! Cache of data used to transpose a parallel data object.
  integer, parameter, public :: TRANSPOSE_II = 1, TRANSPOSE_JJ = 2
  type transpose_cache
    ! Target decomposition type
    integer :: transpose_type
    ! Target decomposition
    type(decomposition), pointer :: decomp
    ! Whether the target decomposition was allocated by the transpose_cache
    ! allocation
    logical :: decomp_allocated
    ! The number of processes to which a given process sends / from which a
    ! given process receives
    integer :: m
    ! The lowest rank process to which a given process sends
    integer :: p_0
    ! Send and receive MPI types for local data matrix (owned), including
    ! boundary nodes
    integer, dimension(:), pointer :: lb_send_types, lb_receive_types
    ! Send and receive MPI types for full data matrix (owned + halo), excluding
    ! halo nodes and including boundary nodes
    integer, dimension(:), pointer :: fb_send_types, fb_receive_types
    ! Start indices in local data matrix (owned) for sends / receives, including
    ! boundary nodes
    integer, dimension(:), pointer :: send_lb_i, send_lb_j
    integer, dimension(:), pointer :: receive_lb_i, receive_lb_j
    ! Start indices in full data matrix (owned + halo) for sends / receives,
    ! excluding halo nodes and including boundary nodes
    integer, dimension(:), pointer :: send_fb_i, send_fb_j
    integer, dimension(:), pointer :: receive_fb_i, receive_fb_j
    
    ! Communicator
    integer :: comm
    ! MPI requests
    integer, dimension(:), pointer :: requests
    ! MPI statuses
    integer, dimension(:), pointer :: statuses
  end type transpose_cache
  
  ! Measures the time spent waiting for transpose communication to complete.
  real(kind = real_8), save, public :: transpose_wait_time = 0.0D0
  real(kind = real_8), save, public :: halo_wait_time = 0.0D0
  ! Used to measure the start time of communication.
  real(kind = real_8), save :: wtime0
  
  ! Decomposition data.
  type decomposition
    ! MPI data
    integer :: n           ! The number of processes
    integer :: p           ! The process rank (indexed from zero)
    
    logical :: periodic_j  ! Periodic in the j-direction?
    
    ! Domain decomposition data
    integer :: cart_comm   ! Topological communicator
    integer :: p_ii, p_jj  ! The number of processes in each dimension
    integer :: p_i, p_j    ! The process logical location in the domain
                           ! decomposition (indexed from one)
    integer :: p_a, p_b, p_c, p_d      ! The process rank of the logical
    integer :: p_ab, p_bc, p_cd, p_da  ! neighbours (indexed from zero)
    integer, dimension(:), pointer :: iis, jjs  ! The local domain sizes for all
                                                ! processes in each row / column
    integer, dimension(:), pointer :: ii2s, jj2s  ! The local domain sizes for
                                                  ! all processes in each row /
                                                  ! column excluding boundary
                                                  ! nodes
    integer :: layers        ! Number of layers
    integer :: ii, jj        ! The local domain size, including boundary nodes
    integer :: ii2, jj2      ! The local domain size, excluding boundary nodes
    integer :: ii_0, jj_0    ! Local index bases
    integer :: ii2_0, jj2_0
    integer :: fa_ii, fa_jj  ! Domain size, including boundary nodes and
                             ! non-owned (halo) nodes
                               
    ! Global domain size
    integer :: g_ii, g_jj          ! Global number of nodes, including boundary
                                   ! nodes
    integer :: g_ii2, g_jj2        ! Global number of nodes, excluding boundary
                                   ! nodes
    integer :: g_ii_int, g_jj_int  ! Global number of node intervals
    
    ! Index ranges for parallel data
    integer :: l_i_min, l_i_max  ! Non-boundary node index ranges in local
    integer :: l_j_min, l_j_max  ! matrix (owned)
    integer :: f_i_min, f_i_max  ! Owned, non-boundary node index ranges in full
    integer :: f_j_min, f_j_max  ! data matrix (owned + halo)
    integer :: fa_i_min, fa_i_max  ! Full node index ranges in full data matrix
    integer :: fa_j_min, fa_j_max  ! (owned + halo), including boundary nodes
                                   ! and non-owned (halo) nodes
    integer :: fb_i_min, fb_i_max  ! Node index ranges in full data matrix
    integer :: fb_j_min, fb_j_max  ! (owned + halo), excluding halo nodes and
                                   ! including boundary nodes
    ! Striding MPI types
    integer :: ts_jj_iip2   ! Vector of jj elements striding by ii + 2.
                            ! Used to update parallel_data objects.
    integer :: ts_jj2_iip3  ! Vector of jj2 elements striding by ii + 3
                            ! Used to update fluxes stored on cell faces.
                            
    ! Decomposition of a parallel data object transposed to a slab decomposition
    ! contiguous in i. This is null if p_ii == 1.
    type(decomposition), pointer :: trans_i_decomp
    ! Cache of data used to perform the transpose. This is unallocated if
    ! p_ii == 1.
    type(transpose_cache) :: trans_i_c

    ! Decomposition of a parallel data object transposed to a slab decomposition
    ! contiguous in j. This is null if p_jj == 1.
    type(decomposition), pointer :: trans_j_decomp
    ! Cache of data used to perform the transpose. This is unallocated if
    ! p_jj == 1.
    type(transpose_cache) :: trans_j_c    
  end type decomposition 
  
  ! Parallel data.
  type parallel_data
    ! Domain decomposition data
    type(decomposition), pointer :: decomp
    ! Send / receive request cache
    type(request_cache), pointer :: requests
    
    ! Communication data
    real(kind = real_8), dimension(:, :, :), pointer :: data   ! Full data matrix (owned + halo)
    real(kind = real_8), dimension(:, :, :), pointer :: sdata  ! Local data matrix (owned)
  end type parallel_data
  
  ! Parallel flux data.
  integer, parameter, public :: FLUX_II = 1, FLUX_JJ = 2
  type parallel_flux_data
    ! Domain decomposition data
    type(decomposition), pointer :: decomp
    ! Send / receive cache
    type(request_cache), pointer :: requests
    ! Stagger direction
    integer :: flux_dir
    
    ! Communication data
    real(kind = real_8), dimension(:, :, :), pointer :: data  ! Full data matrix (owned + halo)
  end type parallel_flux_data
  
  interface allocate
    module procedure allocate_transpose_cache, allocate_requests, &
      & allocate_decomp, allocate_pdata_rank0, allocate_pdata_rank1, &
      & allocate_fpdata
  end interface allocate
  
  interface deallocate
    module procedure deallocate_transpose_cache, deallocate_requests, &
      & deallocate_decomp, deallocate_pdata, deallocate_pdata_rank1, &
      & deallocate_fpdata
  end interface deallocate
  
  interface zero
    module procedure zero_pdata, zero_pdata_rank1, zero_fpdata
  end interface zero
  
  interface set
    module procedure set_data_data, set_data_pdata, set_pdata_real, &
      & set_pdata_pdata, set_pdata_pdata_rank1, set_fpdata_real
  end interface set
  
  interface add
    module procedure add_pdata_pdata
  end interface add
  
  interface subtract
    module procedure subtract_pdata_pdata
  end interface subtract
  
  interface init_rand
    module procedure rand_pdata_rank0, rand_pdata_rank1
  end interface init_rand
  
  interface scatter
    module procedure scatter_pdata_data_rank0_sender, &
      & scatter_pdata_data_rank1_sender, scatter_pdata_data_rank0_receiver, &
      & scatter_pdata_data_rank1_receiver
  end interface scatter

  interface gather
    module procedure gather_pdata_data_rank0_sender, &
      & gather_pdata_data_rank1_sender, gather_pdata_data_rank0_receiver, &
      & gather_pdata_data_rank1_receiver
  end interface gather
  
  interface print
    module procedure print_decomp
  end interface print
  
  interface update
    module procedure update_pdata, update_pdata_rank1, update_fpdata
  end interface update
  
  interface ptranspose
    module procedure transpose_pdata, transpose_data_layer, &
      & transpose_data_layers
  end interface ptranspose
  interface transpose_data
    module procedure transpose_data_layer, transpose_data_layers
  end interface transpose_data
  interface transpose_sdata_i
    module procedure transpose_sdata_i_layer, transpose_sdata_i_layers
  end interface transpose_sdata_i
  interface transpose_data_i
    module procedure transpose_data_i_layer, transpose_data_i_layers
  end interface transpose_data_i
  interface transpose_sdata_j
    module procedure transpose_sdata_j_layer, transpose_sdata_j_layers
  end interface transpose_sdata_j
  interface transpose_data_j
    module procedure transpose_data_j_layer, transpose_data_j_layers
  end interface transpose_data_j
  
  interface pinverse_transpose
    module procedure inverse_transpose_pdata, inverse_transpose_data_layer, &
      & inverse_transpose_data_layers
  end interface pinverse_transpose
  interface inverse_transpose_data
    module procedure inverse_transpose_data_layer, inverse_transpose_data_layers
  end interface inverse_transpose_data
  interface inverse_transpose_sdata_i
    module procedure inverse_transpose_sdata_i_layer, &
      & inverse_transpose_sdata_i_layers
  end interface inverse_transpose_sdata_i
  interface inverse_transpose_data_i
    module procedure inverse_transpose_data_i_layer, &
      & inverse_transpose_data_i_layers
  end interface inverse_transpose_data_i
  interface inverse_transpose_sdata_j
    module procedure inverse_transpose_sdata_j_layer, &
      & inverse_transpose_sdata_j_layers
  end interface inverse_transpose_sdata_j
  interface inverse_transpose_data_j
    module procedure inverse_transpose_data_j_layer, &
      & inverse_transpose_data_j_layers
  end interface inverse_transpose_data_j
  
  interface waitall
    module procedure waitall_requests
  end interface waitall
  
  interface flush
    module procedure flush_transpose_cache, flush_ptranspose_decomp, &
      & flush_ptranspose_pdata, flush_requests, flush_pdata_rank0, &
      & flush_pdata_rank1, flush_pfdata
  end interface flush
  
  interface next_request
    module procedure next_request_requests
  end interface next_request
  
  interface increment_tag
    module procedure increment_tag_requests
  end interface increment_tag
  
contains

  ! Return a NaN.
  function nan()
    real(kind = real_8) :: nan
    
    real(kind = real_8), save :: r1 = 0.0D0, r2 = 0.0D0
    
    nan = r1 / r2
    
  end function nan
  
  ! Return whether the supplied optional is both present and true.
  pure function present_and_true(l) result(p)
    logical, optional, intent(in) :: l
    
    logical :: p
    
    if(present(l)) then
      p = l
    else
      p = .false.
    end if
  
  end function present_and_true
  
  ! Seed the random number generator using an integer.
  subroutine rseed(i)
    integer, intent(in) :: i
    
    integer :: rsize
    integer, dimension(:), allocatable :: put
    
    call random_seed(size = rsize)
    allocate(put(rsize))
    put = i
    call random_seed(put = put)
    deallocate(put)
    
  end subroutine rseed

  ! The length of the supplied integer when written to a string.
  pure function int2str_len(i) result(length)    
    integer, intent(in) :: i
    
    integer :: length
    
    if(i == 0) then
      length = 1
    else
      length = int(dlog10(dabs(real(i, kind = real_8)))) + 1
      if(i < 0) length = length + 1
    end if
  
  end function int2str_len

  ! Write the supplied integer to a string.
  function int2str(i) result(s)  
    integer, intent(in) :: i
    
    character(len = int2str_len(i)) :: s
        
    write(s, "(i0)") i
  
  end function int2str 

  ! Abort with an error message, filename, and file line number.
  subroutine pabort(error, filename, line)  
    character(len = *), intent(in) :: error
    character(len = *), intent(in) :: filename
    integer, intent(in) :: line
    
    integer :: ierr, init
   
    call mpi_initialized(init, ierr)
    if(ierr /= MPI_SUCCESS) then
      write(0, "(a)") "Failed to determine MPI initialisation status"
    else if(init == 1) then
      write(0, "(a,i0,a,i0)") "Error on process ", procno(), ", line " // trim(filename) // ":", line
      write(0, "(a)") trim(error)
      call mpi_abort(MPI_COMM_WORLD, MPI_ERR_OTHER, ierr)
    else
      write(0, "(a,i0)") "Error on line " // trim(filename) // ":", line
      write(0, "(a)") trim(error)
    end if
    stop 1
    
  end subroutine pabort
  
  ! Parallel initialisation.
  subroutine parallel_init()

    integer :: ierr

#ifdef PETSC
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
#else
    call mpi_init(ierr)
#endif

    call allocate(rcache, comm = MPI_COMM_WORLD)

  end subroutine parallel_init
  
  ! Parallel finalisation.
  subroutine parallel_finalise()

    integer :: ierr

    call deallocate(rcache)
  
#ifdef PETSC
    call PetscFinalize(ierr)
#else
    call mpi_finalize(ierr)
#endif
  
  end subroutine parallel_finalise
  
  ! The number of processes.
  function nprocs()  
    integer :: nprocs
    
    integer :: ierr
    
    call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr);  assert(ierr == MPI_SUCCESS)
    
  end function nprocs
  
  ! The process rank.
  function rank()  
    integer :: rank
    
    integer :: ierr
    
    call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr);  assert(ierr == MPI_SUCCESS)
  
  end function rank
  
  ! The process number (rank + 1).
  function procno()  
    integer :: procno
    
    integer :: ierr
    
    call mpi_comm_rank(MPI_COMM_WORLD, procno, ierr);  assert(ierr == MPI_SUCCESS)
    procno = procno + 1
  
  end function procno
  
  ! Collective logical and operation.
  subroutine alland(b)
    logical, intent(inout) :: b
    
    integer :: lb
    
    if(b) then
      lb = 1
    else
      lb = 0
    end if
    
    call allmin(lb)
    
    b = lb == 1
    
  end subroutine alland
  
  ! Collective logical or operation.
  subroutine allor(b)
    logical, intent(inout) :: b
    
    integer :: lb
    
    if(b) then
      lb = 1
    else
      lb = 0
    end if
    
    call allmax(lb)
    
    b = lb == 1
    
  end subroutine allor
  
  ! Collective max operation.
  subroutine allmax_integer(i)
    integer, intent(inout) :: i
    
    integer :: ierr, li
    
    li = i
    call mpi_allreduce(li, i, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
  
  end subroutine allmax_integer
    
  ! Collective max operation.
  subroutine allmax_real(r)
    real(kind = real_8), intent(inout) :: r
    
    integer :: ierr
    real(kind = real_8) :: lr
    
    lr = r
    call mpi_allreduce(lr, r, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
  
  end subroutine allmax_real
  
  ! Collective min operation.
  subroutine allmin_integer(i)
    integer, intent(inout) :: i
    
    integer :: ierr, li
    
    li = i
    call mpi_allreduce(li, i, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
  
  end subroutine allmin_integer
  
  ! Collective min operation.
  subroutine allmin_real(r)
    real(kind = real_8), intent(inout) :: r
    
    integer :: ierr
    real(kind = real_8) :: lr
    
    lr = r
    call mpi_allreduce(lr, r, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
  
  end subroutine allmin_real
  
  ! Collective sum operation.
  subroutine allsum_integer(i)
    integer, intent(inout) :: i
    
    integer :: ierr, li
    
    li = i
    call mpi_allreduce(li, i, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
  
  end subroutine allsum_integer
  
  ! Collective sum operation.
  subroutine allsum_real(r)
    real(kind = real_8), intent(inout) :: r
    
    integer :: ierr
    real(kind = real_8) :: lr
    
    lr = r
    call mpi_allreduce(lr, r, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
  
  end subroutine allsum_real
  
  ! Collective sum operation.
  subroutine allsum_real_rank1(r)
    real(kind = real_8), dimension(:), intent(inout) :: r
    
    integer :: ierr
    real(kind = real_8), dimension(size(r)) :: lr
    
    lr = r
    
    call mpi_allreduce(lr, r, size(r), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
  
  end subroutine allsum_real_rank1
  
  ! Collective sum operation.
  subroutine allsum_real_rank2(r)
    real(kind = real_8), dimension(:, :), intent(inout) :: r
    
    integer :: ierr
    real(kind = real_8), dimension(size(r, 1), size(r, 2)) :: lr
    
    lr = r
    
    call mpi_allreduce(lr, r, size(r, 1) * size(r, 2), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
  
  end subroutine allsum_real_rank2

  ! The length of a parallel filename for this process.
  pure function parallel_filename_procno_len(base) result(length)
    character(len = *), intent(in) :: base
    
    integer :: length
    
    length = len_trim(base) + int2str_len(huge(0))
    
  end function parallel_filename_procno_len
  
  ! A parallel filename for this process. Note that the returned string must be
  ! trimmed.
  function parallel_filename_procno(base) result(filename)
    character(len = *), intent(in) :: base
    
    character(len = parallel_filename_len(base)) :: filename

    if(nprocs() > 1) then
      write(filename, "(a,i0)") trim(base) // ".", procno() - 1
    else
      filename = trim(base)
    end if
  
  end function parallel_filename_procno

  ! The length of a parallel filename for this process.
  pure function parallel_filename_decomp_len(base, decomp) result(length)  
    character(len = *), intent(in) :: base
    type(decomposition), intent(in) :: decomp
    
    integer :: length
    
    if(decomp%n > 1) then
      length = len_trim(base) + 1 + int2str_len(decomp%p)
    else
      length = len_trim(base)
    end if
    
  end function parallel_filename_decomp_len
  
  ! A parallel filename for this process.
  function parallel_filename_decomp(base, decomp) result(filename)  
    character(len = *), intent(in) :: base
    type(decomposition), intent(in) :: decomp
    
    character(len = parallel_filename_len(base, decomp)) :: filename
    
    if(decomp%n > 1) then
      write(filename, "(a,i0)") trim(base) // ".", decomp%p
    else
      filename = trim(base)
    end if
  
  end function parallel_filename_decomp
  
  ! Allocate a requests cache.
  subroutine allocate_requests(requests, comm)
    type(request_cache), intent(out) :: requests
    integer, optional, intent(in) :: comm
    
    integer :: ierr
    integer, parameter :: min_requests = 1024
    
    if(present(comm)) then
      call mpi_comm_dup(comm, requests%comm, ierr);  assert(ierr == MPI_SUCCESS)
    else
      call mpi_comm_dup(MPI_COMM_WORLD, requests%comm, ierr);  assert(ierr == MPI_SUCCESS)
    end if
    
    allocate(requests%requests(min_requests), requests%statuses(min_requests * MPI_STATUS_SIZE))
    requests%requests = MPI_REQUEST_NULL
    requests%statuses = 0
    requests%max_nrequests = min_requests
    
    requests%nrequests = 0
    requests%tag_base = 1
  
  end subroutine allocate_requests
  
  subroutine deallocate_requests(requests)
    type(request_cache), intent(inout) :: requests
    
    integer :: ierr
    
    call mpi_comm_free(requests%comm, ierr);  assert(ierr == MPI_SUCCESS)
    deallocate(requests%requests, requests%statuses)
    
  end subroutine deallocate_requests
  
  ! Increase maximum number of available requests.
  subroutine increase_requests(requests)
    type(request_cache), intent(inout) :: requests
    
    integer :: new_max_nrequests, old_max_nrequests
    integer, dimension(:), allocatable :: old_requests, old_statuses

    old_max_nrequests = requests%max_nrequests
    new_max_nrequests = old_max_nrequests * 2
    assert(new_max_nrequests > old_max_nrequests)
    
    allocate(old_requests(old_max_nrequests), old_statuses(old_max_nrequests * MPI_STATUS_SIZE))
    old_requests = requests%requests
    old_statuses = requests%statuses

    deallocate(requests%requests, requests%statuses)
    allocate(requests%requests(new_max_nrequests), requests%statuses(new_max_nrequests * MPI_STATUS_SIZE))
    requests%requests(:old_max_nrequests) = old_requests
    requests%requests(old_max_nrequests + 1:) = MPI_REQUEST_NULL
    requests%statuses(:old_max_nrequests * MPI_STATUS_SIZE) = old_statuses
    requests%statuses(old_max_nrequests * MPI_STATUS_SIZE + 1:) = 0
    requests%max_nrequests = new_max_nrequests
    
    deallocate(old_requests, old_statuses)
    
  end subroutine increase_requests
  
  ! Wait for all non-blocking send / receives to complete.
  subroutine waitall_requests(requests)  
    type(request_cache), intent(inout) :: requests
    
    integer :: ierr
    
    ! Prevents Electric Fence reporting zero sized allocation in the mpi_waitall
    if(requests%nrequests == 0) return
    
!    write(6, "(a,i0,a)") "Waiting for ", requests%nrequests, " send / receives to complete"
    call mpi_waitall(requests%nrequests, requests%requests, requests%statuses, ierr);  assert(ierr == MPI_SUCCESS)
!    write(6, "(a)") "Communication complete"
      
  end subroutine waitall_requests
  
  ! Flush a non-blocking send / receive cache (including waiting for all send /
  ! receives to complete).
  subroutine flush_requests(requests)  
    type(request_cache), target, optional, intent(inout) :: requests
    
    type(request_cache), pointer :: lrequests
    
    if(present(requests)) then
      lrequests => requests
    else
      lrequests => rcache
    end if
    
    wtime0 = mpi_wtime()
    call waitall(lrequests)
    halo_wait_time = halo_wait_time + mpi_wtime() - wtime0
    lrequests%nrequests = 0
    lrequests%tag_base = 1
    
  end subroutine flush_requests
  
  ! Return the next request from the send / receive cache.
  function next_request_requests(requests) result(request)    
    type(request_cache), target, intent(inout) :: requests
    
    integer, pointer :: request
    
    requests%nrequests = requests%nrequests + 1
    if(requests%nrequests > requests%max_nrequests) call increase_requests(requests)
    request => requests%requests(requests%nrequests)
   
  end function next_request_requests
   
  ! Increment the MPI tag base for the send / receive cache.
  subroutine increment_tag_requests(requests, increment)
    type(request_cache), intent(inout) :: requests
    integer, optional, intent(in) :: increment
    
    if(present(increment)) then
      requests%tag_base = requests%tag_base + increment
    else
      requests%tag_base = requests%tag_base + 1
    end if
    
  end subroutine increment_tag_requests
  
  ! Allocate a decomain decomposition.
  recursive subroutine allocate_decomp(decomp, iis, jjs, layers, periodic_j &
    & , decomp_i, decomp_j)  
    type(decomposition), target, intent(out) :: decomp
    integer, dimension(:), intent(in) :: iis
    integer, dimension(:), intent(in) :: jjs
    integer, intent(in) :: layers
    logical, optional, intent(in) :: periodic_j
    type(decomposition), optional, intent(in) :: decomp_i
    type(decomposition), optional, intent(in) :: decomp_j
    
    integer :: ierr
    integer, dimension(2) :: cart_coord, cart_periods
#ifdef DDEBUG
    integer :: init
    
    call mpi_initialized(init, ierr);  assert(ierr == MPI_SUCCESS)
    assert(init /= 0)
#endif
    
    ! Periodic?
    if(present(periodic_j)) then
      decomp%periodic_j = periodic_j
    else
      decomp%periodic_j = .false.
    end if
    
    ! Number of layers
    decomp%layers = layers
    
    ! Domain decomposition
    assert(size(iis) > 0 .and. size(jjs) > 0)
    if(size(iis) * size(jjs) /= nprocs()) then
      write(0, "(a,i0,a,i0,a,i0)") "p_ii, p_jj, n = ", size(iis), ", ", size(jjs), ", ", nprocs()
      labort("Invalid decomposition")
    end if
    decomp%p_ii = size(iis)
    decomp%p_jj = size(jjs)
    
    ! Global domain size
    decomp%g_ii = sum(iis)
    decomp%g_jj = sum(jjs)
    decomp%g_ii2 = decomp%g_ii - 2
    decomp%g_ii_int = decomp%g_ii - 1
    if(decomp%periodic_j) then
      decomp%g_jj2 = decomp%g_jj
      decomp%g_jj_int = decomp%g_jj
    else
      decomp%g_jj2 = decomp%g_jj - 2
      decomp%g_jj_int = decomp%g_jj - 1
    end if
    
    ! Cartesian topology
    if(decomp%periodic_j) then
      cart_periods = (/0, 1/)
    else
      cart_periods = (/0, 0/)
    end if
    call mpi_cart_create(MPI_COMM_WORLD, 2, (/decomp%p_ii, decomp%p_jj/), cart_periods, 1, decomp%cart_comm, ierr);  assert(ierr == MPI_SUCCESS)
    call mpi_comm_rank(decomp%cart_comm, decomp%p, ierr);  assert(ierr == MPI_SUCCESS)
    call mpi_comm_size(decomp%cart_comm, decomp%n, ierr);  assert(ierr == MPI_SUCCESS)
    call mpi_cart_coords(decomp%cart_comm, decomp%p, 2, cart_coord, ierr)
    decomp%p_i = cart_coord(1) + 1
    decomp%p_j = cart_coord(2) + 1
    
    ! Process neighbours
    decomp%p_a = -1;  decomp%p_b = -1;  decomp%p_c = -1;  decomp%p_d = -1
    decomp%p_ab = -1;  decomp%p_bc = -1;  decomp%p_cd = -1;  decomp%p_da = -1
    if(decomp%p_j < decomp%p_jj .or. decomp%periodic_j) then
      call mpi_cart_rank(decomp%cart_comm, (/cart_coord(1), cart_coord(2) + 1/), decomp%p_a, ierr);  assert(ierr == MPI_SUCCESS)
      if(decomp%p_i > 1) call mpi_cart_rank(decomp%cart_comm, (/cart_coord(1) - 1, cart_coord(2) + 1/), decomp%p_ab, ierr);  assert(ierr == MPI_SUCCESS)
      if(decomp%p_i < decomp%p_ii) call mpi_cart_rank(decomp%cart_comm, (/cart_coord(1) + 1, cart_coord(2) + 1/), decomp%p_da, ierr);  assert(ierr == MPI_SUCCESS)
    end if
    if(decomp%p_j > 1 .or. decomp%periodic_j) then
      call mpi_cart_rank(decomp%cart_comm, (/cart_coord(1), cart_coord(2) - 1/), decomp%p_c, ierr);  assert(ierr == MPI_SUCCESS)
      if(decomp%p_i > 1) call mpi_cart_rank(decomp%cart_comm, (/cart_coord(1) - 1, cart_coord(2) - 1/), decomp%p_bc, ierr);  assert(ierr == MPI_SUCCESS)
      if(decomp%p_i < decomp%p_ii) call mpi_cart_rank(decomp%cart_comm, (/cart_coord(1) + 1, cart_coord(2) - 1/), decomp%p_cd, ierr);  assert(ierr == MPI_SUCCESS)
    end if
    if(decomp%p_i > 1) call mpi_cart_rank(decomp%cart_comm, (/cart_coord(1) - 1, cart_coord(2)/), decomp%p_b, ierr);  assert(ierr == MPI_SUCCESS)
    if(decomp%p_i < decomp%p_ii) call mpi_cart_rank(decomp%cart_comm, (/cart_coord(1) + 1, cart_coord(2)/), decomp%p_d, ierr);  assert(ierr == MPI_SUCCESS)
    
    ! Local domain size, including and excluding boundary nodes
    allocate(decomp%iis(decomp%p_ii));  decomp%iis = iis
    allocate(decomp%jjs(decomp%p_jj));  decomp%jjs = jjs
    
    decomp%ii = iis(decomp%p_i);  decomp%ii2 = decomp%ii
    decomp%jj = jjs(decomp%p_j);  decomp%jj2 = decomp%jj 
    if(decomp%p_b == -1) decomp%ii2 = decomp%ii2 - 1
    if(decomp%p_d == -1) decomp%ii2 = decomp%ii2 - 1
    if(decomp%p_c == -1) decomp%jj2 = decomp%jj2 - 1
    if(decomp%p_a == -1) decomp%jj2 = decomp%jj2 - 1
    
    allocate(decomp%ii2s(decomp%p_ii))
    allocate(decomp%jj2s(decomp%p_jj))
    decomp%ii2s = decomp%iis
    decomp%ii2s(1) = decomp%ii2s(1) - 1
    decomp%ii2s(decomp%p_ii) = decomp%ii2s(decomp%p_ii) - 1
    decomp%jj2s = decomp%jjs
    if(.not. decomp%periodic_j) then
      decomp%jj2s(1) = decomp%jj2s(1) - 1
      decomp%jj2s(decomp%p_jj) = decomp%jj2s(decomp%p_jj) - 1
    end if
    
    ! Index bases, including and excluding boundary nodes
    decomp%ii_0 = sum(iis(:decomp%p_i - 1))
    decomp%jj_0 = sum(jjs(:decomp%p_j - 1))
    decomp%ii2_0 = max(decomp%ii_0 - 1, 0)
    if(decomp%periodic_j) then
      decomp%jj2_0 = decomp%jj_0
    else
      decomp%jj2_0 = max(decomp%jj_0 - 1, 0)
    end if
    
    ! Domain size, including boundary and halo nodes
    decomp%fa_ii = decomp%ii2 + 2
    decomp%fa_jj = decomp%jj2 + 2
    
    ! Owned, non-halo and non-boundary index ranges in local and full data
    ! matrices
    if(decomp%p_b >= 0) then
      decomp%f_i_min = 2
    else
      decomp%f_i_min = 3
    end if
    if(decomp%p_d >= 0) then
      decomp%f_i_max = decomp%ii + 1
    else
      decomp%f_i_max = decomp%ii
    end if
    if(decomp%p_c >= 0) then
      decomp%f_j_min = 2
    else
      decomp%f_j_min = 3
    end if
    if(decomp%p_a >= 0) then
      decomp%f_j_max = decomp%jj + 1
    else
      decomp%f_j_max = decomp%jj
    end if    
    decomp%l_i_min = decomp%f_i_min - 1;  decomp%l_i_max = decomp%l_i_min + decomp%ii2 - 1
    decomp%l_j_min = decomp%f_j_min - 1;  decomp%l_j_max = decomp%l_j_min + decomp%jj2 - 1
    
    ! Index ranges in full data matrix, including boundary and halo nodes
    decomp%fa_i_min = decomp%f_i_min - 1;  decomp%fa_i_max = decomp%f_i_max + 1
    decomp%fa_j_min = decomp%f_j_min - 1;  decomp%fa_j_max = decomp%f_j_max + 1
    
    ! Index ranges in full data matrix, including boundary nodes and excluding
    ! halo nodes
    if(decomp%p_b >= 0) then
      decomp%fb_i_min = decomp%f_i_min
    else
      decomp%fb_i_min = decomp%f_i_min - 1
    end if
    if(decomp%p_d >= 0) then
      decomp%fb_i_max = decomp%f_i_max
    else
      decomp%fb_i_max = decomp%f_i_max + 1
    end if
    if(decomp%p_c >= 0) then
      decomp%fb_j_min = decomp%f_j_min
    else
      decomp%fb_j_min = decomp%f_j_min - 1
    end if
    if(decomp%p_a >= 0) then
      decomp%fb_j_max = decomp%f_j_max
    else
      decomp%fb_j_max = decomp%f_j_max + 1
    end if
    
    ! Striding MPI type for halo updates in the j-direction
    call mpi_type_vector(decomp%jj, 1, decomp%ii + 2, MPI_DOUBLE_PRECISION, decomp%ts_jj_iip2, ierr);  assert(ierr == MPI_SUCCESS)
    call mpi_type_commit(decomp%ts_jj_iip2, ierr);  assert(ierr == MPI_SUCCESS)
    call mpi_type_vector(decomp%jj2, 1, decomp%ii + 3, MPI_DOUBLE_PRECISION, decomp%ts_jj2_iip3, ierr);  assert(ierr == MPI_SUCCESS)
    call mpi_type_commit(decomp%ts_jj2_iip3, ierr);  assert(ierr == MPI_SUCCESS)
  
    ! Data used to perform parallel data transposes
    if(decomp%p_ii > 1) then
      call allocate(decomp%trans_i_c, decomp, TRANSPOSE_II, target_decomp = decomp_i)
      decomp%trans_i_decomp => decomp%trans_i_c%decomp
    else
      nullify(decomp%trans_i_decomp)
    end if
    if(decomp%p_jj > 1) then
      call allocate(decomp%trans_j_c, decomp, TRANSPOSE_JJ, target_decomp = decomp_j)
      decomp%trans_j_decomp => decomp%trans_j_c%decomp
    else
      nullify(decomp%trans_j_decomp)
    end if
  
  end subroutine allocate_decomp
  
  ! Deallocate a domain decomposition.
  recursive subroutine deallocate_decomp(decomp)    
    type(decomposition), intent(inout) :: decomp
    
    integer :: ierr
        
    deallocate(decomp%iis)
    deallocate(decomp%jjs)
    deallocate(decomp%ii2s)
    deallocate(decomp%jj2s)
    
    call mpi_type_free(decomp%ts_jj_iip2, ierr);  assert(ierr == MPI_SUCCESS)
    call mpi_type_free(decomp%ts_jj2_iip3, ierr);  assert(ierr == MPI_SUCCESS)  
    
    if(decomp%p_ii > 1) call deallocate(decomp%trans_i_c)
    if(decomp%p_jj > 1) call deallocate(decomp%trans_j_c)
    
  end subroutine deallocate_decomp
  
  ! Print decomposition information.
  subroutine print_decomp(decomp)  
    type(decomposition), intent(in) :: decomp
    
    write(6, "(a)") "Decomposition:"
    write(6, "(a,i0,a,i0,a,i0)") "Global domain size: ", decomp%g_ii, " x ", decomp%g_jj, " x ", decomp%layers
    write(6, "(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)") "Local domain extents: ", decomp%ii_0 + 1, " - ", decomp%ii_0 + decomp%ii, " x ", decomp%jj_0 + 1, " - ", decomp%jj_0 + decomp%jj, " x ", 1, " - ", decomp%layers
    write(6, "(a,i0)") "Process: ", decomp%p + 1
    write(6, "(a,i0)") "Number of processes: ", decomp%n
    write(6, "(a,i0,a,i0)") "Logical coordinate: ", decomp%p_i, ", ", decomp%p_j
  
  end subroutine print_decomp
  
  ! Allocate a parallel data transpose cache.
  recursive subroutine allocate_transpose_cache(trans_c, donor_decomp, transpose_type, target_decomp)
    type(transpose_cache), intent(out) :: trans_c
    type(decomposition), intent(in) :: donor_decomp
    integer, intent(in) :: transpose_type
    type(decomposition), optional, intent(in) :: target_decomp
    
    select case(transpose_type)
      case(TRANSPOSE_II)
        call allocate_transpose_i_cache(trans_c, donor_decomp, target_decomp = target_decomp)
      case(TRANSPOSE_JJ)
        call allocate_transpose_j_cache(trans_c, donor_decomp, target_decomp = target_decomp)
      case default
        write(6, "(a,i0)") "Transpose type: ", transpose_type
        labort("Unrecognised transpose type")
    end select
    
  end subroutine allocate_transpose_cache
  
  ! Allocate a parallel data i-transpose cache.
  recursive subroutine allocate_transpose_i_cache(trans_c, donor_decomp, target_decomp)
    type(transpose_cache), intent(out) :: trans_c
    type(decomposition), intent(in) :: donor_decomp
    type(decomposition), target, optional, intent(in) :: target_decomp
    
    integer :: donor_ii, donor_jj, i, ierr, p_j, target_ii, target_jj
    integer, dimension(:), pointer :: donor_iis, donor_jjs, target_jjs
    type(decomposition), pointer :: ltarget_decomp
        
    trans_c%transpose_type = TRANSPOSE_II
    
    ! The number of processes to which we will send / from which we will receive
    trans_c%m = donor_decomp%p_ii
    if(trans_c%m == 1) then
      labort("Unable to allocate i-transpose data")
    end if
    ! The lowest rank process to which we will send
    trans_c%p_0 = (donor_decomp%p_j - 1) * trans_c%m
    
    ! Donor domain decomposition definition
    donor_iis => donor_decomp%iis
    donor_jjs => donor_decomp%jjs
    donor_ii = donor_decomp%ii
    donor_jj = donor_decomp%jj
    
    ! The target decomposition, decomposed with a slab decomposition contiguous
    ! in i
    if(present(target_decomp)) then
      assert(target_decomp%g_ii == donor_decomp%g_ii)
      assert(target_decomp%g_jj == donor_decomp%g_jj)
      assert(target_decomp%p_ii == 1)
      assert(target_decomp%p_jj == donor_decomp%n)
      trans_c%decomp => target_decomp
      
      trans_c%decomp_allocated = .false.
    else
      allocate(trans_c%decomp)
      ! Target domain decomposition definition
      allocate(target_jjs(donor_decomp%n))
      do p_j = 1, donor_decomp%p_jj
        target_jjs((p_j - 1) * trans_c%m + 1:p_j * trans_c%m) = donor_jjs(p_j) / trans_c%m
        do i = 1, donor_jjs(p_j) - target_jjs((p_j - 1) * trans_c%m + 1) * trans_c%m
          target_jjs(p_j * trans_c%m - i + 1) = target_jjs(p_j * trans_c%m - i + 1) + 1
        end do
      end do
      assert(sum(target_jjs) == donor_decomp%g_jj)
      
      if(donor_decomp%p_jj == 1) then
        call allocate(trans_c%decomp, &
          & iis = (/donor_decomp%g_ii/), jjs = target_jjs, &
          & layers = donor_decomp%layers, periodic_j = donor_decomp%periodic_j, &
          & decomp_j = donor_decomp)
      else
        call allocate(trans_c%decomp, &
          & iis = (/donor_decomp%g_ii/), jjs = target_jjs, &
          & layers = donor_decomp%layers, periodic_j = donor_decomp%periodic_j)
      end if
      
      deallocate(target_jjs)
      
      trans_c%decomp_allocated = .true.
    end if
    ltarget_decomp => trans_c%decomp
    target_jjs => ltarget_decomp%jjs
    target_ii = ltarget_decomp%ii
    target_jj = ltarget_decomp%jj
    
    ! Set up send / receive MPI types
    allocate(trans_c%lb_send_types(trans_c%m))
    allocate(trans_c%lb_receive_types(trans_c%m))
    allocate(trans_c%fb_send_types(trans_c%m))
    allocate(trans_c%fb_receive_types(trans_c%m))
    do p_j = 1, trans_c%m
      call mpi_type_vector(target_jjs(trans_c%p_0 + p_j), donor_ii, donor_ii, MPI_DOUBLE_PRECISION, trans_c%lb_send_types(p_j), ierr);  assert(ierr == MPI_SUCCESS)
      call mpi_type_commit(trans_c%lb_send_types(p_j), ierr);  assert(ierr == MPI_SUCCESS)
     
      call mpi_type_vector(target_jj, donor_iis(p_j), target_ii, MPI_DOUBLE_PRECISION, trans_c%lb_receive_types(p_j), ierr);  assert(ierr == MPI_SUCCESS)
      call mpi_type_commit(trans_c%lb_receive_types(p_j), ierr);  assert(ierr == MPI_SUCCESS)
      
      call mpi_type_vector(target_jjs(trans_c%p_0 + p_j), donor_ii, donor_ii + 2, MPI_DOUBLE_PRECISION, trans_c%fb_send_types(p_j), ierr);  assert(ierr == MPI_SUCCESS)
      call mpi_type_commit(trans_c%fb_send_types(p_j), ierr);  assert(ierr == MPI_SUCCESS)
     
      call mpi_type_vector(target_jj, donor_iis(p_j), target_ii + 2, MPI_DOUBLE_PRECISION, trans_c%fb_receive_types(p_j), ierr);  assert(ierr == MPI_SUCCESS)
      call mpi_type_commit(trans_c%fb_receive_types(p_j), ierr);  assert(ierr == MPI_SUCCESS)
    end do
    
    ! Find the starting indices in local and  full parallel data matrix
    allocate(trans_c%send_lb_i(1))
    allocate(trans_c%send_lb_j(trans_c%m))
    allocate(trans_c%receive_lb_i(trans_c%m))
    allocate(trans_c%receive_lb_j(1))
    allocate(trans_c%send_fb_i(1))
    allocate(trans_c%send_fb_j(trans_c%m))
    allocate(trans_c%receive_fb_i(trans_c%m))
    allocate(trans_c%receive_fb_j(1))
    trans_c%send_lb_i(1) = 1
    trans_c%send_lb_j(1) = 1
    trans_c%receive_lb_i(1) = 1
    trans_c%receive_lb_j(1) = 1
    do p_j = 2, trans_c%m
      trans_c%send_lb_j(p_j) = trans_c%send_lb_j(p_j - 1) + target_jjs(trans_c%p_0 + p_j - 1)
      trans_c%receive_lb_i(p_j) = trans_c%receive_lb_i(p_j - 1) + donor_iis(p_j - 1)
    end do
    trans_c%send_fb_i = trans_c%send_lb_i + 1
    trans_c%send_fb_j = trans_c%send_lb_j + 1
    trans_c%receive_fb_i = trans_c%receive_lb_i + 1
    trans_c%receive_fb_j = trans_c%receive_lb_j + 1
    
    ! Communicator
    call mpi_comm_dup(trans_c%decomp%cart_comm, trans_c%comm, ierr);  assert(ierr == MPI_SUCCESS)
    ! MPI requests / statuses
    allocate(trans_c%requests(2 * donor_decomp%layers * trans_c%m))
    trans_c%requests = MPI_REQUEST_NULL
    allocate(trans_c%statuses(2 * donor_decomp%layers * trans_c%m * MPI_STATUS_SIZE))
    
  end subroutine allocate_transpose_i_cache
  
  ! Allocate a parallel data j-transpose cache.
  recursive subroutine allocate_transpose_j_cache(trans_c, donor_decomp, target_decomp)
    type(transpose_cache), intent(out) :: trans_c
    type(decomposition), intent(in) :: donor_decomp
    type(decomposition), target, optional, intent(in) :: target_decomp
    
    integer :: donor_ii, donor_jj, i, ierr, p_i, target_ii
    integer, dimension(:), pointer :: donor_iis, donor_jjs, target_iis
    type(decomposition), pointer :: ltarget_decomp
        
    trans_c%transpose_type = TRANSPOSE_JJ
                
    ! The number of processes to which we will send / from which we will receive
    trans_c%m = donor_decomp%p_jj          
    if(trans_c%m == 1) then
      labort("Unable to allocate j-transpose data")
    end if
    ! The lowest rank process to which we will send / from which we will receive   
    trans_c%p_0 = (donor_decomp%p_i - 1) * trans_c%m  

    ! Donor domain decomposition definition
    donor_iis => donor_decomp%iis
    donor_jjs => donor_decomp%jjs
    donor_ii = donor_decomp%ii
    donor_jj = donor_decomp%jj
    
    ! The target decomposition, decomposed with a slab decomposition contiguous
    ! in j
    if(present(target_decomp)) then      
      assert(target_decomp%g_ii == donor_decomp%g_ii)
      assert(target_decomp%g_jj == donor_decomp%g_jj)
      assert(target_decomp%p_ii == donor_decomp%n)
      assert(target_decomp%p_jj == 1)
      trans_c%decomp => target_decomp
      
      trans_c%decomp_allocated = .false.
    else
      allocate(trans_c%decomp)
      ! Target domain decomposition definition
      allocate(target_iis(donor_decomp%n))
      do p_i = 1, donor_decomp%p_ii
        target_iis((p_i - 1) * trans_c%m + 1:p_i * trans_c%m) = donor_iis(p_i) / trans_c%m
        do i = 1, donor_iis(p_i) - target_iis((p_i - 1) * trans_c%m + 1) * trans_c%m
          target_iis(p_i * trans_c%m - i + 1) = target_iis(p_i * trans_c%m - i + 1) + 1
        end do
      end do
      assert(sum(target_iis) == donor_decomp%g_ii)
      
      if(donor_decomp%p_ii == 1) then
        call allocate(trans_c%decomp, &
          & iis = target_iis, jjs = (/donor_decomp%g_jj/), &
          & layers = donor_decomp%layers, periodic_j = donor_decomp%periodic_j, &
          & decomp_i = donor_decomp)
      else
        call allocate(trans_c%decomp, &
          & iis = target_iis, jjs = (/donor_decomp%g_jj/), &
          & layers = donor_decomp%layers, periodic_j = donor_decomp%periodic_j)
      end if
      
      deallocate(target_iis)
      
      trans_c%decomp_allocated = .true.
    end if
    ltarget_decomp => trans_c%decomp
    target_iis => ltarget_decomp%iis
    target_ii = ltarget_decomp%ii
     
    ! Set up send / receive MPI types
    allocate(trans_c%lb_send_types(trans_c%m))
    allocate(trans_c%lb_receive_types(trans_c%m))
    allocate(trans_c%fb_send_types(trans_c%m))
    allocate(trans_c%fb_receive_types(trans_c%m))
    do p_i = 1, trans_c%m
      call mpi_type_vector(donor_jj, target_iis(trans_c%p_0 + p_i), donor_ii, MPI_DOUBLE_PRECISION, trans_c%lb_send_types(p_i), ierr);  assert(ierr == MPI_SUCCESS)
      call mpi_type_commit(trans_c%lb_send_types(p_i), ierr);  assert(ierr == MPI_SUCCESS)
     
      call mpi_type_vector(donor_jjs(p_i), target_ii, target_ii, MPI_DOUBLE_PRECISION, trans_c%lb_receive_types(p_i), ierr);  assert(ierr == MPI_SUCCESS)
      call mpi_type_commit(trans_c%lb_receive_types(p_i), ierr);  assert(ierr == MPI_SUCCESS)
    
      call mpi_type_vector(donor_jj, target_iis(trans_c%p_0 + p_i), donor_ii + 2, MPI_DOUBLE_PRECISION, trans_c%fb_send_types(p_i), ierr);  assert(ierr == MPI_SUCCESS)
      call mpi_type_commit(trans_c%fb_send_types(p_i), ierr);  assert(ierr == MPI_SUCCESS)
     
      call mpi_type_vector(donor_jjs(p_i), target_ii, target_ii + 2, MPI_DOUBLE_PRECISION, trans_c%fb_receive_types(p_i), ierr);  assert(ierr == MPI_SUCCESS)
      call mpi_type_commit(trans_c%fb_receive_types(p_i), ierr);  assert(ierr == MPI_SUCCESS)
    end do
 
    ! Find the starting indices in local and full parallel data matrix
    allocate(trans_c%send_lb_i(trans_c%m))
    allocate(trans_c%send_lb_j(1))
    allocate(trans_c%receive_lb_i(1))
    allocate(trans_c%receive_lb_j(trans_c%m))
    allocate(trans_c%send_fb_i(trans_c%m))
    allocate(trans_c%send_fb_j(1))
    allocate(trans_c%receive_fb_i(1))
    allocate(trans_c%receive_fb_j(trans_c%m))
    trans_c%send_lb_i(1) = 1
    trans_c%send_lb_j(1) = 1
    trans_c%receive_lb_i(1) = 1
    trans_c%receive_lb_j(1) = 1
    do p_i = 2, trans_c%m
      trans_c%send_lb_i(p_i) = trans_c%send_lb_i(p_i - 1) + target_iis(trans_c%p_0 + p_i - 1)
      trans_c%receive_lb_j(p_i) = trans_c%receive_lb_j(p_i - 1) + donor_jjs(p_i - 1)
    end do
    trans_c%send_fb_i = trans_c%send_lb_i + 1
    trans_c%send_fb_j = trans_c%send_lb_j + 1
    trans_c%receive_fb_i = trans_c%receive_lb_i + 1
    trans_c%receive_fb_j = trans_c%receive_lb_j + 1
    
    ! Communicator
    call mpi_comm_dup(trans_c%decomp%cart_comm, trans_c%comm, ierr);  assert(ierr == MPI_SUCCESS)
    ! MPI requests / statuses
    allocate(trans_c%requests(2 * donor_decomp%layers * trans_c%m))
    trans_c%requests = MPI_REQUEST_NULL
    allocate(trans_c%statuses(2 * donor_decomp%layers * trans_c%m * MPI_STATUS_SIZE))
    
  end subroutine allocate_transpose_j_cache
  
  ! Deallocate a parallel data transpose cache.
  recursive subroutine deallocate_transpose_cache(trans_c)
    type(transpose_cache), intent(inout) :: trans_c
    
    integer :: i, ierr
    
    if(trans_c%decomp_allocated) then
      call deallocate(trans_c%decomp)
      deallocate(trans_c%decomp)
    end if
    
    do i = 1, size(trans_c%fb_send_types)
      call mpi_type_free(trans_c%lb_send_types(i), ierr);  assert(ierr == MPI_SUCCESS)
      call mpi_type_free(trans_c%lb_receive_types(i), ierr);  assert(ierr == MPI_SUCCESS)
      call mpi_type_free(trans_c%fb_send_types(i), ierr);  assert(ierr == MPI_SUCCESS)
      call mpi_type_free(trans_c%fb_receive_types(i), ierr);  assert(ierr == MPI_SUCCESS)
    end do
    deallocate(trans_c%lb_send_types)
    deallocate(trans_c%lb_receive_types)
    deallocate(trans_c%fb_send_types)
    deallocate(trans_c%fb_receive_types)
    
    deallocate(trans_c%send_lb_i)
    deallocate(trans_c%send_lb_j)
    deallocate(trans_c%receive_lb_i)
    deallocate(trans_c%receive_lb_j)
    deallocate(trans_c%send_fb_i)
    deallocate(trans_c%send_fb_j)
    deallocate(trans_c%receive_fb_i)
    deallocate(trans_c%receive_fb_j)
  
    call mpi_comm_free(trans_c%comm, ierr);  assert(ierr == MPI_SUCCESS)
    deallocate(trans_c%requests)
    deallocate(trans_c%statuses)
  
  end subroutine deallocate_transpose_cache
  
  ! Complete transpose communications.
  subroutine flush_transpose_cache(trans_c, layer)
    type(transpose_cache), intent(inout) :: trans_c
    integer, optional, intent(in) :: layer
    
    integer :: ierr
    
    if(present(layer)) then
      wtime0 = mpi_wtime()
      call mpi_waitall(2 * trans_c%m, trans_c%requests(2 * (layer - 1) * trans_c%m + 1), trans_c%statuses, ierr);  assert(ierr == MPI_SUCCESS)
      trans_c%requests(2 * (layer - 1) * trans_c%m + 1:2 * layer * trans_c%m) = MPI_REQUEST_NULL
      transpose_wait_time = transpose_wait_time + mpi_wtime() - wtime0
    else
      wtime0 = mpi_wtime()
      call mpi_waitall(size(trans_c%requests), trans_c%requests, trans_c%statuses, ierr);  assert(ierr == MPI_SUCCESS)
      trans_c%requests = MPI_REQUEST_NULL
      transpose_wait_time = transpose_wait_time + mpi_wtime() - wtime0
    end if
    
  end subroutine flush_transpose_cache
  
  ! Allocate a parallel data object.
  subroutine allocate_pdata_rank0(pdata, decomp, requests)
    type(parallel_data), intent(out) :: pdata
    type(decomposition), target, intent(in) :: decomp
    type(request_cache), target, optional, intent(in) :: requests
        
    pdata%decomp => decomp
    if(present(requests)) then
      pdata%requests => requests
    else
      pdata%requests => rcache
    end if
    
    allocate(pdata%data(decomp%ii + 2, decomp%jj + 2, decomp%layers))
    pdata%sdata => pdata%data(2:decomp%ii + 1, 2:decomp%jj + 1, :)  
    
  end subroutine allocate_pdata_rank0
  
  ! Allocate parallel data objects.
  subroutine allocate_pdata_rank1(pdata, decomp, requests)
    type(parallel_data), dimension(:), intent(out) :: pdata
    type(decomposition), intent(in) :: decomp
    type(request_cache), optional, intent(in) :: requests
    
    integer :: i
    
    do i = 1, size(pdata)
      call allocate(pdata(i), decomp, requests = requests)
    end do
    
  end subroutine allocate_pdata_rank1
  
  ! Deallocate a parallel data object.
  subroutine deallocate_pdata(pdata)
    type(parallel_data), intent(inout) :: pdata
    
    deallocate(pdata%data)
    
  end subroutine deallocate_pdata
  
  ! Deallocate parallel data objects.
  subroutine deallocate_pdata_rank1(pdata)
    type(parallel_data), dimension(:), intent(inout) :: pdata
    
    integer :: i
    
    do i = 1, size(pdata)
      call deallocate(pdata(i))
    end do
    
  end subroutine deallocate_pdata_rank1
  
  ! Zero data in a parallel data object.
  subroutine zero_pdata(pdata)
    type(parallel_data), intent(inout) :: pdata
    
    pdata%data = 0.0D0
    
  end subroutine zero_pdata
  
  ! Zero data in parallel data objects.
  subroutine zero_pdata_rank1(pdata)
    type(parallel_data), dimension(:), intent(inout) :: pdata
    
    integer :: i
    
    do i = 1, size(pdata)
      call zero(pdata(i))
    end do
    
  end subroutine zero_pdata_rank1
  
  ! Copy decomposed data.
  subroutine set_data_data(decomp, target, donor)
    type(decomposition), intent(in) :: decomp
    real(kind = real_8), dimension(:, :, :), intent(out) :: target
    real(kind = real_8), dimension(:, :, :), intent(in) :: donor
    
    integer :: i, j, layer
    
    assert(size(target, 1) >= decomp%ii + 2)
    assert(size(target, 2) >= decomp%jj + 2)
    assert(size(target, 3) == decomp%layers)
    assert(size(donor, 1) >= decomp%ii + 2)
    assert(size(donor, 2) >= decomp%jj + 2)
    assert(size(donor, 3) == decomp%layers)
    
    do layer = 1, decomp%layers
      do j = decomp%fa_j_min, decomp%fa_j_max
        do i = decomp%fa_i_min, decomp%fa_i_max
          target(i, j, layer) = donor(i, j, layer)
        end do
      end do
    end do
    
  end subroutine set_data_data
  
  ! Copy decomposed data.
  subroutine set_data_pdata(target, donor)
    type(parallel_data), intent(in) :: donor
    real(kind = real_8), dimension(:, :, :), intent(out) :: target
    
    integer :: i, j, layer
    type(decomposition), pointer :: decomp
    
    assert(size(target, 1) >= donor%decomp%ii + 2)
    assert(size(target, 2) >= donor%decomp%jj + 2)
    assert(size(target, 3) == donor%decomp%layers)
    
    decomp => donor%decomp
    
    do layer = 1, decomp%layers
      do j = decomp%fa_j_min, decomp%fa_j_max
        do i = decomp%fa_i_min, decomp%fa_i_max
          target(i, j, layer) = donor%data(i, j, layer)
        end do
      end do
    end do
    
  end subroutine set_data_pdata
  
  ! Set data in one parallel data object.
  subroutine set_pdata_real(target, donor, layer)
    type(parallel_data), intent(inout) :: target
    real(kind = real_8), intent(in) :: donor
    integer, optional, intent(in) :: layer
    
    if(present(layer)) then
      target%data(:, :, layer) = donor
    else
      target%data = donor
    end if
  
  end subroutine set_pdata_real
  
  ! Set data in one parallel data object.
  subroutine set_pdata_pdata(target, donor)
    type(parallel_data), intent(inout) :: target
    type(parallel_data), intent(in) :: donor
    
    assert(associated(target%decomp, donor%decomp))
    
    target%data = donor%data
  
  end subroutine set_pdata_pdata
  
  ! Set data in parallel data objects.
  subroutine set_pdata_pdata_rank1(target, donor)
    type(parallel_data), dimension(:), intent(inout) :: target
    type(parallel_data), dimension(size(target)), intent(in) :: donor
    
    integer :: i
    
    do i = 1, size(target)
      call set(target(i), donor(i))
    end do
  
  end subroutine set_pdata_pdata_rank1
  
  ! Add parallel data.
  subroutine add_pdata_pdata(target, donor)
    type(parallel_data), intent(inout) :: target
    type(parallel_data), intent(in) :: donor
    
    target%data = target%data + donor%data
  
  end subroutine add_pdata_pdata
  
  ! Subtract parallel data.
  subroutine subtract_pdata_pdata(target, donor)
    type(parallel_data), intent(inout) :: target
    type(parallel_data), intent(in) :: donor
    
    target%data = target%data - donor%data
  
  end subroutine subtract_pdata_pdata
  
  ! Initialise the supplied parallel data object with random data. The random
  ! data is generated only on the root process and distributed to all other
  ! processes (thereby enabling reproducibility between serial and parallel
  ! simulations). Note that this flushes the request cache associated with the
  ! parallel data object.
  subroutine rand_pdata_rank0(pdata, rscale, root)
    type(parallel_data), intent(inout) :: pdata
    real(kind = real_8), intent(in) :: rscale
    integer, optional, intent(in) :: root
    
    integer :: i, j, layer, lroot
    real(kind = real_8) :: r
    real(kind = real_8), dimension(:, :, :), allocatable :: data
    type(decomposition), pointer :: decomp
    
    if(present(root)) then
      assert(root >= 1 .and. root <= nprocs())
      lroot = root
    else
      lroot = 1
    end if
    
    decomp => pdata%decomp
          
    if(decomp%p + 1 == lroot) then
      allocate(data(decomp%g_ii, decomp%g_jj, decomp%layers))
      do layer = 1, decomp%layers
        do j = 1, decomp%g_jj
          do i = 1, decomp%g_ii
            call random_number(r)
            data(i, j, layer) = r * 2.0D0 * rscale - rscale
          end do
        end do
      end do
      
      call scatter(pdata, data)
      
      deallocate(data)
    else
      call scatter(pdata, lroot)
    end if    
    
  end subroutine rand_pdata_rank0
  
  ! Initialise the supplied parallel data objects with random data. The random
  ! data is generated only on the root process and distributed to all other
  ! processes (thereby enabling reproducibility between serial and parallel
  ! simulations). Note that this flushes the request caches associated with the
  ! parallel data objects.
  subroutine rand_pdata_rank1(pdata, rscale, root)
    type(parallel_data), dimension(:), intent(inout) :: pdata
    real(kind = real_8), intent(in) :: rscale
    integer, optional, intent(in) :: root
    
    integer :: i, j, layer, lroot, k
    real(kind = real_8) :: r
    real(kind = real_8), dimension(:, :, :, :), allocatable :: data
    type(decomposition), pointer :: decomp
    
#ifdef DDEBUG
    do i = 2, size(pdata)
      assert(associated(pdata(i)%decomp, pdata(1)%decomp))
    end do
#endif

    if(present(root)) then
      assert(root >= 1 .and. root <= nprocs())
      lroot = root
    else
      lroot = 1
    end if

    decomp => pdata(1)%decomp
    
    if(decomp%p + 1 == lroot) then
      allocate(data(decomp%g_ii, decomp%g_jj, decomp%layers, size(pdata)))
      
      do k = 1, size(pdata)
        do layer = 1, decomp%layers
          do j = 1, decomp%g_jj
            do i = 1, decomp%g_ii
              call random_number(r)
              data(i, j, layer, k) = r * 2.0D0 * rscale - rscale
            end do
          end do
        end do
      end do
      call scatter(pdata, data)
      
      deallocate(data)
    else
      call scatter(pdata, lroot)
    end if
    
  end subroutine rand_pdata_rank1
  
  ! Distribute a global array across all processes. Note that this flushes the
  ! request cache associated with the parallel data object.
  subroutine scatter_pdata_data_rank0_sender(pdata, data)
    type(parallel_data), intent(in) :: pdata
    real(kind = real_8), dimension(:, :, :), intent(in) :: data
    
    type(parallel_data), dimension(1) :: lpdata
    
    lpdata(1) = pdata
    call scatter(lpdata, spread(data, 4, 1))
    
  end subroutine scatter_pdata_data_rank0_sender
  
  ! Distribute a global array across all processes. Note that this flushes the
  ! request caches associated with the parallel data objects.
  subroutine scatter_pdata_data_rank1_sender(pdata, data)
    type(parallel_data), dimension(:), intent(inout) :: pdata
    real(kind = real_8), dimension(:, :, :, :), intent(in) :: data
    
    integer :: i, ierr, j, layer, p, receive_type, tag_base
    integer, dimension(2) :: cart_coord
    integer, dimension(pdata(1)%decomp%p_ii) :: ii_0s
    integer, dimension(pdata(1)%decomp%p_jj) :: jj_0s
    integer, dimension(pdata(1)%decomp%n) :: g_i_min, g_j_min, iis, jjs, &
      & send_types
    type(decomposition), pointer :: decomp

#ifdef DDEBUG
    do i = 2, size(pdata)
      assert(associated(pdata(i)%decomp, pdata(1)%decomp))
    end do
    assert(size(data, 1) == pdata(1)%decomp%g_ii)
    assert(size(data, 2) == pdata(1)%decomp%g_jj)
    assert(size(data, 3) == pdata(1)%decomp%layers)
    assert(size(data, 4) == size(pdata))
#endif
    
    decomp => pdata(1)%decomp
        
    ! The following could be cached in the decomposition. However, the global
    ! distribution is likely to dominate the cost here.
    ii_0s(1) = 0
    do i = 2, decomp%p_ii
      ii_0s(i) = ii_0s(i - 1) + decomp%iis(i - 1)
    end do
    jj_0s(1) = 0
    do j = 2, decomp%p_jj
      jj_0s(j) = jj_0s(j - 1) + decomp%jjs(j - 1)
    end do
    
    do i = 1, decomp%p_ii
      do j = 1, decomp%p_jj
        cart_coord = (/i - 1, j - 1/)
        call mpi_cart_rank(decomp%cart_comm, cart_coord, p, ierr);  assert(ierr == MPI_SUCCESS)
        p = p + 1
        
        iis(p) = decomp%iis(i)        
        jjs(p) = decomp%jjs(j)
        g_i_min(p) = ii_0s(i) + 1
        g_j_min(p) = jj_0s(j) + 1
      end do
    end do
  
    do i = 1, decomp%n
      call mpi_type_vector(jjs(i), iis(i), decomp%g_ii, MPI_DOUBLE_PRECISION, send_types(i), ierr);  assert(ierr == MPI_SUCCESS)
      call mpi_type_commit(send_types(i), ierr);  assert(ierr == MPI_SUCCESS)
    end do
    
    call mpi_type_vector(decomp%jj, decomp%ii, decomp%ii + 2, MPI_DOUBLE_PRECISION, receive_type, ierr);  assert(ierr == MPI_SUCCESS)
    call mpi_type_commit(receive_type, ierr);  assert(ierr == MPI_SUCCESS)
    
    do i = 1, size(pdata)
      do layer = 1, decomp%layers
        tag_base = pdata(i)%requests%tag_base
        do j = 1, decomp%n
          call mpi_isend(data(g_i_min(j), g_j_min(j), layer, i), 1, send_types(j), j - 1, tag_base, pdata(i)%requests%comm, next_request(pdata(i)%requests), ierr)
            assert(ierr == MPI_SUCCESS)
        end do
        call mpi_irecv(pdata(i)%data(2, 2, layer), 1, receive_type, decomp%p, tag_base, pdata(i)%requests%comm, next_request(pdata(i)%requests), ierr);  assert(ierr == MPI_SUCCESS)
        call increment_tag(pdata(i)%requests)
      end do
    end do    
    call flush(pdata)
    
    do i = 1, decomp%n
      call mpi_type_free(send_types(i), ierr);  assert(ierr == MPI_SUCCESS)
    end do
    call mpi_type_free(receive_type, ierr);  assert(ierr == MPI_SUCCESS)
    
    call update(pdata)
    call flush(pdata)
      
  end subroutine scatter_pdata_data_rank1_sender
  
  ! Distribute a global array across all processes. Note that this flushes the
  ! request cache associated with the parallel data object.
  subroutine scatter_pdata_data_rank0_receiver(pdata, root)
    type(parallel_data), intent(inout) :: pdata
    integer, intent(in) :: root
    
    type(parallel_data), dimension(1) :: lpdata
    
    lpdata(1) = pdata
    call scatter(lpdata, root)
    
  end subroutine scatter_pdata_data_rank0_receiver
  
  ! Distribute a global array across all processes. Note that this flushes the
  ! request caches associated with the parallel data objects.
  subroutine scatter_pdata_data_rank1_receiver(pdata, root)
    type(parallel_data), dimension(:), intent(inout) :: pdata
    integer, intent(in) :: root
    
    integer :: i, ierr, layer, receive_type, tag_base
    type(decomposition), pointer :: decomp
    
#ifdef DDEBUG
    do i = 2, size(pdata)
      assert(associated(pdata(i)%decomp, pdata(1)%decomp))
    end do
    assert(root >= 1 .and. root <= pdata(1)%decomp%n)
    assert(root /= pdata(1)%decomp%p + 1)
#endif
    
    decomp => pdata(1)%decomp
        
    ! The following could be cached in the decomposition. However, the global
    ! distribution is likely to dominate the cost here.
    call mpi_type_vector(decomp%jj, decomp%ii, decomp%ii + 2, MPI_DOUBLE_PRECISION, receive_type, ierr);  assert(ierr == MPI_SUCCESS)
    call mpi_type_commit(receive_type, ierr);  assert(ierr == MPI_SUCCESS)
    
    do i = 1, size(pdata)
      do layer = 1, decomp%layers
        tag_base = pdata(i)%requests%tag_base
        call mpi_irecv(pdata(i)%data(2, 2, layer), 1, receive_type, root - 1, tag_base, pdata(i)%requests%comm, next_request(pdata(i)%requests), ierr);  assert(ierr == MPI_SUCCESS)
        call increment_tag(pdata(i)%requests)
      end do
    end do    
    call flush(pdata)
    
    call mpi_type_free(receive_type, ierr);  assert(ierr == MPI_SUCCESS)
    
    call update(pdata)
    call flush(pdata)
      
  end subroutine scatter_pdata_data_rank1_receiver

  ! Gather a parallel data object onto one process. This flushes the supplied
  ! request cache (or rcache if no request cache is supplied).
  subroutine gather_pdata_data_rank0_receiver(pdata, data, requests)
    type(parallel_data), intent(in) :: pdata
    real(kind = real_8), dimension(:, :, :), intent(out) :: data
    type(request_cache), optional, intent(inout) :: requests

    real(kind = real_8), dimension(:, :, :, :), allocatable :: ldata

    allocate(ldata(size(data, 1), size(data, 2), size(data, 3), 1))
    call gather((/pdata/), ldata, requests = requests)
    data = ldata(:, :, :, 1)

  end subroutine gather_pdata_data_rank0_receiver

  ! Gather a parallel data object onto one process. This flushes the supplied
  ! request cache (or rcache if no request cache is supplied).
  subroutine gather_pdata_data_rank1_receiver(pdata, data, requests)
    type(parallel_data), dimension(:), intent(in) :: pdata
    real(kind = real_8), dimension(:, :, :, :), intent(out) :: data
    type(request_cache), target, optional, intent(inout) :: requests

    integer :: i, ierr, j, layer, p, send_type, tag_base
    integer, dimension(2) :: cart_coord
    integer, dimension(pdata(1)%decomp%p_ii) :: ii_0s
    integer, dimension(pdata(1)%decomp%p_jj) :: jj_0s
    integer, dimension(pdata(1)%decomp%n) :: g_i_min, g_j_min, iis, jjs, &
      & receive_types
    type(decomposition), pointer :: decomp
    type(request_cache), pointer :: lrequests

#ifdef DDEBUG
    do i = 2, size(pdata)
      assert(associated(pdata(i)%decomp, pdata(1)%decomp))
    end do
    assert(size(data, 1) == pdata(1)%decomp%g_ii)
    assert(size(data, 2) == pdata(1)%decomp%g_jj)
    assert(size(data, 3) == pdata(1)%decomp%layers)
    assert(size(data, 4) == size(pdata))
#endif
    
    if(present(requests)) then
      lrequests => requests
    else
      lrequests => rcache
    end if
    
    decomp => pdata(1)%decomp

    ! The following could be cached in the decomposition. However, the global
    ! communication is likely to dominate the cost here.
    ii_0s(1) = 0
    do i = 2, decomp%p_ii
      ii_0s(i) = ii_0s(i - 1) + decomp%iis(i - 1)
    end do
    jj_0s(1) = 0
    do j = 2, decomp%p_jj
      jj_0s(j) = jj_0s(j - 1) + decomp%jjs(j - 1)
    end do
    
    do i = 1, decomp%p_ii
      do j = 1, decomp%p_jj
        cart_coord = (/i - 1, j - 1/)
        call mpi_cart_rank(decomp%cart_comm, cart_coord, p, ierr);  assert(ierr == MPI_SUCCESS)
        p = p + 1
        
        iis(p) = decomp%iis(i)        
        jjs(p) = decomp%jjs(j)
        g_i_min(p) = ii_0s(i) + 1
        g_j_min(p) = jj_0s(j) + 1
      end do
    end do
  
    do i = 1, decomp%n
      call mpi_type_vector(jjs(i), iis(i), decomp%g_ii, MPI_DOUBLE_PRECISION, receive_types(i), ierr);  assert(ierr == MPI_SUCCESS)
      call mpi_type_commit(receive_types(i), ierr);  assert(ierr == MPI_SUCCESS)
    end do  
    call mpi_type_vector(decomp%jj, decomp%ii, decomp%ii + 2, MPI_DOUBLE_PRECISION, send_type, ierr);  assert(ierr == MPI_SUCCESS)
    call mpi_type_commit(send_type, ierr);  assert(ierr == MPI_SUCCESS)

    do i = 1, size(pdata)
      do layer = 1, decomp%layers
        tag_base = lrequests%tag_base
        do j = 1, decomp%n
          call mpi_irecv(data(g_i_min(j), g_j_min(j), layer, i), 1, receive_types(j), j - 1, tag_base, lrequests%comm, next_request(lrequests), ierr);  assert(ierr == MPI_SUCCESS)
        end do
        call mpi_isend(pdata(i)%data(2, 2, layer), 1, send_type, decomp%p, tag_base, lrequests%comm, next_request(lrequests), ierr);  assert(ierr == MPI_SUCCESS)
        call increment_tag(lrequests)
      end do
    end do    
    call flush(lrequests)

    do i = 1, decomp%n
      call mpi_type_free(receive_types(i), ierr);  assert(ierr == MPI_SUCCESS)
    end do
    call mpi_type_free(send_type, ierr);  assert(ierr == MPI_SUCCESS)

  end subroutine gather_pdata_data_rank1_receiver

  ! Gather a parallel data object onto one process. This flushes the supplied
  ! request cache (or rcache if no request cache is supplied).
  subroutine gather_pdata_data_rank0_sender(pdata, root, requests)
    type(parallel_data), intent(in) :: pdata
    integer, intent(in) :: root
    type(request_cache), optional, intent(inout) :: requests

    call gather((/pdata/), root, requests = requests)

  end subroutine gather_pdata_data_rank0_sender

  ! Gather a parallel data object onto one process. This flushes the supplied
  ! request cache (or rcache if no request cache is supplied).
  subroutine gather_pdata_data_rank1_sender(pdata, root, requests)
    type(parallel_data), dimension(:), intent(in) :: pdata
    integer, intent(in) :: root
    type(request_cache), target, optional, intent(inout) :: requests

    integer :: i, ierr, layer, send_type, tag_base
    type(decomposition), pointer :: decomp
    type(request_cache), pointer :: lrequests

#ifdef DDEBUG
    do i = 2, size(pdata)
      assert(associated(pdata(i)%decomp, pdata(1)%decomp))
    end do
    assert(root >= 1 .and. root <= pdata(1)%decomp%n)
    assert(root /= pdata(1)%decomp%p + 1)
#endif
    
    if(present(requests)) then
      lrequests => requests
    else
      lrequests => rcache
    end if
    
    decomp => pdata(1)%decomp

    ! The following could be cached in the decomposition. However, the global
    ! communication is likely to dominate the cost here.
    call mpi_type_vector(decomp%jj, decomp%ii, decomp%ii + 2, MPI_DOUBLE_PRECISION, send_type, ierr);  assert(ierr == MPI_SUCCESS)
    call mpi_type_commit(send_type, ierr);  assert(ierr == MPI_SUCCESS)

    do i = 1, size(pdata)
      do layer = 1, decomp%layers
        tag_base = lrequests%tag_base
        call mpi_isend(pdata(i)%data(2, 2, layer), 1, send_type, root - 1, tag_base, lrequests%comm, next_request(lrequests), ierr);  assert(ierr == MPI_SUCCESS)
        call increment_tag(lrequests)
      end do
    end do    
    call flush(lrequests)

    call mpi_type_free(send_type, ierr);  assert(ierr == MPI_SUCCESS)

  end subroutine gather_pdata_data_rank1_sender
    
  ! Update a parallel data object via communication. Note that this only sets up
  ! the non-blocking send / receives - the relevant request cache must be
  ! flushed for the data update to complete.
  subroutine update_pdata(pdata)  
    type(parallel_data), intent(inout) :: pdata
    
    integer :: ierr, layer, tag_base
    type(decomposition), pointer :: decomp
    type(request_cache), pointer :: requests
    
    decomp => pdata%decomp
    if(decomp%n == 1 .and. .not. decomp%periodic_j) return
    requests => pdata%requests
    
    do layer = 1, decomp%layers
      tag_base = requests%tag_base
    
      if(decomp%p_c >= 0) then
        ! Send / receive with logical left
        call mpi_isend(pdata%data(2, 2, layer), decomp%ii, MPI_DOUBLE_PRECISION, decomp%p_c, tag_base, requests%comm, next_request(requests), ierr);  assert(ierr == MPI_SUCCESS)
        call mpi_irecv(pdata%data(2, 1, layer), decomp%ii, MPI_DOUBLE_PRECISION, decomp%p_c, tag_base + 1, requests%comm, next_request(requests), ierr); assert(ierr == MPI_SUCCESS)
      
        if(decomp%p_cd >= 0) then
          ! Send / receive with logical lower left
          call mpi_isend(pdata%data(decomp%ii + 1, 2, layer), 1, MPI_DOUBLE_PRECISION, decomp%p_cd, tag_base + 4, requests%comm, next_request(requests), ierr)
            assert(ierr == MPI_SUCCESS)
          call mpi_irecv(pdata%data(decomp%ii + 2, 1, layer), 1, MPI_DOUBLE_PRECISION, decomp%p_cd, tag_base + 5, requests%comm, next_request(requests), ierr)
            assert(ierr == MPI_SUCCESS)
        end if
        if(decomp%p_bc >= 0) then
          ! Send / receive with logical upper left
          call mpi_isend(pdata%data(2, 2, layer), 1, MPI_DOUBLE_PRECISION, decomp%p_bc, tag_base + 6, requests%comm, next_request(requests), ierr);  assert(ierr == MPI_SUCCESS)
          call mpi_irecv(pdata%data(1, 1, layer), 1, MPI_DOUBLE_PRECISION, decomp%p_bc, tag_base + 7, requests%comm, next_request(requests), ierr);  assert(ierr == MPI_SUCCESS)
        end if
      end if
      if(decomp%p_a >= 0) then
        ! Send / receive with logical right
        call mpi_isend(pdata%data(2, decomp%jj + 1, layer), decomp%ii, MPI_DOUBLE_PRECISION, decomp%p_a, tag_base + 1, requests%comm, next_request(requests), ierr)
          assert(ierr == MPI_SUCCESS)
        call mpi_irecv(pdata%data(2, decomp%jj + 2, layer), decomp%ii, MPI_DOUBLE_PRECISION, decomp%p_a, tag_base, requests%comm, next_request(requests), ierr)
          assert(ierr == MPI_SUCCESS)
              
        if(decomp%p_da >= 0) then
          ! Send / receive with logical lower right
          call mpi_isend(pdata%data(decomp%ii + 1, decomp%jj + 1, layer), 1, MPI_DOUBLE_PRECISION, decomp%p_da, tag_base + 7, requests%comm, next_request(requests), ierr)
            assert(ierr == MPI_SUCCESS)
          call mpi_irecv(pdata%data(decomp%ii + 2, decomp%jj + 2, layer), 1, MPI_DOUBLE_PRECISION, decomp%p_da, tag_base + 6, requests%comm, next_request(requests), ierr)
            assert(ierr == MPI_SUCCESS)
        end if
        
        if(decomp%p_ab >= 0) then
          ! Send / receive with logical upper right
          call mpi_isend(pdata%data(2, decomp%jj + 1, layer), 1, MPI_DOUBLE_PRECISION, decomp%p_ab, tag_base + 5, requests%comm, next_request(requests), ierr)
            assert(ierr == MPI_SUCCESS)
          call mpi_irecv(pdata%data(1, decomp%jj + 2, layer), 1, MPI_DOUBLE_PRECISION, decomp%p_ab, tag_base + 4, requests%comm, next_request(requests), ierr)
            assert(ierr == MPI_SUCCESS)
        end if
      end if    
      if(decomp%p_d >= 0) then
        ! Send / receive with logical below
        call mpi_isend(pdata%data(decomp%ii + 1, 2, layer), 1, decomp%ts_jj_iip2, decomp%p_d, tag_base + 2, requests%comm, next_request(requests), ierr);  assert(ierr == MPI_SUCCESS)
        call mpi_irecv(pdata%data(decomp%ii + 2, 2, layer), 1, decomp%ts_jj_iip2, decomp%p_d, tag_base + 3, requests%comm, next_request(requests), ierr);  assert(ierr == MPI_SUCCESS)
      end if
      if(decomp%p_b >= 0) then
        ! Send / receive with logical above
        call mpi_isend(pdata%data(2, 2, layer), 1, decomp%ts_jj_iip2, decomp%p_b, tag_base + 3, requests%comm, next_request(requests), ierr);  assert(ierr == MPI_SUCCESS)
        call mpi_irecv(pdata%data(1, 2, layer), 1, decomp%ts_jj_iip2, decomp%p_b, tag_base + 2, requests%comm, next_request(requests), ierr);  assert(ierr == MPI_SUCCESS)
      end if
    
      call increment_tag(requests, increment = 8)
    end do
    
  end subroutine update_pdata
  
  ! Update parallel data objects via communication. Note that this only sets up
  ! the non-blocking send / receives - the relevant request cache must be
  ! flushed for the data update to complete.
  subroutine update_pdata_rank1(pdata)
    type(parallel_data), dimension(:), intent(inout) :: pdata
    
    integer :: i
    
    do i = 1, size(pdata)
      call update(pdata(i))
    end do
  
  end subroutine update_pdata_rank1

  ! Transpose a parallel data object. Note that this only sets up the
  ! non-blocking send / receives - the relevant request cache (attached to the
  ! *target*) must be flushed for the data update to complete.
  subroutine transpose_pdata(donor, target, transpose_type)
    type(parallel_data), intent(in) :: donor
    type(parallel_data), intent(inout) :: target
    integer, intent(in) :: transpose_type
   
    select case(transpose_type)
      case(TRANSPOSE_II)
        call transpose_pdata_i(donor, target)
      case(TRANSPOSE_JJ)
        call transpose_pdata_j(donor, target)
      case default
        write(6, "(a,i0)") "Transpose type: ", transpose_type
        labort("Unrecognised transpose type")
    end select
    
  end subroutine transpose_pdata
  
  ! Transpose parallel data. Note that this only sets up the non-blocking send
  ! / receives - the relevant request cache must be flushed for the data update
  ! to complete.
  subroutine transpose_data_layer(donor_decomp, donor, target, transpose_type, layer)
    type(decomposition), intent(in) :: donor_decomp
    real(kind = real_8), dimension(:, :, :), intent(in) :: donor
    real(kind = real_8), dimension(:, :, :), intent(out) :: target
    integer, intent(in) :: transpose_type
    integer, intent(in) :: layer
    
    select case(transpose_type)
      case(TRANSPOSE_II)
        if(size(donor, 1) == donor_decomp%ii) then
          assert(size(donor, 2) == donor_decomp%jj)
          assert(size(donor, 3) == donor_decomp%layers)
          assert(size(target, 1) == donor_decomp%trans_i_decomp%ii)
          assert(size(target, 2) == donor_decomp%trans_i_decomp%jj)
          assert(size(target, 3) == donor_decomp%trans_i_decomp%layers)
          call transpose_sdata_i(donor_decomp, donor, target, layer = layer)
        else
          assert(size(donor, 2) == donor_decomp%jj + 2)
          assert(size(donor, 3) == donor_decomp%layers)
          assert(size(target, 1) == donor_decomp%trans_i_decomp%ii + 2)
          assert(size(target, 2) == donor_decomp%trans_i_decomp%jj + 2)
          assert(size(target, 3) == donor_decomp%trans_i_decomp%layers)
          call transpose_data_i(donor_decomp, donor, target, layer = layer)
        end if
      case(TRANSPOSE_JJ)
        if(size(donor, 1) == donor_decomp%ii) then
          assert(size(donor, 2) == donor_decomp%jj)
          assert(size(donor, 3) == donor_decomp%layers)
          assert(size(target, 1) == donor_decomp%trans_j_decomp%ii)
          assert(size(target, 2) == donor_decomp%trans_j_decomp%jj)
          assert(size(target, 3) == donor_decomp%trans_j_decomp%layers)
          call transpose_sdata_j(donor_decomp, donor, target, layer = layer)
        else
          assert(size(donor, 2) == donor_decomp%jj + 2)
          assert(size(donor, 3) == donor_decomp%layers)
          assert(size(target, 1) == donor_decomp%trans_j_decomp%ii + 2)
          assert(size(target, 2) == donor_decomp%trans_j_decomp%jj + 2)
          assert(size(target, 3) == donor_decomp%trans_j_decomp%layers)
          call transpose_data_j(donor_decomp, donor, target, layer = layer)
        end if
      case default
        write(6, "(a,i0)") "Transpose type: ", transpose_type
        labort("Unrecognised transpose type")
    end select
    
  end subroutine transpose_data_layer
  
  ! Transpose parallel data. Note that this only sets up the non-blocking send
  ! / receives - the relevant requests must be flushed for the data update to
  ! complete.
  subroutine transpose_data_layers(donor_decomp, donor, target, transpose_type)
    type(decomposition), intent(in) :: donor_decomp
    real(kind = real_8), dimension(:, :, :), intent(in) :: donor
    real(kind = real_8), dimension(:, :, :), intent(out) :: target
    integer, intent(in) :: transpose_type
    
    select case(transpose_type)
      case(TRANSPOSE_II)
        if(size(donor, 1) == donor_decomp%ii) then
          assert(size(donor, 2) == donor_decomp%jj)
          assert(size(donor, 3) == donor_decomp%layers)
          assert(size(target, 1) == donor_decomp%trans_i_decomp%ii)
          assert(size(target, 2) == donor_decomp%trans_i_decomp%jj)
          assert(size(target, 3) == donor_decomp%trans_i_decomp%layers)
          call transpose_sdata_i(donor_decomp, donor, target)
        else
          assert(size(donor, 2) == donor_decomp%jj + 2)
          assert(size(donor, 3) == donor_decomp%layers)
          assert(size(target, 1) == donor_decomp%trans_i_decomp%ii + 2)
          assert(size(target, 2) == donor_decomp%trans_i_decomp%jj + 2)
          assert(size(target, 3) == donor_decomp%trans_i_decomp%layers)
          call transpose_data_i(donor_decomp, donor, target)
        end if
      case(TRANSPOSE_JJ)
        if(size(donor, 1) == donor_decomp%ii) then
          assert(size(donor, 2) == donor_decomp%jj)
          assert(size(donor, 3) == donor_decomp%layers)
          assert(size(target, 1) == donor_decomp%trans_j_decomp%ii)
          assert(size(target, 2) == donor_decomp%trans_j_decomp%jj)
          assert(size(target, 3) == donor_decomp%trans_j_decomp%layers)
          call transpose_sdata_j(donor_decomp, donor, target)
        else
          assert(size(donor, 2) == donor_decomp%jj + 2)
          assert(size(donor, 3) == donor_decomp%layers)
          assert(size(target, 1) == donor_decomp%trans_j_decomp%ii + 2)
          assert(size(target, 2) == donor_decomp%trans_j_decomp%jj + 2)
          assert(size(target, 3) == donor_decomp%trans_j_decomp%layers)
          call transpose_data_j(donor_decomp, donor, target)
        end if
      case default
        write(6, "(a,i0)") "Transpose type: ", transpose_type
        labort("Unrecognised transpose type")
    end select
    
  end subroutine transpose_data_layers
  
  ! Transpose decomposed data to a slab decomposition contiguous in i.
  subroutine transpose_pdata_i(donor, target)
    type(parallel_data), intent(in) :: donor
    type(parallel_data), intent(inout) :: target
    
    assert(associated(donor%decomp%trans_i_decomp))
    assert(associated(target%decomp, donor%decomp%trans_i_decomp))
    
    call transpose_data_i(donor%decomp, donor%data, target%data)
    
  end subroutine transpose_pdata_i
  
  ! Transpose a local (owned) data matrix to a slab decomposition contiguous in
  ! i.
  subroutine transpose_sdata_i_layer(donor_decomp, donor, target, layer)
    type(decomposition), target, intent(in) :: donor_decomp
    real(kind = real_8), dimension(donor_decomp%ii, donor_decomp%jj, donor_decomp%layers), intent(in) :: donor
    real(kind = real_8), dimension(donor_decomp%trans_i_decomp%ii, donor_decomp%trans_i_decomp%jj, donor_decomp%trans_i_decomp%layers), intent(out) :: target
    integer, intent(in) :: layer
    
    integer :: ierr, p, p_j, tag
    integer, dimension(:), pointer :: requests
    type(decomposition), pointer :: target_decomp
    type(transpose_cache), pointer :: trans_c
    
    assert(layer >= 1 .and. layer <= donor_decomp%layers)
    
    trans_c => donor_decomp%trans_i_c
    target_decomp => donor_decomp%trans_i_decomp
    
    do p_j = 1, trans_c%m
      requests => trans_c%requests(2 * (p_j + (layer - 1) * trans_c%m) - 1:2 * (p_j + (layer - 1) * trans_c%m))
      assert(all(requests == MPI_REQUEST_NULL))
      tag = (layer - 1)
      
      p = trans_c%p_0 + p_j - 1
      call mpi_isend(donor(1, trans_c%send_lb_j(p_j), layer), 1, trans_c%lb_send_types(p_j), p, tag, trans_c%comm, requests(1), ierr);  assert(ierr == MPI_SUCCESS)
      p = ((target_decomp%p_j - 1) / trans_c%m) + (p_j - 1) * donor_decomp%p_jj
      call mpi_irecv(target(trans_c%receive_lb_i(p_j), 1, layer), 1, trans_c%lb_receive_types(p_j), p, tag, trans_c%comm, requests(2), ierr);  assert(ierr == MPI_SUCCESS)
    end do
  
  end subroutine transpose_sdata_i_layer
  
  ! Transpose a local (owned) data matrix to a slab decomposition contiguous in
  ! i.
  subroutine transpose_sdata_i_layers(donor_decomp, donor, target)
    type(decomposition), target, intent(in) :: donor_decomp
    real(kind = real_8), dimension(donor_decomp%ii, donor_decomp%jj, donor_decomp%layers), intent(in) :: donor
    real(kind = real_8), dimension(donor_decomp%trans_i_decomp%ii, donor_decomp%trans_i_decomp%jj, donor_decomp%trans_i_decomp%layers), intent(out) :: target
    
    integer :: layer
    
    do layer = 1, donor_decomp%layers
      call transpose_sdata_i(donor_decomp, donor, target, layer = layer)
    end do
  
  end subroutine transpose_sdata_i_layers
  
  ! Transpose a full (owned + halo) data matrix to a slab decomposition
  ! contiguous in i.
  subroutine transpose_data_i_layer(donor_decomp, donor, target, layer)
    type(decomposition), target, intent(in) :: donor_decomp
    real(kind = real_8), dimension(donor_decomp%ii + 2, donor_decomp%jj + 2, donor_decomp%layers), intent(in) :: donor
    real(kind = real_8), dimension(donor_decomp%trans_i_decomp%ii + 2, donor_decomp%trans_i_decomp%jj + 2, donor_decomp%trans_i_decomp%layers), intent(out) :: target
    integer, intent(in) :: layer
    
    integer :: ierr, p, p_j, tag
    integer, dimension(:), pointer :: requests
    type(decomposition), pointer :: target_decomp
    type(transpose_cache), pointer :: trans_c
    
    trans_c => donor_decomp%trans_i_c
    target_decomp => donor_decomp%trans_i_decomp
    
    do p_j = 1, trans_c%m
      requests => trans_c%requests(2 * (p_j + (layer - 1) * trans_c%m) - 1:2 * (p_j + (layer - 1) * trans_c%m))
      assert(all(requests == MPI_REQUEST_NULL))
      tag = (layer - 1) + target_decomp%layers
    
      p = trans_c%p_0 + p_j - 1
      call mpi_isend(donor(2, trans_c%send_fb_j(p_j), layer), 1, trans_c%fb_send_types(p_j), p, tag, trans_c%comm, requests(1), ierr);  assert(ierr == MPI_SUCCESS)
      p = ((target_decomp%p_j - 1) / trans_c%m) + (p_j - 1) * donor_decomp%p_jj
      call mpi_irecv(target(trans_c%receive_fb_i(p_j), 2, layer), 1, trans_c%fb_receive_types(p_j), p, tag, trans_c%comm, requests(2), ierr);  assert(ierr == MPI_SUCCESS)
    end do
  
  end subroutine transpose_data_i_layer
  
  ! Transpose a full (owned + halo) data matrix to a slab decomposition
  ! contiguous in i.
  subroutine transpose_data_i_layers(donor_decomp, donor, target)
    type(decomposition), target, intent(in) :: donor_decomp
    real(kind = real_8), dimension(donor_decomp%ii + 2, donor_decomp%jj + 2, donor_decomp%layers), intent(in) :: donor
    real(kind = real_8), dimension(donor_decomp%trans_i_decomp%ii + 2, donor_decomp%trans_i_decomp%jj + 2, donor_decomp%trans_i_decomp%layers), intent(out) :: target
    
    integer :: layer
    
    do layer = 1, donor_decomp%layers
      call transpose_data_i(donor_decomp, donor, target, layer = layer)
    end do
  
  end subroutine transpose_data_i_layers
  
  ! Transpose decomposed data to a slab decomposition contiguous in j.
  subroutine transpose_pdata_j(donor, target)
    type(parallel_data), intent(in) :: donor
    type(parallel_data), intent(inout) :: target
    
    assert(associated(donor%decomp%trans_j_decomp))
    assert(associated(target%decomp, donor%decomp%trans_j_decomp))
    
    call transpose_data_j(donor%decomp, donor%data, target%data)
    
  end subroutine transpose_pdata_j
  
  ! Transpose a local (owned) data matrix to a slab decomposition contiguous in
  ! j.
  subroutine transpose_sdata_j_layer(donor_decomp, donor, target, layer)
    type(decomposition), target, intent(in) :: donor_decomp
    real(kind = real_8), dimension(donor_decomp%ii, donor_decomp%jj, donor_decomp%layers), intent(in) :: donor
    real(kind = real_8), dimension(donor_decomp%trans_j_decomp%ii, donor_decomp%trans_j_decomp%jj, donor_decomp%trans_j_decomp%layers), intent(out) :: target
    integer, intent(in) :: layer
  
    integer :: ierr, p, p_i, tag
    integer, dimension(:), pointer :: requests
    type(decomposition), pointer :: target_decomp
    type(transpose_cache), pointer :: trans_c
    
    assert(layer >= 1 .and. layer <= donor_decomp%layers)
    
    trans_c => donor_decomp%trans_j_c
    target_decomp => donor_decomp%trans_j_decomp
      
    do p_i = 1, trans_c%m
      requests => trans_c%requests(2 * (p_i + (layer - 1) * trans_c%m) - 1:2 * (p_i + (layer - 1) * trans_c%m))
      assert(all(requests == MPI_REQUEST_NULL))
      tag = (layer - 1) + 2 * target_decomp%layers
      
      p = trans_c%p_0 + p_i - 1
      call mpi_isend(donor(trans_c%send_lb_i(p_i), 1, layer), 1, trans_c%lb_send_types(p_i), p, tag, trans_c%comm, requests(1), ierr);  assert(ierr == MPI_SUCCESS)
      call mpi_irecv(target(1, trans_c%receive_lb_j(p_i), layer), 1, trans_c%lb_receive_types(p_i), p, tag, trans_c%comm, requests(2), ierr);  assert(ierr == MPI_SUCCESS)
    end do
    
  end subroutine transpose_sdata_j_layer
  
  ! Transpose a local (owned) data matrix to a slab decomposition contiguous in
  ! j.
  subroutine transpose_sdata_j_layers(donor_decomp, donor, target)
    type(decomposition), target, intent(in) :: donor_decomp
    real(kind = real_8), dimension(donor_decomp%ii, donor_decomp%jj, donor_decomp%layers), intent(in) :: donor
    real(kind = real_8), dimension(donor_decomp%trans_j_decomp%ii, donor_decomp%trans_j_decomp%jj, donor_decomp%trans_j_decomp%layers), intent(out) :: target
  
    integer :: layer
    
    do layer = 1, donor_decomp%layers
      call transpose_sdata_j(donor_decomp, donor, target, layer = layer)
    end do
    
  end subroutine transpose_sdata_j_layers
  
  ! Transpose a full (owned + halo) data matrix to a slab decomposition
  ! contiguous in j.
  subroutine transpose_data_j_layer(donor_decomp, donor, target, layer)
    type(decomposition), target, intent(in) :: donor_decomp
    real(kind = real_8), dimension(donor_decomp%ii + 2, donor_decomp%jj + 2, donor_decomp%layers), intent(in) :: donor
    real(kind = real_8), dimension(donor_decomp%trans_j_decomp%ii + 2, donor_decomp%trans_j_decomp%jj + 2, donor_decomp%trans_j_decomp%layers), intent(out) :: target
    integer, intent(in) :: layer
  
    integer :: ierr, p, p_i, tag
    integer, dimension(:), pointer :: requests
    type(decomposition), pointer :: target_decomp
    type(transpose_cache), pointer :: trans_c
    
    trans_c => donor_decomp%trans_j_c
    target_decomp => donor_decomp%trans_j_decomp
      
    do p_i = 1, trans_c%m
      requests => trans_c%requests(2 * (p_i + (layer - 1) * trans_c%m) - 1:2 * (p_i + (layer - 1) * trans_c%m))
      assert(all(requests == MPI_REQUEST_NULL))
      tag = (layer - 1) + 3 * target_decomp%layers
    
      p = trans_c%p_0 + p_i - 1
      call mpi_isend(donor(trans_c%send_fb_i(p_i), 2, layer), 1, trans_c%fb_send_types(p_i), p, tag, trans_c%comm, requests(1), ierr);  assert(ierr == MPI_SUCCESS)
      call mpi_irecv(target(2, trans_c%receive_fb_j(p_i), layer), 1, trans_c%fb_receive_types(p_i), p, tag, trans_c%comm, requests(2), ierr);  assert(ierr == MPI_SUCCESS)
    end do
    
  end subroutine transpose_data_j_layer
  
  ! Transpose a full (owned + halo) data matrix to a slab decomposition
  ! contiguous in j.
  subroutine transpose_data_j_layers(donor_decomp, donor, target)
    type(decomposition), target, intent(in) :: donor_decomp
    real(kind = real_8), dimension(donor_decomp%ii + 2, donor_decomp%jj + 2, donor_decomp%layers), intent(in) :: donor
    real(kind = real_8), dimension(donor_decomp%trans_j_decomp%ii + 2, donor_decomp%trans_j_decomp%jj + 2, donor_decomp%trans_j_decomp%layers), intent(out) :: target
  
    integer :: layer
    type(decomposition), pointer :: target_decomp
    type(transpose_cache), pointer :: trans_c
    
    trans_c => donor_decomp%trans_j_c
    target_decomp => donor_decomp%trans_j_decomp
    
    do layer = 1, donor_decomp%layers    
      call transpose_data_j(donor_decomp, donor, target, layer = layer)
    end do
    
  end subroutine transpose_data_j_layers
  
  ! Perform the inverse operation to transpose_pdata.
  subroutine inverse_transpose_pdata(donor, target, transpose_type)
    type(parallel_data), intent(in) :: donor
    type(parallel_data), intent(inout) :: target
    integer, intent(in) :: transpose_type
    
    select case(transpose_type)
      case(TRANSPOSE_II)
        call inverse_transpose_pdata_i(donor, target)
      case(TRANSPOSE_JJ)
        call inverse_transpose_pdata_j(donor, target)
      case default
        write(6, "(a,i0)") "Transpose type: ", transpose_type
        labort("Unrecognised transpose type")
    end select

  end subroutine inverse_transpose_pdata
  
  ! Perform the inverse operation to transpose_data.
  subroutine inverse_transpose_data_layer(target_decomp, donor, target, transpose_type, layer)
    type(decomposition), intent(in) :: target_decomp
    real(kind = real_8), dimension(:, :, :), intent(in) :: donor
    real(kind = real_8), dimension(:, :, :), intent(out) :: target
    integer, intent(in) :: transpose_type
    integer, intent(in) :: layer
    
    select case(transpose_type)
      case(TRANSPOSE_II)
        if(size(target, 1) == target_decomp%ii) then
          assert(size(target, 2) == target_decomp%jj)
          assert(size(target, 3) == target_decomp%layers)
          assert(size(donor, 1) == target_decomp%trans_i_decomp%ii)
          assert(size(donor, 2) == target_decomp%trans_i_decomp%jj)
          assert(size(donor, 3) == target_decomp%trans_i_decomp%layers)
          call inverse_transpose_sdata_i(target_decomp, donor, target, layer = layer)
        else
          assert(size(target, 2) == target_decomp%jj + 2)
          assert(size(target, 3) == target_decomp%layers)
          assert(size(donor, 1) == target_decomp%trans_i_decomp%ii + 2)
          assert(size(donor, 2) == target_decomp%trans_i_decomp%jj + 2)
          assert(size(donor, 3) == target_decomp%trans_i_decomp%layers)
          call inverse_transpose_data_i(target_decomp, donor, target, layer = layer)
        end if
      case(TRANSPOSE_JJ)
        if(size(target, 1) == target_decomp%ii) then
          assert(size(target, 2) == target_decomp%jj)
          assert(size(target, 3) == target_decomp%layers)
          assert(size(donor, 1) == target_decomp%trans_j_decomp%ii)
          assert(size(donor, 2) == target_decomp%trans_j_decomp%jj)
          assert(size(donor, 3) == target_decomp%trans_j_decomp%layers)
          call inverse_transpose_sdata_j(target_decomp, donor, target, layer = layer)
        else
          assert(size(target, 2) == target_decomp%jj + 2)
          assert(size(target, 3) == target_decomp%layers)
          assert(size(donor, 1) == target_decomp%trans_j_decomp%ii + 2)
          assert(size(donor, 2) == target_decomp%trans_j_decomp%jj + 2)
          assert(size(donor, 3) == target_decomp%trans_j_decomp%layers)
          call inverse_transpose_data_j(target_decomp, donor, target, layer = layer)
        end if
      case default
        write(6, "(a,i0)") "Transpose type: ", transpose_type
        labort("Unrecognised transpose type")
    end select
    
  end subroutine inverse_transpose_data_layer
  
  ! Perform the inverse operation to transpose_data.
  subroutine inverse_transpose_data_layers(target_decomp, donor, target, transpose_type)
    type(decomposition), intent(in) :: target_decomp
    real(kind = real_8), dimension(:, :, :), intent(in) :: donor
    real(kind = real_8), dimension(:, :, :), intent(out) :: target
    integer, intent(in) :: transpose_type
    
    select case(transpose_type)
      case(TRANSPOSE_II)
        if(size(target, 1) == target_decomp%ii) then
          assert(size(target, 2) == target_decomp%jj)
          assert(size(target, 3) == target_decomp%layers)
          assert(size(donor, 1) == target_decomp%trans_i_decomp%ii)
          assert(size(donor, 2) == target_decomp%trans_i_decomp%jj)
          assert(size(donor, 3) == target_decomp%trans_i_decomp%layers)
          call inverse_transpose_sdata_i(target_decomp, donor, target)
        else
          assert(size(target, 2) == target_decomp%jj + 2)
          assert(size(target, 3) == target_decomp%layers)
          assert(size(donor, 1) == target_decomp%trans_i_decomp%ii + 2)
          assert(size(donor, 2) == target_decomp%trans_i_decomp%jj + 2)
          assert(size(donor, 3) == target_decomp%trans_i_decomp%layers)
          call inverse_transpose_data_i(target_decomp, donor, target)
        end if
      case(TRANSPOSE_JJ)
        if(size(target, 1) == target_decomp%ii) then
          assert(size(target, 2) == target_decomp%jj)
          assert(size(target, 3) == target_decomp%layers)
          assert(size(donor, 1) == target_decomp%trans_j_decomp%ii)
          assert(size(donor, 2) == target_decomp%trans_j_decomp%jj)
          assert(size(donor, 3) == target_decomp%trans_j_decomp%layers)
          call inverse_transpose_sdata_j(target_decomp, donor, target)
        else
          assert(size(target, 2) == target_decomp%jj + 2)
          assert(size(target, 3) == target_decomp%layers)
          assert(size(donor, 1) == target_decomp%trans_j_decomp%ii + 2)
          assert(size(donor, 2) == target_decomp%trans_j_decomp%jj + 2)
          assert(size(donor, 3) == target_decomp%trans_j_decomp%layers)
          call inverse_transpose_data_j(target_decomp, donor, target)
        end if
      case default
        write(6, "(a,i0)") "Transpose type: ", transpose_type
        labort("Unrecognised transpose type")
    end select
    
  end subroutine inverse_transpose_data_layers
  
  ! Perform the inverse operation to transpose_pdata_i.
  subroutine inverse_transpose_pdata_i(donor, target)
    type(parallel_data), intent(in) :: donor
    type(parallel_data), intent(inout) :: target
    
    assert(associated(target%decomp%trans_i_decomp))
    assert(associated(donor%decomp, target%decomp%trans_i_decomp))
    
    call inverse_transpose_data_i(target%decomp, donor%data, target%data)
    
  end subroutine inverse_transpose_pdata_i
  
  ! Perform the inverse operation to transpose_sdata_i.
  subroutine inverse_transpose_sdata_i_layer(target_decomp, donor, target, layer)
    type(decomposition), target, intent(in) :: target_decomp
    real(kind = real_8), dimension(target_decomp%trans_i_decomp%ii, target_decomp%trans_i_decomp%jj, target_decomp%trans_i_decomp%layers), intent(in) :: donor
    real(kind = real_8), dimension(target_decomp%ii, target_decomp%jj, target_decomp%layers), intent(out) :: target
    integer, intent(in) :: layer

    integer :: ierr, p, p_j, tag
    integer, dimension(:), pointer :: requests
    type(decomposition), pointer :: donor_decomp
    type(transpose_cache), pointer :: trans_c

    trans_c => target_decomp%trans_i_c
    donor_decomp => target_decomp%trans_i_decomp
    
    do p_j = 1, trans_c%m
      requests => trans_c%requests(2 * (p_j + (layer - 1) * trans_c%m) - 1:2 * (p_j + (layer - 1) * trans_c%m))
      assert(all(requests == MPI_REQUEST_NULL))
      tag = (layer - 1) + 4 * target_decomp%layers
    
      p = ((donor_decomp%p_j - 1) / trans_c%m) + (p_j - 1) * target_decomp%p_jj
      call mpi_isend(donor(trans_c%receive_lb_i(p_j), 1, layer), 1, trans_c%lb_receive_types(p_j), p, tag, trans_c%comm, requests(1), ierr);  assert(ierr == MPI_SUCCESS)
      p = trans_c%p_0 + p_j - 1
      call mpi_irecv(target(1, trans_c%send_lb_j(p_j), layer), 1, trans_c%lb_send_types(p_j), p, tag, trans_c%comm, requests(2), ierr);  assert(ierr == MPI_SUCCESS)
    end do

  end subroutine inverse_transpose_sdata_i_layer
  
  ! Perform the inverse operation to transpose_sdata_i.
  subroutine inverse_transpose_sdata_i_layers(target_decomp, donor, target)
    type(decomposition), target, intent(in) :: target_decomp
    real(kind = real_8), dimension(target_decomp%trans_i_decomp%ii, target_decomp%trans_i_decomp%jj, target_decomp%trans_i_decomp%layers), intent(in) :: donor
    real(kind = real_8), dimension(target_decomp%ii, target_decomp%jj, target_decomp%layers), intent(out) :: target

    integer :: layer
    
    do layer = 1, target_decomp%layers
      call inverse_transpose_sdata_i(target_decomp, donor, target, layer = layer)
    end do

  end subroutine inverse_transpose_sdata_i_layers
  
  ! Perform the inverse operation to transpose_data_i.
  subroutine inverse_transpose_data_i_layer(target_decomp, donor, target, layer)
    type(decomposition), target, intent(in) :: target_decomp
    real(kind = real_8), dimension(target_decomp%trans_i_decomp%ii + 2, target_decomp%trans_i_decomp%jj + 2, target_decomp%trans_i_decomp%layers), intent(in) :: donor
    real(kind = real_8), dimension(target_decomp%ii + 2, target_decomp%jj + 2, target_decomp%layers), intent(out) :: target
    integer, intent(in) :: layer

    integer :: ierr, p, p_j, tag
    integer, dimension(:), pointer :: requests
    type(decomposition), pointer :: donor_decomp
    type(transpose_cache), pointer :: trans_c

    trans_c => target_decomp%trans_i_c
    donor_decomp => target_decomp%trans_i_decomp
    
    do p_j = 1, trans_c%m
      requests => trans_c%requests(2 * (p_j + (layer - 1) * trans_c%m) - 1:2 * (p_j + (layer - 1) * trans_c%m))
      assert(all(requests == MPI_REQUEST_NULL))
      tag = (layer - 1) + 5 * target_decomp%layers
    
      p = ((donor_decomp%p_j - 1) / trans_c%m) + (p_j - 1) * target_decomp%p_jj
      call mpi_isend(donor(trans_c%receive_fb_i(p_j), 2, layer), 1, trans_c%fb_receive_types(p_j), p, tag, trans_c%comm, requests(1), ierr);  assert(ierr == MPI_SUCCESS)
      p = trans_c%p_0 + p_j - 1
      call mpi_irecv(target(2, trans_c%send_fb_j(p_j), layer), 1, trans_c%fb_send_types(p_j), p, tag, trans_c%comm, requests(2), ierr);  assert(ierr == MPI_SUCCESS)
    end do

  end subroutine inverse_transpose_data_i_layer
  
  ! Perform the inverse operation to transpose_data_i.
  subroutine inverse_transpose_data_i_layers(target_decomp, donor, target)
    type(decomposition), target, intent(in) :: target_decomp
    real(kind = real_8), dimension(target_decomp%trans_i_decomp%ii + 2, target_decomp%trans_i_decomp%jj + 2, target_decomp%trans_i_decomp%layers), intent(in) :: donor
    real(kind = real_8), dimension(target_decomp%ii + 2, target_decomp%jj + 2, target_decomp%layers), intent(out) :: target

    integer :: layer
    
    do layer = 1, target_decomp%layers      
      call inverse_transpose_data_i(target_decomp, donor, target, layer = layer)
    end do

  end subroutine inverse_transpose_data_i_layers
  
  ! Perform the inverse operation to transpose_pdata_j.
  subroutine inverse_transpose_pdata_j(donor, target)
    type(parallel_data), intent(in) :: donor
    type(parallel_data), intent(inout) :: target
    
    assert(associated(target%decomp%trans_j_decomp))
    assert(associated(donor%decomp, target%decomp%trans_j_decomp))
    
    call inverse_transpose_data_j(target%decomp, donor%data, target%data)
    
  end subroutine inverse_transpose_pdata_j
  
  ! Perform the inverse operation to transpose_sdata_j.
  subroutine inverse_transpose_sdata_j_layer(target_decomp, donor, target, layer)
    type(decomposition), target, intent(in) :: target_decomp
    real(kind = real_8), dimension(target_decomp%trans_j_decomp%ii, target_decomp%trans_j_decomp%jj, target_decomp%trans_j_decomp%layers), intent(in) :: donor
    real(kind = real_8), dimension(target_decomp%ii, target_decomp%jj, target_decomp%layers), intent(out) :: target
    integer, intent(in) :: layer

    integer :: ierr, p, p_i, tag
    integer, dimension(:), pointer :: requests
    type(decomposition), pointer :: donor_decomp
    type(transpose_cache), pointer :: trans_c

    trans_c => target_decomp%trans_j_c
    donor_decomp => target_decomp%trans_j_decomp
    
    do p_i = 1, trans_c%m
      requests => trans_c%requests(2 * (p_i + (layer - 1) * trans_c%m) - 1:2 * (p_i + (layer - 1) * trans_c%m))
      assert(all(requests == MPI_REQUEST_NULL))
      tag = (layer - 1) + 6 * target_decomp%layers
    
      p = trans_c%p_0 + p_i - 1
      call mpi_isend(donor(1, trans_c%receive_lb_j(p_i), layer), 1, trans_c%lb_receive_types(p_i), p, tag, trans_c%comm, requests(1), ierr);  assert(ierr == MPI_SUCCESS)
      call mpi_irecv(target(trans_c%send_lb_i(p_i), 1, layer), 1, trans_c%lb_send_types(p_i), p, tag, trans_c%comm, requests(2), ierr);  assert(ierr == MPI_SUCCESS)
    end do

  end subroutine inverse_transpose_sdata_j_layer
  
  ! Perform the inverse operation to transpose_sdata_j.
  subroutine inverse_transpose_sdata_j_layers(target_decomp, donor, target)
    type(decomposition), target, intent(in) :: target_decomp
    real(kind = real_8), dimension(target_decomp%trans_j_decomp%ii, target_decomp%trans_j_decomp%jj, target_decomp%trans_j_decomp%layers), intent(in) :: donor
    real(kind = real_8), dimension(target_decomp%ii, target_decomp%jj, target_decomp%layers), intent(out) :: target

    integer :: layer
    
    do layer = 1, target_decomp%layers
      call inverse_transpose_sdata_j(target_decomp, donor, target, layer = layer)
    end do

  end subroutine inverse_transpose_sdata_j_layers
  
  ! Perform the inverse operation to transpose_data_j.
  subroutine inverse_transpose_data_j_layer(target_decomp, donor, target, layer)
    type(decomposition), target, intent(in) :: target_decomp
    real(kind = real_8), dimension(target_decomp%trans_j_decomp%ii + 2, target_decomp%trans_j_decomp%jj + 2, target_decomp%trans_j_decomp%layers), intent(in) :: donor
    real(kind = real_8), dimension(target_decomp%ii + 2, target_decomp%jj + 2, target_decomp%layers), intent(out) :: target
    integer, intent(in) :: layer

    integer :: ierr, p, p_i, tag
    integer, dimension(:), pointer :: requests
    type(decomposition), pointer :: donor_decomp
    type(transpose_cache), pointer :: trans_c

    trans_c => target_decomp%trans_j_c
    donor_decomp => target_decomp%trans_j_decomp
    
    do p_i = 1, trans_c%m
      requests => trans_c%requests(2 * (p_i + (layer - 1) * trans_c%m) - 1:2 * (p_i + (layer - 1) * trans_c%m))
      assert(all(requests == MPI_REQUEST_NULL))
      tag = (layer - 1) + 7 * target_decomp%layers
    
      p = trans_c%p_0 + p_i - 1
      call mpi_isend(donor(2, trans_c%receive_fb_j(p_i), layer), 1, trans_c%fb_receive_types(p_i), p, tag, trans_c%comm, requests(1), ierr);  assert(ierr == MPI_SUCCESS)
      call mpi_irecv(target(trans_c%send_fb_i(p_i), 2, layer), 1, trans_c%fb_send_types(p_i), p, tag, trans_c%comm, requests(2), ierr);  assert(ierr == MPI_SUCCESS)
    end do

  end subroutine inverse_transpose_data_j_layer
  
  ! Perform the inverse operation to transpose_data_j.
  subroutine inverse_transpose_data_j_layers(target_decomp, donor, target)
    type(decomposition), target, intent(in) :: target_decomp
    real(kind = real_8), dimension(target_decomp%trans_j_decomp%ii + 2, target_decomp%trans_j_decomp%jj + 2, target_decomp%trans_j_decomp%layers), intent(in) :: donor
    real(kind = real_8), dimension(target_decomp%ii + 2, target_decomp%jj + 2, target_decomp%layers), intent(out) :: target

    integer :: layer
    
    do layer = 1, target_decomp%layers      
      call inverse_transpose_data_j(target_decomp, donor, target, layer = layer)
    end do

  end subroutine inverse_transpose_data_j_layers
  
  ! Wait for transpose communication to complete.
  subroutine flush_ptranspose_decomp(decomp, transpose_type, layer)
    type(decomposition), intent(inout) :: decomp
    integer, intent(in) :: transpose_type
    integer, optional, intent(in) :: layer
  
    select case(transpose_type)
      case(TRANSPOSE_II)
        call flush(decomp%trans_i_c, layer = layer)
      case(TRANSPOSE_JJ)
        call flush(decomp%trans_j_c, layer = layer)
      case default
        write(6, "(a,i0)") "Transpose type: ", transpose_type
        labort("Unrecognised transpose type")
    end select
    
  end subroutine flush_ptranspose_decomp
  
  ! Wait for transpose communication to complete.
  subroutine flush_ptranspose_pdata(pdata, transpose_type, layer)
    type(parallel_data), intent(inout) :: pdata
    integer, intent(in) :: transpose_type
    integer, optional, intent(in) :: layer
  
    call flush(pdata%decomp, transpose_type, layer = layer)
    
  end subroutine flush_ptranspose_pdata
  
  ! Flush a non-blocking send / receive cache (including waiting for all send /
  ! receives to complete).
  subroutine flush_pdata_rank0(pdata)
    type(parallel_data), intent(inout) :: pdata
    
    call flush(requests = pdata%requests)
    
  end subroutine flush_pdata_rank0
  
  ! Flush a non-blocking send / receive caches (including waiting for all send /
  ! receives to complete).
  subroutine flush_pdata_rank1(pdata)
    type(parallel_data), dimension(:), intent(inout) :: pdata
    
    integer :: i
    
    do i = 1, size(pdata)
      call flush(pdata(i))
    end do
    
  end subroutine flush_pdata_rank1
  
  ! Test halo communication for a parallel data object, based upon a supplied
  ! decomposition.
  subroutine test_pdata(decomp, stat)
    type(decomposition), intent(in) :: decomp
    ! Note that this stat applies for this process *only*.
    integer, optional, intent(out) :: stat
    
    type(parallel_data) :: pdata
    
    integer :: i, j
    type(request_cache) :: rcache
    
    if(present(stat)) stat = 0
    
    call allocate(rcache, comm = decomp%cart_comm)
    call allocate(pdata, decomp, requests = rcache)
    
    if(decomp%p_jj > 1 .or. decomp%periodic_j) then
      call zero(pdata)
      do i = 1, decomp%ii
        pdata%sdata(i, 1, 1) = -real(i, kind = real_8)
        pdata%sdata(i, decomp%jj, 1) = real(i, kind = real_8)
      end do
      
      call update(pdata)
      call flush(pdata)
      
      if(decomp%p_c >= 0) then
        do i = decomp%fb_i_min, decomp%fb_i_max
          if(fne(pdata%data(i, 1, 1), real(i - decomp%fb_i_min + 1, kind = real_8))) then
            if(present(stat)) then
              stat = 1
            else
              labort("test_pdata: left receive failure")
            end if
          end if
        end do
      end if      
      if(decomp%p_a >= 0) then
        do i = decomp%fb_i_min, decomp%fb_i_max
          if(fne(pdata%data(i, decomp%jj + 2, 1), -real(i - decomp%fb_i_min + 1, kind = real_8))) then
            if(present(stat)) then
              stat = 2
            else
              labort("test_pdata: right receive failure")
            end if
          end if
        end do
      end if
      
      if(decomp%p_ii > 1) then
        call zero(pdata)
        pdata%sdata(1, 1, 1) = 1.0D0
        pdata%sdata(decomp%ii, 1, 1) = 2.0D0
        pdata%sdata(1, decomp%jj, 1) = 3.0D0
        pdata%sdata(decomp%ii, decomp%jj, 1) = 4.0D0
        
        call update(pdata)
        call flush(pdata)
        
        if(decomp%p_c >= 0) then
          if(decomp%p_cd >= 0) then
            if(fne(pdata%data(decomp%ii + 2, 1, 1), 3.0D0)) then
              if(present(stat)) then
                stat = 3
              else
                labort("test_pdata: bottom left receive failure")
              end if
            end if
          end if
          if(decomp%p_bc >= 0) then
            if(fne(pdata%data(1, 1, 1), 4.0D0)) then
              if(present(stat)) then
                stat = 4
              else
                labort("test_pdata: top left receive failure")
              end if
            end if
          end if
        else
          if(decomp%p_da >= 0) then
            if(fne(pdata%data(decomp%ii + 2, decomp%jj + 2, 1), 1.0D0)) then
              if(present(stat)) then
                stat = 5
              else
                labort("test_pdata: bottom right receive failure")
              end if
            end if
          end if
          if(decomp%p_ab >= 0) then
            if(fne(pdata%data(1, decomp%jj + 2, 1), 2.0D0)) then
              if(present(stat)) then
                stat = 6
              else
                labort("test_pdata: top right receive failure")
              end if
            end if
          end if
        end if
      end if
    end if    
    if(decomp%p_ii > 1) then
      call zero(pdata)
      do j = 1, decomp%jj
        pdata%sdata(1, j, 1) = -real(j, kind = real_8)
        pdata%sdata(decomp%ii, j, 1) = real(j, kind = real_8)
      end do
      
      call update(pdata)
      call flush(pdata)
      
      if(decomp%p_b >= 0) then
        do j = decomp%fb_j_min, decomp%fb_j_max
          if(fne(pdata%data(1, j, 1), real(j - decomp%fb_j_min + 1, kind = real_8))) then
            if(present(stat)) then
              stat = 7
            else
              labort("test_pdata: top receive failure")
            end if
          end if
        end do
      end if
      if(decomp%p_d >= 0) then
        do j = decomp%fb_j_min, decomp%fb_j_max
          if(fne(pdata%data(decomp%ii + 2, j, 1), -real(j - decomp%fb_j_min + 1, kind = real_8))) then
            if(present(stat)) then
              stat = 8
            else
              labort("test_pdata: bottom receive failure")
            end if
          end if
        end do
      end if
    end if
    
    call deallocate(pdata)
    call deallocate(rcache)
      
    if(.not. present(stat)) write(6, "(a)") "test_pdata: success"
    
  end subroutine test_pdata
  
  ! Allocate a parallel flux data object.
  subroutine allocate_fpdata(fpdata, decomp, flux_dir, requests)
    type(parallel_flux_data), intent(out) :: fpdata
    type(decomposition), target, intent(in) :: decomp
    integer, intent(in) :: flux_dir
    type(request_cache), target, optional, intent(in) :: requests
        
    fpdata%decomp => decomp
    if(present(requests)) then
      fpdata%requests => requests
    else
      fpdata%requests => rcache
    end if
    
    fpdata%flux_dir = flux_dir
    select case(flux_dir)
      case(FLUX_II)
        assert(decomp%ii > 2)
        allocate(fpdata%data(decomp%ii + 3, decomp%jj + 2, decomp%layers))
      case(FLUX_JJ)
        assert(decomp%jj > 2)
        allocate(fpdata%data(decomp%ii + 2, decomp%jj + 3, decomp%layers))
      case default
        write(6, "(a,i0)") "Flux direction: ", flux_dir
        labort("Unrecognised flux direction")
    end select
  
  end subroutine allocate_fpdata
  
  ! Deallocate a parallel flux data object.
  subroutine deallocate_fpdata(fpdata)
    type(parallel_flux_data), intent(inout) :: fpdata
    
    deallocate(fpdata%data)
    
  end subroutine deallocate_fpdata
  
  ! Zero a parallel flux data object.
  subroutine zero_fpdata(fpdata)
    type(parallel_flux_data), intent(inout) :: fpdata
    
    fpdata%data = 0.0D0
  
  end subroutine zero_fpdata
  
  ! Set data in one parallel flux data object.
  subroutine set_fpdata_real(target, donor, layer)
    type(parallel_flux_data), intent(inout) :: target
    real(kind = real_8), intent(in) :: donor
    integer, optional, intent(in) :: layer
    
    if(present(layer)) then
      target%data(:, :, layer) = donor
    else
      target%data = donor
    end if
    
  end subroutine set_fpdata_real
  
  ! Update a parallel flux data object via communication. Note that this only
  ! sets up the non-blocking send / receives - the relevant request cache must
  ! be flushed for the data update to complete.
  subroutine update_fpdata(fpdata, level_1)
    type(parallel_flux_data), intent(inout) :: fpdata
    logical, optional, intent(in) :: level_1
    
    select case(fpdata%flux_dir)
      case(FLUX_II)
        call update_flux_i(fpdata%data, fpdata%decomp, requests = fpdata%requests, level_1 = level_1)
      case(FLUX_JJ)
        call update_flux_j(fpdata%data, fpdata%decomp, requests = fpdata%requests, level_1 = level_1)
      case default
        write(6, "(a,i0)") "Flux direction: ", fpdata%flux_dir
        labort("Unrecognised flux direction")       
    end select
    
  end subroutine update_fpdata
  
  ! Update a cell face flux variable, stored at the mid-point of two nodes
  ! separated in the i-direction (meridional). Note that this only sets up the
  ! non-blocking send / receives - the relevant request cache must be flushed
  ! for the data update to complete.
  !          |
  !  --- zeta_i,j   ---
  !          |
  !      flux_i+1,j
  !          |
  !  --- zeta_i+1,j ---
  subroutine update_flux_i(flux_i, decomp, requests, level_1)
    type(decomposition), intent(in) :: decomp
    real(kind = real_8), dimension(decomp%ii + 3, decomp%jj + 2, decomp%layers), intent(inout) :: flux_i
    type(request_cache), target, optional, intent(inout) :: requests
    logical, optional, intent(in) :: level_1
    
    integer :: ierr, layer, tag_base
    type(request_cache), pointer :: lrequests
    
    if(decomp%p_d == -1 .and. decomp%p_b == -1) return
    
    if(present(requests)) then
      lrequests => requests
    else
      lrequests => rcache
    end if  
    
    do layer = 1, decomp%layers
      tag_base = lrequests%tag_base

      if(decomp%p_d >= 0) then
        ! Send / receive with logical below
        call mpi_isend(flux_i(decomp%fa_i_max - 1, decomp%fa_j_min + 1, layer), 1, decomp%ts_jj2_iip3, decomp%p_d, tag_base, lrequests%comm, next_request(lrequests), ierr)
          assert(ierr == MPI_SUCCESS)
        call mpi_irecv(flux_i(decomp%fa_i_max + 1, decomp%fa_j_min + 1, layer), 1, decomp%ts_jj2_iip3, decomp%p_d, tag_base + 1, lrequests%comm, next_request(lrequests), ierr)
          assert(ierr == MPI_SUCCESS)
      end if
      if(decomp%p_b >= 0) then
        ! Send / receive with logical above
        call mpi_isend(flux_i(decomp%fa_i_min + 2, decomp%fa_j_min + 1, layer), 1, decomp%ts_jj2_iip3, decomp%p_b, tag_base + 1, lrequests%comm, next_request(lrequests), ierr)
          assert(ierr == MPI_SUCCESS)
        call mpi_irecv(flux_i(decomp%fa_i_min, decomp%fa_j_min + 1, layer), 1, decomp%ts_jj2_iip3, decomp%p_b, tag_base, lrequests%comm, next_request(lrequests), ierr)
          assert(ierr == MPI_SUCCESS)        
      end if

      call increment_tag(lrequests, increment = 2)
    end do
    
    if(present_and_true(level_1)) then
      do layer = 1, decomp%layers
        tag_base = lrequests%tag_base

        if(decomp%p_d >= 0) then
          ! Send / receive with logical below
          call mpi_irecv(flux_i(decomp%fa_i_max, decomp%fa_j_min + 1, layer), 1, decomp%ts_jj2_iip3, decomp%p_d, tag_base + 1, lrequests%comm, next_request(lrequests), ierr)
            assert(ierr == MPI_SUCCESS)
        end if
        if(decomp%p_b >= 0) then
          ! Send / receive with logical above
          call mpi_isend(flux_i(decomp%fa_i_min + 1, decomp%fa_j_min + 1, layer), 1, decomp%ts_jj2_iip3, decomp%p_b, tag_base + 1, lrequests%comm, next_request(lrequests), ierr)
            assert(ierr == MPI_SUCCESS)   
        end if

        call increment_tag(lrequests)
      end do
    end if
    
  end subroutine update_flux_i
  
  ! Update a cell face flux variable, stored at the mid-point of two nodes
  ! separated in the j-direction (zonal). Note that this only sets up the
  ! non-blocking send / receives - the relevant request cache must be flushed
  ! for the data update to complete.
  !        |                           |
  !    zeta_i,j --- flux_i,j+1 --- zeta_i,j+1
  !        |                           |
  subroutine update_flux_j(flux_j, decomp, requests, level_1)
    type(decomposition), intent(in) :: decomp
    real(kind = real_8), dimension(decomp%ii + 2, decomp%jj + 3, decomp%layers), intent(inout) :: flux_j
    type(request_cache), target, optional, intent(inout) :: requests
    logical, optional, intent(in) :: level_1
    
    integer :: ierr, layer, tag_base
    type(request_cache), pointer :: lrequests
    
    if(decomp%p_c == -1 .and. decomp%p_a == -1) return
   
    if(present(requests)) then
      lrequests => requests
    else
      lrequests => rcache
    end if   
    
    do layer = 1, decomp%layers
      tag_base = lrequests%tag_base
            
      if(decomp%p_c >= 0) then
        ! Send / receive with logical left
        call mpi_isend(flux_j(decomp%fa_i_min + 1, decomp%fa_j_min + 2, layer), decomp%ii2, MPI_DOUBLE_PRECISION, decomp%p_c, tag_base, lrequests%comm, next_request(lrequests), ierr)
          assert(ierr == MPI_SUCCESS)
        call mpi_irecv(flux_j(decomp%fa_i_min + 1, decomp%fa_j_min, layer), decomp%ii2, MPI_DOUBLE_PRECISION, decomp%p_c, tag_base + 1, lrequests%comm, next_request(lrequests), ierr)
          assert(ierr == MPI_SUCCESS)
      end if
      if(decomp%p_a >= 0) then
        ! Send / receive with logical right
        call mpi_isend(flux_j(decomp%fa_i_min + 1, decomp%fa_j_max - 1, layer), decomp%ii2, MPI_DOUBLE_PRECISION, decomp%p_a, tag_base + 1, lrequests%comm, next_request(lrequests), ierr)
          assert(ierr == MPI_SUCCESS)
        call mpi_irecv(flux_j(decomp%fa_i_min + 1, decomp%fa_j_max + 1, layer), decomp%ii2, MPI_DOUBLE_PRECISION, decomp%p_a, tag_base, lrequests%comm, next_request(lrequests), ierr)
          assert(ierr == MPI_SUCCESS)
      end if

      call increment_tag(lrequests, increment = 2)
    end do
    
    if(present_and_true(level_1)) then
      do layer = 1, decomp%layers
        tag_base = lrequests%tag_base
              
        if(decomp%p_c >= 0) then
          ! Send / receive with logical left
          call mpi_isend(flux_j(decomp%fa_i_min + 1, decomp%fa_j_min + 1, layer), decomp%ii2, MPI_DOUBLE_PRECISION, decomp%p_c, tag_base, lrequests%comm, next_request(lrequests), ierr)
            assert(ierr == MPI_SUCCESS)
        end if
        if(decomp%p_a >= 0) then
          ! Send / receive with logical right
          call mpi_irecv(flux_j(decomp%fa_i_min + 1, decomp%fa_j_max, layer), decomp%ii2, MPI_DOUBLE_PRECISION, decomp%p_a, tag_base, lrequests%comm, next_request(lrequests), ierr)
            assert(ierr == MPI_SUCCESS)
        end if

        call increment_tag(lrequests)
      end do
    end if
    
  end subroutine update_flux_j
  
  ! Flush a non-blocking send / receive cache (including waiting for all send /
  ! receives to complete).
  subroutine flush_pfdata(pfdata)
    type(parallel_flux_data), intent(inout) :: pfdata
    
    call flush(requests = pfdata%requests)
    
  end subroutine flush_pfdata
  
end module qg_parallel
