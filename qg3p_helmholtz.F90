#include "config.h"

! Helmholtz solver routines for the 3D QG model. The FFT code here is based on
! FFTPACK.
module qg_helmholtz

  use lapack_interfaces
  use qg_parallel
  use qg_parameters, only : hsqt2, INFINITY, pi, real_8, sqrt3, two_pi
  
  implicit none
  
  include "mpif.h"
  
  private
  
  public :: helm
  public :: helm_cache
  public :: allocate, deallocate
  
  ! Maximum number of factors in the rfft factorisation.
  integer, parameter :: MAX_FACTORS = 1024
  ! Cache of rfft data.
  type rfft_cache
    ! The size of the data
    integer :: n
    ! Working memory
    real(kind = real_8), dimension(:), pointer :: wa
    ! Number of factors of n
    integer :: nf
    ! The factors of n
    integer, dimension(MAX_FACTORS) :: ifac
  end type rfft_cache
  
  ! Cache of sint data.
  type sint_cache
    ! The size of the data
    integer :: n
    ! Working memory
    real(kind = real_8), dimension(:), pointer :: cc
    real(kind = real_8), dimension(:), pointer :: was
    ! Cache of rfft data.
    type(rfft_cache) :: rfft_c
  end type
  
  ! Cache of helm data.
  type helm_cache
    ! Decomposition
    type(decomposition), pointer :: decomp
    ! Cache of sint data
    type(sint_cache) :: sint_c
    ! Helmholtz operator coefficients
    real(kind = real_8), dimension(:), pointer :: ss
#ifdef PH1D
    ! Parallel 1D Helmholtz solver cache
    type(h1d_cache), dimension(:), pointer :: h1d_c
#else
    ! Coefficients used to perform Gaussian elimination of the tri-diagonal
    ! (non-periodic) 1D Helmholtz matrices
    real(kind = real_8), dimension(:, :, :), pointer :: eig
    ! Coefficients used to perform supplementary Gaussian elimination operations
    ! for the periodic 1D Helmholtz matrices
    real(kind = real_8), dimension(:, :, :), pointer :: periodic_b, periodic_c
#endif
    ! Cache used to store the result of a transpose to a decomposition
    ! contiguous in i
    real(kind = real_8), dimension(:, :, :), pointer :: zeta_i
    ! Cache used to store the result of a transpose to a decomposition
    ! contiguous in j
    real(kind = real_8), dimension(:, :, :), pointer :: zeta_j
  end type helm_cache
    
  interface allocate
    module procedure allocate_h1d_c, allocate_rfft_c, allocate_sint_c, &
      & allocate_helm_c
  end interface allocate
  
  interface deallocate
    module procedure deallocate_h1d_c, deallocate_rfft_c, deallocate_sint_c, &
      & deallocate_helm_c
  end interface deallocate
  
  ! Cache for parallel 1D Helmholtz solves.
  type h1d_cache
    ! The number of 1D problems owned by this process
    integer :: ii2
    ! Maximum j index
    integer :: max_j
    ! Decomposition contiguous in i
    type(decomposition), pointer :: decomp_i
    ! Decomposition contiguous in j
    type(decomposition), pointer :: decomp_j
    ! Coefficients used to perform the local Gaussian elimination
    real(kind = real_8), dimension(:, :), pointer :: ge_fact
    ! Local boundary contributions
    real(kind = real_8), dimension(:, :), pointer :: psi1, psi2
    ! Coefficients used to perform the global Gaussian elimination
    real(kind = real_8), dimension(:, :, :), pointer :: pge_fact
    ! The RHS vectors for the boundary coefficients 
    real(kind = real_8), dimension(:, :), pointer :: rhs
    ! The solution for the boundary coefficients
    real(kind = real_8), dimension(:, :), pointer :: coeff
    ! Cache for received data
    real(kind = real_8), dimension(:, :, :), pointer :: recv
    ! mpi_alltoallv arguments
    integer :: comm
    integer, dimension(:), pointer :: rcounts1, rcounts2, rindices, scounts1, &
      & scounts2, sindices
#ifdef PH1D_BUFFER
    ! Buffers for mpi_alltoall
    integer :: npsi_buffer
    real(kind = real_8), dimension(:, :, :, :), pointer :: psi_buffer
    integer :: nzeta_buffer
    real(kind = real_8), dimension(:, :, :), pointer :: zeta_buffer
    integer :: ncoeff_buffer
    real(kind = real_8), dimension(:, :, :, :), pointer :: coeff_buffer
#endif     
  end type h1d_cache
  
  ! Time spent in the mpi_alltoallv when performing the parallel 1D Helmholtz
  ! solves.
  real(kind = real_8), save, public :: h1d_time = 0.0D0
  ! Used to measure the start time of communication.
  real(kind = real_8), save :: wtime0
  
contains

  ! Allocate an rfft cache object.
  subroutine allocate_rfft_c(rfft_c, n)
    type(rfft_cache), intent(out) :: rfft_c
    integer, intent(in) :: n
    
    rfft_c%n = n
    allocate(rfft_c%wa(n))
  
  end subroutine allocate_rfft_c
    
  ! Allocate a sint cache object.
  subroutine allocate_sint_c(sint_c, n)
    type(sint_cache), intent(out) :: sint_c
    integer, intent(in) :: n
    
    sint_c%n = n
    allocate(sint_c%cc(n + 1))
    allocate(sint_c%was(n / 2))
    
  end subroutine allocate_sint_c
  
  ! Deallocate an rfft cache object.
  subroutine deallocate_rfft_c(rfft_c)
    type(rfft_cache), intent(inout) :: rfft_c
    
    deallocate(rfft_c%wa)
    
  end subroutine deallocate_rfft_c
  
  ! Deallocate a sint cache object.
  subroutine deallocate_sint_c(sint_c)
    type(sint_cache), intent(inout) :: sint_c
    
    deallocate(sint_c%cc)
    deallocate(sint_c%was)
    call deallocate(sint_c%rfft_c)
    
  end subroutine deallocate_sint_c

  ! Initialise a sint cache object.
  subroutine sinti(sint_c, n)
    type(sint_cache), intent(out) :: sint_c
    integer, intent(in) :: n
    
    integer :: i, np1, ns2
    real(kind = real_8) :: dt
    
    call allocate(sint_c, n)
    
    ns2 = n / 2
    np1 = n + 1
    
    dt = pi / real(np1, kind = real_8)
    do i = 1, ns2
       sint_c%was(i) = 2.0D0 * dsin(real(i, kind = real_8) * dt)
    end do
    call rffti(sint_c%rfft_c, np1)
      
  end subroutine sinti

  ! Initialise an rfft cache object.
  subroutine rffti(rfft_c, n)
    type(rfft_cache), intent(out) :: rfft_c
    integer, intent(in) :: n

    call allocate(rfft_c, n)
    
    if(n == 0) return

    call rffti1(rfft_c, n)
    
  end subroutine rffti

  ! Low-level rfft initialisation routine.
  subroutine rffti1(rfft_c, n)
    type(rfft_cache), intent(inout) :: rfft_c
    integer, intent(in) :: n

    integer :: i, ido, ii, ip, is, j, k1, l1, l2, ld, nf, nl, ntry
    integer, parameter, dimension(9) :: ntryh = (/4, 2, 3, 5, 7, 11, 13, 17, 19/)
    ! This disables the factor four transform routine
    !integer, parameter, dimension(8) :: ntryh = (/2, 3, 5, 7, 11, 13, 17, 19/)
    real(kind = real_8) :: arg, argh, argld, fi
  
    ! The following seeks factors ntry for nl. It first searched through some
    ! specific factors (in order), 4, 2, 3, 5, before searching through all
    ! remaining odd numbers in increasing order. This is not a prime
    ! factorisation, as factors of 4 are searched for before factors of 2.
    nl = n
    nf = 0
    j = 1
    ntry = ntryh(1)
    do while(nl /= 1)
      do while(mod(nl, ntry) /= 0)
        j = j + 1
        if(j <= size(ntryh)) then
          ntry = ntryh(j)
        else
          ntry = ntry + 2
        end if
      end do
      ! We've found a factor
      nf = nf + 1
      if(nf > MAX_FACTORS) then
        labort("Insufficient ifac matrix size (too many factors)")
      end if
      ! Store the factor
      rfft_c%ifac(nf) = ntry
      ! This is what's left
      nl = nl / ntry

      ! Move all the factors of two to the head of the factors list. Hence the
      ! list will contain, in order: factors of 2, 4, 3, 5, 7, 9, etc.
      if(ntry == 2 .and. nf > 1) then
        do i = nf, 2, -1
          rfft_c%ifac(i) = rfft_c%ifac(i - 1)
        end do  
        rfft_c%ifac(1) = 2
      end if

    end do
    ! Store some additional meta-data
    ! The number of factors
    rfft_c%nf = nf
    
    if(rfft_c%nf == 0) then
      ! We have no factors
      return
    end if
    
    argh = two_pi / real(n, kind = real_8)
    is = 0
    l1 = 1
    do k1 = 1, rfft_c%nf - 1
      ip = rfft_c%ifac(k1)
      
      ! l1 = product of factors, excluding the current processes factor
      ! l2 = product of factors, including the currently processed factor
      l2 = l1 * ip
      ! ido, l2 = factor pair
      ido = n / l2
      
      ld = 0
      do j = 1, ip - 1
        ld = ld + l1  ! = j * l1
        argld = real(ld, kind = real_8) * argh
        fi = 0.0D0
        do ii = 2, ido - 1, 2
           i = is + ii
           fi = fi + 1.0D0  ! = ii / 2.0
           
           arg = fi * argld
           rfft_c%wa(i - 1) = dcos(arg)
           rfft_c%wa(i) = dsin(arg)
        end do
        is = is + ido
     end do
     
     l1 = l2
    end do
  
  end subroutine rffti1
  
  ! Allocate a helm cache.
  subroutine allocate_helm_c(helm_c, decomp, ss)
    type(helm_cache), intent(out) :: helm_c
    type(decomposition), target, intent(in) :: decomp
    real(kind = real_8), dimension(decomp%layers), intent(in) :: ss
        
    logical :: trans_i, trans_j
    type(decomposition), pointer :: decomp_i, decomp_j
#ifdef PH1D
    integer :: layer
#endif
    
    helm_c%decomp => decomp
    
    call sinti(helm_c%sint_c, decomp%g_ii2)
    
    allocate(helm_c%ss(decomp%layers))
    helm_c%ss = ss
    
#ifdef INTERMEDIATE_TRANSPOSES
    trans_i = decomp%p_ii > 1
    trans_j = decomp%p_jj > 1
    if(trans_i) then
      decomp_i => decomp%trans_i_decomp
    else
      decomp_i => decomp
    end if   
    if(trans_j) then
      decomp_j => decomp%trans_j_decomp
    else
      decomp_j => decomp
    end if
#else
    trans_i = decomp%p_ii > 1
    trans_j = trans_i .or. decomp%p_jj > 1
    if(trans_i) then
      decomp_i => decomp%trans_i_decomp
    else
      decomp_i => decomp
    end if   
    if(trans_j) then
      decomp_j => decomp_i%trans_j_decomp
    else
      decomp_j => decomp_i
    end if
#endif
    
#ifdef PH1D
    allocate(helm_c%h1d_c(decomp_i%layers))
    do layer = 1, decomp_i%layers
      call allocate(helm_c%h1d_c(layer), decomp_i, decomp_j, ss(layer))
    end do
#else
    allocate(helm_c%eig(decomp_j%ii, decomp_j%jj, decomp_j%layers))
    call poinit_j(decomp_j, ss, helm_c%eig)
    
    if(decomp%periodic_j) then
      allocate(helm_c%periodic_b(decomp_j%ii, decomp_j%jj, decomp_j%layers))
      allocate(helm_c%periodic_c(decomp_j%ii, decomp_j%jj, decomp_j%layers))
      call periodic_coefficients(decomp_j, ss, helm_c%eig, helm_c%periodic_b, helm_c%periodic_c)
    else
      nullify(helm_c%periodic_b)
      nullify(helm_c%periodic_c)
    end if
#endif

#ifdef INTERMEDIATE_TRANSPOSES
    if(trans_i) then
      allocate(helm_c%zeta_i(decomp_i%ii, decomp_i%jj, decomp_i%layers))
#ifdef DDEBUG
      helm_c%zeta_i = INFINITY
#endif
    end if
    if(trans_j) then
      allocate(helm_c%zeta_j(decomp_j%ii, decomp_j%jj, decomp_j%layers))
#ifdef DDEBUG
      helm_c%zeta_j = INFINITY
#endif
    end if
#else
    if(trans_i) then
      allocate(helm_c%zeta_i(decomp_i%ii, decomp_i%jj, decomp_i%layers))
      allocate(helm_c%zeta_j(decomp_j%ii, decomp_j%jj, decomp_j%layers))
#ifdef DDEBUG
      helm_c%zeta_i = INFINITY
      helm_c%zeta_j = INFINITY
#endif
    else if(trans_j) then
      allocate(helm_c%zeta_j(decomp_j%ii, decomp_j%jj, decomp_j%layers))
#ifdef DDEBUG
      helm_c%zeta_j = INFINITY
#endif
    end if
#endif
    
  end subroutine allocate_helm_c
  
  ! Deallocate a helm cache.
  subroutine deallocate_helm_c(helm_c)
    type(helm_cache), intent(inout) :: helm_c
    
#ifdef PH1D
    integer :: layer
#endif
    
    call deallocate(helm_c%sint_c)
    deallocate(helm_c%ss)
#ifdef PH1D
    do layer = 1, size(helm_c%h1d_c)
      call deallocate(helm_c%h1d_c(layer))
    end do
    deallocate(helm_c%h1d_c)
#else
    deallocate(helm_c%eig)
    if(helm_c%decomp%periodic_j) then
      deallocate(helm_c%periodic_b)
      deallocate(helm_c%periodic_c)
    end if
#endif
#ifdef INTERMEDIATE_TRANSPOSES
    if(helm_c%decomp%p_ii > 1) then
      deallocate(helm_c%zeta_i)
    end if
    if(helm_c%decomp%p_jj > 1) then
      deallocate(helm_c%zeta_j)
    end if
#else
    if(helm_c%decomp%p_ii > 1) then
      deallocate(helm_c%zeta_i)
      deallocate(helm_c%zeta_j)
    else if(helm_c%decomp%p_jj > 1) then
      deallocate(helm_c%zeta_j)
    end if
#endif
    
  end subroutine deallocate_helm_c

  ! Initialise coefficients used to invert the (non-periodic) Helmholtz matrices.
  subroutine poinit_j(decomp, ss, eig)
    type(decomposition), intent(in) :: decomp
    real(kind = real_8), dimension(decomp%layers), intent(in) :: ss
    real(kind = real_8), dimension(:, :, :), intent(out) :: eig

    integer :: i, j, layer
    real(kind = real_8) :: cff  
    
    assert(decomp%p_jj == 1)
    assert(size(eig, 1) == decomp%ii)
    assert(size(eig, 2) == decomp%jj)
    assert(size(eig, 3) == decomp%layers)
    
#ifdef DDEBUG
    eig = INFINITY
#endif
    
    cff = pi / (2.0D0 * real(decomp%g_ii_int, kind = real_8))
    do layer = 1, decomp%layers
      do i = decomp%l_i_min, decomp%l_i_max
        eig(i, decomp%l_j_min:decomp%l_j_max, layer) = -2.0D0 - ss(layer) - 4.0D0 * dsin(real(decomp%ii_0 + i - 1, kind = real_8) * cff) ** 2
      end do
    end do
    
    do layer = 1, decomp%layers        
      do j = 3, decomp%jj - 1
        do i = decomp%l_i_min, decomp%l_i_max
          eig(i, j, layer) = eig(i, j, layer) - 1.0D0 / eig(i, j - 1, layer)
        end do
      end do
    end do
        
    do layer = 1, decomp%layers
      do j = 2, decomp%jj - 1
        do i = decomp%l_i_min, decomp%l_i_max
          eig(i, j, layer) = 1.0D0 / eig(i, j, layer)
        end do
      end do
    end do
    
  end subroutine poinit_j
  
  ! Compute coefficients required to perform supplementary Gaussian elimination
  ! operations for inversion of the periodic Helmholtz matrices.
  subroutine periodic_coefficients(decomp, ss, eig, b, c)
    type(decomposition), intent(in) :: decomp
    real(kind = real_8), dimension(decomp%layers), intent(in) :: ss
    real(kind = real_8), dimension(:, :, :), intent(in) :: eig
    real(kind = real_8), dimension(:, :, :), intent(out) :: b
    real(kind = real_8), dimension(:, :, :), intent(out) :: c
    
    integer :: i, j, layer
    real(kind = real_8) :: a, norm_fact
    
    assert(decomp%periodic_j)
    assert(decomp%p_jj == 1)
    assert(decomp%jj >= 3)
    assert(size(eig, 1) == decomp%ii)
    assert(size(eig, 2) == decomp%jj)
    assert(size(eig, 3) == decomp%layers)
    assert(size(b, 1) == decomp%ii)
    assert(size(b, 2) == decomp%jj)
    assert(size(b, 3) == decomp%layers)
    assert(size(c, 1) == decomp%ii)
    assert(size(c, 2) == decomp%jj)
    assert(size(c, 3) == decomp%layers)
    
    ! Assemble the first and last column of the Helmholtz matrices
    norm_fact = pi / (2.0D0 * real(decomp%g_ii_int, kind = real_8))
    do layer = 1, decomp%layers
      do i = decomp%l_i_min, decomp%l_i_max
        ! Diagonal entries of the Helmholtz matrix
        a = -2.0D0 - ss(layer) - 4.0D0 * dsin(real(decomp%ii_0 + i - 1, kind = real_8) * norm_fact) ** 2
      
        ! Entries in the first column of the Helmholtz matrix
        b(i, 1, layer) = a
        b(i, 2, layer) = 1.0D0
        b(i, 3:decomp%jj - 1, layer) = 0.0D0
        b(i, decomp%jj, layer) = 1.0D0
        
        ! Entries in the last column of the Helmholtz matrix
        c(i, 1, layer) = 1.0D0
        c(i, 2:decomp%jj - 2, layer) = 0.0D0
        c(i, decomp%jj - 1, layer) = 1.0D0
        c(i, decomp%jj, layer) = a
      end do
    
      ! Apply Gaussian elimination operations for the non-periodic system to the
      ! first and last columns of the Helmholtz matrix
      do j = 3, decomp%jj - 1   
        do i = decomp%l_i_min, decomp%l_i_max
          b(i, j, layer) = b(i, j, layer) - b(i, j - 1, layer) * eig(i, j - 1, layer)
          c(i, j, layer) = c(i, j, layer) - c(i, j - 1, layer) * eig(i, j - 1, layer)
        end do
      end do
      do i = decomp%l_i_min, decomp%l_i_max
        b(i, decomp%jj - 1, layer) = b(i, decomp%jj - 1, layer) * eig(i, decomp%jj - 1, layer)
        c(i, decomp%jj - 1, layer) = c(i, decomp%jj - 1, layer) * eig(i, decomp%jj - 1, layer)
      end do
      do j = decomp%jj - 2, 2, -1
        do i = decomp%l_i_min, decomp%l_i_max 
          b(i, j, layer) = (b(i, j, layer) - b(i, j + 1, layer)) * eig(i, j, layer)
          c(i, j, layer) = (c(i, j, layer) - c(i, j + 1, layer)) * eig(i, j, layer)
        end do
      end do

      ! Apply Gaussian elimination operations for the first and last rows to the
      ! first and last columns of the Helmholtz matrix
      do i = decomp%l_i_min, decomp%l_i_max
        b(i, 1, layer) = b(i, 1, layer) - b(i, 2, layer)
        c(i, 1, layer) = c(i, 1, layer) - c(i, 2, layer)
      end do      
      do i = decomp%l_i_min, decomp%l_i_max
        b(i, decomp%jj, layer) = b(i, decomp%jj, layer) - b(i, decomp%jj - 1, layer)
        c(i, decomp%jj, layer) = c(i, decomp%jj, layer) - c(i, decomp%jj - 1, layer)
      end do

      ! Apply Gaussian elimination operations for the first column to the last
      ! column of the Helmholtz matrix
      do i = decomp%l_i_min, decomp%l_i_max
        c(i, 1, layer) = c(i, 1, layer) / b(i, 1, layer)
      end do

      do j = 2, decomp%jj
        do i = decomp%l_i_min, decomp%l_i_max
          c(i, j, layer) = c(i, j, layer) - c(i, 1, layer) * b(i, j, layer)
        end do
      end do
    end do    
    
  end subroutine periodic_coefficients
  
  ! Perform supplementary Gaussian elimination operations required to
  ! invert the periodic Helmholtz matrices. This should be applied *after* the
  ! Gaussian elimination operations required to invert the non-periodic
  ! system.
  subroutine periodic_ge(decomp, b, c, zeta)
    type(decomposition), intent(in) :: decomp
    real(kind = real_8), dimension(:, :), intent(in) :: b
    real(kind = real_8), dimension(:, :), intent(in) :: c
    real(kind = real_8), dimension(:, :), intent(inout) :: zeta
  
    integer :: i, j
  
    assert(decomp%periodic_j)
    assert(decomp%p_jj == 1)
    assert(decomp%jj >= 2)
    assert(size(b, 1) == decomp%ii)
    assert(size(b, 2) == decomp%jj)
    assert(size(c, 1) == decomp%ii)
    assert(size(c, 2) == decomp%jj)
    assert(size(zeta, 1) == decomp%ii)
    assert(size(zeta, 2) == decomp%jj)
  
    ! Gaussian elimination for the first and last rows
    do i = decomp%l_i_min, decomp%l_i_max
      zeta(i, 1) = zeta(i, 1) - zeta(i, 2)
    end do    
    do i = decomp%l_i_min, decomp%l_i_max
      zeta(i, decomp%jj) = zeta(i, decomp%jj) - zeta(i, decomp%jj - 1)
    end do
    
    ! Gaussian elimination for the first and last columns
    do i = decomp%l_i_min, decomp%l_i_max
      zeta(i, 1) = zeta(i, 1) / b(i, 1)
    end do
    do j = 2, decomp%jj
      do i = decomp%l_i_min, decomp%l_i_max
        zeta(i, j) = zeta(i, j) - zeta(i, 1) * b(i, j)
      end do
    end do
    
    do i = decomp%l_i_min, decomp%l_i_max
      zeta(i, decomp%jj) = zeta(i, decomp%jj) / c(i, decomp%jj)
    end do
    do j = decomp%jj - 1, 1, -1
      do i = decomp%l_i_min, decomp%l_i_max
        zeta(i, j) = zeta(i, j) - zeta(i, decomp%jj) * c(i, j)
      end do
    end do
    
  end subroutine periodic_ge
  
  ! Parallel Helmholtz solver. On input zeta is the RHS, on output zeta is the
  ! solution vector.
  subroutine helm(helm_c, zeta)
    type(helm_cache), intent(inout) :: helm_c
    real(kind = real_8), dimension(:, :, :), target, intent(inout) :: zeta

    integer :: i, j, layer
    real(kind = real_8) :: cff
    
    type(decomposition), pointer :: decomp
    logical :: trans_i, trans_j
    real(kind = real_8), dimension(:, :, :), pointer :: zeta_i, zeta_j
    type(decomposition), pointer :: decomp_i, decomp_j

    assert(size(zeta, 1) == helm_c%decomp%ii)
    assert(size(zeta, 2) == helm_c%decomp%jj)
    assert(size(zeta, 3) == helm_c%decomp%layers)
    
    decomp => helm_c%decomp
    
#ifdef INTERMEDIATE_TRANSPOSES
    trans_i = decomp%p_ii > 1
    trans_j = decomp%p_jj > 1
    if(trans_i) then
      zeta_i => helm_c%zeta_i
      decomp_i => decomp%trans_i_decomp
    else
      zeta_i => zeta
      decomp_i => decomp
    end if   
    if(trans_j) then
      zeta_j => helm_c%zeta_j
      decomp_j => decomp%trans_j_decomp
    else
      zeta_j => zeta
      decomp_j => decomp
    end if
#else
    trans_i = decomp%p_ii > 1
    trans_j = trans_i .or. decomp%p_jj > 1
    if(trans_i) then
      zeta_i => helm_c%zeta_i
      decomp_i => decomp%trans_i_decomp
    else
      zeta_i => zeta
      decomp_i => decomp
    end if   
    if(trans_j) then
      zeta_j => helm_c%zeta_j
      decomp_j => decomp_i%trans_j_decomp
    else
      zeta_j => zeta_i
      decomp_j => decomp_i
    end if
#endif

    ! Transpose to a decomposition contiguous in i
    if(trans_i) call ptranspose(decomp, zeta, helm_c%zeta_i, TRANSPOSE_II)
    
    ! Forward transform
    cff = 0.5D0 / real(decomp_i%g_ii_int, kind = real_8)
    do layer = 1, decomp_i%layers
      if(trans_i) call flush(decomp, TRANSPOSE_II, layer = layer)
    
      do j = decomp_i%l_j_min, decomp_i%l_j_max
        do i = decomp_i%l_i_min, decomp_i%l_i_max 
          zeta_i(i, j, layer) = cff * zeta_i(i, j, layer)
        end do
        call sint(helm_c%sint_c, zeta_i(decomp_i%l_i_min:decomp_i%l_i_max + 1, j, layer))
      end do
      
#ifdef PH1D
#ifdef INTERMEDIATE_TRANSPOSES
      labort("Intermediate transposes not supported with parallel 1D Helmholtz solves")
#endif

      call h1d(helm_c%h1d_c(layer), zeta_i(:, :, layer))
#else
#ifdef INTERMEDIATE_TRANSPOSES
      ! Transpose back to the original decomposition
      if(trans_i) call pinverse_transpose(decomp, zeta_i, zeta, TRANSPOSE_II, layer = layer)
    end do 
    
    do layer = 1, decomp%layers
      if(trans_i) call flush(decomp, TRANSPOSE_II, layer = layer)
      
      ! Transpose to a decomposition contiguous in j
      if(trans_j) call ptranspose(decomp, zeta, zeta_j, TRANSPOSE_JJ, layer = layer)
    end do
    
    do layer = 1, decomp_j%layers
      if(trans_j) call flush(decomp, TRANSPOSE_JJ, layer = layer)
#else
      ! Transpose to a decomposition contiguous in j
      if(trans_j) call ptranspose(decomp_i, zeta_i, zeta_j, TRANSPOSE_JJ, layer = layer)
    end do
    
    do layer = 1, decomp_j%layers
      if(trans_j) call flush(decomp_i, TRANSPOSE_JJ, layer = layer)
#endif
    
      ! Inversion of the non-periodic 1D Helmholtz matrices:
      ! Lower triangular Gaussian elimination
      do j = 3, decomp_j%jj - 1  
        do i = decomp_j%l_i_min, decomp_j%l_i_max
          zeta_j(i, j, layer) = zeta_j(i, j, layer) - zeta_j(i, j - 1, layer) * helm_c%eig(i, j - 1, layer)
        end do
      end do
      ! Upper triangular Gaussian elimination, combined with normalisation
      do i = decomp_j%l_i_min, decomp_j%l_i_max
        zeta_j(i, decomp_j%jj - 1, layer) = zeta_j(i, decomp_j%jj - 1, layer) * helm_c%eig(i, decomp_j%jj - 1, layer)
      end do
      do j = decomp_j%jj - 2, 2, -1
        do i = decomp_j%l_i_min, decomp_j%l_i_max 
          zeta_j(i, j, layer) = (zeta_j(i, j, layer) - zeta_j(i, j + 1, layer)) * helm_c%eig(i, j, layer)
        end do
      end do
      
      ! Additional operations required for inversion of the periodic
      ! 1D Helmholtz matrices
      if(decomp%periodic_j) call periodic_ge(decomp_j, helm_c%periodic_b(:, :, layer), helm_c%periodic_c(:, :, layer), zeta_j(:, :, layer))
      
#ifdef INTERMEDIATE_TRANSPOSES
      ! Transpose back to the original decomposition
      if(trans_j) call pinverse_transpose(decomp, zeta_j, zeta, TRANSPOSE_JJ, layer = layer)
    end do
    
    do layer = 1, decomp%layers
      if(trans_j) call flush(decomp, TRANSPOSE_JJ, layer = layer)
      
      ! Transpose to a decomposition contiguous in i
      if(trans_i) call ptranspose(decomp, zeta, zeta_i, TRANSPOSE_II, layer = layer)
    end do
    
    do layer = 1, decomp_i%layers
      if(trans_i) call flush(decomp, TRANSPOSE_II, layer = layer)
#else
      ! Transpose back to the decomposition contiguous in j
      if(trans_j) call pinverse_transpose(decomp_i, zeta_j, zeta_i, TRANSPOSE_JJ, layer = layer)
    end do
    
    do layer = 1, decomp_i%layers
      if(trans_j) call flush(decomp_i, TRANSPOSE_JJ, layer = layer)
#endif
#endif
      ! Inverse transform
      do j = decomp_i%l_j_min, decomp_i%l_j_max
        call sint(helm_c%sint_c, zeta_i(decomp_i%l_i_min:decomp_i%l_i_max + 1, j, layer))
        zeta_i(decomp_i%ii, j, layer) = 0.0D0
      end do
      
      ! Transpose back to the original decomposition
      if(trans_i) call pinverse_transpose(decomp, zeta_i, zeta, TRANSPOSE_II, layer = layer)
    end do
    
    if(trans_i) call flush(decomp, TRANSPOSE_II)

  end subroutine helm
  
  ! Allocate a parallel 1D Helmholtz solver cache.
  subroutine allocate_h1d_c(h1d_c, decomp_i, decomp_j, ss)
    type(h1d_cache), intent(out) :: h1d_c
    type(decomposition), target, intent(in) :: decomp_i
    type(decomposition), target, intent(in) :: decomp_j
    real(kind = real_8), intent(in) :: ss
    
    integer :: i, ierr, ii_0, j
    real(kind = real_8) :: norm_fact
    real(kind = real_8), dimension(:), allocatable :: diag
    real(kind = real_8), dimension(:, :, :), allocatable :: d
    
    assert(decomp_i%p_ii == 1)
    assert(decomp_j%p_jj == 1)
    
    h1d_c%decomp_i => decomp_i
    h1d_c%decomp_j => decomp_j
       
    h1d_c%ii2 = decomp_j%ii2   ! Number of problems owned by this process
    h1d_c%max_j = decomp_i%jj  ! Maximum j index for the local problems
    if(decomp_i%p_a < 0) h1d_c%max_j = h1d_c%max_j - 1
    
    ! Compute factors used for the local Gaussian elimination
    allocate(diag(decomp_i%ii), h1d_c%ge_fact(decomp_i%ii, decomp_i%jj))
    norm_fact = pi / (2.0D0 * real(decomp_i%ii - 1, kind = real_8))
    do i = 2, decomp_i%ii - 1
      diag(i) = -2.0D0 - ss - 4.0D0 * dsin(real(i - 1, kind = real_8) * norm_fact) ** 2  ! Cache the diagonal of the system as this will be needed later
      h1d_c%ge_fact(i, 2) = diag(i)
    end do
    do j = 3, h1d_c%max_j
      do i = 2, decomp_i%ii - 1
        h1d_c%ge_fact(i, j) = diag(i) - 1.0D0 / h1d_c%ge_fact(i, j - 1)
      end do
    end do
    do j = 2, h1d_c%max_j
      do i = 2, decomp_i%ii - 1
        h1d_c%ge_fact(i, j) = 1.0D0 / h1d_c%ge_fact(i, j)
      end do
    end do
    
    ! Compute boundary contributions
    allocate(h1d_c%psi1(decomp_i%ii, decomp_i%jj))  ! Bc equal to one on the *left*, zero on the right
    do i = 2, decomp_i%ii - 1
      h1d_c%psi1(i, 2) = -1.0D0
      do j = 3, h1d_c%max_j
        h1d_c%psi1(i, j) = 0.0D0
      end do
    end do
    call ge_local(h1d_c, h1d_c%psi1)
    do i = 2, decomp_i%ii - 1
      h1d_c%psi1(i, 1) = 1.0D0
    end do
    
    allocate(h1d_c%psi2(decomp_i%ii, decomp_i%jj))  ! Bc equal to one on the *right*, zero on the left
    do i = 2, decomp_i%ii - 1
      do j = 2, h1d_c%max_j - 1
        h1d_c%psi2(i, j) = 0.0D0
      end do
      h1d_c%psi2(i, h1d_c%max_j) = -1.0D0
    end do
    call ge_local(h1d_c, h1d_c%psi2)
    do i = 2, decomp_i%ii - 1
      h1d_c%psi2(i, 1) = 0.0D0
    end do
    
    ! Set up communication data
    allocate(h1d_c%recv(h1d_c%ii2, decomp_j%n, 3), h1d_c%coeff(decomp_i%ii, 2))
    allocate(h1d_c%scounts1(decomp_j%n), h1d_c%scounts2(decomp_j%n), &
      & h1d_c%sindices(decomp_j%n), h1d_c%rcounts1(decomp_j%n), &
      & h1d_c%rcounts2(decomp_j%n), h1d_c%rindices(decomp_j%n))
    h1d_c%sindices(1) = 1
    h1d_c%rindices(1) = 0
    do i = 1, decomp_j%n - 1
      h1d_c%sindices(i + 1) = h1d_c%sindices(i) + decomp_j%ii2s(i)
      h1d_c%rindices(i + 1) = h1d_c%rindices(i) + h1d_c%ii2
    end do
    if(decomp_i%p_c < 0) then
      h1d_c%scounts1 = 0
    else
      h1d_c%scounts1 = decomp_j%ii2s
    end if
    if(decomp_i%p_a < 0) then
      h1d_c%scounts2 = 0
    else
      h1d_c%scounts2 = decomp_j%ii2s
    end if
    if(decomp_j%periodic_j) then
      h1d_c%rcounts1(1) = h1d_c%ii2
      h1d_c%rcounts2(decomp_j%n) = h1d_c%ii2
    else
      h1d_c%rcounts1(1) = 0
      h1d_c%rcounts2(decomp_j%n) = 0
    end if
    h1d_c%rcounts1(2:) = h1d_c%ii2
    h1d_c%rcounts2(:decomp_j%n - 1) = h1d_c%ii2
    call mpi_comm_dup(decomp_j%cart_comm, h1d_c%comm, ierr);  assert(ierr == MPI_SUCCESS)
#ifdef PH1D_BUFFER
    h1d_c%npsi_buffer = max(maxval(h1d_c%scounts1), maxval(h1d_c%scounts2))
    call allmax(h1d_c%npsi_buffer)
    allocate(h1d_c%psi_buffer(h1d_c%npsi_buffer, 2, decomp_j%n, 2))
    h1d_c%nzeta_buffer = maxval(h1d_c%scounts1)
    call allmax(h1d_c%nzeta_buffer)
    allocate(h1d_c%zeta_buffer(h1d_c%nzeta_buffer, decomp_j%n, 2))
    h1d_c%ncoeff_buffer = max(maxval(h1d_c%rcounts1), maxval(h1d_c%rcounts2))
    call allmax(h1d_c%ncoeff_buffer)
    allocate(h1d_c%coeff_buffer(h1d_c%ncoeff_buffer, 2, decomp_j%n, 2)) 

    ! Ensure that the buffers are initialised as the whole buffer will be
    ! communicated
    h1d_c%psi_buffer = INFINITY    
    h1d_c%zeta_buffer = INFINITY
    h1d_c%coeff_buffer = INFINITY
#endif
    
    ! Form global Gaussian elimination coefficients
    allocate(d(3, decomp_j%n, h1d_c%ii2))
    d = 0.0D0
        
    call h1d_update_psi(h1d_c, h1d_c%psi1)
    ii_0 = decomp_j%ii_0
    if(decomp_j%p_b < 0) ii_0 = ii_0 + 1
    do i = 1, h1d_c%ii2
      do j = 1, decomp_j%n
        if(j == 1 .and. .not. decomp_j%periodic_j) then
          d(2, j, i) = 1.0D0
        else
          if(j > 1) then
            d(1, j - 1, i) = d(1, j - 1, i) + h1d_c%recv(i, j - 1, 2)
          else if(decomp_j%periodic_j) then
            d(1, decomp_j%n, i) = d(1, decomp_j%n, i) + h1d_c%recv(i, decomp_j%n, 2)
          end if
          d(2, j, i) = d(2, j, i) + diag(i + ii_0) + h1d_c%recv(i, j, 1)
        end if
      end do
    end do
    call h1d_update_psi(h1d_c, h1d_c%psi2)
    do i = 1, h1d_c%ii2
      do j = 1, decomp_j%n
        if(j > 1 .or. decomp_j%periodic_j) then
          if(j > 1) then
            d(2, j, i) = d(2, j, i) + h1d_c%recv(i, j - 1, 2)
          else if(decomp_j%periodic_j) then
            d(2, j, i) = d(2, j, i) + h1d_c%recv(i, decomp_j%n, 2)
          end if
          if(j < decomp_j%n) then
            d(3, j + 1, i) = d(3, j + 1, i) + h1d_c%recv(i, j, 1)
          else
            d(3, 1, i) = d(3, 1, i) + h1d_c%recv(i, j, 1)
          end if
        end if
      end do
    end do
    deallocate(diag)
    
    if(decomp_j%periodic_j .and. decomp_j%n > 2) then
      allocate(h1d_c%pge_fact(5, decomp_j%n, h1d_c%ii2))
    else
      allocate(h1d_c%pge_fact(3, decomp_j%n, h1d_c%ii2))
    end if
    do i = 1, h1d_c%ii2
      if(decomp_j%periodic_j) then
        select case(decomp_j%n)
          case(1)
            d(2, 1, i) = d(1, 1, i) + d(2, 1, i) + d(3, 1, i)
          case(2)
            d(1, 1, i) = d(1, 1, i) + d(3, 1, i)
            d(3, 2, i) = d(1, 2, i) + d(3, 2, i)
          case default
        end select
      end if
    
      do j = 1, decomp_j%n
        h1d_c%pge_fact(2, j, i) = d(2, j, i)
      end do
      do j = 2, decomp_j%n
        h1d_c%pge_fact(1, j, i) = -d(1, j - 1, i) / h1d_c%pge_fact(2, j - 1, i)
        h1d_c%pge_fact(2, j, i) = h1d_c%pge_fact(2, j, i) + h1d_c%pge_fact(1, j, i) * d(3, j, i)
      end do
      do j = 1, decomp_j%n
        h1d_c%pge_fact(2, j, i) = 1.0D0 / h1d_c%pge_fact(2, j, i)
      end do
      do j = decomp_j%n - 1, 1, -1
        h1d_c%pge_fact(3, j, i) = -d(3, j + 1, i)
      end do
      
      if(decomp_j%periodic_j .and. decomp_j%n > 2) then
        h1d_c%pge_fact(4, 1, i) = d(2, 1, i)
        h1d_c%pge_fact(4, 2, i) = d(1, 1, i)
        h1d_c%pge_fact(4, 3:decomp_j%n - 1, i) = 0.0D0
        h1d_c%pge_fact(4, decomp_j%n, i) = d(3, 1, i)
        
        do j = 2, decomp_j%n
          h1d_c%pge_fact(4, j, i) = h1d_c%pge_fact(4, j, i) + h1d_c%pge_fact(4, j - 1, i) * h1d_c%pge_fact(1, j, i)
        end do
        h1d_c%pge_fact(4, decomp_j%n, i) = h1d_c%pge_fact(4, decomp_j%n, i) * h1d_c%pge_fact(2, decomp_j%n, i)
        do j = decomp_j%n - 1, 1, -1
          h1d_c%pge_fact(4, j, i) = (h1d_c%pge_fact(4, j, i) + h1d_c%pge_fact(3, j, i) * h1d_c%pge_fact(4, j + 1, i)) * h1d_c%pge_fact(2, j, i)
        end do   
        
        h1d_c%pge_fact(4, 1, i) = 1.0D0 / h1d_c%pge_fact(4, 1, i)
        do j = 2, decomp_j%n
          h1d_c%pge_fact(4, j, i) = -h1d_c%pge_fact(4, j, i)
        end do  
        
        h1d_c%pge_fact(5, 1, i) = d(1, decomp_j%n, i)
        h1d_c%pge_fact(5, 2:decomp_j%n - 2, i) = 0.0D0
        h1d_c%pge_fact(5, decomp_j%n - 1, i) = d(3, decomp_j%n, i)
        h1d_c%pge_fact(5, decomp_j%n, i) = d(2, decomp_j%n, i)
        
        do j = 2, decomp_j%n
          h1d_c%pge_fact(5, j, i) = h1d_c%pge_fact(5, j, i) + h1d_c%pge_fact(5, j - 1, i) * h1d_c%pge_fact(1, j, i)
        end do
        h1d_c%pge_fact(5, decomp_j%n, i) = h1d_c%pge_fact(5, decomp_j%n, i) * h1d_c%pge_fact(2, decomp_j%n, i)
        do j = decomp_j%n - 1, 1, -1
          h1d_c%pge_fact(5, j, i) = (h1d_c%pge_fact(5, j, i) + h1d_c%pge_fact(3, j, i) * h1d_c%pge_fact(5, j + 1, i)) * h1d_c%pge_fact(2, j, i)
        end do     
        
        h1d_c%pge_fact(5, 1, i) = h1d_c%pge_fact(4, 1, i) * h1d_c%pge_fact(5, 1, i)
        do j = 2, decomp_j%n
          h1d_c%pge_fact(5, j, i) = h1d_c%pge_fact(5, j, i) + h1d_c%pge_fact(4, j, i) * h1d_c%pge_fact(5, 1, i)
        end do
        
        h1d_c%pge_fact(5, decomp_j%n, i) = 1.0D0 / h1d_c%pge_fact(5, decomp_j%n, i)
        do j = decomp_j%n - 1, 1, -1
          h1d_c%pge_fact(5, j, i) = -h1d_c%pge_fact(5, j, i)
        end do
      end if
    end do    
    deallocate(d)
    
    ! Allocate a RHS
    allocate(h1d_c%rhs(h1d_c%ii2, decomp_j%n + 1))
    
  end subroutine allocate_h1d_c
    
  ! Deallocate a parallel 1D Helmholtz solver cache.
  subroutine deallocate_h1d_c(h1d_c)
    type(h1d_cache), intent(inout) :: h1d_c
  
    integer :: ierr
  
    deallocate(h1d_c%ge_fact)
    deallocate(h1d_c%psi1)
    deallocate(h1d_c%psi2)
    deallocate(h1d_c%recv, h1d_c%coeff)
    deallocate(h1d_c%rcounts1, h1d_c%rcounts2, h1d_c%rindices, h1d_c%scounts1, &
      & h1d_c%scounts2, h1d_c%sindices)
#ifdef PH1D_BUFFER
    deallocate(h1d_c%psi_buffer, h1d_c%zeta_buffer, h1d_c%coeff_buffer)
#endif
    call mpi_comm_free(h1d_c%comm, ierr);  assert(ierr == MPI_SUCCESS)
    deallocate(h1d_c%pge_fact)
    deallocate(h1d_c%rhs)
  
  end subroutine deallocate_h1d_c
  
  ! Perform local Gaussian elimination operations.
  subroutine ge_local(h1d_c, psi)
    type(h1d_cache), intent(in) :: h1d_c
    real(kind = real_8), dimension(:, :), intent(inout) :: psi
    
    integer :: i, j
    type(decomposition), pointer :: decomp_i
    
    assert(size(psi, 1) == h1d_c%decomp_i%ii)
    assert(size(psi, 2) == h1d_c%decomp_i%jj)
    
    decomp_i => h1d_c%decomp_i
    
    ! Lower triangular Gaussian elimination
    do j = 3, h1d_c%max_j
      do i = 2, decomp_i%ii - 1
        psi(i, j) = psi(i, j) - psi(i, j - 1) * h1d_c%ge_fact(i, j - 1)
      end do
    end do
    ! Upper triangular Gaussian elimination, combined with normalisation
    do i = 2, decomp_i%ii - 1
      psi(i, h1d_c%max_j) = psi(i, h1d_c%max_j) * h1d_c%ge_fact(i, h1d_c%max_j)
    end do
    do j = h1d_c%max_j - 1, 2, -1
      do i = 2, decomp_i%ii - 1
        psi(i, j) = (psi(i, j) - psi(i, j + 1)) * h1d_c%ge_fact(i, j)
      end do
    end do
    do i = 2, decomp_i%ii - 1
      psi(i, 1) = 0.0D0
    end do
    
  end subroutine ge_local
  
  ! Perform global 1D Helmholtz solves. The algorithm is as described in 
  ! Schumann and Strietzel 1995 (JSC). Each problem (in the j-direction) is
  ! decomposed into a series of local homogeneous and particular solutions, and
  ! the resulting boundary component coefficients for each problem are
  ! determined by the solution of a global problem on one process. The global
  ! problems are themselves decomposed (in the i-direction).
  subroutine h1d(h1d_c, zeta)
    type(h1d_cache), intent(inout) :: h1d_c
    real(kind = real_8), dimension(:, :), intent(inout) :: zeta
    
    integer :: i, j, j1, j2
    type(decomposition), pointer :: decomp_i, decomp_j
    
    assert(size(zeta, 1) == h1d_c%decomp_i%ii)
    assert(size(zeta, 2) == h1d_c%decomp_i%jj)
        
    decomp_i => h1d_c%decomp_i
    decomp_j => h1d_c%decomp_j
        
    ! Communicate global GE RHS
    call h1d_update_zeta(h1d_c, zeta)
    ! Local GE
    call ge_local(h1d_c, zeta)
    ! Communicate global GE RHS
    call h1d_update_psi(h1d_c, zeta)
    
    do i = 1, h1d_c%ii2
      ! Form global GE RHS
      do j2 = 1, decomp_j%n
        j1 = j2 - 1
        if(j1 == 0) j1 = decomp_j%n
        if(j2 == 1 .and. .not. decomp_j%periodic_j) then
          h1d_c%rhs(i, j2) = 0.0D0
        else
          h1d_c%rhs(i, j2) = h1d_c%recv(i, j2, 3) - h1d_c%recv(i, j1, 2) - h1d_c%recv(i, j2, 1)
        end if
      end do
      
      ! Perform global GE
      if(decomp_j%periodic_j) then
        ! Periodic case. Two additional sweeps of GE are required to handle the
        ! periodic elements of the matrix.        
        
        do j = 2, decomp_j%n
          h1d_c%rhs(i, j) = h1d_c%rhs(i, j) + h1d_c%rhs(i, j - 1) * h1d_c%pge_fact(1, j, i)
        end do
        h1d_c%rhs(i, decomp_j%n) = h1d_c%rhs(i, decomp_j%n) * h1d_c%pge_fact(2, decomp_j%n, i)
        do j = decomp_j%n - 1, 1, -1
          h1d_c%rhs(i, j) = (h1d_c%rhs(i, j) + h1d_c%pge_fact(3, j, i) * h1d_c%rhs(i, j + 1)) * h1d_c%pge_fact(2, j, i)
        end do
        
        if(decomp_j%n > 2) then
          h1d_c%rhs(i, 1) = h1d_c%rhs(i, 1) * h1d_c%pge_fact(4, 1, i)
          do j = 2, decomp_j%n
            h1d_c%rhs(i, j) = h1d_c%rhs(i, j) + h1d_c%pge_fact(4, j, i) * h1d_c%rhs(i, 1)
          end do
          
          h1d_c%rhs(i, decomp_j%n) = h1d_c%rhs(i, decomp_j%n) * h1d_c%pge_fact(5, decomp_j%n, i)
          do j = decomp_j%n - 1, 1, -1
            h1d_c%rhs(i, j) = h1d_c%rhs(i, j) + h1d_c%pge_fact(5, j, i) * h1d_c%rhs(i, decomp_j%n)
          end do
        end if
      else
        ! Non-periodic case
      
        do j = 3, decomp_j%n
          h1d_c%rhs(i, j) = h1d_c%rhs(i, j) + h1d_c%rhs(i, j - 1) * h1d_c%pge_fact(1, j, i)
        end do
        h1d_c%rhs(i, decomp_j%n) = h1d_c%rhs(i, decomp_j%n) * h1d_c%pge_fact(2, decomp_j%n, i)
        do j = decomp_j%n - 1, 2, -1
          h1d_c%rhs(i, j) = (h1d_c%rhs(i, j) + h1d_c%pge_fact(3, j, i) * h1d_c%rhs(i, j + 1)) * h1d_c%pge_fact(2, j, i)
        end do
      end if
    end do
    
    if(decomp_j%periodic_j .and. decomp_j%n > 1) then
      ! One additional periodic coefficient requires communication in the
      ! parallel periodic case.
      do i = 1, h1d_c%ii2
        h1d_c%rhs(i, decomp_j%n + 1) = h1d_c%rhs(i, 1)
      end do
    end if
    
    ! Communicate the solution of the global GE.
    call h1d_update_coeff(h1d_c, h1d_c%rhs)
    
    ! Add on the boundary contributions.
    if(decomp_i%p_c >= 0) then
      if(decomp_i%p_a >= 0) then
        do j = decomp_i%l_j_min, decomp_i%l_j_max
          do i = decomp_i%l_i_min, decomp_i%l_i_max
            zeta(i, j) = zeta(i, j) + h1d_c%coeff(i, 1) * h1d_c%psi1(i, j) + h1d_c%coeff(i, 2) * h1d_c%psi2(i, j)
          end do
        end do
      else
        do j = decomp_i%l_j_min, decomp_i%l_j_max
          do i = decomp_i%l_i_min, decomp_i%l_i_max
            zeta(i, j) = zeta(i, j) + h1d_c%coeff(i, 1) * h1d_c%psi1(i, j)
          end do
        end do
      end if
    else if(decomp_i%p_a >= 0) then
      do j = decomp_i%l_j_min, decomp_i%l_j_max
        do i = decomp_i%l_i_min, decomp_i%l_i_max
          zeta(i, j) = zeta(i, j) + h1d_c%coeff(i, 2) * h1d_c%psi2(i, j)
        end do
      end do
    end if        
        
  end subroutine h1d
  
  ! Perform parallel 1D Helmholtz solver communications. The communicated data
  ! is received on the h1d_cache.
  subroutine h1d_update_psi(h1d_c, psi)
    type(h1d_cache), target, intent(inout) :: h1d_c
    real(kind = real_8), dimension(:, :), intent(in) :: psi
    
    integer :: ierr
    type(decomposition), pointer :: decomp_i, decomp_j
#ifdef PH1D_BUFFER
    integer :: i, j
#endif
    
    assert(size(psi, 1) == h1d_c%decomp_i%ii)
    assert(size(psi, 2) == h1d_c%decomp_i%jj)
    
    decomp_i => h1d_c%decomp_i
    decomp_j => h1d_c%decomp_j

#ifdef DDEBUG
    h1d_c%recv(:, :, 1:2) = INFINITY
#ifdef PH1D_BUFFER
    h1d_c%psi_buffer = INFINITY
#endif
#endif

    if(decomp_j%n == 1) then
      if(decomp_j%periodic_j) then
        h1d_c%recv(:, 1, 1) = psi(2:decomp_i%ii - 1, 2)
        h1d_c%recv(:, 1, 2) = psi(2:decomp_i%ii - 1, h1d_c%max_j)
      end if
    else
#ifdef PH1D_BUFFER
      do i = 1, decomp_j%n
        do j = 1, h1d_c%scounts1(i)
          h1d_c%psi_buffer(j, 1, i, 1) = psi(h1d_c%sindices(i) + j, 2)
        end do
        do j = 1, h1d_c%scounts2(i)
          h1d_c%psi_buffer(j, 2, i, 1) = psi(h1d_c%sindices(i) + j, h1d_c%max_j)
        end do
      end do
      
      wtime0 = mpi_wtime()
      call mpi_alltoall(h1d_c%psi_buffer(1, 1, 1, 1), 2 * h1d_c%npsi_buffer, &
        & MPI_DOUBLE_PRECISION, h1d_c%psi_buffer(1, 1, 1, 2), &
        & 2 * h1d_c%npsi_buffer, MPI_DOUBLE_PRECISION, h1d_c%comm, ierr);  assert(ierr == MPI_SUCCESS)
      h1d_time = h1d_time + mpi_wtime() - wtime0
      
      do i = 1, decomp_j%n
        do j = 1, h1d_c%rcounts1(i)
          h1d_c%recv(j, i, 1) = h1d_c%psi_buffer(j, 1, i, 2)
        end do
        do j = 1, h1d_c%rcounts2(i)
          h1d_c%recv(j, i, 2) = h1d_c%psi_buffer(j, 2, i, 2)
        end do
      end do
#else
      wtime0 = mpi_wtime()
      call mpi_alltoallv(psi(1, 2), h1d_c%scounts1, h1d_c%sindices, &
        & MPI_DOUBLE_PRECISION, h1d_c%recv(1, 1, 1), h1d_c%rcounts1, &
        & h1d_c%rindices, MPI_DOUBLE_PRECISION, h1d_c%comm, ierr);  assert(ierr == MPI_SUCCESS)
      call mpi_alltoallv(psi(1, h1d_c%max_j), h1d_c%scounts2, h1d_c%sindices, &
        & MPI_DOUBLE_PRECISION, h1d_c%recv(1, 1, 2), h1d_c%rcounts2, &
        & h1d_c%rindices, MPI_DOUBLE_PRECISION, h1d_c%comm, ierr);  assert(ierr == MPI_SUCCESS)
      h1d_time = h1d_time + mpi_wtime() - wtime0
#endif
    end if
    
  end subroutine h1d_update_psi
  
  ! Perform parallel 1D Helmholtz solver communications. The communicated data
  ! is received on the h1d_cache.
  subroutine h1d_update_zeta(h1d_c, zeta)
    type(h1d_cache), target, intent(inout) :: h1d_c
    real(kind = real_8), dimension(:, :), intent(in) :: zeta
    
    integer :: ierr
    type(decomposition), pointer :: decomp_i, decomp_j
#ifdef PH1D_BUFFER
    integer :: i, j
#endif    

    assert(size(zeta, 1) == h1d_c%decomp_i%ii)
    assert(size(zeta, 2) == h1d_c%decomp_i%jj)
    
    decomp_i => h1d_c%decomp_i
    decomp_j => h1d_c%decomp_j

#ifdef DDEBUG
    h1d_c%recv(:, :, 3) = INFINITY
#ifdef PH1D_BUFFER
    h1d_c%zeta_buffer = INFINITY
#endif
#endif

    if(decomp_j%n == 1) then
      if(decomp_j%periodic_j) then
        h1d_c%recv(:, 1, 3) = zeta(2:decomp_i%ii - 1, 1)
      end if
    else
#ifdef PH1D_BUFFER
      do i = 1, decomp_j%n
        do j = 1, h1d_c%scounts1(i)
          h1d_c%zeta_buffer(j, i, 1) = zeta(h1d_c%sindices(i) + j, 1)
        end do
      end do

      wtime0 = mpi_wtime()
      call mpi_alltoall(h1d_c%zeta_buffer(1, 1, 1), h1d_c%nzeta_buffer, &
        & MPI_DOUBLE_PRECISION, h1d_c%zeta_buffer(1, 1, 2), &
        & h1d_c%nzeta_buffer, MPI_DOUBLE_PRECISION, h1d_c%comm, ierr);  assert(ierr == MPI_SUCCESS)
      h1d_time = h1d_time + mpi_wtime() - wtime0

      do i = 1, decomp_j%n
        do j = 1, h1d_c%rcounts1(i)
          h1d_c%recv(j, i, 3) = h1d_c%zeta_buffer(j, i, 2)
        end do
      end do
#else
      wtime0 = mpi_wtime()
      call mpi_alltoallv(zeta, h1d_c%scounts1, h1d_c%sindices, &
        & MPI_DOUBLE_PRECISION, h1d_c%recv(1, 1, 3), h1d_c%rcounts1, &
        & h1d_c%rindices, MPI_DOUBLE_PRECISION, h1d_c%comm, ierr);  assert(ierr == MPI_SUCCESS)
      h1d_time = h1d_time + mpi_wtime() - wtime0
#endif
    end if
    
  end subroutine h1d_update_zeta
  
  ! Perform parallel 1D Helmholtz solver communications. The communicated data
  ! is received on the h1d_cache.
  subroutine h1d_update_coeff(h1d_c, c)
    type(h1d_cache), target, intent(inout) :: h1d_c
    real(kind = real_8), dimension(:, :), intent(in) :: c
    
    integer :: ierr
    type(decomposition), pointer :: decomp_i, decomp_j
#ifdef PH1D_BUFFER    
    integer :: i, j
#endif

    assert(size(c, 1) == h1d_c%ii2)
    assert(size(c, 2) == h1d_c%decomp_j%n + 1)
    
    decomp_i => h1d_c%decomp_i
    decomp_j => h1d_c%decomp_j

#ifdef DDEBUG
    h1d_c%coeff = INFINITY
#ifdef PH1D_BUFFER
    h1d_c%coeff_buffer = INFINITY
#endif
#endif

    if(decomp_j%n == 1) then
      if(decomp_j%periodic_j) then
        h1d_c%coeff(2:decomp_i%ii - 1, 1) = c(:, 1)
        h1d_c%coeff(2:decomp_i%ii - 1, 2) = c(:, 1)
      end if
    else
#ifdef PH1D_BUFFER
      do i = 1, decomp_j%n
        do j = 1, h1d_c%rcounts1(i)
          h1d_c%coeff_buffer(j, 1, i, 1) = c(j, i)
        end do
        do j = 1, h1d_c%rcounts2(i)
          h1d_c%coeff_buffer(j, 2, i, 1) = c(j, i + 1)
        end do
      end do

      wtime0 = mpi_wtime()
      call mpi_alltoall(h1d_c%coeff_buffer(1, 1, 1, 1), &
        & 2 * h1d_c%ncoeff_buffer, MPI_DOUBLE_PRECISION, &
        & h1d_c%coeff_buffer(1, 1, 1, 2), 2 * h1d_c%ncoeff_buffer, &
        & MPI_DOUBLE_PRECISION, h1d_c%comm, ierr);  assert(ierr == MPI_SUCCESS)
      h1d_time = h1d_time + mpi_wtime() - wtime0
      
      do i = 1, decomp_j%n
        do j = 1, h1d_c%scounts1(i)
          h1d_c%coeff(h1d_c%sindices(i) + j, 1) = h1d_c%coeff_buffer(j, 1, i, 2)
        end do
        do j = 1, h1d_c%scounts2(i)
          h1d_c%coeff(h1d_c%sindices(i) + j, 2) = h1d_c%coeff_buffer(j, 2, i, 2)
        end do
      end do
#else
      wtime0 = mpi_wtime()
      call mpi_alltoallv(c(1, 1), h1d_c%rcounts1, h1d_c%rindices, &
        & MPI_DOUBLE_PRECISION, h1d_c%coeff(1, 1), h1d_c%scounts1, &
        & h1d_c%sindices, MPI_DOUBLE_PRECISION, h1d_c%comm, ierr);  assert(ierr == MPI_SUCCESS)
      call mpi_alltoallv(c(1, 2), h1d_c%rcounts2, h1d_c%rindices, &
        & MPI_DOUBLE_PRECISION, h1d_c%coeff(1, 2), h1d_c%scounts2, &
        & h1d_c%sindices, MPI_DOUBLE_PRECISION, h1d_c%comm, ierr);  assert(ierr == MPI_SUCCESS)
      h1d_time = h1d_time + mpi_wtime() - wtime0
#endif
    end if

  end subroutine h1d_update_coeff

  subroutine sint(sint_c, x)
    type(sint_cache), intent(inout) :: sint_c
    real(kind = real_8), dimension(:), intent(inout) :: x

    assert(size(x) == sint_c%n + 1)

    call sint1(sint_c, x)

  end subroutine sint

  subroutine sint1(sint_c, x)
    type(sint_cache), intent(inout) :: sint_c
    real(kind = real_8), dimension(:), intent(inout) :: x

    integer :: i, k, kc, modn, n, np1, ns2
    real(kind = real_8) :: xhold, t1, t2

    ! Optimisations, based on testing using gcc 4.4, -O3 -ffast-math
    !   1. Some line re-ordering

    assert(size(x) == sint_c%n + 1)
          
    n = sint_c%n
    assert(n >= 2)
    if(n < 2) then
      x(1) = x(1) + x(1)
    else if(n == 2) then
      xhold = sqrt3 * (x(1) + x(2))
      x(2) = sqrt3 * (x(1) - x(2))
      x(1) = xhold
    else    
      np1 = n + 1
      ns2 = n / 2

      sint_c%cc(1) = 0.0D0
      do k = 1, ns2
        kc = np1 - k
        t1 = x(k) - x(kc)
        t2 = sint_c%was(k) * (x(k) + x(kc))
        sint_c%cc(k + 1) = t1 + t2
        sint_c%cc(kc + 1) = t2 - t1
      end do
      modn = mod(n, 2)
      if(modn /= 0) sint_c%cc(ns2 + 2) = 4.0D0 * x(ns2 + 1)
      
      call rfftf1(sint_c%rfft_c, sint_c%cc, x)

      x(1) = 0.5D0 * sint_c%cc(1)
      do i = 3, n, 2
        x(i) = x(i - 2) + sint_c%cc(i - 1)
        x(i - 1) = -sint_c%cc(i)
      end do

      if(modn == 0) x(n) = -sint_c%cc(n + 1)
    end if

  end subroutine sint1
      
  subroutine rfftf1(rfft_c, cc, ch)
    type(rfft_cache), intent(inout) :: rfft_c
    real(kind = real_8), dimension(:), intent(inout) :: cc
    real(kind = real_8), dimension(:), intent(inout) :: ch
    
    integer :: idl1, ido, ip, iw, ix2, ix3, i, kh, l1, l2, n, nf

    assert(size(cc) == rfft_c%n)
    assert(size(ch) == rfft_c%n)

    n = rfft_c%n
    nf = rfft_c%nf
    l2 = n
    iw = n
    do i = 1, nf
      kh = nf - i
      ip = rfft_c%ifac(kh + 1)
      l1 = l2 / ip
      ido = n / l2
      idl1 = ido * l1
      iw = iw - (ip - 1) * ido
      select case(ip)
        case(4)
          ix2 = iw + ido
          ix3 = ix2 + ido
          if(mod(i, 2) == 1) then
            call radf4(ido, l1, cc, ch, rfft_c%wa(iw:ix2 - 1), rfft_c%wa(ix2:ix3 - 1), rfft_c%wa(ix3:))
          else
            call radf4(ido, l1, ch, cc, rfft_c%wa(iw:ix2 - 1), rfft_c%wa(ix2:ix3 - 1), rfft_c%wa(ix3:))
          end if
        case(2)
          if(mod(i, 2) == 1) then
            call radf2(ido, l1, cc, ch, rfft_c%wa(iw:))
          else
            call radf2(ido, l1, ch, cc, rfft_c%wa(iw:))
          end if
        case default
         labort('Missing "radf"... subroutines')
      end select
      
      l2 = l1
    end do
    
    if(mod(nf, 2) == 0) return

    cc = ch

  end subroutine rfftf1

  subroutine radf2(ido, l1, cc, ch, wa1)
    integer, intent(in) :: ido
    integer, intent(in) :: l1
    real(kind = real_8), dimension(ido, l1, 2), intent(in) :: cc
    real(kind = real_8), dimension(ido, 2, l1), intent(out) :: ch
    real(kind = real_8), dimension(*), intent(in) :: wa1

    ! Optimisations, based on testing using gcc 4.4, -O3 -ffast-math
    !   1. dimension(*)

    integer :: i, ic, idp2, k
    real(kind = real_8) :: ti2, tr2
        
    do k = 1, l1
       ch(1, 1, k) = cc(1, k, 1) + cc(1, k, 2)
       ch(ido, 2, k) = cc(1, k, 1) - cc(1, k, 2)
       ch(ido, 2, k) = cc(1, k, 1) - cc(1, k, 2)
    end do

    if(ido > 2) then
      idp2 = ido + 2

      do k = 1, l1
         do i = 3, ido, 2
           ic = idp2 - i
           tr2 = wa1(i - 2) * cc(i - 1, k, 2) + wa1(i - 1) * cc(i, k, 2)
           ti2 = wa1(i - 2) * cc(i, k, 2) - wa1(i - 1) * cc(i - 1, k, 2)
           ch(i, 1, k) = cc(i, k, 1) + ti2
           ch(ic, 2, k) = ti2 - cc(i, k, 1)
           ch(i - 1, 1, k) = cc(i - 1, k, 1) + tr2
           ch(ic - 1, 2, k) = cc(i - 1, k, 1) - tr2
         end do
      end do
    end if

    if(mod(ido, 2) == 0) then
      do k = 1, l1
         ch(1, 2, k) = -cc(ido, k, 2)
         ch(ido, 1, k) = cc(ido, k, 1)
      end do
    end if

  end subroutine radf2

  subroutine radf4(ido, l1, cc, ch, wa1, wa2, wa3)
    integer, intent(in) :: ido
    integer, intent(in) :: l1
    real(kind = real_8), dimension(ido, l1, 4), intent(in) :: cc
    real(kind = real_8), dimension(ido, 4, l1), intent(out) :: ch
    real(kind = real_8), dimension(*), intent(in) :: wa1
    real(kind = real_8), dimension(*), intent(in) :: wa2
    real(kind = real_8), dimension(*), intent(in) :: wa3

    ! Optimisations, based on testing using gcc 4.4, -O3 -ffast-math
    !   1. dimension(*)

    integer :: i, ic, idp2, k
    real(kind = real_8) :: ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, &
      & tr1, tr2, tr3, tr4     

    do k = 1, l1
       tr1 = cc(1, k, 2) + cc(1, k, 4)
       tr2 = cc(1, k, 1) + cc(1, k, 3)
       ch(1, 1, k) = tr1 + tr2
       ch(ido, 4, k) = tr2 - tr1
       ch(ido, 2, k) = cc(1, k, 1) - cc(1, k, 3)
       ch(1, 3, k) = cc(1, k, 4) - cc(1, k, 2)
    end do

    if(ido > 2) then
      idp2 = ido + 2

      do k = 1, l1
        do i = 3, ido, 2
          ic = idp2 - i
          cr2 = wa1(i - 2) * cc(i - 1, k, 2) + wa1(i - 1) * cc(i, k, 2)
          ci2 = wa1(i - 2) * cc(i, k, 2) - wa1(i - 1) * cc(i - 1, k, 2)
          cr2 = wa1(i - 2) * cc(i - 1, k, 2) + wa1(i - 1) * cc(i, k, 2)
          ci2 = wa1(i - 2) * cc(i, k, 2) - wa1(i - 1) * cc(i - 1, k, 2)
          cr3 = wa2(i - 2) * cc(i - 1, k, 3) + wa2(i - 1) * cc(i, k, 3)
          ci3 = wa2(i - 2) * cc(i, k, 3) - wa2(i - 1) * cc(i - 1, k, 3)
          cr4 = wa3(i - 2) * cc(i - 1, k, 4) + wa3(i - 1) * cc(i, k, 4)
          ci4 = wa3(i - 2) * cc(i, k, 4) - wa3(i - 1) * cc(i - 1, k, 4)
          tr1 = cr2 + cr4
          tr4 = cr4 - cr2
          ti1 = ci2 + ci4
          ti4 = ci2 - ci4
          ti2 = cc(i, k, 1) + ci3
          ti3 = cc(i, k, 1) - ci3
          tr2 = cc(i - 1, k, 1) + cr3
          tr3 = cc(i - 1, k, 1) - cr3
          ch(i - 1, 1, k) = tr1 + tr2
          ch(ic - 1, 4, k) = tr2 - tr1
          ch(i, 1, k) = ti1 + ti2
          ch(ic, 4, k) = ti1 - ti2
          ch(i - 1, 3, k) = ti4 + tr3
          ch(ic - 1, 2, k) = tr3 - ti4
          ch(i, 3, k) = tr4 + ti3
          ch(ic, 2, k) = tr4 - ti3
        end do
      end do
    end if

    if(mod(ido, 2) == 0) then
      do k = 1, l1
         ti1 = -hsqt2 * (cc(ido, k, 2) + cc(ido, k, 4))
         tr1 =  hsqt2 * (cc(ido, k, 2) - cc(ido, k, 4))
         ch(ido, 1, k) = tr1 + cc(ido, k, 1)
         ch(ido, 3, k) = cc(ido, k, 1) - tr1
         ch(1, 2, k) = ti1 - cc(ido, k, 3)
         ch(1, 4, k) = ti1 + cc(ido, k, 3)
      end do
   end if

  end subroutine radf4
 
end module qg_helmholtz
