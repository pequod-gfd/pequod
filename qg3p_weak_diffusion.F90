! Experimental barotropic Crank--Nicolson source substep with no-slip walls.
module qg_weak_diffusion
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use lapack_interfaces, only: dpotrf, dpotrs, dpocon
  use qg_irregular_transport, only: source_momentum_budget
  use qg_cmm, only: cmm_cache, setup_cmm, cmm_solve, release_cmm
  use qg_geometry, only: fluid_node, wall_node
  use qg_irregular_parallel, only: distributed_geometry, broadcast_error
  use qg_irregular_operators, only: viscous_face_plan, setup_viscous_faces, viscous_open_flux, &
                                    stiffness_flux, sum_open_flux, mean_to_pv_nodes
  use qg_parallel, only: decomposition, parallel_data, request_cache, allocate, deallocate, update, flush
  use qg_parameters, only: real_8
  implicit none
  include 'mpif.h'
  private
  public :: weak_diffusion_cache, setup_weak_diffusion, advance_weak_diffusion, release_weak_diffusion

  type weak_diffusion_cache
    private
    type(distributed_geometry), pointer :: geom => null()
    type(cmm_cache) :: cmm
    type(viscous_face_plan) :: viscous
    type(parallel_data) :: halo
    type(request_cache), pointer :: requests => null()
    real(real_8), allocatable :: flux(:,:),harmonic(:,:,:,:),factor(:,:)
  end type weak_diffusion_cache

contains

  subroutine setup_weak_diffusion(cache,decomp,geom,error)
    type(weak_diffusion_cache), intent(inout) :: cache
    type(decomposition), target, intent(in) :: decomp
    type(distributed_geometry), target, intent(in) :: geom
    character(*), intent(out) :: error
    real(real_8), allocatable :: nodal(:,:,:),zero(:,:,:),lap(:,:,:),trace(:,:),work(:)
    real(real_8) :: scale,rcond,asymmetry
    integer, allocatable :: iwork(:)
    integer :: islands,k,ierr,info

    error = ''
    if (associated(cache%geom)) error = 'weak diffusion cache is already active'
    if (decomp%layers /= 1 .or. decomp%periodic_j) &
      error = 'experimental weak diffusion requires one layer and a closed basin'
    call broadcast_error(geom%comm,geom%rank,geom%ranks,error)
    if (len_trim(error) /= 0) return
    call setup_viscous_faces(geom,cache%viscous,error)
    if (len_trim(error) /= 0) return
    call setup_cmm(cache%cmm,decomp,geom,[0.0D0],error)
    if (len_trim(error) /= 0) return
    cache%geom => geom
    allocate(cache%requests,cache%flux(1,size(geom%faces)))
    call allocate(cache%requests,comm=decomp%cart_comm)
    call allocate(cache%halo,decomp,requests=cache%requests)

    islands = ubound(geom%global%loops,1)
    allocate(cache%harmonic(geom%ni,geom%nj,1,islands),cache%factor(islands,islands))
    if (islands == 0) return
    allocate(nodal(geom%ni,geom%nj,1),zero(geom%ni,geom%nj,1),lap(geom%ni,geom%nj,1),trace(0:islands,1))
    zero = 0.0D0
    do k = 1,islands
      trace = 0.0D0
      trace(k,1) = 1.0D0
      call cmm_solve(cache%cmm,zero,trace,nodal,error)
      if (len_trim(error) /= 0) exit
      call impose_trace(geom,nodal,trace(:,1),cache%harmonic(:,:,:,k))
      call evaluate_stiffness(cache,cache%harmonic(:,:,:,k),lap,error)
      if (len_trim(error) /= 0) exit
      call island_rows(geom,lap,cache%factor(:,k))
    end do
    if (geom%rank == 0 .and. len_trim(error) == 0) then
      asymmetry = maxval(abs(cache%factor-transpose(cache%factor)))
      scale = maxval(sum(abs(cache%factor),dim=1))
      if (asymmetry > 1.0D-11*scale) then
        error = 'no-slip island stiffness influence is asymmetric'
      else
        cache%factor = -cache%factor
        call dpotrf('L',islands,cache%factor,islands,info)
        if (info == 0) then
          allocate(work(3*islands),iwork(islands))
          call dpocon('L',islands,cache%factor,islands,scale,rcond,work,iwork,info)
          if (info == 0) then
            if (rcond < 1.0D-12) info = 1
          end if
        end if
        if (info /= 0) error = 'no-slip island stiffness influence is singular or ill-conditioned'
      end if
    end if
    call mpi_bcast(error,len(error),MPI_CHARACTER,0,geom%comm,ierr)
    if (len_trim(error) /= 0) call release_weak_diffusion(cache)
  end subroutine setup_weak_diffusion

  subroutine impose_trace(geom,nodal,trace,field)
    type(distributed_geometry), intent(in) :: geom
    real(real_8), intent(in) :: nodal(:,:,:),trace(0:)
    real(real_8), intent(out) :: field(:,:,:)
    integer :: i,j,k

    call mean_to_pv_nodes(geom,nodal,field)
    do j = 1,geom%nj
      do i = 1,geom%ni
        k = geom%global%component(geom%i0+i,geom%j0+j)
        if (k >= 0) field(i,j,1) = trace(k)
      end do
    end do
  end subroutine impose_trace

  subroutine evaluate_stiffness(cache,psi,lap,error)
    type(weak_diffusion_cache), intent(inout) :: cache
    real(real_8), intent(in) :: psi(:,:,:)
    real(real_8), intent(out) :: lap(:,:,:)
    character(*), intent(out) :: error

    cache%halo%sdata = psi
    call update(cache%halo)
    call flush(cache%halo)
    call stiffness_flux(cache%geom,cache%halo%data,cache%flux,error)
    if (len_trim(error) == 0) call sum_open_flux(cache%geom,cache%flux,lap)
  end subroutine evaluate_stiffness

  subroutine island_rows(geom,field,rows)
    type(distributed_geometry), intent(in) :: geom
    real(real_8), intent(in) :: field(:,:,:)
    real(real_8), intent(out) :: rows(:)
    real(real_8) :: local(size(rows))
    integer :: i,j,k,ierr

    local = 0.0D0
    do j = 1,geom%nj
      do i = 1,geom%ni
        k = geom%global%component(geom%i0+i,geom%j0+j)
        if (k > 0) local(k) = local(k)+field(i,j,1)
      end do
    end do
    call mpi_allreduce(local,rows,size(rows),MPI_DOUBLE_PRECISION,MPI_SUM,geom%comm,ierr)
  end subroutine island_rows

  ! Poisson closure sets the no-slip weak wall derivative to zero.
  ! Consequently the compatible vorticity mean is B psi/A on every cell.
  subroutine evaluate_diffusion(cache,psi,lap,omega,open_divergence,error)
    type(weak_diffusion_cache), intent(inout) :: cache
    real(real_8), intent(in) :: psi(:,:,:)
    real(real_8), intent(out) :: lap(:,:,:),omega(:,:,:),open_divergence(:,:,:)
    character(*), intent(out) :: error
    integer :: i,j

    associate(g => cache%geom)
      call evaluate_stiffness(cache,psi,lap,error)
      if (len_trim(error) /= 0) return
      omega = 0.0D0
      do j = 1,g%nj
        do i = 1,g%ni
          if (g%area(i,j)>0.0D0) omega(i,j,1) = lap(i,j,1)*g%h**2/g%area(i,j)
        end do
      end do
      call viscous_open_flux(g,cache%viscous,omega,cache%flux,error)
      if (len_trim(error) /= 0) return
      call sum_open_flux(g,cache%flux,open_divergence)
    end associate
  end subroutine evaluate_diffusion

  ! A single source substep in native units, with the outer gauge fixed to
  ! zero. source contains time-averaged integrated curl over each wet cell.
  ! CMM preconditions the coupled strict equations. The returned wall flux
  ! is the algebraic value needed by the very same cell-mean update.
  subroutine advance_weak_diffusion(cache,comm,psi_old,source,dt,viscosity,psi,omega,wall_flux,iterations,error,parent_torque)
    type(weak_diffusion_cache), intent(inout) :: cache
    integer, intent(in) :: comm
    real(real_8), intent(in) :: psi_old(:,:,:),source(:,:,:),dt,viscosity
    real(real_8), intent(out) :: psi(:,:,:),omega(:,:,:),wall_flux(:,:,:)
    integer, intent(out) :: iterations
    character(*), intent(out) :: error
    real(real_8), optional, intent(in) :: parent_torque(0:)
    integer, parameter :: krylov=24,max_restart=8
    real(real_8), allocatable :: old_omega(:,:,:),old_diff(:,:,:),lap(:,:,:),diff(:,:,:),rhs(:,:,:),target(:,:,:), &
                                 vector(:,:,:),mapped(:,:,:),basis(:,:,:,:)
    real(real_8), allocatable :: prior(:,:,:),left(:,:,:),trace(:),torque(:),boundary_rhs(:,:), &
                                 trace_local(:,:),trace_extrema(:,:)
    real(real_8) :: source_residual(1)
    real(real_8) :: hessenberg(krylov+1,krylov),c(krylov),s(krylov),g(krylov+1),coeff(krylov), &
                    theta,settings(2),root(2),size_v,rho,tmp,tolerance,volume
    integer :: ni,nj,k,l,r,i,j,ierr,invalid,any_invalid,comparison,nk,info,component
    logical :: converged

    error = ''
    iterations = 0
    invalid = 0
    if (associated(cache%geom)) then
      ni = cache%geom%ni
      nj = cache%geom%nj
      nk = ubound(cache%geom%global%loops,1)
      if (present(parent_torque)) then
        if (size(parent_torque) /= nk+1) invalid = 1
      else if (nk > 0) then
        invalid = 1
      end if
      call mpi_comm_compare(comm,cache%geom%comm,comparison,ierr)
      if (comparison == MPI_UNEQUAL .or. any(shape(psi_old) /= [ni,nj,1]) .or. &
          any(shape(source) /= [ni,nj,1]) .or. any(shape(psi) /= [ni,nj,1]) .or. &
          any(shape(omega) /= [ni,nj,1]) .or. any(shape(wall_flux) /= [ni,nj,1])) invalid = 1
    else
      invalid = 1
    end if
    if (.not. all(ieee_is_finite([dt,viscosity])) .or. min(dt,viscosity)<=0.0D0) invalid = 1
    call mpi_allreduce(invalid,any_invalid,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (any_invalid /= 0) then
      error = 'invalid coupled weak diffusion cache, shape or time step'
      return
    end if
    settings = [dt,viscosity]
    root = settings
    call mpi_bcast(root,2,MPI_DOUBLE_PRECISION,0,comm,ierr)
    if (any(root/=settings)) error = 'weak diffusion settings differ across ranks'
    associate(geom => cache%geom)
      allocate(trace(0:nk),torque(0:nk),trace_local(2,0:nk),trace_extrema(2,0:nk),boundary_rhs(nk,1))
      trace_local = -huge(0.0D0)
      do j = 1,nj
        do i = 1,ni
          k = geom%global%kind(geom%i0+i,geom%j0+j)
          if (geom%area(i,j)>0.0D0 .and. .not. ieee_is_finite(source(i,j,1))) error = 'nonfinite wet-cell source'
          if (k==fluid_node .and. .not. ieee_is_finite(psi_old(i,j,1))) error = 'nonfinite strict streamfunction'
          component = geom%global%component(geom%i0+i,geom%j0+j)
          if (component > 0) then
            volume = psi_old(i,j,1)
            if (.not. ieee_is_finite(volume)) error = 'nonfinite old island trace'
            trace_local(:,component) = max(trace_local(:,component),[volume,-volume])
          end if
        end do
      end do
      call broadcast_error(geom%comm,geom%rank,geom%ranks,error)
      if (len_trim(error) /= 0) return
      call mpi_allreduce(trace_local,trace_extrema,size(trace_local),MPI_DOUBLE_PRECISION,MPI_MAX,geom%comm,ierr)
      if (any(trace_extrema(1,1:) /= -trace_extrema(2,1:))) then
        error = 'old streamfunction is not constant on an island'
        return
      end if
      if (nk > 0) then
        torque = parent_torque
        call source_momentum_budget(geom,source,reshape(torque,[nk+1,1]),source_residual,error)
        if (len_trim(error) /= 0) return
      end if
      trace = 0.0D0
      trace(1:) = trace_extrema(1,1:)
      allocate(old_omega(ni,nj,1),old_diff(ni,nj,1),lap(ni,nj,1),diff(ni,nj,1),rhs(ni,nj,1),target(ni,nj,1), &
               vector(ni,nj,1),mapped(ni,nj,1),basis(ni,nj,1,krylov+1),prior(ni,nj,1),left(ni,nj,1))
      theta = dt*viscosity/2
      call impose_trace(geom,psi_old,trace,prior)
      call evaluate_diffusion(cache,prior,lap,old_omega,old_diff,error)
      if (len_trim(error) /= 0) return
      rhs = lap+theta*old_diff
      do j = 1,nj
        do i = 1,ni
          if (geom%area(i,j)>0.0D0) rhs(i,j,1) = rhs(i,j,1)+dt*source(i,j,1)
        end do
      end do
      call inverse_constraint(rhs,target,include_torque=.true.)
      if (len_trim(error) /= 0) return
      psi = prior
      tolerance = 2.0D-11*max(tiny(theta),sqrt(inner_product(target,target)))
      converged = .false.
      do r = 1,max_restart
        call apply_operator(psi,vector)
        if (len_trim(error) /= 0) return
        vector = target-vector
        rho = sqrt(inner_product(vector,vector))
        if (rho <= tolerance) then
          converged = .true.
          exit
        end if
        basis(:,:,:,1) = vector/rho
        hessenberg = 0.0D0
        g = 0.0D0
        g(1) = rho
        do k = 1,krylov
          call apply_operator(basis(:,:,:,k),vector)
          if (len_trim(error) /= 0) return
          do l = 1,k
            hessenberg(l,k) = inner_product(vector,basis(:,:,:,l))
            vector = vector-hessenberg(l,k)*basis(:,:,:,l)
          end do
          size_v = sqrt(inner_product(vector,vector))
          hessenberg(k+1,k) = size_v
          do l = 1,k-1
            tmp = c(l)*hessenberg(l,k)+s(l)*hessenberg(l+1,k)
            hessenberg(l+1,k) = -s(l)*hessenberg(l,k)+c(l)*hessenberg(l+1,k)
            hessenberg(l,k) = tmp
          end do
          tmp = hypot(hessenberg(k,k),size_v)
          if (tmp <= tiny(tmp)) then
            error = 'coupled weak diffusion Krylov basis broke down'
            return
          end if
          c(k) = hessenberg(k,k)/tmp
          s(k) = size_v/tmp
          hessenberg(k,k) = tmp
          g(k+1) = -s(k)*g(k)
          g(k) = c(k)*g(k)
          iterations = iterations+1
          if (abs(g(k+1))<=tolerance .or. size_v<=tiny(size_v)) exit
          basis(:,:,:,k+1) = vector/size_v
        end do
        k = min(k,krylov)
        do l = k,1,-1
          coeff(l) = (g(l)-dot_product(hessenberg(l,l+1:k),coeff(l+1:k)))/hessenberg(l,l)
          psi = psi+coeff(l)*basis(:,:,:,l)
        end do
      end do
      if (.not. converged) then
        call apply_operator(psi,vector)
        if (len_trim(error) /= 0) return
        vector = target-vector
        if (sqrt(inner_product(vector,vector)) > tolerance) then
          error = 'coupled weak diffusion did not reach the Krylov tolerance'
          return
        end if
      end if

      call evaluate_diffusion(cache,psi,lap,omega,diff,error)
      if (len_trim(error) /= 0) return
      wall_flux = 0.0D0
      do j = 1,nj
        do i = 1,ni
          if (geom%global%kind(geom%i0+i,geom%j0+j) /= wall_node) cycle
          volume = geom%area(i,j)/geom%h**2
          wall_flux(i,j,1) = (volume*(omega(i,j,1)-old_omega(i,j,1))/dt-source(i,j,1))/viscosity &
                             -0.5D0*(old_diff(i,j,1)+diff(i,j,1))
        end do
      end do
      if (nk > 0) then
        call island_rows(geom,wall_flux,boundary_rhs(:,1))
        if (any(abs(viscosity*boundary_rhs(:,1)+torque(1:)) > &
            1.0D-9*max(1.0D0,abs(torque(1:)),abs(viscosity*boundary_rhs(:,1))))) &
          error = 'no-slip component torque failed the coupled update'
      end if
    end associate

  contains

    subroutine inverse_constraint(pv,field,include_torque)
      real(real_8), intent(in) :: pv(:,:,:)
      real(real_8), intent(out) :: field(:,:,:)
      logical, intent(in) :: include_torque
      integer :: island

      trace = 0.0D0
      call mean_to_pv_nodes(cache%geom,pv,left)
      call cmm_solve(cache%cmm,left,reshape(trace,[nk+1,1]),mapped,error)
      if (len_trim(error) /= 0) return
      call impose_trace(cache%geom,mapped,trace,field)
      if (nk == 0) return
      call evaluate_stiffness(cache,field,lap,error)
      if (len_trim(error) /= 0) return
      call island_rows(cache%geom,lap-pv,boundary_rhs(:,1))
      if (include_torque) boundary_rhs(:,1) = boundary_rhs(:,1)+dt*torque(1:)
      if (cache%geom%rank == 0) then
        call dpotrs('L',nk,1,cache%factor,nk,boundary_rhs,nk,info)
        if (info /= 0) error = 'no-slip island stiffness factor solve failed'
      end if
      call mpi_bcast(error,len(error),MPI_CHARACTER,0,cache%geom%comm,ierr)
      if (len_trim(error) /= 0) return
      call mpi_bcast(boundary_rhs,nk,MPI_DOUBLE_PRECISION,0,cache%geom%comm,ierr)
      do island = 1,nk
        field = field+boundary_rhs(island,1)*cache%harmonic(:,:,:,island)
      end do
    end subroutine inverse_constraint

    subroutine apply_operator(field,result)
      real(real_8), intent(in) :: field(:,:,:)
      real(real_8), intent(out) :: result(:,:,:)

      call evaluate_diffusion(cache,field,lap,omega,diff,error)
      if (len_trim(error) /= 0) return
      rhs = theta*diff
      call inverse_constraint(rhs,result,include_torque=.false.)
      if (len_trim(error) == 0) result = field-result
    end subroutine apply_operator

    function inner_product(first,second) result(product)
      real(real_8), intent(in) :: first(:,:,:),second(:,:,:)
      real(real_8) :: product,local
      integer :: status

      local = sum(first*second)
      call mpi_allreduce(local,product,1,MPI_DOUBLE_PRECISION,MPI_SUM,cache%geom%comm,status)
    end function inner_product

  end subroutine advance_weak_diffusion

  subroutine release_weak_diffusion(cache)
    type(weak_diffusion_cache), intent(inout) :: cache

    if (.not. associated(cache%geom)) return
    call release_cmm(cache%cmm)
    call deallocate(cache%halo)
    call deallocate(cache%requests)
    deallocate(cache%requests,cache%flux,cache%harmonic,cache%factor)
    nullify(cache%geom)
  end subroutine release_weak_diffusion

end module qg_weak_diffusion
