! Viscous source and integrated Navier-slip closure for a single axial coastline.
module qg_cape_slip
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use mpi, only: mpi_allreduce, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_MAX
  use qg_parameters, only: real_8
  use qg_geometry, only: solid_node, wall_node
  use qg_irregular_parallel, only: distributed_geometry, half_coordinates
  use qg_irregular_operators, only: wall_mean_plan, setup_wall_means, reconstruct_wall_means, &
                                   stiffness_flux, sum_open_flux, mean_to_pv_nodes
  use qg_parallel, only: decomposition, parallel_data, request_cache, allocate, deallocate, update, flush
  use qg_diagnostics, only: stratification
  use qg_constraints, only: basin_mass_cache, invert_basin_mass
  implicit none
  private
  public :: cape_slip_cache, setup_cape_slip, release_cape_slip, cape_viscous_state, cape_implicit_step

  type cape_slip_cache
    private
    type(distributed_geometry), pointer :: geom
    type(stratification), pointer :: strat
    type(basin_mass_cache), pointer :: inverse
    type(wall_mean_plan) :: wall_mean
    type(parallel_data) :: halo
    type(request_cache), pointer :: requests
    real(real_8), allocatable :: slip_weight(:,:), flux(:,:), omega(:,:,:), mean_psi(:,:,:)
    real(real_8) :: viscosity
  end type cape_slip_cache

contains

  subroutine setup_cape_slip(cache,decomp,geom,strat,inverse,viscosity,slip_length,error)
    type(cape_slip_cache), intent(out) :: cache
    type(decomposition), target, intent(in) :: decomp
    type(distributed_geometry), target, intent(in) :: geom
    type(stratification), target, intent(in) :: strat
    type(basin_mass_cache), target, intent(inout) :: inverse
    real(real_8), intent(in) :: viscosity,slip_length ! native units
    character(*), intent(out) :: error
    integer :: k,i,j,edge,previous,ni,point(2),incoming(2),outgoing(2)

    error=''
    if (size(geom%global%loops)/=1 .or. .not. all(ieee_is_finite([viscosity,slip_length])) .or. &
        min(viscosity,slip_length)<=0.0D0) then
      error='cape partial slip requires one coastline and positive finite viscosity and slip length'
      return
    end if
    associate(loop=>geom%global%loops(0)%ij)
      do edge=1,size(loop,2)
        previous=modulo(edge-2,size(loop,2))+1
        if (sum(abs(loop(:,edge)-loop(:,previous)))/=1) then
          error='cape partial slip requires an axial coastline'
          return
        end if
      end do
    end associate
    call setup_wall_means(geom,cache%wall_mean,error)
    if (len_trim(error)/=0) return
    cache%geom=>geom
    cache%strat=>strat
    cache%inverse=>inverse
    cache%viscosity=viscosity
    allocate(cache%requests)
    call allocate(cache%requests,comm=geom%comm)
    call allocate(cache%halo,decomp,requests=cache%requests)
    allocate(cache%flux(strat%layers,size(geom%faces)),cache%omega(geom%ni,geom%nj,strat%layers), &
             cache%mean_psi(geom%ni,geom%nj,strat%layers),cache%slip_weight(geom%ni,geom%nj))
    cache%slip_weight=0.0D0
    ni=size(geom%global%kind,1)
    do k=1,size(geom%walls)
      edge=geom%walls(k)%edge
      associate(loop=>geom%global%loops(0)%ij)
        previous=modulo(edge-2,size(loop,2))+1
        i=loop(1,edge)-geom%i0
        j=loop(2,edge)-geom%j0
        point=geom%walls(k)%endpoints(:,1)
        incoming=point-half_coordinates(loop(1,previous),loop(2,previous),ni)
        outgoing=geom%walls(k)%endpoints(:,2)-point
      end associate
      ! Half of each adjacent wall segment, expressed in native grid units.
      cache%slip_weight(i,j)=slip_length*(0.25D0*sum(abs(incoming)+abs(outgoing)))
    end do
  end subroutine setup_cape_slip

  ! The strict cell values are the independent PV state. At the wall,
  ! Gamma = A*omega - B*psi = -ell*L*omega, with a lumped wall trace of
  ! length L. Use the same symmetric open-face stiffness B on omega:
  ! the viscous contribution to energy is -nu * (B*psi)'*omega.
  subroutine cape_viscous_state(cache,q,psi,wall,rate,slip_defect,error)
    type(cape_slip_cache), intent(inout) :: cache
    real(real_8), intent(inout) :: q(:,:,:)
    real(real_8), intent(out) :: psi(:,:,:),wall(0:,:),rate(:,:,:),slip_defect
    character(*), intent(out) :: error
    real(real_8) :: target(size(q,3)-1),area,weight,defect
    integer :: i,j,kind,ierr

    associate(g=>cache%geom,s=>cache%strat)
      target=0.0D0
      call invert_basin_mass(cache%inverse,g%comm,q,target,psi,wall(0,:),error)
      if (len_trim(error)/=0) return
      call reconstruct_wall_means(cache%wall_mean,psi,wall,cache%mean_psi)
      where (spread(g%area==0.0D0,3,s%layers)) psi=0.0D0
      cache%halo%sdata=psi
      call update(cache%halo)
      call flush(cache%halo)
      call stiffness_flux(g,cache%halo%data,cache%flux,error)
      if (len_trim(error)/=0) return
      call sum_open_flux(g,cache%flux,rate)

      defect=0.0D0
      cache%omega=0.0D0
      do j=1,g%nj
        do i=1,g%ni
          kind=g%global%kind(g%i0+i,g%j0+j)
          select case(kind)
            case(solid_node)
              q(i,j,:)=0.0D0
            case(wall_node)
              area=g%area(i,j)/g%h**2
              weight=cache%slip_weight(i,j)
              cache%omega(i,j,:)=rate(i,j,:)/(area+weight)
              q(i,j,:)=cache%omega(i,j,:)+matmul(s%s_mat,cache%mean_psi(i,j,:))
              defect=max(defect,maxval(abs((area+weight)*cache%omega(i,j,:)-rate(i,j,:))))
            case default
              cache%omega(i,j,:)=q(i,j,:)-matmul(s%s_mat,psi(i,j,:))
          end select
        end do
      end do
      call mpi_allreduce(defect,slip_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,g%comm,ierr)

      cache%halo%sdata=cache%omega
      call update(cache%halo)
      call flush(cache%halo)
      call stiffness_flux(g,cache%halo%data,cache%flux,error)
      if (len_trim(error)/=0) return
      call sum_open_flux(g,cache%flux,rate)
      do j=1,g%nj
        do i=1,g%ni
          area=g%area(i,j)/g%h**2
          if (area>0.0D0) rate(i,j,:)=cache%viscosity*rate(i,j,:)/area
        end do
      end do
    end associate
  end subroutine cape_viscous_state

  ! Solve (I-step*D) on strict PV with restarted GMRES. The viscosity
  ! operator D includes the modal inverse and the algebraic wall solve.
  subroutine cape_implicit_step(cache,step,source,q,error,relative_tolerance)
    type(cape_slip_cache), intent(inout) :: cache
    real(real_8), intent(in) :: step,source(:,:,:)
    real(real_8), intent(out) :: q(:,:,:)
    character(*), intent(out) :: error
    real(real_8), optional, intent(in) :: relative_tolerance
    integer, parameter :: krylov=24,restarts=8
    integer :: ni,nj,nl,r,k,l
    real(real_8), allocatable :: basis(:,:,:,:), rhs(:,:,:), residual(:,:,:), work(:,:,:), psi(:,:,:), rate(:,:,:)
    real(real_8) :: h(krylov+1,krylov),cs(krylov),sn(krylov),g(krylov+1),coeff(krylov)
    real(real_8) :: wall(0:0,size(q,3)),defect,tolerance,rho,length,tmp

    error=''
    ni=size(q,1)
    nj=size(q,2)
    nl=size(q,3)
    allocate(rhs(ni,nj,nl),residual(ni,nj,nl),work(ni,nj,nl),psi(ni,nj,nl),rate(ni,nj,nl))
    call mean_to_pv_nodes(cache%geom,source,rhs)
    q=rhs
    rho=sqrt(inner_product(rhs,rhs))
    if (rho==0.0D0) return
    allocate(basis(ni,nj,nl,krylov+1))
    tolerance=1.0D-12*rho
    if (present(relative_tolerance)) tolerance=relative_tolerance*rho
    do r=1,restarts
      call apply(q,residual)
      if (len_trim(error)/=0) return
      residual=rhs-residual
      rho=sqrt(inner_product(residual,residual))
      if (rho<=tolerance) return
      basis(:,:,:,1)=residual/rho
      h=0.0D0
      g=0.0D0
      g(1)=rho
      do k=1,krylov
        call apply(basis(:,:,:,k),residual)
        if (len_trim(error)/=0) return
        do l=1,k
          h(l,k)=inner_product(residual,basis(:,:,:,l))
          residual=residual-h(l,k)*basis(:,:,:,l)
        end do
        length=sqrt(inner_product(residual,residual))
        do l=1,k-1
          tmp=cs(l)*h(l,k)+sn(l)*h(l+1,k)
          h(l+1,k)=-sn(l)*h(l,k)+cs(l)*h(l+1,k)
          h(l,k)=tmp
        end do
        tmp=hypot(h(k,k),length)
        if (tmp<=tiny(tmp)) then
          error='cape implicit viscosity Krylov basis broke down'
          return
        end if
        cs(k)=h(k,k)/tmp
        sn(k)=length/tmp
        h(k,k)=tmp
        g(k+1)=-sn(k)*g(k)
        g(k)=cs(k)*g(k)
        if (abs(g(k+1))<=tolerance .or. length<=tiny(length)) exit
        basis(:,:,:,k+1)=residual/length
      end do
      k=min(k,krylov)
      do l=k,1,-1
        coeff(l)=(g(l)-dot_product(h(l,l+1:k),coeff(l+1:k)))/h(l,l)
        q=q+coeff(l)*basis(:,:,:,l)
      end do
    end do
    call apply(q,residual)
    if (len_trim(error)/=0) return
    residual=rhs-residual
    if (sqrt(inner_product(residual,residual))>tolerance) error='cape implicit viscosity did not converge'

  contains
    subroutine apply(field,result)
      real(real_8), intent(in) :: field(:,:,:)
      real(real_8), intent(out) :: result(:,:,:)

      work=field
      call cape_viscous_state(cache,work,psi,wall,rate,defect,error)
      if (len_trim(error)/=0) return
      call mean_to_pv_nodes(cache%geom,rate,result)
      result=field-step*result
    end subroutine apply

    function inner_product(a,b) result(value)
      real(real_8), intent(in) :: a(:,:,:),b(:,:,:)
      real(real_8) :: value,local
      integer :: ierr

      local=sum(a*b)
      call mpi_allreduce(local,value,1,MPI_DOUBLE_PRECISION,MPI_SUM,cache%geom%comm,ierr)
    end function inner_product
  end subroutine cape_implicit_step

  subroutine release_cape_slip(cache)
    type(cape_slip_cache), intent(inout) :: cache

    call deallocate(cache%halo)
    call deallocate(cache%requests)
    deallocate(cache%requests)
  end subroutine release_cape_slip

end module qg_cape_slip
