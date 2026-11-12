! Discrete incidence, weak Poisson circulation and viscous flux balance.
module qg_irregular_operators
  use qg_geometry, only: fluid_node, wall_node
  use qg_irregular_parallel, only: distributed_geometry, share_face_flux, broadcast_error, half_coordinates
  use qg_moments, only: n_moments, x_power, y_power, monomial_index, segment_moments
  use qg_reconstruction, only: reconstruction_plan, cubic_terms, build_reconstruction, evaluate_reconstruction
  use qg_parameters, only: real_8
  implicit none
  include 'mpif.h'
  private
  public :: mean_to_pv_nodes, stiffness_flux, sum_open_flux, weak_circulation, diffusion_balance
  public :: viscous_face_plan, setup_viscous_faces, viscous_open_flux
  public :: wall_mean_plan, setup_wall_means, reconstruct_wall_means
  public :: wall_trace_plan, setup_wall_trace, reconstruct_wall_trace

  type viscous_face_plan
    private
    type(reconstruction_plan) :: reconstruction
    integer, allocatable :: owned_faces(:)
  end type viscous_face_plan

  type wall_mean_plan
    private
    type(reconstruction_plan) :: reconstruction
    integer, allocatable :: ij(:,:)
  end type wall_mean_plan

  type wall_trace_plan
    private
    type(reconstruction_plan) :: reconstruction
    integer, allocatable :: ij(:,:)
  end type wall_trace_plan

contains

  ! Fit a cubic to visible cell means and integrate its normal derivative
  ! over the actual shared segment. Only the north/west owner builds a row.
  subroutine setup_viscous_faces(geom, plan, error)
    type(distributed_geometry), intent(in) :: geom
    type(viscous_face_plan), intent(out) :: plan
    character(*), intent(out) :: error
    integer, allocatable :: origins(:,:)
    real(real_8), allocatable :: target(:,:)
    integer :: n,f,k,m,a,b
    real(real_8) :: half_segment(2), normal(2), moment(n_moments)

    n = count(geom%faces%owner == geom%rank)
    allocate(plan%owned_faces(n),origins(2,n),target(cubic_terms,n))
    target = 0.0D0
    k = 0
    do f = 1,size(geom%faces)
      associate(face => geom%faces(f))
        if (face%owner /= geom%rank) cycle
        k = k+1
        plan%owned_faces(k) = f
        origins(:,k) = sum(face%endpoints,dim=2)
        half_segment = 0.25D0*(face%endpoints(:,2)-face%endpoints(:,1))
      end associate
      moment = segment_moments(-half_segment,half_segment)
      normal = 2.0D0*[half_segment(2),-half_segment(1)]
      do m = 2,cubic_terms
        a = x_power(m)
        b = y_power(m)
        if (a > 0) target(m,k) = a*normal(1)*moment(monomial_index(a-1,b))
        if (b > 0) target(m,k) = target(m,k)+b*normal(2)*moment(monomial_index(a,b-1))
      end do
    end do
    call build_reconstruction(geom,origins,target,.true.,plan%reconstruction,error)
  end subroutine setup_viscous_faces

  subroutine viscous_open_flux(geom, plan, omega_mean, flux, error)
    type(distributed_geometry), intent(in) :: geom
    type(viscous_face_plan), intent(in) :: plan
    real(real_8), intent(in) :: omega_mean(:,:,:)
    real(real_8), intent(out) :: flux(:,:)
    character(*), intent(out) :: error
    real(real_8), allocatable :: values(:,:)
    integer :: layers

    error = ''
    layers = size(omega_mean,3)
    if (any(shape(omega_mean) /= [geom%ni,geom%nj,layers]) .or. any(shape(flux) /= [layers,size(geom%faces)])) &
      error = 'viscous face arrays have incompatible shapes'
    call broadcast_error(geom%comm,geom%rank,geom%ranks,error)
    if (len_trim(error) /= 0) return
    allocate(values(layers,size(plan%owned_faces)))
    call evaluate_reconstruction(plan%reconstruction,omega_mean,values)
    flux = 0.0D0
    flux(:,plan%owned_faces) = values
    call share_face_flux(geom,flux,error)
  end subroutine viscous_open_flux

  ! Map owned streamfunction point data to means only on cut cells.
  ! The wall point is prescribed, while the average is evaluated by moments.
  subroutine setup_wall_means(geom, plan, error)
    type(distributed_geometry), intent(in) :: geom
    type(wall_mean_plan), intent(out) :: plan
    character(*), intent(out) :: error
    integer, allocatable :: origin(:,:)
    real(real_8), allocatable :: target(:,:)
    integer :: n,k,i,j,point(2)

    n = size(geom%walls)
    allocate(plan%ij(2,n),origin(2,n),target(cubic_terms,n))
    do k = 1,n
      point = geom%walls(k)%endpoints(:,1)
      i = size(geom%owner,1)-point(2)/2-geom%i0
      j = 1+point(1)/2-geom%j0
      plan%ij(:,k) = [i,j]
      origin(:,k) = 2*point
      target(:,k) = geom%moments(:cubic_terms,i,j)
    end do
    call build_reconstruction(geom,origin,target,.false.,plan%reconstruction,error)
  end subroutine setup_wall_means

  subroutine reconstruct_wall_means(plan, psi, wall, cut_mean)
    type(wall_mean_plan), intent(in) :: plan
    real(real_8), intent(in) :: psi(:,:,:),wall(0:,:)
    real(real_8), intent(out) :: cut_mean(:,:,:)
    real(real_8) :: values(size(psi,3),size(plan%ij,2))
    integer :: k

    cut_mean = 0.0D0
    call evaluate_reconstruction(plan%reconstruction,psi,values,wall)
    do k = 1,size(plan%ij,2)
      cut_mean(plan%ij(1,k),plan%ij(2,k),:) = values(:,k)
    end do
  end subroutine reconstruct_wall_means

  ! For finite slip use the integral of omega along the two actual wall
  ! half-segments. The fit consumes cell means and is a smooth functional.
  ! Its singular trace convergence requires a coupled discrete test.
  subroutine setup_wall_trace(geom, plan, error)
    type(distributed_geometry), intent(in) :: geom
    type(wall_trace_plan), intent(out) :: plan
    character(*), intent(out) :: error
    integer, allocatable :: origins(:,:)
    real(real_8), allocatable :: target(:,:)
    integer :: n,k,gi,edge,prev,start(2),finish(2),incoming(2),a(2,2),b(2,2),side,row,component
    real(real_8) :: half_step(2),moment(n_moments)

    n = size(geom%walls)
    gi = size(geom%owner,1)
    allocate(plan%ij(2,n),origins(2,2*n),target(cubic_terms,2*n))
    do k = 1,n
      start = geom%walls(k)%endpoints(:,1)
      finish = geom%walls(k)%endpoints(:,2)
      edge = geom%walls(k)%edge
      component = geom%walls(k)%component
      associate(loop => geom%global%loops(component)%ij)
        prev = modulo(edge-2,size(loop,2))+1
        incoming = half_coordinates(loop(1,prev),loop(2,prev),gi)
        plan%ij(:,k) = loop(:,edge)-[geom%i0,geom%j0]
      end associate
      a(:,1) = (incoming+start)/2
      b(:,1) = start
      a(:,2) = start
      b(:,2) = (start+finish)/2
      do side = 1,2
        row = 2*(k-1)+side
        origins(:,row) = a(:,side)+b(:,side)
        half_step = 0.25D0*(b(:,side)-a(:,side))
        moment = segment_moments(-half_step,half_step)
        target(:,row) = 2*norm2(half_step)*moment(:cubic_terms)
      end do
    end do
    call build_reconstruction(geom,origins,target,.true.,plan%reconstruction,error)
  end subroutine setup_wall_trace

  subroutine reconstruct_wall_trace(plan, omega_mean, integral)
    type(wall_trace_plan), intent(in) :: plan
    real(real_8), intent(in) :: omega_mean(:,:,:)
    real(real_8), intent(out) :: integral(:,:,:)
    real(real_8) :: values(size(omega_mean,3),2*size(plan%ij,2))
    integer :: k

    integral = 0.0D0
    call evaluate_reconstruction(plan%reconstruction,omega_mean,values)
    do k = 1,size(plan%ij,2)
      integral(plan%ij(1,k),plan%ij(2,k),:) = values(:,2*k-1)+values(:,2*k)
    end do
  end subroutine reconstruct_wall_trace

  ! Strict nodes have a full centred square. The mean-to-point identity is
  ! second order there. Boundary-cell means do not enter Dirichlet rows.
  subroutine mean_to_pv_nodes(geom, cell_mean, nodal)
    type(distributed_geometry), intent(in) :: geom
    real(real_8), intent(in) :: cell_mean(:,:,:)
    real(real_8), intent(out) :: nodal(:,:,:)
    integer :: i,j

    nodal = 0.0D0
    do j = 1,geom%nj
      do i = 1,geom%ni
        if (geom%global%kind(geom%i0+i,geom%j0+j) == fluid_node) nodal(i,j,:) = cell_mean(i,j,:)
      end do
    end do
  end subroutine mean_to_pv_nodes

  ! The right-triangle P1 stiffness is the axial difference times the open
  ! face length in native grid units. This is exact on strict Poisson rows.
  ! Input contains one populated halo ring, including physical wall values.
  subroutine stiffness_flux(geom, nodal, flux, error)
    type(distributed_geometry), intent(in) :: geom
    real(real_8), intent(in) :: nodal(:,:,:)
    real(real_8), intent(out) :: flux(:,:)
    character(*), intent(out) :: error
    integer :: f, i(2), j(2), layers
    real(real_8) :: length

    error = ''
    layers = size(nodal,3)
    if (any(shape(nodal) /= [geom%ni+2,geom%nj+2,layers]) .or. any(shape(flux) /= [layers,size(geom%faces)])) &
      error = 'stiffness flux has incompatible array shape'
    call broadcast_error(geom%comm,geom%rank,geom%ranks,error)
    if (len_trim(error) /= 0) return

    flux = 0.0D0
    do f = 1,size(geom%faces)
      associate(face => geom%faces(f))
        if (face%owner /= geom%rank) cycle
        i = face%cells(1,:)-geom%i0+1
        j = face%cells(2,:)-geom%j0+1
        length = 0.5D0*sum(abs(face%endpoints(:,2)-face%endpoints(:,1)))
        flux(:,f) = length*(nodal(i(2),j(2),:)-nodal(i(1),j(1),:))
      end associate
    end do
    call share_face_flux(geom,flux,error)
  end subroutine stiffness_flux

  ! Flux is already owner-synchronized. Return integrated open divergence,
  ! without area division. A positive value points out of the cell.
  subroutine sum_open_flux(geom, flux, integrated)
    type(distributed_geometry), intent(in) :: geom
    real(real_8), intent(in) :: flux(:,:)
    real(real_8), intent(out) :: integrated(:,:,:)
    integer :: i,j,side,f,sign

    integrated = 0.0D0
    do j = 1,geom%nj
      do i = 1,geom%ni
        do side = 1,4
          f = geom%face_index(side,i,j)
          if (f == 0) cycle
          sign = merge(1,-1,side<=2)
          integrated(i,j,:) = integrated(i,j,:)+sign*flux(:,f)
        end do
      end do
    end do
  end subroutine sum_open_flux

  ! Discrete Poisson remainder. At each wall cell this is the weak wall
  ! normal derivative integral. Strict-cell remainders are reported for
  ! consistency checks. omega_mean must use the stage's stretching map.
  subroutine weak_circulation(geom, omega_mean, psi_flux, circulation, strict_defect, remainder, error)
    type(distributed_geometry), intent(in) :: geom
    real(real_8), intent(in) :: omega_mean(:,:,:), psi_flux(:,:)
    real(real_8), intent(out) :: circulation(0:,:), strict_defect, remainder(:,:,:)
    character(*), intent(out) :: error
    real(real_8) :: local(0:ubound(geom%global%loops,1),size(omega_mean,3)), local_defect
    integer :: i,j,component,kind,ierr

    call validate_balance_arrays(geom,omega_mean,psi_flux,shape(circulation),shape(remainder),error)
    if (len_trim(error) /= 0) return
    call sum_open_flux(geom,psi_flux,remainder)
    local = 0.0D0
    local_defect = 0.0D0
    do j = 1,geom%nj
      do i = 1,geom%ni
        if (geom%area(i,j) == 0.0D0) cycle
        remainder(i,j,:) = (geom%area(i,j)/geom%h**2)*omega_mean(i,j,:)-remainder(i,j,:)
        kind = geom%global%kind(geom%i0+i,geom%j0+j)
        if (kind == wall_node) then
          component = geom%global%component(geom%i0+i,geom%j0+j)
          local(component,:) = local(component,:)+remainder(i,j,:)
        else
          local_defect = max(local_defect,maxval(abs(remainder(i,j,:))))
        end if
      end do
    end do
    call mpi_allreduce(local,circulation,size(local),MPI_DOUBLE_PRECISION,MPI_SUM,geom%comm,ierr)
    call mpi_allreduce(local_defect,strict_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,geom%comm,ierr)
  end subroutine weak_circulation

  ! The same outward wall fluxes determine the source and the component
  ! torque. Multiply both by the layer viscosity at the source stage.
  subroutine diffusion_balance(geom, open_flux, wall_flux, laplacian, torque, error)
    type(distributed_geometry), intent(in) :: geom
    real(real_8), intent(in) :: open_flux(:,:), wall_flux(:,:,:)
    real(real_8), intent(out) :: laplacian(:,:,:), torque(0:,:)
    character(*), intent(out) :: error
    real(real_8) :: local(0:ubound(geom%global%loops,1),size(wall_flux,3))
    integer :: i,j,component,ierr

    call validate_balance_arrays(geom,wall_flux,open_flux,shape(torque),shape(laplacian),error)
    if (len_trim(error) /= 0) return
    call sum_open_flux(geom,open_flux,laplacian)
    local = 0.0D0
    do j = 1,geom%nj
      do i = 1,geom%ni
        if (geom%area(i,j) == 0.0D0) cycle
        component = geom%global%component(geom%i0+i,geom%j0+j)
        if (component >= 0) then
          laplacian(i,j,:) = laplacian(i,j,:)+wall_flux(i,j,:)
          local(component,:) = local(component,:)+wall_flux(i,j,:)
        end if
        laplacian(i,j,:) = laplacian(i,j,:)/(geom%area(i,j)/geom%h**2)
      end do
    end do
    call mpi_allreduce(local,torque,size(local),MPI_DOUBLE_PRECISION,MPI_SUM,geom%comm,ierr)
  end subroutine diffusion_balance

  subroutine validate_balance_arrays(geom, cell, face, contour_shape, result_shape, error)
    type(distributed_geometry), intent(in) :: geom
    real(real_8), intent(in) :: cell(:,:,:), face(:,:)
    integer, intent(in) :: contour_shape(2), result_shape(3)
    character(*), intent(out) :: error
    integer :: layers, counts(2), ierr

    error = ''
    layers = size(cell,3)
    if (any(shape(cell) /= [geom%ni,geom%nj,layers]) .or. any(result_shape /= shape(cell)) .or. &
        any(shape(face) /= [layers,size(geom%faces)]) .or. any(contour_shape /= [size(geom%global%loops),layers])) &
      error = 'irregular flux balance has incompatible array shape'
    call broadcast_error(geom%comm,geom%rank,geom%ranks,error)
    if (len_trim(error) /= 0) return
    call mpi_allreduce([layers,-layers],counts,2,MPI_INTEGER,MPI_MAX,geom%comm,ierr)
    if (counts(1) /= -counts(2)) error = 'irregular flux balance layer counts differ across ranks'
  end subroutine validate_balance_arrays

end module qg_irregular_operators
