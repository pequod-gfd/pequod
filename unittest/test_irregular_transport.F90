program test_irregular_transport
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
  use qg_irregular_parallel, only: distributed_geometry, load_distributed_geometry, gather_second_face_trace
  use qg_irregular_bc, only: vertex_plan, setup_vertex_plan, instantaneous_transport
  use qg_irregular_transport, only: transport_fluxes, wet_cell_rate, staggered_trace, affine_cabaret_trace, integrate_stage_rate
  use qg_irregular_transport, only: source_substep, source_momentum_budget, cut_cabaret_plan, &
                                    setup_cut_cabaret, cut_cabaret_traces
  use qg_parallel, only: parallel_init, parallel_finalise
  use qg_parameters, only: real_8
  implicit none
  include 'mpif.h'

  integer, parameter :: n=33,layers=2
  type(distributed_geometry) :: geom
  type(vertex_plan) :: vertices
  type(cut_cabaret_plan) :: cut_plan
  integer :: comm,rank,ranks,ierr,box(4),grid(2),row,column,i,j,f,l,k,side,cell,bad_rank,candidate,gi,gj,fallback
  real(real_8), allocatable :: psi(:,:,:),wall(:,:),velocity(:,:,:),trace(:,:,:,:),flux(:,:)
  real(real_8), allocatable :: torque(:,:)
  real(real_8), allocatable :: old_trace(:,:,:,:),cell_velocity(:,:,:,:),source_in(:,:,:),source_out(:,:,:)
  real(real_8), allocatable :: reference(:,:),outgoing(:,:,:),source(:,:,:),rate(:,:,:)
  real(real_8) :: point(2),half(2),exact,actual,volume,local(layers),global(layers),limits(0:2)
  real(real_8) :: z,uq,qq,peak,reference_peak,length,difference,nan,min_area,local_min
  real(real_8) :: t(4),integral,target,coarse,fine
  real(real_8) :: centroid(2),samples(2,3),times(3),sources(3),observations(3),vel(2),ends(2,2),values(2)
  real(real_8) :: gradient(2),expected(2),corrected(3),stage_dt,source_rate,face_error
  character(1024) :: path,error

  if (command_argument_count() /= 1) then
    print '(a)', 'Irregular transport fixture requires a .dat path'
    stop 2
  end if
  call get_command_argument(1,path)
  call parallel_init()
  call mpi_comm_rank(MPI_COMM_WORLD,rank,ierr)
  call mpi_comm_size(MPI_COMM_WORLD,ranks,ierr)
  call mpi_comm_split(MPI_COMM_WORLD,0,ranks-1-rank,comm,ierr)
  grid = 0
  call mpi_dims_create(ranks,2,grid,ierr)
  row = rank/grid(2)
  column = mod(rank,grid(2))
  box = [row*n/grid(1),column*n/grid(2),(row+1)*n/grid(1)-row*n/grid(1),(column+1)*n/grid(2)-column*n/grid(2)]
  call load_distributed_geometry(path,n,n,box,comm,2.5D0,geom,error)
  call require(len_trim(error)==0,error)
  f = size(geom%faces)
  allocate(velocity(2,layers,f),trace(2,2,layers,f),flux(layers,f),reference(layers,f),outgoing(2,layers,f))
  allocate(source(geom%ni,geom%nj,layers),rate(geom%ni,geom%nj,layers))
  allocate(psi(geom%ni,geom%nj,layers),wall(0:ubound(geom%global%loops,1),layers))
  allocate(torque(0:ubound(geom%global%loops,1),layers))

  ! Independent Gauss integration of sign-changing affine face velocities
  ! and distinct affine traces. Reverse signs on alternate faces and layers.
  nan = ieee_value(0.0D0,ieee_quiet_nan)
  velocity = nan
  trace = nan
  do f = 1,size(geom%faces)
    length = 0.5D0*sum(abs(geom%faces(f)%endpoints(:,2)-geom%faces(f)%endpoints(:,1)))
    limits = [0.0D0,0.25D0,1.0D0]
    do l = 1,layers
      exact = (-1.0D0)**(geom%faces(f)%id+l)
      if (geom%faces(f)%owner == geom%rank) then
        velocity(:,l,f) = exact*[-1.0D0,3.0D0]
        trace(:,:,l,f) = reshape([1.0D0,4.0D0,6.0D0,4.0D0],[2,2])
      end if
      reference(l,f) = 0.0D0
      outgoing(:,l,f) = 0.0D0
      do k = 1,2
        cell = merge(1,2,exact*(k-1.5D0) > 0)
        do i = -1,1,2
          z = 0.5D0*(limits(k-1)+limits(k))+i*(limits(k)-limits(k-1))/(2*sqrt(3.0D0))
          uq = exact*(4*z-1)
          qq = merge(1+3*z,6-2*z,cell==1)
          volume = 0.5D0*length*(limits(k)-limits(k-1))
          reference(l,f) = reference(l,f)+volume*uq*qq
          outgoing(cell,l,f) = outgoing(cell,l,f)+volume*abs(uq)
        end do
      end do
    end do
  end do
  call transport_fluxes(geom,velocity,trace,flux,peak,error)
  call require(len_trim(error)==0,error)
  exact = max(0.0D0,maxval(abs(flux-reference)))
  call mpi_allreduce(exact,difference,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
  call require(difference < 3.0D-15,'shared upwind flux disagrees with independent quadrature')
  reference_peak = 0.0D0
  do j = 1,geom%nj
    do i = 1,geom%ni
      if (geom%area(i,j) == 0.0D0) cycle
      local = 0.0D0
      do side = 1,4
        f = geom%face_index(side,i,j)
        if (f == 0) cycle
        cell = merge(1,2,side<=2)
        local = local+outgoing(cell,:,f)
      end do
      reference_peak = max(reference_peak,maxval(local)*geom%h**2/geom%area(i,j))
    end do
  end do
  call mpi_allreduce(reference_peak,actual,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
  call require(abs(peak-actual)<1.0D-13,'CFL omitted one side of a reversing face')
  reference_peak = actual

  ! Integrated finite-volume source and owner-face cancellation, including
  ! the smallest cut cells. Dry inputs are deliberately poisoned.
  source = nan
  local = 0.0D0
  do j = 1,geom%nj
    do i = 1,geom%ni
      if (geom%area(i,j) == 0.0D0) cycle
      source(i,j,:) = geom%area(i,j)/geom%h**2*[1.0D0,-2.0D0]
      local = local+source(i,j,:)
    end do
  end do
  call wet_cell_rate(geom,flux,source,rate,error)
  call require(len_trim(error)==0,error)
  do j = 1,geom%nj
    do i = 1,geom%ni
      local = local-geom%area(i,j)/geom%h**2*rate(i,j,:)
    end do
  end do
  call mpi_allreduce(local,global,layers,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  call require(maxval(abs(global))<2.0D-12,'finite-volume rate did not close the source budget')

  local = 0.0D0
  do j = 1,geom%nj
    do i = 1,geom%ni
      if (geom%area(i,j)>0) local = local+source(i,j,:)
    end do
  end do
  call mpi_allreduce(local,global,layers,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  do k = 1,ubound(torque,1)
    torque(k,:) = k*[0.125D0,-0.25D0]
  end do
  torque(0,:) = global-sum(torque(1:,:),dim=1)
  call source_momentum_budget(geom,source,torque,global,error)
  call require(len_trim(error)==0,error)
  call require(maxval(abs(global))<2.0D-12,'curl and torque did not share a global budget')
  torque(0,1) = torque(0,1)+0.125D0
  call source_momentum_budget(geom,source,torque,global,error)
  call require(error=='volume source and parent boundary torque are incompatible','incompatible curl and torque accepted')

  ! Construct divergence-free integrated volume from actual vertex values,
  ! then add a zero-mean slope along each face. A constant PV must be
  ! preserved even if upwinding switches within the face.
  call setup_vertex_plan(geom,vertices,error)
  call require(len_trim(error)==0,error)
  wall = 0.0D0
  do j = 1,geom%nj
    do i = 1,geom%ni
      point = [real(geom%j0+j-1,real_8),real(n-geom%i0-i,real_8)]
      psi(i,j,:) = sin(point(1)/6)*cos(point(2)/5)*[1.0D0,-2.0D0]
    end do
  end do
  call instantaneous_transport(geom,vertices,psi,wall,reference,error)
  call require(len_trim(error)==0,error)
  do f = 1,size(geom%faces)
    if (geom%faces(f)%owner /= geom%rank) cycle
    half = 0.5D0*(geom%faces(f)%endpoints(:,2)-geom%faces(f)%endpoints(:,1))
    do l = 1,layers
      velocity(:,l,f) = reference(l,f)/sum(abs(half))+[-0.3D0,0.3D0]
    end do
    trace(:,:,1,f) = 3.0D0
    trace(:,:,2,f) = -5.0D0
  end do
  call transport_fluxes(geom,velocity,trace,flux,peak,error)
  call require(len_trim(error)==0,error)
  source = 0.0D0
  call wet_cell_rate(geom,flux,source,rate,error)
  call require(len_trim(error)==0,error)
  call require(maxval(abs(rate))<2.0D-13,'divergence-free transport did not preserve constant PV')

  if (ranks > 1) then
    candidate = ranks
    if (geom%rank > 0 .and. any(geom%faces%owner==geom%rank)) candidate = geom%rank
    call mpi_allreduce(candidate,bad_rank,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (bad_rank < ranks) then
      if (geom%rank == bad_rank) then
        f = findloc(geom%faces%owner,geom%rank,dim=1)
        trace(1,1,1,f) = nan
      end if
      call transport_fluxes(geom,velocity,trace,flux,peak,error)
      call require(error=='nonfinite owned cut-cell face velocity or trace','rank-local NaN was accepted')
    end if
  end if

  ! Affine source-adjusted extrapolation at a non-midpoint time, with a
  ! quadratic accumulated source and independent limiter case.
  actual = staggered_trace(1.0D0,1.734375D0,0.375D0,0.234375D0,2.0D0,-10.0D0,10.0D0)
  call require(abs(actual-13.0D0/3.0D0)<1.0D-14,'affine trace extrapolation failed')
  actual = staggered_trace(1.0D0,9.0D0,0.375D0,0.234375D0,2.0D0,0.75D0,1.5D0)
  call require(actual==3.5D0,'limiter discarded the net physical source')

  ! Verify all positive-area cell centroids, including asymmetric cut
  ! cells. Observations at unequal times are affine face/cell means.
  samples = reshape([-0.5D0,0.25D0, 0.0D0,-0.5D0, -0.875D0,0.375D0],[2,3])
  ends = reshape([0.5D0,0.0D0, 0.5D0,0.5D0],[2,2])
  times = [-0.23D0,-0.37D0,-0.11D0]
  sources = times+times**2
  gradient = [1.5D0,-2.0D0]
  vel = [0.31D0,-0.17D0]
  do k = 1,3
    observations(k) = 0.4D0+dot_product(gradient,samples(:,k)-times(k)*vel)+sources(k)
  end do
  do k = 1,2
    expected(k) = 0.4D0+dot_product(gradient,ends(:,k)-0.17D0*vel)+0.1989D0
  end do
  do j = 1,geom%nj
    do i = 1,geom%ni
      if (geom%area(i,j) == 0.0D0) cycle
      centroid = geom%moments(2:3,i,j)
      exact = 0.4D0+dot_product(gradient,centroid)
      call affine_cabaret_trace(centroid,exact,samples,times,sources,observations,vel,ends,0.17D0, &
                                0.1989D0,.false.,values,error)
      call require(len_trim(error)==0,error)
      call require(maxval(abs(values-expected))<2.0D-14,'affine staggered moment recovery failed')
      call affine_cabaret_trace(centroid,exact,samples,times,sources,observations,vel,ends,0.17D0, &
                                0.1989D0,.true.,values,error)
      call require(len_trim(error)==0,error)
      corrected = observations-sources
      call require(maxval(values-0.1989D0)<=max(exact,maxval(corrected))+1.0D-14 .and. &
                   minval(values-0.1989D0)>=min(exact,minval(corrected))-1.0D-14,'source-aware bounds failed')
    end do
  end do
  samples = 0.0D0
  call affine_cabaret_trace([0.0D0,0.0D0],0.4D0,samples,times,sources,observations,vel,ends,0.17D0, &
                            0.1989D0,.false.,values,error)
  call require(error=='rank-deficient staggered characteristic moments','rank-deficient trace was accepted')

  ! Use actual open-face and wet-cell centroids at unequal stage times.
  ! A reversal supplies a different inflow segment on each side of the
  ! same face. MPI neighbors must exchange the cell-2 traces back to owners.
  allocate(old_trace, mold=trace)
  allocate(cell_velocity(2,geom%ni,geom%nj,layers),source_in(geom%ni,geom%nj,layers), &
           source_out(geom%ni,geom%nj,layers))
  call setup_cut_cabaret(geom,cut_plan,error)
  call require(len_trim(error)==0,error)
  stage_dt = 0.375D0
  source_rate = 0.75D0+dot_product(gradient,vel)
  source_in = -0.5D0*stage_dt*source_rate
  source_out = -source_in
  do j = 1,geom%nj
    do i = 1,geom%ni
      point = [real(geom%j0+j-1,real_8),real(n-geom%i0-i,real_8)]+geom%moments(2:3,i,j)
      if (geom%area(i,j) == 0.0D0) then
        psi(i,j,:) = nan
      else
        psi(i,j,:) = 0.4D0+dot_product(gradient,point)
      end if
      do l = 1,layers
        cell_velocity(:,i,j,l) = vel
      end do
    end do
  end do
  trace = nan
  do f = 1,size(geom%faces)
    ends = 0.5D0*geom%faces(f)%endpoints
    do l = 1,layers
      velocity(:,l,f) = [-1.0D0,3.0D0]*(-1.0D0)**l
    end do
    do cell = 1,2
      do k = 1,2
        exact = 0.4D0+dot_product(gradient,ends(:,k))
        old_trace(k,cell,:,f) = exact-0.375D0*stage_dt
        gi = geom%faces(f)%cells(1,cell)
        gj = geom%faces(f)%cells(2,cell)
        if (geom%owner(gi,gj)==geom%rank .and. geom%global%wedge(gi,gj)==255) &
          trace(k,cell,:,f) = exact+0.375D0*stage_dt
      end do
    end do
  end do
  call cut_cabaret_traces(geom,cut_plan,psi,old_trace,velocity,cell_velocity,source_in,source_out, &
                          stage_dt,.false.,trace,fallback,error)
  call require(len_trim(error)==0,error)
  call require(fallback==0,'affine reversing-face fixture lost a characteristic direction')
  call gather_second_face_trace(geom,trace,error)
  call require(len_trim(error)==0,error)
  face_error = 0.0D0
  do f = 1,size(geom%faces)
    if (geom%faces(f)%owner/=geom%rank) cycle
    ends = 0.5D0*geom%faces(f)%endpoints
    do k = 1,2
      exact = 0.4D0+dot_product(gradient,ends(:,k))+0.375D0*stage_dt
      call require(all(ieee_is_finite(trace(k,:,:,f))),'face owner did not receive both traces')
      face_error = max(face_error,maxval(abs(trace(k,:,:,f)-exact)))
    end do
  end do
  call mpi_allreduce(face_error,exact,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
  call require(exact<1.0D-12,'geometry-aware affine face trace failed')
  call transport_fluxes(geom,velocity,trace,flux,peak,error)
  call require(len_trim(error)==0,error)

  if (ranks>1) then
    call cut_cabaret_traces(geom,cut_plan,psi,old_trace,velocity,cell_velocity,source_in,source_out, &
                            merge(stage_dt,2*stage_dt,geom%rank/=1),.false.,trace,fallback,error)
    call require(error=='cut CABARET stage differs across ranks','rank-local stage time was accepted')
  end if
  psi = 1.0D0
  velocity = 0.0D0
  old_trace = 1.0D0
  trace = 1.0D0+source_out(1,1,1)
  call cut_cabaret_traces(geom,cut_plan,psi,old_trace,velocity,cell_velocity,source_in,source_out, &
                          stage_dt,.true.,trace,fallback,error)
  call require(len_trim(error)==0,error)
  local_min = minval(geom%area,mask=geom%area>0)
  call mpi_allreduce(local_min,min_area,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
  if (min_area < geom%h**2/4) call require(fallback>0,'one-open-face degeneracy was not counted')
  call require(max(0.0D0,maxval(abs(trace-1.0D0-source_out(1,1,1))))<1.0D-13,'fallback lost constant PV')

  ! Real half-stage spacings from three unequal integer step lengths.
  t = [0.0D0,0.125D0,0.5D0,0.875D0]
  target = 0.5D0*(t(3)+t(4))
  z = 0.5D0*(t(2)+t(3))
  uq = 0.5D0*(t(1)+t(2))
  integral = integrate_stage_rate(2*uq-1,2*z-1,z-uq,target-z)
  call require(abs(integral-(target*target-target-z*z+z))<1.0D-14,'unequal-step affine forcing failed')
  volume = 0.5D0*(t(3)-t(2))
  exact = 0.5D0*(t(4)-t(3))
  actual = source_substep(2*uq-1,2*z-1,z-uq,0.0D0,volume)+ &
           source_substep(2*uq-1,2*z-1,z-uq,volume,exact)
  call require(abs(actual-integral)<1.0D-14,'predictor and corrector did not integrate the same stage rate')
  coarse = temporal_error(16)
  fine = temporal_error(32)
  call require(coarse/fine > 3.75D0,'stage quadrature failed smooth time refinement')
  call require(drag_error(16)/drag_error(32)>3.5D0,'explicit drag stage integration lost temporal order')


  local_min = minval(geom%area,mask=geom%area>0)
  call mpi_allreduce(local_min,min_area,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
  if (rank == 0) print '(a,4es14.6,a)', 'cut-transport errors, CFL, area, time order: ', &
    difference,reference_peak,min_area/geom%h**2,log(coarse/fine)/log(2.0D0),' PASS'
  call mpi_comm_free(comm,ierr)
  call parallel_finalise()

contains

  function temporal_error(steps) result(err)
    integer, intent(in) :: steps
    real(real_8) :: err,current,prev,forward_dt,previous_dt,total
    integer :: step

    prev = -0.8D0/steps
    current = 0.0D0
    total = 0.0D0
    do step = 1,steps
      previous_dt = current-prev
      forward_dt = merge(0.8D0,1.2D0,mod(step,2)==0)/steps
      total = total+integrate_stage_rate(exp(prev),exp(current),previous_dt,forward_dt)
      prev = current
      current = current+forward_dt
    end do
    err = abs(total-(exp(current)-1))
  end function temporal_error

  function drag_error(steps) result(err)
    integer, intent(in) :: steps
    real(real_8) :: err,previous_dt,forward_dt,t,current,previous,next
    integer :: step

    ! Gamma_t=-r Gamma, with the same AB2 source used for the PV curl
    ! of a linear drag. A single startup Euler interval is O(dt**2).
    current = 1.0D0
    previous_dt = 0.8D0/steps
    t = previous_dt
    previous = current
    current = current-previous_dt*0.7D0*current
    do step = 2,steps
      forward_dt = merge(1.2D0,0.8D0,mod(step,2)==0)/steps
      next = current+integrate_stage_rate(-0.7D0*previous,-0.7D0*current,previous_dt,forward_dt)
      previous = current
      current = next
      t = t+forward_dt
      previous_dt = forward_dt
    end do
    err = abs(current-exp(-0.7D0*t))
  end function drag_error

  subroutine require(test,message)
    logical, intent(in) :: test
    character(*), intent(in) :: message
    integer :: status

    if (test) return
    print '(a)', trim(message)
    call mpi_abort(MPI_COMM_WORLD,1,status)
  end subroutine require

end program test_irregular_transport
