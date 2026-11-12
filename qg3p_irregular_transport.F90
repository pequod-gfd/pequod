! Spatial flux and time-quadrature kernels for experimental cut-cell CABARET.
module qg_irregular_transport
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use qg_irregular_parallel, only: distributed_geometry, share_face_flux, broadcast_error, donor_plan, &
                                   build_donor_plan, gather_donors, half_coordinates, make_sector_moments, &
                                   cell_mean_moments, gather_second_face_trace
  use qg_moments, only: n_moments
  use qg_irregular_operators, only: sum_open_flux
  use qg_parameters, only: real_8
  implicit none
  include 'mpif.h'
  private
  public :: split_face_flux, transport_fluxes, wet_cell_rate
  public :: staggered_trace, affine_cabaret_trace, integrate_stage_rate, source_substep
  public :: source_momentum_budget
  public :: cut_cabaret_plan, setup_cut_cabaret, cut_cabaret_traces
  public :: cabaret_face_traces, cabaret_cell_velocity

  type cut_cell_stencil
    integer :: ij(2),face(4),neighbor(4)
    real(real_8) :: centroid(2),neighbor_centroid(2,4)
  end type cut_cell_stencil

  type cut_cabaret_plan
    private
    type(cut_cell_stencil), allocatable :: cells(:)
    type(donor_plan) :: donors
    integer :: donor_count
    logical :: includes_full = .false.
  end type cut_cabaret_plan

contains

  ! Endpoint velocities use the same normal as the face incidence. Both traces
  ! are affine in the normalized face coordinate and indexed by adjacent cell.
  ! Positive outflow leaves cell 1. Negative velocity leaves cell 2.
  pure subroutine split_face_flux(length, velocity, trace, flux, outflow)
    real(real_8), intent(in) :: length,velocity(2),trace(2,2)
    real(real_8), intent(out) :: flux,outflow(2)
    real(real_8) :: limits(0:2),a,b,width,u_a,u_b,q_a,q_b
    integer :: intervals,interval,side

    limits = [0.0D0,1.0D0,1.0D0]
    intervals = 1
    if (minval(velocity) < 0.0D0 .and. maxval(velocity) > 0.0D0) then
      limits(1) = -velocity(1)/(velocity(2)-velocity(1))
      intervals = 2
    end if

    flux = 0.0D0
    outflow = 0.0D0
    do interval = 1,intervals
      a = limits(interval-1)
      b = limits(interval)
      width = length*(b-a)
      u_a = (1.0D0-a)*velocity(1)+a*velocity(2)
      u_b = (1.0D0-b)*velocity(1)+b*velocity(2)
      side = merge(1,2,u_a+u_b >= 0.0D0)
      q_a = (1.0D0-a)*trace(1,side)+a*trace(2,side)
      q_b = (1.0D0-b)*trace(1,side)+b*trace(2,side)
      flux = flux+width*(2*u_a*q_a+u_a*q_b+u_b*q_a+2*u_b*q_b)/6.0D0
      outflow(side) = outflow(side)+0.5D0*width*abs(u_a+u_b)
    end do
  end subroutine split_face_flux

  ! Only the face owner supplies endpoint velocities and both reconstructed
  ! traces. Share the flux and both directed volumes together. A reversing
  ! face can therefore contribute to the CFL restriction on both cells.
  subroutine transport_fluxes(geom, velocity, trace, flux, max_outgoing_rate, error)
    type(distributed_geometry), intent(in) :: geom
    real(real_8), intent(in) :: velocity(:,:,:),trace(:,:,:,:)
    real(real_8), intent(out) :: flux(:,:),max_outgoing_rate
    character(*), intent(out) :: error
    real(real_8), allocatable :: packet(:,:)
    real(real_8) :: length,volume(2),local_rate
    integer :: f,layer,layers,n_faces,i,j,side,cell,ierr

    error = ''
    layers = size(flux,1)
    n_faces = size(geom%faces)
    if (any(shape(velocity) /= [2,layers,n_faces]) .or. any(shape(trace) /= [2,2,layers,n_faces]) .or. &
        size(flux,2) /= n_faces .or. layers < 1) error = 'incompatible cut-cell transport arrays'
    call broadcast_error(geom%comm,geom%rank,geom%ranks,error)
    if (len_trim(error) /= 0) return

    allocate(packet(3*layers,n_faces))
    packet = 0.0D0
    do f = 1,n_faces
      if (geom%faces(f)%owner /= geom%rank) cycle
      if (.not. all(ieee_is_finite(velocity(:,:,f))) .or. .not. all(ieee_is_finite(trace(:,:,:,f)))) then
        error = 'nonfinite owned cut-cell face velocity or trace'
        exit
      end if
      length = 0.5D0*sum(abs(geom%faces(f)%endpoints(:,2)-geom%faces(f)%endpoints(:,1)))
      do layer = 1,layers
        call split_face_flux(length,velocity(:,layer,f),trace(:,:,layer,f),packet(layer,f),volume)
        packet(layers+layer,f) = volume(1)
        packet(2*layers+layer,f) = volume(2)
      end do
    end do
    call broadcast_error(geom%comm,geom%rank,geom%ranks,error)
    if (len_trim(error) /= 0) return
    call share_face_flux(geom,packet,error)
    if (len_trim(error) /= 0) return
    flux = packet(:layers,:)

    local_rate = 0.0D0
    do j = 1,geom%nj
      do i = 1,geom%ni
        if (geom%area(i,j) == 0.0D0) cycle
        do layer = 1,layers
          length = 0.0D0
          do side = 1,4
            f = geom%face_index(side,i,j)
            if (f == 0) cycle
            cell = merge(1,2,side <= 2)
            length = length+packet(cell*layers+layer,f)
          end do
          local_rate = max(local_rate,length*geom%h**2/geom%area(i,j))
        end do
      end do
    end do
    call mpi_allreduce(local_rate,max_outgoing_rate,1,MPI_DOUBLE_PRECISION,MPI_MAX,geom%comm,ierr)
  end subroutine transport_fluxes

  ! The source input is an integral over each wet cell in native grid units.
  ! The same operator can supply the predictor and corrector residuals.
  subroutine wet_cell_rate(geom, flux, volume_source, rate, error)
    type(distributed_geometry), intent(in) :: geom
    real(real_8), intent(in) :: flux(:,:),volume_source(:,:,:)
    real(real_8), intent(out) :: rate(:,:,:)
    character(*), intent(out) :: error
    integer :: i,j,layers

    error = ''
    layers = size(rate,3)
    if (any(shape(rate) /= [geom%ni,geom%nj,layers]) .or. any(shape(volume_source) /= shape(rate)) .or. &
        any(shape(flux) /= [layers,size(geom%faces)]) .or. layers < 1) error = 'incompatible cut-cell residual arrays'
    call broadcast_error(geom%comm,geom%rank,geom%ranks,error)
    if (len_trim(error) /= 0) return
    call sum_open_flux(geom,flux,rate)
    do j = 1,geom%nj
      do i = 1,geom%ni
        if (geom%area(i,j) == 0.0D0) cycle
        rate(i,j,:) = (volume_source(i,j,:)-rate(i,j,:))*geom%h**2/geom%area(i,j)
      end do
    end do
  end subroutine wet_cell_rate

  ! Extrapolate in source-adjusted variables. fraction is the time fraction
  ! at the intermediate observation, strictly between zero and one. The
  ! caller aligns sample with the characteristic using the spatial gradient
  ! and supplies bounds at the entry time in source-adjusted variables.
  pure elemental function staggered_trace(entry, sample, fraction, source_sample, source_exit, &
                                           lower, upper) result(exit_trace)
    real(real_8), intent(in) :: entry,sample,fraction,source_sample,source_exit,lower,upper
    real(real_8) :: exit_trace,adjusted

    adjusted = entry+(sample-source_sample-entry)/fraction
    exit_trace = min(max(adjusted,lower),upper)+source_exit
  end function staggered_trace

  ! Affine space-time polynomial constrained by q_t+u.grad(q)=S(t).
  ! Time offsets are relative to the half-stage cell observation. Correct
  ! each datum by the source integral relative to that same time. The input
  ! positions are actual cell or face centroids in native grid coordinates.
  subroutine affine_cabaret_trace(centroid,mean,positions,times,sources,observations,velocity, &
                                   endpoints,exit_time,exit_source,limit,trace,error,reference_gradient)
    real(real_8), intent(in) :: centroid(2),mean,positions(:,:),times(:),sources(:),observations(:)
    real(real_8), intent(in) :: velocity(2),endpoints(2,2),exit_time,exit_source
    logical, intent(in) :: limit
    real(real_8), intent(out) :: trace(2)
    character(*), intent(out) :: error
    real(real_8), optional, intent(in) :: reference_gradient(2)
    real(real_8) :: moment(2,2),rhs(2),offset(2),corrected,gradient(2),det,scale,lower,upper
    integer :: k,n

    error = ''
    n = size(times)
    if (any(shape(positions) /= [2,n]) .or. size(sources) /= n .or. size(observations) /= n) then
      error = 'invalid staggered affine reconstruction observations'
      return
    end if
    if (.not. all(ieee_is_finite([centroid,mean,velocity,exit_time,exit_source])) .or. &
        .not. all(ieee_is_finite(endpoints)) .or. .not. all(ieee_is_finite(positions)) .or. &
        .not. all(ieee_is_finite(times)) .or. .not. all(ieee_is_finite(sources)) .or. &
        .not. all(ieee_is_finite(observations))) then
      error = 'nonfinite staggered affine reconstruction observations'
      return
    end if

    moment = 0.0D0
    rhs = 0.0D0
    lower = mean
    upper = mean
    do k = 1,n
      offset = positions(:,k)-centroid-times(k)*velocity
      corrected = observations(k)-sources(k)
      moment(:,1) = moment(:,1)+offset*offset(1)
      moment(:,2) = moment(:,2)+offset*offset(2)
      rhs = rhs+offset*(corrected-mean)
      lower = min(lower,corrected)
      upper = max(upper,corrected)
    end do
    scale = moment(1,1)+moment(2,2)
    det = moment(1,1)*moment(2,2)-moment(1,2)*moment(2,1)
    if (det <= 1.0D-10*scale**2 .or. scale < 1.0D-12) then
      error = 'rank-deficient staggered characteristic moments'
      return
    end if
    gradient(1) = (moment(2,2)*rhs(1)-moment(1,2)*rhs(2))/det
    gradient(2) = (moment(1,1)*rhs(2)-moment(2,1)*rhs(1))/det
    do k = 1,2
      offset = endpoints(:,k)-centroid-exit_time*velocity
      corrected = mean+dot_product(offset,gradient)
      trace(k) = corrected
    end do
    if (limit) then
      if (present(reference_gradient)) then
        call limit_characteristic_trace(centroid,mean,positions,times,sources,observations,velocity, &
                                        endpoints,exit_time,reference_gradient,trace)
      else
        trace = min(max(trace,lower),upper)
      end if
    end if
    trace = trace+exit_source
  end subroutine affine_cabaret_trace

  ! Geometry rows for the staggered trace. Each cell owns its half-level
  ! observation. Neighbor means come only from shared open faces.
  subroutine setup_cut_cabaret(geom,plan,error,include_full)
    type(distributed_geometry), intent(in) :: geom
    type(cut_cabaret_plan), intent(out) :: plan
    character(*), intent(out) :: error
    logical, optional, intent(in) :: include_full
    integer, allocatable :: donor(:,:)
    integer :: cells,k,i,j,side,f,cell,n,ij(2),mask,ni
    real(real_8) :: sector(n_moments,0:7),moments(n_moments)

    ni = size(geom%owner,1)
    cells = count(geom%area > 0.0D0 .and. geom%area < geom%h**2)
    plan%includes_full = .false.
    if (present(include_full)) plan%includes_full = include_full
    if (plan%includes_full) cells = count(geom%area > 0.0D0)
    allocate(plan%cells(cells),donor(2,4*cells))
    sector = make_sector_moments()
    k = 0
    n = 0
    do j = 1,geom%nj
      do i = 1,geom%ni
        if (geom%area(i,j) <= 0.0D0) cycle
        if (.not. plan%includes_full .and. geom%area(i,j) >= geom%h**2) cycle
        k = k+1
        associate(stencil => plan%cells(k))
          stencil%ij = [i,j]
          stencil%face = geom%face_index(:,i,j)
          stencil%centroid = 0.5D0*half_coordinates(geom%i0+i,geom%j0+j,ni)+geom%moments(2:3,i,j)
          stencil%neighbor = 0
          do side = 1,4
            f = stencil%face(side)
            if (f == 0) cycle
            cell = merge(2,1,side<=2)
            ij = geom%faces(f)%cells(:,cell)
            mask = geom%global%wedge(ij(1),ij(2))
            moments = cell_mean_moments(mask,sector)
            n = n+1
            donor(:,n) = ij
            stencil%neighbor(side) = n
            stencil%neighbor_centroid(:,side) = 0.5D0*half_coordinates(ij(1),ij(2),ni)+moments(2:3)
          end do
        end associate
      end do
    end do
    plan%donor_count = n
    call build_donor_plan(geom,donor(:,:n),plan%donors,error,cell_means=.true.)
  end subroutine setup_cut_cabaret

  ! Select incoming history on the part of each previous face that entered
  ! the cut cell. A sign reversal changes the observation's centroid and
  ! mean. Full cells use opposing-face reflection, when included in the
  ! plan. This kernel fills locally owned slots for the shared exchange.
  subroutine cut_cabaret_traces(geom,plan,mean_half,old_trace,old_speed,cell_velocity,source_in,source_out, &
                                dt,limit,trace,fallback_count,error)
    type(distributed_geometry), intent(in) :: geom
    type(cut_cabaret_plan), intent(in) :: plan
    real(real_8), intent(in) :: mean_half(:,:,:),old_trace(:,:,:,:),old_speed(:,:,:),cell_velocity(:,:,:,:)
    real(real_8), intent(in) :: source_in(:,:,:),source_out(:,:,:),dt
    logical, intent(in) :: limit
    real(real_8), intent(inout) :: trace(:,:,:,:)
    integer, intent(out) :: fallback_count
    character(*), intent(out) :: error
    real(real_8), allocatable :: neighbors(:,:)
    real(real_8) :: positions(2,8),times(8),sources(8),observations(8),endpoint(2,2),speed(2),centre
    real(real_8) :: stage(3),root_stage(3),gradient(2),adjacent(2),mean
    integer :: i,j,k,layer,nl,n,f,side,cell,nd,spatial,local_count,ierr,other
    logical :: full
    character(80) :: fit_error

    error = ''
    nl = size(mean_half,3)
    n = size(geom%faces)
    if (nl < 1 .or. any(shape(mean_half) /= [geom%ni,geom%nj,nl]) .or. &
        any(shape(old_speed) /= [2,nl,n]) .or. any(shape(old_trace) /= [2,2,nl,n]) .or. &
        any(shape(trace) /= [2,2,nl,n]) .or. any(shape(cell_velocity) /= [2,geom%ni,geom%nj,nl]) .or. &
        any(shape(source_in) /= shape(mean_half)) .or. any(shape(source_out) /= shape(mean_half))) &
      error = 'incompatible cut CABARET arrays'
    if (.not. allocated(plan%cells)) error = 'cut CABARET plan is uninitialized'
    if (.not. ieee_is_finite(dt) .or. dt <= 0.0D0) error = 'invalid cut CABARET step'
    call broadcast_error(geom%comm,geom%rank,geom%ranks,error)
    if (len_trim(error) /= 0) return
    stage = [dt,merge(1.0D0,0.0D0,limit),merge(1.0D0,0.0D0,plan%includes_full)]
    root_stage = stage
    call mpi_bcast(root_stage,size(root_stage),MPI_DOUBLE_PRECISION,0,geom%comm,ierr)
    if (any(root_stage /= stage)) error = 'cut CABARET stage differs across ranks'
    call broadcast_error(geom%comm,geom%rank,geom%ranks,error)
    if (len_trim(error) /= 0) return
    allocate(neighbors(nl,plan%donor_count))
    call gather_donors(plan%donors,mean_half,neighbors)

    local_count = 0
    do k = 1,size(plan%cells)
      associate(stencil => plan%cells(k))
        i = stencil%ij(1)
        j = stencil%ij(2)
        full = geom%area(i,j)==geom%h**2
        do layer = 1,nl
          mean = mean_half(i,j,layer)
          nd = 0
          do side = 1,4
            f = stencil%face(side)
            if (f == 0) cycle
            n = stencil%neighbor(side)
            nd = nd+1
            positions(:,nd) = stencil%neighbor_centroid(:,side)
            times(nd) = 0.0D0
            sources(nd) = 0.0D0
            observations(nd) = neighbors(layer,n)
          end do
          spatial = nd
          if (limit .and. plan%includes_full) call spatial_gradient(stencil%centroid,mean, &
            positions(:,:spatial),observations(:spatial),gradient)
          do side = 1,4
            if (full) exit
            f = stencil%face(side)
            if (f == 0) cycle
            ! The shared normal points outward from side 1 and inward to 2.
            cell = merge(1,2,side<=2)
            speed = old_speed(:,layer,f)*merge(1.0D0,-1.0D0,cell==1)
            if (.not. all(ieee_is_finite(speed))) then
              error = 'nonfinite cut CABARET incoming face speed'
              cycle
            end if
            if (minval(speed) >= 0.0D0) cycle
            centre = 0.5D0
            if (maxval(speed) > 0.0D0) then
              centre = 0.5D0*(merge(0.0D0,1.0D0,speed(1)<0.0D0)-speed(1)/(speed(2)-speed(1)))
            end if
            endpoint = 0.5D0*geom%faces(f)%endpoints
            nd = nd+1
            positions(:,nd) = (1-centre)*endpoint(:,1)+centre*endpoint(:,2)
            times(nd) = -0.5D0*dt
            sources(nd) = source_in(i,j,layer)
            observations(nd) = dot_product([1-centre,centre],old_trace(:,3-cell,layer,f))
          end do
          do side = 1,4
            f = stencil%face(side)
            if (f == 0) cycle
            cell = merge(1,2,side<=2)
            endpoint = 0.5D0*geom%faces(f)%endpoints
            if (full) then
              other = stencil%face(modulo(side+1,4)+1)
              adjacent = upwind_history(old_trace(:,:,layer,other),old_speed(:,layer,other))
              trace(:,cell,layer,f) = 2*mean-(adjacent(2:1:-1)-source_in(i,j,layer))
              if (limit) then
                positions(:,spatial+1:spatial+2) = 0.5D0*geom%faces(other)%endpoints
                positions(:,spatial+3:spatial+4) = endpoint
                times(spatial+1:spatial+4) = -0.5D0*dt
                sources(spatial+1:spatial+4) = source_in(i,j,layer)
                observations(spatial+1:spatial+2) = adjacent
                observations(spatial+3:spatial+4) = upwind_history(old_trace(:,:,layer,f),old_speed(:,layer,f))
                call limit_characteristic_trace(stencil%centroid,mean,positions(:,:spatial+4),times(:spatial+4), &
                  sources(:spatial+4),observations(:spatial+4),cell_velocity(:,i,j,layer),endpoint,0.5D0*dt, &
                  gradient,trace(:,cell,layer,f))
              end if
              trace(:,cell,layer,f) = trace(:,cell,layer,f)+source_out(i,j,layer)
              cycle
            end if
            if (limit .and. plan%includes_full) then
              call affine_cabaret_trace(stencil%centroid,mean,positions(:,:nd),times(:nd),sources(:nd), &
                  observations(:nd),cell_velocity(:,i,j,layer),endpoint,0.5D0*dt,source_out(i,j,layer), &
                  limit,trace(:,cell,layer,f),fit_error,reference_gradient=gradient)
            else
              call affine_cabaret_trace(stencil%centroid,mean,positions(:,:nd),times(:nd),sources(:nd), &
                  observations(:nd),cell_velocity(:,i,j,layer),endpoint,0.5D0*dt,source_out(i,j,layer), &
                  limit,trace(:,cell,layer,f),fit_error)
            end if
            select case(trim(fit_error))
              case('')
              case('rank-deficient staggered characteristic moments')
                trace(:,cell,layer,f) = mean_half(i,j,layer)+source_out(i,j,layer)
                local_count = local_count+1
              case default
                error = fit_error
            end select
          end do
        end do
      end associate
    end do
    call broadcast_error(geom%comm,geom%rank,geom%ranks,error)
    if (len_trim(error) /= 0) return
    call mpi_allreduce(local_count,fallback_count,1,MPI_INTEGER,MPI_SUM,geom%comm,ierr)
  end subroutine cut_cabaret_traces

  pure function upwind_history(trace,speed) result(value)
    real(real_8), intent(in) :: trace(2,2),speed(2)
    real(real_8) :: value(2)

    value = merge(trace(:,1),trace(:,2),speed>=0.0D0)
  end function upwind_history

  ! Fit a reference slope using only the simultaneous neighboring means.
  ! Restrict its predictions at those centroids to their measured range.
  pure subroutine spatial_gradient(centre,mean,positions,observations,gradient)
    real(real_8), intent(in) :: centre(2),mean,positions(:,:),observations(:)
    real(real_8), intent(out) :: gradient(2)
    real(real_8) :: moment(2,2),rhs(2),offset(2),det,scale,lo,hi,change
    integer :: k

    moment = 0.0D0
    rhs = 0.0D0
    gradient = 0.0D0
    do k = 1,size(observations)
      offset = positions(:,k)-centre
      moment(:,1) = moment(:,1)+offset*offset(1)
      moment(:,2) = moment(:,2)+offset*offset(2)
      rhs = rhs+offset*(observations(k)-mean)
    end do
    scale = moment(1,1)+moment(2,2)
    det = moment(1,1)*moment(2,2)-moment(1,2)*moment(2,1)
    if (det <= 1.0D-10*scale**2 .or. scale < 1.0D-12) return
    gradient = [moment(2,2)*rhs(1)-moment(1,2)*rhs(2),moment(1,1)*rhs(2)-moment(1,2)*rhs(1)]/det
    lo = min(0.0D0,minval(observations)-mean)
    hi = max(0.0D0,maxval(observations)-mean)
    scale = 1.0D0
    do k = 1,size(observations)
      change = dot_product(gradient,positions(:,k)-centre)
      if (change < lo) scale = min(scale,lo/change)
      if (change > hi) scale = min(scale,hi/change)
    end do
    gradient = scale*gradient
  end subroutine spatial_gradient

  ! Apply the bounds to departures from the reference affine characteristic.
  ! For smooth fields these departures are O(h**2) when dt scales with h.
  ! This preserves an affine solution even at a face outside the hull of
  ! the available wet centroids. A zero reference gives scalar bounds.
  pure subroutine limit_characteristic_trace(centre,mean,positions,times,sources,observations,velocity, &
                                              endpoints,exit_time,gradient,trace)
    real(real_8), intent(in) :: centre(2),mean,positions(:,:),times(:),sources(:),observations(:),velocity(2)
    real(real_8), intent(in) :: endpoints(2,2),exit_time,gradient(2)
    real(real_8), intent(inout) :: trace(2)
    real(real_8) :: lo,hi,reference,residual
    integer :: k

    lo = 0.0D0
    hi = 0.0D0
    do k = 1,size(observations)
      reference = mean+dot_product(gradient,positions(:,k)-centre-times(k)*velocity)
      residual = observations(k)-sources(k)-reference
      lo = min(lo,residual)
      hi = max(hi,residual)
    end do
    do k = 1,2
      reference = mean+dot_product(gradient,endpoints(:,k)-centre-exit_time*velocity)
      trace(k) = reference+min(max(trace(k)-reference,lo),hi)
    end do
  end subroutine limit_characteristic_trace

  ! Mean velocity from the first moments of the normal flux. Impermeable
  ! walls contribute zero. For a solenoidal velocity, div(x*u) = u_x
  ! and div(y*u) = u_y, including on the clipped control volumes.
  subroutine cabaret_cell_velocity(geom,speed,cell_velocity)
    type(distributed_geometry), intent(in) :: geom
    real(real_8), intent(in) :: speed(:,:,:)
    real(real_8), intent(out) :: cell_velocity(:,:,:,:)
    real(real_8) :: centre(2),a(2),b(2),length,sign
    integer :: i,j,side,f,k,ni

    ni = size(geom%owner,1)
    cell_velocity = 0.0D0
    do j = 1,geom%nj
      do i = 1,geom%ni
        if (geom%area(i,j) == 0.0D0) cycle
        centre = 0.5D0*half_coordinates(geom%i0+i,geom%j0+j,ni)+geom%moments(2:3,i,j)
        do side = 1,4
          f = geom%face_index(side,i,j)
          if (f == 0) cycle
          sign = merge(1.0D0,-1.0D0,side<=2)
          a = 0.5D0*geom%faces(f)%endpoints(:,1)-centre
          b = 0.5D0*geom%faces(f)%endpoints(:,2)-centre
          length = sum(abs(b-a))
          do k = 1,size(speed,2)
            cell_velocity(:,i,j,k) = cell_velocity(:,i,j,k)+sign*length/6 * &
              ((2*a+b)*speed(1,k,f)+(a+2*b)*speed(2,k,f))
          end do
        end do
        cell_velocity(:,i,j,:) = cell_velocity(:,i,j,:)*(geom%h**2/geom%area(i,j))
      end do
    end do
  end subroutine cabaret_cell_velocity

  subroutine cabaret_face_traces(geom,plan,mean_half,old_trace,old_speed,cell_velocity,source_in,source_out, &
                                 dt,limit,trace,fallback_count,error)
    type(distributed_geometry), intent(in) :: geom
    type(cut_cabaret_plan), intent(in) :: plan
    real(real_8), intent(in) :: mean_half(:,:,:),old_trace(:,:,:,:),old_speed(:,:,:),cell_velocity(:,:,:,:)
    real(real_8), intent(in) :: source_in(:,:,:),source_out(:,:,:),dt
    logical, intent(in) :: limit
    real(real_8), intent(out) :: trace(:,:,:,:)
    integer, intent(out) :: fallback_count
    character(*), intent(out) :: error
    real(real_8), allocatable :: packet(:,:)
    integer :: nl,nf

    nl = size(mean_half,3)
    nf = size(geom%faces)
    error = ''
    if (.not. plan%includes_full) error = 'CABARET face history requires a plan including full cells'
    call broadcast_error(geom%comm,geom%rank,geom%ranks,error)
    if (len_trim(error) /= 0) return
    trace = 0.0D0
    call cut_cabaret_traces(geom,plan,mean_half,old_trace,old_speed,cell_velocity,source_in,source_out, &
                            dt,limit,trace,fallback_count,error)
    if (len_trim(error) /= 0) return
    call gather_second_face_trace(geom,trace,error)
    if (len_trim(error) /= 0) return
    ! Both cells need both traces as history on the next step.
    packet = reshape(trace,[4*nl,nf])
    call share_face_flux(geom,packet,error)
    if (len_trim(error) == 0) trace = reshape(packet,shape(trace))
  end subroutine cabaret_face_traces

  ! Check global curl/torque compatibility for a stage with a supplied
  ! momentum force. The parent torque includes the outer boundary.
  ! This cannot recover individual torques from a supplied curl.
  subroutine source_momentum_budget(geom, source, torque, residual, error)
    type(distributed_geometry), intent(in) :: geom
    real(real_8), intent(in) :: source(:,:,:),torque(0:,:)
    real(real_8), intent(out) :: residual(:)
    character(*), intent(out) :: error
    real(real_8) :: local(2,size(residual)),total(2,size(residual)),root_torque(size(torque,1),size(torque,2))
    integer :: i,j,nl,counts(2),ierr

    error = ''
    nl = size(residual)
    if (any(shape(source) /= [geom%ni,geom%nj,nl]) .or. &
        any(shape(torque) /= [size(geom%global%loops),nl]) .or. nl < 1) error = 'invalid curl/torque budget shapes'
    if (.not. all(ieee_is_finite(torque))) error = 'nonfinite parent torque'
    call broadcast_error(geom%comm,geom%rank,geom%ranks,error)
    if (len_trim(error) /= 0) return
    call mpi_allreduce([nl,-nl],counts,2,MPI_INTEGER,MPI_MAX,geom%comm,ierr)
    if (counts(1) /= -counts(2)) then
      error = 'curl/torque layer counts differ across ranks'
      return
    end if
    root_torque = torque
    call mpi_bcast(root_torque,size(root_torque),MPI_DOUBLE_PRECISION,0,geom%comm,ierr)
    if (any(torque /= root_torque)) error = 'parent torque differs across ranks'
    local = 0.0D0
    do j = 1,geom%nj
      do i = 1,geom%ni
        if (geom%area(i,j) == 0.0D0) cycle
        if (.not. all(ieee_is_finite(source(i,j,:)))) error = 'nonfinite wet-cell curl source'
        local(1,:) = local(1,:)+source(i,j,:)
        local(2,:) = local(2,:)+abs(source(i,j,:))
      end do
    end do
    call broadcast_error(geom%comm,geom%rank,geom%ranks,error)
    if (len_trim(error) /= 0) return
    call mpi_allreduce(local,total,size(local),MPI_DOUBLE_PRECISION,MPI_SUM,geom%comm,ierr)
    residual = total(1,:)-sum(torque,dim=1)
    if (any(abs(residual) > 1.0D-11*max(tiny(0.0D0),total(2,:)+sum(abs(torque),dim=1)))) &
      error = 'volume source and parent boundary torque are incompatible'
  end subroutine source_momentum_budget

  ! For an interval of length dt starting offset after the last inversion,
  ! integrate the same linear rate history in PV and momentum. Predictor
  ! and corrector use distinct offsets but identical history.
  pure elemental function source_substep(previous,current,previous_dt,offset,dt) result(increment)
    real(real_8), intent(in) :: previous,current,previous_dt,offset,dt
    real(real_8) :: increment

    increment = dt*(current+((offset+0.5D0*dt)/previous_dt)*(current-previous))
  end function source_substep

  ! Integrate a linearly extrapolated stage rate over the forward interval.
  ! previous_dt is the distance between previous and current stage times.
  ! Use identical quadrature for parent momentum and integrated PV sources.
  pure elemental function integrate_stage_rate(previous, current, previous_dt, forward_dt) result(increment)
    real(real_8), intent(in) :: previous,current,previous_dt,forward_dt
    real(real_8) :: increment

    increment = source_substep(previous,current,previous_dt,0.0D0,forward_dt)
  end function integrate_stage_rate

end module qg_irregular_transport
