! Standalone, wind-driven axial cape experiment using the qg3 modal inversion.
program qg3_cape
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: output_unit
  use mpi, only: mpi_comm_rank, mpi_comm_size, mpi_bcast, mpi_dims_create, mpi_allreduce, &
                 MPI_COMM_WORLD, MPI_INTEGER, MPI_DOUBLE_PRECISION, MPI_CHARACTER, MPI_SUM, MPI_MAX
  use qg_parameters, only: real_8, day_s, basinscale, beta, u_scale, rd_bcl1, rd_bcl2, depths => h, visc, alpha, io_unit => unit
  use qg_geometry, only: write_rectangular_cape, solid_node, wall_node, fluid_node
  use qg_parallel, only: parallel_init, parallel_finalise, pabort, decomposition, parallel_data, request_cache, &
                         allocate, deallocate, gather, int2str
  use qg_diagnostics, only: stratification, allocate, deallocate, z_to_modes
  use qg_irregular_parallel, only: distributed_geometry, load_distributed_geometry, broadcast_error, &
                                   donor_plan, build_donor_plan, gather_donors, share_face_flux
  use qg_constraints, only: basin_mass_cache, setup_basin_mass, invert_basin_mass, release_basin_mass, fluid_integrals
  use qg_irregular_bc, only: vertex_plan, setup_vertex_plan, instantaneous_transport
  use qg_irregular_forcing, only: wind_curl_torque, cape_wind_curl_torque
  use qg_irregular_operators, only: sum_open_flux
  use qg_irregular_transport, only: transport_fluxes, wet_cell_rate, source_momentum_budget, cut_cabaret_plan, &
                                    setup_cut_cabaret, cabaret_face_traces, cabaret_cell_velocity
  use qg_cape_slip, only: cape_slip_cache, setup_cape_slip, release_cape_slip, cape_viscous_state, cape_implicit_step
  use qg_io, only: write_s_header, open_parallel_data, write_unformatted
  implicit none

  integer, parameter :: layers = 3
  real(real_8), parameter :: cfl_target = 0.45D0
  character(*), parameter :: geometry_path = 'cape_geometry.dat', state_path = 'cape_state.csv'
  integer :: nx=65,ny=65,cape_start=24,cape_width=16,cape_height=16
  real(real_8) :: basin_km=basinscale/1.0D5,run_days=5.0D0,max_dt_seconds=3600.0D0,beta_gradient=beta
  real(real_8) :: wind_stress=0.8D0,wind_shift=0.0D0,wind_tilt=0.2D0,wind_asymmetry=0.9D0
  character(32) :: wind_profile='cape_single'
  character(32) :: advection='upwind'
  real(real_8) :: viscosity_m2_s=visc/1.0D4,slip_length_km=alpha/1.0D5,snapshot_days=1.0D0
  namelist /cape/ nx,ny,cape_start,cape_width,cape_height,basin_km,run_days,max_dt_seconds,beta_gradient, &
                   wind_profile,wind_stress,wind_shift,wind_tilt,wind_asymmetry,viscosity_m2_s,slip_length_km,snapshot_days, &
                   advection

  type(decomposition), target :: decomp
  type(distributed_geometry), target :: geom
  type(stratification), target :: strat
  type(basin_mass_cache), target :: inverse
  type(cape_slip_cache) :: slip
  type(vertex_plan) :: vertices
  type(donor_plan) :: face_cells
  type(cut_cabaret_plan) :: cabaret
  integer :: ierr,rank,ranks,grid(2),settings(5),i,j,k,f,n,steps,retries,unit,status,snapshot_index
  integer :: fallback,peak_fallback,rejected_steps
  integer, allocatable :: iis(:),jjs(:),owned_faces(:),requested(:,:)
  real(real_8) :: real_settings(11),h_grid,time_scale,end_time,max_dt,time,dt,rate0,rate1,peak_cfl
  real(real_8) :: next_output,remaining,slip_defect,peak_slip_defect,wall_pv(layers),local_wall_pv(layers)
  real(real_8) :: target_mass(layers-1),mass(layers),wall_values(0:0,layers),torque(0:0,layers),budget(layers)
  real(real_8) :: wet_area,local_area,strength,wind_west,wind_south,beta_grid,length
  real(real_8) :: wind_total(layers),pv_total(layers),pv_scale(layers)
  real(real_8) :: solver_pv(layers),local_solver_pv(layers),budget_scale(layers)
  real(real_8) :: max_psi,local_psi,divergence_error,local_divergence
  real(real_8) :: previous_dt
  real(real_8), allocatable :: q(:,:,:),predicted(:,:,:),rhs0(:,:,:),rhs1(:,:,:),psi(:,:,:),wind(:,:,:), &
                               planetary(:,:,:),pv_at_faces(:,:),volume(:,:),velocity(:,:,:),trace(:,:,:,:),flux(:,:)
  real(real_8), allocatable :: diff_old(:,:,:),diff_new(:,:,:),step_rhs(:,:,:),next_state(:,:,:)
  real(real_8), allocatable :: wind_rate(:,:,:),diff_half(:,:,:),source_in(:,:,:),source_out(:,:,:),cell_velocity(:,:,:,:), &
                               speed_half(:,:,:),previous_half(:,:,:),speed_next(:,:,:),trace_next(:,:,:,:)
  character(1024) :: config_path,error
  logical :: already_exists,viscous,use_cabaret

  call parallel_init()
  call mpi_comm_rank(MPI_COMM_WORLD,rank,ierr)
  call mpi_comm_size(MPI_COMM_WORLD,ranks,ierr)
  error = ''
  if (rank==0) then
    if (command_argument_count()/=1) then
      error = 'usage: qg3_cape PATH_TO_CAPE_NAMELIST'
    else
      call get_command_argument(1,config_path,status=status)
      if (status==0) then
        open(newunit=unit,file=config_path,status='old',action='read',iostat=status,iomsg=error)
        if (status==0) then
          read(unit,nml=cape,iostat=status,iomsg=error)
          close(unit)
        end if
      else
        error = 'could not read the cape namelist path'
      end if
    end if
    if (len_trim(error)==0) then
      inquire(file=geometry_path,exist=already_exists)
      if (already_exists) error = 'cape_geometry.dat already exists. Use a fresh run directory.'
      inquire(file=state_path,exist=already_exists)
      if (already_exists) error = 'cape_state.csv already exists. Use a fresh run directory.'
    end if
  end if
  call broadcast_error(MPI_COMM_WORLD,rank,ranks,error)
  call check(error)
  settings = [nx,ny,cape_start,cape_width,cape_height]
  real_settings = [basin_km,run_days,max_dt_seconds,beta_gradient,wind_stress,wind_shift,wind_tilt,wind_asymmetry, &
                   viscosity_m2_s,slip_length_km,snapshot_days]
  call mpi_bcast(settings,size(settings),MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(real_settings,size(real_settings),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(wind_profile,len(wind_profile),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(advection,len(advection),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  nx=settings(1)
  ny=settings(2)
  cape_start=settings(3)
  cape_width=settings(4)
  cape_height=settings(5)
  basin_km=real_settings(1)
  run_days=real_settings(2)
  max_dt_seconds=real_settings(3)
  beta_gradient=real_settings(4)
  wind_stress=real_settings(5)
  wind_shift=real_settings(6)
  wind_tilt=real_settings(7)
  wind_asymmetry=real_settings(8)
  viscosity_m2_s=real_settings(9)
  slip_length_km=real_settings(10)
  snapshot_days=real_settings(11)

  if (.not. all(ieee_is_finite(real_settings)) .or. min(basin_km,run_days,max_dt_seconds)<=0.0D0) &
    call check('basin_km, run_days and max_dt_seconds must be finite and positive, and all parameters finite')
  select case(wind_profile)
    case('cape_single')
      ! Bounds are set from the cape geometry below.
    case('double_gyre')
      if (abs(wind_shift)+0.5D0*abs(wind_tilt)>=0.5D0 .or. wind_asymmetry<=0.0D0) &
        call check('double_gyre requires abs(wind_shift) + abs(wind_tilt)/2 < 0.5 and wind_asymmetry > 0')
    case default
      call check("wind_profile must be 'cape_single' or 'double_gyre'")
  end select
  if (viscosity_m2_s<0.0D0 .or. (viscosity_m2_s>0.0D0 .and. slip_length_km<=0.0D0)) &
    call check('viscosity_m2_s must be nonnegative and viscous runs require slip_length_km > 0')
  if (snapshot_days<=0.0D0 .or. snapshot_days<run_days/real(huge(snapshot_index)-1,real_8)) &
    call check('snapshot_days must be positive and the snapshot count must fit in a default integer')
  viscous=viscosity_m2_s>0.0D0
  if (advection/='upwind' .and. advection/='cabaret') call check("advection must be 'upwind' or 'cabaret'")
  use_cabaret=advection=='cabaret'
  grid=0
  call mpi_dims_create(ranks,2,grid,ierr)
  if (ny/grid(1)<3 .or. nx/grid(2)<3) call check('cape decomposition requires at least three nodes per rank direction')
  if (cape_start<1 .or. cape_start>=nx-2) call check('cape_start must be an interior offset from the west wall')
  h_grid=basin_km*1.0D5/(ny-1)
  time_scale=h_grid/u_scale
  end_time=run_days*day_s/time_scale
  max_dt=max_dt_seconds/time_scale
  beta_grid=beta_gradient*h_grid*time_scale
  if (.not. all(ieee_is_finite([h_grid,time_scale,end_time,max_dt,beta_grid])) .or. &
      min(h_grid,time_scale,end_time,max_dt)<=0.0D0) call check('invalid derived cape scales')
  if (rank==0) call write_rectangular_cape(geometry_path,ny,nx,1+cape_start,cape_width,cape_height,h_grid,error)
  call broadcast_error(MPI_COMM_WORLD,rank,ranks,error)
  call check(error)
  iis = [(ny/grid(1)+merge(1,0,i>grid(1)-mod(ny,grid(1))),i=1,grid(1))]
  jjs = [(nx/grid(2)+merge(1,0,j>grid(2)-mod(nx,grid(2))),j=1,grid(2))]
  call allocate(decomp,iis=iis,jjs=jjs,layers=layers)
  call load_distributed_geometry(geometry_path,ny,nx,[decomp%ii_0,decomp%jj_0,decomp%ii,decomp%jj], &
                                 decomp%cart_comm,h_grid,geom,error)
  call check(error)
  call allocate(strat,rd_bcl1,rd_bcl2,depths(1),depths(2),depths(3),l_scale=h_grid)
  call setup_basin_mass(inverse,decomp,geom,strat,error)
  call check(error)
  if (viscous) then
    call setup_cape_slip(slip,decomp,geom,strat,inverse,viscosity_m2_s*(1.0D4/(h_grid*u_scale)), &
                         slip_length_km*(1.0D5/h_grid),error)
    call check(error)
  end if
  if (viscous .or. use_cabaret) then
    allocate(diff_old(geom%ni,geom%nj,layers),diff_new(geom%ni,geom%nj,layers), &
             step_rhs(geom%ni,geom%nj,layers),next_state(geom%ni,geom%nj,layers))
    diff_new=0.0D0
  end if
  call setup_vertex_plan(geom,vertices,error,normal_profile=use_cabaret)
  call check(error)

  if (use_cabaret) then
    call setup_cut_cabaret(geom,cabaret,error,include_full=.true.)
    call check(error)
  else
    n=count(geom%faces%owner==geom%rank)
    allocate(owned_faces(n),requested(2,2*n),pv_at_faces(layers,2*n))
    k=0
    do f=1,size(geom%faces)
      if (geom%faces(f)%owner/=geom%rank) cycle
      k=k+1
      owned_faces(k)=f
      requested(:,2*k-1:2*k)=geom%faces(f)%cells
    end do
    call build_donor_plan(geom,requested,face_cells,error,cell_means=.true.)
    call check(error)
  end if
  allocate(q(geom%ni,geom%nj,layers),predicted(geom%ni,geom%nj,layers),rhs0(geom%ni,geom%nj,layers), &
           rhs1(geom%ni,geom%nj,layers),psi(geom%ni,geom%nj,layers),wind(geom%ni,geom%nj,layers), &
           planetary(geom%ni,geom%nj,layers),volume(layers,size(geom%faces)), &
           velocity(2,layers,size(geom%faces)),trace(2,2,layers,size(geom%faces)),flux(layers,size(geom%faces)))
  if (wind_profile=='cape_single') then
    wind_west=real(cape_start+cape_width,real_8)
    wind_south=real(ny-1-cape_height,real_8)
    strength=wind_stress/(u_scale**2*strat%h(1))
    call cape_wind_curl_torque(geom,strength,wind_west,wind_south,wind,torque,error)
  else
    strength=wind_stress/(u_scale**2*(ny-1)*strat%h(1))
    call wind_curl_torque(geom,strength,wind_shift,wind_tilt,wind_asymmetry,wind,torque,error)
  end if
  call check(error)
  if (geom%rank==0) then
    print '(a,1x,a)', 'Advection',trim(advection)
    print '(a,1x,a)', 'Wind profile',trim(wind_profile)
    if (wind_profile=='cape_single') print '(a,4(1x,f12.3))', 'Wind patch west east south north (km)', &
      [wind_west,real(nx-1,real_8),wind_south,real(ny-1,real_8)]*h_grid/1.0D5
  end if
  call source_momentum_budget(geom,wind,torque,budget,error)
  call check(error)
  call mpi_allreduce(sum(sum(wind,dim=1),dim=1),wind_total,layers,MPI_DOUBLE_PRECISION,MPI_SUM,geom%comm,ierr)

  planetary=0.0D0
  do j=1,geom%nj
    do i=1,geom%ni
      if (geom%area(i,j)>0.0D0) &
        planetary(i,j,:)=beta_grid*(0.5D0*(ny+1)-geom%i0-i+geom%moments(3,i,j))
    end do
  end do
  local_area=sum(geom%area)
  call mpi_allreduce(local_area,wet_area,1,MPI_DOUBLE_PRECISION,MPI_SUM,geom%comm,ierr)

  q=0.0D0
  target_mass=0.0D0
  time=0.0D0
  peak_cfl=0.0D0
  peak_slip_defect=0.0D0
  local_wall_pv=0.0D0
  local_solver_pv=0.0D0
  divergence_error=0.0D0
  steps=0
  previous_dt=0.0D0
  peak_fallback=0
  rejected_steps=0
  snapshot_index=0
  call evaluate_state(q)
  if (use_cabaret) call initialize_cabaret()
  call write_snapshot()
  next_output=min(run_days,snapshot_days)*day_s/time_scale
  do while (time<end_time)
    if (viscous .or. use_cabaret) diff_old=diff_new
    if (use_cabaret) then
      call face_tendency(velocity,trace,rhs0,rate0)
    else
      call tendency(q,rhs0,rate0)
    end if
    remaining=next_output-time
    dt=min(max_dt,remaining)
    if (remaining-dt<=8*spacing(next_output)) dt=remaining
    if (rate0>0.0D0) dt=min(dt,cfl_target/rate0)
    retries=0
    do
      if (time+dt<=time) call check('cape time step underflow')
      if (use_cabaret) then
        call predict_cabaret(dt)
      else
        if (viscous) then
          step_rhs=q+dt*rhs0
          call cape_implicit_step(slip,dt,step_rhs,predicted,error)
          call check(error)
        else
          predicted=q+dt*rhs0
        end if
        call evaluate_state(predicted)
        call tendency(predicted,rhs1,rate1)
      end if
      if (dt*rate1<=cfl_target) exit
      dt=0.9D0*cfl_target/rate1
      retries=retries+1
      if (retries>12) call check('cape corrector CFL could not be satisfied')
    end do
    rejected_steps=rejected_steps+retries
    if (use_cabaret) then
      call correct_cabaret(dt)
    else if (viscous) then
      ! Heun advection with a backward-Euler predictor and Crank-Nicolson diffusion corrector.
      step_rhs=q+0.5D0*dt*(rhs0+rhs1+diff_old)
      call cape_implicit_step(slip,0.5D0*dt,step_rhs,next_state,error)
      call check(error)
      call evaluate_state(next_state)
      call record_pv_residual(next_state-q-0.5D0*dt*(rhs0+rhs1+diff_old+diff_new))
      q=next_state
    else
      q=q+0.5D0*dt*(rhs0+rhs1)
      call evaluate_state(q)
    end if
    steps=steps+1
    peak_cfl=max(peak_cfl,dt*max(rate0,rate1))
    if (dt>=remaining) then
      time=next_output
      snapshot_index=snapshot_index+1
      call write_snapshot()
      if (time>=end_time) exit
      next_output=min(run_days,real(snapshot_index+1,real_8)*snapshot_days)*day_s/time_scale
      if (end_time-next_output<=8*spacing(end_time)) next_output=end_time
    else
      time=time+dt
    end if
  end do

  call z_to_modes(strat,psi,predicted)
  call fluid_integrals(geom,predicted,mass)
  call fluid_integrals(geom,q,pv_total)
  call fluid_integrals(geom,abs(q),pv_scale)
  call mpi_allreduce(local_wall_pv,wall_pv,layers,MPI_DOUBLE_PRECISION,MPI_SUM,geom%comm,ierr)
  budget=pv_total/h_grid**2-time*wind_total-wall_pv
  budget_scale=max(1.0D-20,pv_scale/h_grid**2,time*abs(wind_total),abs(wall_pv))
  if (viscous) call mpi_allreduce(local_solver_pv,solver_pv,layers,MPI_DOUBLE_PRECISION,MPI_SUM,geom%comm,ierr)
  if (any(abs(budget)>2.0D-11*budget_scale)) then
    if (geom%rank==0) then
      print '(a,3(1x,es14.6))', 'Integrated PV budget residual by layer',budget
      print '(a,3(1x,es14.6))', 'Integrated PV budget scales by layer',budget_scale
      if (viscous) print '(a,3(1x,es14.6))', 'Accumulated strict-cell implicit residual by layer',solver_pv
      flush(output_unit)
    end if
    call check('cape integrated PV budget failed')
  end if
  local_psi=maxval(abs(psi),mask=spread(geom%area>0.0D0,3,layers))
  call mpi_allreduce(local_psi,max_psi,1,MPI_DOUBLE_PRECISION,MPI_MAX,geom%comm,ierr)
  call write_state()
  if (geom%rank==0) then
    print '(a,1x,a)', 'Geometry SHA256',geom%global%digest
    print '(a,1x,f12.3)', 'Wet area (km^2)',wet_area/1.0D10
    print '(a,1x,i0,1x,a,1x,f10.5)', 'Completed steps',steps,'days',time*time_scale/day_s
    print '(a,1x,es14.6)', 'Max outgoing CFL',peak_cfl
    if (use_cabaret) then
      print '(a,1x,i0)', 'Rejected CABARET steps',rejected_steps
      print '(a,1x,i0)', 'Peak CABARET rank-deficient cut traces per step',peak_fallback
    end if
    print '(a,1x,es14.6)', 'Max open-face transport divergence (native)',divergence_error
    print '(a,1x,es14.6)', 'Max baroclinic area-mean residual (native)',maxval(abs(mass(2:)))/wet_area
    print '(a,1x,es14.6)', 'Max integrated PV budget residual (native)',maxval(abs(budget))
    if (viscous) then
      print '(a,1x,es14.6)', 'Max discrete Navier wall residual (native)',peak_slip_defect
      print '(a,3(1x,es14.6))', 'Diagnosed wall PV input by layer (native)',wall_pv
      if (use_cabaret) print '(a,3(1x,es14.6))', 'Accumulated strict-cell implicit residual by layer (native)',solver_pv
    end if
    print '(a,1x,es14.6)', 'Max wet streamfunction (m^2/s)',max_psi*h_grid*u_scale/1.0D4
  end if
  if (viscous) call release_cape_slip(slip)
  call release_basin_mass(inverse)
  call deallocate(strat)
  call deallocate(decomp)
  call parallel_finalise()

contains

  subroutine evaluate_state(state)
    real(real_8), intent(inout) :: state(:,:,:)

    if (viscous) then
      call cape_viscous_state(slip,state,psi,wall_values,diff_new,slip_defect,error)
      call check(error)
      peak_slip_defect=max(peak_slip_defect,slip_defect)
    else
      call invert_basin_mass(inverse,geom%comm,state,target_mass,psi,wall_values(0,:),error)
      call check(error)
      where (spread(geom%area==0.0D0,3,layers)) psi=0.0D0
    end if
  end subroutine evaluate_state

  subroutine initialize_cabaret()
    integer :: face,endpoint

    allocate(wind_rate,mold=q)
    allocate(diff_half,source_in,source_out,mold=q)
    allocate(speed_half,previous_half,speed_next,mold=velocity)
    allocate(trace_next,mold=trace)
    allocate(cell_velocity(2,geom%ni,geom%nj,layers))
    flux=0.0D0
    call wet_cell_rate(geom,flux,wind,wind_rate,error)
    call check(error)
    velocity=0.0D0
    ! Relative PV and velocity are initially zero. Initialize total PV
    ! at each endpoint, including the exact affine planetary term.
    do face=1,size(geom%faces)
      do endpoint=1,2
        trace(endpoint,:,:,face)=beta_grid*(0.5D0*geom%faces(face)%endpoints(2,endpoint)-0.5D0*(ny-1))
      end do
    end do
  end subroutine initialize_cabaret

  subroutine cabaret_stage(state,stage_dt,forcing,diffusion,result)
    real(real_8), intent(in) :: state(:,:,:),stage_dt,forcing(:,:,:),diffusion(:,:,:)
    real(real_8), intent(out) :: result(:,:,:)

    if (viscous) then
      step_rhs=state+stage_dt*(forcing+0.5D0*diffusion)
      ! Two implicit solves per step can accumulate an integral bias over
      ! long runs. Use a tighter residual for the CABARET half steps.
      call cape_implicit_step(slip,0.5D0*stage_dt,step_rhs,result,error,relative_tolerance=2.0D-14)
      call check(error)
    else
      result=state+stage_dt*forcing
    end if
    call evaluate_state(result)
  end subroutine cabaret_stage

  subroutine predict_cabaret(step)
    real(real_8), intent(in) :: step
    real(real_8) :: half_step,weight
    real(real_8), allocatable :: packet(:,:)

    half_step=0.5D0*step
    call cabaret_stage(q,half_step,rhs0,diff_old,predicted)
    diff_half=diff_new
    call instantaneous_transport(geom,vertices,psi,wall_values,volume,error,speed_half)
    call check(error)
    call check_transport()
    packet=reshape(speed_half,[2*layers,size(geom%faces)])
    call share_face_flux(geom,packet,error)
    call check(error)
    speed_half=reshape(packet,shape(speed_half))
    call cabaret_cell_velocity(geom,speed_half,cell_velocity)
    if (previous_dt>0.0D0) then
      weight=step/(previous_dt+step)
      speed_next=speed_half+weight*(speed_half-previous_half)
    else
      speed_next=2*speed_half-velocity
    end if

    ! Relative to the midpoint, subtract the first half's sources from
    ! incoming traces. The inferred cut-cell source includes the slip
    ! projection. Extrapolate the smooth strict-cell diffusive rate to
    ! the outgoing time and hold the diagnosed cut-cell rate for a half step.
    source_in=q+half_step*(rhs0-wind_rate)-predicted
    source_out=-source_in
    if (viscous) then
      where (spread(geom%area==h_grid**2,3,layers)) source_out=source_out+half_step*(diff_half-diff_old)
    end if
    call cabaret_face_traces(geom,cabaret,predicted+planetary,trace,velocity,cell_velocity,source_in,source_out, &
                             step,.true.,trace_next,fallback,error)
    call check(error)
    call face_tendency(speed_next,trace_next,rhs1,rate1)
  end subroutine predict_cabaret

  subroutine correct_cabaret(step)
    real(real_8), intent(in) :: step
    real(real_8) :: half_step

    half_step=0.5D0*step
    call cabaret_stage(predicted,half_step,rhs1,diff_half,next_state)
    if (viscous) call record_pv_residual(next_state-q-half_step*(rhs0+rhs1+diff_half+0.5D0*(diff_old+diff_new)))
    q=next_state
    velocity=speed_next
    trace=trace_next
    previous_half=speed_half
    previous_dt=step
    peak_fallback=max(peak_fallback,fallback)
  end subroutine correct_cabaret

  subroutine record_pv_residual(residual)
    real(real_8), intent(in) :: residual(:,:,:)
    integer :: row,col

    do col=1,geom%nj
      do row=1,geom%ni
        select case(geom%global%kind(geom%i0+row,geom%j0+col))
          case(wall_node)
            local_wall_pv=local_wall_pv+(geom%area(row,col)/h_grid**2)*residual(row,col,:)
          case(fluid_node)
            local_solver_pv=local_solver_pv+residual(row,col,:)
        end select
      end do
    end do
  end subroutine record_pv_residual

  subroutine face_tendency(speed,pv_trace,rate,outgoing)
    real(real_8), intent(in) :: speed(:,:,:),pv_trace(:,:,:,:)
    real(real_8), intent(out) :: rate(:,:,:),outgoing

    call transport_fluxes(geom,speed,pv_trace,flux,outgoing,error)
    call check(error)
    call wet_cell_rate(geom,flux,wind,rate,error)
    call check(error)
  end subroutine face_tendency

  subroutine check_transport()
    real(real_8) :: divergence(geom%ni,geom%nj,layers),stage_divergence

    call sum_open_flux(geom,volume,divergence)
    local_divergence=maxval(abs(divergence))
    call mpi_allreduce(local_divergence,stage_divergence,1,MPI_DOUBLE_PRECISION,MPI_MAX,geom%comm,ierr)
    divergence_error=max(divergence_error,stage_divergence)
    if (.not. all(ieee_is_finite(volume))) call check('nonfinite cape face transport')
  end subroutine check_transport

  subroutine tendency(state,rate,outgoing)
    real(real_8), intent(in) :: state(:,:,:)
    real(real_8), intent(out) :: rate(:,:,:),outgoing
    integer :: face,slot

    call instantaneous_transport(geom,vertices,psi,wall_values,volume,error)
    call check(error)
    call check_transport()
    call gather_donors(face_cells,state+planetary,pv_at_faces)
    velocity=0.0D0
    trace=0.0D0
    do slot=1,size(owned_faces)
      face=owned_faces(slot)
      length=0.5D0*sum(abs(geom%faces(face)%endpoints(:,2)-geom%faces(face)%endpoints(:,1)))
      velocity(1,:,face)=volume(:,face)/length
      velocity(2,:,face)=velocity(1,:,face)
      trace(1,1,:,face)=pv_at_faces(:,2*slot-1)
      trace(1,2,:,face)=pv_at_faces(:,2*slot)
      trace(2,:,:,face)=trace(1,:,:,face)
    end do
    call face_tendency(velocity,trace,rate,outgoing)
  end subroutine tendency

  subroutine write_snapshot()
    real(real_8) :: day
    integer :: layer

    day=time*time_scale/day_s
    if (time==end_time) day=run_days
    call write_s_header(decomp,snapshot_index,h_grid,time_scale,single_level=.true.,overwrite=.false.)
    call open_parallel_data(decomp,io_unit,'qg3_'//int2str(ny)//'_s',output_index=snapshot_index,overwrite=.false.)
    call write_unformatted(io_unit,day)
    call write_unformatted(io_unit,day)
    do layer=1,layers
      call write_unformatted(io_unit,psi(:,:,layer))
    end do
    do layer=1,layers
      call write_unformatted(io_unit,q(:,:,layer))
    end do
    close(io_unit)
  end subroutine write_snapshot

  subroutine write_state()
    type(parallel_data) :: field
    type(request_cache), target :: requests
    real(real_8), allocatable :: stream_global(:,:,:),q_global(:,:,:)
    integer :: row,col,layer

    call allocate(requests,comm=geom%comm)
    call allocate(field,decomp,requests=requests)
    field%sdata=psi
    if (geom%rank==0) then
      allocate(stream_global(ny,nx,layers),q_global(ny,nx,layers))
      call gather(field,stream_global,requests=requests)
    else
      call gather(field,1,requests=requests)
    end if
    field%sdata=q
    if (geom%rank==0) then
      call gather(field,q_global,requests=requests)
    else
      call gather(field,1,requests=requests)
    end if
    call deallocate(field)
    call deallocate(requests)

    if (geom%rank/=0) return
    open(newunit=unit,file=state_path,status='new',action='write',iostat=status,iomsg=error)
    call check(error)
    write(unit,'(a)',iostat=status,iomsg=error) 'i,j,layer,kind,x_km,y_km,wet_fraction,psi_m2_s,relative_pv_s_1'
    call check(error)
    do layer=1,layers
      do row=1,ny
        do col=1,nx
          if (geom%global%kind(row,col)==solid_node) cycle
          write(unit,'(3(i0,","),i0,5(",",es24.16e3))',iostat=status,iomsg=error) &
            row,col,layer,geom%global%kind(row,col), &
            (col-1)*h_grid/1.0D5,(ny-row)*h_grid/1.0D5,popcnt(geom%global%wedge(row,col))/8.0D0, &
            stream_global(row,col,layer)*h_grid*u_scale/1.0D4,q_global(row,col,layer)/time_scale
          call check(error)
        end do
      end do
    end do
    close(unit,iostat=status,iomsg=error)
    call check(error)
  end subroutine write_state

  subroutine check(message)
    character(*), intent(in) :: message

    if (len_trim(message)/=0) call pabort(trim(message),__FILE__,__LINE__)
  end subroutine check
end program qg3_cape
