#include "config.h"
#include "leapfrog.h"

! Multi-layer QG model
program qg3

  use qg_aux
  use qg_diagnostics
  use qg_helmholtz
  use qg_io
  use qg_laplace_bc
  use qg_parallel
  use qg_parameters, only : alpha_dim => alpha, beta_dim => beta, &
    & dt_dim => dt, max_spin_dim => max_spin, max_time_dim => max_time, &
    & start_time_override_dim => start_time_override, &
    & visc_bot_dim => visc_bot, viscosity_dim => visc, &
    & zeta_perturb_dim => zeta_perturb
  use qg_parameters, only : basinscale, day_s, h, ii_m, iis, jj_m, jjs, &
    & k_cp, k_dump, k_out, k_save, layers, partial_slip, periodic, pi, quad, &
    & rd_bcl1, rd_bcl2, real_8, real_8_chr, save_log, set_start_time, &
    & set_strat, start, START_PERTURB, START_PICKUP, START_TIME_AUTO, &
    & START_TIME_MANUAL, START_ZERO, STRAT_RD, STRAT_INPUT, u_scale, unit
#ifdef WIND_FORCING
  use qg_parameters, only : tau_dim => tau
  use qg_parameters, only : wind_asym, wind_shift, wind_tilt
#endif
#ifdef BACKGROUND_VSHEAR
  use qg_parameters, only : u_background_v_dim => u
  use qg_parameters, only : set_shear, SHEAR_U, SHEAR_INPUT
#endif
#ifdef CABARET
  use qg_parameters, only : max_dt_dim => max_dt
  use qg_parameters, only : adapt_dt, max_dt_increase, target_cfl
#endif

  implicit none

  include "mpif.h"

  integer :: i, j, layer
  integer :: day_number, dump_count, last_diagnostics_day, last_dump_day, &
    & s_index
  logical :: exit_timestep_loop
  real(kind = real_8) :: diagnostic_time, model_time, model_time_day, &
    & start_time, start_time_day

  real(kind = real_8) :: cff, l_scale, t_scale
  real(kind = real_8) :: alpha, dt, max_spin, max_time, visc_bot, viscosity, &
    &  zeta_perturb
  real(kind = real_8), dimension(layers) :: beta
#ifdef WIND_FORCING
  real(kind = real_8) :: tau
#endif
#ifdef BACKGROUND_VSHEAR
  real(kind = real_8), dimension(layers) :: u_background_v
#endif
#ifdef CABARET
  real(kind = real_8) :: cfl, max_dt
#endif

  real(kind = real_8), dimension(ii_m + 2, layers) :: p_beta_y
  type(stratification) :: strat
  type(velocity_boundary_condition) :: vbc
#ifdef WIND_FORCING
  real(kind = real_8) :: a_wind, b_wind, x, y, y_mid
  real(kind = real_8), dimension(ii_m, jj_m) :: force_wind
#endif

#ifdef CABARET
  !!!!!! CABARET data
  ! Time integration weights
  real(kind = real_8) :: tweight_a, tweight_u
  ! Cell faces
  type(parallel_flux_data) :: zeta_fi, zeta_fj
  type(parallel_flux_data) :: ui, uj
  type(parallel_flux_data) :: uimo, ujmo
  type(parallel_flux_data) :: flux_i, flux_j
  ! Cell residuals
  real(kind = real_8), dimension(ii_m, jj_m, layers) :: res, res_beta
#else
  !!!!!! Leapfrog data
  integer :: leapfrog
  real(kind = real_8) :: T
#endif
  real(kind = real_8), dimension(ii_m + 2, jj_m + 2, layers TNLEAPFROG) :: zeta_old

  !!!!!! Diagnostics
  ! Energy
  real(kind = real_8) :: epot
  real(kind = real_8), dimension(layers) :: ekin

  ! Time averages
  real(kind = real_8), dimension(ii_m, jj_m, layers) :: psi_av
  ! Covariances
  real(kind = real_8), dimension(:, :, :), allocatable :: u, v, u_i, v_i, &
    & tness
  type(eddy_covariance) :: uu_av, vv_av, uu_i_av, vv_i_av, uv_av, ur_av, &
    & vr_av, rr_av, pr_av, uz_av, vz_av, zz_av, zr_av, pp_av, pz_av, ut_av, &
    & vt_av, rt_av, tt_av, sc_ub_av, sc_vb_av, sc_bb_av, sc_bb_i_av
  ! Timers
  real(kind = real_8) :: dtime, htime, itime, lwtime_0, wtime, wtime_0

  !!!!!! Helmholtz solver data
  type(helmholtz_cache) :: hcache

  !!!!!! Parallel data
  type(decomposition) :: decomp
  type(parallel_data) DIMLEAPFROG, target :: psi, rel, zeta
#ifdef BACKGROUND_HSHEAR
  type(parallel_data) :: psi_b
  type(parallel_data), target :: rel_b
  type(parallel_data), pointer :: zeta_b
  type(parallel_data) DIMLEAPFROG :: rel_p
#else
  type(parallel_data) DIMSLEAPFROG, pointer :: rel_p
#endif

  !!!!!! Parallel initialisation

  call parallel_init()
  if(save_log) call log_output("qg" // int2str(layers))
  call print_config()
  call check_decomp()
  call allocate(decomp, iis = iis, jjs = jjs, layers = layers, periodic_j = periodic)
  call print(decomp)
  if(decomp%n > 1 .or. decomp%periodic_j) call test_pdata(decomp)

  !!!!!! Seed the random number generator

  call rseed(42)

  !!!!!! Allocate prognostic variables

  DOLEAPFROG
  call allocate(psi BLEAPFROG, decomp)
  call allocate(rel BLEAPFROG, decomp)
  call allocate(zeta BLEAPFROG, decomp)
  ENDDOLEAPFROG
#ifdef CABARET
  call allocate(zeta_fi, decomp, FLUX_II)
  call allocate(zeta_fj, decomp, FLUX_JJ)
  call allocate(ui, decomp, FLUX_II)
  call allocate(uj, decomp, FLUX_JJ)
  call allocate(uimo, decomp, FLUX_II)
  call allocate(ujmo, decomp, FLUX_JJ)
  call allocate(flux_i, decomp, FLUX_II)
  call allocate(flux_j, decomp, FLUX_JJ)
#endif

  !!!!!! Non-dimensionalisation

  l_scale = basinscale / real(decomp%g_ii_int, kind = real_8)
  t_scale = l_scale / u_scale

  dt = dt_dim / t_scale
  max_spin = max_spin_dim * day_s / t_scale
  max_time = max_time_dim * day_s / t_scale
#ifdef CABARET
  if(adapt_dt) max_dt = max_dt_dim / t_scale
#endif
  alpha = l_scale / alpha_dim
  beta = beta_dim * l_scale * l_scale / u_scale
#ifdef WIND_FORCING
  tau = tau_dim / (u_scale * u_scale * real(decomp%g_ii_int, kind = real_8))
#endif
  viscosity = viscosity_dim / (l_scale * u_scale)
  visc_bot = visc_bot_dim * l_scale / u_scale
  zeta_perturb = zeta_perturb_dim * t_scale
  
#ifdef BACKGROUND_VSHEAR
  !!!!!! Background vertical shear

  select case(set_shear)
    case(SHEAR_U)
      u_background_v = u_background_v_dim / u_scale
    case(SHEAR_INPUT)
      call read_shear_v(u_background_v, u_scale)
    case default
      write(6, "(a,i0)") "For shear option: ", set_shear
      labort("Unrecognised shear option")
  end select
  write(6, "(a)") "Background velocity (dimensional)"
  do layer = 1, layers
    write(6, "(a,i0,a,"//real_8_chr//")") "Layer ", layer, ": ", u_background_v(layer) * u_scale
  end do
  write(6, "(a)") "Background velocity (non-dimensional)"
  do layer = 1, layers
    write(6, "(a,i0,a,"//real_8_chr//")") "Layer ", layer, ": ", u_background_v(layer)
  end do
  call write_shear_v(u_background_v, u_scale)
#endif

#ifdef BACKGROUND_HSHEAR
  !!!!!! Background horizontal shear
  
  call allocate(psi_b, decomp)
  call allocate(rel_b, decomp)
  zeta_b => rel_b
  call initialise_background(psi_b, rel_b)
  call update(rel_b)
  call flush()
  DOLEAPFROG
  call allocate(rel_p BLEAPFROG, decomp)
  ENDDOLEAPFROG
#else
  rel_p => rel
#endif

  !!!!!! Stratification

  select case(set_strat)
    case(STRAT_RD)
      select case(layers)
        case(3)
          call allocate(strat, rd_bcl1, rd_bcl2, h(1), h(2), h(3), l_scale = l_scale)
        case(2)
          call allocate(strat, rd_bcl1, h(1), h(2), l_scale = l_scale)
        case(1)
          call allocate(strat, h(1), l_scale = l_scale)
        case default
          labort("Unable to allocate stratification with " // int2str(layers) // " layers")
      end select
    case(STRAT_INPUT)
      call read_stratification(strat, l_scale = l_scale)
      if(strat%layers /= layers) then
        write(6, "(a,i0)") "Number of layers: ", strat%layers
        labort("Invalid number of layers")
      end if
    case default
      write(6, "(a,i0)") "For stratification option: ", set_strat
      labort("Unrecognised stratification option")
  end select
  call print(strat)
  call write_stratification(strat)
#ifdef WIND_FORCING
  tau = tau / strat%h(1)
#endif

  !!!!!! Boundary condition

  if(partial_slip) then
    if(alpha_dim >= 1.0D100) then
      call initialise(vbc, BC_FREE_SLIP)
    else
      call initialise(vbc, BC_PARTIAL_SLIP, alpha = alpha)
    end if
  else
    call initialise(vbc, BC_NO_SLIP)
  end if

  !!!!!! Wind stress

#ifdef WIND_FORCING
  a_wind = -1.0D0 * wind_asym
  b_wind = 1.0D0 / wind_asym

  do j = 1, decomp%jj
     do i = 1, decomp%ii
        x = real(decomp%jj_0 + j - 1, kind = real_8) / real(decomp%g_jj - 1, kind = real_8) - 0.5D0
        y = -real(decomp%ii_0 + i - 1, kind = real_8) / real(decomp%g_ii - 1, kind = real_8) + 0.5D0

        y_mid = wind_shift + wind_tilt * x

        if(y < y_mid)then
          force_wind(i, j) = a_wind * dsin(pi * (0.5D0 + y) / (0.5D0 + y_mid))
        else
          force_wind(i, j) = b_wind * dsin(pi * (y - y_mid) / (0.5D0 - y_mid))
        end if
        force_wind(i, j) =  tau * 2.0D0 * pi * force_wind(i, j)
     end do
  end do
!  call open_parallel_header(decomp, unit, filename = "wind_stress")
!  call write_parallel_header(decomp, unit, l_scale, t_scale)
!  write(unit, "(a)") "force_wind  ii jj"
!  close(unit)
!  open(unit, file = "wind_stress.d", access = "stream", form = "unformatted", status = "replace")
!  call write_unformatted(unit, force_wind)
!  close(unit)
#endif

  !!!!!! Background PV gradient

#ifdef BACKGROUND_VSHEAR
  if(layers > 1) then
    beta(1) = beta(1) + strat%s2(1) * (u_background_v(1) - u_background_v(2))
    do layer = 2, layers - 1
      beta(layer) = beta(layer) + &
        & strat%s1(layer) * (u_background_v(layer) - u_background_v(layer - 1)) + &
        & strat%s2(layer) * (u_background_v(layer) - u_background_v(layer + 1))
    end do
    beta(layers) = beta(layers) + &
      & strat%s1(layers) * (u_background_v(layers) - u_background_v(layers - 1))
  end if
#endif
  write(6, "(a)") "Beta (dimensional):"
  do layer = 1, layers
    write(6, "(a,i0,a,"//real_8_chr//")") "Layer ", layer, ": ", beta(layer) / (l_scale * t_scale)
  end do
  write(6, "(a)") "Beta (non-dimensional):"
  do layer = 1, layers
    write(6, "(a,i0,a,"//real_8_chr//")") "Layer ", layer, ": ", beta(layer)
  end do
  do layer = 1, layers
    do i = 0, decomp%ii + 1
      p_beta_y(i + 1, layer) = beta(layer) * ((decomp%g_ii + 1.0D0) / 2.0D0 - (decomp%ii_0 + i))
    end do
  end do

  !!!!!! Elliptic solver initialisation

  call allocate(hcache, decomp, strat, vbc)

  !!!!!! Prognostic variable initialisation

  select case(start)
    case(START_PICKUP)
      write(6, "(a)") "Restarting from checkpoint"
#ifdef CABARET
      ! Read the stream function and relative potential vorticity, and other
      ! non-recoverable data
      call read_parallel_checkpoint(start_time_day, psi, zeta, dt, ui, uj, zeta_fi, zeta_fj, res, res_beta)
      start_time = start_time_day * day_s / t_scale
      last_dump_day = int(start_time_day + 1.0D-3)
      last_diagnostics_day = last_dump_day
      call update(psi)
      !call update(zeta)
      call update(zeta_fi)
      call update(zeta_fj)
      call flush()
      ! Update flux variables
      call cabaret_cell_boundaries(psi, uimo, ujmo &
#ifdef BACKGROUND_VSHEAR
        & , u_background_v &
#endif
        & )
      call cabaret_flux(ui, uj, zeta_fi, zeta_fj, flux_i, flux_j)
      ! Take one half timestep forwards
      call cabaret_corrector(dt, flux_i, flux_j, res, zeta)
      call update(zeta)
      call flush()
#else
      ! Read the stream function and relative potential vorticity
      call read_parallel_checkpoint(start_time_day, psi, zeta)
      select case(set_start_time)
        case(START_TIME_AUTO)
        case(START_TIME_MANUAL)
          start_time_day = start_time_override_dim
        case default
          write(6, "(a,i0)") "For start time option: ", set_start_time
          labort("Unrecognised start time option")
      end select
      start_time = start_time_day * day_s / t_scale
      last_dump_day = int(start_time_day + 1.0D-3)
      last_diagnostics_day = last_dump_day
      call update(psi)
      call update(zeta)
      call flush()
      ! Update the relative vorticity
      DOLEAPFROG
      call rel_from_psi(vbc, psi BLEAPFROG, rel BLEAPFROG)
      call update(rel BLEAPFROG)
      ENDDOLEAPFROG
      call flush()
#endif
    case(START_PERTURB)
      write(6, "(a)") "Starting with perturbation"
      select case(set_start_time)
        case(START_TIME_AUTO)
          start_time_day = 0.0D0
          start_time = 0.0D0
          last_dump_day = 0
          last_diagnostics_day = 0
        case(START_TIME_MANUAL)
          start_time_day = start_time_override_dim
          start_time = start_time_day * day_s / t_scale
          last_dump_day = int(start_time_day + 1.0D-3)
          last_diagnostics_day = last_dump_day
        case default
          write(6, "(a,i0)") "For start time option: ", set_start_time
          labort("Unrecognised start time option")
      end select
      call init_rand(zeta, rscale = zeta_perturb, root = 1)
      DOLEAPFROG
      call apply_dirichlet_bc(0.0D0, zeta BLEAPFROG)
      call update(zeta BLEAPFROG)
      call flush()
      call pv_inversion(strat, vbc, zeta BLEAPFROG, hcache, psi BLEAPFROG)
      call update(psi BLEAPFROG)
      call flush()
      call rel_from_zeta_and_psi(strat, vbc, zeta BLEAPFROG, psi BLEAPFROG, rel BLEAPFROG)
      call update_zeta_bc(strat, rel BLEAPFROG, psi BLEAPFROG, zeta BLEAPFROG)
      ENDDOLEAPFROG
#ifdef CABARET
      call zero(zeta_fi)
      call zero(zeta_fj)
      call zero(ui)
      call zero(uj)
      call zero(uimo)
      call zero(ujmo)
      call zero(flux_i)
      call zero(flux_j)
      res(:decomp%ii, :decomp%jj, :) = 0.0D0
      res_beta(:decomp%ii, :decomp%jj, :) = 0.0D0
#endif
    case(START_ZERO)
      write(6, "(a)") "Starting without perturbation"
      select case(set_start_time)
        case(START_TIME_AUTO)
          start_time_day = 0.0D0
          start_time = 0.0D0
          last_dump_day = 0
          last_diagnostics_day = 0
        case(START_TIME_MANUAL)
          start_time_day = start_time_override_dim
          start_time = start_time_day * day_s / t_scale
          last_dump_day = int(start_time_day + 1.0D-3)
          last_diagnostics_day = last_dump_day
        case default
          write(6, "(a,i0)") "For start time option: ", set_start_time
          labort("Unrecognised start time option")
      end select
      call zero(psi)
      call zero(rel)
      call zero(zeta)
#ifdef CABARET
      call zero(zeta_fi)
      call zero(zeta_fj)
      call zero(ui)
      call zero(uj)
      call zero(uimo)
      call zero(ujmo)
      call zero(flux_i)
      call zero(flux_j)
      res(:decomp%ii, :decomp%jj, :) = 0.0D0
      res_beta(:decomp%ii, :decomp%jj, :) = 0.0D0
#endif
    case default
      write(6, "(a,i0)") "For start option: ", start
      labort("Unrecognised start option")
  end select
#ifdef BACKGROUND_HSHEAR
  DOLEAPFROG
  call add(psi BLEAPFROG, psi_b)
  call set(rel_p BLEAPFROG, rel BLEAPFROG)
  call add(rel BLEAPFROG, rel_b)
  call add(zeta BLEAPFROG, zeta_b)
  ENDDOLEAPFROG
#endif
  ! Initialise the old value of the relative potential vorticity
  DOLEAPFROG
  call set(zeta_old(:, :, : TLEAPFROG), zeta BLEAPFROG)
  ENDDOLEAPFROG

  !!!!!! Timestep initialisation

#ifdef CABARET
  model_time = dt / 2.0D0
  tweight_a = 0.5D0
  tweight_u = 0.5D0
#else
  model_time = 0.0D0
  T = 2.0D0 * dt
#ifdef FORWARD_EULER_START
  select case(start)
    case(START_PICKUP)
    case(START_PERTURB, START_ZERO)
      call zeta_forward_euler(dt, strat, p_beta_y, &
#ifdef WIND_FORCING
        & force_wind, &
#endif
#ifdef BACKGROUND_VSHEAR
        & u_background_v, &
#endif
        & viscosity, visc_bot, psi(1), rel(1), rel_p(1), zeta(1), zeta_np1 = zeta(2))
      call update(zeta(2))
      call flush()
#ifdef BACKGROUND_HSHEAR
      call subtract(zeta(2), rel_b)
#endif
      call pv_inversion(strat, vbc, zeta(2), hcache, psi(2))
      call update(psi(2))
      call flush()
      call rel_from_zeta_and_psi(strat, vbc, zeta(2), psi(2), rel(2))
      call update_zeta_bc(strat, rel(2), psi(2), zeta(2))
#ifdef BACKGROUND_HSHEAR
      call add(psi(2), psi_b)
      call set(rel_p(2), rel(2))
      call add(rel(2), rel_b)
      call add(zeta(2), zeta_b)
#endif
      call set(zeta_old(:, :, :, 2), zeta(2))
      model_time = model_time + dt
    case default
      write(6, "(a,i0)") "For start option: ", start
      labort("Unrecognised start option")
  end select
#endif
#endif

  !!!!!! Diagnostic initialisation

  diagnostic_time = 0.0D0
  psi_av(:decomp%ii, :decomp%jj, :) = 0.0D0

  if(quad) then
    allocate(u(decomp%ii, decomp%jj, decomp%layers))
    allocate(v(decomp%ii, decomp%jj, decomp%layers))
    allocate(u_i(decomp%ii, decomp%jj, decomp%layers - 1))
    allocate(v_i(decomp%ii, decomp%jj, decomp%layers - 1))
    allocate(tness(decomp%ii, decomp%jj, decomp%layers))
    u = 0.0D0
    v = 0.0D0
    u_i = 0.0D0
    v_i = 0.0D0
    tness = 0.0D0
    call allocate(uu_av, decomp)
    call allocate(vv_av, decomp)
    call allocate(uu_i_av, decomp, ec_type = EC_INTERFACES)
    call allocate(vv_i_av, decomp, ec_type = EC_INTERFACES)
    call allocate(uv_av, decomp)
    call allocate(ur_av, decomp)
    call allocate(vr_av, decomp)
    call allocate(rr_av, decomp)
    call allocate(pr_av, decomp)
    call allocate(uz_av, decomp)
    call allocate(vz_av, decomp)
    call allocate(zz_av, decomp)
    call allocate(zr_av, decomp)
    call allocate(pp_av, decomp)
    call allocate(pz_av, decomp)
    call allocate(ut_av, decomp)
    call allocate(vt_av, decomp)
    call allocate(rt_av, decomp)
    call allocate(tt_av, decomp)
    call allocate(sc_ub_av, decomp, ec_type = EC_INTERFACES)
    call allocate(sc_vb_av, decomp, ec_type = EC_INTERFACES)
    call allocate(sc_bb_av, decomp)
    call allocate(sc_bb_i_av, decomp, ec_type = EC_INTERFACES)
    call zero(uu_av)
    call zero(vv_av)
    call zero(uu_i_av)
    call zero(vv_i_av)
    call zero(uv_av)
    call zero(ur_av)
    call zero(vr_av)
    call zero(rr_av)
    call zero(pr_av)
    call zero(uz_av)
    call zero(vz_av)
    call zero(zz_av)
    call zero(zr_av)
    call zero(pp_av)
    call zero(pz_av)
    call zero(sc_ub_av)
    call zero(sc_vb_av)
    call zero(sc_bb_av)
    call zero(sc_bb_i_av)
  end if

  !!!!!! Output initialisation

  dump_count = 0
  s_index = -1
  call write_s_av_header(decomp, l_scale, t_scale)
  if(quad) call write_s_quad_av_header(decomp, l_scale, t_scale)

  ! ###################################
  ! ############ Time-loop ############
  ! ###################################

  dtime = 0.0D0
  htime = 0.0D0
  itime = 0.0D0
  ! Zero communication times, as only time spent in the time loop should be
  ! measured (i.e. startup times are not measured)
  transpose_wait_time = 0.0D0
#ifdef PH1D
  h1d_time = 0.0D0
#endif
  halo_wait_time = 0.0D0
  wtime_0 = mpi_wtime()
  write(6, "(a)") "Beginning time stepping"
  write(6, "(a,"//real_8_chr//")") "Initial model time (days)  = ", start_time_day
  write(6, "(a,"//real_8_chr//")") "Initial model time (years) = ", start_time_day / 365.0D0
  write(6, "(a,"//real_8_chr//")") "Model run time (days)      = ", (max_spin + max_time) * t_scale / day_s
  write(6, "(a,"//real_8_chr//")") "Model run time (years)     = ", (max_spin + max_time) * t_scale / (day_s * 365.0D0)
  flush(0)
  flush(6)
  exit_timestep_loop = .false.
  timestep_loop: do
  DOLEAPFROG

#ifdef CABARET

    !!!!!! CABARET adaptive timestep

    if(adapt_dt) then
      call cabaret_adapt_dt(target_cfl, max_dt, max_dt_increase, ui, uj, dt, tweight_a, tweight_u)  
    end if

    !!!!!! CABARET predictor step

    lwtime_0 = mpi_wtime()
    call cabaret_predictor(dt, zeta_old, res, zeta)

    !!!!!! Beta term

    call cabaret_beta(dt, tweight_a, beta, ui, zeta, res_beta)
    call update(zeta)
    call flush()
    itime = itime + mpi_wtime() - lwtime_0
#else

    !!!!!! Centred in space and time integration

    lwtime_0 = mpi_wtime()
    call zeta_centred(leapfrog, T, strat, p_beta_y, &
#ifdef WIND_FORCING
      & force_wind, &
#endif
#ifdef BACKGROUND_VSHEAR
      & u_background_v, &
#endif
      & viscosity, visc_bot, psi, rel, rel_p, zeta)
    call update(zeta BLEAPFROG)
    call flush()

    !!!!!! Robert-Asselin-Williams time filter

    call robert_asselin_filter(leapfrog, zeta, zeta_old)
    itime = itime + mpi_wtime() - lwtime_0
#endif

    !!!!!! PV inversion (elliptic solves)

    lwtime_0 = mpi_wtime()
#ifdef BACKGROUND_HSHEAR
    call subtract(zeta BLEAPFROG, zeta_b)
#endif
    call pv_inversion(strat, vbc, zeta BLEAPFROG, hcache, psi BLEAPFROG)
    call update(psi BLEAPFROG)
    call flush()

    call rel_from_zeta_and_psi(strat, vbc, zeta BLEAPFROG, psi BLEAPFROG, rel BLEAPFROG)
    call update_zeta_bc(strat, rel BLEAPFROG, psi BLEAPFROG, zeta BLEAPFROG)
#ifdef BACKGROUND_HSHEAR
    call add(psi BLEAPFROG, psi_b)
    call set(rel_p BLEAPFROG, rel BLEAPFROG)
    call add(rel BLEAPFROG, rel_b)
    call add(zeta BLEAPFROG, zeta_b)
#endif
    htime = htime + mpi_wtime() - lwtime_0

#ifdef CABARET
    ! With CABARET zeta has moved one half timestep forwards
    model_time = model_time + 0.5D0 * dt
#else
    ! With the centred scheme the model has moved one timestep forwards
    model_time = model_time + dt
#endif
    model_time_day = model_time * t_scale / day_s
    day_number = int(start_time_day + model_time_day + 1.0D-3)

    !!!!!! Diagnostics
    
#ifdef BACKGROUND_HSHEAR
    call subtract(psi BLEAPFROG, psi_b)
    call subtract(rel BLEAPFROG, rel_b)
    call subtract(zeta BLEAPFROG, zeta_b)
#endif

    if(start_time + model_time > max_spin) then
      lwtime_0 = mpi_wtime()

      do layer = 1, layers
        do j = 1, decomp%jj
          do i = 1, decomp%ii
            psi_av(i, j, layer) = psi_av(i, j, layer) + dt * psi BLEAPFROG %sdata(i, j, layer)
          end do
        end do
      end do

      if(quad) then
        call zonal_velocity(vbc, psi BLEAPFROG, rel BLEAPFROG, u)
        call meridional_velocity(vbc, psi BLEAPFROG, rel BLEAPFROG, v)
        call interface_velocity(decomp, u, v, u_i, v_i)
        call thickness(psi BLEAPFROG, strat, tness)

        ! Velocity - velocity
        call add(uu_av, u, dt = dt)
        call add(vv_av, v, dt = dt)
        call add(uu_i_av, u_i, dt = dt)
        call add(vv_i_av, v_i, dt = dt)
        call add(uv_av, u, v, dt = dt)
        ! Velocity - vorticity
        call add(ur_av, u, rel BLEAPFROG, dt = dt)
        call add(vr_av, v, rel BLEAPFROG, dt = dt)
        ! Velocity - relative potential vorticity
        call add(uz_av, u, zeta BLEAPFROG, dt = dt)
        call add(vz_av, v, zeta BLEAPFROG, dt = dt)
        ! Vorticity - vorticity
        call add(rr_av, rel BLEAPFROG, dt = dt)
        ! Stream function - vorticity
        call add(pr_av, psi BLEAPFROG, rel BLEAPFROG, dt = dt)
        ! Relative potential vorticity - relative potential vorticity
        call add(zz_av, zeta BLEAPFROG, dt = dt)
        ! Relative potential vorticity - vorticity
        call add(zr_av, zeta BLEAPFROG, rel BLEAPFROG, dt = dt)
        ! Stream function - stream function
        call add(pp_av, psi BLEAPFROG, dt = dt)
        ! Stream function - relative potential vorticity
        call add(pz_av, psi BLEAPFROG, zeta BLEAPFROG, dt = dt)
        ! Velocity - thickness
        call add(ut_av, u, tness, dt = dt)
        ! Velocity - thickness
        call add(vt_av, v, tness, dt = dt)
        ! Vorticity - thickness
        call add(rt_av, tness, rel BLEAPFROG, dt = dt)
        ! Thickness - thickness
        call add(tt_av, tness, dt = dt)
        ! Scaled buoyancy covariances
        call add_scaled_buoyancy_corr(sc_ub_av, sc_vb_av, sc_bb_av, sc_bb_i_av, strat, u, v, psi BLEAPFROG, dt)
      end if

      diagnostic_time = diagnostic_time + dt

      dtime = dtime + mpi_wtime() - lwtime_0
    end if

    !!!!!! Diagnostic output

    if( &
#ifndef CABARET
      & leapfrog == 2 .and. &
#endif
      & mod(day_number, k_out) == 0 .and. last_diagnostics_day < day_number) then
      last_diagnostics_day = day_number

      wtime = mpi_wtime() - wtime_0
      lwtime_0 = wtime - itime - htime - dtime
#ifdef CABARET
      cfl = cabaret_cfl(dt, ui, uj)
#endif
      call energy(strat, psi BNLEAPFROG, ekin, epot)

      if(decomp%p == 0) then
        write(6, "(a)") "Writing diagnostics"
        write(6, "(a,"//real_8_chr//")") "Time (days)       = ", start_time_day + model_time_day
        write(6, "(a,"//real_8_chr//")") "Time (years)      = ", (start_time_day + model_time_day) / 365.0D0
#ifdef CABARET
        write(6, "(a,"//real_8_chr//")") "dt (hours)        = ", dt * t_scale / 3600.0D0
        write(6, "(a,"//real_8_chr//")") "Max. Courant no.  = ", cfl
#endif

        write(6, "(a,"//real_8_chr//")") "Total             = ", wtime
        write(6, "(a,"//real_8_chr//",a,i0,a)") "Integration       = ", itime, " (", int((itime / wtime) * 100.0D0), "%)"
        write(6, "(a,"//real_8_chr//",a,i0,a)") "Helmholtz         = ", htime, " (", int((htime / wtime) * 100.0D0), "%)"
        write(6, "(a,"//real_8_chr//",a,i0,a)") "Diagnostics       = ", dtime, " (", int((dtime / wtime) * 100.0D0), "%)"
        write(6, "(a,"//real_8_chr//",a,i0,a)") "Other             = ", lwtime_0, " (", int((lwtime_0 / wtime) * 100.0D0), "%)"
        if(decomp%n > 1) then
#ifndef PETSC
          write(6, "(a,"//real_8_chr//",a,i0,a)") "Transpose wait    = ", transpose_wait_time, &
            & " (", int((transpose_wait_time / wtime) * 100.0D0), "%)"
#ifdef PH1D
          write(6, "(a,"//real_8_chr//",a,i0,a)") "1D Helmholtz comm = ", h1d_time, &
            & " (", int((h1d_time / wtime) * 100.0D0), "%)"
#endif
#endif
          write(6, "(a,"//real_8_chr//",a,i0,a)") "Halo wait         = ", halo_wait_time, &
            & " (", int((halo_wait_time / wtime) * 100.0D0), "%)"
        end if
        call print(ekin, epot)
      end if

    end if

#ifdef CABARET
    lwtime_0 = mpi_wtime()

    !!!!!! Viscosity, wind forcing, and bottom friction

    call cabaret_forcing(dt, &
#ifdef WIND_FORCING
      & force_wind, &
#endif
      & viscosity, visc_bot, rel_p, zeta)
    call update(zeta)
    call flush()

    !!!!!! Flux variables of the CABARET scheme

    call cabaret_cell_boundaries(tweight_u, psi, ui, uj, uimo, ujmo &
#ifdef BACKGROUND_VSHEAR
      & , u_background_v &
#endif
      & )
    call update(ui, level_1 = .true.)
    call update(uj, level_1 = .true.)
    call flush()

    !!!!!! CABARET extrapolator step

    call cabaret_extrapolator(dt, zeta, zeta_old, ui, uj, zeta_fi, zeta_fj)
    call update(zeta_fi, level_1 = .true.)
    call update(zeta_fj, level_1 = .true.)
    call flush()
    call cabaret_flux(ui, uj, zeta_fi, zeta_fj, flux_i, flux_j)
    itime = itime + mpi_wtime() - lwtime_0
#endif

    ! #################################################
    ! ############ Time loop exit criteria ############
    ! #################################################

    if( &
#ifndef CABARET
      & leapfrog == 2 .and. &
#endif
      & start_time + model_time >= max_spin + max_time) exit_timestep_loop = .true.

    !!!!!! Data output

    if(exit_timestep_loop .or. ( &
#ifndef CABARET
      & leapfrog == 2 .and. &
#endif
      & (mod(day_number, k_save) == 0 .and. last_dump_day < day_number))) then
      write(6, "(a)") "Writing output"
      last_dump_day = day_number
      dump_count = dump_count + 1

      if(mod(dump_count, k_dump) == 1 .or. k_dump == 1) then
        s_index = s_index + 1
        call write_s_header(decomp, s_index, l_scale, t_scale)
      end if
      call open_parallel_data(decomp, unit, "qg" // int2str(layers) // "_" // int2str(decomp%g_ii) // "_s", &
        & output_index = s_index)
      call write_unformatted(unit, model_time_day)
      call write_unformatted(unit, start_time_day + model_time_day)
      do layer = 1, layers
#ifdef CABARET
        call write_unformatted(unit, psi%sdata(:, :, layer))
#else
        call write_unformatted(unit, psi(1)%sdata(:, :, layer), psi(2)%sdata(:, :, layer))
#endif
      end do
      do layer = 1, layers
#ifdef CABARET
        call write_unformatted(unit, zeta%sdata(:, :, layer))
#else
        call write_unformatted(unit, zeta(1)%sdata(:, :, layer), zeta(2)%sdata(:, :, layer))
#endif
      end do
      close(unit)

#ifdef CABARET
      if(exit_timestep_loop .or. mod(dump_count, k_cp) == 0) then
        call write_s_supp_header(decomp, s_index, l_scale, t_scale)
        call open_parallel_data(decomp, unit, "qg" // int2str(layers) // "_" // int2str(decomp%g_ii) // "_s_supp", &
          & output_index = s_index)
        call write_unformatted(unit, model_time_day)
        call write_unformatted(unit, start_time_day + model_time_day)
        call write_unformatted(unit, dt)
        do layer = 1, layers
          call write_unformatted(unit, ui%data(decomp%fa_i_min + 1:decomp%fa_i_max, decomp%fa_j_min + 1:decomp%fa_j_max - 1, layer))
        end do
        do layer = 1, layers
          call write_unformatted(unit, uj%data(decomp%fa_i_min + 1:decomp%fa_i_max - 1, decomp%fa_j_min + 1:decomp%fa_j_max, layer))
        end do    
        do layer = 1, layers
          call write_unformatted(unit, zeta_fi%data(decomp%fa_i_min + 1:decomp%fa_i_max, decomp%fa_j_min + 1:decomp%fa_j_max - 1, layer))
        end do
        do layer = 1, layers
          call write_unformatted(unit, zeta_fj%data(decomp%fa_i_min + 1:decomp%fa_i_max - 1, decomp%fa_j_min + 1:decomp%fa_j_max, layer))
        end do      
        do layer = 1, layers
          call write_unformatted(unit, res(:decomp%ii, :decomp%jj, layer))
        end do       
        do layer = 1, layers
          call write_unformatted(unit, res_beta(:decomp%ii, :decomp%jj, layer))
        end do
        close(unit)
      end if
#endif

      if(exit_timestep_loop .or. (start_time + model_time > max_spin &
#ifdef CABARET
        & .and. mod(dump_count, k_cp) == 0 &
#endif
        & )) then
        cff = 1.0D0 / diagnostic_time

        call open_parallel_data(decomp, unit, "qg" // int2str(layers) // "_" // int2str(decomp%g_ii) // "_s_av")
        call write_unformatted(unit, model_time_day)
        call write_unformatted(unit, start_time_day + model_time_day)
        do layer = 1, layers
          call write_unformatted(unit, psi_av(:decomp%ii, :decomp%jj, layer) * cff)
        end do
        close(unit)

        if(quad) then
          call open_parallel_data(decomp, unit, "qg" // int2str(layers) // "_" // int2str(decomp%g_ii) // "_quad_av")
          call write_unformatted(unit, model_time_day)
          call write_unformatted(unit, start_time_day + model_time_day)
          call write_unformatted(unit, uu_av%data * cff)
          call write_unformatted(unit, vv_av%data * cff)
          call write_unformatted(unit, uu_i_av%data * cff)
          call write_unformatted(unit, vv_i_av%data * cff)
          call write_unformatted(unit, uv_av%data * cff)
          call write_unformatted(unit, ur_av%data * cff)
          call write_unformatted(unit, vr_av%data * cff)
          call write_unformatted(unit, rr_av%data * cff)
          call write_unformatted(unit, pr_av%data * cff)
          call write_unformatted(unit, uz_av%data * cff)
          call write_unformatted(unit, vz_av%data * cff)
          call write_unformatted(unit, zz_av%data * cff)
          call write_unformatted(unit, zr_av%data * cff)
          call write_unformatted(unit, pp_av%data * cff)
          call write_unformatted(unit, pz_av%data * cff)
          call write_unformatted(unit, ut_av%data * cff)
          call write_unformatted(unit, vt_av%data * cff)
          call write_unformatted(unit, rt_av%data * cff)
          call write_unformatted(unit, tt_av%data * cff)
          call write_unformatted(unit, sc_ub_av%data * cff)
          call write_unformatted(unit, sc_vb_av%data * cff)
          call write_unformatted(unit, sc_bb_av%data * cff)
          call write_unformatted(unit, sc_bb_i_av%data * cff)
          close(unit)
        end if
      end if

      flush(0)
      flush(6)
    end if
    
#ifdef BACKGROUND_HSHEAR
    call add(psi BLEAPFROG, psi_b)
    call add(rel BLEAPFROG, rel_b)
    call add(zeta BLEAPFROG, zeta_b)
#endif

    ! ##############################################
    ! ############ Time loop exit point ############
    ! ##############################################

    if(exit_timestep_loop) exit timestep_loop

#ifdef CABARET
    !!!!!! CABARET corrector step

    lwtime_0 = mpi_wtime()
    call cabaret_corrector(dt, flux_i, flux_j, res, zeta)
    call update(zeta)
    call flush()

    ! Store old relative potential vorticity values for use on the next timestep
    call set(zeta_old, zeta)

    ! zeta has moved one half timestep forwards
    model_time = model_time + 0.5D0 * dt
    !model_time_day = time * t_scale / day_s
    !day_number = int(start_time_day + model_time_day + 1.0D-3)
    itime = itime + mpi_wtime() - lwtime_0
#endif

  ENDDOLEAPFROG
  end do timestep_loop
  write(6, "(a)") "Completed time stepping"
  write(6, "(a,"//real_8_chr//")") "Final model time (days)  = ", start_time_day + model_time_day
  write(6, "(a,"//real_8_chr//")") "Final model time (years) = ", (start_time_day + model_time_day) / 365.0D0

  ! ##########################################
  ! ############ End of time-loop ############
  ! ##########################################

  !!!!!! Finalisation

  ! Ensure all communications have completed prior to deallocation
  call flush()
  if(quad) then
    deallocate(u)
    deallocate(v)
    deallocate(u_i)
    deallocate(v_i)
    deallocate(tness)
    call deallocate(uu_av)
    call deallocate(vv_av)
    call deallocate(uu_i_av)
    call deallocate(vv_i_av)
    call deallocate(uv_av)
    call deallocate(ur_av)
    call deallocate(vr_av)
    call deallocate(rr_av)
    call deallocate(pr_av)
    call deallocate(uz_av)
    call deallocate(vz_av)
    call deallocate(zz_av)
    call deallocate(zr_av)
    call deallocate(pp_av)
    call deallocate(pz_av)
    call deallocate(ut_av)
    call deallocate(vt_av)
    call deallocate(rt_av)
    call deallocate(tt_av)
    call deallocate(sc_ub_av)
    call deallocate(sc_vb_av)
    call deallocate(sc_bb_av)
    call deallocate(sc_bb_i_av)
  end if
  call deallocate(hcache)
  call deallocate(strat)
  call deallocate(psi)
  call deallocate(rel)
  call deallocate(zeta)
#ifdef CABARET
  call deallocate(zeta_fi)
  call deallocate(zeta_fj)
  call deallocate(ui)
  call deallocate(uj)
  call deallocate(uimo)
  call deallocate(ujmo)
  call deallocate(flux_i)
  call deallocate(flux_j)
#endif
#ifdef BACKGROUND_HSHEAR
  call deallocate(psi_b)
  call deallocate(rel_b)
  call deallocate(rel_p)
#endif
  call deallocate(decomp)
  call parallel_finalise()      

end program qg3
