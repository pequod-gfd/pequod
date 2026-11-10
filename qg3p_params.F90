#include "config.h"

! Parameters for the 3D QG model
module qg_parameters

  implicit none
      
  integer, parameter :: FILENAME_LEN = 1024
  integer, parameter :: complex_8 = kind((1.0D0, 1.0D0))
  integer, parameter :: real_8 = kind(1.0D0)
  integer, parameter :: real_8_bit_size = 8
  character(len = *), parameter :: real_8_chr = "e26.17e3"
  
  ! Grid
  integer, parameter :: layers = 3
  integer, parameter, dimension(2) :: iis = (/128, 129/)
  integer, parameter, dimension(2) :: jjs = (/128, 129/)
  ! The maximum grid size over all processes
  integer, parameter :: ii_m = 129, jj_m = 129
  
  ! Periodic (in longitude) channel geometry?
  logical, parameter :: periodic = .false.
  
  ! Time stepping  
  real(kind = real_8), parameter :: dt = 1800.0D0
#ifdef CABARET
  logical, parameter :: adapt_dt = .true.
  real(kind = real_8), parameter :: max_dt = 86400.0D0, max_dt_increase = 1.5D0
  real(kind = real_8), parameter :: target_cfl = 0.4D0
#else
  ! Leapfrog filter weights
  real(kind = real_8), parameter :: ra_nu = 0.02D0
#ifdef RA_DOUBLE_FILTER
  ! Williams 2009 (MWR) modified filter parameter
  real(kind = real_8), parameter :: ra_alpha = 0.5D0
#endif
#endif
  ! How to set the simulation start time
  integer, parameter :: START_TIME_AUTO = 1, START_TIME_MANUAL = 2
  integer, parameter :: set_start_time = START_TIME_AUTO
  ! Simulation start time in days when using START_TIME_MANUAL
  real(kind = real_8) :: start_time_override = huge(0.0D0)
  ! Time in days for spinup (no output) and run time
  real(kind = real_8), parameter :: max_spin = 0.0D0, max_time = 1.0D0 * 10.0D0
  
  ! Physical parameters
  real(kind = real_8), parameter :: basinscale = 3840.0D5
  real(kind = real_8), parameter :: beta = 2.0D-13
  ! How to set the stratification
  integer, parameter :: STRAT_RD = 1, STRAT_INPUT = 2
  integer, parameter :: set_strat = STRAT_RD
  ! Stratification parameters set the stratification with STRAT_RD
  real(kind = real_8), dimension(layers), save :: H = (/0.25D5, 0.75D5, 3.0D5/)
  real(kind = real_8), parameter :: Rd_bcl1 = 40.0D5, Rd_bcl2 = 23.0D5
#ifdef WIND_FORCING
  ! Wind forcing
  real(kind = real_8), parameter :: tau = 0.8D0, wind_shift = 0.0D0, wind_tilt = 0.2D0, wind_asym = 0.9D0
#endif
#ifdef BACKGROUND_VSHEAR
  ! How to set the shear
  integer, parameter :: SHEAR_U = 1, SHEAR_INPUT = 2
  integer, parameter :: set_shear = SHEAR_U
  real(kind = real_8), dimension(layers) :: U = (/6.0D0, 3.0D0, 0.0D0/)
#endif
  ! Viscosity and bottom drag
  real(kind = real_8), parameter :: visc = 100.0D4, visc_bot = 4.0D-8
  
  ! Boundary conditions
  ! Use partial slip? Otherwise, use no slip.
  logical, parameter :: partial_slip = .true.
  ! Partial slip parameter (length scale)
  real(kind = real_8), parameter :: alpha = 120.0D5
  
  ! Input / Output
  ! Time in days for energy diagnostic output
  integer, parameter :: k_out = 1
  ! Time in days for data output. The data is saved every k_save days, a new
  ! dump file is written to every k_dump * k_save days, and a full checkpoint
  ! is written every k_cp * k_save days. The latter applies only when using
  ! CABARET.
  integer, parameter :: k_save = 1, k_dump = 1, k_cp = 10
  ! Enable quadratic term diagnostics
  logical, parameter :: quad = .false.
  ! Input / output unit
  integer, parameter :: unit = 11
  ! How to start the simulation
  integer, parameter :: START_ZERO = 1, START_PERTURB = 2, START_PICKUP = 3
  integer, parameter :: start = START_ZERO
  ! Perturbation applied to the relative potential vorticity when starting with
  ! START_PERTURB
  real(kind = real_8) :: zeta_perturb = 1.0D-7
  ! Log standard output and standard error?
  logical, parameter :: save_log = .false.
  
  ! Non-dimensional velocity scale
  real(kind = real_8), parameter :: u_scale = 1.0D0
  
  ! Constants  
  real(kind = real_8), parameter :: INFINITY = huge(0.0D0) * epsilon(0.0D0)
  real(kind = real_8), parameter :: pi = 3.14159265358979323846D0
  real(kind = real_8), parameter :: two_pi = 6.28318530717959D0
  real(kind = real_8), parameter :: sqrt3 = 1.73205080756888D0
  real(kind = real_8), parameter :: hsqt2 = 0.7071067811865475D0
  real(kind = real_8), parameter :: day_s = 86400.0D0

#ifdef BACKGROUND_HSHEAR
contains

  ! Background velocity profile.
  function u_background_h(y) result(u)
    real(kind = real_8), intent(in) :: y
    
    real(kind = real_8) :: u

    u = 0.0D0
    
  end function u_background_h
#endif

end module qg_parameters
