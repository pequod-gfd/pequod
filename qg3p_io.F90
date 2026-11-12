#include "config.h"
#include "leapfrog.h"

! QG model I/O routines.
module qg_io

  use qg_diagnostics
  use qg_parallel
  use qg_parameters, only : FILENAME_LEN, real_8, real_8_bit_size, real_8_chr
  use test_tools

  implicit none
  
  include "mpif.h"
  
  private
  
  public :: read_stratification, write_stratification, read_shear_v, &
    & write_shear_v, read_parallel_checkpoint, write_unformatted
  public :: print_config, log_output, open_parallel_header, &
    & open_parallel_data, write_parallel_header, write_s_header, &
    & write_s_av_header, write_s_quad_av_header
#ifdef CABARET
  public :: write_s_supp_header
#endif
  
  interface read_unformatted
    module procedure read_unformatted_rank0, read_unformatted_rank1, &
      & read_unformatted_rank2, read_unformatted_rank2_rank2, &
      & read_unformatted_rank3
  end interface read_unformatted
  
  interface write_unformatted
    module procedure write_unformatted_rank0, write_unformatted_rank1, &
      & write_unformatted_rank2, write_unformatted_rank2_rank2, &
      & write_unformatted_rank3
  end interface write_unformatted
  
contains

  ! Check that a given file exists, and exit with an error otherwise.
  subroutine assert_file_exists(filename)
    character(len = *), intent(in) :: filename
    
    logical :: file_exists
    
    inquire(file = trim(filename), exist = file_exists)
    if(.not. file_exists) then
      labort("File " // trim(filename) // " not found")
    end if

  end subroutine assert_file_exists

  ! Read stratification parameters from a stratification parameter file. This
  ! allocates the stratification.
  subroutine read_stratification(strat, l_scale)
    use qg_parameters, only : unit
    
    type(stratification), intent(out) :: strat
    real(kind = real_8), optional, intent(in) :: l_scale
    
    character(len = *), parameter :: filename = "input/stratification"
    integer :: ierr, layer, layers, stat
    real(kind = real_8), dimension(:), allocatable :: h, s1, s2
    
    if(procno() == 1) then    
      call assert_file_exists(filename)      
      open(unit, file = trim(filename), iostat = stat)
      if(stat /= 0) then
        labort("Failed to open input file")
      end if
      
      read(unit, *) layers
      if(layers <= 0) then
        write(6, "(a,i0)") "Number of layers: ", layers
        labort("Invalid number of layers")
      end if
      call mpi_bcast(layers, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
      allocate(h(layers), s1(layers - 1), s2(layers - 1))
      if(layers == 1) then
        read(unit, *) h(1)
      else
        read(unit, *) h(1), s2(1)
        do layer = 2, layers - 1
          read(unit, *) h(layer), s1(layer - 1), s2(layer)
        end do
        read(unit, *) h(layers), s1(layers - 1)
      end if
      close(unit)      
    else
      call mpi_bcast(layers, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
      allocate(h(layers), s1(layers - 1), s2(layers - 1))
    end if    
    
    call mpi_bcast(h(1), size(h), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
    call mpi_bcast(s1(1), size(s1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
    call mpi_bcast(s2(1), size(s2), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
    
    call allocate(strat, h, s1, s2, l_scale = l_scale)
    deallocate(h, s1, s2)
    
  end subroutine read_stratification
  
  ! Read vertical shear parameters from a vertical shear parameter file.
  subroutine read_shear_v(u, u_scale)
    use qg_parameters, only : unit
    
    real(kind = real_8), dimension(:), intent(out) :: u
    real(kind = real_8), intent(in) :: u_scale
    
    character(len = *), parameter :: filename = "input/u"
    integer :: i, ierr, layers, stat
    
    if(procno() == 1) then
      call assert_file_exists(filename)
      open(unit, file = trim(filename), iostat = stat)
      if(stat /= 0) then
        labort("Failed to open input file")
      end if
      
      read(unit, *) layers
      if(layers /= size(u)) then
        write(6, "(a,i0)") "Number of layers: ", layers
        labort("Invalid number of layers")
      end if
      do i = 1, layers
        read(unit, *) u(i)
      end do
      close(unit)
    end if
    
    call mpi_bcast(u(1), size(u), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
    
    u = u / u_scale
      
  end subroutine read_shear_v
  
  ! Write (dimensional) stratification parameters to a stratification parameter
  ! file.
  subroutine write_stratification(strat)
    use qg_parameters, only : unit
  
    character(len = *), parameter :: filename = "stratification"
    integer :: layer, stat
    type(stratification), intent(in) :: strat
    
    if(procno() == 1) then
      open(unit, file = trim(filename), iostat = stat)
      if(stat /= 0) then
        labort("Failed to open output file")
      end if
      
      write(unit, "(i0)") strat%layers
      if(strat%layers == 1) then
        write(unit, "("//real_8_chr//")") strat%h(1) * strat%l_scale
      else
        write(unit, "("//real_8_chr//","//real_8_chr//")") &
          & strat%h(1) * strat%l_scale, strat%s2(1) / (strat%l_scale * strat%l_scale)
        do layer = 2, strat%layers - 1
          write(unit, "("//real_8_chr//","//real_8_chr//","//real_8_chr//")") &
            & strat%h(layer) * strat%l_scale, strat%s1(layer) / (strat%l_scale * strat%l_scale), strat%s2(layer) / (strat%l_scale * strat%l_scale)
        end do
        write(unit, "("//real_8_chr//","//real_8_chr//")") &
          & strat%h(strat%layers) * strat%l_scale, strat%s1(strat%layers) / (strat%l_scale * strat%l_scale)
      end if
      
      close(unit)
    end if
    
  end subroutine write_stratification
  
  ! Write (dimensional) vertical shear parameters to a vertical shear parameter
  ! file.
  subroutine write_shear_v(u, u_scale)
    use qg_parameters, only : unit
  
    real(kind = real_8), dimension(:), intent(in) :: u
    real(kind = real_8), intent(in) :: u_scale
    
    character(len = *), parameter :: filename = "u"
    integer :: layer, stat
    
    if(procno() == 1) then
      open(unit, file = trim(filename), iostat = stat)
      if(stat /= 0) then
        labort("Failed to open output file")
      end if
      
      write(unit, "(i0)") size(u)
      do layer = 1, size(u)
        write(unit, "("//real_8_chr//")") u(layer) * u_scale
      end do
      
      close(unit)
    end if
    
  end subroutine write_shear_v
  
  ! Read non-recoverable data from parallel data files, in order to pickup a
  ! simulation. This does *not* update the halo regions (only level 1 flux halos
  ! are updated for flux variables).
  subroutine read_parallel_checkpoint(total_time_day, psi, zeta &
#ifdef CABARET
    & , dt, ui, uj, zeta_fi, zeta_fj, res, res_beta &
#endif
    & )
    use qg_parameters, only : unit
  
    real(kind = real_8), intent(out) :: total_time_day
    type(parallel_data) DIMSLEAPFROG, intent(inout) :: psi
    type(parallel_data) DIMSLEAPFROG, intent(inout) :: zeta
#ifdef CABARET
    real(kind = real_8), intent(out) :: dt
    type(parallel_flux_data), intent(inout) :: ui
    type(parallel_flux_data), intent(inout) :: uj
    type(parallel_flux_data), intent(inout) :: zeta_fi
    type(parallel_flux_data), intent(inout) :: zeta_fj
    real(kind = real_8), dimension(:, :, :), intent(out) :: res
    real(kind = real_8), dimension(:, :, :), intent(out) :: res_beta
#endif
    
    character(len = FILENAME_LEN) :: filename
    integer :: layer, layers, stat, marker_size, scalar_size, field_size
    integer(kind = selected_int_kind(18)) :: file_size, snapshot_size
    real(kind = real_8) :: model_time_day
    type(decomposition), pointer :: decomp
#ifdef CABARET
    real(kind = real_8) :: model_time_day_supp, total_time_day_supp
#endif
    
#ifdef CABARET
    assert(associated(zeta%decomp, psi%decomp))
    assert(associated(zeta_fi%decomp, psi%decomp))
    assert(zeta_fi%flux_dir == FLUX_II)
    assert(associated(zeta_fj%decomp, psi%decomp))
    assert(zeta_fj%flux_dir == FLUX_JJ)
    assert(associated(ui%decomp, psi%decomp))
    assert(ui%flux_dir == FLUX_II)
    assert(associated(uj%decomp, psi%decomp))
    assert(uj%flux_dir == FLUX_JJ)
    assert(size(res, 1) >= psi%decomp%ii)
    assert(size(res, 2) >= psi%decomp%jj)
    assert(size(res, 3) == psi%decomp%layers)
    assert(size(res_beta, 1) >= psi%decomp%ii)
    assert(size(res_beta, 2) >= psi%decomp%jj)
    assert(size(res_beta, 3) == psi%decomp%layers)
#else
    assert(associated(psi(1)%decomp, psi(2)%decomp))
    assert(associated(zeta(1)%decomp, psi(2)%decomp))
    assert(associated(zeta(2)%decomp, psi(2)%decomp))
#endif
    
    decomp => psi BNLEAPFROG %decomp
    layers = decomp%layers
    
    filename = trim(parallel_filename("input/qg" // int2str(layers) // "_" // int2str(decomp%g_ii) // "_s.d"))
    call assert_file_exists(filename)    
    open(unit, file = trim(filename), access = "stream", form = "unformatted", iostat = stat)
    if(stat /= 0) then
      labort("Failed to open input file")
    end if
    
    ! Each dump can contain several snapshots. Pick up from the last one.
    inquire(iolength = marker_size) 0, 0
    inquire(iolength = scalar_size) model_time_day
#ifdef CABARET
    inquire(iolength = field_size) psi%sdata(:, :, 1)
#else
    inquire(iolength = field_size) psi(1)%sdata(:, :, 1), psi(2)%sdata(:, :, 1)
#endif
    snapshot_size = 2 * (scalar_size + marker_size) + 2 * int(layers, kind(snapshot_size)) * (field_size + marker_size)
    inquire(unit, size = file_size)
    if(file_size < snapshot_size .or. mod(file_size, snapshot_size) /= 0) then
      labort("Invalid snapshot file size")
    end if
    read(unit, pos = file_size - snapshot_size + 1)
    call read_unformatted(unit, model_time_day)
    call read_unformatted(unit, total_time_day)
#ifdef CABARET
    do layer = 1, layers
      call read_unformatted(unit, psi%sdata(:, :, layer))
    end do
    do layer = 1, layers
      call read_unformatted(unit, zeta%sdata(:, :, layer))
    end do
#else
    do layer = 1, layers
      call read_unformatted(unit, psi(1)%sdata(:, :, layer), psi(2)%sdata(:, :, layer))
    end do
    do layer = 1, layers
      call read_unformatted(unit, zeta(1)%sdata(:, :, layer), zeta(2)%sdata(:, :, layer))
    end do
#endif
    close(unit)
    
#ifdef CABARET
    filename = parallel_filename("input/qg" // int2str(layers) // "_" // int2str(decomp%g_ii) // "_s_supp.d")
    call assert_file_exists(filename)    
    open(unit, file = trim(filename), access = "stream", form = "unformatted", iostat = stat)
    if(stat /= 0) then
      labort("Failed to open input file")
    end if
        
    call read_unformatted(unit, model_time_day_supp)
    call read_unformatted(unit, total_time_day_supp)
    if(fne(model_time_day_supp, model_time_day) .or. &
     & fne(total_time_day_supp, total_time_day)) then
      labort("Input times do not match in CABARET checkpoint data")
    end if
    call read_unformatted(unit, dt)
    do layer = 1, layers
      call read_unformatted(unit, ui%data(decomp%fa_i_min + 1:decomp%fa_i_max, decomp%fa_j_min + 1:decomp%fa_j_max - 1, layer))
    end do
    do layer = 1, layers
      call read_unformatted(unit, uj%data(decomp%fa_i_min + 1:decomp%fa_i_max - 1, decomp%fa_j_min + 1:decomp%fa_j_max, layer))
    end do
    do layer = 1, layers
      call read_unformatted(unit, zeta_fi%data(decomp%fa_i_min + 1:decomp%fa_i_max, decomp%fa_j_min + 1:decomp%fa_j_max - 1, layer))
    end do
    do layer = 1, layers
      call read_unformatted(unit, zeta_fj%data(decomp%fa_i_min + 1:decomp%fa_i_max - 1, decomp%fa_j_min + 1:decomp%fa_j_max, layer))
    end do
    do layer = 1, layers
      call read_unformatted(unit, res(:decomp%ii, :decomp%jj, layer))
    end do
    do layer = 1, layers
      call read_unformatted(unit, res_beta(:decomp%ii, :decomp%jj, layer))
    end do
    close(unit)
#endif
    
  end subroutine  read_parallel_checkpoint
  
  ! Print preprocessor defined model configuration options.
  subroutine print_config()
    PRINT_PEQUOD_VERSION
    write(6, "(a)") "Advection scheme: " // &
#ifdef CABARET
#ifdef CABARET_LIMITER
      & "CABARET (with flux limiter)"
#else
      & "CABARET (with *no* flux limiter)"
#endif
#else
      & "Centred leapfrog " // &
#ifdef RA_DOUBLE_FILTER
      & "with Robert-Asselin-Williams filter"
#else
      & "with Robert-Asselin filter"
#endif
#ifdef CENTRED_BOTTOM_FRICTION
    write(6, "(a)") "WARNING: Centred bottom friction"
#endif
#ifndef FORWARD_EULER_START
    write(6, "(a)") "WARNING: No forward Euler timestep at start"
#endif
#endif

    write(6, "(a)") "Linear solver: " // &
#ifdef PETSC
      & "PETSc"
#else
      & "Native"
#endif
  
  end subroutine print_config
    
  ! Redirect standard error and output to log files.
  subroutine log_output(base)
    character(len = *), intent(in) :: base
    
    integer :: ierr
    
#ifdef __INTEL_COMPILER
    labort("log_output does not work with the Intel Fortran compiler")
#else
    close(0, iostat = ierr)
    if(ierr /= 0) then
      labort("Failed to close standard error")
    end if
    open(0, file = trim(parallel_filename(base // ".err")), iostat = ierr)
    if(ierr /= 0) then
      labort("Failed to open standard error")
    end if
    
    close(6, iostat = ierr)
    if(ierr /= 0) then
      labort("Failed to close standard output")
    end if
    open(6, file = trim(parallel_filename(base // ".log")), iostat = ierr)
    if(ierr /= 0) then
      labort("Failed to open standard output")
    end if
#endif
    
  end subroutine log_output
  
  ! Open a parallel data header file.
  subroutine open_parallel_header(decomp, unit, filename, output_index, overwrite)
    type(decomposition), intent(in) :: decomp
    integer, intent(in) :: unit
    character(len = *), intent(in) :: filename
    integer, optional, intent(in) :: output_index
    logical, optional, intent(in) :: overwrite
  
    character(len = filename_len) :: lfilename
    character(7) :: open_status
    integer :: ierr
  
    lfilename = filename
    if(present(output_index)) lfilename = trim(lfilename) // "." // int2str(output_index)
    lfilename = parallel_filename(lfilename, decomp)
    
    open_status = 'replace'
    if (present(overwrite)) then
      if (.not. overwrite) open_status = 'new'
    end if
    open(unit, file = trim(lfilename), form = "formatted", status = open_status, iostat = ierr)
    if(ierr /= 0) then
      labort("Failed to open header file: " // trim(lfilename))
    end if
  
  end subroutine open_parallel_header
  
  ! Open a parallel data file.
  subroutine open_parallel_data(decomp, unit, filename, output_index, append, overwrite)
    type(decomposition), intent(in) :: decomp
    integer, intent(in) :: unit
    character(len = *), intent(in) :: filename
    integer, optional, intent(in) :: output_index
    logical, optional, intent(in) :: append
    logical, optional, intent(in) :: overwrite
    
    character(len = filename_len) :: lfilename
    character(7) :: open_status
    integer :: ierr
    logical :: append_data
    
    lfilename = filename
    if(present(output_index)) lfilename = trim(lfilename) // "." // int2str(output_index)
    lfilename = parallel_filename(trim(lfilename) // ".d", decomp)
    
    append_data = .false.
    if(present(append)) append_data = append
    if(append_data) then
      open(unit, file = trim(lfilename), access = "stream", form = "unformatted", status = "old", &
        & position = "append", iostat = ierr)
    else
      open_status = 'replace'
      if (present(overwrite)) then
        if (.not. overwrite) open_status = 'new'
      end if
      open(unit, file = trim(lfilename), access = "stream", form = "unformatted", status = open_status, iostat = ierr)
    end if
    if(ierr /= 0) then
      labort("Failed to open output data file: " // trim(lfilename))
    end if
    
  end subroutine open_parallel_data
    
  ! Write generic header information to an open parallel header file.
  subroutine write_parallel_header(decomp, unit, l_scale, t_scale)
    type(decomposition), intent(in) :: decomp
    integer, intent(in) :: unit
    real(kind = real_8), intent(in) :: l_scale
    real(kind = real_8), intent(in) :: t_scale
    
    write(unit, "(i0,a,i0,a,i0,a,i0,a,i0,a,i0)") decomp%n, " ", decomp%p_jj, " ", decomp%p_ii, " ", decomp%p, " ", decomp%p_j, " ", decomp%p_ii - decomp%p_i + 1
    write(unit, "(i0,a,i0,a,i0)") decomp%g_ii, " ", decomp%g_jj, " ", decomp%layers
    write(unit, "(i0,a,i0,a,i0,a,i0)") decomp%ii_0, " ", decomp%jj_0, " ", decomp%ii, " ", decomp%jj
    write(unit, "("//real_8_chr//",a,"//real_8_chr//")") l_scale, " ", t_scale
   
  end subroutine write_parallel_header
  
  ! Write header information to a snapshot header file.
  subroutine write_s_header(decomp, output_index, l_scale, t_scale, single_level, overwrite)
    use qg_parameters, only : unit
  
    type(decomposition), intent(in) :: decomp
    integer, intent(in) :: output_index
    real(kind = real_8), intent(in) :: l_scale
    real(kind = real_8), intent(in) :: t_scale
    logical, optional, intent(in) :: single_level,overwrite
    
    integer :: layer
    logical :: one_level
    character(16) :: dimensions
    
#ifdef CABARET
    one_level = .true.
#else
    one_level = .false.
#endif
    if (present(single_level)) one_level = single_level
    dimensions = '  ii jj'
    if (.not. one_level) dimensions = '  ii jj leapfrog'
    call open_parallel_header(decomp, unit, filename = "qg" // int2str(decomp%layers) // "_" // int2str(decomp%g_ii) // "_s", &
                              output_index = output_index, overwrite = overwrite)
    call write_parallel_header(decomp, unit, l_scale, t_scale)
    write(unit, "(a)") "model_time_day"
    write(unit, "(a)") "total_time_day"
    do layer = 1, decomp%layers
      write(unit, "(a,i0,a)") "psi", layer, trim(dimensions)
    end do
    do layer = 1, decomp%layers
      write(unit, "(a,i0,a)") "zeta", layer, trim(dimensions)
    end do
    close(unit)
    
  end subroutine write_s_header
  
#ifdef CABARET
  ! Write header information to a CABARET snapshot supplementary data file.
  subroutine write_s_supp_header(decomp, output_index, l_scale, t_scale)
    use qg_parameters, only : unit
  
    type(decomposition), intent(in) :: decomp
    integer, intent(in) :: output_index
    real(kind = real_8), intent(in) :: l_scale
    real(kind = real_8), intent(in) :: t_scale
    
    integer :: layer
    
    call open_parallel_header(decomp, unit, filename = "qg" // int2str(decomp%layers) // "_" // int2str(decomp%g_ii) // "_s_supp", output_index = output_index)
    call write_parallel_header(decomp, unit, l_scale, t_scale)
    write(unit, "(a)") "model_time_day"
    write(unit, "(a)") "total_time_day"
    write(unit, "(a)") "dt"
    do layer = 1, decomp%layers
      write(unit, "(a,i0,a)") "ui", layer, "  flux_i_ii flux_i_jj"
    end do
    do layer = 1, decomp%layers
      write(unit, "(a,i0,a)") "uj", layer, "  flux_j_ii flux_j_jj"
    end do
    do layer = 1, decomp%layers
      write(unit, "(a,i0,a)") "zeta_fi", layer, "  flux_i_ii flux_i_jj"
    end do
    do layer = 1, decomp%layers
      write(unit, "(a,i0,a)") "zeta_fj", layer, "  flux_j_ii flux_j_jj"
    end do
    do layer = 1, decomp%layers
      write(unit, "(a,i0,a)") "res", layer, "  ii jj"
    end do
    do layer = 1, decomp%layers
      write(unit, "(a,i0,a)") "res_beta", layer, "  ii jj"
    end do
    close(unit)
    
  end subroutine write_s_supp_header
#endif
  
  ! Write header information to a time-averaged header file.
  subroutine write_s_av_header(decomp, l_scale, t_scale)
    use qg_parameters, only : unit
  
    type(decomposition), intent(in) :: decomp
    real(kind = real_8), intent(in) :: l_scale
    real(kind = real_8), intent(in) :: t_scale
    
    integer :: layer
    
    call open_parallel_header(decomp, unit, filename = "qg" // int2str(decomp%layers) // "_" // int2str(decomp%g_ii) // "_s_av")
    call write_parallel_header(decomp, unit, l_scale, t_scale)
    write(unit, "(a)") "model_time_day"
    write(unit, "(a)") "total_time_day"
    do layer = 1, decomp%layers
      write(unit, "(a,i0,a)") "psi", layer, "_av  ii jj"
    end do
    close(unit)
    
  end subroutine write_s_av_header 
  
  ! Write header information to a time-averaged quadratic terms header file.
  subroutine write_s_quad_av_header(decomp, l_scale, t_scale)
    use qg_parameters, only : unit
  
    type(decomposition), intent(in) :: decomp
    real(kind = real_8), intent(in) :: l_scale
    real(kind = real_8), intent(in) :: t_scale
    
    call open_parallel_header(decomp, unit, filename = "qg" // int2str(decomp%layers) // "_" // int2str(decomp%g_ii) // "_quad_av")
    call write_parallel_header(decomp, unit, l_scale, t_scale)
    write(unit, "(a)") "model_time_day"
    write(unit, "(a)") "total_time_day"
    write(unit, "(a)") "uu_av       ii jj layer"
    write(unit, "(a)") "vv_av       ii jj layer"
    write(unit, "(a)") "uu_i_av     ii jj interface"
    write(unit, "(a)") "vv_i_av     ii jj interface"
    write(unit, "(a)") "uv_av       ii jj layer"
    write(unit, "(a)") "ur_av       ii jj layer"
    write(unit, "(a)") "vr_av       ii jj layer"
    write(unit, "(a)") "rr_av       ii jj layer"
    write(unit, "(a)") "pr_av       ii jj layer"
    write(unit, "(a)") "uz_av       ii jj layer"
    write(unit, "(a)") "vz_av       ii jj layer"
    write(unit, "(a)") "zz_av       ii jj layer"
    write(unit, "(a)") "zr_av       ii jj layer"
    write(unit, "(a)") "pp_av       ii jj layer"
    write(unit, "(a)") "pz_av       ii jj layer"
    write(unit, "(a)") "ut_av       ii jj layer"
    write(unit, "(a)") "vt_av       ii jj layer"
    write(unit, "(a)") "rt_av       ii jj layer"
    write(unit, "(a)") "tt_av       ii jj layer"
    write(unit, "(a)") "sc_ub_av    ii jj interface"
    write(unit, "(a)") "sc_vb_av    ii jj interface"
    write(unit, "(a)") "sc_bb_av    ii jj layer"
    write(unit, "(a)") "sc_bb_i_av  ii jj interface"
    close(unit)    
    
  end subroutine write_s_quad_av_header

  ! Read from a file open for stream input in a manner consistent with usual
  ! Fortran unformatted input.
  subroutine read_unformatted_rank0(unit, data)
    integer, intent(in) :: unit
    real(kind = real_8), intent(out) :: data
    
    integer :: s1, s2
    
    read(unit) s1
    assert(s1 == real_8_bit_size)
    read(unit) data
    read(unit) s2
    assert(s2 == s1)
  
  end subroutine read_unformatted_rank0
  
  ! Read from a file open for stream input in a manner consistent with usual
  ! Fortran unformatted input.
  subroutine read_unformatted_rank1(unit, data)
    integer, intent(in) :: unit
    real(kind = real_8), dimension(:), intent(out) :: data
    
    integer :: s1, s2
    
    read(unit) s1
    assert(s1 == size(data) * real_8_bit_size)
    read(unit) data
    read(unit) s2
    assert(s2 == s1)
  
  end subroutine read_unformatted_rank1

  ! Read from a file open for stream input in a manner consistent with usual
  ! Fortran unformatted input.
  subroutine read_unformatted_rank2(unit, data)
    integer, intent(in) :: unit
    real(kind = real_8), dimension(:, :), intent(out) :: data
    
    integer :: s1, s2
    
    read(unit) s1
    assert(s1 == size(data) * real_8_bit_size)
    read(unit) data
    read(unit) s2
    assert(s2 == s1)
  
  end subroutine read_unformatted_rank2
  
  ! Read from a file open for stream input in a manner consistent with usual
  ! Fortran unformatted input.
  subroutine read_unformatted_rank2_rank2(unit, data1, data2)
    integer, intent(in) :: unit
    real(kind = real_8), dimension(:, :), intent(out) :: data1
    real(kind = real_8), dimension(:, :), intent(out) :: data2

    integer :: s1, s2
    
    read(unit) s1
    assert(s1 == (size(data1) + size(data2)) * real_8_bit_size)
    read(unit) data1, data2
    read(unit) s2
    assert(s1 == s2)
  
  end subroutine read_unformatted_rank2_rank2
  
  ! Read from a file open for stream input in a manner consistent with usual
  ! Fortran unformatted input.
  subroutine read_unformatted_rank3(unit, data)
    integer, intent(in) :: unit
    real(kind = real_8), dimension(:, :, :), intent(out) :: data
    
    integer :: s1, s2
    
    read(unit) s1
    assert(s1 == size(data) * real_8_bit_size)
    read(unit) data
    read(unit) s2
    assert(s2 == s1)
  
  end subroutine read_unformatted_rank3

  ! Write to a file open for stream output in a manner consistent with usual
  ! Fortran unformatted output.
  subroutine write_unformatted_rank0(unit, data)
    integer, intent(in) :: unit
    real(kind = real_8), intent(in) :: data
    
    integer, parameter :: s = real_8_bit_size
    
    write(unit) s
    write(unit) data
    write(unit) s
    
  end subroutine write_unformatted_rank0

  ! Write to a file open for stream output in a manner consistent with usual
  ! Fortran unformatted output.
  subroutine write_unformatted_rank1(unit, data)
    integer, intent(in) :: unit
    real(kind = real_8), dimension(:), intent(in) :: data
    
    integer :: s
    
    s = size(data) * real_8_bit_size
    
    write(unit) s
    write(unit) data
    write(unit) s
    
  end subroutine write_unformatted_rank1
  
  ! Write to a file open for stream output in a manner consistent with usual
  ! Fortran unformatted output.
  subroutine write_unformatted_rank2(unit, data)
    integer, intent(in) :: unit
    real(kind = real_8), dimension(:, :), intent(in) :: data
    
    integer :: s
    
    s = size(data) * real_8_bit_size
    
    write(unit) s
    write(unit) data
    write(unit) s
    
  end subroutine write_unformatted_rank2
  
  ! Write to a file open for stream output in a manner consistent with usual
  ! Fortran unformatted output.
  subroutine write_unformatted_rank2_rank2(unit, data1, data2)
    integer, intent(in) :: unit
    real(kind = real_8), dimension(:, :), intent(in) :: data1
    real(kind = real_8), dimension(:, :), intent(in) :: data2
    
    integer :: s
    
    s = (size(data1) + size(data2)) * real_8_bit_size
    
    write(unit) s
    write(unit) data1, data2
    write(unit) s
    
  end subroutine write_unformatted_rank2_rank2
  
  ! Write to a file open for stream output in a manner consistent with usual
  ! Fortran unformatted output.
  subroutine write_unformatted_rank3(unit, data)
    integer, intent(in) :: unit
    real(kind = real_8), dimension(:, :, :), intent(in) :: data
    
    integer :: s
    
    s = size(data) * real_8_bit_size
    
    write(unit) s
    write(unit) data
    write(unit) s
    
  end subroutine write_unformatted_rank3

end module qg_io
