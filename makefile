include makefile.compiler

DEFINES =

OBJS = lapack_interfaces.o qg3p_aux.o qg3p_base.o qg3p_diagnostics.o \
  qg3p_helmholtz.o qg3p_helmholtz_petsc.o qg3p_io.o qg3p_laplace_bc.o \
  qg3p_params.o qg3p_geometry.o qg3p_irregular_parallel.o qg3p_sha256.o qg3p_cmm.o test_tools.o
OBJS += qg3p_cmm_dense.o qg3p_constraints.o qg3p_irregular_diagnostics.o qg3p_irregular_bc.o
OBJS += qg3p_moments.o qg3p_irregular_operators.o qg3p_reconstruction.o qg3p_smooth_corners.o
OBJS += qg3p_weak_diffusion.o qg3p_island_constraints.o qg3p_irregular_transport.o qg3p_irregular_forcing.o
MODS = lapack_interfaces.mod qg_diagnostics.mod qg_aux.mod qg_helmholtz.mod \
  qg_helmholtz_petsc.mod qg_io.mod qg_laplace_bc.mod qg_parallel.mod \
  qg_parameters.mod qg_geometry.mod qg_irregular_parallel.mod qg_sha256.mod qg_cmm.mod test_tools.mod
MODS += qg_cmm_dense.mod qg_constraints.mod qg_irregular_diagnostics.mod qg_irregular_bc.mod
MODS += qg_moments.mod qg_irregular_operators.mod qg_reconstruction.mod qg_smooth_corners.mod
MODS += qg_weak_diffusion.mod qg_island_constraints.mod qg_irregular_transport.mod qg_irregular_forcing.mod
CONFIG_HEADERS = config.h version.h

CAPE_OBJECTS = qg3p_params.o test_tools.o qg3p_base.o lapack_interfaces.o qg3p_diagnostics.o \
  qg3p_helmholtz.o qg3p_helmholtz_petsc.o qg3p_laplace_bc.o qg3p_sha256.o qg3p_geometry.o \
  qg3p_moments.o qg3p_irregular_parallel.o qg3p_cmm_dense.o qg3p_cmm.o qg3p_constraints.o \
  qg3p_reconstruction.o qg3p_irregular_bc.o qg3p_irregular_operators.o qg3p_irregular_forcing.o \
  qg3p_irregular_transport.o qg3p_io.o qg3p_cape_slip.o
CAPE_BUILD ?= build/cape-native
CAPE_FC ?= mpif90
CAPE_FCFLAGS ?= -O2 -g -std=gnu -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -fallow-argument-mismatch
CAPE_CONFIG ?= examples/cape.nml
CAPE_RUN_DIR ?= run-cape
CAPE_RANKS ?= 4
MPIEXEC ?= mpirun
  
default: qg3
qg3: $(OBJS) qg3.o
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)
clean: clean-bench clean-unittest clean-manual
	rm -fv qg3 $(MODS) $(OBJS) qg3.o
    
.IGNORE: clean-bench clean-unittest clean-manual
.PHONY: bench manual test unittest cmm-check cape run-cape
PYTHON ?= python3
CMM_CHECK_ARGS ?=
%.o: %.F90
	$(FC) $(FCFLAGS) -c $< -o $@

lapack_interfaces.o: $(CONFIG_HEADERS) debug.h qg3p_base.o qg3p_params.o
qg3p_aux.o: $(CONFIG_HEADERS) debug.h leapfrog.h lapack_interfaces.o \
  qg3p_base.o qg3p_diagnostics.o qg3p_helmholtz.o qg3p_helmholtz_petsc.o \
  qg3p_laplace_bc.o qg3p_params.o test_tools.o
qg3p_base.o: $(CONFIG_HEADERS) debug.h qg3p_params.o test_tools.o
qg3p_diagnostics.o: $(CONFIG_HEADERS) debug.h lapack_interfaces.o qg3p_base.o \
  qg3p_params.o test_tools.o
qg3p_helmholtz.o: $(CONFIG_HEADERS) debug.h lapack_interfaces.o qg3p_base.o \
  qg3p_params.o
qg3p_helmholtz_petsc.o: $(CONFIG_HEADERS) debug.h qg3p_base.o \
  qg3p_laplace_bc.o qg3p_params.o
qg3p_io.o: $(CONFIG_HEADERS) debug.h qg3p_base.o qg3p_diagnostics.o \
  qg3p_params.o test_tools.o
qg3p_laplace_bc.o: $(CONFIG_HEADERS) debug.h qg3p_base.o qg3p_params.o
qg3p_params.o: $(CONFIG_HEADERS)
qg3p_geometry.o: qg3p_params.o qg3p_sha256.o
qg3p_irregular_parallel.o: $(CONFIG_HEADERS) debug.h qg3p_geometry.o qg3p_base.o qg3p_moments.o
qg3p_moments.o: qg3p_params.o
qg3p_irregular_operators.o: qg3p_irregular_parallel.o qg3p_reconstruction.o
qg3p_reconstruction.o: qg3p_irregular_parallel.o qg3p_moments.o lapack_interfaces.o
qg3p_smooth_corners.o: qg3p_reconstruction.o qg3p_diagnostics.o
qg3p_irregular_transport.o qg3p_irregular_forcing.o: qg3p_irregular_operators.o
qg3p_cape_slip.o: qg3p_constraints.o qg3p_irregular_operators.o
qg3p_weak_diffusion.o: qg3p_cmm.o qg3p_irregular_operators.o qg3p_irregular_transport.o
qg3p_island_constraints.o: qg3p_constraints.o qg3p_irregular_operators.o lapack_interfaces.o
qg3p_cmm_dense.o: $(CONFIG_HEADERS) qg3p_base.o
qg3p_cmm.o: $(CONFIG_HEADERS) qg3p_cmm_dense.o lapack_interfaces.o qg3p_irregular_parallel.o \
  qg3p_helmholtz.o qg3p_helmholtz_petsc.o qg3p_laplace_bc.o
qg3p_constraints.o: qg3p_cmm.o qg3p_diagnostics.o
qg3p_irregular_diagnostics.o: $(CONFIG_HEADERS) debug.h qg3p_irregular_parallel.o qg3p_diagnostics.o
qg3p_irregular_bc.o: $(CONFIG_HEADERS) debug.h qg3p_irregular_parallel.o qg3p_diagnostics.o qg3p_reconstruction.o
qg3.o: $(CONFIG_HEADERS) debug.h leapfrog.h qg3p_aux.o qg3p_base.o \
  qg3p_diagnostics.o qg3p_io.o qg3p_helmholtz.o qg3p_params.o
qg3_cape.o: $(CAPE_OBJECTS)
qg3_cape: $(CAPE_OBJECTS) qg3_cape.o
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)
test_tools.o: test_tools.F90 $(CONFIG_HEADERS) debug.h qg3p_params.o

test:
	test/run_tests
			
bench: $(OBJS)
	test -d bench && make -C bench -I..
clean-bench:
	test -d bench && make clean -C bench -I..
			
unittest-bin: $(OBJS)
	test -d unittest && $(MAKE) default -C unittest -I.. MODEL_OBJECTS="$(addprefix ../,$(OBJS))"
unittest: unittest-bin
	test -d unittest && make unittest -C unittest -I..
cmm-check:
	$(PYTHON) unittest/run_cmm_checks.py $(CMM_CHECK_ARGS)

# Copy sources into a native-only build so a cape experiment can coexist
# with the normal qg3 objects, compiler settings and simulation files.
cape:
	mkdir -p "$(CAPE_BUILD)"
	set -e; for src in $(CAPE_OBJECTS:.o=.F90) qg3_cape.F90 version.h debug.h leapfrog.h \
	    makefile makefile.compiler makefile.gfortran; do \
	  cmp -s "$$src" "$(CAPE_BUILD)/$$src" || cp "$$src" "$(CAPE_BUILD)/$$src"; \
	done
	sed 's/^#define PETSC$$/#undef PETSC/' config.h > "$(CAPE_BUILD)/config.h.tmp"
	if cmp -s "$(CAPE_BUILD)/config.h.tmp" "$(CAPE_BUILD)/config.h"; then \
	  rm "$(CAPE_BUILD)/config.h.tmp"; \
	else mv "$(CAPE_BUILD)/config.h.tmp" "$(CAPE_BUILD)/config.h"; fi
	$(MAKE) -C "$(CAPE_BUILD)" FC="$(CAPE_FC)" LD="$(CAPE_FC)" FCFLAGS="$(CAPE_FCFLAGS)" \
	  LDFLAGS= LIBS='-llapack -lblas' qg3_cape

run-cape: cape
	mkdir -p "$(CAPE_RUN_DIR)"
	cd "$(CAPE_RUN_DIR)" && $(MPIEXEC) -np $(CAPE_RANKS) "$(abspath $(CAPE_BUILD))/qg3_cape" "$(abspath $(CAPE_CONFIG))"

clean-unittest:
	test -d unittest && make clean -C unittest -I..

manual:
	make -C manual
clean-manual:
	test -d manual && make clean -C manual
