include makefile.compiler

DEFINES =

OBJS = lapack_interfaces.o qg3p_aux.o qg3p_base.o qg3p_diagnostics.o \
  qg3p_helmholtz.o qg3p_helmholtz_petsc.o qg3p_io.o qg3p_laplace_bc.o \
  qg3p_params.o qg3p_geometry.o qg3p_irregular_parallel.o qg3p_sha256.o qg3p_cmm.o test_tools.o
OBJS += qg3p_constraints.o qg3p_irregular_diagnostics.o qg3p_irregular_bc.o
MODS = lapack_interfaces.mod qg_diagnostics.mod qg_aux.mod qg_helmholtz.mod \
  qg_helmholtz_petsc.mod qg_io.mod qg_laplace_bc.mod qg_parallel.mod \
  qg_parameters.mod qg_geometry.mod qg_irregular_parallel.mod qg_sha256.mod qg_cmm.mod test_tools.mod
MODS += qg_constraints.mod qg_irregular_diagnostics.mod qg_irregular_bc.mod
CONFIG_HEADERS = config.h version.h
  
default: qg3
qg3: $(OBJS) qg3.o
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)
clean: clean-bench clean-unittest clean-manual
	rm -fv qg3 $(MODS) $(OBJS) qg3.o
    
.IGNORE: clean-bench clean-unittest clean-manual
.PHONY: bench manual test unittest
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
qg3p_irregular_parallel.o: $(CONFIG_HEADERS) debug.h qg3p_geometry.o qg3p_base.o
qg3p_cmm.o: $(CONFIG_HEADERS) lapack_interfaces.o qg3p_irregular_parallel.o \
  qg3p_helmholtz.o qg3p_helmholtz_petsc.o qg3p_laplace_bc.o
qg3p_constraints.o: qg3p_cmm.o qg3p_diagnostics.o
qg3p_irregular_diagnostics.o: $(CONFIG_HEADERS) debug.h qg3p_irregular_parallel.o qg3p_diagnostics.o
qg3p_irregular_bc.o: $(CONFIG_HEADERS) debug.h qg3p_irregular_parallel.o qg3p_diagnostics.o
qg3.o: $(CONFIG_HEADERS) debug.h leapfrog.h qg3p_aux.o qg3p_base.o \
  qg3p_diagnostics.o qg3p_io.o qg3p_helmholtz.o qg3p_params.o
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
clean-unittest:
	test -d unittest && make clean -C unittest -I..

manual:
	make -C manual
clean-manual:
	test -d manual && make clean -C manual
