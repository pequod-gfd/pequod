include ../makefile.compiler

UNITTESTS = $(basename $(wildcard *.F90))
ARG_FIXTURES = test_cmm_dirichlet test_irregular_mass test_irregular_mpi_geometry test_irregular_straight_wall test_irregular_operators
ARG_FIXTURES += test_irregular_vertex_order test_irregular_smooth_corners
ARG_FIXTURES += test_cmm_scalapack test_cmm_profile test_weak_diffusion test_irregular_island_constraints test_irregular_transport test_irregular_forcing
ARG_FIXTURES += test_cape_cabaret
MODEL_OBJECTS ?= $(filter-out ../qg3.o, $(wildcard ../*.o))

default: $(UNITTESTS)
clean:
	rm -fv $(UNITTESTS) $(addsuffix .o, $(UNITTESTS))

unittest:
	@echo "SKIP $(ARG_FIXTURES) (run with make cmm-check)"
	./run_tests $(filter-out $(ARG_FIXTURES),$(UNITTESTS))

$(UNITTESTS): %: %.F90 $(MODEL_OBJECTS)
	$(FC) $(FCFLAGS) -c -o $@.o $< -I..
	$(LD) $(LDFLAGS) -o $@ $@.o $(MODEL_OBJECTS) $(LIBS)
