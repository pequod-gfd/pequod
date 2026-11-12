include ../makefile.compiler

UNITTESTS = $(basename $(wildcard *.F90))
MODEL_OBJECTS ?= $(filter-out ../qg3.o, $(wildcard ../*.o))

default: $(UNITTESTS)
clean:
	rm -fv $(UNITTESTS) $(addsuffix .o, $(UNITTESTS))

unittest:
	./run_tests $(UNITTESTS)

$(UNITTESTS): %: %.F90 $(MODEL_OBJECTS)
	$(FC) $(FCFLAGS) -c -o $@.o $< -I..
	$(LD) $(LDFLAGS) -o $@ $@.o $(MODEL_OBJECTS) $(LIBS)
