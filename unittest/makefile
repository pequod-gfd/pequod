include ../makefile.compiler

UNITTESTS = $(basename $(wildcard *.F90))

default: $(UNITTESTS)
clean:
	rm -fv $(UNITTESTS) $(addsuffix .o, $(UNITTESTS))

unittest:
	./run_tests $(UNITTESTS)

%: %.F90
	$(FC) $(FCFLAGS) -c -o $@.o $< -I..
	$(LD) $(LDFLAGS) -o $@ $@.o $(filter-out ../qg3.o, $(wildcard ../*.o)) $(LIBS)
