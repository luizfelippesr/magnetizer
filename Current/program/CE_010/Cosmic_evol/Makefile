default=fmpi
DIR=Cosmic_evol

FC=mpif90
FCFLAGS=-lfgsl -I/usr/local/include/fgsl -I/usr/include/ -fbacktrace  -ffpe-trap=zero,invalid -g
DEP=
COMMON_OBJ= bessel_functions.o constants.o grid.o global_input_parameters.o input_parameters.o random.o profiles.o gutsdynamo.o IO_simple.o output.o dynamo.o

%.o: %.f90 $(DEP)
	$(FC) $(FCFLAGS) -c $^

fmpi: $(COMMON_OBJ) mpicalldynamo.f90
	$(FC) -o fmpicall $(COMMON_OBJ) mpicalldynamo.f90 $(FCFLAGS)

serial: $(COMMON_OBJ) calldynamo.f90
	$(FC)  -o fcall $(COMMON_OBJ) calldynamo.f90 $(FCFLAGS)

# pympi needs to be tested...
pympi: serial
	f2py -c -m ${FCFLAGS)} pycall constants.f90 grid.f90 global_input_parameters.f90 profiles.f90 input_parameters.f90 dynamo.f90 gutsdynamo.f90 random.f90

runserial: serial
	./fcall

clean:
	rm -fv *.mod *.o
	rm -fv fcall fmpicall

