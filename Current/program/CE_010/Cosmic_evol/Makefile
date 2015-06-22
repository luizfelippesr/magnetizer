default=fmpi
DIR=Cosmic_evol

FC=mpif90
FCFLAGS=
DEP=
COMMON_OBJ=constants.o parameters.o random.o gutsdynamo.o dynamo.o

%.o: %.f90 $(DEP)
	$(FC) $(FCFLAGS) -c $^

fmpi: $(COMMON_OBJ) mpicalldynamo.f90
	$(FC) -o fmpicall $(COMMON_OBJ) mpicalldynamo.f90

serial: $(COMMON_OBJ) calldynamo.f90
	$(FC) -o fcall $(COMMON_OBJ) calldynamo.f90

pympi: serial
	f2py -c -m pycall constants.f90 parameters.f90 dynamo.f90 gutsdynamo.f90 random.f90

runserial: serial
	./fcall

clean:
	rm -fv *.mod *.o
	rm -fv fcall fmpicall

