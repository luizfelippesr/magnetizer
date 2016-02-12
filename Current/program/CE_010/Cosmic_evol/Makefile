default=fmpi
# IO=IO_hdf5
IO=IO_simple
# IO options: IO_hdf5, IO_simple 
FC=mpif90
# FC=h5pfc
FCFLAGS=-lfgsl -I/usr/local/include/fgsl -I/usr/include/ -fbacktrace  -ffpe-trap=zero,invalid -g
DEP=
COMMON_OBJ= bessel_functions.o constants.o grid.o global_input_parameters.o input_parameters.o outflow.o random.o pressureEquilibrium.o profiles.o gutsdynamo.o ts_arrays.o $(IO).o output.o dynamo.o

%.o: %.f90 $(DEP)
	$(FC) $(FCFLAGS) -c $^

fmpi: $(COMMON_OBJ) mpicalldynamo.f90
	$(FC) -o fmpicall $(COMMON_OBJ) mpicalldynamo.f90 $(FCFLAGS)

serial: $(COMMON_OBJ) calldynamo.f90
	$(FC)  -o fcall $(COMMON_OBJ) calldynamo.f90 $(FCFLAGS)

# pympi needs to be tested...
pympi: serial
	f2py -c -m ${FCFLAGS)} pycall constants.f90 grid.f90 global_input_parameters.f90 profiles.f90 input_parameters.f90 dynamo.f90 gutsdynamo.f90 random.f90

clean:
	rm -fv *.mod *.o
	rm -fv fcall fmpicall

