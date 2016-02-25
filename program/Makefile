default=mpi
# IO options: IO_hdf5, IO_simple
# IO=IO_simple
IO=IO_hdf5
# FC=mpif90
FC=h5pfc
FCFLAGS=-lfgsl -I/usr/local/include/fgsl -I/usr/include/ -fbacktrace  -ffpe-trap=zero,invalid,overflow -g
DEP=

OBJ= bessel_functions.o root_finder.o constants.o grid.o global_input_parameters.o pressureEquilibrium.o outflow.o random.o  input_parameters.o $(IO).o profiles.o gutsdynamo.o ts_arrays.o  output.o dynamo.o

%.o: %.f90 $(DEP)
	$(FC) $(FCFLAGS) -c $^

pressureEquilibrium.f90: root_finder.o constants.o global_input_parameters.o
outflow.f90: input_parameters.o
grid.f90: constants.o
input_parameters.f90: grid.o
gutsdynamo.f90: pressureEquilibrium.o outflow.o

common: $(OBJ)
	$(FC) $(FCFLAGS) -c $^

mpi: $(OBJ) mpicalldynamo.f90
	$(FC) -o magnetize_galform $(OBJ) mpicalldynamo.f90 $(FCFLAGS)

serial: $(OBJ) calldynamo.f90
	$(FC)  -o magnetize_galform_serial $(OBJ) calldynamo.f90 $(FCFLAGS)

clean:
	rm -fv *.mod *.o
	rm -fv magnetize_galform_serial magnetize_galform
	rm -fv tests/*.exe

testProfiles: pressureEquilibrium.o outflow.o tests/printProfiles.f90
	$(FC) -o printProfiles.exe pressureEquilibrium.o outflow.o tests/printProfiles.f90 $(FCFLAGS)

testRoots: root_finder.o tests/testRoot.f90
	$(FC) -o tests/testRoot.exe $(OBJ) tests/testRoot.f90 $(FCFLAGS)
