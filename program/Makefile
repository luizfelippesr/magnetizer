default=mpi
# IO options: IO_hdf5, IO_simple
# IO=IO_simple
IO=IO_hdf5
# FC=mpif90
FC=h5pfc
FCFLAGS=-lfgsl -I/usr/local/include/fgsl -I/usr/include/ -fbacktrace  -ffpe-trap=zero,invalid,overflow -g
DEP=
COMMON_OBJ= bessel_functions.o constants.o grid.o global_input_parameters.o input_parameters.o outflow.o random.o root_finder.o pressureEquilibrium.o profiles.o gutsdynamo.o ts_arrays.o $(IO).o output.o dynamo.o

%.o: %.f90 $(DEP)
	$(FC) $(FCFLAGS) -c $^

mpi: $(COMMON_OBJ) mpicalldynamo.f90
	$(FC) -o magnetize_galform $(COMMON_OBJ) mpicalldynamo.f90 $(FCFLAGS)

serial: $(COMMON_OBJ) calldynamo.f90
	$(FC)  -o magnetize_galform_serial $(COMMON_OBJ) calldynamo.f90 $(FCFLAGS)

clean:
	rm -fv *.mod *.o
	rm -fv magnetize_galform_serial magnetize_galform
	rm -fv tests/*.exe

testProfiles: $(COMMON_OBJ) tests/printProfiles.f90
	$(FC) -o tests/printProfiles.exe $(COMMON_OBJ) tests/printProfiles.f90 $(FCFLAGS)

testRoots: $(COMMON_OBJ) tests/testRoot.f90
	$(FC) -o tests/testRoot.exe $(COMMON_OBJ) tests/testRoot.f90 $(FCFLAGS)
