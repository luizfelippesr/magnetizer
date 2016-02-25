default=mpi
# IO options: IO_hdf5, IO_simple
# IO=IO_simple
IO=IO_hdf5
FC=h5pfc
FC_nonMPI=h5fc
srcdir=program
builddir=build
FCFLAGS=-I./${srcdir}/ -J./${builddir}/ -I./${builddir}/ -lfgsl -I/usr/local/include/fgsl -I/usr/include/ -fbacktrace  -ffpe-trap=zero,invalid,overflow -g

OBJ= bessel_functions.o root_finder.o constants.o grid.o global_input_parameters.o pressureEquilibrium.o outflow.o random.o  input_parameters.o $(IO).o profiles.o gutsdynamo.o ts_arrays.o  output.o dynamo.o

# Builds parallel version
mpi: $(OBJ) ${srcdir}/mpicalldynamo.f90
	$(FC) ${builddir}/$< $(FCFLAGS) -o magnetize_galform $(OBJ) mpicalldynamo.f90

# Builds serial version
serial: $(OBJ) calldynamo.f90
	$(FC)  -o magnetize_galform_serial $(OBJ) calldynamo.f90 $(FCFLAGS)

# Builds all objects/modules
%.o : ${srcdir}/%.f90 $(DEP)
	$(FC) -I$(srcdir) $(FCFLAGS) -c $^ -o ${builddir}/$@

# Explicit dependencies between files
$(srcdir)/pressureEquilibrium.f90: root_finder.o constants.o global_input_parameters.o
$(srcdir)/outflow.f90: input_parameters.o
$(srcdir)/grid.f90: constants.o
$(srcdir)/input_parameters.f90: grid.o
$(srcdir)/gutsdynamo.f90: pressureEquilibrium.o outflow.o profiles.o
$(srcdir)/ts_arrays.f90: gutsdynamo.o
$(srcdir)/dynamo.f90: ts_arrays.o
$(srcdir)/output.f90: $(IO).o gutsdynamo.o
$(srcdir)/$(IO).f90: grid.o
$(srcdir)/profiles.f90: pressureEquilibrium.o outflow.o

# Tides up
clean:
	rm -fv ${builddir}/*.mod ${builddir}/*.o
	rm -fv magnetize_galform_serial magnetize_galform
	rm -fv *.exe

# Test programs
testProfiles: pressureEquilibrium.o outflow.o ${srcdir}/tests/printProfiles.f90
	$(FC) -o printProfiles.exe ${builddir}/pressureEquilibrium.o ${builddir}/outflow.o tests/printProfiles.f90 $(FCFLAGS)
testRoots: root_finder.o ${srcdir}/tests/testRoot.f90
	$(FC) -o testRoot.exe ${builddir}/root_finder.o ${srcdir}/tests/testRoot.f90 $(FCFLAGS)
