default=mpi
# IO options: IO_hdf5, IO_simple
# IO=IO_simple
IO=IO_hdf5
FC=h5pfc
FC_nonMPI=h5fc
srcdir=program
builddir=build
FCFLAGS=-I. -I./${srcdir}/ -J./${builddir}/ -fintrinsic-modules-path ./${builddir} -I./${builddir}/ -lfgsl -I/usr/local/include/fgsl -I/usr/include/ -fbacktrace  -ffpe-trap=zero,invalid,overflow -g
VPATH=${srcdir}:${builddir}
_OBJ= bessel_functions.o root_finder.o constants.o grid.o global_input_parameters.o pressureEquilibrium.o outflow.o random.o  input_parameters.o $(IO).o profiles.o gutsdynamo.o ts_arrays.o  output.o dynamo.o
OBJ = $(patsubst %,$(builddir)/%,$(_OBJ))

# h5pfc bessel_functions.o -I. -I./program/ -J./build/ -fintrinsic-modules-path ./build -I./build/ -lfgsl -I/usr/local/include/fgsl -I/usr/include/ -fbacktrace  -ffpe-trap=zero,invalid,overflow -g -o magnetize_galform.exe  program/mpicalldynamo.f90

# Builds parallel version
mpi: $(OBJ) ${srcdir}/mpicalldynamo.f90
	$(FC) $^ $(FCFLAGS) -o magnetize_galform.exe

# Builds serial version
serial: $(OBJ) ${srcdir}/calldynamo.f90
	$(FC_nonMPI) $< $(FCFLAGS) -o magnetize_galform_serial.exe $(OBJ) ${srcdir}/calldynamo.f90

# Test programs
testProfiles: pressureEquilibrium.o tests.printProfiles.o global_input_parameters.o root_finder.o
	$(FC) $(FCFLAGS) ${builddir}/pressureEquilibrium.o ${builddir}/tests.printProfiles.o ${builddir}/global_input_parameters.o ${builddir}/root_finder.o  -o printProfiles.exe
testRoots: root_finder.o tests.testRoots.o
	$(FC) ${builddir}/root_finder.o ${builddir}/tests.testRoots.o $(FCFLAGS) -o testRoots.exe
# All programs
all: testRoots testProfiles mpi

# Builds all objects/modules following
${builddir}/%.o : ${srcdir}/%.f90 $(DEP)
	$(FC)  $(FCFLAGS) -c $^ -o $@

# Explicit dependencies between files
$(srcdir)/pressureEquilibrium.f90: ${builddir}/root_finder.o ${builddir}/constants.o ${builddir}/global_input_parameters.o
$(srcdir)/outflow.f90: ${builddir}/input_parameters.o
$(srcdir)/grid.f90: ${builddir}/constants.o
$(srcdir)/input_parameters.f90: ${builddir}/grid.o
$(srcdir)/gutsdynamo.f90: ${builddir}/pressureEquilibrium.o ${builddir}/outflow.o ${builddir}/profiles.o
$(srcdir)/ts_arrays.f90: ${builddir}/gutsdynamo.o
$(srcdir)/dynamo.f90: ${builddir}/ts_arrays.o
$(srcdir)/output.f90: ${builddir}/$(IO).o ${builddir}/gutsdynamo.o
$(srcdir)/$(IO).f90: ${builddir}/grid.o
$(srcdir)/profiles.f90: ${builddir}/pressureEquilibrium.o ${builddir}/outflow.o

# Tides up
clean:
	rm -fv ${builddir}/*.mod ${builddir}/*.o
	rm -fv magnetize_galform_serial magnetize_galform
	rm -fv *.exe
