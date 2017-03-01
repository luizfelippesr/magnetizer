default=mpi
# IO options: IO_hdf5, IO_simple
# IO=IO_simple
IO=IO_hdf5
FC=h5pfc
FC_nonMPI=h5fc
srcdir=source
builddir=build

_OBJ= bessel_functions.o root_finder.o constants.o grid.o global_input_parameters.o pressureEquilibrium.o outflow.o random.o  input_parameters.o $(IO).o profiles.o gutsdynamo.o ts_arrays.o  output.o dynamo.o rotationCurves.o deriv.o messages.o interpolation.o integration.o
OBJ = $(patsubst %,$(builddir)/%,$(_OBJ))

# Special configurations
ifeq ($(HOSTNAME),topsy-login-1-0.local)
	FCFLAGS_special= -fintrinsic-modules-path /share/apps/fgsl/fgsl-1.0.0/include/fgsl -I/share/apps/fgsl/fgsl-1.0.0/include/fgsl
endif

FCFLAGS+=-I. -I./${srcdir}/ -J./${builddir}/ -fintrinsic-modules-path ./${builddir} -I./${builddir}/ -lfgsl  -I/usr/include/ ${FCFLAGS_special} -fbacktrace  -ffpe-trap=zero,invalid,overflow -fbounds-check -g -Wall -Wuninitialized

# Builds parallel version
mpi: $(OBJ) ${builddir}/mpicalldynamo.o
	$(FC) $^ $(FCFLAGS) -o magnetize_galform.exe

# Builds serial version (currently not working!)
serial: $(OBJ) ${builddir}/calldynamo.o
	$(FC_nonMPI) $^ $(FCFLAGS) -o magnetize_galform_serial.exe

# Test programs
testProfiles: ${builddir}/pressureEquilibrium.o ${builddir}/tests.printProfiles.o ${builddir}/global_input_parameters.o ${builddir}/root_finder.o
	$(FC) $(FCFLAGS) ${builddir}/pressureEquilibrium.o ${builddir}/tests.printProfiles.o ${builddir}/global_input_parameters.o ${builddir}/root_finder.o  -o printProfiles.exe
testInterpolation: ${builddir}/interpolation.o ${builddir}/tests.interpolation.o
	$(FC) $(FCFLAGS) ${builddir}/interpolation.o ${builddir}/tests.interpolation.o -o testInterpolation.exe
testRoots: ${builddir}/root_finder.o ${builddir}/tests.testRoots.o
	$(FC) ${builddir}/root_finder.o ${builddir}/tests.testRoots.o $(FCFLAGS) -o testRoots.exe
testPressureEquilibrium: ${builddir}/constants.o ${builddir}/messages.o ${builddir}/global_input_parameters.o ${builddir}/pressureEquilibrium.o ${builddir}/tests.pressureEquilibrium.o
	$(FC) ${builddir}/constants.o ${builddir}/messages.o ${builddir}/global_input_parameters.o ${builddir}/root_finder.o ${builddir}/pressureEquilibrium.o ${builddir}/tests.pressureEquilibrium.o $(FCFLAGS) -o testPressureEquilibrium.exe
testIntegration: ${builddir}/integration.o ${builddir}/tests.integration.o
	$(FC) $(FCFLAGS) ${builddir}/integration.o ${builddir}/tests.integration.o -o testIntegration.exe
# All programs
all: testRoots testProfiles mpi

# Builds all objects/modules following
${builddir}/%.o : ${srcdir}/%.f90
	echo ${HOSTNAME}
	$(FC)  $(FCFLAGS) -c $^ -o $@

# Explicit dependencies between files
$(srcdir)/pressureEquilibrium.f90: ${builddir}/root_finder.o ${builddir}/constants.o ${builddir}/global_input_parameters.o ${builddir}/input_parameters.o ${builddir}/integration.o
$(srcdir)/outflow.f90: ${builddir}/input_parameters.o
$(srcdir)/grid.f90: ${builddir}/constants.o ${builddir}/global_input_parameters.o ${builddir}/messages.o ${builddir}/interpolation.o
$(srcdir)/deriv.f90: ${builddir}/grid.o
$(srcdir)/input_parameters.f90: ${builddir}/grid.o ${builddir}/$(IO).o
$(srcdir)/gutsdynamo.f90: ${builddir}/pressureEquilibrium.o ${builddir}/outflow.o ${builddir}/profiles.o ${builddir}/deriv.o
$(srcdir)/ts_arrays.f90: ${builddir}/gutsdynamo.o ${builddir}/interpolation.o
$(srcdir)/dynamo.f90: ${builddir}/output.o ${builddir}/messages.o ${builddir}/ts_arrays.o ${builddir}/interpolation.o
$(srcdir)/output.f90: ${builddir}/$(IO).o ${builddir}/gutsdynamo.o ${builddir}/ts_arrays.o
$(srcdir)/$(IO).f90: ${builddir}/grid.o ${builddir}/messages.o
$(srcdir)/profiles.f90: ${builddir}/pressureEquilibrium.o ${builddir}/outflow.o ${builddir}/rotationCurves.o ${builddir}/input_parameters.o ${builddir}/grid.o
$(srcdir)/mpicalldynamo.f90: ${builddir}/dynamo.o ${builddir}/grid.o ${builddir}/messages.o
$(srcdir)/rotationCurves.f90: ${builddir}/bessel_functions.o ${builddir}/deriv.o

# Tides up
clean:
	rm -fv ${builddir}/*.mod ${builddir}/*.o

cleanall: clean
	rm -fv scripts/*.pyc
	rm -fv scripts/*/*.pyc
	rm -fv *.exe
