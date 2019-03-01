default=prod
FC=h5pfc
srcdir=source
builddir=build

_OBJ= tsDataObj.o bessel_functions.o root_finder.o constants.o grid.o floor_field.o global_input_parameters.o pressureEquilibrium.o outflow.o random.o  input_parameters.o IO_hdf5.o profiles.o gutsdynamo.o ts_arrays.o  output.o dynamo.o rotationCurves.o deriv.o messages.o interpolation.o integration.o seed_field.o
OBJ = $(patsubst %,$(builddir)/%,$(_OBJ))

FCFLAGS+= -lfgsl -I. -I./${srcdir}/ -J./${builddir}/ -fintrinsic-modules-path ./${builddir} -I./${builddir}/ -I/usr/include/  -fbacktrace  -ffpe-trap=zero,invalid,overflow -fbounds-check

FCFLAGS_TEST=-g -Wall
FCFLAGS_PROD=-O2

help:
	@echo '---------------------'
	@echo ' Magnetizer Makefile'
	@echo '---------------------'
	@echo '  make prod -> Builds main program for a production run.'
	@echo '  make test -> Builds main program for a test run (debugging and backtracing enabled).'
	@echo '  make clean -> Removes all object files.'
	@echo '  make cleanall -> Removes all object files, executables files and python bytecode.'
	@echo '  make all -> Builds all test programs and main program (test settings).'
	@echo '  make help -> Displays this help'

prod: PRD main

test: TST main

PRD:
	@echo ''
	@echo 'Building Magnetizer production run'
	@echo ''
	$(eval FCFLAGS += $(FCFLAGS_PROD))
TST:
	@echo ''
	@echo 'Building Magnetizer test/debug run'
	@echo ''
	$(eval FCFLAGS += $(FCFLAGS_TEST))

main: $(OBJ) ${builddir}/Magnetizer.o
	$(FC) $^ $(FCFLAGS) -o Magnetizer.exe


integrate: $(OBJ) ${builddir}/path_integrate.o
	$(FC) $^ $(FCFLAGS) -o path_integrate.exe

# Other programs
LoSintegrate: TST ${builddir}/constants.o ${builddir}/LoSintegrate_aux.o ${builddir}/messages.o ${builddir}/global_input_parameters.o ${builddir}/interpolation.o ${builddir}/grid.o  ${builddir}/IO_hdf5.o ${builddir}/LoSintegrate.o
	$(FC) $(FCFLAGS) ${builddir}/constants.o ${builddir}/LoSintegrate_aux.o ${builddir}/messages.o ${builddir}/global_input_parameters.o ${builddir}/interpolation.o ${builddir}/grid.o  ${builddir}/IO_hdf5.o ${builddir}/LoSintegrate.o  -o LoSintegrate.exe

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
all: test testRoots testProfiles

# Tides up
clean:
	rm -fv ${builddir}/*.mod ${builddir}/*.o Magnetizer.exe

cleanall: clean
	rm -fv python/*.pyc
	rm -fv python/*/*.pyc
	rm -fv python/*/*/*.pyc
	rm -fv *.exe

# Builds all objects/modules following
${builddir}/%.o : ${srcdir}/%.f90
	$(FC)  $(FCFLAGS) -c $^ -o $@

# Explicit dependencies between files
$(srcdir)/pressureEquilibrium.f90: ${builddir}/root_finder.o ${builddir}/constants.o ${builddir}/global_input_parameters.o ${builddir}/input_parameters.o ${builddir}/integration.o
$(srcdir)/outflow.f90: ${builddir}/input_parameters.o
$(srcdir)/grid.f90: ${builddir}/constants.o ${builddir}/global_input_parameters.o ${builddir}/messages.o ${builddir}/interpolation.o
$(srcdir)/deriv.f90: ${builddir}/grid.o
$(srcdir)/LoSintegrate_aux.f90: ${builddir}/global_input_parameters.o ${builddir}/interpolation.o ${builddir}/messages.o
$(srcdir)/input_parameters.f90: ${builddir}/grid.o ${builddir}/IO_hdf5.o
$(srcdir)/gutsdynamo.f90: ${builddir}/pressureEquilibrium.o ${builddir}/outflow.o ${builddir}/profiles.o ${builddir}/deriv.o ${builddir}/floor_field.o ${builddir}/seed_field.o
$(srcdir)/ts_arrays.f90: ${builddir}/gutsdynamo.o ${builddir}/interpolation.o ${builddir}/tsDataObj.o
$(srcdir)/dynamo.f90: ${builddir}/output.o ${builddir}/messages.o ${builddir}/ts_arrays.o ${builddir}/interpolation.o ${builddir}/floor_field.o
$(srcdir)/output.f90: ${builddir}/IO_hdf5.o ${builddir}/gutsdynamo.o ${builddir}/ts_arrays.o ${builddir}/tsDataObj.o
$(srcdir)/IO_hdf5.f90: ${builddir}/grid.o ${builddir}/messages.o
$(srcdir)/profiles.f90: ${builddir}/pressureEquilibrium.o ${builddir}/outflow.o ${builddir}/rotationCurves.o ${builddir}/input_parameters.o ${builddir}/grid.o
$(srcdir)/Magnetizer.f90: ${builddir}/dynamo.o ${builddir}/grid.o ${builddir}/messages.o
$(srcdir)/rotationCurves.f90: ${builddir}/bessel_functions.o ${builddir}/deriv.o
$(srcdir)/floor_field.f90: ${builddir}/grid.o ${builddir}/global_input_parameters.o
$(srcdir)/seed_field.f90: ${builddir}/grid.o ${builddir}/profiles.o ${builddir}/global_input_parameters.o
