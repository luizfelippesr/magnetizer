Instructions for compiling and running code

Using fortran (serial only)
In bash shell, from the directory "program/CE_###/" type 

  bash fcompile
  bash frun

or to remove '.mod' and '.o' files before compiling, type

  bash fcompileclean
  bash frun

Using python to wrap fortran code (mpi using mpi4py)
From the directory "program/CE_###/" type 

  mpicompile
  mpirun

List of files

In directory "program/CE_###/":

  fcompile		#compiles serial code
  frun                  #runs serial code
  fcompileclean         #compiles serial code with clean
  mpicompile            #compiles mpi code
  mpirun                #runs mpi code
  mpidynamo.py		#top-level python code that calls "Cosmic_evol/pycall.so", which wraps fortran routines
  diagnostic.out        #diagnostic file containing info from run
  r.pro			#idl routine that reads in simulation output
  p.pro			#idl routine that plots simulation output
  r.py			#python routine that reads in simulation output
  p.py			#python routine that plots simulation output

In directory "program/CE_###/Cosmic_evol/":

  __init__.py		#tells python to look in this directory
  __init__.pyc		#tells python to look in this directory
  fcall 		#executable created by 'fcompile' or 'fcompileclean'
  calldynamo.f90	#top-level fortran code that calls simulation, looping over the galaxies (this file is not necessary for mpi4py runs)
  dynamo.f90		#middle-level fortran code that runs the simulation for a given galaxy
  gutsdynamo.f90	#bottom-level fortran code that contains most of the necessary subroutines for dynamo.f90
  pycall.so		#executable created by f2py fortran wrapper for running fortran routines
  random.f90		#module for seeding variables with random noise
  noise.pro		#for testing the random.f90 module
  setup.sh		#was used for testing python
  mpiExample.py		#was used for testing mpi4py
