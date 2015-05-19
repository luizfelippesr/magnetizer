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

  bash mpicompile
  bash mpirun

Using fortran with mpi
From the directory "program/CE_###/" type 

  bash fmpicompile
  bash fmpirun

List of files:

In directory "program/":

  readme.txt		#this file
  rc.py			#generalizes program/CE_###/r.py to read output from multiple runs
  pc.py			#generalizes program/CE_###/p.py to plot output from multiple runs for comparison purposes
  rcm.py		#generalizes program/rc.py to read output from multiple galaxies from each run
  pcm.py		#generalizes program/pc.py to plot output from multiple galaxies from each run
  CE_###/

In directory "program/CE_###/":

  fcompile		#compiles serial code
  frun                  #runs serial code
  fcompileclean         #compiles serial code with clean
  mpicompile            #compiles python mpi code
  fmpicompile           #compiles fortran mpi code
  mpirun                #runs mpi code using python wrapper
  fmpirun               #runs mpi code completely within fortran
  mpidynamo.py		#top-level python code that calls "Cosmic_evol/pycall.so", which wraps fortran routines
  diagnostic.out        #diagnostic file containing info from run
  r.pro			#idl routine that reads in simulation output
  p.pro			#idl routine that plots simulation output
  r.py			#python routine that reads in simulation output
  p.py			#python routine that plots simulation output
  rm.py			#generalizes program/CE_###/r.py to read output from multiple galaxies
  pm.py			#generalizes program/CE_###/p.py to plot output from multiple galaxies
  Cosmic_evol/

In directory "program/CE_###/Cosmic_evol/":

  __init__.py		#tells python to look in this directory
  __init__.pyc		#tells python to look in this directory
  fcall 		#executable created by 'fcompile' or 'fcompileclean'
  calldynamo.f90	#top-level fortran code that calls simulation for serial runs without mpi
  mpicalldynamo.f90	#top-level fortran code that calls simulation for mpi runs within fortran
  dynamo.f90		#middle-level fortran code that runs the simulation for a given galaxy
  gutsdynamo.f90	#bottom-level fortran code that contains most of the necessary subroutines for dynamo.f90
  pycall.so		#executable created by f2py fortran wrapper for running fortran routines
  random.f90		#module for seeding variables with random noise
  noise.pro		#for testing the random.f90 module
  setup.sh		#was used for testing python
  mpiExample.py		#was used for testing mpi4py
