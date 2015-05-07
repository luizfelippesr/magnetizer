#!/bin/bash

# Compile module files
gfortran -c random.f90 dynamoguts.f90 topdynamo.f90
# Compile F code
gfortran -o call calltopdynamo.f90 topdynamo.o dynamoguts.o random.o

# Compile python module
f2py -c -m pycall calltopdynamo.f90 topdynamo.f90 dynamoguts.f90 random.f90

# Make directory behave like a module
touch __init__.py

# Go up one
cd ..

# Test import
python -c 'import cosmic_evol.pycall as pycall'
