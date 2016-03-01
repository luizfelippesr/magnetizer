#!/usr/bin/env python

"""
Run as:
mpirun -np NPROCS ./mpiExample.py
"""
import timeit

start = timeit.default_timer()

import sys
from mpi4py import MPI
import Cosmic_evol.pycall as pycall

world=MPI.COMM_WORLD
rank=world.rank
size=world.size

#gals=range(1,101)
gals=range(1,9)
my_gals=gals[rank::size]

# intent(in)
#  info:   how much info to print during run
#          options are 0=min, 1=standard, 2=max
#  gal_id: id number of galaxy
#          ranges from 1 to up to 99,999,999
#
# intent(out)
#  flag:   status of run
#          results are -1=run successful, 1=, 2=, 3=
#
# Note: calculated physical variables (ts_t, ts_Br, ts_Bp, ts_Bzmod, ts_alp_m) are written to files, not passed
#
# Fortran format:  call dynamo_run(info, gal_id, flag)
#  Standard parameters: (0, igal, flag), where igal is looped over within python

# Loop over galaxies
PY_INFO=0
for igal in my_gals:
    print "Processor %d handling file %s" % (rank,igal)
    # Do some work
    flag=pycall.dynamo.dynamo_run(PY_INFO,igal)
    print 'flag=',flag,'for galaxy',igal
    pass

world.Barrier()

if rank==0:
  stop = timeit.default_timer()
  print "Total simulation time="
  print stop-start

sys.exit(0)
