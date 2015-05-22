#!/usr/bin/env python

"""
Run as:
mpirun -np NPROCS ./mpiExample.py
"""

import sys
from mpi4py import MPI
#import cosmic_evol.pycall as pycall

world=MPI.COMM_WORLD
rank=world.rank
size=world.size

gals=range(1,101)
my_gals=gals[rank::size]

# Loop over galaxies
#PYTCH_INFO=0
for igal in my_gals:
    print "Processor %d handling file %s" % (rank,igal)
    # Do some work
    #flag=pycall.topdynamo.dynamo_run(PYTCH_INFO,igal)
    #print flag
    pass

sys.exit(0)
