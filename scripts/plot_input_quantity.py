#!/usr/bin/env python
import h5py
import sys
import numpy as N
import pylab as P
import argparse

#filename = sys.argv[1]
#igals = N.array(sys.argv[2:]).astype(int)

#print igals




if __name__ == "__main__"  :
    parser = argparse.ArgumentParser(
             description='This script plots the inputs of specific galaxies.')
    parser.add_argument("INPUT_FILENAME", help="Filename of the input file." )
    parser.add_argument('IGALS', help='Galaxy number(s) (i.e. igals). '
                        'If using list of igals, it must be comma separated'
                        ' (with no spaces).')
    parser.add_argument('-q', "--quantities",
                   help='List of quantities to be plotted (comma separated).',
                   default=None )
    parser.add_argument('-z', "--use_redshift", help="Keeps galaxies.hdf5 files.", action="store_true")

    args = parser.parse_args()

    input_filename = args.INPUT_FILENAME
    igals = N.array(args.IGALS.split(',')).astype(int)
    quantities = args.quantities
    use_redshift = args.use_redshift

    f = h5py.File(input_filename,'r')
    Input = f['Input']

    if quantities == None:
        quantities = Input.keys()
        quantities.remove('t')
        if 'r_max' in quantities:
          quantities.remove('r_max')
        quantities.remove('z')
    else:
        quantities = quantities.split(',')

    if use_redshift:
        t = Input['z']
    else:
        t = Input['t']

    for q in quantities:
        if q not in Input:
            print q, 'not available'
            print Input.keys()
            continue

        P.figure()
        for igal in igals:
            data = N.array(Input[q][igal,:])

            ok = data>-1e-5
            P.plot(t[ok], data[ok],marker='.')
        if use_redshift:
            P.xlabel(r'$z$')
        else:
            P.xlabel(r'$t\,[{\rm Gyr}]$')
        P.ylabel(q)

    P.show()
