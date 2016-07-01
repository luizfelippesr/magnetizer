#!/usr/bin/env python
import h5py
import sys
import numpy as N
import pylab as P
import argparse
import portfolio as pf
from extra_quantities import compute_extra_quantity

if __name__ == "__main__"  :
    parser = argparse.ArgumentParser(
             description='This script plots the outputs of specific galaxies.')
    parser.add_argument("INPUT_FILENAME", help="Filename of the input file." )
    parser.add_argument('IGALS', help='Galaxy number(s) (i.e. igals). '
                        'If using list of igals, it must be comma separated'
                        ' (with no spaces).')
    parser.add_argument('-q', "--quantities", default=None,
                   help='List of quantities to be plotted (comma separated).')
    parser.add_argument('-z', "--use_redshift", action="store_true",
                        help="Uses redshift instead of time.")

    args = parser.parse_args()

    input_filename = args.INPUT_FILENAME
    igals = N.array(args.IGALS.split(',')).astype(int)
    quantities = args.quantities
    use_redshift = args.use_redshift

    f = h5py.File(input_filename,'r')
    Output = f['Output']
    Input = f['Input']

    data_dict = {'t' : Input['t'],
                'r' : Output['r'],
                'Omega' : Output['Omega']}

    if quantities == None:
        exit('erro')
    else:
        quantities = quantities.split(',')
        for q in quantities:
            if q in Output:
                data_dict[q] = Output[q]
            else:
                data_dict[q] = compute_extra_quantity(q, Output)

    if use_redshift:
        t = Input['z']
    else:
        t = Input['t']

    ts =t[:]
    #ts = t[(t[:]>=6.)*(t[:]<6.5)]
    #print ts
    #print t[:]
    #exit()
    for j, quantity in enumerate(quantities):
        for i, igal in enumerate(igals):
            if j+i==0:
                fig = P.figure(1)
            else:
                fig = P.figure()
            ax = fig.add_subplot(111)
            P.title('Galaxy {0}'.format(igal+1))
            pf.plot_quantity(igal, quantity, data_dict, cmap=P.cm.viridis,
                          ax=ax, ts=ts)

            if use_redshift:
                vmax=t[:].min(); vmin=t[:].max()
            else:
                vmin=ts[:].min(); vmax=ts[:].max()
            norm = P.mpl.colors.Normalize(vmin=vmin, vmax=vmax)
            ax = fig.add_axes([.065, .04, .9, .01])
            x = P.mpl.colorbar.ColorbarBase(ax, cmap=P.cm.viridis, norm=norm,
                                            orientation='horizontal')


    P.show()
