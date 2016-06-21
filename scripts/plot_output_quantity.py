#!/usr/bin/env python
import h5py
import sys
import numpy as N
import pylab as P
import argparse
import portfolio as pf

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
            data_dict[q] = Output[q]

    if use_redshift:
        t = Input['z']
    else:
        t = Input['t']

    for quantity in quantities:
        for i, igal in enumerate(igals):
            if i==0:
                fig = P.figure(1)
            else:
                fig = P.figure()
            ax = fig.add_subplot(111)
            P.title('Galaxy {0}'.format(igal+1))
            pf.plot_quantity(igal, quantity, data_dict, cmap=P.cm.viridis,
                          ax=ax)

            norm = P.mpl.colors.Normalize(vmin=data_dict['t'][:].min(),
                                          vmax=data_dict['t'][:].max())
            ax = fig.add_axes([.065, .04, .9, .01])
            x = P.mpl.colorbar.ColorbarBase(ax, cmap=P.cm.YlGnBu, norm=norm,
                                            orientation='horizontal')


    P.show()
