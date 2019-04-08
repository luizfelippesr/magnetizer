#!/usr/bin/env python
# Copyright (C) 2018,2019 Luiz Felippe S. Rodrigues, Luke Chamandy
#
# This file is part of Magnetizer.
#
# Magnetizer is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Magnetizer is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Magnetizer.  If not, see <http://www.gnu.org/licenses/>.
#
import numpy as N
import h5py
from hdf5_util import add_dataset
# Warning: order is important
columns = ('mcold', 'mstars_disk', 'rdisk', 'vdisk','mstardot', 'vbulge', 'rbulge', 'vhalo', 'mhalo', 'strc')

def read_parameters(filename):
    """ Reads the galaxy parameters from file
        Input: filename
        Output: dictionary of dictionaries containing all data
                example:   Mgas = d['M31']['Mgas'] """

    names = N.genfromtxt(filename, usecols=(0), dtype=('S5'))
    data = N.genfromtxt(filename)
    data = data[:,1:] # Removes names


    d = {}
    for i, name in enumerate(names):
        d[name] = dict()
        for j, c in enumerate(columns):
            d[name][c] = data[i,j]
    return d

def make_galaxies_file(data_dict, t0, tf, nt, output_filename):
    """ Prepares a galaxies.hdf5 file from a dictionary of galaxy
        data. Will assume a constant SFR and compute the evolution
        of Mgas and Mstars accordingly.

        Input: data_dict -> dictionary with structure:
                            d[name] = {'Mgas':Mgas,
                                      'Mstars':Mstars,
                                      'r_disk':r_disk,
                                      'SFR':SFR}
               t0 -> initial snapshot time (earliest), in Gyr
               tf -> final snapshot time (present), in Gyr
               nt -> number of snapshots
               optional, output_filename -> name of the output filename
    """

    f = h5py.File(output_filename,'w')

    tout_array = N.linspace(tf, t0, nt)
    zout_array = tout_array*0.0 

    params = f.create_group('Parameters')
    params['h0'] = 1.0

    Output_Times = f.create_group('Output_Times')
    Output_Times["zout"] = zout_array
    Output_Times["tout"] = tout_array
    z_indices = N.arange(len(zout_array))

    gals = sorted(data_dict.keys())
    gals_ID = N.arange(len(gals))

    for i, zidx in enumerate(z_indices[:]):
        output_id = "Output%.3i" % (zidx+1) # outputs labeled Fortran style.
        output = f.create_group(output_id)

        for ID, gal in zip(gals_ID, gals):
            # Saves disk stellar and disk gas masses
            mcold_f = data_dict[gal]['mcold']
            mstars_f = data_dict[gal]['mstars_disk']
            sfr = data_dict[gal]['mstardot']*1e9 # Converts 1/yr -> 1/Gyr
            mstars = mstars_f - sfr*(tf - tout_array[i])
            mcold = mcold_f + sfr*(tf - tout_array[i])
            add_dataset(output, 'mstardot', N.array([sfr,]))
            add_dataset(output, 'mcold', N.array([mcold,]))
            add_dataset(output, 'mstars_disk', N.array([mstars,]))

            # Radii information has to have their units tweaked
            for rname in ('rdisk','rbulge'): #,'halo_r_virial'):
                r = data_dict[gal][rname]*1e-3 # kpc -> Gpc
                add_dataset(output, rname, N.array([r,]))

            # Saves all other physical properties
            for c in columns:
                if c in ('mcold', 'mstars_disk', 'mstardot','rdisk','rbulge','halo_r_virial'):
                    continue
                add_dataset(output, c, N.array([data_dict[gal][c],]))

            # Saves all ID-related fields
            for c in ('FirstProgenitorID','EndMainBranchID','DescendantID',
                      'GalaxyID','LastProgenitorID'):
                add_dataset(output, c, N.array([ID,]))
            add_dataset(output, 'names', N.array([gal,]))

            # Manual adjustments
            add_dataset(output, 'mstardot_average', N.array([data_dict[gal]['mstardot'],]))
            add_dataset(output, 'mhot', N.array([0.,]))
            add_dataset(output, 'mcold_burst', N.array([0.,]))
            add_dataset(output, 'mstars_bulge', N.array([0.,]))
            add_dataset(output, 'galaxy_weight', N.array([1.,]))
            add_dataset(output, 'is_central', N.array([1.,]))
    f.close()


if __name__ == "__main__"  :
    import argparse

    parser = argparse.ArgumentParser(description='Prepares an fake Galform '
                            'output from a text table of galaxy properties')

    parser.add_argument("PROPERTIES_TABLE", help="File path to the properties"
                        " table.")

    parser.add_argument("FAKE_SAM_OUTPUT",help="File path of the fake Galform "
                        'HDF5 output file.')

    parser.add_argument('-ft', '--final_time', default=13.5,
                        help='Final comoving time for the fake cosmological'
                        ' run in Gyr. Default: 13.5.')
    parser.add_argument('-it', '--initial_time', default=10.0,
                        help='Initial comoving time for the fake cosmological'
                        ' run in Gyr. Default: 10.0.')
    parser.add_argument('-n', '--number_of_snapshots', default=10.0,
                        help='Number of snapshots in the fake cosmological'
                        ' run. Default: 10.')

    args = parser.parse_args()

    param_dict = read_parameters(args.PROPERTIES_TABLE)
    make_galaxies_file(param_dict, args.initial_time, args.final_time,
                       args.number_of_snapshots, args.FAKE_SAM_OUTPUT )
