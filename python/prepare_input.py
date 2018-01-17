#!/usr/bin/env python
# Copyright (C) 2018  Luiz Felippe S. Rodrigues, Luke Chamandy
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
""" Contains routines which prepare galaxy input files for the magnetic field
    evolution code. """
from read_data import read_time_data
import numpy as np

# Some constants
Mpc_to_kpc = 1000
G_SI=6.67259e-11
MSOLAR=1.9891e30 # The mass of the Sun in kg
KM2M=1.0e3 # km to m
MPC2M=3.0856775807e22 # The number of metres in a Mpc
G=G_SI*MSOLAR/MPC2M/(KM2M**2) # The gravitational constant
                              # in units of (km/s)^2 Mpc/Msun.
description_dictionary = {
        'IDs': 'Galaxy ID or name.',
        'z':'Redshift',
        't':'Time since the Big Bang.',
        'SFR' : 'Star formation rate',
        'Mstars_disk':'Stelar mass of the galaxy disc',
        'Mstars_bulge':'Stelar mass of the galaxy bulge',
        'Mhalo': 'Halo mass',
        'Mgas_disk': 'Total gas mass of the galaxy disc',
        'r_disk': 'Half mass radius of the galaxy disc',
        'v_disk': 'Circular velocity at the half-mass radius of the galaxy disc',
        'r_bulge': 'Half mass radius of the galaxy bulge',
        'v_bulge': 'Circular velocity at the half-mass radius of the galaxy bulge',
        'r_halo': 'Virial radius of the dark matter halo',
        'v_halo': 'Circular velocity at the virial radius of the dark matter halo',
        'nfw_cs1': 'Inverse of the NFW concentration parameter',
        'weight': 'The density of this type of galaxy in the simulation.',
        'central': 'Whether the galaxy is a central (1) or a satellite (0).',
        'GalaxyID':'Together with the redshift, this can be used '
                    'to locate a galaxy within Galform\'s output.',
        'last_burst':'Time since the last burst of star formation (Gyr)',
        'sample_z':'Redshift used for sampling the galaxies.'
        }

dataset_names = {'rdisk' : 'r_disk',
                 'vdisk' : 'v_disk',
                 'rbulge' : 'r_bulge',
                 'vbulge' : 'v_bulge' ,
                 'halo_r_virial' : 'r_halo',
                 'vhalo' : 'v_halo',
                 'strc' : 'nfw_cs1',
                 'mcold' : 'Mgas_disk',
                 'mchalo' : 'Mhalo',
                 'mstars_disk' : 'Mstars_disk',
                 'mstars_bulge' : 'Mstars_bulge',
                 'mstardot_average' : 'SFR',
                 'galaxy_weight' : 'weight',
                 'is_central' : 'central',
                 'GalaxyID' :'GalaxyID',
                 'tburst':'last_burst'
                   }

units_dictionary = {
        't' :'Gyr',
        'SFR': 'Msun/yr',
        'Mstars_disk': 'Msun',
        'Mstars_bulge': 'Msun',
        'Mgas_disk': 'Msun',
        'Mhalo': 'Msun',
        'r_disk' : 'kpc',
        'v_disk' : 'km/s',
        'r_bulge': 'kpc',
        'v_bulge': 'km/s',
        'r_halo' : 'kpc',
        'v_halo' : 'km/s',
        'weight': 'Mpc^-3',
        'last_burst' :'Gyr'
        }

def prepares_hdf5_input(data_dict, output_file):
    from hdf5_util import add_dataset
    import h5py

    h5file = h5py.File(output_file)
    input_grp = h5file.create_group('Input')
    galform_grp = h5file.create_group('Galform Parameters')
    for param in data_dict['Galform_Parameters']:
        galform_grp[param] = data_dict['Galform_Parameters'][param]

    h0 = data_dict['h0']

    IDs = data_dict['IDs']

    for dataset in dataset_names:
        name = dataset_names[dataset]
        dset = []
        for ID in IDs:
            if dataset not in data_dict[ID]: continue

            dset.append(data_dict[ID][dataset])
        if len(dset)==0: continue
        dset = np.vstack(dset)
        # Loads what is necessary (and removes the little-h dependence)
        if name in ('Mstars_disk','Mgas_disk','Mstars_bulge', 'weight','Mhalo'):
            dset = dset/h0
        elif name == 'SFR':
            dset = dset*1e-9/h0 # Msun/Gyr/h -> Msun/yr
        elif name in ('r_disk', 'r_bulge', 'r_halo'):
            dset = dset*Mpc_to_kpc/h0
        # Adds it to the
        add_dataset(input_grp, name, dset)

        input_grp[name].attrs['Description'] = description_dictionary[name]
        if name in units_dictionary:
            input_grp[name].attrs['Units'] = units_dictionary[name]

    if 'r_halo' not in input_grp:
        r_halo = G*input_grp['Mhalo'][...]/(input_grp['v_halo'][...]**2)
        r_halo *= Mpc_to_kpc
        add_dataset(input_grp, 'r_halo', r_halo)
        input_grp['r_halo'].attrs['Description'] = description_dictionary['r_halo']
        if name in units_dictionary:
            input_grp['r_halo'].attrs['Units'] = units_dictionary['r_halo']

    for name in ('t', 'z', 'sample_z', 'IDs'):
        add_dataset(input_grp, name, data_dict[name])
        input_grp[name].attrs['Description'] = description_dictionary[name]
        if name in units_dictionary:
            input_grp[name].attrs['Units'] = units_dictionary[name]


    ngals, nzs = input_grp['r_disk'].shape
    h5file.attrs['Number of galaxies'] = ngals
    h5file.attrs['Number of snapshots'] = nzs

    h5file.close()


if __name__ == "__main__"  :
    import time
    import argparse

    parser = argparse.ArgumentParser(
             description='Prepares an input file for the Magnetizer.')

    parser.add_argument("SAM_OUTPUT",help="HDF5 output file of a Galform run.")

    parser.add_argument("MAGNETIZER_INPUT", help="Name of the Magnetizer input"
                        " file to be prepared.")

    parser.add_argument('-n', '--number_of_galaxies', default=1e10,
                        help='Approximate *maximum* number of galaxies to '
                        'extract from the SAM_OUTPUT file at z=0. '
                        'The evolution of these galaxies will then be followed.'
                        'Default: 1e10.')

    parser.add_argument('-nz', '--number_of_extra_galaxies_per_z',
                        default=None,
                        help='Approximate *maximum* number of galaxies to '
                        'extract from the SAM_OUTPUT file for each redshift'
                        '(except z=0). Default: same as number_of_galaxies.')

    parser.add_argument('-BoT', "--maximum_B_over_T", default=0.75,
                        help='Maximum bulge to total mass ratio.'
                        ' Default: 0.5.')

    parser.add_argument('-ms', "--minimum_stellar_mass", default=1e7,
                        help="Minimum disk stellar mass at z=0 (in Msun)."
                        " Default: 1e7.")

    parser.add_argument('-mg', "--minimum_gas_mass", default=1e6,
                        help="Minimum disk gas mass at z=0 (in Msun)."
                        " Default: 1e6.")

    parser.add_argument('-Ms', "--maximum_stellar_mass", default=1e14,
                        help="Maximum disk stellar mass at z=0 (in Msun)."
                        " Default: 1e14.")

    parser.add_argument('-r', "--minimum_disk_size", default=0.5,
                        help="Minimum disk half mass radius at z=0 (in kpc)."
                        " Default: 0.5.")

    parser.add_argument('-z', "--max_redshift", default=4,
                        help="Maximum redshift to use. Default: 6.")

    parser.add_argument('-naz', "--do_not_sample_all_z", action="store_true",
                        help="If present, galaxies will be sampled only at "
                        "z=0 and their ancestors will be followed to high z "
                        "(this is equivalent to -nz 0). "
                        "If absent, extra galaxies will me sampled for each redshift following -nz.")


    args = parser.parse_args()

    start = time.time()

    nz = args.number_of_extra_galaxies_per_z
    if nz is None:
        nz = int(float(args.number_of_galaxies))
    else:
        nz = int(float(nz))

    data_dict = read_time_data(args.SAM_OUTPUT,
                               max_z = float(args.max_redshift),
                               maximum_final_B_over_T=float(args.maximum_B_over_T),
                               minimum_final_stellar_mass=float(args.minimum_stellar_mass),
                               maximum_final_stellar_mass=float(args.maximum_stellar_mass),
                               minimum_final_gas_mass=float(args.minimum_gas_mass),
                               number_of_galaxies=int(float(args.number_of_galaxies)),
                               minimum_final_disk_size=1e-3*float(args.minimum_disk_size),
                               empirical_disks=False,
                               sample_all_z=not args.do_not_sample_all_z,
                               number_of_galaxies_high_z=nz
                               )
    end = time.time()
    prepares_hdf5_input(data_dict, args.MAGNETIZER_INPUT)
    print 'preparation time', end-start,'s'
