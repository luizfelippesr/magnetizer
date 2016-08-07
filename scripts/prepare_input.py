#!/usr/bin/env python
""" Contains routines which prepare galaxy input files for the magnetic field
    evolution code. """
from read_data import read_time_data
import numpy as N

# Some constants
Mpc_to_kpc = 1000
G_SI=6.67259e-11
MSOLAR=1.9891e30 # The mass of the Sun in kg
KM2M=1.0e3 # km to m
MPC2M=3.0856775807e22 # The number of metres in a Mpc
G=G_SI*MSOLAR/MPC2M/(KM2M**2) # The gravitational constant
                              # in units of (km/s)^2 Mpc/Msun.
description_dictionary = {
        'ID': 'Galaxy ID or name.',
        'z':'Redshift',
        't':'Time since the Big Bang.',
        'SFR' : 'Star formation rate',
        'Mstars_disk':'Stelar mass of the galaxy disc',
        'Mstars_bulge':'Stelar mass of the galaxy bulge',
        'Mgas_disk': 'Total gas mass of the galaxy disc',
        'r_disk': 'Half mass radius of the galaxy disc',
        'v_disk': 'Circular velocity at the half-mass radius of the galaxy disc',
        'r_bulge': 'Half mass radius of the galaxy bulge',
        'v_bulge': 'Circular velocity at the half-mass radius of the galaxy bulge',
        'r_halo': 'Virial radius of the dark matter halo',
        'v_halo': 'Circular velocity at the virial radius of the dark matter halo',
        'nfw_cs1': 'Inverse of the NFW concentration parameter',
        'r_max': 'Maximum radius considered'
        }

units_dictionary = {
        't' :'Gyr',
        'SFR': 'Msun/yr',
        'Mstars_disk': 'Msun',
        'Mstars_bulge': 'Msun',
        'Mgas_disk': 'Msun',
        'r_disk' : 'kpc',
        'v_disk' : 'km/s',
        'r_bulge': 'kpc',
        'v_bulge': 'km/s',
        'r_halo' : 'kpc',
        'v_halo' : 'km/s',
        'r_max' : 'kpc'
        }

def prepares_hdf5_input(data_dict, output_file):
    from hdf5_util import add_dataset
    import h5py

    number_of_r50 = 2.5

    h5file = h5py.File(output_file)
    input_grp = h5file.create_group('Input')

    galform_grp = h5file.create_group('Galform Parameters')
    for param in data_dict['Galform Parameters']:
        galform_grp[param] = data_dict['Galform Parameters'][param]

    ts = data_dict['tout']
    zs = data_dict['zout']
    h0 = data_dict['h0']

    IDs = data_dict[ts[0]]['ID']

    datasets =('r_disk',
               'v_disk',
               'r_bulge',
               'v_bulge',
               'r_halo',
               'v_halo',
               'nfw_cs1',
               'Mgas_disk',
               'Mstars_disk',
               'Mstars_bulge',
               'SFR')

    if 'names' in data_dict:
        add_dataset(input_grp, 'ID', data_dict['names'])
    else:
        add_dataset(input_grp, 'ID', IDs)
    for i in N.argsort(ts):
        add_dataset(input_grp, 't', [ts[i],])
        add_dataset(input_grp, 'z', [zs[i],])

    for k in ('t','z','ID'):
        input_grp[k].attrs['Description'] = description_dictionary[k]
    input_grp['t'].attrs['Units'] = units_dictionary['t']

    for i, ID in enumerate(IDs):
        # Constructs a temporary dictionary to store the time series
        # extracted from galform
        tmp = dict()
        r_disk_max = 0
        # Initializes all arrays with the "INVALID" mark: -999999
        for d in datasets:
            tmp[d] = N.ones_like(ts)*(-999999)

        for j, t in enumerate(sorted(ts)):
            select = data_dict[t]['ID'] == ID
            # Skips missing times..
            if not select.any():
                continue

            # Loads what is necessary (and removes the little-h dependence)
            tmp['Mstars_disk'][j] = data_dict[t]['mstars_disk'][select][0]/h0

            tmp['Mgas_disk'][j] = data_dict[t]['mcold'][select][0]/h0
            tmp['Mstars_bulge'][j] = data_dict[t]['mstars_bulge'][select][0]/h0
            # NB mstars_bulge and mcold_burst not used

            tmp['SFR'][j] = data_dict[t]['mstardot'][select][0]
            tmp['SFR'][j] *= 1e-9 /h0 # Msun/Gyr/h -> Msun/yr

            r_disk = data_dict[t]['rdisk'][select][0]
            r_disk *= Mpc_to_kpc/h0

            tmp['r_disk'][j] = r_disk

            tmp['v_disk'][j] = data_dict[t]['vdisk'][select][0]

            tmp['r_bulge'][j] = data_dict[t]['rbulge'][select][0]
            tmp['r_bulge'][j] *= Mpc_to_kpc/h0

            tmp['v_bulge'][j] = data_dict[t]['vbulge'][select][0]

            tmp['v_halo'][j] = data_dict[t]['vhalo'][select][0]
            if 'halo_r_virial' in data_dict[t]:
                r_halo  = data_dict[t]['halo_r_virial'][select][0]
            else:
                mhalo = data_dict[t]['mhalo'][select][0]
                r_halo = G*mhalo/(tmp['v_halo'][j]**2)
            r_halo *= Mpc_to_kpc/h0
            tmp['r_halo'][j] = r_halo

            tmp['nfw_cs1'][j] = data_dict[t]['strc'][select][0]

        for dataset_name in tmp:
            add_dataset(input_grp, dataset_name, [tmp[dataset_name],])


        for dataset_name in tmp:
            input_grp[dataset_name].attrs['Description'] = \
                                          description_dictionary[dataset_name]
            if dataset_name in units_dictionary:
                input_grp[dataset_name].attrs['Units'] = \
                                                units_dictionary[dataset_name]
        ngals, nzs = input_grp['r_disk'].shape
        h5file.attrs['Number of galaxies'] = ngals
        h5file.attrs['Number of snapshots'] = nzs



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
                        'extract from the SAM_OUTPUT file. Default: 1e10.')

    parser.add_argument('-BoT', "--maximum_B_over_T", default=0.5,
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

    parser.add_argument('-z', "--max_redshift", default=6,
                        help="Maximum redshift to use. Default: 6.")

    args = parser.parse_args()

    start = time.time()
    data_dict = read_time_data( args.SAM_OUTPUT,
                                max_z = float(args.max_redshift),
                                maximum_final_B_over_T=float(args.maximum_B_over_T),
                                minimum_final_stellar_mass=float(args.minimum_stellar_mass),
                                maximum_final_stellar_mass=float(args.maximum_stellar_mass),
                                minimum_final_gas_mass=float(args.minimum_gas_mass),
                                number_of_galaxies=int(args.number_of_galaxies),
                                minimum_final_disk_size=1e-3*float(args.minimum_disk_size),
                                empirical_disks=False)
    middle = time.time()
    prepares_hdf5_input(data_dict, args.MAGNETIZER_INPUT)
    end = time.time()

    print 'read_time_data', middle-start,'s'
    print 'prepares_hdf5_input', end-middle,'s'
