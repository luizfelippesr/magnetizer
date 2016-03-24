""" Contains routines which prepare galaxy input files for the magnetic field
    evolution code. """
from read_data import read_time_data
import numpy as N
Mpc_to_kpc = 1000


description_dictionary = {
        't':'Time since the Big Bang.',
        'SFR' : 'Star formation rate',
        'Mstars_disk':'Stelar mass of the galaxy disc',
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
        'Mgas_disk': 'Msun',
        'r_disk' : 'kpc',
        'v_disk' : 'km/s',
        'r_bulge': 'kpc',
        'v_bulge': 'km/s',
        'r_halo' : 'kpc',
        'v_halo' : 'km/s',
        'r_max' : 'kpc'
        }

def prepares_text_input(model_dir, odir):

    number_of_r50 = 2.5

    data_dict = read_time_data( model_dir,
                                max_z = 10,
                                maximum_final_B_over_T=0.5,
                                minimum_final_stellar_mass=2e9,
                                maximum_final_stellar_mass=1e14,
                                minimum_final_gas_mass=1e7,
                                number_of_galaxies=100,
                                minimum_final_disk_size=5e-4, # Gpc, i.e. 0.5 kpc
                                empirical_disks=False,
                                ivol_dir='ivol0')
    ts = data_dict['tout']

    IDs = data_dict[ts[0]]['ID']
    names = data_dict[ts[0]]['names'].astype(str)

    for i, ID in enumerate(IDs):
        f_gal_id = '{0:8d}'.format(i+1).replace(' ','0')
        name_input_file = '{1}/name_{0}.in'.format(f_gal_id, odir)
        with open(name_input_file, 'w+') as f:
            f.write(names[i])


        t_dep_input_file = '{1}/time_dep_params_{0}.in'.format(f_gal_id, odir)
        t_indep_input_file='{1}/time_indep_params_{0}.in'.format(f_gal_id, odir)
        header =('time        |'
                 'r_disk_kpc  |'
                 'v_disk_kpc  |'
                 'r_bulge_kpc |'
                 'v_bulge_kpc |'
                 'r_halo_kpc  |'
                 'v_halo_kpc  |'
                 'nfw_cs1     |'
                 'Mgas_disk   |'
                 'Mstars_disk |'
                 'SFR\n')

        with open(t_dep_input_file, 'w+') as fdep:
            t_dep_input = header
            r_disk_max = 0
            for t in sorted(ts):
                #print 't = {0} Gyr'.format(t)
                select = data_dict[t]['ID'] == ID
                # Skips missing times..
                if not select.any():
                    continue

                # Loads what is necessary (and removes the little-h dependence)
                M_stars =data_dict[t]['mstars_disk'][select][0]/data_dict['h0']

                M_gas = data_dict[t]['mcold'][select][0]/data_dict['h0']

                SFR = data_dict[t]['mstardot'][select][0]
                SFR *= 1e-9 /data_dict['h0'] # Msun/Gyr/h -> Msun/yr

                r_disk = data_dict[t]['rdisk'][select][0]/data_dict['h0']
                r_disk *= Mpc_to_kpc/data_dict['h0']

                v_disk = data_dict[t]['vdisk'][select][0]
                # NB mstars_bulge and mcold_burst not considered
                #    they would need different profiles...

                r_disk_max = max(r_disk, r_disk_max)

                r_bulge = data_dict[t]['rbulge'][select][0]
                r_bulge *= Mpc_to_kpc/data_dict['h0']

                v_bulge = data_dict[t]['vbulge'][select][0]

                v_halo = data_dict[t]['vhalo'][select][0]
                if 'halo_r_virial' in data_dict[t]:
                    r_halo  = data_dict[t]['halo_r_virial'][select][0]
                else:
                    G_SI=6.67259e-11
                    MSOLAR=1.9891e30 # The mass of the Sun in kg
                    KM2M=1.0e3 # km to m
                    MPC2M=3.0856775807e22 # The number of metres in a Mpc
                    G=G_SI*MSOLAR/MPC2M/(KM2M**2) # The gravitational constant
                                                  # in units of (km/s)^2 Mpc/Msun.
                    mhalo = data_dict[t]['mhalo'][select][0]
                    r_halo = G*mhalo/(v_halo**2)
                r_halo *= Mpc_to_kpc/data_dict['h0']


                nfw_cs1 = data_dict[t]['strc'][select][0]


                t_dep_input += '{0:.3e}    '.format(t)
                t_dep_input += '{0:.3e}    '.format(r_disk)
                t_dep_input += '{0:.3e}    '.format(v_disk)
                t_dep_input += '{0:.3e}    '.format(r_bulge)
                t_dep_input += '{0:.3e}    '.format(v_bulge)
                t_dep_input += '{0:.3e}    '.format(r_halo)
                t_dep_input += '{0:.3e}    '.format(v_halo)
                t_dep_input += '{0:.3e}    '.format(nfw_cs1)
                t_dep_input += '{0:.3e}    '.format(M_gas)
                t_dep_input += '{0:.3e}    '.format(M_stars)
                t_dep_input += '{0:.3e}    '.format(SFR)
                t_dep_input += '\n'

            # Writes the file
            fdep.write(t_dep_input)
        with open(t_indep_input_file, 'w+') as findep:
            t_indep_input_file = 'r_max_kpc\n'
            t_indep_input_file += '{0:.3e}\n'.format(
                                             number_of_r50*r_disk_max)

            findep.write(t_indep_input_file)


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

    IDs = data_dict[ts[0]]['ID']
    names = data_dict[ts[0]]['names'].astype(str)

    datasets =(
               'r_disk',
               'v_disk',
               'r_bulge',
               'v_bulge',
               'r_halo',
               'v_halo',
               'nfw_cs1',
               'Mgas_disk',
               'Mstars_disk',
               'SFR')

    add_dataset(input_grp, 't', [sorted(ts),])
    input_grp['t'].attrs['Description'] = description_dictionary['t']
    input_grp['t'].attrs['Units'] = units_dictionary['t']


    for i, ID in enumerate(IDs):

        # Constructs a temporary dictionary to store the time series
        # extracted from galform
        tmp = dict()
        r_disk_max = 0
        for d in datasets:
            tmp[d] = N.ones_like(ts)*(-999999)
        #tmp['t']= N.ones_like(ts)*N.NaN

        for j, t in enumerate(sorted(ts)):
            select = data_dict[t]['ID'] == ID
            # Skips missing times..
            if not select.any():
                continue

            # Loads what is necessary (and removes the little-h dependence)
            tmp['Mstars_disk'][j] = data_dict[t]['mstars_disk'][select][0]/data_dict['h0']

            tmp['Mgas_disk'][j] = data_dict[t]['mcold'][select][0]/data_dict['h0']

            tmp['SFR'][j] = data_dict[t]['mstardot'][select][0]
            tmp['SFR'][j] *= 1e-9 /data_dict['h0'] # Msun/Gyr/h -> Msun/yr

            r_disk = data_dict[t]['rdisk'][select][0]/data_dict['h0']
            r_disk *= Mpc_to_kpc/data_dict['h0']

            tmp['r_disk'][j] = r_disk

            r_disk_max = max(r_disk, r_disk_max)

            tmp['v_disk'][j] = data_dict[t]['vdisk'][select][0]
            # NB mstars_bulge and mcold_burst not considered
            #    they would need different profiles...


            tmp['r_bulge'][j] = data_dict[t]['rbulge'][select][0]
            tmp['r_bulge'][j] *= Mpc_to_kpc/data_dict['h0']

            tmp['v_bulge'][j] = data_dict[t]['vbulge'][select][0]

            tmp['v_halo'][j] = data_dict[t]['vhalo'][select][0]
            if 'halo_r_virial' in data_dict[t]:
                r_halo  = data_dict[t]['halo_r_virial'][select][0]
            else:
                G_SI=6.67259e-11
                MSOLAR=1.9891e30 # The mass of the Sun in kg
                KM2M=1.0e3 # km to m
                MPC2M=3.0856775807e22 # The number of metres in a Mpc
                G=G_SI*MSOLAR/MPC2M/(KM2M**2) # The gravitational constant
                                              # in units of (km/s)^2 Mpc/Msun.
                mhalo = data_dict[t]['mhalo'][select][0]
                r_halo = G*mhalo/(tmp['v_halo'][j]**2)
            r_halo *= Mpc_to_kpc/data_dict['h0']
            tmp['r_halo'][j] = r_halo

            tmp['nfw_cs1'][j] = data_dict[t]['strc'][select][0]

        tmp['r_max'] = r_disk_max


        for dataset_name in tmp:
            add_dataset(input_grp, dataset_name, [tmp[dataset_name],])


        for dataset_name in tmp:
            input_grp[dataset_name].attrs['Description'] = \
                                          description_dictionary[dataset_name]
            if dataset_name in units_dictionary:
                input_grp[dataset_name].attrs['Units'] = \
                                                units_dictionary[dataset_name]



if __name__ == "__main__"  :
    import argparse

    model_dir = 'test'
    odir = '../input/real'

    model_dir = '/home/nlfsr/galform_models/GON'
    odir = '../input/GON/ivol0'

    data_dict = read_time_data( model_dir,
                                max_z = 10,
                                maximum_final_B_over_T=0.5,
                                minimum_final_stellar_mass=2e9,
                                maximum_final_stellar_mass=1e14,
                                minimum_final_gas_mass=1e7,
                                number_of_galaxies=10,
                                minimum_final_disk_size=5e-4, # Gpc, i.e. 0.5 kpc
                                empirical_disks=False,
                                ivol_dir='ivol1')

    #prepares_text_input(data_dict, odir)
    prepares_hdf5_input(data_dict, 'test.hdf5')
