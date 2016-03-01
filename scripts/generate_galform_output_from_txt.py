import numpy as N
import h5py

columns = ('mcold', 'mstars_disk', 'rdisk', 'vdisk','mstardot', 'vbulge', 'rbulge', 'vchalo', 'halo_r_virial', 'strc')

def add_dataset(storage, dataset_name, dataset, compression=False):
    """ Concatenates a matrix (dataset) called dataset_name to an openned HDF5 file (storage), creating if necessary """

    # If the dataset doesn't exist, create it!
    if dataset_name not in storage:
        shape = list(N.shape(dataset))
        shape[0]=None
        if compression:
            storage.create_dataset(dataset_name, data=dataset, maxshape=shape, compression='gzip')
        else:
            storage.create_dataset(dataset_name, data=dataset, maxshape=shape)
    # Otherwise, reshapes the dataset in the HDF5 storage and add the new data
    else:
        old_length = storage[dataset_name].shape[0]
        storage[dataset_name].resize(old_length+len(dataset), axis=0)

        storage[dataset_name][old_length:] = dataset

    return


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

def make_galaxies_file(data_dict, t0, tf, nt,
                       output_filename='test/ivol0/galaxies.hdf5'):
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
    zout_array = tout_array*N.NaN # This is not really being used...

    params = f.create_group('Parameters')
    params['h0'] = 1.0

    Output_Times = f.create_group('Output_Times')
    Output_Times["zout"] = zout_array
    Output_Times["tout"] = tout_array
    z_indices = N.arange(len(zout_array))

    gals = sorted(data_dict.keys())
    gals_ID = N.arange(len(gals))
    for gal in gals:
      add_dataset(f, 'Names', N.array([gal,]))

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

            # Saves all other physical properties
            for c in columns:
                if c in ('mcold', 'mstars_disk', 'mstardot'):
                    continue
                add_dataset(output, c, N.array([data_dict[gal][c],]))
            # Saves the galaxy name on all ID-related fields
            for c in ('FirstProgenitorID','EndMainBranchID','DescendantID',
                      'GalaxyID','LastProgenitorID'):
                add_dataset(output, c, N.array([ID,]))

            # Manual adjustments
            add_dataset(output, 'mstardot_average', N.array([data_dict[gal]['mstardot'],]))
            add_dataset(output, 'mhot', N.array([0.,]))
            add_dataset(output, 'mcold_burst', N.array([0.,]))
            add_dataset(output, 'mstars_bulge', N.array([0.,]))
            add_dataset(output, 'galaxy_weight', N.array([1.,]))
            add_dataset(output, 'is_central', N.array([1.,]))
    f.close()

a = read_parameters('example_galaxies.txt')
make_galaxies_file(a, 10.0, 13.5, 10)
