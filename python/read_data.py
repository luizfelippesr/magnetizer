""" Contains routines which allow extracting galform data from an output hdf5
    file and repack it in a more usable way. """
import numpy as np
import h5py
import os
from numpy.random import random, shuffle

def add_or_create(idx, quantity_name, dictionary, quantity, size,
                  fill_value=np.NaN):
    if quantity_name not in dictionary:
        dictionary[quantity_name] = np.full(size, fill_value)

    dictionary[quantity_name][idx] = quantity

    return

def read_time_data(sam_output_filepath, maximum_final_B_over_T=0.5,
                   max_z = 1000, minimum_final_stellar_mass=1e6,
                   maximum_final_stellar_mass=1e14, minimum_final_gas_mass=1e5,
                   number_of_galaxies=100, empirical_disks=True,
                   minimum_final_disk_size=0.0, datasets=None,
                   sample_all_z=True, fill_value=-9999999.,
                   number_of_galaxies_high_z=100):
    """
    Reads data from the galaxies.hdf5 file inside sam_output_filepath.
    The option number_of_galaxies sets the _approximate_ number of galaxies
    in the output (which are drawn, randomly, from the sample, taking into
    account halo Press-Schechter weights.

    Returns a dictionary with the structure:
        { Galaxy_ID : { dataset_name : time_series_array }}
    for all galaxies/datasets.

    Missing items in the time series are replaced by the value fill_value.
    """

    if datasets == None:
        # Which data sets will be read from the output
        datasets = (
                    # Essential physical properties
                    'vdisk', 'mstars_disk', 'mstars_bulge', 'mcold', 'mhot',
                    'mcold_burst','rdisk', 'rbulge', 'mstardot', 'mchalo',
                    'is_central', 'mstardot_average',
                    'vbulge','vhalo',
                    'halo_r_virial',
                    'strc',
                    # New style indexing and extras
                    #'FirstProgenitorID',
                    #'EndMainBranchID','DescendantID',
                    'GalaxyID',
                    #'LastProgenitorID',
                    'tburst',
                    'galaxy_weight'
                   )

    # Opens the hdf5 file
    f = h5py.File(sam_output_filepath)

    # Reads the little-h
    h0 = f['Parameters']['h0'][()]

    # Gets the redshifts and times
    zout_array = f["Output_Times/zout"][:]
    tout_array = f["Output_Times/tout"][:]
    z_indices = np.arange(len(zout_array[zout_array<=max_z]))
    nzs = z_indices.size

    # Initializes dictionary and copy galform run parameters
    data_dict = {'Galform_Parameters': {}, 'h0' : h0}
    for param in f['Parameters']:
        data_dict['Galform_Parameters'][param] = f['Parameters'][param][()]

    target_IDs = dict()
    IDs = []

    # Loops over redshifs
    for zidx in z_indices[:]:
        save_idx = nzs-zidx-1 # Stores earlier times first
        t = tout_array[zidx]
        add_or_create(save_idx, 't', data_dict, t, nzs, fill_value=fill_value)
        add_or_create(save_idx, 'z', data_dict, zout_array[zidx], nzs,
                      fill_value=fill_value)

        output_id = "Output%.3i" % (zidx+1)   # outputs labeled Fortran style.
        fout = f[output_id]
        print('Working on redshift {0}, time {1}'.format(zout_array[zidx],t))


        # If galaxy weights are not in the file, include them
        if "galaxy_weight" not in f[output_id]:
            print('Adding galaxy weights (this may take some time).')
            # Weights are associated with haloes.
            # Need to assign a weight to each galaxy.
            # The same goes with jm and jtree
            halo_weight = f[output_id+"/Trees/weight"]
            halo_ngals = f[output_id+"/Trees/ngals"]
            #halo_jm = f[output_id+"/Trees/jm"]
            #halo_jtree = f[output_id+"/Trees/jtree"]

            if ((np.sum(halo_ngals) != fout['mstars_disk'].size) or
                 (len(halo_weight) != len(halo_ngals))):
                print "ERROR: not enough galaxies in the file"
                exit()

            # count galaxies, assiging weight from their halo weight.
            # initial weights are -1
            galaxy_weight = np.zeros( len(f[output_id+"/mstars_disk"][:]) )-1.
            #jtree = np.zeros( len(f[output_id+"/mstars_disk"][:]) ) -1.
            #jm = np.zeros( len(f[output_id+"/mstars_disk"][:]) ) -1.

            igal = 0
            for ihalo in range( len(halo_ngals) ):
                if (halo_ngals[ihalo] > 0):
                    galaxy_weight[igal:igal+halo_ngals[ihalo]] = halo_weight[ihalo]
                    #jtree[igal:igal+halo_ngals[ihalo]] = halo_jtree[ihalo]
                    #jm[igal:igal+halo_ngals[ihalo]] = halo_jm[ihalo]
                    igal = igal+halo_ngals[ihalo]

            fout["galaxy_weight"] = galaxy_weight
            #f[output_id+"/galaxy_jtree"] = jtree
            #f[output_id+"/galaxy_jm"] = jm
            f.flush()

        # First goes through all the previously selected galaxies
        # and stores the indices
        old_indices = []
        GalaxyIDs = f[output_id]['GalaxyID'][:]

        IDs_with_target_IDs = target_IDs.keys()

        for ID in IDs_with_target_IDs:
            old = np.where(GalaxyIDs==target_IDs[ID])[0]
            if old.size!=0:
                # If the target galaxy was found in the present redshift
                idx = old[0]
                # Saves the index for later
                old_indices.append(idx)
                # Stores the ID of the galaxy in the next earlier redshift
                # (actually, the "most massive progenitor")
                target_IDs[ID] = fout['FirstProgenitorID'][idx]
                # Adds data to dictionary
                for dataset in datasets:
                    if dataset not in fout: continue
                    add_or_create(idx=save_idx,
                                  quantity_name=dataset,
                                  dictionary=data_dict[ID],
                                  quantity=fout[dataset][idx],
                                  size=nzs,
                                  fill_value=fill_value)
            else:
                # Otherwise, remove galaxy from target dictionary
                del target_IDs[ID]

        if (zidx == 0) or (sample_all_z and zidx!=z_indices[-1]):
            if zidx != 0 :
                number_of_galaxies = number_of_galaxies_high_z
            print number_of_galaxies
            # Pre-loads parts of the HDF5 file to RAM
            # (this converts the 'HDF5 dataset' into a numpy array)
            mstars_disk = fout['mstars_disk'][:]
            mcold = fout['mcold'][:]
            mstars_bulge = fout['mstars_bulge'][:]
            mcold_burst = fout['mcold_burst'][:]
            # Computes bulge to total mass ratio
            BoT = (mstars_bulge + mcold_burst)/(mstars_disk + mcold
                                                + mstars_bulge + mcold_burst)
            # Selects galaxies with a minimum stellar mass, gass mass
            # and the correct bulge to total mass ratio
            ok = mstars_disk > minimum_final_stellar_mass
            ok *= mstars_disk < maximum_final_stellar_mass
            ok *= mcold > minimum_final_gas_mass
            ok *= BoT < maximum_final_B_over_T

            # Constructs an array of valid indices
            new_indices= np.where(ok)[0]

            # If a maximum number of galaxies was specified
            if number_of_galaxies:
                weight = fout['galaxy_weight'][ok]
                randm = random(weight.size)
                # Selects a sample of galaxies based on weight
                rfac = number_of_galaxies / weight.sum()

                new_indices= new_indices[randm/weight<rfac]

                if new_indices.size>number_of_galaxies:
                    shuffle(new_indices)
                    new_indices= new_indices[:number_of_galaxies]

            for idx in new_indices:
                # If a pre-selected index was selected, one has already dealt
                # with it. Thus, skips it.
                if idx in old_indices:
                    continue

                ID = GalaxyIDs[idx]
                IDs.append(ID)

                # Stores the ID of the galaxy in the next earlier redshift
                # (actually, the "most massive progenitor")
                target_IDs[ID] = fout['FirstProgenitorID'][idx]

                # Initializes dictionary
                data_dict[ID] = dict()

                # Adds data to dictionary
                for dataset in datasets:
                    if dataset not in fout: continue
                    add_or_create(idx=save_idx,
                                  quantity_name=dataset,
                                  dictionary=data_dict[ID],
                                  quantity=fout[dataset][idx],
                                  size=nzs,
                                  fill_value=fill_value)
    f.close()
    data_dict['IDs'] = sorted(IDs)
    return data_dict



def plot_mass_evolution(sam_output_filepath, gtype='all'):
    """ Tests the module ploting the time evolution of galaxy mass.
        gtype can be either: all, central/cen or sat/satellite."""
    import pylab as P

    # Will color them by mass
    cmap = P.cm.get_cmap('YlGnBu',200)

    data_dict = read_time_data(sam_output_filepath,
                               maximum_final_B_over_T=1.0,
                               minimum_final_stellar_mass=1e8,
                               minimum_final_gas_mass=1e7,
                               number_of_galaxies=20,
                               empirical_disks=False,
                               datasets = ('mstars_disk', 'mstars_bulge',
                                           'mcold', 'mhot',
                                           'mcold_burst'),
                               number_of_galaxies_high_z = 5
                               )
    max_mass = None
    for ID in data_dict:
        if ID in ['h0','IDs','t','z','Galform_Parameters']:
            continue

        mass  = (data_dict[ID]['mstars_disk'] + data_dict[ID]['mstars_bulge'])

        this_max_mass = mass[mass>0].max()
        if max_mass is None:
            max_mass = this_max_mass
            min_mass = this_max_mass
        max_mass = max(this_max_mass, max_mass)
        min_mass = min(this_max_mass, min_mass)

        c=cmap((this_max_mass-min_mass)/(max_mass-min_mass) )
        P.plot(data_dict['t'],np.array(mass)/1e10, marker='.', color=c)
    P.xlabel(r'$t\,[{{\rm Gyr}}]$')
    P.ylabel(r'$M\,/\, (10^{10}\,{{\rm M}}_\odot)$')
    P.show()


if __name__ == "__main__"  :
    model_dir = 'scripts/test_SAM_output/galaxies.hdf5'
    model_dir = '/data/nlfsr/galform_models/GON/GON-partial-ivolume2/galaxies.hdf5'

    plot_mass_evolution(model_dir, gtype='central')
