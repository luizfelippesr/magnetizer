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
""" Contains routines which allow extracting galform data from an output hdf5
    file and repack it in a more usable way. """
import numpy as N
import h5py
import os
from numpy.random import random, shuffle
#from empirical_disks import empirical_disk # Not yet implemented
from data_dict_util import filter_dictionary_inplace, filter_dict

def read_time_data(sam_output_filepath, maximum_final_B_over_T=0.5,
                   max_z = 1000.0, minimum_final_stellar_mass=1e6,
                   maximum_final_stellar_mass=1e14, minimum_final_gas_mass=1e5,
                   number_of_galaxies=100, empirical_disks=True,
                   minimum_final_disk_size=0.0, datasets=None):
    """
    Reads data from the galaxies.hdf5 file inside sam_output_filepath.
    The option number_of_galaxies sets the _approximate_ number of galaxies
    in the output (which are drawn, randomly, from the sample, taking into
    account halo press-schechter weights. Other options allow to select
    by minimum mass of the final galaxy.

    Returns a dictionary of dictionaries containing times as keys in the
    first level, dataset names in the second level and the arrays
    associated with each dataset in the last level.
    Redshift and time information can be obtained from the special keys: zout
    and tout.
    """

    if datasets == None:
        # Which data sets will be read from the output
        datasets = (
                    # Essential physical properties
                    'vdisk', 'mstars_disk', 'mstars_bulge', 'mcold', 'mhot',
                    'mcold_burst','rdisk', 'rbulge', 'mstardot', 'mhalo',
                    'is_central', 'mstardot_average',
                    'vbulge','vhalo','halo_r_virial','strc',
                    # New style indexing and extras
                    'FirstProgenitorID',
                    'EndMainBranchID','DescendantID',
                    'GalaxyID', 'LastProgenitorID',
                    'galaxy_weight'
                   )

    # Opens the hdf5 file
    f = h5py.File(sam_output_filepath)

    # Reads the little-h
    h0 = f['Parameters']['h0'][()]

    # Gets the redshifts and times
    zout_array = f["Output_Times/zout"][:]
    tout_array = f["Output_Times/tout"][:]
    z_indices = N.arange(len(zout_array[zout_array<=max_z]))

    # Initializes dictionary and copy galform run parameters
    data_dict = {'Galform Parameters': dict()}
    for param in f['Parameters']:
        data_dict['Galform Parameters'][param] = f['Parameters'][param][()]

    # Loops over redshifs
    for i, zidx in enumerate(z_indices[:]):
        t = tout_array[zidx]

        print('Working on redshift {0}, time {1}'.format(zout_array[zidx],t))
        # Initializes dictionary
        data_dict[t] = {}
        output_id = "Output%.3i" % (zidx+1)   # outputs labeled Fortran style.

        # Reads the data from the hdf5 file.
        for k in datasets:
            if k in f[output_id]:
                data_dict[t][k] = f[output_id][k]

        # If galaxy weights are not in the file, include them
        if "galaxy_weight" not in f[output_id]:
            print('Adding weights, jm and jtree (this may take some time).')
            # Weights are associated with haloes.
            # Need to assign a weight to each galaxy.
            # The same goes with jm and jtree
            halo_weight = f[output_id+"/Trees/weight"]
            halo_ngals = f[output_id+"/Trees/ngals"]
            #halo_jm = f[output_id+"/Trees/jm"]
            #halo_jtree = f[output_id+"/Trees/jtree"]

            if ((N.sum(halo_ngals) != len(data_dict[t]['mstars_disk']) or
                 len(halo_weight) != len(halo_ngals))):
                print "ERROR: not enough galaxies in the file"
                exit()

            # count galaxies, assiging weight from their halo weight.
            # initial weights are -1
            galaxy_weight = N.zeros( len(f[output_id+"/mstars_disk"][:]) )-1.
            #jtree = N.zeros( len(f[output_id+"/mstars_disk"][:]) ) -1.
            #jm = N.zeros( len(f[output_id+"/mstars_disk"][:]) ) -1.

            igal = 0
            for ihalo in range( len(halo_ngals) ):
                if (halo_ngals[ihalo] > 0):
                    galaxy_weight[igal:igal+halo_ngals[ihalo]] = halo_weight[ihalo]
                    #jtree[igal:igal+halo_ngals[ihalo]] = halo_jtree[ihalo]
                    #jm[igal:igal+halo_ngals[ihalo]] = halo_jm[ihalo]
                    igal = igal+halo_ngals[ihalo]

            f[output_id+"/galaxy_weight"] = galaxy_weight
            #f[output_id+"/galaxy_jtree"] = jtree
            #f[output_id+"/galaxy_jm"] = jm
            f.flush()
        else:
            galaxy_weight = f[output_id+"/galaxy_weight"]

        data_dict[t]['weight'] = galaxy_weight

        print 'Galaxies read:', len(galaxy_weight)
        if i==0:
            # Pre-loads parts of the HDF5 file to RAM (allowing a ~20% speedup)
            # (this converts the 'HDF5 dataset' into a numpy array)
            data_dict[t]['mstars_disk'] = data_dict[t]['mstars_disk'][:]
            data_dict[t]['mcold'] = data_dict[t]['mcold'][:]
            data_dict[t]['mstars_bulge'] = data_dict[t]['mstars_bulge'][:]
            data_dict[t]['mcold_burst'] = data_dict[t]['mcold_burst'][:]

            # Computes bulge to total mass ratio
            data_dict[t]['BoT'] = (data_dict[t]['mstars_bulge'] +
                                   data_dict[t]['mcold_burst']
                                  ) / (data_dict[t]['mstars_disk'] +
                                       data_dict[t]['mcold'] +
                                       data_dict[t]['mstars_bulge'] +
                                       data_dict[t]['mcold_burst']  )

            # Selects galaxies with a minimum stellar mass, gass mass
            # and the correct bulge to total mass ratio
            ok = data_dict[t]['mstars_disk'] > minimum_final_stellar_mass
            ok *= data_dict[t]['mstars_disk'] < maximum_final_stellar_mass
            ok *= data_dict[t]['mcold'] > minimum_final_gas_mass
            ok *= data_dict[t]['rdisk'][:] > minimum_final_disk_size
            ok *= data_dict[t]['BoT'] < maximum_final_B_over_T
            # Filters the dictionary
            filter_dictionary_inplace(ok,data_dict[t])

            print('Number of galaxies after initial filtering: {0}'.format(
                                                  len(data_dict[t]['weight'])))

            # If a maximum number of galaxies was specified
            if number_of_galaxies:
                print("Sampling up to {0} galaxies based on weight".format(
                                                           number_of_galaxies))
                randm = random( len(data_dict[t]['weight']) )
                # Selects a sample of galaxies based on weight
                rfac = number_of_galaxies / N.sum(data_dict[t]['weight'])
                ok = ( randm / data_dict[t]['weight'] < rfac )
                print('Actual number of selected galaxies: {0}'.format(
                                                                  len(ok[ok])))
                filter_dictionary_inplace(ok,data_dict[t])

            # Stores indexing information and reference ID
            GalaxyID_target    = data_dict[t]['FirstProgenitorID']
            data_dict[t]['ID']= data_dict[t]['GalaxyID']
            previous_ID = data_dict[t]['ID']
            previous_t = t
        else:
            # Selects only the most massive ascendents of the gals in
            # the previous t
            filt = N.zeros(data_dict[t]['GalaxyID'].shape, dtype=bool)
            data_dict[t]['ID'] = N.empty(data_dict[t]['GalaxyID'].shape)*N.NaN

            for target_gid, ID in zip(GalaxyID_target,previous_ID):
                # Finds the index of a target galaxy
                match = target_gid == data_dict[t]['GalaxyID']

                if match.any():
                    # If there is a match, updates filter
                    filt +=  match
                    # Includes the z=0 ID in the dictionary
                    data_dict[t]['ID'][match] = ID

            # Applies the filter
            filter_dictionary_inplace(filt,data_dict[t])

            print('Number of galaxies after sampling: {0}'.format(
                                                      len(data_dict[t]['ID'])))

            # Creates a list of missing galaxies

            GalaxyID_target = data_dict[t]['FirstProgenitorID'][:]
            previous_ID = data_dict[t]['ID']
            previous_t = t

    if 'names' in f[output_id]:
        # If there is a list of galaxy names (as in the fake output) use it
        data_dict['names'] = f[output_id+'/names']

    if empirical_disks: # TODO not working yet
        raise NotImplementedError

    data_dict['tout'] = tout_array[zout_array<max_z]
    data_dict['zout'] = zout_array[zout_array<max_z]
    data_dict['h0'] = h0

    return data_dict


def plot_mass_evolution(sam_output_filepath, gtype='all'):
    """ Tests the module ploting the time evolution of galaxy mass.
        gtype can be either: all, central/cen or sat/satellite."""
    import pylab as P

    # Will color them by mass
    cmap = P.cm.get_cmap('YlGnBu',200)

    data_dict = read_time_data(sam_output_filepath,
                               maximum_final_B_over_T=1.0,
                               minimum_final_stellar_mass=1e10,
                               minimum_final_gas_mass=1e7,
                               number_of_galaxies=100,
                               empirical_disks=False)

    ts = sorted(data_dict['tout'])[::-1]
    IDs = data_dict[ts[0]]['ID']
    if gtype == 'all':
        pass
    elif 'sat' in gtype:
        IDs = IDs[data_dict[ts[0]]['is_central']==0]
    elif 'cen' in gtype:
        IDs = IDs[data_dict[ts[0]]['is_central']==1]
    else:
        exit('plot_mass_evolution: Error! Unrecongnized gtype.')

    shuffle(IDs)
    IDs = IDs[:10]

    gals_dict = dict()
    max_mass=None

    for ID in sorted(IDs):
        gals_dict[ID] = []
        for t in ts:
            select = data_dict[t]['ID'] == ID
            mass  = (
                    #data_dict[t]['mcold'][select] +
                    data_dict[t]['mstars_disk'][select] +
                    #data_dict[t]['mhot'][select] +
                      #data_dict[t]['mstars_bulge'][select]+
                    data_dict[t]['mcold_burst'][select])

            if len(mass)==1:
                if max_mass==None:
                    max_mass = mass[0]
                    min_mass = mass[0]
                max_mass = max(mass[0], max_mass)
                min_mass = min(mass[0], min_mass)
                gals_dict[ID].append(mass[0])
            else:
                gals_dict[ID].append(N.NaN)

    for ID in IDs:
        c=cmap( (gals_dict[ID][0]-min_mass)/(max_mass-min_mass) )
        P.plot(ts,N.array(gals_dict[ID])/1e10, marker='.', color=c)
    P.xlabel(r'$t\,[{{\rm Gyr}}]$')
    P.ylabel(r'$M\,/\, (10^{10}\,{{\rm M}}_\odot)$')
    P.show()


if __name__ == "__main__"  :
    model_dir = 'scripts/test_SAM_output/galaxies.hdf5'
    plot_mass_evolution(model_dir, gtype='cen')
